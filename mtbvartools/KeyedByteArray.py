import zlib, os, pickle, tqdm
from collections.abc import Iterable
import numpy as np
from .dasktools import findClient

# Require Python Version >= 3.7 for ordered dictionaries

class _iLocKeyedByteArray:
    def __init__(self, kba):
        self.kba = kba

    def __getitem__(self, index):
        def parseItem(index):
            if isinstance(index, int):
                return self.kba.read(self.kba.row[index])
            elif isinstance(index, slice):
                start = index.start
                stop = index.stop
                if start == None:
                    start = 0
                if stop == None:
                    stop = len(self.kba.row) + 1
                return np.asarray([self.kba.read(r) for r in self.kba.row[start:stop]])
            elif isinstance(index, Iterable):
                return np.asarray([self.kba.read(r) for r in self.kba.row[index]])
            raise IndexError(f'Unhandled indexing input {index}.')
        if isinstance(index, tuple):
            if len(index) == 2:
                return parseItem(index[0])[:, index[1]]
            else:
                raise NotImplemented('Only len == 2 tuples allowed for indexing')
        else:
            return parseItem(index)

class _locKeyedByteArray:
    def __init__(self, kba):
        self.kba = kba

    def __getitem__(self, index):
        if isinstance(index, list):  # assume list of keys
            return np.asarray([self.kba.read(k) for k in index])
        else:  # assume key
            return self.kba.read(index)


class KeyedByteArray():
    def write(self, key, value):
        if self.file.mode == 'wb':
            self._checkKey(key)
            raw_bytes = bytearray(  # enforce dtype and write to byte array
                np.asarray(value, dtype=self.dtype))
            compressed_bytes = self.compress(raw_bytes, **self.compress_kwargs)
            start_pointer = self.file.tell()
            array_len = self.file.write(compressed_bytes)
            self.pointers[key] = (start_pointer, array_len)
            self.bytes += len(raw_bytes)
            self.bytes_written += array_len
        else:
            raise ValueError(f'Function not supported for mode: {self.file.mode}.')

    def writebinary(self, keys, inputbytes, raw_bytes=0):
        if self.file.mode == 'wb':
            start_pointer = self.file.tell()
            for key, value in keys.items():
                self._checkKey(key)
                self.pointers[key] = (value[0] + start_pointer, value[1])
            # no conversion directly writes binaryvalue
            array_len = self.file.write(inputbytes)
            self.bytes += raw_bytes  # must be user provided, otherwise do not change
            self.bytes_written += array_len
        else:
            raise ValueError(f'Function not supported for mode: {self.file.mode}.')

    @property
    def iloc(self):
        return _iLocKeyedByteArray(self)

    @property
    def loc(self):
        return _locKeyedByteArray(self)
        
    def readFile(self, key):
        pointer, streamlen = self.row_dict[key]
        self.file.seek(pointer)
        return np.frombuffer(
            self.decompress(self.file.read(streamlen)), dtype=self.dtype)
    
    def readMem(self, key):
        pointer, streamlen = self.row_dict[key]
        return np.frombuffer(
            self.decompress(self.stream[pointer:pointer + streamlen]), dtype=self.dtype)
        
    def close(self):
        if self.file.mode == 'wb':
            self._writeIndex()
            del self.index, self.pointers
        elif self.file.mode == 'rb':
            del self.index, self.pointers, self.row, self.row_dict, self.col, self.col_dict
        self.file.close()

    def rechunk(self, output_path, jobs=100):
        # define job
        def rechunkJob(origin_path, tmp_path, col_idxs):
            # open origin and tmp kba
            origin_kba = KeyedByteArray(origin_path, mode='r')
            new_rows = origin_kba.col[col_idxs[0]:col_idxs[1]]
            tmp_kba = KeyedByteArray(
                tmp_path, mode='w', dtype=origin_kba.dtype, columns=origin_kba.row,
                compression=origin_kba.index['compression']['compression_type'],
                compression_kwargs=origin_kba.index['compression']['kwargs'])
            # store part of the full array uncompressed in memory
            tmp_arr = np.empty(
                shape=(len(origin_kba.row), len(new_rows)),
                dtype=origin_kba.dtype)
            for i, row_key in enumerate(origin_kba.row):
                tmp_arr[i, :] = origin_kba.read(row_key)[col_idxs[0]:col_idxs[1]]
            origin_kba.close()
            # write to the tmp kba
            for i, row in enumerate(new_rows):
                tmp_kba.write(row, tmp_arr[:, i])
            tmp_kba.close()
            return tmp_path
        # prepare to do jobs
        client = findClient()
        job_idxs = list(zip(
            np.linspace(0, len(self.col) + 1, jobs + 1).astype(int)[:-1],
            np.linspace(0, len(self.col) + 1, jobs + 1).astype(int)[1:]))
        # submit futures to the client
        futures = []
        for i, jidxs in enumerate(job_idxs):
            tmp_filepath = f'{os.path.dirname(self.file.name)}/tmp.{i}.kba'
            futures.append(client.submit(
                rechunkJob, self.file.name, tmp_filepath, jidxs,
                priority=len(job_idxs)-i))
        # read full files, write to output
        print('rechunking compression (for off axis random access)....')
        output_kba = KeyedByteArray(
                output_path, mode='w', dtype=self.dtype, columns=self.row,
                compression=self.index['compression']['compression_type'],
                compression_kwargs=self.index['compression']['kwargs'])
        for f in tqdm.tqdm(futures):
            tmp_path = client.gather(f)
            tmp_kba = KeyedByteArray(  # read entire stream into memory
                tmp_path, mode='r', mem=True)
            output_kba.writebinary(
                tmp_kba.row_dict, tmp_kba.stream, raw_bytes=tmp_kba.index['compression']['bytes'])
            tmp_kba.close()
            client.cancel(f)
            del tmp_kba
            os.remove(tmp_path), os.remove(f'{tmp_path}.kbi')
        output_kba.close()

    def _checkKey(self, key):
        try:
            # key will not be unique, close file and return error
            existing_pointer = self.pointers[key]
            pointer = self.file.tell()
            self.close()
            raise ValueError(f'Cannot add at {pointer}!\nKey {key} already exists, pointing to {existing_pointer}.')
        except KeyError:
            pass  # keys should be unique, do nothing

    def _writeIndex(self):
        if self.file.mode == 'wb':
            self.index['compression']['dtype'] = self.dtype
            self.index['compression']['bytes'] = self.bytes
            self.index['compression']['bytes_written'] = self.bytes_written
            rows = []
            pointers = []
            for k, v in self.pointers.items():
                rows.append(k)
                pointers.append(v)
            self.index['indexing']['rows'] = rows
            self.index['indexing']['pointers'] = pointers
            with open(f'{self.file.name}.kbi', 'wb') as index_file:
                index_file.write(
                    zlib.compress(pickle.dumps(self.index)))
        else:
            raise ValueError(f'Function not supported for mode: {self.file.mode}.')

    def _initCompression(self, compression, **kwargs):
        # define compression and decompression functions
        self.compress_kwargs = kwargs
        if compression == 'zlib':
            self.compress = zlib.compress
            self.decompress = zlib.decompress
        elif compression == 'None':
            self.compress = lambda x: x
            self.decompress = lambda x: x
        else:
            raise NotImplementedError(f"Compression function {compression} not implemented.")

    def __init__(self, filepath, mode='r', columns=None, dtype=None, compression='zlib', compression_kwargs={}, mem=False):
        if mode == 'w':
            if columns is None:
                raise ValueError('Must define column names.')
            if dtype is None:
                raise ValueError('Must define dtype.')
            self.file = open(filepath, mode='wb')
            self._initCompression(
                compression, **compression_kwargs)
            # initialize write-specific attributes
            self.dtype = dtype
            self.bytes = 0
            self.bytes_written = 0
            self.pointers = {}
            self.index = {
                'compression': {
                    'compression_type': compression,
                    'kwargs': compression_kwargs},
                'indexing': {
                    'columns': columns}}
        elif mode == 'r':
            self.file = open(filepath, mode='rb')
            # load index
            with open(f'{self.file.name}.kbi', 'rb') as index_file:
                self.index = pickle.loads(zlib.decompress((index_file.read())))
            self._initCompression(
                self.index['compression']['compression_type'], **self.index['compression']['kwargs'])
            self.dtype = self.index['compression']['dtype']
            # initialize read-specific attributes
            self.pointers = np.asarray(
                self.index['indexing']['pointers'], dtype=int)
            self.row = self.index['indexing']['rows']
            self.row_dict = {  # maps row names to stream pointers
                r: p for r,p in zip(self.row, self.pointers)}
            self.col = self.index['indexing']['columns']
            self.col_dict = {  # maps column names to numerical indexes
                v: i for i,v in enumerate(self.col)}
            # define read function
            if mem == True:
                self.stream = self.file.read()  # read entire stream into memory
                self.read = self.readMem
            else:
                self.read = self.readFile
        else:
            raise NotImplementedError("Supported modes are 'r' and 'w'")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()