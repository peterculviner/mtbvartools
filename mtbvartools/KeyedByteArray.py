import zlib, os, pickle, tqdm, numbers
from collections.abc import Iterable
import numpy as np
from .dasktools import findClient, subproc

# Require Python Version >= 3.7 for ordered dictionaries

class _iLocKeyedByteArray:
    def __init__(self, kba):
        self.kba = kba

    def __getitem__(self, index):
        def parseItem(index):
            if isinstance(index, numbers.Integral):
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
        
    def preparebinary(self, keys, values):
        raw_byte_count = 0
        pointer_i = 0
        bytestream = b''
        pointers = {}
        for key, value in zip(keys, values):
             # save bytestream info
            uncompressed_bytes = bytearray(value)
            compressed_bytes = self.compress(uncompressed_bytes, **self.compress_kwargs)
            bytestream += compressed_bytes
            # save pointers
            pointers[key] = (pointer_i, len(compressed_bytes))
            # tick forward counters
            raw_byte_count += len(uncompressed_bytes)
            pointer_i += len(compressed_bytes)
        return pointers, bytestream, raw_byte_count
    
    def subset(self, new_path, mask=None, index=None):
        if (mask is None and index is None) or (mask is not None and index is not None):
            raise ValueError('Define exactly one of mask and index.')
        # initialize new KBA
        new_kba = KeyedByteArray(
                new_path, mode='w', columns=self.col, dtype=self.dtype,
                compression=self.index['compression']['compression_type'],
                compression_kwargs=self.index['compression']['kwargs'])
        if mask is not None:  # write by mask positions
            if mask.dtype != bool or len(mask.shape) != 1:
                raise ValueError('Mask should be a 1D bool numpy array')
            for i in np.where(mask == True)[0]:
                new_kba.write(self.row[i], self.iloc[i])
        elif index is not None:  # write by explicit keys
            for key in index:
                new_kba.write(key, self.loc[key])
        new_kba.close()

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
        @subproc
        def rechunkJob(origin_path, tmp_path, col_idxs):
            # open origin kba
            origin_kba = KeyedByteArray(origin_path, mode='r')
            new_rows = origin_kba.col[  # new rows are origin columns
                col_idxs[0]:col_idxs[1]]
            # write to tmp numpy array in memory
            tmp_arr = np.empty(
                shape=(len(origin_kba.row), len(new_rows)),
                dtype=origin_kba.dtype)
            for i, row_key in enumerate(origin_kba.row):
                tmp_arr[i, :] = origin_kba.read(row_key)[col_idxs[0]:col_idxs[1]]
            # compress and output bytes
            pointers, compressed_bytes, raw_byte_count = origin_kba.preparebinary(
                new_rows, tmp_arr.T)
            origin_kba.close()
            return pointers, compressed_bytes, raw_byte_count
        
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
            pointers, compressed_bytes, raw_byte_count = client.gather(f)
            output_kba.writebinary(
                pointers, compressed_bytes, raw_bytes=raw_byte_count)
            f.release()
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
            with open(f'{filepath}.kbi', 'rb') as index_file:
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