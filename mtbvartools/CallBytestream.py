import numpy as np
import os
from tqdm import tqdm
from .KeyedByteArray import KeyedByteArray

class CallBytestream():
    def __init__(self, data_path, init_calls=True, init_nodes=True):
        self.path = data_path
        if init_calls:
            # load in data indexed by call
            self.calls = KeyedByteArray(
                f'{self.path}/by_variant.kba')
        if init_nodes:
            # load in data indexed by node
            self.nodes = KeyedByteArray(
                f'{self.path}/by_node.kba')

    def close(self):
        try: self.calls.close()
        except AttributeError: pass
        try: self.nodes.close()
        except AttributeError: pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def subsampleCallBytestream(input_path, output_path, call_names=None, node_names=None):
    input_vcb = CallBytestream(  # load in target vcb
        input_path)
    os.makedirs(output_path, exist_ok=True)
    # prepare indexes
    print('Finding call indexes....')
    if call_names == None:
        call_names = input_vcb.calls.row
        call_i = np.arange(len(input_vcb.calls.row))
    else:
        call_i = np.asarray(
            [input_vcb.nodes.col_dict[name] for name in call_names])
    print(f'Writing {len(call_i)}/{len(input_vcb.calls.row)}')
    print('Finding node indexes....')
    if node_names == None:
        node_names = input_vcb.nodes.row
        node_i = np.arange(len(input_vcb.nodes.row))
    else:
        node_i = np.asarray(
            [input_vcb.calls.col_dict[name] for name in call_names])
    print(f'Writing {len(node_i)}/{len(input_vcb.nodes.row)}')
    # writing by call
    print('Writing call KBA....')
    by_call_kba = KeyedByteArray(
        f'{output_path}/by_variant.kba', mode='w', columns=node_names, dtype='uint8',
        compression=input_vcb.calls.index['compression']['compression_type'],
        compression_kwargs=input_vcb.calls.index['compression']['kwargs'])
    for key, index in tqdm(list(zip(call_names, call_i))):
        by_call_kba.write(key, input_vcb.calls.iloc[index][node_i])
    by_call_kba.close()
    # writing by node
    print('Writing node KBA....')
    by_node_kba = KeyedByteArray(
        f'{output_path}/by_node.kba', mode='w', columns=call_names, dtype='uint8',
        compression=input_vcb.nodes.index['compression']['compression_type'],
        compression_kwargs=input_vcb.nodes.index['compression']['kwargs'])
    for key, index in tqdm(list(zip(node_names, node_i))):
        by_node_kba.write(key, input_vcb.nodes.iloc[index][call_i])
    by_node_kba.close()