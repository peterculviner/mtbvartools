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