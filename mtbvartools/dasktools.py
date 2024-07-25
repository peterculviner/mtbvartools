from dask.distributed import get_client, LocalCluster
import multiprocess as mp
from functools import wraps

def subproc(func):
    """
    Launches a python function as a separate process. Use as a decorator for jobs or functions that cause memory leaks.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        in_subproc = kwargs.pop('_subproc', False)
        queue = kwargs.pop('_queue', None)
        if in_subproc:  # run function, push results to queue
            output = func(*args, **kwargs)
            queue.put(output)
        else:  # spawn subprocess
            queue = mp.Queue()
            process = mp.Process(  # pass queue and this function
                target=wrapper, args=args, kwargs={'_subproc': True, '_queue': queue, **kwargs})
            process.start()
            output = queue.get()
            process.join()
            process.close()
            return output
    return wrapper

def findClient(n_workers=None, threads_per_worker=1, **kwargs):
    try:
        get_client()
    except ValueError:
        LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker, **kwargs).get_client()
    return get_client()


