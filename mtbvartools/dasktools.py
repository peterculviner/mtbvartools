import os, re
from datetime import datetime
from dask.distributed import get_client, LocalCluster
from dask_jobqueue import SLURMCluster
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


def timedsubproc(timeout):
    """
    Launches a python function as a separate process. Use as a decorator for jobs or functions that cause memory leaks.
    Requires provision of a time in seconds for timeout. Will kill subproc and issue TimeoutError at timeout.
    """
    def decorator(func):
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
                process.join(timeout=timeout)
                if process.is_alive():  # waited for process to complete, terminating instead
                    process.terminate()
                    process.join()
                    process.close()
                    raise TimeoutError(f'timeout expired at {timeout}')
                else:
                    output = queue.get()
                    process.close()
                    return output
        return wrapper
    return decorator


def findClient(n_workers=None, threads_per_worker=1, **kwargs):
    try:
        get_client()
    except ValueError:
        LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker, **kwargs).get_client()
    return get_client()


def startClient(
        n_workers=1,
        use_local=True,
        use_slurm=False,
        log_dir='.',
        queue='sapphire',
        process_per_node=1,
        cores_per_process=1,
        memory_per_process='4GB',
        walltime='1:00:00',
        dashboard_address=':10000',
        job_script_prologue=[],
        **kwargs):
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    if sum([use_local, use_slurm]) != 1:
        raise ValueError('Choose exactly 1 of {use_local, use_slurm}.')
    if use_local:
        print(f'{ts} - Launching LocalCluster client ({n_workers} workers)')
        return LocalCluster(
            n_workers=n_workers,
            threads_per_worker=1,
            dashboard_address=dashboard_address).get_client()
    if use_slurm:
        # prepare outfile directory
        outfile_dir = f'{log_dir}/' + datetime.now().strftime('outfiles_%y%m%d_%Hh%Mm%Ss')
        os.makedirs(outfile_dir)
        os.makedirs(f'{outfile_dir}/workers')
        cores_per_node = int(cores_per_process) * int(process_per_node)
        split_memory = re.split(r"(?<=\D)(?=\d)|(?<=\d)(?=\D)", memory_per_process)
        memory_per_node = str(int(split_memory[0]) * int(process_per_node)) + ''.join(split_memory[1:])
        print(f'{ts} - Launching SLURMCluster client ({n_workers} workers x {cores_per_process} CPU x {memory_per_process} x {walltime} @ {queue} with {process_per_node} workers / node)')
        cluster = SLURMCluster(
            processes=int(process_per_node),
            cores=cores_per_node,  # per node
            memory=memory_per_node,
            queue=queue,
            walltime=walltime,
            worker_extra_args=[
                '--resources "jobs=1"'],
            job_extra_directives=[
                f'-o {outfile_dir}/workers/worker_node-%j.out',
                '--open-mode=append'],
            scheduler_options={
                'dashboard_address': dashboard_address},
            job_script_prologue=job_script_prologue)
        cluster.scale(n_workers)
        return cluster.get_client()