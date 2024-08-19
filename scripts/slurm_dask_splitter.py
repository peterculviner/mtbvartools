#!/usr/bin/env -S python -u

import sys, argparse, time, os, re
import pandas as pd
import numpy as np
from datetime import datetime
import mtbvartools as vt
from dask_jobqueue import SLURMCluster
from dask.distributed import fire_and_forget
from dask import config

config.set(  # jobs expensive and text is cheap to send over network
    {'distributed.scheduler.worker-saturation': 1,
     'distributed.scheduler.unknown-task-duration': '1000s'})


def getIdle(client):
    return set([
        w for w, ws in client.cluster.scheduler.workers.items() if len(ws.processing) == 0])

def incrementalUpscale(client, start_scale, workers_per_tick, maximum_workers, wait=1):
    cluster = client.cluster
    current_scale = start_scale
    while current_scale <= maximum_workers:
        n_tasks = len(cluster.scheduler.tasks)
        n_workers = len(cluster.scheduler.workers)
        if n_tasks <= n_workers:
            print(f'Number of workers, {n_workers}, greater than number of tasks, {n_tasks} - stop upscaling.')
            break
        time.sleep(wait)
        cluster.scale(current_scale)
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f'{ts} - Working on {n_tasks} tasks, scaling {n_workers} (current) -> {current_scale} (set).')
        current_scale += workers_per_tick

def incrementalDownscale(maximum_workers, client, wait=60):
    idle_last_tick = set()
    n_workers = len(client.cluster.scheduler.workers)
    n_tasks = len(client.cluster.scheduler.tasks)
    cluster = client.cluster
    while n_tasks > 0:
        ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if n_tasks >= n_workers:  # work aplenty
            target_workers = min(maximum_workers, n_tasks)
            if target_workers > n_workers:
                print(f'{ts} - {n_tasks} tasks / {n_workers} active workers. Rescaling to {target_workers}')
                cluster.scale(target_workers)
            else:
                print(f'{ts} - {n_tasks} tasks / {n_workers} active workers.')
        else:  # start retiring idle workers
            idle_this_tick = getIdle(client)
            to_close = idle_last_tick.intersection(idle_this_tick)
            client.retire_workers(to_close)
            idle_last_tick = idle_this_tick
            print(f'{ts} - {n_tasks} tasks / {n_workers} active workers. Closing {len(to_close)} idle workers.')
        # update stats
        time.sleep(wait)
        n_workers = len(client.cluster.scheduler.workers)
        n_tasks = len(client.cluster.scheduler.tasks)
    print('No tasks left! Closing the cluster.')
    client.shutdown()


# handle arguments
parser = argparse.ArgumentParser(
    description="""Run a process in individual jobs on SLURM using a spreadsheet for arguments. argsheet specifications:
    required columns:
       "label": job label to be integrated into outfile
    format specifications:
    - positional arguments columns should be called "pos" and will be handled in order given
    - non-positional arguments columns should be called "opt <flag> <value>"
    - to NOT call an argument in a row, set the value to NA (see pandas NA values in function pd.read_csv)
    - for non-positional arguments without a value, columns should be called "flg <flag>" if the value is True, it will be called""",
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    'process', type=str, help='path to process to run')
parser.add_argument(
    'argsheet', type=str, 
    help='path to csv containing arguments')
parser.add_argument(
    '--log-dir', default='.', type=str, help='location for outfile logging')
parser.add_argument(
    '--max-processes', default=10, type=int, help='maximum number of processes to start')
parser.add_argument(
    '--queue', default='sapphire', type=str, help='slurm queue to submit to')
parser.add_argument(
    '--process-per-node', default=1, type=int, help='n processes per node')
parser.add_argument(
    '--cores-per-process', default=1, type=int, help='n cores per process')
parser.add_argument(
    '--memory-per-process', default='4GB', type=str, help='memory per process')
parser.add_argument(
    '--walltime', default='1:00:00', type=str, help='memory per process')
parser.add_argument(
    '--dashboard-address', default=':10000', type=str)
parser.add_argument(
    '--upscale-wait', default=1, type=float, help='wait time (s) between scale ticks (default: 1s)')
parser.add_argument(
    '--downscale-wait', default=30, type=float, help='wait time (s) between downscale ticks (default: 30s)')
parser.add_argument(
    '--initial-nodes', default=50, type=int, help='maximum number of nodes to request to start.')
args = parser.parse_args()

# prepare outfile directory
outfile_dir = f'{args.log_dir}/' + datetime.now().strftime('outfiles_%y%m%d_%Hh%Mm%Ss')
os.makedirs(outfile_dir)
os.makedirs(f'{outfile_dir}/jobs')
os.makedirs(f'{outfile_dir}/workers')

ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
print(f'{ts} - Launching client....')

cores_per_node = int(args.cores_per_process) * int(args.process_per_node)
split_memory = re.split(r"(?<=\D)(?=\d)|(?<=\d)(?=\D)", args.memory_per_process)
memory_per_node = str(int(split_memory[0]) * int(args.process_per_node)) + ''.join(split_memory[1:])

cluster = SLURMCluster(
    processes=int(args.process_per_node),
    cores=cores_per_node,  # per node
    memory=memory_per_node,
    queue=args.queue,
    walltime=args.walltime,
    worker_extra_args=[
        '--resources "jobs=1"'],
    job_extra_directives=[
        f'-o {outfile_dir}/workers/worker_node-%j.out',
        '--open-mode=append'],
    scheduler_options={
        'dashboard_address': args.dashboard_address})

# use initial nodes setting to set initial node number
starting_scale = min(args.max_processes, int(args.initial_nodes) * int(args.process_per_node))
cluster.scale(starting_scale)  # scale refers to the number of processes to launch
print(f'Starting with scale of {starting_scale} processes....') 
client = cluster.get_client()

print('Submitting tasks....')
process_df = pd.read_csv(args.argsheet)
for i, (label, row) in enumerate(process_df.iterrows()):
    label = ''
    sbatch = ''
    pos_args = []
    flag_args = []
    # parse argument csv
    for flag, val in zip(process_df.columns, row):
        if not isinstance(val, str) and np.isnan(val):
            continue  # skip argument for this call
        if flag == 'label':
            label = str(val)
        elif flag[:4] == 'pos ':
            pos_args.append(str(val))
        elif flag[:4] == 'opt ':
            flag_args.append(flag[4:])
            flag_args.append(str(val))
        elif flag[:4] == 'flg ':
            if val:
                flag_args.append(flag[4:])
        else:
            raise ValueError(f'Unhandled argument column type {flag} / {val}')
    # format command and submit
    pos_str = ' '.join(pos_args)
    flag_str = ' '.join(flag_args)
    cmd = f'{args.process} {pos_str} {flag_str} >> {outfile_dir}/jobs/{label}.out'
    fire_and_forget(
        client.submit(vt.contShell, cmd, resources={'jobs': 1}))

# let the scheduler get its bearings on tasks    
time.sleep(60)

print('\n\nStarting incremental upscaling....')
incrementalUpscale(
    client, starting_scale, args.process_per_node, args.max_processes,
    wait=args.upscale_wait)

# let the scheduler get its bearings on tasks    
time.sleep(60)

print('\n\nStarting incremental downscaling....')
incrementalDownscale(
    args.max_processes, client, wait=args.downscale_wait)

sys.exit(0)