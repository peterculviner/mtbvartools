#!/usr/bin/env -S python -u

import subprocess, sys, argparse, time, os, dask
import pandas as pd
import numpy as np
from datetime import datetime
from dask_jobqueue import SLURMCluster

# run on command on shell while updating stdout
def contShell(cmd, return_output=False):
    def yieldcmd(cmd):
        popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)
    if return_output:
        return [line for line in yieldcmd(cmd)]
    else:
        for line in yieldcmd(cmd):
            print(line[:-1])


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
    '-p', '--n-processes', default=10, type=int, help='number of processes to start')
parser.add_argument(
    '--job-spec', default='shared,1,1,4GB,1:00:00', type=str, help='queue,processes,cores_per_process,memory_per_process,walltime')
parser.add_argument(
    '--dashboard-address', default=':10000', type=str)
parser.add_argument(
    '--wait', default=2, type=float, help='wait time (s) between submissions (default: 2s)')
args = parser.parse_args()

# prepare outfile directory
outfile_dir = datetime.now().strftime('outfiles_%y%m%d_%Hh%Mm%Ss')
os.makedirs(outfile_dir)
os.makedirs(f'{outfile_dir}/jobs')
os.makedirs(f'{outfile_dir}/workers')

print('launching client....')
queue, processes_per_node, cores_per_node, memory, walltime = args.job_spec.split(',')

cluster = SLURMCluster(
    processes=int(processes_per_node),
    cores=int(cores_per_node),
    memory=memory,
    queue=queue,
    walltime=walltime,
    worker_extra_args=[
        '--resources "jobs=1"'],
    job_extra_directives=[
        f'-o {outfile_dir}/workers/worker_node-%j.out'],
    scheduler_options={
        'dashboard_address': args.dashboard_address})

# start an initial maximum of 10 nodes, this is few enough to not upset the SLURM scheduler
current_scale = min(args.n_processes, 10 * int(processes_per_node))
cluster.scale(current_scale)  # scale refers to the number of processes to launch
print(f'requesting scale of {current_scale} processes....') 
client = cluster.get_client()


print('submitting work....')
futures = []

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
    futures.append(client.submit(contShell, cmd, resources={'jobs': 1}))


while current_scale < args.n_processes:
    time.sleep(args.wait)
    current_scale += int(processes_per_node)
    cluster.adapt(minimum=0, maximum=current_scale)
    print(f'  scaling to {current_scale}....')

# wait on completion of futures before closing
print('waiting for jobs to complete....')
client.gather(futures, errors='skip')
print('shutting down....')

# shutdown and exit
client.shutdown()
sys.exit(0)