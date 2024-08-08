#!/usr/bin/env -S python -u

import argparse, os, sys
import numpy as np
import pandas as pd
import mtbvartools as vt
from shutil import rmtree

# functions for parsing outputs
def parseBAMStats(filename):
    try:
        with open(filename, 'r') as file:
            for line in file:
                if line.startswith('SN\treads mapped:\t'):
                    mapped = int(line.replace('SN\treads mapped:\t', ''))
                if line.startswith('SN\tsequences:\t'):
                    total = int(line.replace('SN\tsequences:\t', ''))
                if line.startswith('SN\terror rate:\t'):
                    error_rate = float(line.replace('SN\terror rate:\t', '').split('\t')[0])
                if line.startswith('FFQ'):
                    break
        return pd.Series({
            'mapping_rate': mapped/total,
            'error_rate': error_rate})
    except:
        return pd.Series({
            'mapping_rate': np.nan,
            'error_rate': np.nan})
    
coverage_keys = ['coverage', 'meandepth']

def getBAMCoverage(filename, target_values):
    try:
        full_series = pd.read_csv(
            filename, delimiter='\t')
        return full_series.loc[0, target_values]
    except FileNotFoundError:
        return pd.Series({key: np.nan for key in target_values})

# argument handling
parser = argparse.ArgumentParser(
    description='Maps with BWA using default settings for PE and SE data.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-i', '--in-fastq', type=str, required=True, 
    help='Path to fastq file stub.')
parser.add_argument(
    '-f', '--fasta', type=str, required=True, help='fasta file, needs to be pre-indexed')
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path')
parser.add_argument(
    '-t', '--threads', type=int, default=1, help='number of threads to use')
parser.add_argument(
    '--tmp-path', type=str, default='/tmp/')
parser.add_argument(
    '--keep-tmp', action='store_true', help='keep step temporary files')
parser.add_argument(
    '--overwrite', action='store_true', help='ignore result files and overwrite')

args, _ = parser.parse_known_args()

# check if paired end / single end / undetermined
fin_path = f'{args.in_fastq}.fastq'
fin_path1 = f'{args.in_fastq}_1.fastq'
fin_path2 = f'{args.in_fastq}_2.fastq'
if os.path.exists(fin_path1) and os.path.exists(fin_path2):
    is_paired = True
elif os.path.exists(fin_path):
    is_paired = False
else:
    raise ValueError('No fastqs matching expected PE or SE patterns.\n' + '\n'.join([fin_path, 'or...', fin_path1, fin_path2]))

tmp_dir = f'{args.tmp_path}/bwa_map'
os.makedirs(tmp_dir, exist_ok=True)
base_dir = f'{args.dir}/{args.output}/bwa_map'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)


if is_paired:
    print(f'Running bwa on PE inputs....\n{fin_path1}\n{fin_path2}')
    cmd = f"\
        bwa mem -t {args.threads} \
        -R '@RG\\tID:{args.output}\\tSM:{args.output}' \
        -o {tmp_dir}/{args.output}.sam \
        {args.fasta} {fin_path1} {fin_path2} 2>&1"

else:  # single end
    print(f'Running bwa on SE input....\n{fin_path}')
    cmd = f"\
        bwa mem -t {args.threads} \
        -R '@RG\\tID:{args.output}\\tSM:{args.output}' \
        -o {tmp_dir}/{args.output}.sam \
        {args.fasta} {fin_path} 2>&1"
vt.contShell(cmd)

# produce a sorted, indexed bam
print('sorting and indexing bam file....')
cmd = f'\
    samtools sort -@ {args.threads} \
    -o {base_dir}/{args.output}.bam \
    {tmp_dir}/{args.output}.sam 2>&1 \
    && samtools index -@ {args.threads} {base_dir}/{args.output}.bam 2>&1'
vt.contShell(cmd)

# calculate statistics
cmd = f'\
    samtools coverage {base_dir}/{args.output}.bam > {base_dir}/{args.output}.bam.coverage'
vt.contShell(cmd)
cmd = f'\
    samtools stats {base_dir}/{args.output}.bam > {base_dir}/{args.output}.bam.stats'
vt.contShell(cmd)

# parse values, write results and exit
results_df = pd.DataFrame(
    data=pd.concat([
        parseBAMStats(f'{base_dir}/{args.output}.bam.stats'),
        getBAMCoverage(f'{base_dir}/{args.output}.bam.coverage', coverage_keys)], axis=0),
    columns=[args.output]).T
results_df.to_csv(
    f'{base_dir}/{args.output}.results.csv')

if args.keep_tmp is not True:
    rmtree(tmp_dir, ignore_errors=True)

sys.exit(0)