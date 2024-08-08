#!/usr/bin/env -S python -u

import argparse, os, sys
import pandas as pd
import mtbvartools as vt
from mtbvartools.misc import getJSONValues

# values to pull from fastp json output
fastp_lookup = [
    ['fastp_type', ('summary', 'sequencing')],
    ['reads_pre_filter', ('summary', 'before_filtering', 'total_reads')],
    ['reads_post_filter', ('summary', 'after_filtering', 'total_reads')],
    ['read1_len_post_filter', ('summary', 'after_filtering', 'read1_mean_length')],
    ['read2_len_post_filter', ('summary', 'after_filtering', 'read2_mean_length')],
    ['adapter_trimmed_reads', ('adapter_cutting', 'adapter_trimmed_reads')],
    ['adapter_trimmed_bases', ('adapter_cutting', 'adapter_trimmed_bases')],
    ['duplicate_rate', ('duplication', 'rate')],]


# argument handling
parser = argparse.ArgumentParser(
    description='Removes duplicates for PE and SE data and trims adapters (if found).',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-i', '--in-fastq', type=str, required=True, 
    help='Path to fastq file stub.')
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path')
parser.add_argument(
    '-t', '--threads', type=int, default=1, help='number of threads to use')
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

base_dir = f'{args.dir}/{args.output}/fastp_trimming'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

if is_paired:
    print(f'Running fastp on PE inputs....\n{fin_path1}\n{fin_path2}')
    cmd = f'\
        fastp -D -w {args.threads} \
        -i {fin_path1} \
        -I {fin_path2} \
        -o {base_dir}/{args.output}_1.fastq \
        -O {base_dir}/{args.output}_2.fastq \
        -j {base_dir}/fastp.json -h {base_dir}/fastp.html 2>&1'

else:  # single end
    print(f'Running fastp on SE input....\n{fin_path}')
    cmd = f'\
        fastp -D -w {args.threads} \
        -i {fin_path} \
        -o {base_dir}/{args.output}.fastq \
        -j {base_dir}/fastp.json -h {base_dir}/fastp.html 2>&1'
vt.contShell(cmd)
    
# parse values, write results and exit
results_df = pd.DataFrame(getJSONValues(f'{base_dir}/fastp.json', fastp_lookup)).T
results_df.index = [args.output]
results_df.to_csv(
    f'{base_dir}/{args.output}.results.csv')
sys.exit(0)