#!/usr/bin/env -S python -u

import argparse, os, sys, timeit
import pandas as pd
import mtbvartools as vt
import pysam
from shutil import copy2, rmtree

def mergeResults(current_df, target_path):
    if os.path.exists(target_path):
        return pd.merge(
            current_df, pd.read_csv(target_path, index_col=0),
            left_index=True, right_index=True)
    else:
        print(f'Results file: "{target_path}" does not exist.')
        sys.exit(1)


def politeExit(args, results_df, exit_type):
    results_df.loc[args.output, 'exit_type'] = exit_type
    if args.keep_intermediates is False:
        for name in os.listdir(f'{args.dir}/{args.output}'):
            if name != 'results':
                rmtree(f'{args.dir}/{args.output}/{name}')
    results_df.to_csv(f'{args.dir}/{args.output}/results/{args.output}.results.csv')
    sys.exit(0)
    

# argument handling
parser = argparse.ArgumentParser(
    description='Runs TBProfiler for BAM file, parses outputs.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# INPUTS
parser.add_argument(
    '--sra', type=str, required=True, 
    help='SRA accession number.')
parser.add_argument(
    '--genbank', type=str, required=True, 
    help='Complete path to genbank reference.')
parser.add_argument(
    '--fasta', type=str, required=True,
    help='Complete path to fasta reference. Before running, index (bwa index X.fasta) and dict (gatk CreateSequenceDictionary -R X.fasta).')

# COMMON INPUTS
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle.')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path.')
parser.add_argument(
    '-t', '--threads', type=int, default=1, help='number of threads to use.')
parser.add_argument(
    '-m', '--memory', type=str, default='6000m', help='amount of memory to allocate.')
parser.add_argument(
    '--overwrite', action='store_true', help='ignore result files and overwrite')
parser.add_argument(
    '--keep-intermediates', action='store_true', help='keep step intermediates')


## EXIT CONDITIONS
parser.add_argument(
    '--target-depth', type=int, default=0,
    help='Estimated depth to target for downsampling, 0 does no downsampling.')
parser.add_argument(
    '--min-depth', type=int, default=20,
    help='minimum mean depth required to proceed to variant calling')
parser.add_argument(
    '--min-aligned', type=float, default=0.8,
    help='minimum aligned fraction of reads by BWA MEM to proceed')
parser.add_argument(
    '--allowed-lineages', type=str, default='lineage1,lineage2,lineage3,lineage4,lineage5,lineage6,lineage7,lineage8,lineage9',
    help='comma separated list of allowed tbprofiler lineages. Set to string "any" to allow any lineage.')

## TB PROFILER INPUTS 
parser.add_argument(
    '--lineage-snp-threshold', type=float, default=90,
    help='minimum required threshold value (lineage snp mean in percent) for lineage snps to make a lineage call.')
parser.add_argument(
    '--lineage-snp-count', type=float, default=5,
    help='number of tbprofiler snps required to define a lineage')


args = parser.parse_args()

results_dir = f'{args.dir}/{args.output}/results/'
os.makedirs(results_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{results_dir}/{args.output}.results.csv'):
    print(f'FINAL SUMMARY FILE: {results_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

# start results DF
results_df = pd.DataFrame(index=[args.output])

# start global timer
global_start = timeit.default_timer()


#### PRE-SCRIPT CALCULATIONS ####
# common args
common_args = f'-o {args.output} -d {args.dir}'
if args.overwrite:
    common_args += ' --overwrite'
# fasta reference length
ref_len = pysam.FastaFile(
    args.fasta).lengths[0]
####


#### SRA DOWNLOAD ####
start = timeit.default_timer()
cmd = f'\
    pipeline_steps/sra_download.py \
    --sra {args.sra} {common_args}'
vt.contShell(cmd)
results_df = mergeResults(  # get results
    results_df, f'{args.dir}/{args.output}/sra_download/{args.output}.results.csv')
## move key files
## decide to proceed or stop
print(f'COMPLETE sra_download.py in {round(timeit.default_timer() - start)} seconds.\n\n')
####


#### DOWNSAMPLE FASTQ ####
start = timeit.default_timer()
cmd = f'\
    pipeline_steps/downsample_fastq.py \
    --in-fastq {args.dir}/{args.output}/sra_download/{args.output} \
    --reference-length {ref_len} --target-depth {args.target_depth} \
    {common_args}'
vt.contShell(cmd)
results_df = mergeResults(  # get results
    results_df, f'{args.dir}/{args.output}/downsample_fastq/{args.output}.results.csv')
## move key files
## decide to proceed or stop
print(f'COMPLETE downsample_fastq.py in {round(timeit.default_timer() - start)} seconds.\n\n')
####


#### BWA MAP & STATS ####
start = timeit.default_timer()
cmd = f'\
    pipeline_steps/bwa_map.py \
    --in-fastq {args.dir}/{args.output}/downsample_fastq/{args.output} \
    --fasta {args.fasta} \
    --threads {args.threads} {common_args}'
vt.contShell(cmd)
results_df = mergeResults(  # get results
    results_df, f'{args.dir}/{args.output}/bwa_map/{args.output}.results.csv')
## move key files
[copy2(file, f'{args.dir}/{args.output}/results') for file in [
    f'{args.dir}/{args.output}/bwa_map/{args.output}.bam.coverage',
    f'{args.dir}/{args.output}/bwa_map/{args.output}.bam.stats',
    ]]
## check stopping conditions
if results_df.loc[args.output, 'meandepth'] < args.min_depth:
    print('Failed minimum mean depth cutoff.')
    politeExit(args, results_df, 'fail_depth')
if results_df.loc[args.output, 'mapping_rate'] < args.min_aligned:
    print('Failed minimum mapping rate cutoff.')
    politeExit(args, results_df, 'fail_mapping_rate')
print(f'COMPLETE bwa_map.py in {round(timeit.default_timer() - start)} seconds.\n\n')
####


#### TBPROFILER ####
start = timeit.default_timer()
cmd = f'\
    pipeline_steps/tb_profiler.py \
    --in-bam {args.dir}/{args.output}/bwa_map/{args.output}.bam \
    --lineage-snp-threshold {args.lineage_snp_threshold} \
    --lineage-snp-count {args.lineage_snp_count} \
    --threads {args.threads} {common_args}'
vt.contShell(cmd)
results_df = mergeResults(  # get results
    results_df, f'{args.dir}/{args.output}/tb_profiler/{args.output}.results.csv')
## move key files
[copy2(file, f'{args.dir}/{args.output}/results') for file in [
    f'{args.dir}/{args.output}/tb_profiler/{args.output}.tbprofiler.json',
    ]]
## check stopping conditions
if results_df.loc[args.output, 'call'] not in args.allowed_lineages.split(','):
    print('Failed to find lineage call in allowed lineages.')
    politeExit(args, results_df, 'fail_lineage')
print(f'COMPLETE tb_profiler.py in {round(timeit.default_timer() - start)} seconds.\n\n')
####


#### MUTECT2 ####
start = timeit.default_timer()
cmd = f'\
    pipeline_steps/mutect2.py \
    --in-bam {args.dir}/{args.output}/bwa_map/{args.output}.bam \
    -r {args.fasta} \
    --memory {args.memory} --threads {args.threads} {common_args}'
vt.contShell(cmd)
results_df = mergeResults(  # get results
    results_df, f'{args.dir}/{args.output}/mutect2/{args.output}.results.csv')
## move key files
[copy2(file, f'{args.dir}/{args.output}/results') for file in [
    f'{args.dir}/{args.output}/mutect2/{args.output}.mutect2.pass.vcf',
    f'{args.dir}/{args.output}/mutect2/{args.output}.mutect2.fail.vcf',
    ]]
## check stopping conditions

print(f'COMPLETE mutect2.py in {round(timeit.default_timer() - start)} seconds.\n\n')
####


#### BRESEQ CONSENSUS ####
start = timeit.default_timer()
cmd = f'\
    pipeline_steps/breseq_consensus.py \
    --in-fastq {args.dir}/{args.output}/downsample_fastq/{args.output} \
    --genbank {args.genbank} \
    --threads {args.threads} {common_args}'
vt.contShell(cmd)
results_df = mergeResults(  # get results
    results_df, f'{args.dir}/{args.output}/breseq_consensus/{args.output}.results.csv')
## move key files
[copy2(file, f'{args.dir}/{args.output}/results') for file in [
    f'{args.dir}/{args.output}/breseq_consensus/{args.output}.breseq.vcf',
    ]]
## check stopping conditions
print(f'COMPLETE breseq_consensus.py in {round(timeit.default_timer() - start)} seconds.\n\n')
####


print(f'COMPLETE SCRIPT in {round(timeit.default_timer() - global_start)} seconds / {round((timeit.default_timer() - global_start)/60)} minutes.')
politeExit(args, results_df, 'complete')