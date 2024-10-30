#!/usr/bin/env -S python -u

import argparse, os, sys, pysam
import pandas as pd
import mtbvartools as vt
from mtbvartools.misc import StopWatch
from shutil import copy2, rmtree, copytree

def mergeResults(current_df, target_path):
    if os.path.exists(target_path):
        return pd.merge(
            current_df, pd.read_csv(target_path, index_col=0),
            left_index=True, right_index=True)
    else:
        print(f'Results file: "{target_path}" does not exist.')
        sys.exit(1)

def errorExit(args, results_df):
    if args.keep_tmp is not True:
        rmtree(tmp_path)
    results_df.loc[args.output, 'exit_type'] = 'ERROR'
    sw.end('pipeline', f'{results_dir}/times.txt', 'ERROR')
    sw.report(f'{results_dir}/times.csv')
    results_df.to_csv(f'{args.dir}/{args.output}/results/{args.output}.error.results.csv')
    vt.contShell(f'touch {args.dir}/{args.output}/results/{args.output}.error')


def politeExit(args, results_df, exit_type):
    if args.keep_tmp is not True:
        rmtree(tmp_path)
    results_df.loc[args.output, 'exit_type'] = exit_type
    sw.end('pipeline', f'{results_dir}/times.txt', exit_type)
    sw.report(f'{results_dir}/times.csv')
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
    '--sra', type=str, required=False, default='nan',
    help='SRA accession number.')
parser.add_argument(
    '--fastq-path', type=str, required=False, default='nan',
    help='FASTQ file(s), comma separated if paired (alternative to SRA accession number).')
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
    '--tmp-path', type=str, default='/tmp/', help='set to FALSE to keep tmp files in base directory.')
parser.add_argument(
    '--keep-tmp', action='store_true', help='keep step temporary files')
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
parser.add_argument(
    '--tbprofiler-fastq',
    action='store_true', help='run tbprofiler on fastq files instead of bwa output')


args = parser.parse_args()

# prepare results path
results_dir = f'{args.dir}/{args.output}/results/'
os.makedirs(results_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{results_dir}/{args.output}.results.csv'):
    print(f'FINAL SUMMARY FILE: {results_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

# prepare temporary path
if args.tmp_path == 'FALSE':
    tmp_path = f'{args.dir}/{args.output}/tmp/'
else:
    tmp_path = f'{args.tmp_path}/{args.output}/'
os.makedirs(tmp_path, exist_ok=True)

# start results DF
results_df = pd.DataFrame(index=[args.output])

# start global timer
sw = StopWatch(dict_path=f'{results_dir}/timings_dict.bin')
sw.start('pipeline', f'{results_dir}/times.txt')

try:  # global error handling
    #### PRE-SCRIPT CALCULATIONS ####
    # common args
    common_args = f'-o {args.output} -d {args.dir} --tmp-path {tmp_path}'
    if args.overwrite:
        common_args += ' --overwrite'
    if args.keep_tmp:
        common_args += ' --keep-tmp'
    # fasta reference length
    ref_len = pysam.FastaFile(
        args.fasta).lengths[0]
    ####

    
    #### SRA DOWNLOAD / FASTQ HANDLING ####
    if str(args.sra) != 'nan':
        sw.start('sra_download.py', f'{results_dir}/times.txt')
        cmd = f'\
            sra_download.py \
            --sra {args.sra} {common_args}'
        vt.contShell(cmd)
        sw.end('sra_download.py', f'{results_dir}/times.txt')
        results_df = mergeResults(  # get results
            results_df, f'{args.dir}/{args.output}/sra_download/{args.output}.results.csv')
    elif str(args.fastq_path) == 'nan':
        raise ValueError('Either --sra or --fastq-path must be defined.')
    ## move key files
    ## decide to proceed or stop
    ####


    #### DOWNSAMPLE FASTQ ####
    sw.start('downsample_fastq.py', f'{results_dir}/times.txt')
    if str(args.sra) != 'nan':
        cmd = f'\
            downsample_fastq.py \
            --in-stub {args.dir}/{args.output}/sra_download/{args.output} \
            --reference-length {ref_len} --target-depth {args.target_depth} \
            {common_args}'
    else:
        cmd = f'\
            downsample_fastq.py \
            --in-fastq {args.fastq_path} \
            --reference-length {ref_len} --target-depth {args.target_depth} \
            {common_args}'
    vt.contShell(cmd)
    sw.end('downsample_fastq.py', f'{results_dir}/times.txt')
    results_df = mergeResults(  # get results
        results_df, f'{args.dir}/{args.output}/downsample_fastq/{args.output}.results.csv')
    ## move key files
    ## decide to proceed or stop
    ####


    #### RUN FASTP ####
    sw.start('fastp_trimming.py', f'{results_dir}/times.txt')
    cmd = f'\
        fastp_trimming.py \
        --in-fastq {args.dir}/{args.output}/downsample_fastq/{args.output} \
        {common_args}'
    vt.contShell(cmd)
    sw.end('fastp_trimming.py', f'{results_dir}/times.txt')
    results_df = mergeResults(  # get results
        results_df, f'{args.dir}/{args.output}/fastp_trimming/{args.output}.results.csv')
    ## move key files
    [copy2(file, f'{args.dir}/{args.output}/results') for file in [
        f'{args.dir}/{args.output}/fastp_trimming/fastp.json',
        ]]
    ## decide to proceed or stop
    ####


    #### BWA MAP & STATS ####
    sw.start('bwa_map.py', f'{results_dir}/times.txt')
    cmd = f'\
        bwa_map.py \
        --in-fastq {args.dir}/{args.output}/fastp_trimming/{args.output} \
        --fasta {args.fasta} \
        --threads {args.threads} {common_args}'
    vt.contShell(cmd)
    sw.end('bwa_map.py', f'{results_dir}/times.txt')
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
    ####


    #### TBPROFILER ####
    sw.start('tb_profiler.py', f'{results_dir}/times.txt')
    if args.tbprofiler_fastq:
        cmd = f'\
            tb_profiler.py \
            --in-fastq {args.dir}/{args.output}/fastp_trimming/{args.output} \
            --lineage-snp-threshold {args.lineage_snp_threshold} \
            --lineage-snp-count {args.lineage_snp_count} \
            --threads {args.threads} {common_args}'
    else:
        cmd = f'\
            tb_profiler.py \
            --in-bam {args.dir}/{args.output}/bwa_map/{args.output}.bam \
            --lineage-snp-threshold {args.lineage_snp_threshold} \
            --lineage-snp-count {args.lineage_snp_count} \
            --threads {args.threads} {common_args}'
    vt.contShell(cmd)
    sw.end('tb_profiler.py', f'{results_dir}/times.txt')
    results_df = mergeResults(  # get results
        results_df, f'{args.dir}/{args.output}/tb_profiler/{args.output}.results.csv')
    ## move key files
    [copy2(file, f'{args.dir}/{args.output}/results') for file in [
        f'{args.dir}/{args.output}/tb_profiler/{args.output}.tbprofiler.json',
        ]]
    ## check stopping conditions
    if args.allowed_lineages != 'any' and results_df.loc[args.output, 'call'] not in args.allowed_lineages.split(','):
        print('Failed to find lineage call in allowed lineages.')
        politeExit(args, results_df, 'fail_lineage')
    ####


    #### MUTECT2 ####
    sw.start('mutect2.py', f'{results_dir}/times.txt')
    cmd = f'\
        mutect2.py \
        --in-bam {args.dir}/{args.output}/bwa_map/{args.output}.bam \
        -r {args.fasta} \
        --memory {args.memory} --threads {args.threads} {common_args}'
    vt.contShell(cmd)
    sw.end('mutect2.py', f'{results_dir}/times.txt')
    results_df = mergeResults(  # get results
        results_df, f'{args.dir}/{args.output}/mutect2/{args.output}.results.csv')
    ## move key files
    [copy2(file, f'{args.dir}/{args.output}/results') for file in [
        f'{args.dir}/{args.output}/mutect2/{args.output}.mutect2.pass.vcf',
        f'{args.dir}/{args.output}/mutect2/{args.output}.mutect2.fail.vcf',
        ]]
    ## check stopping conditions
    ####


    #### BRESEQ CONSENSUS ####
    sw.start('breseq_consensus.py', f'{results_dir}/times.txt')
    cmd = f'\
        breseq_consensus.py \
        --in-fastq {args.dir}/{args.output}/fastp_trimming/{args.output} \
        --genbank {args.genbank} \
        --threads {args.threads} {common_args}'
    vt.contShell(cmd)
    sw.end('breseq_consensus.py', f'{results_dir}/times.txt')
    results_df = mergeResults(  # get results
        results_df, f'{args.dir}/{args.output}/breseq_consensus/{args.output}.results.csv')
    ## move key files
    [copy2(file, f'{args.dir}/{args.output}/results') for file in [
        f'{args.dir}/{args.output}/breseq_consensus/{args.output}.breseq.vcf',
        ]]
    copytree(
        f'{args.dir}/{args.output}/breseq_consensus/{args.output}.miss.breseq.zarr',
        f'{args.dir}/{args.output}/results/{args.output}.miss.breseq.zarr',
        dirs_exist_ok=True)
    ## check stopping conditions
    ####

    politeExit(args, results_df, 'complete')

except Exception as original_exception:
    errorExit(args, results_df)
    raise original_exception