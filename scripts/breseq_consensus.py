#!/usr/bin/env -S python -u

import argparse, os, sys, zarr, pysam, re
import numpy as np
import pandas as pd
import mtbvartools as vt

def reformatBreseqVCF(input_vcf, output_vcf, label):
    with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
        for line in fin:
            outline = line  # default
            if line == '##fileDate\n':  # remove filedate line
                continue
            if line[:2] == '#C':
                outline = f'{line[:-1]}\tFORMAT\t{label}\n'
            if line[0] != '#':  # is a record
                outline = f'{line[:-1]}\tGT\t1/1\n'
                # check line integrity
                split_line = outline.split('\t')
                if (len(split_line[3]) < 1 or
                    len(split_line[4]) < 1 or
                    re.search('[^ATGC]+', split_line[3]) is not None or
                    re.search('[^ATGC]+', split_line[4]) is not None):
                    print(f'dropping poorly formed record in {label}: {outline}')
                    continue
            fout.write(outline)


def storeBreseqMissing(gd_path, zarray_path):
    zarr_miss = zarr.open(zarray_path, mode='a')
    # prepare tmp array
    tmp_array = np.zeros(zarr_miss.shape[0]).astype(bool)
    with open(gd_path) as f:
        for line in f:
            if line[:2] == 'UN':
                fields = line.split('\t')
                # both values are 1-indexed, end is inclusive
                start, end = int(fields[-2]), int(fields[-1])
                tmp_array[start - 1: end] = True
    # store in zarr array
    zarr_miss[:] = tmp_array


def bcftoolsStats(bcftools_keys, filename):
    output_data = {}
    for key, filterstr in bcftools_keys:
        if os.path.isfile(filename):
            output_data[key] = float(vt.contShell(
                f'bcftools view -H {filterstr} {filename} | wc -l', is_return=True).strip())
        else:
            output_data[key] = np.nan
    return pd.Series(output_data)


def getBreseqMiss(filename):
    try:
        return pd.Series({'breseq_miss_count': zarr.open(filename, 'r')[:].sum()})
    except:
        return pd.Series({'breseq_miss_count': np.nan})

# argument handling
parser = argparse.ArgumentParser(
    description='Runs Breseq in consensus (default) mode.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-i', '--in-fastq', type=str, required=True, 
    help='Path to fastq file stub.')
parser.add_argument(
    '-g', '--genbank', type=str, required=True, help='fasta file, needs to be pre-indexed')
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

base_dir = f'{args.dir}/{args.output}/breseq_consensus'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)


 # run the mapping
if is_paired:
    cmd = f'\
        breseq --brief-html-output \
        -n {args.output} \
        -o {base_dir} \
        -j {args.threads} \
        -r {args.genbank} \
        {fin_path1} {fin_path2} 2>&1'
else:
    cmd = f'\
        breseq --brief-html-output \
        -n {args.output} \
        -o {base_dir} \
        -j {args.threads} \
        -r {args.genbank} \
        {fin_path} 2>&1'
vt.contShell(cmd)

# reprocess breseq VCF
reformatBreseqVCF(
    f'{base_dir}/output/output.vcf',
    f'{base_dir}/{args.output}.breseq.vcf',
    args.output)

# store breseq miss locations
ref_len = pysam.FastaFile(f'{base_dir}/data/reference.fasta').lengths[0]
zarr_miss = zarr.open(
    mode='w', store=f'{base_dir}/{args.output}.miss.breseq.zarr', dtype='bool',
    shape=(ref_len, ), chunks=(-1, ),)
storeBreseqMissing(
    f'{base_dir}/output/output.gd',
    f'{base_dir}/{args.output}.miss.breseq.zarr')

# store breseq metrics and output results file
stats_keys = [
    ['breseq_consensus_snp_count', '-m2 -M2 -v snps'],]

results_df = pd.DataFrame(
    data=pd.concat([
        bcftoolsStats(stats_keys, f'{base_dir}/{args.output}.breseq.vcf'),
        getBreseqMiss(f'{base_dir}/{args.output}.miss.breseq.zarr')], axis=0),
    columns=[args.output]).T
results_df.to_csv(
    f'{base_dir}/{args.output}.results.csv')
sys.exit(0)