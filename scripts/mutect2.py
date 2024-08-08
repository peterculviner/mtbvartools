#!/usr/bin/env -S python -u

import argparse, os, sys, pysam
import mtbvartools as vt
import numpy as np
import pandas as pd
from shutil import rmtree


filter_functions = [
    lambda r: 'PASS' in r.filter.keys()]

def editMutect2VCF(in_path, out_path, fail_path=None, filter_functions=None):
    n_fail = 0
    n_pass = 0
    with pysam.VariantFile(in_path, 'r') as fin, pysam.VariantFile(out_path, 'w', header=fin.header) as passout:
        if fail_path is not None and filter_functions is not None:
            failout = pysam.VariantFile(fail_path, 'w', header=fin.header)
        sample = fin.header.samples[0]
        for record in fin.fetch():
            if len(record.alts) != 1:  # check for multiple ALTs
                raise ValueError(f'Expected records with exactly one ALT, found {len(record.alts)} at\n{record}')
            if set(record.samples[sample]['GT']) != set((0, 1)):
                raise ValueError(f'Unexpected GT value for Mutect2 in record\n{record}')
            # recalculate AF
            AD = record.samples[sample]['AD']
            try:
                AF = AD[1] / (AD[0] + AD[1])
            except ZeroDivisionError:
                AF = record.samples[sample]['AF']
                print(f'Zero AD found writing to fail VCF:\n{record}')
                failout.write(record)
                continue
            record.samples[sample]['AF'] = AF
            # reset GT based on haploid AD support
            record.samples[sample].phased = False  # remove phasing
            if AF == 1:  # no evidence for mixture
                record.samples[sample]['GT'] = (1, 1)
            elif AF > 0:  # possible mixed support
                record.samples[sample]['GT'] = (0, 1)
            else:  # no evidence for ALT, record fails
                failout.write(record)
                continue
            # write to fail or pass depending on filters
            if fail_path is not None and filter_functions is not None and not np.all([f(record) for f in filter_functions]):
                failout.write(record)
                n_fail += 1
            else:
                passout.write(record)
                n_pass += 1
        return pd.Series({'m2_pass': n_pass, 'm2_fail': n_fail})

def getMutect2SNPStats(target_vcf):
    snp_counts = {
        'm2_fix': 0,
        'm2_near_fix': 0,
        'm2_mix_call': 0,
        'm2_low_alt': 0,}
    unhandled = 0
    af_list = []
    with pysam.VariantFile(target_vcf, 'r') as fin:
        sample = fin.header.samples[0]
        for record in fin.fetch():
            if len(record.ref) == 1 and len(record.alts[0]) == 1:  # only consider SNPs
                AF = record.samples[sample]['AF'][0]
                af_list.append(AF)
                if AF == 1:
                    snp_counts['m2_fix'] += 1
                elif AF >= 0.9:
                    snp_counts['m2_near_fix'] += 1
                elif AF >= 0.1:
                    snp_counts['m2_mix_call'] += 1
                elif AF > 0:
                    snp_counts['m2_low_alt'] += 1
                else:
                    unhandled += 1
    return pd.Series(snp_counts), unhandled, af_list


# argument handling
parser = argparse.ArgumentParser(
    description='Runs TBProfiler for BAM file, parses outputs.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-i', '--in-bam', type=str, required=True, 
    help='Complete path to bam file; should be de-duplicated either by fastp or picard.')
parser.add_argument(
    '-r', '--reference', type=str, required=True,
    help='Complete path to fasta reference.')
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path')
parser.add_argument(
    '-t', '--threads', type=int, default=1, help='number of threads to use')
parser.add_argument(
    '-m', '--memory', type=str, default='6000m', help='amount of memory to allocate.')
parser.add_argument(
    '--tmp-path', type=str, default='/tmp/')
parser.add_argument(
    '--keep-tmp', action='store_true', help='keep step temporary files')
parser.add_argument(
    '--overwrite', action='store_true', help='ignore result files and overwrite')

args, _ = parser.parse_known_args()

# check if files exist
if not os.path.exists(args.in_bam):
    raise ValueError(f'BAM file {args.in_bam} not found.')

tmp_dir = f'{args.tmp_path}/mutect2'
os.makedirs(tmp_dir, exist_ok=True)
base_dir = f'{args.dir}/{args.output}/mutect2'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

print('\nRunning Mutect2....')

# --af-of-alleles-not-in-resource 0.05 to 0.0004 - using ~avg n variants / strain across L1-4 rel to H37Rv
cmd = f'export OMP_NUM_THREADS={args.threads} && \
    gatk --java-options "-Xmx{args.memory}" Mutect2 \
    -R {args.reference} \
    -I {args.in_bam} \
    -O {tmp_dir}/{args.output}.vcf  \
    --annotation StrandBiasBySample \
    --read-filter NotSupplementaryAlignmentReadFilter \
    --num-matching-bases-in-dangling-end-to-recover 1 \
    --min-base-quality-score 20 --dont-use-soft-clipped-bases \
    --bam-output {tmp_dir}/{args.output}.mt2.bam \
    --min-pruning 2  \
    --f1r2-tar-gz {tmp_dir}/{args.output}.f1r2.tar.gz \
    --af-of-alleles-not-in-resource 0.0004 \
    --max-reads-per-alignment-start 0 \
    2>&1'
vt.contShell(cmd)  # run GATK mutect2

cmd = f'gatk LearnReadOrientationModel \
    -I {tmp_dir}/{args.output}.f1r2.tar.gz \
    -O {tmp_dir}/{args.output}.readorientationmodel.tar.gz 2>&1'
vt.contShell(cmd)  # learn read orientations

# min reads per strand down from 5 -> 1, we can do additional downstream filtering
# remove hard filtering of AF (--min-allele-fraction 0.05), conduct own AF filtering
cmd = f'export OMP_NUM_THREADS={args.threads} && \
    gatk --java-options "-Xmx{args.memory}" FilterMutectCalls \
    -V {tmp_dir}/{args.output}.vcf \
    -R {args.reference} \
    -O {tmp_dir}/{args.output}.tmp.vcf \
	--microbial-mode \
	--min-reads-per-strand 1 \
	--max-events-in-region 1 \
    --ob-priors {tmp_dir}/{args.output}.readorientationmodel.tar.gz 2>&1'
vt.contShell(cmd)  # filter VCF using read orientation model

# split multiallele records
cmd = f'bcftools norm -m - {tmp_dir}/{args.output}.tmp.vcf -o {tmp_dir}/{args.output}.split.tmp.vcf 2>&1'
vt.contShell(cmd)

# recalculate AF and apply filters
pass_stats = editMutect2VCF(
    f'{tmp_dir}/{args.output}.split.tmp.vcf',
    f'{base_dir}/{args.output}.mutect2.pass.vcf',
    fail_path=f'{base_dir}/{args.output}.mutect2.fail.vcf',
    filter_functions=filter_functions)

# get SNP stats
snp_stats, _, _ = getMutect2SNPStats(
    f'{base_dir}/{args.output}.mutect2.pass.vcf')

# output results
results_df = pd.DataFrame(
    data=pd.concat([pass_stats, snp_stats]),
    columns=[args.output]).T

results_df.to_csv(
    f'{base_dir}/{args.output}.results.csv')

if args.keep_tmp is not True:
    rmtree(tmp_dir, ignore_errors=True)

sys.exit(0)