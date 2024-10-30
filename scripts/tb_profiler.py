#!/usr/bin/env -S python -u

import argparse, os, sys
import numpy as np
import pandas as pd
import mtbvartools as vt
import json
from shutil import rmtree


def parseLineageCalls(path, lower=5, upper=95):
    strain_info = json.load(open(path, 'r'))
    # main_call = strain_info['main_lineage']
    # sub_call = strain_info['sub_lineage']
    evidence_dict = {}
    for lin_evidence in strain_info['lineage']:
        lin_dict = {}
        evidence_dict[lin_evidence['lineage']] = lin_dict
        lin_dict['low'], lin_dict['high'], lin_dict['mixed'] = 0, 0, 0
        support_fraction = []
        for variant in lin_evidence['support']:
            support_fraction.append(variant['target_allele_percent'])
            if variant['target_allele_percent'] < lower:
                lin_dict['low'] += 1
            elif variant['target_allele_percent'] > upper:
                lin_dict['high'] += 1
            else:
                lin_dict['mixed'] += 1
        lin_dict['mean'] = np.mean(support_fraction)
        lin_dict['median'] = np.median(support_fraction)
        lin_dict['std'] = np.std(support_fraction)
        lin_dict['raw'] = [int(np.round(s, 0)) for s in support_fraction]
    return evidence_dict


def getLineageCall(evidence_dict, target_level=0, min_snps=2, min_fraction=95):
    all_keys = np.asarray(
        [k for k in evidence_dict.keys()])
    key_levels = np.asarray(
        [k.count('.') for k in evidence_dict.keys()])
    target_keys = all_keys[key_levels == target_level]
    output_lineages = []
    for key in target_keys:
        if evidence_dict[key]['mean'] > min_fraction and len(evidence_dict[key]['raw']) > min_snps:
            output_lineages.append(key)
    return pd.Series(
        index=['call', 'n_calls'],
        data=[','.join(output_lineages), len(output_lineages)])


def getDeepestCall(evidence_dict, min_snps=2, min_fraction=95):
    evidence = pd.Series(evidence_dict)
    try:
        if len(evidence) == 0:
            return 'NONE'
        levels = evidence.index.str.count('\.')
        target_level = max(levels)
        while target_level != -1:
            target_evidence = evidence[levels == target_level]
            passing_evidence = []
            for i, e in target_evidence.items():
                if e['mean'] > min_fraction and len(e['raw']) > min_snps:
                    passing_evidence.append(i)
            if len(passing_evidence) > 1:
                return 'CONFLICT'
            if len(passing_evidence) == 1:
                return passing_evidence[0]
            target_level -= 1
        return 'NONE'
    except:
        return 'NONE'


def getLineageConflicts(evidence_dict):
        visited_keys = []  # record visited nodes
        all_keys = np.asarray(
            [k for k in evidence_dict.keys()])
        key_levels = np.asarray(
            [k.count('.') for k in evidence_dict.keys()])
        # step through levels
        visited_keys = []
        conflict_data = []
        for level in np.sort(np.unique(key_levels)):
            # check for conflicts on this level
            found_lineages = np.unique(
                ['.'.join(k.split('.')[:level + 1]) for k in all_keys if k.count('.') >= level])
            if len(found_lineages) < 2:
                visited_keys += list(all_keys[key_levels == level])
                continue  # no conflicts
            else:
                # get all sublineages of conflicting keys in a single list
                conflicting_keys = all_keys[np.any(
                    [[fndk in k for fndk in found_lineages] for k in all_keys], axis=1)]
                visited_keys += list(conflicting_keys)  # add to visited keys
                # record data on the conflict
                n_mixed, n_high = 0, 0
                proportion_dict = {
                    fndk: {'count': [], 'mean': []} for fndk in found_lineages}
                for k in conflicting_keys:
                    n_high += evidence_dict[k]['high']
                    n_mixed += evidence_dict[k]['mixed']
                    # populate the proportion dict
                    n_found = 0
                    for fndk in proportion_dict.keys():
                        if fndk in k:
                            proportion_dict[fndk]['count'].append(
                                evidence_dict[k]['high'] + evidence_dict[k]['mixed'] + evidence_dict[k]['low'])
                            proportion_dict[fndk]['mean'].append(evidence_dict[k]['mean'])
                            n_found += 1
                    if n_found > 1:
                        raise ValueError(f'multiple matches for {k} in {list(proportion_dict.keys())}!')
                # prepare output strings
                conflict_string = ', '.join(found_lineages)
                average_proportions = []
                n_snps = []
                for fndk in found_lineages:
                    average_proportions.append(
                        str(int(np.round(np.average(
                            proportion_dict[fndk]['mean'],
                            weights=proportion_dict[fndk]['count']), 0))))
                    n_snps.append(
                        str(np.sum(proportion_dict[fndk]['count'])))
                count_string = ', '.join(n_snps)
                proportion_string = ', '.join(average_proportions)
                conflict_data.append([
                conflict_string, count_string, proportion_string, n_high, n_mixed, (n_mixed / (n_mixed + n_high))])
                break
            # # stop iterating when all keys have been visited
            # if np.all(np.isin(visited_keys, all_keys)):
            #     break
        if len(conflict_data) == 0:
            conflict_data =[['', '', '', np.nan, np.nan, np.nan]]
        return pd.DataFrame(
            data=conflict_data,
            columns=['conflict_lineages', 'conflict_snps', 'conflict_avgs', 'n_high', 'n_mixed', 'mix_freq']).sort_values(
                'mix_freq', ascending=False)


def parseConflictData(filename, min_snps=5, min_fraction=95):
        try:
            edict = parseLineageCalls(filename)
            lineage_call = getLineageCall(edict, target_level=0, min_snps=min_snps, min_fraction=min_fraction)
            deepest_call = pd.Series(
                    index=['deepest_call'],
                    data=[getDeepestCall(edict, min_snps=min_snps, min_fraction=min_fraction)])
            conflict_data = getLineageConflicts(edict).iloc[0]
            return pd.concat([lineage_call, deepest_call, conflict_data])
        except FileNotFoundError:
            return pd.Series()


# argument handling
parser = argparse.ArgumentParser(
    description='Runs TBProfiler for BAM file, parses outputs.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '--in-bam', type=str, required=False, default='None',
    help='Complete path to bam file.')
parser.add_argument(
    '--in-fastq', type=str, required=False, default='None',
    help='Path to fastq stubs.')
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path')
parser.add_argument(
    '-t', '--threads', type=int, default=1, help='number of threads to use')
parser.add_argument(
    '--lineage-snp-threshold', type=float, default=95,
    help='minimum required threshold value (lineage snp mean in percent) for lineage snps to make a lineage call.')
parser.add_argument(
    '--lineage-snp-count', type=float, default=5,
    help='number of tbprofiler snps required to define a lineage')
parser.add_argument(
    '--tmp-path', type=str, default='/tmp/')
parser.add_argument(
    '--keep-tmp', action='store_true', help='keep step temporary files')
parser.add_argument(
    '--overwrite', action='store_true', help='ignore result files and overwrite')

args, _ = parser.parse_known_args()

tmp_dir = f'{args.tmp_path}/tb_profiler'
os.makedirs(tmp_dir, exist_ok=True)
base_dir = f'{args.dir}/{args.output}/tb_profiler'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)


# run TBProfiler
print('\nRunning TBProfiler for lineage calling....')
 # copy tb-profiler database to a temporary folder to prevent file access traffic jams
os.makedirs(f'{tmp_dir}/database', exist_ok=True)
cmd = f'\
    db=$(tb-profiler list_db | cut -f5); cp $db* {tmp_dir}/database'
vt.contShell(cmd)
if args.in_bam != 'None':
    print(f'running on BWA output {args.in_bam}')
    # check if files exist
    if not os.path.exists(args.in_bam):
        raise ValueError(f'BAM file {args.in_bam} not found.')
    # write a fresh bam file with "Chromosome" as reference name (required for proper execution of tb-profiler)
    cmd = f"\
        samtools view -h {args.in_bam} | sed 's/NC_000962\.3/Chromosome/g' | samtools view -b -o {tmp_dir}/{args.output}.bam"
    vt.contShell(cmd)
    # run tb-profiler, save results file
    cmd = f'\
        cd {tmp_dir} && tb-profiler profile --threads {args.threads} --external_db database/tbdb --bam {args.output}.bam --no_delly 2>&1'
    vt.contShell(cmd)
    cmd = f'\
        mv {tmp_dir}/results/tbprofiler.results.json {base_dir}/{args.output}.tbprofiler.json'
    vt.contShell(cmd)
elif args.in_fastq != 'None':
    print(f'running on fastq input {args.in_fastq}')
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
    # prepare tbprofiler run
    if is_paired:
        cmd = f'\
            cd {tmp_dir} && tb-profiler profile --threads {args.threads} --external_db database/tbdb --no_delly \
            -1 {os.path.abspath(fin_path1)} -2 {os.path.abspath(fin_path2)} 2>&1'
    else:
        cmd = f'\
            cd {tmp_dir} && tb-profiler profile --threads {args.threads} --external_db database/tbdb --no_delly \
            -1 {os.path.abspath(fin_path)} 2>&1'
    # run tbprofiler, save results file
    vt.contShell(cmd)
    cmd = f'\
        mv {tmp_dir}/results/tbprofiler.results.json {base_dir}/{args.output}.tbprofiler.json'
    vt.contShell(cmd)
else:
    raise ValueError('Provide one of --in-bam or --in-fastq')

# get lineage and conflict calls, write to results file
results_df = pd.DataFrame(
    data=parseConflictData(f'{base_dir}/{args.output}.tbprofiler.json', min_snps=args.lineage_snp_count, min_fraction=args.lineage_snp_threshold),
    columns=[args.output]).T
results_df.to_csv(
    f'{base_dir}/{args.output}.results.csv')

if args.keep_tmp is not True:
    rmtree(tmp_dir, ignore_errors=True)

sys.exit(0)