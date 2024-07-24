#!/usr/bin/env -S python -u

import argparse, os, sys
import numpy as np
import pandas as pd
import mtbvartools as vt
import json


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
            conflict_data = getLineageConflicts(edict).iloc[0]
            return pd.concat([lineage_call, conflict_data])
        except FileNotFoundError:
            return pd.Series()


# argument handling
parser = argparse.ArgumentParser(
    description='Runs TBProfiler for BAM file, parses outputs.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-i', '--in-bam', type=str, required=True, 
    help='Complete path to bam file.')
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
    '--overwrite', action='store_true', help='ignore result files and overwrite')

args = parser.parse_args()

# check if files exist
if not os.path.exists(args.in_bam):
    raise ValueError(f'BAM file {args.in_bam} not found.')

base_dir = f'{args.dir}/{args.output}/tbprofiler/'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

# run TBProfiler
print('\nRunning TBProfiler for lineage calling....')
# write a fresh bam file with "Chromosome" as reference name (required for proper execution of tb-profiler)
cmd = f"\
    samtools view -h {args.in_bam} | sed 's/NC_000962\.3/Chromosome/g' | samtools view -b -o {base_dir}/{args.output}.bam"
vt.contShell(cmd)
# copy tb-profiler database to a temporary folder to prevent file access traffic jams
os.makedirs(f'{base_dir}/database', exist_ok=True)
cmd = f'\
    db=$(tb-profiler list_db | cut -f5); cp $db* {base_dir}/database'
vt.contShell(cmd)
# run tb-profiler, save results file
cmd = f'\
    cd {base_dir} && tb-profiler profile --threads {args.threads} --external_db database/tbdb --bam {args.output}.bam --no_delly 2>&1'
vt.contShell(cmd)
cmd = f'\
    mv {base_dir}/results/tbprofiler.results.json {base_dir}/{args.output}.tbprofiler.json'
vt.contShell(cmd)

# get lineage and conflict calls, write to results file
evidence_dict = parseLineageCalls(f'{base_dir}/{args.output}.tbprofiler.json')
lineage_calls = getLineageCall(
    evidence_dict, target_level=0, min_snps=args.lineage_snp_count, min_fraction=args.lineage_snp_threshold)
results_df = pd.DataFrame(
    data=parseConflictData(f'{base_dir}/{args.output}.tbprofiler.json', min_snps=args.lineage_snp_count, min_fraction=args.lineage_snp_threshold),
    columns=[args.output]).T
results_df.to_csv(
    f'{base_dir}/{args.output}.results.csv')
sys.exit(0)