#!/usr/bin/env -S python -u

import argparse, pysam, os, sys
import numpy as np
import pandas as pd
import mtbvartools as vt



# argument handling
parser = argparse.ArgumentParser(
    description='',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '--in-stub', type=str, required=False, default='None',
    help='Path to fastq file stub.')
parser.add_argument(
    '--in-fastq', type=str, required=False, default='None',
    help='Explicit path to fastq file(s).')
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path')
parser.add_argument(
    '-r', '--reference-length', type=int, required=True)
parser.add_argument(
    '-t', '--target-depth', type=int, required=False, default=0,
    help='target depth of 0 (default) will copy all reads.')
parser.add_argument(
    '--seed', type=int, default=0, help='random seed, if 0 (default), will convert output string into seed for replicable results.')
parser.add_argument(
    '--overwrite', action='store_true', help='ignore result files and overwrite')

args, _ = parser.parse_known_args()

# define random seed
if args.seed == 0:
    np.random.seed(
        int.from_bytes(args.output.encode(), 'big') % 2**30)
else:
    np.random.seed(args.seed)

# handle file path or stub inputs
if args.in_stub != 'None':
    # check if paired end / single end / undetermined
    fin_path = f'{args.in_stub}.fastq'
    fin_path1 = f'{args.in_stub}_1.fastq'
    fin_path2 = f'{args.in_stub}_2.fastq'
    # check that input files exist
    if os.path.exists(fin_path1) and os.path.exists(fin_path2):
        is_paired = True
    elif os.path.exists(fin_path):
        is_paired = False
    else:
        raise ValueError(f'No fastqs found matching expected PE or SE patterns.\nInputs: --in-stub {args.in_stub} --in-fastq {args.in_fastq}')
elif args.in_fastq != 'None':
    if len(args.in_fastq.split(',')) == 1:
        fin_path = args.in_fastq
        is_paired = False
        if not os.path.exists(fin_path):
            raise ValueError(f'{fin_path} not found.')
    elif len(args.in_fastq.split(',')) == 2:
        fin_path1, fin_path2 = args.in_fastq.split(',')
        is_paired = True
        if not os.path.exists(fin_path1) or not os.path.exists(fin_path2):
            raise ValueError(f'{fin_path1} and/or {fin_path2} not found.')

base_dir = f'{args.dir}/{args.output}/downsample_fastq'
os.makedirs(base_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

if is_paired:
    print(f'Downsampling PE inputs....\n{fin_path1}\n{fin_path2}')
    with pysam.FastxFile(fin_path1) as fin1, pysam.FastxFile(fin_path2) as fin2:
        fin1_lens = [
            len(record.sequence) for record in fin1]
        fin2_lens = [
            len(record.sequence) for record in fin2]
    # calculate required fragments
    n_fragments = len(fin1_lens)
    required_fragments = int(args.target_depth * args.reference_length / np.mean(fin1_lens + fin2_lens) / 2)
    
    # write fragments
    fout_path1 = f'{base_dir}/{args.output}_1.fastq'
    fout_path2 = f'{base_dir}/{args.output}_2.fastq'
    # write all fragments (copy file)
    if required_fragments >= n_fragments or args.target_depth == 0:
        chosen_i = np.arange(n_fragments)
        print('copying all fragments')
    else:
        chosen_i = np.sort(
            np.random.choice(np.arange(n_fragments), size=required_fragments, replace=False))
        print('first 200 chosen indexes:\n', str(chosen_i[:200]))
    for fin_path, fout_path in zip([fin_path1, fin_path2], [fout_path1, fout_path2]):
        pointer = 0
        with pysam.FastxFile(fin_path) as fin, open(fout_path, 'w') as fout:
            for i, record in enumerate(fin):
                if chosen_i[pointer] == i:
                    fout.write(str(record) + '\n')
                    pointer += 1
                if pointer == len(chosen_i):
                    break
        n_written = pointer
        print(f'Wrote {n_written}/{n_fragments} reads to {fout_path}')
            
    
    # write results and exit
    pd.DataFrame(
        columns=['avglen_1', 'stdlen_1', 'avglen_2', 'stdlen_2', 'n_required', 'n_written'],
        data=[[np.mean(fin1_lens), np.std(fin1_lens), np.mean(fin2_lens), np.std(fin2_lens), required_fragments, n_written]],
        index=[args.output]).to_csv(
            f'{base_dir}/{args.output}.results.csv')
    sys.exit(0)


else:  # single end
    print(f'Downsampling SE input....\n{fin_path}')
    with pysam.FastxFile(fin_path) as fin:
        fin_lens = [
            len(record.sequence) for record in fin]
    # calculate required fragments
    n_fragments = len(fin_lens)
    required_fragments = int(args.target_depth * args.reference_length / np.mean(fin_lens))

    # write fragments
    fout_path = f'{base_dir}/{args.output}.fastq'
    # write all fragments (copy file)
    if required_fragments >= n_fragments or args.target_depth == 0:
        chosen_i = np.arange(n_fragments)
        print('copying all fragments')
    else:
        chosen_i = np.sort(
            np.random.choice(np.arange(n_fragments), size=required_fragments, replace=False))
        print('first 200 chosen indexes:\n', str(chosen_i[:200]))
    pointer = 0
    with pysam.FastxFile(fin_path) as fin, open(fout_path, 'w') as fout:
        for i, record in enumerate(fin):
            if chosen_i[pointer] == i:
                fout.write(str(record) + '\n')
                pointer += 1
            if pointer == len(chosen_i):
                break
    n_written = pointer
    print(f'Wrote {n_written}/{n_fragments} reads to {fout_path}')

    # write results and exit
    pd.DataFrame(
        columns=['avglen_1', 'stdlen_1', 'n_required',  'n_written'],
        data=[[np.mean(fin_lens), np.std(fin_lens), required_fragments, n_written]],
        index=[args.output]).to_csv(
            f'{base_dir}/{args.output}.results.csv')
    sys.exit(0)