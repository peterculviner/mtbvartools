#!/usr/bin/env -S python -u

import os, argparse, sys, random, subprocess, time
import mtbvartools as vt
import pandas as pd
from shutil import rmtree

# argument handling
parser = argparse.ArgumentParser(
    description='',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument(
    '-i', '--sra', type=str, required=True, 
    help='accession number from short read archive comma separated SRA accessions will be concatenated')
parser.add_argument(
    '-o', '--output', type=str, required=True, help='output file handle')
parser.add_argument(
    '-d', '--dir', type=str, default='.', help='output directory path')
parser.add_argument(
    '--tmp-path', type=str, default='/tmp/')
parser.add_argument(
    '--keep-tmp', action='store_true', help='keep step temporary files')
parser.add_argument(
    '--overwrite', action='store_true', help='ignore result files and overwrite')

args, _ = parser.parse_known_args()

# dealing with failure to connect to SRA
download_attempts = 3

tmp_dir = f'{args.tmp_path}/{args.output}/sra_download'
os.makedirs(tmp_dir, exist_ok=True)
base_dir = f'{args.dir}/{args.output}/sra_download'
os.makedirs(base_dir, exist_ok=True)

sra_dir = f'{tmp_dir}/sra/'
os.makedirs(sra_dir, exist_ok=True)

fq_dir = f'{tmp_dir}/fastq/'
os.makedirs(fq_dir, exist_ok=True)

if not args.overwrite and os.path.exists(f'{base_dir}/{args.output}.results.csv'):
    print(f'{base_dir}/{args.output}.results.csv found, exiting....')
    sys.exit(0)

# parse SRA input(s)
data_inputs = []
iidx = 0
try:
    for sra_input in args.sra.split(','):
        data_inputs.append([sra_input, f'_{iidx}'])
        iidx += 1
except TypeError:
    pass

# in the case of a single input, a file suffix is not needed
if len(data_inputs) == 1:
    # set suffix to empty string
    data_inputs[0][1] = ''

print('handling the following input(s)....')
for sra_input, suffix in data_inputs:
    print(f'{sra_input} > {args.output}{suffix}')


# step through processing and mapping of reads
fastq_outputs = []
fragments = []
is_paired = []
for sra_input, suffix in data_inputs:
    base_output = f'{args.output}{suffix}'
    for ntry in range(download_attempts):
        try:
            vt.contShell(
                f'rm -r -f {sra_dir}/{sra_input}')  # clear partial download(s) of SRA
            vt.contShell(
                f'prefetch -O {sra_dir} {sra_input}')
            break
        except subprocess.CalledProcessError as cpe:
            sleep_time = random.randint(0, 120)
            print(f'failed download on try {ntry + 1} retrying in {sleep_time} seconds.')
            time.sleep(sleep_time)
    # downloaded, now extract
    retval = vt.contShell(
        f'cd {sra_dir} && fasterq-dump -f -3 --skip-technical -O ../fastq {sra_input}', is_return=True)
    for val in retval.split('\n\n'):
        if val.startswith('spots read'):
            fragments.append(int(val.split(': ')[1].replace(',', '')))
            break
    else:
        raise ValueError(f'Could not parse counts for fasterq-dump returned value:\n{retval}')
    # check if paired end / single end / undetermined
    se = f'{fq_dir}/{sra_input}.fastq'
    pe1 = f'{fq_dir}/{sra_input}_1.fastq'
    pe2 = f'{fq_dir}/{sra_input}_2.fastq'
    if os.path.exists(pe1) and os.path.exists(pe2):
        target_fastq = [pe1, pe2]
        is_paired.append(True)
    elif os.path.exists(se):
        target_fastq = [se]
        is_paired.append(False)
    else:
        raise ValueError('SRA fastq does not fit expected PE or SE patterns.')
    fastq_outputs.append(target_fastq)

# verify is_paired integrity
if any(is_paired) != all(is_paired):
    raise NotImplementedError('Mixture of paired and single end SRAs not implemented.')


# combine outputs
if len(fastq_outputs) == 1:  # single fastq
    if len(fastq_outputs[0]) == 2:  # paired end
        vt.contShell(
            f'mv {fastq_outputs[0][0]} {base_dir}/{args.output}_1.fastq')
        vt.contShell(
            f'mv {fastq_outputs[0][1]} {base_dir}/{args.output}_2.fastq')
    elif len(fastq_outputs[0]) == 1:  # single end
        vt.contShell(
            f'mv {fastq_outputs[0][0]} {base_dir}/{args.output}.fastq')
    else:
        raise ValueError('Unexpected number of fastq outputs.')

elif len(fastq_outputs) > 1:  # concatenate SRAs
    if len(fastq_outputs[0]) == 2:  # paired end
        cat_string_1 = ' '.join([fq[0] for fq in fastq_outputs])
        vt.contShell(
            f'cat {cat_string_1} > {base_dir}/{args.output}_1.fastq')
        cat_string_2 = ' '.join([fq[1] for fq in fastq_outputs])
        vt.contShell(
            f'cat {cat_string_2} > {base_dir}/{args.output}_2.fastq')
    elif len(fastq_outputs[1]) == 1:  # single end
        cat_string = ' '.join([fq[0] for fq in fastq_outputs])
        vt.contShell(
            f'cat {cat_string} > {base_dir}/{args.output}.fastq')
    else:
        raise ValueError('Unexpected number of fastq outputs.')

# clean up temporary directories
vt.contShell(
    f'rm -r {fq_dir} {sra_dir}')

# write results file
pd.DataFrame(
    data=[[sum(fragments), all(is_paired)]],
    columns=['sra_fragments', 'is_paired'],
    index=[args.output]).to_csv(f'{base_dir}/{args.output}.results.csv')

print(f'Wrote {sum(fragments)} fragments, is_paired == {all(is_paired)}')

if args.keep_tmp is not True:
    rmtree(tmp_dir, ignore_errors=True)

sys.exit(0)