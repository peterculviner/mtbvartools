#!/usr/bin/env -S python -u

import sys, argparse, os, timeit, pysam, zarr, warnings, subprocess
import pandas as pd
import numpy as np
from dask.array import from_zarr
from dask.distributed import LocalCluster
import mtbvartools as vt
from mtbvartools.dasktools import subproc
from mtbvartools.vcf import bedToMask, maskToBed

warnings.simplefilter(action='ignore', category=FutureWarning)

if __name__ == '__main__':  # required for multiprocessing with dask "process" workers
    start_time = timeit.default_timer()  # timer

    # argument handling
    parser = argparse.ArgumentParser(
        description="""
        Merge VCF files described in sample sheet into a single fasta file for input into tree building software such as RAxML.
        Required columns:
          "output_name": output label for sample (from -o flag of scaled_tbprofiler_breseq.py input row)
          "output_dir": output directory for sample (from -d flag of scaled_tbprofiler_breseq.py input row)""",
        formatter_class=argparse.RawTextHelpFormatter)

    # required/global inputs
    parser.add_argument(
        '--input-csv', type=str, required=True, help='csv sample sheet, must have "label", "vcf_path", and "miss_path" columns.')
    parser.add_argument(
        '--input-fasta', type=str, required=True, help='input fasta')
    parser.add_argument(
        '--inputs-dir', type=str, default='.', help='base directory to look for files in.')
    parser.add_argument(
        '--out-dir', type=str, default='.', help='output directory')
    parser.add_argument(
        '--output', type=str, default='ouptut', help='output file name')
    parser.add_argument(
        '--mask', type=str, help='bed file (tsv) with regions to remove from all genomes')
    parser.add_argument(
        '--miss-threshold', type=float, default=1, help='')
    
    # multithreading options
    parser.add_argument(
            '--local-threads', type=int, default=1)

    args = parser.parse_args()

    # start workers
    # spawn dask processes
    print('Spawning workers....')
    cluster = LocalCluster(
        n_workers=args.local_threads,
        threads_per_worker=1)
    client = cluster.get_client()

    # parse sample sheet for errors
    print('Reading sample sheet....')
    sample_sheet_df = pd.read_csv(args.input_csv)
    if any(~np.isin(['label', 'vcf_path', 'miss_path'], sample_sheet_df.columns)):
        raise ValueError('One or more required columns in sample sheet is missing.')
    do_not_exist = []
    # check for vcf
    for i, rdata in sample_sheet_df.iterrows():
        for f in [
                f'{args.inputs_dir}/{rdata.vcf_path}',
                f'{args.inputs_dir}/{rdata.miss_path}']:
            if not os.path.exists(f):
                do_not_exist.append(f)
    if len(do_not_exist) > 0:
        output_str = '\n  '.join(do_not_exist)
        raise ValueError(f'\nThe following files were not found:\n  {output_str}')

    print('Initializing outputs....')
    # generate directory structure
    tmp_dir = f'{args.out_dir}/tmp/'
    os.makedirs(f'{tmp_dir}', exist_ok=True)

    # initialize zarr array tracking miss locations
    fasta_file = pysam.FastaFile(args.input_fasta)
    genome_len = fasta_file.lengths[0]
    miss_array = zarr.open(
            f'{args.out_dir}/{args.output}.miss.zarr', mode='w',
            shape=(len(sample_sheet_df), genome_len),
            chunks=(1, genome_len),
            dtype=bool)

    # initialize zarr array tracking all var locations
    var_array = zarr.open(
            f'{args.out_dir}/{args.output}.variants.zarr', mode='w',
            shape=(len(sample_sheet_df), genome_len),
            chunks=(1, genome_len),
            dtype=bool)
    
    # initialize results dictionary
    results_dict = {'n_samples': len(sample_sheet_df)}

    print('Generating miss mask....')
    @subproc
    def writeMiss(arr_idx, miss_path):
        miss_indexes = zarr.open(miss_path, mode='r')[:]  # load miss zarr
        # store missing sites in large zarr
        miss_array[arr_idx, :] = miss_indexes

    futures = []
    for arr_idx, rdata in sample_sheet_df.iterrows():
        futures.append(client.submit(
            writeMiss,
            arr_idx,
            f'{args.inputs_dir}/{rdata.miss_path}',
            priority=len(sample_sheet_df)-arr_idx))
    # sequentially print futures as completed
    for i, f in enumerate(futures):
        print(f'{i + 1} / {len(futures)}')
        client.gather(f)
        f.release()
        
    # MISSING SAMPLE HANDLING
    # sum the number of missing positions across all samples
    miss_sum = from_zarr(f'{args.out_dir}/{args.output}.miss.zarr').sum(axis=0)
    fraction_miss = miss_sum.compute() / len(sample_sheet_df)
    miss_mask = fraction_miss >= args.miss_threshold
    results_dict['n_miss_masked'] = np.sum(miss_mask)

    print('Combine masks....')
    # load pre-defined BED mask (if any)
    if args.mask is not None:
        bed_mask = bedToMask(args.mask, genome_len)
    else:
        bed_mask = np.zeros(genome_len).astype(bool)
    # combine with miss mask
    combined_mask = np.any(
        [miss_mask, bed_mask], axis=0)
    maskToBed(combined_mask, f'{args.out_dir}/combined_mask.bed', 'chromosome')
    # write results
    results_dict['n_in_bed_masked'] = np.sum(bed_mask)
    results_dict['total_masked'] = np.sum(combined_mask)
    
    # prepare function to parallelize
    @subproc
    def preparefasta(arr_idx, label, vcf_path, miss_path):
        output = f'{label}:\n'
        
        # write a SNP-only VCF
        output += vt.contShell(f'\
            bcftools view {vcf_path} -m2 -M2 -v snps | \
            awk \'{{gsub(/NC_000962/,"NC_000962.3"); print}}\' - | \
            bcftools view -o {tmp_dir}/{label}.tmp.vcf.gz -O z \
            && tabix {tmp_dir}/{label}.tmp.vcf.gz', is_return=True)
        
        # edit the input fasta with the variants
        result = subprocess.run(f'\
            bcftools consensus -f {args.input_fasta} -H 1 {tmp_dir}/{label}.tmp.vcf.gz',
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        output += result.stderr
        genome_str = ''.join(result.stdout.split('\n')[1:])
        
        # apply individual missing sites to the fasta
        miss_indexes = miss_array[arr_idx, :]
        fasta_array = np.asarray(
            [bp for bp in genome_str])
        fasta_array[miss_indexes] = 'N'

        # apply combined masked sites to the fasta - these sites will be removed completely
        cmask = bedToMask(f'{args.out_dir}/combined_mask.bed', genome_len)
        fasta_array[cmask] = ''

        # write an output fasta file
        with open(f'{tmp_dir}/{label}.fasta', 'w') as f:
            f.write(f'>{label}\n{"".join(fasta_array)}\n')
        
        # store the variant positions - ignore possible variants in a sample's miss regions
        variant_file = pysam.VariantFile(
            f'{tmp_dir}/{label}.tmp.vcf.gz')
        is_variant = np.zeros(genome_len).astype(bool)
        var_pos = [v.pos - 1 for v in variant_file.fetch()]  # vcf files are read in as 1-idx
        is_variant[var_pos] = True
        # ignore miss locations
        is_variant[miss_indexes] = False
        # ignore variants in masked regions
        combined_mask = bedToMask(
            f'{args.out_dir}/combined_mask.bed', genome_len)
        is_variant[combined_mask] = False
        # write to variant array
        var_array[arr_idx, :] = is_variant
        return output

    print('Preparing output fastas....')
    # write separate fasta files in parallel using above function
    # (1) write a SNP only VCF
    # (2) edit the genome using bcftools consensus
    # (3) apply combined BED mask excluding particular regions
    # (4) store which positions were variant in a zarr file for summing
    # (5) write a fasta file for this sample to the temporary folder
    # (6) store variant position information
    futures = []
    for arr_idx, rdata in sample_sheet_df.iterrows():
        futures.append(client.submit(
            preparefasta,
            arr_idx,
            rdata.label,
            f'{args.inputs_dir}/{rdata.vcf_path}',
            f'{args.inputs_dir}/{rdata.miss_path}',
            priority=len(sample_sheet_df)-arr_idx))
    # sequentially print futures as completed
    for i, f in enumerate(futures):
        print(f'{i + 1} / {len(futures)}')
        print(client.gather(f))
        f.release()

    # WRITE SUMMARY INFO
    # run var sum in parallel
    var_sum = from_zarr(f'{args.out_dir}/{args.output}.variants.zarr').sum(axis=0)
    n_alt = var_sum.compute()  # calculate the number of strains with alt call at this position
    # write summary to dictionary
    results_dict['genome_len'] = genome_len
    results_dict['miss_threshold'] = args.miss_threshold
    results_dict['variant_sites'] = np.sum(n_alt != 0)
    
    print('Writing results summaries....')
    pd.Series(
        data=results_dict,
        name=args.output).to_csv(f'{args.out_dir}/{args.output}.results.csv')
    
    value, count = np.unique(
        n_alt, return_counts=True)
    pd.Series(
        index=value, data=count, name='pos_counts').to_csv(
        f'{args.out_dir}/{args.output}.histogram.csv')

    print('Removing tmp files....')
    vt.contShell(f'\
        for file in {tmp_dir}/*.fasta; do cat "$file" >> {args.out_dir}/{args.output}.fasta; done')
    vt.contShell(f'\
        rm -r {tmp_dir}')
    
    elapsed = timeit.default_timer() - start_time
    print(f'finished in {int(elapsed)}s or {int(elapsed/60)}m')
    client.shutdown()
    sys.exit(0)