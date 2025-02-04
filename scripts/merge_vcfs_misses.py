#!/usr/bin/env -S python -u

import argparse, os, pysam, zarr, sys, timeit, glob
import pandas as pd
import numpy as np
import os.path
from rechunker import rechunk
from dask.distributed import as_completed, wait
import mtbvartools as vt
from mtbvartools.dasktools import subproc, startClient

env_bin = ''  # for testing if environment path needed

if __name__ == '__main__':  # required for multiprocessing with dask "process" workers

    start_time = timeit.default_timer()  # timer

    @subproc
    def mergeJob(vcf_list, output_name, chromosome, start, end, merge_level, threads=1):
        # prepare merge list
        merge_string = ' '.join(vcf_list)
        if merge_level == 1:
            if len(vcf_list) > 1:
                # conduct the merge
                vt.contShell(f'\
                    {env_bin}bcftools merge --threads {threads} \
                    -r {chromosome}:{start}-{end} \
                    -i "-" --missing-to-ref -m none \
                    {merge_string} | \
                    {env_bin}bcftools view --threads {threads} \
                    -O z0 -o {output_name} -t {chromosome}:{start}-{end} \
                    && {env_bin}tabix {output_name}')
            else:
                vt.contShell(f'\
                    {env_bin}bcftools view --threads {threads} \
                    -O z0 -o {output_name} \
                    -t {chromosome}:{start}-{end} -r {chromosome}:{start}-{end} \
                    {merge_string} && {env_bin}tabix {output_name} ')
        else:
            if len(vcf_list) > 1:
                # conduct the merge
                vt.contShell(f'\
                    {env_bin}bcftools merge --threads {threads} \
                    -r {chromosome}:{start}-{end} \
                    -i "-" --missing-to-ref -m none \
                    -O z0 -o {output_name} {merge_string} && {env_bin}tabix {output_name}')
            else:
                vt.contShell(f'\
                    mv {merge_string} {output_name} && {env_bin}tabix {output_name}')
        return output_name
    
    @subproc
    def deleteInputs(parent_future, vcf_list):
        pass
        # for path in vcf_list:
        #     os.remove(path)
        #     os.remove(f'{path}.tbi')


    def treeMerge(vcf_list, output_stub, genome_len, default_merge=10, genome_splits=2, threads=1, tmp_dir='tmp', output_dir='.', merge_dict={1: 20}):
        client = vt.findClient()
        final_futures = []
        rm_futures = []
        os.makedirs(tmp_dir, exist_ok=True)
        # get refrence name
        with pysam.VariantFile(vcf_list[0]) as open_vcf:
            chromosome = open_vcf.get_reference_name(0)
        # split the work by genome length
        starts, ends = np.asarray([
            np.linspace(0, genome_len, num=genome_splits + 1)[:-1] + 1,
            np.linspace(0, genome_len, num=genome_splits + 1)[1:]]).astype(int)
        split_stubs = [
            f'{tmp_dir}/{output_stub}.G{i + 1}' for i in range(genome_splits)]
        # iterate by genome splits
        with open(f'{output_dir}/merge_plan.log', 'w') as merge_plan:
            for i, (stub, start, end) in enumerate(zip(split_stubs, starts, ends)):
                print(f'preparing jobs at {start}-{end}')
                merge_plan.write(f'TREE {start} - {end} merges:\n')
                # iterate by merges
                current_paths, current_futures = vcf_list, vcf_list
                merge_it = 1
                merge_level = 1
                while len(current_paths) > 1:
                    next_paths, next_futures = [], []
                    try:
                        merge_size = merge_dict[merge_level]
                    except KeyError:
                        merge_size = default_merge
                    for i in range(0, len(current_paths), merge_size):
                        merge_name = f'{stub}.L{merge_level}.I{merge_it}.vcf.gz'
                        target_paths, target_futures = current_paths[i:i + merge_size], current_futures[i:i + merge_size]
                        # submit merge job
                        current_future = client.submit(
                            mergeJob, target_futures, merge_name, chromosome, start, end, merge_level,
                            threads=threads, priority=10)
                        if merge_level > 1:
                            # submit removal job, do not remove first level files
                            rm_futures.append(client.submit(
                                deleteInputs, current_future, target_paths, priority=0))
                        # write to log
                        merge_plan.write(
                            f'{" ".join(target_paths)} > {merge_name}\n')
                        next_paths.append(merge_name), next_futures.append(current_future)
                        merge_it += 1
                    merge_level += 1
                    current_paths, current_futures = next_paths, next_futures
                print(f'Merging via tree with {merge_it - 1} steps x {merge_level - 1} levels.')
                final_futures.append(current_future)
            wait(final_futures)
            wait(rm_futures)
        return final_futures
    

    @subproc
    def writeRecords(input_paths, output_path, start, end):
        # open miss array
        miss_array = zarr.open(
            f'{args.out_dir}/miss_array_by_pos.zarr', mode='r')
        # open output vcf
        vcf_out = pysam.VariantFile(
            output_path, 'w', header=pysam.VariantFile(input_paths[0]).header)
        written, skipped, outgroup_only = 0, 0, 0
        for in_path in input_paths:
            vcf_in = pysam.VariantFile(in_path)
            # fetch records in required region
            for record in vcf_in.fetch(contig=vcf_in.get_reference_name(0), start=start, end=end):
                # ignore records outside of defined region
                if record.pos - 1 < start or record.pos - 1 >= end:
                    continue
                # get ordered ref and miss array
                sample_ids = np.asarray([sid for sid in record.samples])
                is_ref = np.asarray([record.samples[sid]['GT'] == (0, 0) for sid in sample_ids])
                is_miss = miss_array[:, record.pos - 1]
                # edit miss positions in record
                for miss_idx in np.where(is_miss)[0]:
                    record.samples[sample_ids[miss_idx]]['GT'] = (None, None)
                if has_outgroup:
                    # if no non-miss, non-outgroup samples are variable, don't write
                    has_variants = np.all([
                        ~is_ref,  # must be alt
                        ~np.any([is_miss, is_outgroup], axis=0)],  # must be neither OG or miss
                        axis=0).any()
                    # check if only outgroup was variable
                    if np.sum(~is_ref) == 1 and np.all([~is_ref, is_outgroup], axis=0).any():
                        outgroup_only += 1
                else:
                    # if no non-miss samples are variable, don't write
                    has_variants = np.all([
                        ~is_ref,  # must be alt
                        ~is_miss],  # must not be miss
                        axis=0).any()
                if has_variants:
                    vcf_out.write(record)
                    written += 1
                else:
                    skipped += 1
            vcf_in.close()
        vcf_out.close()
        return written, skipped, outgroup_only

    # argument handling
    parser = argparse.ArgumentParser(
        description=""" """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required/global inputs
    parser.add_argument(
        '--input-csv', type=str, required=True, help='csv sample sheet')
    parser.add_argument(
        '--inputs-dir', type=str, default='.', help='base directory to look for files in.')
    parser.add_argument(
        '--out-dir', default='.', type=str, required=True,
        help='output directory')
    parser.add_argument(
        '--outgroup', type=str, default='', help='outgroup label (if any). multiple allowed with comma separation.')
    parser.add_argument(
        '--split-genome', type=int, default=2, help='number of files to split genome into (for parallelization, best to ~= n_workers).')
    parser.add_argument(
        '--merge-size', type=int, default=10, help='number of vcfs to merge at a time (for parallelization)')
    parser.add_argument(
        '--keep-intermediates', action='store_true')
    
    # multithreading options
    # RECOMMENDED TO RUN WITH AT LEAST 10 THREADS
    parser.add_argument(
        '--local-threads', type=int, default=1)
    parser.add_argument(
        '--n-workers', type=int, default=1)
    parser.add_argument(
        '--use-local', action='store_true')
    parser.add_argument(
        '--use-slurm', action='store_true')
    
    # SLURM options
    parser.add_argument(
        '--queue', type=str, default='sapphire')
    parser.add_argument(
        '--process-per-node', type=int, default=1)
    parser.add_argument(
        '--cores-per-process', type=int, default=1)
    parser.add_argument(
        '--memory-per-process', type=str, default='4GB')
    parser.add_argument(
        '--walltime', type=str, default='1:00:00')
    

    args = parser.parse_args()

    # start parallel workers with dask
    client = startClient(log_dir=args.out_dir, **args.__dict__)

    # parse sample sheet for errors
    print('Reading sample sheet....')
    sample_sheet_df = pd.read_csv(args.input_csv)
    if any(~np.isin(['label', 'vcf_path', 'miss_path'], sample_sheet_df.columns)):
        raise ValueError('One or more required columns in sample sheet is missing.')
    
    # generate directory structure
    tmp_dir = f'{args.out_dir}/tmp/'
    os.makedirs(f'{tmp_dir}', exist_ok=True)

    print(f'Zipping VCF files if not already ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # split records to monoallelics into tmp folder
    @subproc
    def copyVCF(label, file_path):
        vt.contShell(f'\
            bcftools norm {file_path} -m - -O z0 -o {tmp_dir}/{label}.vcf.gz && \
            tabix {tmp_dir}/{label}.vcf.gz')
        
    futures = []
    updated_vcf_path = []
    for i, rdata in sample_sheet_df.iterrows():
        futures.append(client.submit(
            copyVCF, rdata.label, f'{args.inputs_dir}/{rdata.vcf_path}'))
        updated_vcf_path.append(f'{tmp_dir}/{rdata.label}.vcf.gz')
    wait(futures)
    del futures
    sample_sheet_df.loc[:, 'vcf_path'] = updated_vcf_path

    print(f'Merging data for {len(sample_sheet_df.vcf_path)} samples ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # initialize zarr array tracking miss locations
    genome_len = zarr.open(
        f"{args.inputs_dir}/{sample_sheet_df.loc[sample_sheet_df.index[0], 'miss_path']}", mode='r').shape[0]
    miss_array = zarr.open(
            f'{tmp_dir}/miss_array_by_strain.zarr', mode='w',
            shape=(len(sample_sheet_df), genome_len),
            chunks=(1, genome_len),
            dtype=bool)

    # WRITE MISS ARRAY
    # (1) write miss array samples (row) to genome (columns) length arrays
    # (2) rechunk miss array to be read by columns (genomic positions)

    # prepare function to parallelize
    @subproc
    def writeMissArray(arr_idx, miss_path):
        # apply missing sites to the fasta
        miss_array[arr_idx, :] = zarr.open(miss_path, mode='r')[:]

    print(f'Writing miss array sample-wise ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # parallel write
    futures = []
    for i, miss_path in enumerate(sample_sheet_df.miss_path):
        futures.append(
            client.submit(writeMissArray, i, f'{args.inputs_dir}/{miss_path}'))
    wait(futures)
    del futures

    print(f'Rechunking miss array position-wise ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # parallel rechunk
    rechunked = rechunk(
        miss_array,
        target_chunks=(miss_array.shape[0], 1000),
        target_store=f'{args.out_dir}/miss_array_by_pos.zarr',
        max_mem=f'1GB',
        temp_store=f'{tmp_dir}/intermediate.zarr')
    rechunked.execute(
        num_workers=args.n_workers)


    print(f'Merging input VCFs ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # MERGE INPUT VCFS BY SAMPLE
    # (1) split VCFs across genome in (genome_len / genome_splits) sized blocks
    # (2) in order, merge (merge_size) groups of VCFs to generate a tree of merges ending with a single file for each genome split 
    merged_vcfs = treeMerge(
        vcf_list=sample_sheet_df.vcf_path.values,
        output_stub='raw_merge',
        default_merge=args.merge_size,
        genome_len=genome_len,
        genome_splits=args.split_genome,
        tmp_dir=tmp_dir, output_dir=args.out_dir)


    print(f'Editing merged VCF to include miss locations ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # WRITE MISSES AND FILTER VCF
    # for each vcf position, store miss values for strains with missing data in miss array
    # ignore sites that are ONLY variant in the outgroup (or not variant)
    # check for outgroup
    if args.outgroup != '':
        print(f'Searching for outgroup {args.outgroup}....')
        is_outgroup = np.isin(sample_sheet_df.label, args.outgroup.split(','))
        if np.any(is_outgroup) == False:
            raise ValueError('Outgroup not found in label values.')
        else:
            print('outgroup found!')
        has_outgroup = True
    else:
        has_outgroup = False

    # split the work equally across workers
    starts, ends = np.asarray([
        np.linspace(0, genome_len, num=args.n_workers + 1)[:-1],
        np.linspace(0, genome_len, num=args.n_workers + 1)[1:]]).astype(int)

    # execute in parallel (last parallel step, block here)
    futures = []
    files_to_merge = []
    for i, (s, e) in enumerate(zip(starts, ends)):
        output_path = f'{tmp_dir}/added_miss.{i}.vcf'
        futures.append(
            client.submit(writeRecords, merged_vcfs, output_path, s, e))
        files_to_merge.append(output_path)
    dask_output = []
    for f in as_completed(futures):
        dask_output.append(client.gather(f))
        f.release()

    # delete vcf merge futures
    del merged_vcfs

    print('Shutting down dask workers (no more parallel steps).')
    client.shutdown()
    print('\n')


    print(f'Merging files split files into final record ({int(timeit.default_timer() - start_time)}s elapsed)....')
    # merge files
    with open(f'{tmp_dir}/concat_list.txt', 'w') as fout:
        for fp in files_to_merge:
            fout.write(f'{fp}\n')
    vt.contShell(f'\
        {env_bin}bcftools concat --threads {args.local_threads} -f {tmp_dir}/concat_list.txt -O z9 -o {args.out_dir}/merged.filtered.vcf.gz && \
        {env_bin}tabix {args.out_dir}/merged.filtered.vcf.gz')
    if has_outgroup:
        written, skipped, outgroup_only = np.asarray(
            dask_output).sum(axis=0)
        print(f'{written} lines written, {skipped} lines skipped, {outgroup_only} outgroup only (skipped).')
    else:
        written, skipped, _ = np.asarray(
            dask_output).sum(axis=0)
        print(f'{written} lines written, {skipped} lines skipped.')
    if args.keep_intermediates is False:
        vt.contShell(f'rm -r {tmp_dir}')
    
    elapsed = timeit.default_timer() - start_time
    print(f'finished in {int(elapsed)}s or {int(elapsed/60)}m')
    sys.exit(0)