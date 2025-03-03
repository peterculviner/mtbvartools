#!/usr/bin/env -S python -u

import sys, argparse, os, timeit, pysam, zarr, warnings
import pandas as pd
import numpy as np
from dask.array import from_zarr
from dask.distributed import LocalCluster
import mtbvartools as vt
from mtbvartools.dasktools import subproc

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
    parser.add_argument(
        '--min-strain-count', type=int, default=1, help='minimum number of strains required to be variant at a position for position to be included in fasta output.')
    
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

    # prepare function to parallelize
    @subproc
    def preparefasta(arr_idx, label, vcf_path, miss_path):
        output = f'{label}:\n'
        # write a SNP-only VCF
        output += vt.contShell(f'\
            bcftools view {vcf_path} -m2 -M2 -v snps | \
            awk \'{{gsub(/NC_000962/,"NC_000962.3"); print}}\' - | bcftools view -o {tmp_dir}/{label}.tmp.vcf.gz -O z && tabix {tmp_dir}/{label}.tmp.vcf.gz',
            is_return=True)
        # edit the input fasta with the variants
        output += vt.contShell(f'\
            bcftools consensus -f {args.input_fasta} -H 1 {tmp_dir}/{label}.tmp.vcf.gz > {tmp_dir}/{label}.tmp.fasta',
            is_return=True)
        # apply missing sites to the fasta
        miss_indexes = zarr.open(miss_path, mode='r')[:]  # load miss zarr
        fasta_file = pysam.FastaFile(
            f'{tmp_dir}/{label}.tmp.fasta')
        fasta_array = np.asarray(
            [bp for bp in fasta_file.fetch(fasta_file.references[0])])
        fasta_array[miss_indexes] = 'N'
        with open(f'{tmp_dir}/{label}.tmp.fasta', 'w') as f:
            f.write(f'> {fasta_file.references[0]}\n{"".join(fasta_array)}')
        output += vt.contShell(f'\
            rm {tmp_dir}/{label}.tmp.fasta.fai', is_return=True)
        # store missing sites in large zarr
        miss_array[arr_idx, :] = miss_indexes
        # store the variant positions - ignore possible variants in a sample's miss regions
        variant_file = pysam.VariantFile(
            f'{tmp_dir}/{label}.tmp.vcf.gz')
        is_variant = np.zeros(genome_len).astype(bool)
        var_pos = [v.pos - 1 for v in variant_file.fetch()]  # vcf files are read in as 1-idx
        is_variant[var_pos] = True
        is_variant[miss_indexes] = False  # set any variants in miss regions back to False
        var_array[arr_idx, :] = is_variant
        # apply global mask
        if args.mask is not None:
            output += 'Applied mask.\n'
            output += vt.contShell(f'\
                bedtools maskfasta -fi {tmp_dir}/{label}.tmp.fasta -fo {tmp_dir}/{label}.fasta \
                -bed {args.mask} -mc X', is_return=True)
        else:
            output += vt.contShell(f'\
                mv {tmp_dir}/{label}.tmp.fasta {tmp_dir}/{label}.fasta',
                is_return=True)
        return output

    print('Preparing unfiltered fastas & merging miss locations....')
    # write separate fasta files in parallel using above function
    # (1) write a SNP only VCF
    # (2) edit the genome using bcftools consensus
    # (3) apply missing sites using the missing masks (stored in zarr files)
    # (4) merge missing sites into a single zarr file for summing
    # (5) store which positions were variant in a zarr file for summing
    # (6) apply a BED mask excluding particular regions (if given)
    # (7) write a fasta file for this sample to the temporary folder
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

    print('Generating masks....')
    # MISSING SAMPLE HANDLING
    # sum the number of missing positions across all samples
    miss_sum = from_zarr(f'{args.out_dir}/{args.output}.miss.zarr').sum(axis=0)
    fraction_miss = miss_sum.compute() / len(sample_sheet_df)

    # generate a global missing position mask where positions are less than threshold
    pass_miss_mask = fraction_miss < args.miss_threshold
    print(f'{pass_miss_mask.sum()}/{genome_len} sites pass miss threshold of {args.miss_threshold}')
    results_dict['genome_len'] = genome_len
    results_dict['miss_threshold'] = args.miss_threshold
    results_dict['pass_miss_threshold_sites'] = pass_miss_mask.sum()


    # VARIANT POSITION HANDLING
    # run var sum in parallel
    var_sum = from_zarr(f'{args.out_dir}/{args.output}.variants.zarr').sum(axis=0)
    n_alt = var_sum.compute()  # calculate the number of strains with alt call at this position

    # generate masks
    # (1) has any variants in samples -> variant_mask
    # (2) has enough variants to meet minimum variant counts -> variant_output_mask
    # in default value of 1, these two masks will be identical
    variant_mask = np.all([
        n_alt > 0,
        n_alt < len(sample_sheet_df)], axis=0)
    variant_output_mask = np.all([
        n_alt >= args.min_strain_count,
        n_alt <= (len(sample_sheet_df) - args.min_strain_count)], axis=0)
    # print statistics
    print(f'{variant_mask.sum()} sites have variants in at least 1 sample.')
    print(f'{variant_output_mask.sum()} sites have variants in at least {args.min_strain_count} sample(s), and will be considered for output.')
    results_dict['variant_sites'] = variant_mask.sum()
    results_dict['minimum_variant_strains_to_consider'] = args.min_strain_count
    results_dict['considered_variant_sites'] = variant_output_mask.sum()


    # BED MASK HANDLING
    pass_bed_mask = np.ones(genome_len).astype(bool)
    if args.mask is not None:
        # get sites that were in the global mask - mask is the same for all files
        fasta_file = pysam.FastaFile(f'{tmp_dir}/{sample_sheet_df.label[0]}.fasta')
        fasta_str = fasta_file.fetch(fasta_file.references[0])
        # convert masked sites ("X") from indexes to a boolean mask
        bed_mask_idx = np.asarray([gidx for gidx, s in enumerate(fasta_str) if s == 'X'])
        pass_bed_mask[bed_mask_idx] = False
        print(f'Using BED mask, masked {len(bed_mask_idx)} positions.')
        results_dict['bed_mask_sites'] = len(bed_mask_idx)
    else:
        print('No BED mask provided, so no positions masked.')
        results_dict['bed_mask_sites'] = 0


    # PRINT INVARIANT, VARIANT, and OUTPUT COUNTS
    passing_invariant_mask = np.all([
        ~variant_mask,  # not variant
        pass_miss_mask,  # not miss
        pass_bed_mask,  # not BED masked
        ], axis=0)
    print(
        f'{passing_invariant_mask.sum()} sites passing miss and BED masks were invariant.')
    passing_variant_mask = np.all([
        variant_mask,  # is variant
        pass_miss_mask,  # not miss
        pass_bed_mask,  # not BED masked
        ], axis=0)
    print(
        f'{passing_variant_mask.sum()} sites passing miss and BED masks were variant in at least 1 sample.')
    passing_output_mask = np.all([
        variant_output_mask,  # meets required variant counts
        pass_miss_mask,  # not miss
        pass_bed_mask,  # not BED masked
        ], axis=0)
    print(
        f'{passing_output_mask.sum()} sites passing miss and BED masks were variant in at least {args.min_strain_count} sample(s) and will be in output.')
    results_dict['passing_variant_sites'] = passing_variant_mask.sum()
    results_dict['passing_output_sites'] = passing_output_mask.sum()
    results_dict['passing_invariant_sites'] = passing_invariant_mask.sum()
    # categorize invariant sites by base
    fasta_file = pysam.FastaFile(
        args.input_fasta)
    invariant_array = np.asarray(
        [bp for bp in fasta_file.fetch(fasta_file.references[0])])[passing_invariant_mask]
    results_dict['invariant_A'] = np.sum(invariant_array == 'A')
    results_dict['invariant_T'] = np.sum(invariant_array == 'T')
    results_dict['invariant_G'] = np.sum(invariant_array == 'G')
    results_dict['invariant_C'] = np.sum(invariant_array == 'C')

    print('Writing filtered fasta outputs....')
    # WRITE OUTPUT AND RECORD MASK
    # save kept mask to disk for use with future samples
    output_mask_zarr = zarr.open(
        f'{args.out_dir}/{args.output}.mask.zarr', mode='w',
        shape=(genome_len),
        chunks=(genome_len),
        dtype=bool)
    output_mask_zarr[:] = passing_output_mask

    # define a function for trimming fastas
    @subproc
    def trimfasta(label):
        fasta_file = pysam.FastaFile(f'{tmp_dir}/{label}.fasta')
        # convert to numpy array for indexing
        fasta_array = np.fromiter(
            fasta_file.fetch(fasta_file.references[0]), dtype='U1')
        # open output mask
        output_mask = zarr.open(
            f'{args.out_dir}/{args.output}.mask.zarr', mode='r')[:]
        trimmed_string = ''.join(fasta_array[output_mask])
        with open(f'{tmp_dir}/{label}.trimmed.fasta', 'w') as fo:
            fo.write(f'>{label}\n{trimmed_string}\n')
        return label

    # write a trimmed fasta for each sample
    futures = [
        client.submit(trimfasta, label)
        for label in sample_sheet_df.label]
    for f in futures:
        client.gather(f)
        f.release()

    print('Merging fasta outputs and cleaning up....')
    vt.contShell(f'\
        for file in {tmp_dir}/*.trimmed.fasta; do cat "$file" >> {args.out_dir}/{args.output}.fasta; done')
    vt.contShell(f'\
        rm -r {tmp_dir}')
    
    print('Writing results summary....')
    pd.Series(
        data=results_dict,
        name=args.output).to_csv(f'{args.out_dir}/{args.output}.results.csv')
    
    elapsed = timeit.default_timer() - start_time
    print(f'finished in {int(elapsed)}s or {int(elapsed/60)}m')
    client.shutdown()
    sys.exit(0)