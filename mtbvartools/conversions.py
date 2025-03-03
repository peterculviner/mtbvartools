import pysam, shutil, os, tqdm, dendropy, shutil
from pastml.acr import pastml_pipeline
import numpy as np
import pandas as pd
from time import sleep
from .KeyedByteArray import KeyedByteArray
from .CallBytestream import CallBytestream
from .trees import writeLabelledTree
from .dasktools import subproc, timedsubproc, findClient
from .misc import getSlices, getNearbyGenes

def writeVariantCalls(
        target_vcf, output_path,
        sample_list=None, compression='zlib', ckwargs={}, n_jobs=(1000, 50),
        output_lookup={
            (0, 0): 0,
            (1, 1): 1, (0, 1): 1, (1, 0): 1,
            (None, None): 2}):
    # define single job
    def job(target_vcf, sample_list, region_start, region_end, cfunc, ckwargs):
        if sample_list == None:  # if none provided, write all samples in VCF order; otherwise iterate across sample list
            with pysam.VariantFile(target_vcf) as open_vcf:
                sample_list = open_vcf.header.samples
        pointers = {}
        # open vcf and iterate through
        open_vcf = pysam.VariantFile(target_vcf)
        region_string = f'{open_vcf.get_reference_name(0)}:{region_start}-{region_end}'  # 1-indexed, inclusive on both ends
        var_i = 0
        pointer_i = 0
        total_bytes = 0
        total_bytes_written = 0
        bytestream = b''
        for record in open_vcf.fetch(region=region_string):
            if record.pos >= region_start and record.pos <= region_end:
                record_calls = np.zeros(len(sample_list), dtype='uint8')
                for i, name in enumerate(sample_list):
                    try:
                        sample_call = record.samples[name]['GT']
                        record_calls[i] = output_lookup[sample_call]
                    except KeyError:
                        raise KeyError(
                            f'Unhandled call: {sample_call}, for {name} at {record}.')
                # save bytestream info
                uncompressed_bytes = bytearray(record_calls)
                compressed_bytes = cfunc(uncompressed_bytes, **ckwargs)
                bytestream += compressed_bytes
                # save pointers
                pointers[(record.pos, record.ref, record.alts[0])] = (pointer_i, len(compressed_bytes))
                # tick forward counters
                total_bytes += len(uncompressed_bytes)
                total_bytes_written += len(compressed_bytes)
                pointer_i += len(compressed_bytes)
                var_i += 1
                if len(record.alts) > 1:
                    raise ValueError(f'Unhandled multiple alts for record {str(record)}')
        return bytestream, pointers, total_bytes, total_bytes_written
    # prepare directory
    shutil.rmtree(output_path, ignore_errors=True)
    os.makedirs(output_path)
    # handle job
    client = findClient()
    # get variant file information
    with pysam.VariantFile(target_vcf) as open_vcf:
        reference_length = open_vcf.header.contigs[open_vcf.get_reference_name(0)].length
        samples = np.asarray(open_vcf.header.samples)
    # initialize bytestream
    variant_kba = KeyedByteArray(
        f'{output_path}/by_variant.kba', mode='w', columns=samples, dtype='uint8',
        compression=compression, compression_kwargs=ckwargs)
    # split based on jobs
    regions = zip(
        np.linspace(1, reference_length + 1, n_jobs[0] + 1).astype(int)[:-1],
        np.linspace(1, reference_length + 1, n_jobs[0] + 1).astype(int)[1:] - 1)
    # submit futures
    futures = []
    for i, (region_start, region_end) in enumerate(regions):
        futures.append(client.submit(
            job, target_vcf, sample_list, region_start, region_end,
            variant_kba.compress, ckwargs,
            priority=n_jobs[0]-i))
    # gather futures in order submitted
    print('importing calls from VCF (indexed & compressed by variants)....')
    for f in tqdm.tqdm(futures):
        bytestream, pointers, total_bytes, _ = client.gather(f)
        f.release()
        variant_kba.writebinary(
            pointers, bytestream, raw_bytes=total_bytes)
    variant_kba.close() # close and write key to file
    # rechunk
    variant_kba = KeyedByteArray(f'{output_path}/by_variant.kba', mode='r')
    variant_kba.rechunk(f'{output_path}/by_node.kba', jobs=n_jobs[1])


def writeAncestorCalls(
        vcb_path, tree_path, output_path,
        step_size=100, method='DOWNPASS', miss_filter=0.50, rechunk_jobs=None, timeout=None):

    def acrJobHandler(job_slices, vcb_path, output_path, method, miss_filter):
        # try with a timeout
        try:
            return [acrJob(job_slices, vcb_path, output_path, method, miss_filter)]
        except TimeoutError:
            shutil.rmtree(  # remove timed out directory
                f'{output_path}/AR_TMP_{job_slices[0]}_{job_slices[1]}')
            singlet_output = []
            for left in range(job_slices[0], job_slices[1]):
                try:
                    singlet_output.append(
                        acrJob((left, left + 1), vcb_path, output_path, method, miss_filter))
                except TimeoutError:
                    singlet_output.append(
                        RuntimeWarning(f'Handled timeout at index {left}'))
                    shutil.rmtree(  # remove timed out directory
                        f'{output_path}/AR_TMP_{left}_{left + 1}')
            return singlet_output
    
    # define ancestral reconstruction job
    @timedsubproc(timeout)
    def acrJob(job_slices, vcb_path, output_path, method, miss_filter):
        tmp_dir = f'{output_path}/AR_TMP_{job_slices[0]}_{job_slices[1]}'
        os.makedirs(tmp_dir, exist_ok=True)
        # open call byestream to access data
        target_vcb = CallBytestream(vcb_path, init_nodes=False)
        input_calls = target_vcb.calls.iloc[job_slices[0]:job_slices[1]]
        input_rows = target_vcb.calls.row[job_slices[0]:job_slices[1]]
        input_cols = target_vcb.calls.col
        pass_filter = np.sum(  # remove rows with too many misses
            input_calls == 2, axis=1) / len(input_cols) <= miss_filter
        if np.all(pass_filter == False):  # none of the sites passed filter
            shutil.rmtree(tmp_dir)
            return None
        # define variant labels
        var_idxs = np.arange(len(input_rows))[pass_filter]
        input_df = pd.DataFrame(  # define dataframe and write to CSV
            data=input_calls[pass_filter],
            index=var_idxs,
            columns=input_cols,
            dtype=str).T
        input_df.values[input_df == '2'] = ''  # replace '2' with pastml missing character ''
        # input_df.index.name = 'ID'
        input_df.to_csv(f'{tmp_dir}/input_calls.csv')
        pastml_pipeline(  # run pastml
            f'{output_path}/tree.nwk',
            data=f'{tmp_dir}/input_calls.csv',
            prediction_method=method,
            data_sep=',',
            work_dir=tmp_dir,
            threads=1)
        node_states_df = pd.read_csv(  # read pastml
            f'{tmp_dir}/combined_ancestral_states.tab',
            delimiter='\t',
            index_col=0)
        node_states_df.columns = node_states_df.columns.astype(int)
        node_states_df = node_states_df.loc[:, np.sort(node_states_df.columns)]
        # identify ambiguous nodes (shown as more than one row)
        node, counts = np.unique(node_states_df.index, return_counts=True)
        keep_mask = np.ones(len(node_states_df)).astype(bool)
        for ambiguous_node in node[counts > 1]:
            # set ambiguous nodes to miss character
            node_states_df.loc[
                ambiguous_node,
                node_states_df.columns[np.all(np.isnan(node_states_df.loc[ambiguous_node]) == False, axis=0)]] = 2
            # update keep mask to remove second
            keep_mask[np.where(node_states_df.index == ambiguous_node)[0][1]] = False
        node_states_df = node_states_df.loc[keep_mask].T.astype(target_vcb.calls.dtype)
        # store pastml outputs
        output_df = []
        for i in node_states_df.index:
            output_df.append(
                pd.read_csv(
                    f'{tmp_dir}/params.character_{i}.method_{method}.tab',
                    delimiter='\t', index_col=0).value)
        output_df = pd.concat(output_df, axis=1).T
        output_df.index = var_idxs + job_slices[0]
        shutil.rmtree(tmp_dir)
        # prepare outputs
        pointers, compressed_bytes, raw_byte_count = target_vcb.calls.preparebinary(
            [input_rows[i] for i in np.where(pass_filter)[0]],
            node_states_df.values)
        col_values = node_states_df.columns.values
        target_vcb.close()
        return col_values, output_df.to_csv(), pointers, compressed_bytes, raw_byte_count
    
    # prepare directory
    shutil.rmtree(output_path, ignore_errors=True)
    os.makedirs(output_path)
    # prepare client if not already prepared
    client = findClient()

    # write the labelled tree
    writeLabelledTree(
        tree_path,
        f'{output_path}/tree.nwk',
        suppress_internal_node_labels=True, suppress_rooting=True)  # req'd for pastml read-in
    
    # load target vcb to initialize work and outputs
    target_vcb = CallBytestream(vcb_path, init_calls=True, init_nodes=False)
    job_slices = getSlices(step_size, len(target_vcb.calls.row))
    futures = []  # prepare futures
    for i, js in enumerate(job_slices):
        futures.append(client.submit(
            acrJobHandler, js, vcb_path, output_path,
            method, miss_filter, priority=len(job_slices)-i))
        sleep(0.001)
    with open(f'{output_path}/acr_summary.csv', 'w') as summary_file:
        # handle first future
        first_columns, summary_csv, _, _, _ = futures[0].result()[0]
        summary_file.write(
            summary_csv[:summary_csv.find('\n') + 1])
        output_kba = KeyedByteArray(  # initialize output
            f'{output_path}/by_variant.kba', mode='w', columns=first_columns, dtype=target_vcb.calls.dtype,
            compression=target_vcb.calls.index['compression']['compression_type'],
            compression_kwargs=target_vcb.calls.index['compression']['kwargs'])
        # gather outputs as futures complete
        for js, f in tqdm.tqdm(list(zip(job_slices, futures))):
            results = f.result()  # handle list of results in future
            for r in results:
                if r is None:  # handle situation with all missing results
                    continue
                elif isinstance(r, Exception):
                    print(r)
                else:
                    columns, summary_csv, pointers, compressed_bytes, raw_byte_count = r
                    # write the pre-compressed data
                    output_kba.writebinary(pointers, compressed_bytes, raw_byte_count)
                    # append to the summary file
                    summary_file.write(
                        summary_csv[summary_csv.find('\n') + 1:])
                    if np.any(first_columns == columns) is False:
                        raise ValueError('columns do not match')
            f.release()
    output_kba.close()
    if rechunk_jobs != None:
        # rechunk output_kba
        output_kba = KeyedByteArray(f'{output_path}/by_variant.kba', mode='r')
        output_kba.rechunk(f'{output_path}/by_node.kba', jobs=rechunk_jobs)
        output_kba.close()
    # add pos / ref / alt index to summary csv
    summary_df = pd.read_csv(
        f'{output_path}/acr_summary.csv', index_col=0)
    ancestor_calls = CallBytestream(
        output_path, init_calls=True, init_nodes=False)
    summary_df.index = pd.MultiIndex.from_tuples(
        ancestor_calls.calls.row, names=['pos', 'ref', 'alt'])
    summary_df.to_csv(f'{output_path}/acr_summary.csv')


def writeEventTransitions(ancestor_calls, output_path, tree_obj, rechunk_jobs=None):
    os.makedirs(output_path, exist_ok=True)
    by_node_kba = KeyedByteArray(
        f'{output_path}/by_node.kba', mode='w', columns=ancestor_calls.nodes.col, dtype='uint8',
        compression=ancestor_calls.nodes.index['compression']['compression_type'],
        compression_kwargs=ancestor_calls.nodes.index['compression']['kwargs'])
    # initialize arrays to record transitions
    transitions_to_ref = np.zeros(len(ancestor_calls.nodes.col), dtype=int)
    transitions_to_alt = np.zeros(len(ancestor_calls.nodes.col), dtype=int)
    # apply well-defined node transitions
    print('writing transitions....')
    for pnode in tqdm.tqdm(list(tree_obj.postorder_internal_node_iter())):
        pnode_state = ancestor_calls.nodes.loc[pnode.label]
        for cnode in pnode.child_nodes():
            # get transition types
            cnode_state = ancestor_calls.nodes.loc[cnode.label]
            tta_mask = np.all(  # ref parent -> alt child
                [pnode_state == 0, cnode_state == 1], axis=0)
            ttr_mask = np.all(  # alt parent -> ref child
                [pnode_state == 1, cnode_state == 0], axis=0)
            # print(transitions_to_ref.shape, tta_mask.shape)
            transitions_to_alt += tta_mask
            transitions_to_ref += ttr_mask
            # prepare output vector
            output_vector = np.zeros(cnode_state.shape, dtype='uint8')
            output_vector[tta_mask] = 1  # transition to ALT is 1
            output_vector[ttr_mask] = 2  # transition to REF is 2
            by_node_kba.write(cnode.label, output_vector)
    by_node_kba.close()
    if rechunk_jobs != None:
        # rechunk output_kba
        output_kba = KeyedByteArray(f'{output_path}/by_node.kba', mode='r')
        output_kba.rechunk(f'{output_path}/by_variant.kba', jobs=rechunk_jobs)
        output_kba.close()
    return transitions_to_ref, transitions_to_alt


def countReversionEvents(node_kba_path, tree_obj):
    by_node_kba = KeyedByteArray(  # open the array in read mode
        node_kba_path, mode='r')
    n_variants = len(by_node_kba.col)
    # initialize arrays to record reversions
    reversions_to_ref = np.zeros(n_variants, dtype=int)
    reversions_to_alt = np.zeros(n_variants, dtype=int)
    # search for reversions on the tree / event data
    print('counting reversions....')
    for tnode in tqdm.tqdm(list(tree_obj.postorder_node_iter())):
        if tnode == tree_obj.seed_node:
            continue  # skip seed node, by definition cannot have transitions
        tnode_rtr = np.zeros(n_variants, dtype=int)
        tnode_rta = np.zeros(n_variants, dtype=int)
        tnode_state = by_node_kba.loc[tnode.label]  # target node state vector
        qnode = tnode
        while qnode != tree_obj.seed_node:
            # get query node state vector
            qnode_state = by_node_kba.loc[qnode.label]
            # compare to target node state vector
            rtr_mask = np.all(  # rev to ref req: query TTA (1), target TTR (2)
                [qnode_state == 1, tnode_state == 2], axis=0)
            rta_mask = np.all(  # rev to alt req: query TTR (2), target TTA (1)
                [qnode_state == 2, tnode_state == 1], axis=0)
            # record as much as 1 reversion event per variant
            tnode_rtr[rtr_mask] = 1
            tnode_rta[rta_mask] = 1
            qnode = qnode.parent_node
        # add reversions to cumulative sum
        reversions_to_ref += tnode_rtr
        reversions_to_alt += tnode_rta
    by_node_kba.close()
    return reversions_to_ref, reversions_to_alt


def writeEventCalls(acr_path, var_path, output_path, annotations_csv=None, gene_table_csv=None, rechunk_jobs=None):
    # load in call bytestreams
    ancestor_calls = CallBytestream(acr_path)
    # load in tree object and pastml output
    tree_obj = dendropy.Tree.get(
        path=f'{acr_path}/tree.nwk',
        schema='newick',
        preserve_underscores=True)
    # copy leaf taxon labels -> leaf label
    for leaf in tree_obj.leaf_nodes():
        leaf.label = leaf.taxon.label
    acr_output_df = pd.read_csv(
        f'{acr_path}/acr_summary.csv', index_col=[0, 1, 2])
    # initialize VCB at output path
    os.makedirs(output_path, exist_ok=True)
    ttr, tta = writeEventTransitions(  # write event transitions to keyed byte array
        ancestor_calls, output_path, tree_obj, rechunk_jobs=rechunk_jobs)
    rtr, rta = countReversionEvents(  # count reversions
        f'{output_path}/by_node.kba', tree_obj)
    variant_calls = CallBytestream(var_path)
    call_counts = []
    for vidx in ancestor_calls.calls.row:
        call_vector = variant_calls.calls.loc[vidx]
        call_counts.append([
            np.sum(call_vector == 0), np.sum(call_vector == 1), np.sum(call_vector == 2)])
    call_counts = pd.DataFrame(
        data=call_counts, columns=['ref_count', 'alt_count', 'miss_count'])
    variant_calls.close()
    ancestor_calls.close()
    # prepare output df
    acr_output_df.loc[:, 'ref_count'] = call_counts.ref_count.values
    acr_output_df.loc[:, 'transitions_to_ref'] = ttr
    acr_output_df.loc[:, 'reversions_to_ref'] = rtr
    acr_output_df.loc[:, 'alt_count'] = call_counts.alt_count.values
    acr_output_df.loc[:, 'transitions_to_alt'] = tta
    acr_output_df.loc[:, 'reversions_to_alt'] = rta
    acr_output_df.loc[:, 'miss_count'] = call_counts.miss_count.values
    # write files to vcb
    acr_output_df.to_csv(f'{output_path}/event_calls.raw.csv')
    shutil.copy(f'{acr_path}/tree.nwk', f'{output_path}/tree.nwk')
    # write other event calls
    if annotations_csv != None:
        print('merging gene-level annotations....')
        # merge onto vcf_annotations to count vcf annotations in file
        annotated_events = pd.merge(
            acr_output_df,
            pd.read_csv(annotations_csv, index_col=['pos', 'ref', 'alt']),
            how='right',
            left_index=True, right_index=True)
        annotated_events.to_csv(f'{output_path}/event_calls.annotated.csv')
    if gene_table_csv != None:
        print('writing genomic annotations....')
        gene_table = pd.read_csv(gene_table_csv, index_col=0)
        nearby_data = []
        for pos in tqdm.tqdm(acr_output_df.index.get_level_values('pos')):
            nearby_data.append(getNearbyGenes(pos - 1, gene_table))
        nearby_df = pd.DataFrame(
            index=acr_output_df.index,
            data=nearby_data)
        pd.merge(
            acr_output_df, nearby_df,
            left_index=True, right_index=True).to_csv(f'{output_path}/event_calls.genomic.csv')
