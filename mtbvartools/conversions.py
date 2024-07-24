import pysam, shutil, os, tqdm, dendropy
from pastml.acr import pastml_pipeline
import numpy as np
import pandas as pd
from .KeyedByteArray import KeyedByteArray
from .CallBytestream import CallBytestream
from .trees import writeLabelledTree
from .dasktools import subproc, findClient
from .misc import getSlices

def writeVariantCalls(target_vcf, output_path, sample_list=None, compression='zlib', ckwargs={}, n_jobs=(1000, 50)):
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
                    sample_call = record.samples[name]['GT']
                    if sample_call == (0, 0):
                        continue  # default is ref (0), no need to change
                    elif sample_call == (1, 1):
                        record_calls[i] = 1  # alt is saved as 1
                    elif sample_call == (None, None):
                        record_calls[i] = 2  # miss is saved as 2
                    else:
                        raise ValueError(
                            f'Unhandled call, {record["GT"]}, for sample {record.name}.')
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
        step_size=100, method='DOWNPASS', miss_filter=0.80, keep_cols=['num_unresolved_nodes', 'steps'],
        compression='zlib', ckwargs={}):
    # define ancestral reconstruction job
    @subproc
    def job(job_slices, vcb_path, output_path, method, miss_filter, keep_cols):
        tmp_dir = f'{output_path}/AR_TMP_{job_slices[0]}'
        os.makedirs(tmp_dir, exist_ok=True)
        with CallBytestream(vcb_path, init_nodes=False) as target_vcb:
            input_calls = target_vcb.calls.iloc[job_slices[0]:job_slices[1]]
            input_rows = target_vcb.calls.row[job_slices[0]:job_slices[1]]
            input_cols = target_vcb.calls.col
        pass_filter = np.sum(  # remove rows with too many misses
            input_calls == 2, axis=1) / len(input_cols) <= miss_filter
        # define variant labels
        var_idxs = np.arange(len(input_rows))[pass_filter]
        input_df = pd.DataFrame(  # define dataframe and write to CSV
            data=input_calls[pass_filter],
            index=var_idxs,
            columns=input_cols,
            dtype=str).replace('2', '').T
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
        node_states_df = node_states_df.loc[keep_mask].T.astype(int)
        # store pastml outputs
        output_df = []
        for i in node_states_df.index:
            output_df.append(
                pd.read_csv(
                    f'{tmp_dir}/params.character_{i}.method_{method}.tab',
                    delimiter='\t', index_col=0).value)
        output_df = pd.concat(output_df, axis=1).T.loc[:, keep_cols]
        output_df.index = var_idxs + job_slices[0]
        # write a temporary KBA, then read the bytestream and pointers
        tmp_kba = KeyedByteArray(
            f'{output_path}/kba{job_slices[0]}.tmp.kba', mode='w', dtype='uint8', columns=node_states_df.columns)
        for i, data in node_states_df.iterrows():
            tmp_kba.write(input_rows[i], data.values.astype('uint8'))
        tmp_kba.close()
        # read in the bytestream for output
        tmp_kba = KeyedByteArray(
            f'{output_path}/kba{job_slices[0]}.tmp.kba')
        pointers = tmp_kba.row_dict
        bytestream = tmp_kba.file.read()
        raw_bytes = tmp_kba.index['compression']['bytes']
        columns = tmp_kba.col
        tmp_kba.close()
        shutil.rmtree(tmp_dir) # remove temporary files
        os.remove(f'{output_path}/kba{job_slices[0]}.tmp.kba')
        os.remove(f'{output_path}/kba{job_slices[0]}.tmp.kba.kbi')
        return columns, output_df, pointers, bytestream, raw_bytes
    # prepare directory
    shutil.rmtree(output_path, ignore_errors=True)
    os.makedirs(output_path)
    # prepare client if not already prepared
    client = findClient()
    # load target vcb to split work
    with CallBytestream(vcb_path, init_calls=True, init_nodes=False) as target_vcb:
        job_slices = getSlices(step_size, len(target_vcb.calls.row))
    # write the labelled tree
    writeLabelledTree(
        tree_path,
        f'{output_path}/tree.nwk',
        suppress_internal_node_labels=True, suppress_rooting=True)  # req'd for pastml read-in
    futures = []  # submit futures
    for i, js in enumerate(job_slices):
        futures.append(client.submit(
            job, js, vcb_path, output_path,
            method, miss_filter, keep_cols,
            priority=len(job_slices)-i))
    # return futures
    # gather futures and release memory
    # initialize data array
    init_columns, _, _, _, _ = client.gather(futures[0])
    variant_kba = KeyedByteArray(  # initialize output
        f'{output_path}/by_variant.kba', mode='w', columns=np.asarray(init_columns), dtype='uint8',
        compression=compression, compression_kwargs=ckwargs)
    # iterate through futures 
    pastml_outputs = []
    for f in tqdm.tqdm(futures):
        columns, output_df, pointers, bytestream, raw_bytes = client.gather(f)
        if not np.all(init_columns == columns):  # verify columns are the same
            raise ValueError('Columns for datasets do not match!')
        f.release()
        variant_kba.writebinary(
            pointers, bytestream, raw_bytes=raw_bytes)
        pastml_outputs.append(output_df)
    # finish output df
    pastml_outputs = pd.concat(pastml_outputs, axis=0)
    pastml_outputs.index = pd.Index(variant_kba.pointers.keys())
    pastml_outputs.index.names = ('pos', 'ref', 'alt')
    variant_kba.close()
    # rechunk variant_kba
    variant_kba = KeyedByteArray(f'{output_path}/by_variant.kba', mode='r')
    variant_kba.rechunk(f'{output_path}/by_node.kba', jobs=50)
    variant_kba.close()
    pastml_outputs.to_csv(f'{output_path}/pastml_outputs.csv')
    return pastml_outputs


def writeEventTransitions(ancestor_calls, output_path, tree_obj):
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
    # reopen and rechunk to complete vcb file
    by_node_kba = KeyedByteArray(
        f'{output_path}/by_node.kba', mode='r')
    by_node_kba.rechunk(f'{output_path}/by_variant.kba')
    by_node_kba.close()
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


def writeEventCalls(acr_path, var_path, output_path):
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
        f'{acr_path}/pastml_outputs.csv', index_col=0)
    # initialize VCB at output path
    os.makedirs(output_path, exist_ok=True)
    ttr, tta = writeEventTransitions(  # write event transitions to keyed byte array
        ancestor_calls, output_path, tree_obj)
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
    acr_output_df.to_csv(f'{output_path}/event_call_outputs.csv')
    shutil.copy(f'{acr_path}/tree.nwk', f'{output_path}/tree.nwk')
    return acr_output_df