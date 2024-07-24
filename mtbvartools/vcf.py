import numpy as np
import pandas as pd
import pysam

def bedToMask(path, genome_length):
    mask = np.zeros(genome_length).astype(bool)
    bed = pd.read_csv(
        path, delimiter='\t', index_col=0, header=None)
    for _, (l, r) in bed.iterrows():
        mask[l:r] = True
    return mask


def SNPEffToSeries(
        record,
        annotation_fields=['Annotation', 'Gene_ID', 'HGVS.c', 'cDNA.pos', 'HGVS.p', 'CDS.pos']):
    # return empty series if not called
    if 'ANN' not in record.info.keys():
        return pd.Series()
    # column definitions
    snpeff_dict = {
        'Allele': 0,
        'Annotation': 1,
        'Annotation_Impact': 2,
        'Gene_Name': 3,
        'Gene_ID': 4,
        'Feature_Type': 5,
        'Feature_ID': 6,
        'Transcript_BioType': 7,
        'Rank': 8,
        'HGVS.c': 9,
        'HGVS.p': 10,
        'cDNA.pos': 11,
        'cDNA.length': 12,
        'CDS.pos': 13,
        'CDS.length': 14,
        'AA.pos': 15,
        'AA.length': 16,
        'Distance': 17,
        'ERRORS / WARNINGS / INFO': 18}
    annotation_data = [[] for f in annotation_fields]
    # parse annotation string
    for annotation_string in record.info['ANN']:
        split_annotations = annotation_string.split('|')
        for i, field in enumerate(annotation_fields):
            try:
                annotation_data[i].append(split_annotations[snpeff_dict[field]])
            except IndexError:  # empty field if value out of range
                annotation_data[i].append('')
    # generate dataframe
    annotation_df = pd.DataFrame(
        index=pd.Index(annotation_fields, name='anno_field'),
        columns=np.arange(len(record.info['ANN'])) + 1,
        data=annotation_data,)
    # melt to single value
    annotation_df = annotation_df.reset_index().melt(
        id_vars=['anno_field'], var_name='anno_num', value_name='annotation')
    # convert into a series for output
    output_series = annotation_df.set_index(['anno_field', 'anno_num']).loc[:, 'annotation']
    return output_series


def SIFTtoSeries(
        record,
        annotation_fields=['GeneId', 'SIFT_num_seqs', 'SIFT_score', 'SIFT_prediction']):
    # return empty series if not called
    if 'SIFTINFO' not in record.info.keys():
        return pd.Series()
    # column definitions
    sift_dict = {
        'Allele': 0,
        'Transcript': 1,
        'GeneId': 2,
        'GeneName': 3,
        'Region': 4,
        'VariantType': 5,
        'Ref_Amino_Acid/Alt_AminoAcid': 6,
        'Amino_position': 7,
        'SIFT_score': 8,
        'SIFT_median': 9,
        'SIFT_num_seqs': 10,
        'Allele_Type': 11,
        'SIFT_prediction': 12}
    annotation_data = [[] for f in annotation_fields]
    # parse annotation string
    for annotation_string in record.info['SIFTINFO']:
        split_annotations = annotation_string.split('|')
        for i, field in enumerate(annotation_fields):
            try:
                annotation_data[i].append(split_annotations[sift_dict[field]])
            except IndexError:  # empty field if value out of range
                annotation_data[i].append('')
    # generate dataframe
    annotation_df = pd.DataFrame(
        index=pd.Index(annotation_fields, name='anno_field'),
        columns=np.arange(len(record.info['SIFTINFO'])) + 1,
        data=annotation_data,)
    # melt to single value
    annotation_df = annotation_df.reset_index().melt(
        id_vars=['anno_field'], var_name='anno_num', value_name='annotation')
    # # convert into a series for output
    output_series = annotation_df.set_index(['anno_field', 'anno_num']).loc[:, 'annotation']
    return output_series


def filterAnnotation(annotation, field, value):
    if len(annotation) == 0:
        return pd.Series()
    match_filter = annotation.loc[pd.IndexSlice[field, :]] == value
    try:
        return annotation.loc[
            pd.IndexSlice[:, match_filter.index[match_filter][0]]]
    except IndexError:
        return pd.Series()


def fetchGeneVariantAnnotations(target_vcf, geneid, start, stop):
    # start stop are 1-based, inclusive
    # preload records into memory to free up the file for other workers
    record_list = []
    with pysam.VariantFile(target_vcf) as open_vcf:
        region_string = f'{open_vcf.get_reference_name(0)}:{start}-{stop}'
        for record in open_vcf.fetch(region=region_string):
            record_list.append(record)
    if len(record_list) == 0:
        return pd.DataFrame()  # return empty dataframe if no records
    # get annotation data
    gene_annotation_data = []
    for record in record_list:
        # get annotation data
        record_info = pd.Series(
            index=['pos', 'ref', 'alt'],
            data=[record.pos, record.ref, record.alts[0]])
        snpeff_out = filterAnnotation(
            SNPEffToSeries(record), 'Gene_ID', geneid)
        sift_out = filterAnnotation(
            SIFTtoSeries(record), 'GeneId', geneid)
        gene_annotation_data.append(
            pd.concat([record_info, snpeff_out, sift_out], axis=0))
    # convert to DataFrame
    gene_annotation_data = pd.concat(
        gene_annotation_data, axis=1).T
    # return
    return gene_annotation_data