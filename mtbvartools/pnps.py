import numpy as np
from Bio.codonalign.codonseq import _count_site_NG86
from Bio.Data import CodonTable
from Bio import SeqIO

def getCodonMask(global_mask, start, end):
    # expand mask to all codons and apply to same shape as global mask
    if (end - start) % 3 != 0:
        return ValueError('Input region should be divisible by three because codons.')
    codon_rows = global_mask[start:end].reshape(-1, 3).copy()
    codon_mask = np.any(codon_rows, axis=1)
    codon_rows[codon_mask, :] = True
    output_mask = np.ones(global_mask.shape, dtype=bool)
    output_mask[start:end] = codon_rows.reshape(-1)
    return output_mask

def getSynNSynCounts(fasta_path, start, end, codon_mask, k=1, is_rev=False):
    # get target fasta sequence
    fasta_seq = SeqIO.parse(
        open(fasta_path), 'fasta').__next__().seq
    # get target gene seq
    gene_seq = fasta_seq[start:end]  # zero indexed, slice coordinates
    gene_mask = codon_mask[start:end]
    if is_rev:
        gene_seq = gene_seq.reverse_complement()
        gene_mask = gene_mask[::-1]
    sum_syn, sum_ns = 0, 0
    for i in np.arange(0, end - start, 3)[1:-1]:
        if all(gene_mask[i:i+3]):
            if any(gene_mask[i:i+3]) != all(gene_mask[i:i+3]):
                raise ValueError(f'codon mask is not properly formed for {start}:{end}, index {i}')
            continue
        # count syn and ns sites omitting start and end
        codon = gene_seq[i:i+3]
        if codon in ['TAG', 'TGA', 'TAA']:
            print(f'WARNING: prestop found between {start} and {end}')
            continue
        syn_sites, ns_sites = _count_site_NG86(
            [codon],
            CodonTable.standard_dna_table, k=k)
        sum_syn += syn_sites
        sum_ns += ns_sites
    return sum_syn, sum_ns