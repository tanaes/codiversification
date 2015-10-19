#!/usr/bin/env python

import sys
import pandas as pd

#Method to parse taxonomy summaries from PhyloSift

def parse_ps_taxa_summary(taxa_f, ok_taxa = {'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}):

    """
    It's a list of inferred taxonomic ranks, ending with the uncertain placements. 
    Read in line by line, 
    if [rank] in [good ranks] and [probability mass] > 0.5
        make that the taxon for that rank. 

    return pandas DataFrame object
    """

    taxa_dict = {}

    taxa_lines = taxa_f.readlines()

    firstLine = taxa_lines.pop(0)

    for line in taxa_lines:
        taxon_line = line.strip().split('\t')

        if taxon_line[0] not in taxa_dict:
            taxa_dict[taxon_line[0]] = {'superkingdom': None, 'phylum': None, 'class': None, 'order': None, 'family': None, 'genus': None, 'species': None}

        if taxon_line[3] in ok_taxa and float(taxon_line[5]) > 0.5:
            taxa_dict[taxon_line[0]][taxon_line[3]] = taxon_line[4]

    return pd.DataFrame.from_dict(taxa_dict, orient='index')

def parse_lencovgc(lencovgc_f):
    
    lencovgc_df = pd.read_table(lencovgc_f, sep = '\t', header = None, names = list(['file','seq','len','cov','gc']), index_col = False)
    
    df_renamed = lencovgc_df.rename(index = {x[0]: '{0}_{1}_{2:.4f}_{3:.4f}'.format(x[1]['seq'],
                                                                 x[1]['len'],
                                                                 x[1]['cov'],
                                                                 x[1]['gc'])
                                for x in lencovgc_df.iterrows()})

    return df_renamed

def main():

    
    # eat a taxon table

    # eat a list of sample directories (these are sample ids)

    # if specified, eat a len/cov/gc file and add to observation metaxata

    # eat the list of markers (these are pOTUs)

    # for each sample, for each marker, eat the alignment

    # for each alignment, replace name with MARKER.SAMPLE.CONTIG and add cat aln

    # add number of seqs in alignment to pOTU biom

    # add cov (if read) or 1 to cOTU biom

    # add obs metadata to cOTU biom

    return

"""
#Sequence_ID    Hit_Coordinates NCBI_Taxon_ID   Taxon_Rank  Taxon_Name  Cumulative_Probability_Mass Markers_Hit
scaffold_124805_325_3.3754_0.5477   1.270   1224    phylum  PROTEOBACTERIA  1   concat
"""