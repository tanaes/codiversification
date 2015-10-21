#!/usr/bin/env python

import sys
import pandas as pd
from biom import Table

#Method to parse taxonomy summaries from PhyloSift

def parse_ps_taxa_summary(taxa_f, ok_taxa = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):

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

    return pd.DataFrame.from_dict(taxa_dict, orient='index').reindex(columns = ok_taxa)

def parse_lencovgc(lencovgc_f):
    
    lencovgc_df = pd.read_table(lencovgc_f, sep = '\t', header = None, names = list(['file','seq','len','cov','gc']), index_col = False)
    
    lencovgc_df['lencovgc_name'] = lencovgc_df.apply(lambda x:'{0}_{1}_{2:.4f}_{3:.4f}'.format(x['seq'],x['len'],x['cov'],x['gc']),axis=1)
    
    lencovgc_df['split_name'] = lencovgc_df.apply(lambda x: x['lencovgc_name'].split('.')[0],axis=1)
 
    return lencovgc_df


def initialize_ps_biom(observ_ids, sam_ids, **kwargs):
    
    data = np.zeros((len(observ_ids),len(sam_ids)))

    return Table(data, observ_ids, sam_ids, **kwargs)


def initialize_ps_df(observ_ids, sam_ids):
    
    data = np.zeros((len(observ_ids),len(sam_ids)))

    df = pd.DataFrame(data, index=observ_ids, columns=sam_ids)

    return df


def read_fasta(f):
    """from http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python"""

    name, seq = None, []
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def tax_df_to_tax(df, obs_name_col, obs_names, ok_taxa = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
    idx_df = df.set_index(obs_name_col)

    return [{'taxonomy': [x for x in idx_df.loc[y,ok_taxa]]} for y in obs_names]

def pd_df_to_biom(df, metadata_df=None, metadata_func=None, de_na=False, **kwargs):
    observ_ids = list(df.index)
    sample_ids = list(df.columns)

    if de_na:
        df = df.fillna(0)

    data = df.as_matrix()

    observation_metadata = metadata_func(metadata_df, obs_names=observ_ids, **kwargs)

    return Table(data, observ_ids, sample_ids, observation_metadata=observation_metadata)


def write_ps_cOTUs(potu_biom, cOTU_bioms, seq_info_df, aln_dict, output_dir, force=False):
    # phylosift_cOTUs/pOTU/phylosift_aligned_seqs/seqs_rep_set_aligned.fasta
    # phylosift_cOTUs/pOTU/otu_table.biom
    # phylosift_cOTUs/pOTU/rep_set.tre
    # phylosift_cOTUs/phylosift_potu_table.biom
    # phylosift_cOTUs/phylosift_seq_info.txt

    with open(os.path.join(output_dir,'phylosift_potu_table.biom'),'w') as potu_f:
        potu_f.write(potu_biom.to_json('phylosift_to_codiv.py'))

    with open(os.path.join(output_dir,'phylosift_seq_info.txt'),'w') as seq_info_f:
        seq_info_df.set_index('newname').to_csv(seq_info_f, sep='\t')

    for marker in cOTU_bioms:

        if cOTU_bioms[marker].is_empty():
            continue

        # write cOTU biom table:

        try:
            os.makedirs(os.path.join(output_dir,marker))
        except OSError:
            if not force:
                raise

        with open(os.path.join(output_dir,marker,'otu_table.biom'), 'w') as otu_f:
            otu_f.write(cOTU_bioms[marker].to_json('phylosift_to_codiv.py'))


        # write alignment:
        try:
            os.makedirs(os.path.join(output_dir,marker,'phylosift_aligned_seqs'))
        except OSError:
            if not force:
                raise

        with open(os.path.join(output_dir,
                               marker,
                               'phylosift_aligned_seqs',
                               'seqs_rep_set_aligned.fasta'), 'w') as aln_f:
            for header, seq in aln_dict[marker]:
                aln_f.write('>%s\n' % header)
                aln_f.write('%s\n\n' % seq)

    return


"""
for path in ./*; do
    [ -d "${path}" ] || continue # if not a directory, skip
    d_name="$(basename "${path}")"
    echo ${d_name}
    fasttree -nt ${d_name}/phylosift_aligned_seqs/seqs_rep_set_aligned.fasta > ${d_name}/rep_set.tre
done

"""


def main():

    
    try:
        os.makedirs(output_dir)
    except OSError:
        if not force:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            raise OSError("Output directory already exists. Please choose "
                "a different directory, or force overwrite with -f.")
    # parse file of sample information

    sample_info = pd.read_table(sample_info_fp, index_col=False)

    sample_ids = [x for x in sample_info.Sample]

    # eat the list of markers (these are pOTUs)

    with open(markers_fp, 'r') as f:
        markers = [x.strip() for x in f.readlines()]

    if marker_name_fp:
        marker_tax_df = pd.read_table(marker_name_fp, header=None, names=['gene'], index_col=0)
        marker_name_tax = [{'taxonomy': [x[1][0]]} for x in marker_taxa_df.reindex(markers).iterrows()]
    else:
        marker_name_tax = None


    # initialize a pOTU biom table -- samples by markers

    pOTU_df = initialize_ps_df(markers, sample_ids)


    # initialize a dict of cOTU biom tables -- sample by contig

    cOTU_dfs = {}
    aln_dict = {}

    for marker in markers:
        cOTU_dfs[marker] = initialize_ps_df([], sample_ids)
        aln_dict[marker] = []

    seq_info_df = pd.DataFrame()

    for row in sample_info.iterrows():

        sample = row[1].SampleID
        sample_dir = row[1].ps_dir
        lencovgc_fp = row[1].lencovgc_file

        # if specified, eat a len/cov/gc file and add to observation metaxata
        with open(lencovgc_fp, 'r') as lencovgc_f:
            lencovgc_df = parse_lencovgc(lencovgc_f)

            lencovgc_df_split = lencovgc_df.set_index('split_name')

        sample_seq_info_df = pd.DataFrame()

        for marker in markers:

            # for each sample, for each marker, eat the alignment

            try:
                with open(os.path.join(sample_dir, 'alignDir', marker + '.codon.updated.1.fasta')) as aln_f:
                    
                    i = 0
                    
                    for x,y in read_fasta(aln_f):
                        # set some variable names for the fasta header and seq
                        oldname = x[1:]
                        seq = y

                        # rename the sequence by the sample, marker, and i
                        newname = '{0}.{1}.{2}'.format(sample, marker, i)

                        # split the name on a period to remove PS cruft
                        split_name = oldname.split('.')[0]

                        # add all of this info to the provenance tracking df
                        sample_seq_info_df = sample_seq_info_df.append({'sample': sample,
                                                    'marker': marker,
                                                    'oldname': oldname,
                                                    'newname': newname,
                                                    'split_name': split_name}, ignore_index=True)

                        # add the fasta tuple to the list of seqs for this marker
                        aln_dict[marker] += [(newname, seq)]

                        # This is ugly, but:
                        # pull the coverage info out of the df and add it to a
                        # df for that marker
                        cOTU_dfs[marker].loc[newname,sample] = lencovgc_df_split.loc[split_name,'cov'] 

                        i += 1

                    # after adding all the cOTUs, append the taxonomy columns to seq name df


                    # and summarize the number of seqs in a pOTU table
                    pOTU_df.loc[marker,sample] = i
            except IOError as (errno, strerror):
                print "I/O error({0}): {1}".format(errno, strerror)


        # eat a taxon table
        with open(os.path.join(sample_dir, 'sequence_taxa_summary.txt'), 'r') as taxa_f:
            taxa_df = parse_ps_taxa_summary(taxa_f)
            taxa_df = taxa_df.reset_index()
            taxa_df['split_name'] = taxa_df.apply(lambda x: x['index'].split('.')[0],axis=1)

        # add the taxonomy columns to the per-sequence mapping file
        seq_info_df = pd.concat([seq_info_df, 
                                 pd.merge(sample_seq_info_df,
                                          taxa_df,
                                          how='left',
                                          on=['split_name','split_name'],
                                          sort=False)],
                                ignore_index=True)


    # now send these to biom objects
    cOTU_bioms = {}
    for marker in cOTU_dfs:
        cOTU_bioms[marker] = pd_df_to_biom(cOTU_dfs[marker],
                                  metadata_df=seq_info_df,
                                  metadata_func=tax_df_to_tax,
                                  de_na=True,
                                  obs_name_col='newname')

        cotu_biom.to_json('phylosift_to_codiv.py')

    potu_biom = pd_df_to_biom(pOTU_df,
                              metadata_df = marker_tax_df.reset_index(),
                              metadata_func=tax_df_to_tax,
                              de_na=True,
                              obs_name_col='index',
                              ok_taxa=['gene'])


    # output in subclustered directory structure

    force = True

    write_ps_cOTUs(potu_biom, cOTU_bioms, seq_info_df, aln_dict, output_dir, force=True)

