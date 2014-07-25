#!/n/sw/Python-2.6.2/bin/python
# File created on 3 October 2011.
from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2011, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Experimental"

from qiime.util import make_option
import os
import numpy
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 get_options_lookup, parse_otu_table
from qiime.parse import parse_qiime_parameters, parse_taxonomy
from Bio import Phylo


#test_cospeciation.py
qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""
Annotates cOTU tree nodes with fdr support values."""
script_info['script_description']="""
This script reads in a directory of fdr-corrected results, a directory of pOTU
preclusters, and then outputs annotated cOTU trees to an output directory.

NOTE: not a good QIIME script; this requires use of the BioPython Phylo library
for appropriate node annotation. I could modify this to use the name_nodes hack
that QIIME uses (e.g. for the jacknife support values). 

Currently, the default is to place the taxon names in quote-postrophes; but you
can optionally pass the script --before and --after values, for example to pre-
fix or suffix the numeric cOTU label.

"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""
The following steps are performed by the command below:

1. Read in each cOTU tree in the pOTU directory

2. Convert to a Bio.Phylo-compatible Newick format (de numericize taxa labels). 

3. Read in cospeciation support values for each node of cOTU tree

4. Store support value as confidence value

5. Write newick-with-confidence tree to output_dir
""","""
annotate_nodes_with_sig_values.py -p pOTU_dir -r results_dir -o output_dir """))
script_info['output_description']="""
We'll see!"""
script_info['required_options']=[\
 make_option('-c','--pOTU_dir',\
            help='pOTU directory [REQUIRED]'),\
 make_option('-r','--input_tree_dir',\
            help='corrected tree results directory [REQUIRED]'),\
 make_option('-o','--output_dir',\
            help='path to output directory [REQUIRED]'),\
 make_option('-m','--mapping_fp',\
            help='path to metadata mapping file [REQUIRED]'),\
 make_option('-l','--mapping_cat',\
            help='mapping category to add [REQUIRED]'),\
]
script_info['optional_options']=[\
 make_option('--before',\
            help='cOTU taxon label prefix'+\
            '[default: single quote]',\
            default=None),\
 make_option('--after',\
            help='cOTU taxon label suffix'+\
            '[default: single quote]',\
            default=None),\
 make_option('--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
]

script_info['version'] = __version__




def stupid_hack(summarized_otu_table):
    import re
    
    p = re.compile('Consensus Lineage')
    q = re.compile('\t\n')
    summarized_otu_table += '\n'
    summarized_otu_table = p.sub('',summarized_otu_table)
    summarized_otu_table = q.sub('\n',summarized_otu_table)
    return summarized_otu_table


def relabel_tree_by_otu_table(orig_tree_fp,otu_table_fp,mapping_fp,mapping_cat='HostSpecies',before="'",after="'"):
    from cogent import LoadTree
    import StringIO
    from qiime.summarize_otu_by_cat import summarize_by_cat
    from Bio import Phylo
    
    #read in confidence-relabeled OTU tree into Bio.Phylo object
    
    orig_tree = Phylo.read(orig_tree_fp,'newick')
    
    #get names from tree, make as keys to a dictionary
    
    name_dict = {}
    
    for terminal in orig_tree.get_terminals():
        if terminal.name[0] == "'" and terminal.name[-1] == "'":
            terminal.name = terminal.name[1:-1]
        name_dict[terminal.name] = {}
    
    mapping_file = open(mapping_fp, 'Ur')
    otu_file = open(otu_table_fp, 'Ur')
    
    #get otu tables and summarize appropriately
    summarized_otu_table = \
     summarize_by_cat(mapping_file,otu_file,mapping_cat,False)
    summarized_otu_table = stupid_hack(summarized_otu_table)
    otu_obj = StringIO.StringIO(summarized_otu_table)
    mapping_file.close()
    otu_file.close()
    
    #read otu table
    
    #parse OTU table 
    try:
        sample_names, taxon_names, data, lineages = parse_otu_table(otu_obj)
    except:
        print " not an OTU table?"
    
    #populate the terminal name dictionary with a dictionary mapping samples to counts
    for i in range(len(taxon_names)):
        
        for j in range(len(sample_names)):
            if data[i][j]:
                name_dict[before + taxon_names[i] + after][sample_names[j]] = data[i][j]
    
    #rename terminals accordingly:
    
    for terminal in orig_tree.get_terminals():
        otu = terminal.name
        new_name = otu
        for key in name_dict[otu]:
            new_name += '_' + key + '[' + str(name_dict[otu][key]) + ']'
        
        terminal.name = new_name
        
    return orig_tree
    
    

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    """
    orig_tree_fp,otu_table_fp,mapping_fp,mapping_cat,before,after
    """
    
    
    pOTU_dir = opts.pOTU_dir
    input_tree_dir = opts.input_tree_dir
    output_dir = opts.output_dir  
    mapping_fp = opts.mapping_fp
    mapping_cat = opts.mapping_cat
    before = opts.before
    after = opts.after
    force = opts.force
    
    if not (before or after):
        before = "'"
        after = "'"
    else:
        if not before:
            before = ''
        if not after:
            after = ''
    
    
    try:
        os.makedirs(output_dir)
    except OSError:
        if force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)
    
    
    
    """
    orig_tree_fp = '/n/home02/jsanders/labdir/AmpliconNoise/95_preclust/otus/fasttree_99_copy/957_fdr_annotated_cOTU.tre'
    otu_table_fp = '/n/home02/jsanders/labdir/AmpliconNoise/95_preclust/otus/fasttree_99_copy/cOTUs/957_seqs_otu_table.txt'
    mapping_fp = '/n/home02/jsanders/labdir/AmpliconNoise/cephalotes_denoise_map.txt'
    mapping_cat = 'HostSpecies'
    before = 'cOTU_'
    after = ''
    
    tree = relabel_tree_by_otu_table(orig_tree_fp,otu_table_fp,mapping_fp,mapping_cat,before,after)
    
    
    
    otu_table_dir = '/n/home02/jsanders/labdir/AmpliconNoise/95_preclust/otus/fasttree_99_copy/cOTUs'
    
    input_tree_dir = '/n/home02/jsanders/labdir/AmpliconNoise/95_preclust/otus/fasttree_99_copy/hommola_corrected_noNA_trees'
    
    output_dir = '/n/home02/jsanders/labdir/AmpliconNoise/95_preclust/otus/fasttree_99_copy/hommola_corrected_noNA_trees_hostlabeled'
    
    mapping_fp = '/n/home02/jsanders/labdir/AmpliconNoise/cephalotes_denoise_map.txt'
    mapping_cat = 'HostSpecies'
    before = 'cOTU_'
    after = ''
    
    """
    
    for file in os.listdir(input_tree_dir):
        if file.endswith(".tre"):
            #print results_dir + '/' + file
            orig_tree_fp = input_tree_dir + '/' + file
            pOTU = file.split('_')[0]
        else:
            continue
        
        otu_table_fp = pOTU_dir + '/' + pOTU + '_seqs_otu_table.txt'
        
        tree = relabel_tree_by_otu_table(orig_tree_fp,otu_table_fp,mapping_fp,mapping_cat,before,after)
        
        Phylo.write(tree,output_dir + '/' + pOTU + '_fdr_annotated_cOTU_withHosts.tre','newick')
        
        
        

if __name__ == "__main__":
    main()
