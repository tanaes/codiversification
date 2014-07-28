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
 make_option('-r','--results_dir',\
            help='corrected results directory [REQUIRED]'),\
 make_option('-o','--output_dir',\
            help='path to output directory [REQUIRED]'),\
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



def load_de_numericized_newick_tree(tree_in,before="'",after="'",root=False):
    from cogent.core.tree import PhyloNode
    from cogent import LoadTree
    import os.path
    
    if os.path.isfile(tree_in):
        tree = LoadTree(tree_in)
    else:
        tree = LoadTree(treestring=tree_in)
    terminals = tree.getTipNames()
    rename_dict = {}
    for tip in terminals:
        rename_dict[tip] = before + str(tip) + after
    tree.reassignNames(rename_dict)
    if root:
        tree = tree.rootAtMidpoint()
    treestring = tree.getNewick(with_distances=True)
    
    return treestring

#adding significance values to nodes
#two ways to do it

#ONE: just replace confidences with pvals

#following method taken directly from summarize_results.py, modified to return OTU basename
def read_results_file(results_fp, hack=True):
    #This reads in the results file and sticks it in dicts in a list
    results_file = open(results_fp,'Ur')
    
    results_list = []
    results_keys = results_file.next().strip().split('\t')
    
    pOTU = os.path.basename(results_fp).split('_')[0]
    
    if hack:
        results_keys.append('h_span')
    
    for line in results_file:
        entry_dict = dict(zip(results_keys,line.strip().split('\t')))
        entry_dict['pOTU'] = pOTU
        results_list.append(entry_dict)
    
    results_file.close()
    
    return pOTU, results_list
    

    

def annotate_cOTU_tree(cOTU_tree_string,results_list):
    from Bio import Phylo
    from StringIO import StringIO
    
    tree = Phylo.read(StringIO(cOTU_tree_string),'newick',rooted=True)
    
    for node_dict in results_list:
        node_tree = Phylo.read(StringIO(load_de_numericized_newick_tree(node_dict['s_nodes'],before="cOTU_",after="")),'newick',rooted=True)
        
        ###debug###
        #print node_tree
        
        node_ref = []
        for terminal in node_tree.get_terminals():
            node_ref.append({"name": terminal.name})
        
        node = tree.common_ancestor(node_ref)
        
        node.confidence = float(node_dict['fdr_p'])
        
        #print node_dict['fdr_p']
    
    out = StringIO()
    
    Phylo.write(tree,out,'newick')
    
    return out.getvalue()


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    pOTU_dir = opts.pOTU_dir
    results_dir = opts.results_dir
    output_dir = opts.output_dir    
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
    
    
    
    #read in cOTU tree
    for file in os.listdir(results_dir):
        if file.endswith("results.txt"):
            #print results_dir + '/' + file
            pOTU, results_list = read_results_file(results_dir + '/' + file)
        else:
            continue
        
        cOTU_tree_fp = pOTU_dir + '/' + pOTU + '_seqs_rep_set.tre'
        
        cOTU_tree_string = load_de_numericized_newick_tree(cOTU_tree_fp,before,after,root=True)
        
        ###debug###
        #print cOTU_tree_string
        
        annotated_tree_string = annotate_cOTU_tree(cOTU_tree_string,results_list)
        
        outfile = open(output_dir + '/' + pOTU + '_with_fdr_vals.tre','w')
        outfile.write(annotated_tree_string)
        outfile.close()
        
        #read in cOTU corrected results file
        #grab each node in turn from results file; read in treestring and fdr pval
        #
    
    #TWO: add additional pvals
    
    



if __name__ == "__main__":
    main()
