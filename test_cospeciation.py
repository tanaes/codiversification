#!/n/sw/python-2.7.1/bin/python
# File created June 2014.
from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2014, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Experimental"

from qiime.util import make_option
import os, sys
from cospeciation import *

from qiime.util import make_option
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 get_options_lookup
from qiime.parse import parse_qiime_parameters, parse_taxonomy
import os.path
import os

# the following 3 imports added 2013-07-09 by Aaron Behr
# for method cogent_dist_to_qiime_dist 
from StringIO import StringIO
from qiime.parse import parse_distmat

#parallel_test_cospeciation.py
#qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""
A script to test for cospeciation between OTUs and their host organisms"""
script_info['script_description']="""
This script will perform recursive unifrac or distance-based cospeciation 
analysis on 'child' OTUs that have previously been subclustered from wider, 
'parent' OTUs. The relationships among the constituents of the 'parent' OTU -- 
i.e., the relationships among 'child' OTUs -- will be compared to the relation-
ships among their hosts, represented as a phylogenetic tree. The script will 
indicate whether each node of the child OTU tree for each parent OTU is signifi-
cantly correlated with the host tree, and more broadly, which of the parent OTUs
 have appeared to codiversify with the hosts. 
 
As input, you need to provide a folder of subclustered OTUs, i.e. the results of
 otu_subcluster.py, a QIIME taxonomy file indicating the taxonomic affiliation 
of each of those parent OTUs, and a Newick-formatted tree of relationships among
 your host organisms. Alternatively, you may provide a DNA alignment or distance
 matrix in place of a host tree.
 
For the host tree, branch length matters, and the tip labels of the tree must be
 identical to the sample names in your OTU tables. Optionally, you may provide a
 metadata sample mapping file to relate sample names to the tips of the host 
tree.

Note that if you used a summarized OTU table for the OTU subclustering step 
(i.e. you combined samples by species or genus prior to subclustering) you will 
need to pass the appropriate mapping category linking SampleID to the summarized
 names represented in your subclustered pOTU tables and host tree.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""
The following steps are performed by the command below:

1. Perform Hommola et al significance test;

2. Output results to ./output
""","""
test_cospeciation_only_nate_draft_2.py -i Example_output/cOTUs_test 
-p Example_input/otu_table_HostSpecies_rarified_filtered.txt 
-a Example_input/host_tree.tre -o Example_output/hommola_test 
-T hommola -t Example_input/taxonomy.txt -m Example_input/sample_map.txt 
-c HostSpecies"""))
script_info['output_description']="""
This script results in a folder containing a series of files describing the 
outcome of the cospeciation analysis for each parent OTU, plus a summary file 
indicating, for each pOTU, the number of nodes tested, the number significant, 
and the taxonomic assignment of that pOTU. 

Per-pOTU output files may vary by the specific test, but always include the 
following information for each tested node of the pOTU tree:
h_nodes: the hosts subtree represented by cOTUs downstream of that node in the 
symbiont tree.
s_nodes: the cOTU tree descendent from the tested node in the symbiont tree.
p_vals: p-value significance of codiversification for that node
s_tips: the number of cOTUs descendent of that node
h_tips: the number of hosts represented by the above cOTUs
h_span: the number of hosts in the minimum monophyletic clade of the host tree 
that is spanned by the above hosts. 
"""
script_info['required_options']=[\
 make_option('-i','--cotu_table_fp',\
            help='the input OTU table file, or directory of files [REQUIRED]'),\
make_option('-p','--potu_table_fp',\
            help='the input pOTU table file [REQUIRED]'),\
 make_option('-a','--host_tree_fp',\
            help='a newick-formatted tree with samples as tips OR a FASTA sequence alignment OR a distance matrix [REQUIRED]'),\
 make_option('-o','--output_dir',\
            help='path to output directory [REQUIRED]'),\
 make_option('-T','--test',\
            help='desired test (unifrac, hommola, hommola_recursive) [REQUIRED]'),\
 make_option('-t','--taxonomy_fp',\
            help='parent OTU taxonomy filepath [REQUIRED]'),\
 make_option('-m','--mapping_fp',
    help='path to the metadata mapping file [REQUIRED]'),\
]
script_info['optional_options']=[
 make_option('-c','--mapping_category',\
            help='map category for which to pool samples'+\
            '[default: %default]',\
            default='SampleID'),\
 make_option('-s','--significance_level',type='float',\
            help='Desired level of significance for permutation test '+\
            '[default: %default]',\
            default=0.05),\
 make_option('-n','--permutations',type='int',\
            help='Desired level of significance for permutation test '+\
            '[default: %default]',\
            default=1000),\
            
 make_option('--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),options_lookup['jobs_to_start_workflow']]

script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    potu_table_fp = opts.potu_table_fp
    cotu_table_fp = opts.cotu_table_fp
    host_tree_fp = opts.host_tree_fp
    mapping_fp = opts.mapping_fp
    mapping_category = opts.mapping_category
    output_dir = opts.output_dir
    significance_level = float(opts.significance_level)
    test = opts.test
    permutations = int(opts.permutations)
    taxonomy_fp = opts.taxonomy_fp
    test_cospeciation(potu_table_fp,cotu_table_fp,host_tree_fp,mapping_fp,mapping_category,output_dir,significance_level,test,permutations,taxonomy_fp,opts.force)

if __name__ == "__main__":
    main()

