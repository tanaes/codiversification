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
import os
import sys
from cospeciation import recursive_hommola, make_dists_and_tree, reconcile_hosts_symbionts
from biom import load_table
from qiime.util import make_option
from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup
from qiime.parse import parse_qiime_parameters, parse_taxonomy
import os.path
import os
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent import LoadTree, LoadSeqs, DNA

# the following 3 imports added 2013-07-09 by Aaron Behr
# for method cogent_dist_to_qiime_dist
from StringIO import StringIO
from qiime.parse import parse_distmat

# parallel_test_cospeciation.py
#qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info = {}
script_info['brief_description'] = """
A script to test for cospeciation between OTUs and their host organisms"""
script_info['script_description'] = """
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
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""", """
The following steps are performed by the command below:

1. Perform Hommola et al significance test;

2. Output results to ./output
""", """
test_cospeciation_only_nate_draft_2.py -i Example_output/cOTUs_test 
-p Example_input/otu_table_HostSpecies_rarified_filtered.txt 
-a Example_input/host_tree.tre -o Example_output/hommola_test 
-T hommola -t Example_input/taxonomy.txt -m Example_input/sample_map.txt 
-c HostSpecies"""))
script_info['output_description'] = """
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
script_info['required_options'] = [
    make_option('-i', '--cotu_table_fp',
                help='the input OTU table file, or directory of files [REQUIRED]'),
    make_option('-p', '--potu_table_fp',
                help='the input pOTU table file [REQUIRED]'),
    make_option('-o', '--output_dir',
                help='path to output directory [REQUIRED]'),
    make_option('-T', '--test',
                help='desired test (unifrac, hommola, hommola_recursive) [REQUIRED]'),
    make_option('-t', '--taxonomy_fp',
                help='parent OTU taxonomy filepath [REQUIRED]'),
    make_option('-m', '--mapping_fp',
                help='path to the metadata mapping file [REQUIRED]')
]
script_info['optional_options'] = [
    make_option('--host_tree_fp',
              help='a newick-formatted tree with samples as tips [This, a host alignment, or a host distance matrix is required]'),
    make_option('--host_align_fp',
              help='A FASTA sequence alignment of the hosts [This, a host tree, or a host distance matrix is required]'),
    make_option('--host_dist_fp',
              help='A distance matrix specifying some pairwise distance -- not necessarily genetic -- between hosts [This, a host tree, or a host alignment is required]'),
    make_option('-c', '--mapping_category',
                help='map category for which to pool samples' +
                '[default: %default]',
                default='SampleID'),
    make_option('-s', '--significance_level', type='float',
                help='Desired level of significance for permutation test ' +
                '[default: %default]',
                default=0.05),
    make_option('-n', '--permutations', type='int',
                help='Desired level of significance for permutation test ' +
                '[default: %default]',
                default=1000),

    make_option('--force', action='store_true',
                dest='force', help='Force overwrite of existing output directory' +
                ' (note: existing files in output_dir will not be removed)' +
                ' [default: %default]'), options_lookup['jobs_to_start_workflow']]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    potu_table_fp = opts.potu_table_fp
    subcluster_dir = opts.cotu_table_fp
    mapping_fp = opts.mapping_fp
    mapping_category = opts.mapping_category
    output_dir = opts.output_dir
    significance_level = float(opts.significance_level)
    test = opts.test
    permutations = int(opts.permutations)
    taxonomy_fp = opts.taxonomy_fp

    if(opts.host_tree_fp):
        host_fp = opts.host_tree_fp
        host_input_type = "tree"
    elif(opts.host_align_fp):
        host_fp = opts.host_align_fp
        host_input_type = "alignment"
    elif(opts.host_dist_fp):
        host_fp = opts.host_dist_fp
        host_input_type = "distances"
    else:
        sys.exit("Must specify path to a host tree (--host_tree_fp), host alignment (--host_align_fp), or host distance matrix (--host_dist_fp).")

    # Convert inputs to absolute paths
    output_dir = os.path.abspath(output_dir)
    host_fp = os.path.abspath(host_tree_fp)
    mapping_fp = os.path.abspath(mapping_fp)
    potu_table_fp = os.path.abspath(potu_table_fp)
    subcluster_dir = os.path.abspath(subcluster_dir)

    # Check Host Tree
    try:
        with open(host_fp) as f:
            pass

    except IOError as e:
        print 'Host Data could not be opened! Are you sure it is located at ' + host_fp + '  ?'
        exit(1)

    # Check pOTU table
    try:
        with open(potu_table_fp) as f:
            pass

    except IOError as e:
        print 'parent OTU table could not be opened! Are you sure it is located at ' + potu_table_fp + '  ?'
        exit(1)

    try:
        os.makedirs(output_dir)

    except OSError:
        if force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose " +\
                "a different directory, or force overwrite with -f."
            exit(1)

    # get sample names present in potu table
    # sample_names, taxon_names, data, lineages
    potu_table = load_table(potu_table_fp)
    sample_names = potu_table.ids()
    potu_names = potu_table.ids(axis="observation")
    lineages = [lm["taxonomy"] for lm in potu_table.metadata(axis="observation")]
    
    # Process host input (tree/alignment/matrix) and take subtree of host
    # supertree
    host_tree, host_dist = make_dists_and_tree(sample_names, host_fp, host_input_type)

    # At this point, the host tree and host dist matrix have the intersect of
    # the samples in the pOTU table and the input host tree/dm.

    summary_file = open(os.path.join(output_dir,'cospeciation_results_summary.txt'), 'w')
    summary_file.write("sig_nodes\tnum_nodes\tfile\n")

    # Load taxonomic assignments for the pOTUs
    otu_to_taxonomy = parse_taxonomy(open(taxonomy_fp, 'Ur'))

    # test that you have a directory, otherwise exit.
    if os.path.isdir(subcluster_dir):
        os.chdir(subcluster_dir)
        print os.getcwd()
        # run test on cOTU tables in directory.
        # use pOTU table to choose which cOTUs to use.
        for potu in potu_names:
            # ignore comment lines

            print "Analyzing pOTU # %s" % potu

            cotu_table_fp = os.path.join(subcluster_dir,potu,'otu_table.biom')
            # Read in cOTU file
            cotu_table = load_table(cotu_table_fp)

            # Reconcile hosts in host DM and cOTU table
            cotu_table_filtered, host_dist_filtered = reconcile_hosts_symbionts(
                cotu_table, host_dist)

            # Read in reconciled cOTU table
            sample_names_filtered = cotu_table_filtered.ids()
            cotu_names_filtered = cotu_table_filtered.ids(axis="observation")

            # exit loop if less than three hosts or cOTUs
            if len(sample_names_filtered) < 3 or len(cotu_names_filtered) < 3:
                print "Less than 3 hosts or cOTUs in cOTU table!"
                continue

            # Import, filter, and root cOTU tree
            cotu_tree_fp = os.path.join(subcluster_dir,potu,"rep_set.tre")
            cotu_tree_file = open(cotu_tree_fp, 'r')
            cotu_tree_unrooted = DndParser(cotu_tree_file, PhyloNode)
            cotu_tree_file.close()
            cotu_subtree_unrooted = cotu_tree_unrooted.getSubTree(cotu_names_filtered)
            # root at midpoint
            # Consider alternate step to go through and find closest DB seq
            # to root?
            cotu_subtree = cotu_subtree_unrooted.rootAtMidpoint()

            # filter host tree
            host_subtree = host_tree.getSubTree(sample_names_filtered)
            
            align_folder = glob.glob(os.path.join(subcluster_dir,potu,'*aligned_seqs'))[0]
            # Load up and filter cOTU sequences
            aligned_otu_seqs = LoadSeqs(
                os.path.join(align_folder,'seqs_rep_set_aligned.fasta'), moltype=DNA, label_to_name=lambda x: x.split()[0])
            cotu_seqs_filtered = aligned_otu_seqs.takeSeqs(list(cotu_names_filtered))

            result = False

            # Run recursive test on this pOTU:
            # DEBUG:
            # print 'in run_test_cospeciation'

            # get number of hosts and cOTUs
            htips = len(host_subtree.getTipNames())
            stips = len(cotu_subtree.getTipNames())

            # if test == 'unifrac':
            #    print 'calling unifrac test'
            #    results_dict, acc_dict = unifrac_recursive_test(host_subtree, cotu_subtree, sample_names_filtered,
            #                                                   cotu_names_filtered, data, permutations)

            if test == 'hommola_recursive':

                # run recursive hommola test
                results_dict, acc_dict = recursive_hommola(cotu_seqs_filtered, host_subtree, host_dist_filtered, cotu_subtree, cotu_table_filtered, permutations, recurse=True)


            if test == 'hommola':

                # run recursive hommola test
                results_dict, acc_dict = recursive_hommola(cotu_seqs_filtered, host_subtree, host_dist_filtered, cotu_subtree, cotu_table_filtered, permutations, recurse=False)


            sig_nodes = 0

            # Count number of significant nodes
            for pval in results_dict['p_vals']:
                if pval < significance_level:
                    sig_nodes += 1

            num_nodes = write_results(
                results_dict, acc_dict, output_dir, potu, test, host_tree)
            result = True

            if result:
                outline = "{0}\t{1}\t{2}\t{3}".format(
                    sig_nodes, num_nodes, potu, otu_to_taxonomy[potu]) + "\n"
            else:
                outline = "ERROR\t\t" + file + "\n"
            print outline
            summary_file.write(outline)

    else:
        print 'Not a directory.'

    summary_file.close()
if __name__ == "__main__":
    main()
