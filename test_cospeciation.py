#!/n/sw/python-2.7.1/bin/python
# File created on 3 October 2011.
from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2011, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Experimental"

from qiime.util import make_option
import os, sys
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 get_options_lookup, parse_otu_table
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_distmat
from cogent import LoadSeqs, DNA
from qiime.summarize_otu_by_cat import summarize_by_cat
import StringIO
import numpy
from cogent.parse.tree import DndParser
from qiime.filter import filter_samples_from_otu_table
from cogent.core.tree import PhyloNode

# the following 3 imports added 2013-07-09 by Aaron Behr
# for method cogent_dist_to_qiime_dist 
from StringIO import StringIO
from qiime.parse import parse_distmat
from cogent.util.dict2d import Dict2D

"""
from cogent.app.muscle import align_unaligned_seqs as muscle_aln
from cogent.app.fasttree import build_tree_from_alignment as fasttree_build_tree
from cogent.core.tree import PhyloNode
from cogent import LoadTree, LoadSeqs, DNA
from cogent.parse.tree import DndParser
from qiime.summarize_otu_by_cat import summarize_by_cat
import re
import numpy
import StringIO
"""


#test_cospeciation.py
qiime_config = load_qiime_config()
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
            help='desired test (unifrac or hommola) [REQUIRED]'),\
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

def hommola_cospeciation_test(host_dist, par_dist, matrix, permutations):
    """Performs the cospeciation test from Hommola et al recursively over a tree.
    
    Takes numpy matrices of jxj host distances, ixi 'parasite' (OTU) distances, 
    and a binary ixj association matrix. 
    
    test data from Hommola et al MB&E 2009: 
    
    hdist = numpy.array([[0,3,8,8,9],[3,0,7,7,8],[8,7,0,6,7],[8,7,6,0,3],[9,8,7,3,0]])
    pdist = numpy.array([[0,5,8,8,8],[5,0,7,7,7],[8,7,0,4,4],[8,7,4,0,2],[8,7,4,2,0]])
    int = numpy.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,1,1]])
    
    This is basically a direct translation from the R code, and not optimized
    in any way for Python. 
    
    """
    import cogent.maths.stats.test as stats
    from random import shuffle
    import numpy
    
    m = matrix.sum()
    
    hosts = [0]*m
    pars = [0]*m
    
    
    
    
    #Generate lists of host and symbiont edges
    s = 0
    while s < m:
        #print "s: " + str(s)
        for i in range(matrix.shape[0]):
            #print "i: " + str(i)
            for j in range(matrix.shape[1]):
                #print "j: " + str(j)
                if matrix[i,j] == 1:
                    #pause = raw_input("yes")
                    hosts[s] = j
                    pars[s] = i
                    s += 1
    
    #now we have a host and parasite list, where the index of the list
    #represents an edge connecting the host to the parasite
    
    #get a vector of pairwise distances for each interaction edge
    x = get_dist(hosts,host_dist,range(matrix.shape[1]))
    y = get_dist(pars,par_dist,range(matrix.shape[0]))
    
    
    #calculate the observed correlation coefficient for this host/symbionts
    r = stats.correlation(x, y)[0]
    
    
    #now do permutaitons. Initialize index lists of the appropriate size.
    mp = range(par_dist.shape[1])
    mh = range(host_dist.shape[1])
    below = 0
    
    
    for i in range(permutations):
        #Generate a shuffled list of indexes for each permutation. This effectively
        #randomizes which host is association with which symbiont, but maintains
        #the distribution of genetic distances.
        shuffle(mp)
        shuffle(mh)
        
        #Get pairwise distances in shuffled order
        y_p = get_dist(pars,par_dist,mp)
        x_p = get_dist(hosts,host_dist,mh)
        
        #calculate shuffled correlation. 
        #If greater than observed value, iterate counter below.
        r_p = stats.correlation(x_p,y_p)[0]
        if r_p >= r:
            below += 1
        
    
    #DEBUG!
    #print "p: " + str(below) + "/" + str(permutations)
    
    return float(below + 1)/float(permutations + 1)

def get_dist(labels,dists,index):
    """Function for picking a subset of pairwise distances from a distance matrix
    according to a set of (randomizable) indices. Derived from Hommola et al R code"""
    m = len(labels)
    vec = []
        
    for i in range(m-1):
        #print "i: " + str(i)
        k = index[labels[i]]
        for j in range(i+1,m):
            #print "j: " + str(j)
            t = index[labels[j]]
            vec.append(dists[k,t])
    
    return vec



def dist2Dict2D(dist_dict,sample_names):
    """This little routine takes distance calc outputs from Cogent and sticks 
    them into a Dict2D object."""
    
    from cogent.util.dict2d import Dict2D
    
    #print "Before dist dict:"
    #print dist_dict
    
    dist_list = []
    
    #this ugliness is necessary because Dict2D guesses wrong with symmetric data
    #structures of size 3 where you want padded diagonals. 
    if len(sample_names)==3:
        for x in range(3):
            row_list=[]
            for y in range(3):
                if x==y:
                    row_list.append(0)
                else:
                    try:
                        row_list.append(dist_dict[(sample_names[x],sample_names[y])])
                    except:
                        row_list.append(dist_dict[(sample_names[y],sample_names[x])])
            dist_list.append(row_list)
    else:       
        for item in dist_dict.iteritems():
            #DEBUG:
            #print item
            #pause = raw_input("")
            dist_list.append([item[0][0],item[0][1],item[1]])
    
    dict2d = Dict2D(dist_list,sample_names,sample_names,0,True)
    
    #print "After dist dict:"
    #print dist_dict
    #pause = raw_input("pause")
    
    return dict2d


def cogent_dist_to_qiime_dist(dist_tuple_dict):
    """
    This takes a dict with tuple keys and distance values, such as is output
    by the getDistances() method of a PhyloNode object, and converts it to a 
    QIIME-style distance matrix object: an ordered tuple with a list of samples
    in [0] and a numpy array of the distance matrix in [1].
    
    EDITED AND UPDATED 2013-07-09 Aaron Behr
    """
    from StringIO import StringIO
    from qiime.parse import parse_distmat
    from cogent.util.dict2d import Dict2D
    
    headers = []
    dist_dict = {}
    
    # loop through dist_tuple_dict, returning (k1,k2):v tuples simultaneously
    for item in dist_tuple_dict.iteritems():
        if item[0][0] not in headers: # if k1 is not in headers, add it to headers
            headers.append(item[0][0])
            dist_dict[item[0][0]] = {item[0][0]: 0.0} # null self-distance
            
        dist_dict[item[0][0]][item[0][1]] = item[1] # dist_dict[k1][k2] = v
    headers.sort()
    
    
    # Initialize dict2d, with data from dist_dict (dict of dicts).
    # Also, RowOrder and ColOrder are set to the order of the sorted headers.
    # NOTE: no longer using the fromDicts() method to pass dist_dict to dict2d
    dict2d = Dict2D(dist_dict, headers, headers)
    
    # output tab-delimited printable string of the items in dict2d including headers.    
    dist_delim = dict2d.toDelimited() 
    
    # generate and return Qiime distance matrix
    return parse_distmat(StringIO(dist_delim[1:]))



"""
dist_tuple_dict = {('SHAJ', 'SHAK'): 0.10750048520885,
 ('SHAJ', 'SHAM'): 0.10750048520885,
 ('SHAJ', 'SHOA'): 0.0147434146325,
 ('SHAJ', 'SHOG'): 0.0147434146325,
 ('SHAK', 'SHAJ'): 0.10750048520885,
 ('SHAK', 'SHAM'): 0.048024926561999998,
 ('SHAK', 'SHOA'): 0.10750048520885,
 ('SHAK', 'SHOG'): 0.10750048520885,
 ('SHAM', 'SHAJ'): 0.10750048520885,
 ('SHAM', 'SHAK'): 0.048024926561999998,
 ('SHAM', 'SHOA'): 0.10750048520885,
 ('SHAM', 'SHOG'): 0.10750048520885,
 ('SHOA', 'SHAJ'): 0.0147434146325,
 ('SHOA', 'SHAK'): 0.10750048520885,
 ('SHOA', 'SHAM'): 0.10750048520885,
 ('SHOA', 'SHOG'): 0.0,
 ('SHOG', 'SHAJ'): 0.0147434146325,
 ('SHOG', 'SHAK'): 0.10750048520885,
 ('SHOG', 'SHAM'): 0.10750048520885,
 ('SHOG', 'SHOA'): 0.0}
 
 
 qiime_distmat = (['SHOA', 'SHOG', 'SHAJ', 'SHAK', 'SHAM'],
 array([[ 0.01474341,  0.01474341,  0.        ,  0.10750049,  0.10750049],
       [ 0.        ,  0.        ,  0.01474341,  0.10750049,  0.10750049],
       [ 0.10750049,  0.10750049,  0.10750049,  0.04802493,  0.        ],
       [ 0.        ,  0.        ,  0.01474341,  0.10750049,  0.10750049],
       [ 0.10750049,  0.10750049,  0.10750049,  0.        ,  0.04802493]]))
 
"""

def recursive_hommola(aligned_otu_seqs,host_subtree,host_dist,otu_tree,sample_names,
                        taxon_names,otu_data,permutations=1000):
    """
    Applies Hommola et al test of cospeciation recursively to OTU tree.
    
    Conceptually similar to the oraganization of the recursive Unifrac method.
    
    Host distances are calculated from the provided host tree, and OTU distances
    from the MUSCLE alignment using the TN93 model of nucleotide evolution. It 
    would probably be better to do pairwise alignments and calculate distances
    that way. 
    
    It returns a dictionary of the results, and an empty accessory dict.
    """
    
    
    
    from cogent.evolve.pairwise_distance import TN93Pair
    from cogent.maths.unifrac.fast_unifrac import fast_unifrac
    from cogent.parse.tree import DndParser
    from cogent.core.tree import PhyloNode, TreeError
    from cogent.util.dict2d import Dict2D
    import numpy
    
    
    print "Performing recursive Hommola et al cospeciation test..." 
    
    #calculate pairise distances between OTUs
    
    dist_calc = TN93Pair(DNA, alignment=aligned_otu_seqs)
    dist_calc.run()
    
    otu_dists = dist_calc.getPairwiseDistances()
    
    #convert pw distances (and tree distances for hosts) to numpy arrays with same
    #column/row headings as host/OTU positions in OTU table numpy array.
    
    
    hdd = dist2Dict2D(host_dist,sample_names)
    hdn = numpy.array(hdd.toLists())
    sdd = dist2Dict2D(otu_dists,taxon_names)
    sdn = numpy.array(sdd.toLists())
    
    #convert OTU table to binary array, throwing out all OTUs below a given thresh.
    
    
    interaction = otu_data.clip(0,1)
    
    #traverse OTU tree and test each node
    
    #initialize our output lists
    s_nodes, h_nodes, p_vals, s_tips, h_tips = [], [], [], [], []
    
    
    #DEBUG!
    #print hdn
    #print sdn
    
    #iterate over the tree of child OTUs
    for node in otu_tree.traverse(self_before=True, self_after=False):
        
        #get just OTUs in this node
        otu_subset = node.getTipNames()
        s_vec = []
        h_vec = []
        h_names = []
        
        #find positional index (from OTU table) for each cOTU represented in this node:
        for i in range(len(taxon_names)):
            if taxon_names[i] in otu_subset:
                s_vec.append(i)
        #slice symbiont distance matrix down to only cOTUs in this node 
        s_slice = sdn[numpy.ix_(s_vec,s_vec)]
        
        #slice interaction matrix down to only cOTUs in this node
        i_s_slice = interaction[numpy.ix_(s_vec)]
        
        #find positional index (this time from OTU table size) for each sample in this node:
        #sum all values in column for each host, if greater than zero, add that host position to h_vec
        for j in range(i_s_slice.shape[1]):
            if i_s_slice[:,j].sum():
                h_vec.append(j)
                h_names.append(sample_names[j])
        
        
        
        #Make sure we have at least 3 hosts and symbionts represented
        if len(h_vec) > 2 and len(s_vec) > 2:
            #slice interaction matrix
            i_slice = interaction[numpy.ix_(s_vec,h_vec)]
            
            #slice host distance matrix
            h_slice = hdn[numpy.ix_(h_vec,h_vec)]
            
            #append symbiont nodes and host subtrees as tree objects
            s_nodes.append(node)
            h_nodes.append(host_subtree.getSubTree(h_names))
            
            #append number of symbionts and hosts for this node
            s_tips.append(len(s_vec))
            h_tips.append(len(h_vec))
            #calculate pemutation p value for hommola test for this node
            p = hommola_cospeciation_test(h_slice, s_slice, i_slice, permutations)
            #append to results list
            p_vals.append(p)
            
            #print node.asciiArt()
            #print p
        #else:
        #   print "Less than three hosts"
        #   #s_nodes.append(node)
        #   #h_nodes.append(host_subtree.getSubTree(h_names))
        #   #s_tips.append(len(s_vec))
        #   #h_tips.append(len(h_vec))
        #   #p_vals.append('NA')
    
    #DEBUG:
    """
    for i in range(len(p_vals)):
        if p_vals[i] < 0.1:
            print s_nodes[i].asciiArt()
            print h_nodes[i].asciiArt()
            print p_vals[i]
            pause = raw_input("")
    """
    results_dict = {'p_vals':p_vals,'s_tips':s_tips,'h_tips':h_tips,'s_nodes':s_nodes,'h_nodes':h_nodes}
    acc_dict = {}
    return (results_dict, acc_dict)
    
    
    
    
    
    
    
    
    
    

def unifrac_recursive_test(ref_tree,tree,sample_names,
                        taxon_names,data,permutations=1000):#, metric=weighted):
    """Performs UniFrac recursively over a tree.
    
    Specifically, for each node in the tree, performs UniFrac clustering.
    Then compares the UniFrac tree to a reference tree of the same taxa using
    the tip-to-tip distances and the subset distances. Assumption is that if
    the two trees match, the node represents a group in which evolution has
    mirrored the evolution of the reference tree.
    
    tree: contains the tree on which UniFrac will be performed recursively.
    envs: environments for UniFrac clustering (these envs should match the
          taxon labels in the ref_tree)
    ref_tree: reference tree that the clustering is supposed to match.
    metric: metric for UniFrac clustering.
    
    Typically, will want to estimate significance by comparing the actual
    values from ref_tree to values obtained with one or more shuffled versions
    of ref_tree (can make these with permute_tip_labels).
    
    
    Note from Jon: 
    
    I've modified this code a bit to test each node against a set of label-
    permuted host trees, and return some additional information about each node.
    
    It doesn't appear to give sensible results, not sure why. Almost none of the
    resulting permutations yield any other than zero or the number of permuta-
    tions. In other words, every permutation yields either a better or worse 
    match than the true tree. 
    """
    UNIFRAC_CLUST_ENVS = "cluster_envs"
    from cogent.maths.unifrac.fast_unifrac import fast_unifrac
    from cogent.parse.tree import DndParser
    from cogent.core.tree import PhyloNode, TreeError
    from random import shuffle
    from qiime.parse import make_envs_dict
    import numpy
    
    lengths,dists,sets,s_nodes,h_nodes,dist_below,sets_below, h_tips, s_tips = [],[],[],[],[],[],[], [], []
    
    #Permute host tips, store permuted trees in a list of tree strings
    print "Permuting host tree..."
    
    permuted_trees = []
    host_names = ref_tree.getTipNames()
    random_names = ref_tree.getTipNames()
    #for i in range(permutations):
    #   shuffle(random_names)
    #   permute_dict = dict(zip(host_names,random_names))
    #   permuted_subtree = ref_tree.copy()
    #   permuted_subtree.reassignNames(permute_dict)
    #   permuted_trees.append(str(permuted_subtree))
    #   
    #alt:
    for i in range(permutations):
        shuffle(random_names)
        permute_dict = dict(zip(host_names,random_names))
        permuted_subtree = ref_tree.copy()
        permuted_subtree.reassignNames(permute_dict)
        permuted_trees.append(permuted_subtree)
    
    interaction = data.clip(0,1)
    #Parse OTU table data into Unifrac-compatible envs tuple
    
    envs = make_envs_dict(data.T, sample_names, taxon_names)
    
    #Pass host tree, new OTU tree, and envs to recursive unifrac
    print "Performing recursive Unifrac analysis..."
    
    for node in tree.traverse(self_before=True, self_after=False):
        
        #pause = raw_input("pause!")
        #print node
        try:
            result = fast_unifrac(node, envs, weighted=False, modes=set([UNIFRAC_CLUST_ENVS]))
            curr_tree = result[UNIFRAC_CLUST_ENVS]
        except ValueError:
            #hit a single node?
            continue
        except AttributeError:
            #hit a zero branch length
            continue
        if curr_tree is None:
            #hit single node?
            continue
        try:
            l = len(curr_tree.tips())
            d = curr_tree.compareByTipDistances(ref_tree)
            s = curr_tree.compareBySubsets(ref_tree, True)
            
            d_b = 0.0
            s_b = 0.0
            
            #for rand_tree_string in permuted_trees:
            #   rand_tree = DndParser(rand_tree_string)
            #   if d >= curr_tree.compareByTipDistances(rand_tree):
            #       d_b += 1
            #   if s >= curr_tree.compareBySubsets(rand_tree):
            #       s_b += 1
            
            for rand_tree in permuted_trees:
                if d >= curr_tree.compareByTipDistances(rand_tree):
                    d_b += 1
                if s >= curr_tree.compareBySubsets(rand_tree):
                    s_b += 1
            
            d_b = d_b / float(len(permuted_trees))
            s_b = s_b / float(len(permuted_trees))
            
            #The following section generates s_tips and h_tips variables
            #get just OTUs in this node
            otu_subset = node.getTipNames()
            s_tips_tmp = 0
            h_tips_tmp = 0
            s_vec = []
            #find positional index (from OTU table) for each cOTU represented in this node:
            for i in range(len(taxon_names)):
                if taxon_names[i] in otu_subset:
                    s_tips_tmp += 1
                    s_vec.append(i)
            
            #slice interaction matrix down to only cOTUs in this node
            i_s_slice = interaction[numpy.ix_(s_vec)]
            
            #find positional index (this time from OTU table size) for each sample in this node:
            #sum all values in column for each host, if greater than zero, add that host position to h_vec
            for j in range(i_s_slice.shape[1]):
                if i_s_slice[:,j].sum():
                    h_tips_tmp += 1
            
            #want to calculate all values before appending so we can bail out
            #if any of the calculations fails: this ensures that the lists
            #remain synchronized.
            
            
            """
            print curr_tree.asciiArt()
            print ref_tree.asciiArt()
            print l
            print d
            print d_b
            print s
            print s_b
            print node
            
            pause = raw_input("pause!")
            """
            
            if l > 2:
                lengths.append(l)
                dists.append(d)
                sets.append(s)
                s_nodes.append(node)
                h_nodes.append(curr_tree)
                dist_below.append(d_b)
                sets_below.append(s_b)
                h_tips.append(h_tips_tmp)
                s_tips.append(s_tips_tmp)
        except ValueError:
            #no common taxa
            continue
    results_dict = {'p_vals':sets_below,'s_tips':s_tips,'h_tips':h_tips,'s_nodes':s_nodes,'h_nodes':h_nodes}
    
    acc_dict = {'lengths':lengths,'dists':dists,'sets':sets,'dist_below':dist_below}
    
    return (results_dict, acc_dict)
    



#This routine fixes a problem with OTU tables that don't have taxonomic information.
#There's an empty column at the end that you need to remove.
def stupid_hack(summarized_otu_table):
    import re
    p = re.compile('Consensus Lineage')
    q = re.compile('\t\n')
    summarized_otu_table += '\n'
    summarized_otu_table = p.sub('',summarized_otu_table)
    summarized_otu_table = q.sub('\n',summarized_otu_table)
    return summarized_otu_table

#host_subtree, otu_tree, sample_names, taxon_names, data, permutations,host_dist, filtered_seqs

def run_test_cospeciation(basename,
     host_tree,
     host_subtree,
     host_dist,
     otu_subtree,
     sample_names,
     taxon_names,
     data,
     filtered_seqs,
     output_dir,
     test,
     significance_level=0.05,
     permutations=100):
    
    """This function actually tests for cospeciation.
    
    As an input, it takes information about the parent OTU and host phylogeny.
    
    If the parent OTU table was summarized prior to the test, it summarizes the 
    child OTU table according to the same category in the mapping file. 
    
    Then, the set of samples shard by the host tree and parent OTU are retained,
    and any orphaned child OTUs removed. This subsetted host tree and child OTU
    table are passed to the actual test routine. 
    
    The routine expects two dicts 
    to be returned from the test routine: a standard results dict, which is ini-
    tialized in this routine with keys for h_tips, s_tips, p_vals, h_nodes, and 
    s_nodes (host and symbionts [cOTUs] present at each node, significance, and 
    newick-formatted nodes themselves). The values in this dict should be 
    ordered lists of node information. It also accepts an accessory dict with 
    any additional per-node information specific to the test used. These data 
    will be appended to the pOTU results file under column headings determined 
    by dict keys.
    
    Writes a file for each pOTU using these two dicts, with information on each
    node of the cOTU tree on a separate line. 
    
    Returns a boolean True if the test completed, as well as ints for the number
    of nodes tested in that pOTU, and the number of those significant under the 
    given criterea. 
    """
    
    from cogent.app.muscle import align_unaligned_seqs as muscle_aln
    from cogent.app.fasttree import build_tree_from_alignment as fasttree_build_tree
    from cogent.core.tree import PhyloNode
    from cogent import LoadTree, LoadSeqs, DNA
    from cogent.parse.tree import DndParser
    from qiime.summarize_otu_by_cat import summarize_by_cat
    import re
    import numpy
    import StringIO
    
    #DEBUG:
    #print 'in run_test_cospeciation'
    
    
    
    #get number of hosts and cOTUs
    htips = len(host_subtree.getTipNames())
    stips = len(otu_subtree.getTipNames())
    
    #DEBUG:
    #print "got tips" + htips + "_" + stips
    
    #Since we got rid of some stuff, check to see if there are enough hosts and cOTUs
    if htips < 3 or stips < 3:
        results_dict = {'h_tips':htips,'s_tips':stips,'s_nodes':otu_tree,'h_nodes':host_subtree,'p_vals': 1} 
        pvals = 'p_vals'
        sig_nodes = 0
        num_nodes = 0
        print "Host or cOTU tree has less than three tips! "
        
        return True,sig_nodes,num_nodes
    else:
        #DEBUG:
        #print "in test else"
        if test == 'unifrac':
            print 'calling unifrac test'
            results_dict, acc_dict = unifrac_recursive_test(host_subtree,otu_subtree,sample_names,
             taxon_names,data,permutations)
            pvals = 'p_vals'
        
        if test == 'hommola':
            
            
            
            #run recursive hommola test
            results_dict, acc_dict = recursive_hommola(filtered_seqs, host_subtree, host_dist, otu_subtree,sample_names,
             taxon_names,data,permutations)
            
            pvals = 'p_vals'
            
        sig_nodes = 0
        
        
        #Count number of significant nodes
        for pval in results_dict[pvals]:
            if pval < significance_level:
                sig_nodes += 1
        
        results_file = open(output_dir + '/' + basename + '_results.txt', 'w')
        
        keys = results_dict.keys()
        acc_keys = acc_dict.keys()
        for key in keys:
            results_file.write(key + "\t")
        for key in acc_keys:
            results_file.write(key + "\t")
        results_file.write("h_span" + "\t")
        #Write results for each node
        num_nodes = len(results_dict[keys[0]])
        for i in range(num_nodes):
            results_file.write("\n")
            for key in keys:
                results_file.write(str(results_dict[key][i]) + "\t")
            for key in acc_keys:
                results_file.write(str(acc_dict[key][i]) + "\t")
            #calculate 'host span' for each node -- 
            #host span is the minimum number of hosts in the subtree of the original
            #input host tree that is spanned by the hosts included in the cOTU table.
            try:
                h_span = str(len(host_tree.lowestCommonAncestor(results_dict['h_nodes'][i].getTipNames()).getTipNames()))
                results_file.write(h_span)
            except:
                print 'h_span error!'
            
        results_file.close()
        return True, sig_nodes, num_nodes
    

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
    
    
    #Check Host Tree
    try:
        with open(host_tree_fp) as f:
            pass
    except IOError as e:
        print 'Host Data could not be opened! Are you sure it is located at ' + host_tree_fp + '  ?'
        exit(1)
        
    try:
        os.makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)
        
    #Convert inputs to absolute paths
    output_dir=os.path.abspath(output_dir)
    host_tree_fp=os.path.abspath(host_tree_fp)
    mapping_fp=os.path.abspath(mapping_fp)
    potu_table_fp=os.path.abspath(potu_table_fp)
    cotu_table_fp=os.path.abspath(cotu_table_fp)
    
    #Process host input (tree/alignment/matrix) and take subtree of host supertree
    host_tree, host_dist = make_dists_and_tree(potu_table_fp, host_tree_fp)
    
    summary_file = open(output_dir + '/' + 'cospeciation_results_summary.txt', 'w')
    summary_file.write("sig_nodes\tnum_nodes\tfile\n")
    
    #Load taxonomic assignments for the pOTUs
    otu_to_taxonomy = parse_taxonomy(open(taxonomy_fp,'Ur'))
    
    #test that you have a directory, otherwise exit.
    if os.path.isdir(cotu_table_fp):
        os.chdir(cotu_table_fp)
        #run test on all otu tables in directory
        for file in os.listdir('.'):
            if file.endswith("otu_table.txt"):
                #DEBUG:
                #print file
                
                #Infer pOTU name from filename
                cotu_basename = file.split('_')[0];
                print "Analyzing pOTU # " + cotu_basename
                
                result=False
                
                
                cotu_table_fp = cotu_basename + '_seqs_otu_table.txt'
                
                try:
                    cotu_file = open(cotu_table_fp, 'Ur')
                except:
                    print "is this a real file?"
                
                
                basename = cotu_basename + "_" + test
                
                #Import cOTU tree
                otu_tree_fp = cotu_basename + "_seqs_rep_set.tre"
                otu_tree_file = open(otu_tree_fp, 'r')
                otu_tree_unrooted = DndParser(otu_tree_file, PhyloNode)
                otu_tree_file.close()
                
                #root at midpoint
                
                ###Consider alternate step to go through and find closest DB seq to root?
                otu_tree = otu_tree_unrooted.rootAtMidpoint()
                
                #Import host tree
                #host_tree_file = open(host_tree_fp, 'r')
                #host_tree = DndParser(host_tree_file, PhyloNode)
                #host_tree_file.close()
                
                
                #DEBUG:
                #print "basename is:"
                #print basename
                
                ###pool samples per host species (or chosen mapping category)
                
                #DEBUG:
                #print mapping_fp
                
                
                #Open otu table file
                #get 
                
                filtered_cotu_table = filter_samples_from_otu_table(cotu_file,\
                                  host_tree.getTipNames(),\
                                  negate=True)
                
                
                
                sample_names, taxon_names, data, lineages = parse_otu_table(filtered_cotu_table)
                
                if len(sample_names) < 3 or len(taxon_names) < 3:
                    print "Less than 3 hosts or cOTUs in cOTU table!"
                    continue
                
                #make a subtree of the included hosts; exclude hosts not in cOTU table
                host_subtree = host_tree.getSubTree(sample_names)
                otu_subtree = otu_tree.getSubTree(taxon_names)
                
                #Filter the host_dists to match the newly trimmed subtree
                for key in host_dist.keys():
                    if key[0] not in sample_names or key[1] not in sample_names:
                        del host_dist[key]
                #DEBUG:
                #print host_subtree.getTipNames()
                
                #Load up cOTU sequences
                aligned_otu_seqs = LoadSeqs(cotu_basename + '_seqs_rep_set_aligned.fasta', moltype=DNA, label_to_name=lambda x: x.split()[0])
                #Take only ones that made it through host filtering
                filtered_seqs = aligned_otu_seqs.takeSeqs(taxon_names)
                
                #Run recursive test on this pOTU:
                try:
                    #DEBUG:
                    #print 'Trying to run cospeciation test'
                    result, sig_nodes, num_nodes = run_test_cospeciation(
                     basename=basename,
                     host_tree=host_tree,
                     host_subtree=host_subtree,
                     host_dist=host_dist.copy(),
                     otu_subtree=otu_subtree,
                     sample_names=sample_names,
                     taxon_names=taxon_names,
                     data=data,
                     filtered_seqs=filtered_seqs,
                     output_dir=output_dir,
                     test=test,
                     significance_level=significance_level,
                     permutations=permutations)
                except Exception as e:
                    print e
                if result:
                    outline = "{0}\t{1}\t{2}\t{3}".format(sig_nodes,num_nodes,cotu_basename,otu_to_taxonomy[cotu_basename]) + "\n"
                else:
                    outline = "ERROR\t\t" + file + "\n"
                print outline
                summary_file.write(outline)
            
        
    else:
        print 'Not a directory.'
    
    summary_file.close()
    
def make_dists_and_tree(potu_table_fp, host_tree_fp):
    from cogent.app.muscle import align_unaligned_seqs as muscle_aln
    from cogent.app.fasttree import build_tree_from_alignment as fasttree_build_tree
    from cogent.core.tree import PhyloNode
    from cogent import LoadTree, LoadSeqs, DNA
    from cogent.parse.tree import DndParser
    import re
    import numpy
    import StringIO
    from qiime.summarize_otu_by_cat import summarize_by_cat
    from qiime.parse import parse_otu_table
    from cogent.phylo import distance, nj
    from cogent.evolve.models import HKY85
    
    """
    This routine reads in your host information (tree, alignment, or distance 
    matrix) and converts it to a distance matrix and a tree. These are subsetted
    to just the samples present in the pOTU table. The resulting subtree is 
    written to the same directory as the original tree for reference. Both the 
    distance matrix and host subtree are passed back to the main routine for 
    testing.
    """
    
    """
    EDIT:
    
    Make this so that it writes all the files in the output directory for ref.
    """
    
    #Open Parent OTU table
    try:
        potu_file = open(potu_table_fp, 'Ur')
    except:
        print "is this a real file?" + potu_table_fp
    
    #parse pOTU table to grab sample names
    #try:
    sample_names, taxon_names, data, lineages = parse_otu_table(potu_file)
    #except:
    #    print "Problem opening pOTU table. Not a valid OTU table?"
    
    #Initialize host distances dictionary
    host_dist = {}
    
    #Attempt to parse the host tree/alignment/distance matrix
    try:
        #Attempt to load input as tree
        host_supertree = LoadTree(host_tree_fp)
        #remove samples from sample_names that aren't present in supertree
        for x in reversed(range(len(sample_names))):
        #is sample name in supplied host tree?
            if sample_names[x] not in host_supertree.getTipNames():
                del sample_names[x]
        #Trim host supertree
        host_tree = host_supertree.getSubTree(sample_names)
        host_dist = host_tree.getDistances()
        print "Input is a tree"
    
    except Exception as e:
        print "Input is not a tree " 
        print e
        #Attempt to load input as alignment
        try:
            al = LoadSeqs(host_tree_fp)
            d = distance.EstimateDistances(al, submodel= HKY85())
            d.run(show_progress=False)
            host_dist = d.getPairwiseDistances()
            #Delete any distances involving samples not in the pOTU table
            for key in host_dist.keys():
                if key[0] not in sample_names or key[1] not in sample_names:
                    del host_dist[key]
            #generate tree  
            host_tree = nj.nj(host_dist)
            print "Input is a set of sequence alignments"
        except:
            print "Input is not aligned sequences"
            #Attempt to load input as pairwise distances
            try:
                dFile = open(host_tree_fp, 'r')
                dists = dFile.readlines()
                #convert input matrix to array and list
                names, matrix = parse_distmat(dists)
                #Loop through matrix. If both sample names exist, convert to a pairwise distance and add to the dictionary
                for i, item in enumerate(matrix):
                    for j, itemtwo in enumerate(matrix[i]):
                        if i != j:
                            if names[i] in sample_names:
                                if names[j] in sample_names:
                                    host_dist[(names[i], names[j])]  = matrix[i][j]
                #generate tree from distance matrix
                host_tree = nj.nj(host_dist)
                print "Input is a distance matrix"
            except Exception as e:
                print "Input is not a distance matrix"
                print e
                print "Input could not be parsed!"
            
    newPath = host_tree_fp.split('.')[0] + '_filtered.tre'
    
    #Write modified host tree to file. Get rid of single quotes -- is there a
    #better way of doing this?
    try:
    
        host_tree.writeToFile(newPath)
        tmpFile = open(newPath, 'r')        
        tmpString = tmpFile.read()
        tmpFile.close()
        tmpFile = open(newPath, 'w')
        tmpString2 = tmpString.replace("'", "")
        tmpFile.write(tmpString2)
        tmpFile.close()
        #print tmpString2
        
    
    
    except Exception as e:
        print e
    return host_tree, host_dist

if __name__ == "__main__":
    main()

