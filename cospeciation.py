from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jon Sanders", "Nate Bresnick", "Aaron Behr"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Development"

import os
import sys
import re
from StringIO import StringIO
import numpy
from random import shuffle
from operator import add

from qiime.util import load_qiime_config, parse_command_line_parameters, get_options_lookup, make_option
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_distmat, make_envs_dict, fields_to_dict
from qiime.filter import filter_samples_from_otu_table, filter_samples_from_distance_matrix
from qiime.workflow.upstream import run_pick_de_novo_otus

from qiime.workflow.util import call_commands_serially, no_status_updates

from biom import load_table

from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent.phylo import distance, nj
from cogent.evolve.models import HKY85
from cogent.evolve.pairwise_distance import TN93Pair
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent import LoadTree, LoadSeqs, DNA
from cogent.util.dict2d import Dict2D, largest

from skbio.math.stats.evolve.hommola import hommola_cospeciation

def cogent_dist_to_qiime_dist(dist_tuple_dict):
    """
    This takes a dict with tuple keys and distance values, such as is output
    by the getDistances() method of a PhyloNode object, and converts it to a 
    QIIME-style distance matrix object: an ordered tuple with a list of samples
    in [0] and a numpy array of the distance matrix in [1].

    EDITED AND UPDATED 2013-07-09 Aaron Behr
    """

    headers = []
    dist_dict = {}

    # loop through dist_tuple_dict, returning (k1,k2):v tuples simultaneously
    for item in dist_tuple_dict.iteritems():
        # if k1 is not in headers, add it to headers
        if item[0][0] not in headers:
            headers.append(item[0][0])
            dist_dict[item[0][0]] = {item[0][0]: 0.0}  # null self-distance

        dist_dict[item[0][0]][item[0][1]] = item[1]  # dist_dict[k1][k2] = v
    headers.sort()

    # Initialize dict2d, with data from dist_dict (dict of dicts).
    # Also, RowOrder and ColOrder are set to the order of the sorted headers.
    # NOTE: no longer using the fromDicts() method to pass dist_dict to dict2d
    dict2d = Dict2D(dist_dict, headers, headers)

    # reflect dict2d so that it is no longer sparse
    dict2d.reflect(largest)

    # output tab-delimited printable string of the items in dict2d including
    # headers.
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

def cache_tipnames(tree):
    for n in tree.postorder(include_self=True):
        if n.isTip():
            n._tip_names = [n.Name]
        else:
            n._tip_names = reduce(add, [c._tip_names for c in n.Children])

def recursive_hommola(aligned_otu_seqs, host_subtree, host_dm, otu_tree, otu_table, 
                permutations=1000, recurse=False):
    """
    Applies Hommola et al test of cospeciation recursively to OTU tree.

    Conceptually similar to the oraganization of the recursive Unifrac method.

    Host distances are calculated from the provided host tree, and OTU distances
    from the MUSCLE alignment using the TN93 model of nucleotide evolution. It 
    would probably be better to do pairwise alignments and calculate distances
    that way. 

    It returns a dictionary of the results, and an empty accessory dict.
    """

    sample_names = otu_table.ids()
    taxon_names = otu_table.ids(axis="observation")
    presence_absence = otu_table.pa(inplace=False)
    interaction = numpy.asarray(list(presence_absence.iter_data(axis='observation')))
    
    # calculate pairise distances between OTUs

    dist_calc = TN93Pair(DNA, alignment=aligned_otu_seqs)
    dist_calc.run()

    otu_dists = dist_calc.getPairwiseDistances()

    otu_dm = cogent_dist_to_qiime_dist(otu_dists)

    host_dm = sort_dm_by_sample(host_dm, sample_names)
    otu_dm = sort_dm_by_sample(otu_dm, taxon_names)

    # convert OTU table to binary array, throwing out all OTUs below a given
    # thresh.

    # traverse OTU tree and test each node

    # initialize our output lists
    s_nodes = []
    h_nodes = []
    p_vals = []
    s_tips = []
    h_tips = []
    r_vals = []
    r_distro_vals = []
    # print "just before loop"
    # iterate over the tree of child OTUs

    cache_tipnames(otu_tree)

    for node in otu_tree.preorder():

        # get just OTUs in this node
        otu_subset = node._tip_names

        # subset dms and interaction matrix to just this node
        otu_dm_sub, host_dm_sub, interaction_sub = \
            filter_dms(otu_dm, host_dm, interaction, otu_subset)

        # Make sure we have at least 3 hosts and symbionts represented
        if len(host_dm_sub[0]) > 2 and len(otu_dm_sub[0]) > 2 \
                and host_dm_sub[1].sum() != 0 and otu_dm_sub[1].sum() != 0:

            # append symbiont nodes and host subtrees as tree objects
            s_nodes.append(node)
            h_nodes.append(host_subtree.getSubTree(host_dm_sub[0]))

            # append number of symbionts and hosts for this node
            s_tips.append(len(otu_dm_sub[0]))
            h_tips.append(len(host_dm_sub[0]))
            # calculate pemutation p value for hommola test for this node
            p, r, r_distro = hommola_cospeciation(host_dm_sub[1], otu_dm_sub[1],
                                                       interaction_sub, permutations)
            # append to results list
            p_vals.append(p)
            r_vals.append(r)
            r_distro_vals.append(r_distro)

        # If only testing top-level node, break out of tree traverse.
        if not recurse:
            break

    results_dict = {'p_vals': p_vals, 's_tips': s_tips,
                    'h_tips': h_tips, 's_nodes': s_nodes, 'h_nodes': h_nodes}
    acc_dict = {'r_vals': r_vals}
    # suppressed: return the distribution of r values
    # 'r_distro_vals':r_distro_vals
    return (results_dict, acc_dict)

def make_dists_and_tree(sample_names, host_fp, host_input_type):
    """
    This routine reads in your host information (tree, alignment, or distance 
    matrix) and converts it to a distance matrix and a tree. These are subsetted
    to just the samples passed to the routine. The resulting subtree is 
    written to the same directory as the original tree for reference. Both the 
    distance matrix and host subtree are passed back to the main routine for 
    testing.
    """
    with open(host_fp, 'r') as host_f:
        host_str = hostf.read()

    # Attempt to parse the host tree/alignment/distance matrix
    if host_input_type == "tree":
        host_tree, host_dist = processTree(host_str)

    elif host_input_type == "alignment":
        host_tree, host_dist = processAlignment(host_str)

    elif host_input_type == "distances":
        host_tree, host_dist = processMatrix(host_str)

    sample_names = [x for x in sample_names if x in host_tree.getTipNames()]
    # Get host subtree and filter distance matrix so they only include samples
    # present in the pOTU table
    host_tree = host_tree.getSubTree(sample_names)

    host_dist = filter_samples_from_distance_matrix(
        host_dist, sample_names, negate=True)
    return host_tree, host_dist


# This function is copied directly from QIIME except it returns a native
# distance matrix instead of formatting it
def filter_samples_from_distance_matrix(dm, samples_to_discard, negate=False):
    from numpy import array, inf
    """ Remove specified samples from distance matrix 
    
        dm: (sample_ids, dm_data) tuple, as returned from 
         qiime.parse.parse_distmat; or a file handle that can be passed
         to qiime.parse.parse_distmat
    
    """
    try:
        sample_ids, dm_data = dm

    except ValueError:
        # input was provide as a file handle
        sample_ids, dm_data = parse_distmat(dm)

    sample_lookup = {}.fromkeys([e.split()[0] for e in samples_to_discard])
    temp_dm_data = []
    new_dm_data = []
    new_sample_ids = []

    if negate:
        def keep_sample(s):
            return s in sample_lookup

    else:

        def keep_sample(s):
            return s not in sample_lookup

    for row, sample_id in zip(dm_data, sample_ids):

        if keep_sample(sample_id):
            temp_dm_data.append(row)
            new_sample_ids.append(sample_id)

    temp_dm_data = array(temp_dm_data).transpose()

    for col, sample_id in zip(temp_dm_data, sample_ids):

        if keep_sample(sample_id):
            new_dm_data.append(col)

    new_dm_data = array(new_dm_data).transpose()

    return (new_sample_ids, new_dm_data)


def sort_dm_by_sample(dm, sample_names):
    """Sorts a qiime distance matrix tuple in the order of sample names given"""

    dm_names_dict = {x:i for i, x in enumerate(sample_names)}

    name_slice = [dm_names_dict[name] for name in dm[0]]

    sorted_dm = dm[1][numpy.ix_(name_slice, name_slice)]

    return (sample_names, sorted_dm)


def filter_dms(otu_dm, host_dm, interaction, otu_subset):
    """This filters a host dm, symbiont dm, and interaction matrix by a set of
    sybionts (otus) defined by otu_subset, and returns the sliced values.
    Also eliminates any hosts that had no otus present."""
    # input host dm, symbiont dm, and otu data

    # return filtered dms,
    s_vec = []
    h_vec = []
    h_names = []
    s_names = []

    # find positional index (from OTU table) for each cOTU represented in this
    # node:
    for i in range(len(otu_dm[0])):
        if otu_dm[0][i] in otu_subset:
            s_vec.append(i)
            s_names.append(otu_dm[0][i])
    # slice symbiont distance matrix down to only cOTUs in this node

    s_slice = otu_dm[1][numpy.ix_(s_vec, s_vec)]

    # slice interaction matrix down to only cOTUs in this node
    i_s_slice = interaction[numpy.ix_(s_vec)]

    # find positional index (this time from OTU table size) for each sample in this node:
    # sum all values in column for each host, if greater than zero, add that
    # host position to h_vec
    for j in range(i_s_slice.shape[1]):
        if i_s_slice[:, j].any():
            h_vec.append(j)
            h_names.append(host_dm[0][j])

    # check to see that the host vector isn't empty
    if not h_vec:
        return(([], []), ([], []), [])

    i_slice = interaction[numpy.ix_(s_vec, h_vec)]

    # slice host distance matrix
    h_slice = host_dm[1][numpy.ix_(h_vec, h_vec)]

    sliced_host_dm = (h_names, h_slice)
    sliced_otu_dm = (s_names, s_slice)
    sliced_interaction = i_slice

    return(sliced_otu_dm, sliced_host_dm, sliced_interaction)

def processTree(fstr):
    # Attempt to load input as tree
    host_tree = LoadTree(treestring=fstr)
    host_dist = cogent_dist_to_qiime_dist(host_tree.getDistances())
    return host_tree, host_dist


def processAlignment(fstr):
    # load sequences and estimate distance matrix
    al = LoadSeqs(data=fstr)
    d = distance.EstimateDistances(al, submodel=HKY85())
    d.run(show_progress=False)
    host_dist = cogent_dist_to_qiime_dist(d.getPairwiseDistances())

    # generate tree from matrix
    host_tree = distmat_to_tree(host_dist)
    return host_tree, host_dist


def processMatrix(fstr):
    dists = fstr.splitlines()
    # Parse distance matrix and build tree
    host_dist = parse_distmat(dists)
    host_tree = distmat_to_tree(host_dist)
    return host_tree, host_dist


def distmat_to_tree(distmat):
    dist_headers, dist_matrix = distmat
    cogent_host_dist = {}
    # Loop through host distance matrix to create a dictionary of pairwise
    # distances
    for i, item in enumerate(dist_matrix):
        for j, itemtwo in enumerate(dist_matrix[i]):
            if i != j:
                cogent_host_dist[
                    (dist_headers[i], dist_headers[j])] = dist_matrix[i][j]
    # Generate tree from distance matrix

    return nj.nj(cogent_host_dist)

def write_results(results_dict, acc_dict, output_dir, potu, test, host_tree):
    # print results_dict
    # print acc_dict

    results_file = open(os.path.join(output_dir,('%s_%s_results.txt' % (potu, test))), 'w')

    keys = results_dict.keys()
    acc_keys = acc_dict.keys()
    for key in keys:
        results_file.write(key + "\t")

    for key in acc_keys:
        results_file.write(key + "\t")

    results_file.write("h_span" + "\t")
    # Write results for each node
    num_nodes = len(results_dict[keys[0]])

    for i in range(num_nodes):
        results_file.write("\n")

        for key in keys:
            results_file.write(str(results_dict[key][i]) + "\t")

        for key in acc_keys:
            results_file.write(str(acc_dict[key][i]) + "\t")

        h_span = calc_h_span(host_tree, results_dict, i)
        results_file.write(h_span)

    results_file.close()
    return num_nodes


def calc_h_span(host_tree, results_dict, i):
    # calculate 'host span' for each node --
    # host span is the minimum number of hosts in the subtree of the original
    # input host tree that is spanned by the hosts included in the cOTU table.

    try:
        h_span = str(len(host_tree.lowestCommonAncestor(
            results_dict['h_nodes'][i].getTipNames()).getTipNames()))
    except:
        print 'h_span error!'
    return h_span


def reconcile_hosts_symbionts(cotu_table, host_dist):

    shared_hosts = set(cotu_table.ids(axis='sample')).intersection(host_dist[0])

    # filter cOTU table by samples present in host_tree/dm

    cotu_table_filtered = cotu_table.filter(shared_hosts, axis='sample', inplace=False)
    

    # filter cOTU table again to get rid of absent cOTUs. BIOM should do this.

    cotu_table_filtered.filter(lambda val, id_, metadata: 1 <= val.sum(), axis='observation', inplace=True)
    
    # Filter the host_dists to match the newly trimmed subtree
    # Note: this is requiring the modified filter_dist method which
    # returns a native dm tuple rather than a string.

    host_dist_filtered = filter_samples_from_distance_matrix(
        host_dist, cotu_table_filtered.ids(axis = 'sample'), negate=True)

    

    return cotu_table_filtered, host_dist_filtered