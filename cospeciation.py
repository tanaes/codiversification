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


import os
import sys
import re
from StringIO import StringIO
import numpy
from random import shuffle

from qiime.util import load_qiime_config, parse_command_line_parameters, get_options_lookup, make_option
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_distmat, make_envs_dict
from qiime.filter import filter_samples_from_otu_table, filter_samples_from_distance_matrix

from biom import load_table

from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent.phylo import distance, nj
from cogent.evolve.models import HKY85
from cogent.evolve.pairwise_distance import TN93Pair
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent import LoadTree, LoadSeqs, DNA
from cogent.util.dict2d import Dict2D, largest


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

    NOTE: the method return signature is now changed.
    For backwards compatibility purposes - 
    when this method is called, 'result' has changed to 'result[0]'
    """
    import cogent.maths.stats.test as stats
    from random import shuffle
    import numpy
    # for testing
    import math

    m = matrix.sum()

    hosts = [0] * m
    pars = [0] * m

    # Generate lists of host and symbiont edges, such that the index
    # of the lists represents an edge connecting the host to the parasite.
    s = 0
    while s < m:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if matrix[i, j] == 1:
                    hosts[s] = j
                    pars[s] = i
                    s += 1

    # get a vector of pairwise distances for each interaction edge
    x = get_dist(hosts, host_dist, range(matrix.shape[1]))
    y = get_dist(pars, par_dist, range(matrix.shape[0]))

    # calculate the observed correlation coefficient for this host/symbionts
    r = stats.correlation(x, y)[0]

    # now do permutaitons. Initialize index lists of the appropriate size.
    mp = range(par_dist.shape[1])
    mh = range(host_dist.shape[1])
    below = 0

    perm_stats = []  # initialize list of shuffled correlation vals

    for i in range(permutations):
        # Generate a shuffled list of indexes for each permutation. This effectively
        # randomizes which host is associated with which symbiont, but maintains
        # the distribution of genetic distances.
        shuffle(mp)
        shuffle(mh)

        # Get pairwise distances in shuffled order
        y_p = get_dist(pars, par_dist, mp)
        x_p = get_dist(hosts, host_dist, mh)

        # calculate shuffled correlation.
        # If greater than observed value, iterate counter below.
        r_p = stats.correlation(x_p, y_p)[0]
        perm_stats.append(r_p)
        if r_p >= r:
            below += 1

    # print "Below: " + str(below)
    # print "Pemutations: " + str(permutations)

    p_val = float(below + 1) / float(permutations + 1)

    return p_val, r, perm_stats


def get_dist(labels, dists, index):
    """Function for picking a subset of pairwise distances from a distance matrix
    according to a set of (randomizable) indices. Derived from Hommola et al R code"""

    m = len(labels)
    vec = []

    for i in range(m - 1):
        k = index[labels[i]]
        for j in range(i + 1, m):
            t = index[labels[j]]
            vec.append(dists[k, t])

    return vec


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

    # print "Performing recursive Hommola et al cospeciation test..."
    
    sample_names = otu_table.ids()
    taxon_names = otu_table.ids(axis="observation")
    otu_data = numpy.asarray([v for v in otu_table.iter_data(axis='sample')])
    
    # calculate pairise distances between OTUs

    dist_calc = TN93Pair(DNA, alignment=aligned_otu_seqs)
    dist_calc.run()

    otu_dists = dist_calc.getPairwiseDistances()

    otu_dm = cogent_dist_to_qiime_dist(otu_dists)

    # convert pw distances (and tree distances for hosts) to numpy arrays with same
    # column/row headings as host/OTU positions in OTU table numpy array.

    #hdd = dist2Dict2D(host_dist,sample_names)
    #hdn = numpy.array(hdd.toLists())
    #sdd = dist2Dict2D(otu_dists,taxon_names)
    #sdn = numpy.array(sdd.toLists())

    # print "got here"

    # print host_dm
    # print sample_names
    # print otu_dm
    # print taxon_names

    host_dm = sort_dm_by_sample(host_dm, sample_names)
    otu_dm = sort_dm_by_sample(otu_dm, taxon_names)

    # convert OTU table to binary array, throwing out all OTUs below a given
    # thresh.

    interaction = otu_data.clip(0, 1)

    # traverse OTU tree and test each node

    # initialize our output lists
    s_nodes, h_nodes, p_vals, s_tips, h_tips, r_vals, r_distro_vals = [
    ], [], [], [], [], [], []
    # print "just before loop"
    # iterate over the tree of child OTUs
    for node in otu_tree.traverse(self_before=True, self_after=False):

        # get just OTUs in this node
        otu_subset = node.getTipNames()

        # subset dms and interaction matrix to just this node
        otu_dm_sub, host_dm_sub, interaction_sub = \
            filter_dms(otu_dm, host_dm, interaction, otu_subset)

        # Make sure we have at least 3 hosts and symbionts represented
        if len(host_dm_sub[0]) > 2 and len(otu_dm_sub[0]) > 2 \
                and host_dm_sub[1].sum() != 0 and otu_dm_sub[1].sum() != 0:

            # print node.asciiArt()

            # append symbiont nodes and host subtrees as tree objects
            s_nodes.append(node)
            h_nodes.append(host_subtree.getSubTree(host_dm_sub[0]))

            # append number of symbionts and hosts for this node
            s_tips.append(len(otu_dm_sub[0]))
            h_tips.append(len(host_dm_sub[0]))
            # calculate pemutation p value for hommola test for this node
            p, r, r_distro = hommola_cospeciation_test(host_dm_sub[1], otu_dm_sub[1],
                                                       interaction_sub, permutations)
            # append to results list
            p_vals.append(p)
            r_vals.append(r)
            r_distro_vals.append(r_distro)

            # print node.asciiArt()
            # print p

        # If only testing top-level node, break out of tree traverse.
        if not recurse:
            break
        # else:
        #   print "Less than three hosts"
        # s_nodes.append(node)
        # h_nodes.append(host_subtree.getSubTree(h_names))
        # s_tips.append(len(s_vec))
        # h_tips.append(len(h_vec))
        # p_vals.append('NA')

    # DEBUG:
    """
    for i in range(len(p_vals)):
        if p_vals[i] < 0.1:
            print s_nodes[i].asciiArt()
            print h_nodes[i].asciiArt()
            print p_vals[i]
            pause = raw_input("")
    """

    # print "finished recursive Hommola"

    results_dict = {'p_vals': p_vals, 's_tips': s_tips,
                    'h_tips': h_tips, 's_nodes': s_nodes, 'h_nodes': h_nodes}
    acc_dict = {'r_vals': r_vals}
    # suppressed: return the distribution of r values
    # 'r_distro_vals':r_distro_vals
    return (results_dict, acc_dict)


def unifrac_recursive_test(ref_tree, tree, sample_names,
                           taxon_names, data, permutations=1000):  # , metric=weighted):
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

    lengths, dists, sets, s_nodes, h_nodes, dist_below, sets_below, h_tips, s_tips = [
    ], [], [], [], [], [], [], [], []

    # Permute host tips, store permuted trees in a list of tree strings
    # print "Permuting host tree..."

    permuted_trees = []
    host_names = ref_tree.getTipNames()
    random_names = ref_tree.getTipNames()
    # for i in range(permutations):
    #   shuffle(random_names)
    #   permute_dict = dict(zip(host_names,random_names))
    #   permuted_subtree = ref_tree.copy()
    #   permuted_subtree.reassignNames(permute_dict)
    #   permuted_trees.append(str(permuted_subtree))
    #
    # alt:
    for i in range(permutations):
        shuffle(random_names)
        permute_dict = dict(zip(host_names, random_names))
        permuted_subtree = ref_tree.copy()
        permuted_subtree.reassignNames(permute_dict)
        permuted_trees.append(permuted_subtree)

    interaction = data.clip(0, 1)
    # Parse OTU table data into Unifrac-compatible envs tuple

    envs = make_envs_dict(data.T, sample_names, taxon_names)

    # Pass host tree, new OTU tree, and envs to recursive unifrac
    # print "Performing recursive Unifrac analysis..."

    for node in tree.traverse(self_before=True, self_after=False):

        #pause = raw_input("pause!")
        # print node
        try:
            result = fast_unifrac(
                node, envs, weighted=False, modes=set([UNIFRAC_CLUST_ENVS]))
            curr_tree = result[UNIFRAC_CLUST_ENVS]
        except ValueError:
            # hit a single node?
            continue
        except AttributeError:
            # hit a zero branch length
            continue
        if curr_tree is None:
            # hit single node?
            continue
        try:
            l = len(curr_tree.tips())
            d = curr_tree.compareByTipDistances(ref_tree)
            s = curr_tree.compareBySubsets(ref_tree, True)

            d_b = 0.0
            s_b = 0.0

            # for rand_tree_string in permuted_trees:
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

            # The following section generates s_tips and h_tips variables
            # get just OTUs in this node
            otu_subset = node.getTipNames()
            s_tips_tmp = 0
            h_tips_tmp = 0
            s_vec = []
            # find positional index (from OTU table) for each cOTU represented
            # in this node:
            for i in range(len(taxon_names)):
                if taxon_names[i] in otu_subset:
                    s_tips_tmp += 1
                    s_vec.append(i)

            # slice interaction matrix down to only cOTUs in this node
            i_s_slice = interaction[numpy.ix_(s_vec)]

            # find positional index (this time from OTU table size) for each sample in this node:
            # sum all values in column for each host, if greater than zero, add
            # that host position to h_vec
            for j in range(i_s_slice.shape[1]):
                if i_s_slice[:, j].sum():
                    h_tips_tmp += 1

            # want to calculate all values before appending so we can bail out
            # if any of the calculations fails: this ensures that the lists
            # remain synchronized.

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
            # no common taxa
            continue
    results_dict = {'p_vals': sets_below, 's_tips': s_tips,
                    'h_tips': h_tips, 's_nodes': s_nodes, 'h_nodes': h_nodes}

    acc_dict = {'lengths': lengths, 'dists': dists,
                'sets': sets, 'dist_below': dist_below}

    return (results_dict, acc_dict)


def make_dists_and_tree(sample_names, host_fp):
    """
    This routine reads in your host information (tree, alignment, or distance 
    matrix) and converts it to a distance matrix and a tree. These are subsetted
    to just the samples passed to the routine. The resulting subtree is 
    written to the same directory as the original tree for reference. Both the 
    distance matrix and host subtree are passed back to the main routine for 
    testing.
    """
    hostf = open(host_fp, 'r')
    host_str = hostf.read()
    hostf.close()

    # Attempt to parse the host tree/alignment/distance matrix
    if isTree(host_str):
        host_tree, host_dist = processTree(host_str)
        print "Input is tree"

    elif isAlignment(host_str):
        host_tree, host_dist = processAlignment(host_str)
        print "Input is alignment"

    elif isMatrix(host_str):
        host_tree, host_dist = processMatrix(host_str)
        print "Input is distance matrix"

    else:
        print "Host information file could not be parsed"

    # Remove any sample names not in host tree
    sample_names = filter(
        lambda x: x if x in host_tree.getTipNames() else None, sample_names)
    print sample_names
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

    dm_names_dict = dict([(x[1], x[0]) for x in enumerate(sample_names)])
    name_slice = []

    for name in dm[0]:
        name_slice.append(dm_names_dict[name])

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
        if i_s_slice[:, j].sum():
            h_vec.append(j)
            h_names.append(host_dm[0][j])

    # check to see that the host vector isn't empty
    if len(h_vec) < 1:
        return(([], []), ([], []), ([]))

    i_slice = interaction[numpy.ix_(s_vec, h_vec)]

    # slice host distance matrix
    h_slice = host_dm[1][numpy.ix_(h_vec, h_vec)]

    sliced_host_dm = (h_names, h_slice)
    sliced_otu_dm = (s_names, s_slice)
    sliced_interaction = i_slice

    return(sliced_otu_dm, sliced_host_dm, sliced_interaction)


def isTree(fstr):

    try:
        LoadTree(treestring=fstr)
        return True

    except:
        return False


def isAlignment(fstr):

    try:
        LoadSeqs(data=fstr)
        return True

    except:
        return False


def isMatrix(fstr):

    try:
        result = parse_distmat(fstr.splitlines())
        if result[0] == None:
            return False

        else:
            return True

    except:
        return False


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


def filter_otu_table_by_min(sample_names, taxon_names, data, lineages, min=1):
    # loop through and remove every otu present at less than min
    # this should be replaced by native QIIME filter code.
    s_vec = []
    s_names = []
    taxonomies = []

    # create list of OTUs to keep
    for otu in range(data.shape[0]):
        if data[otu, :].sum() >= min:
            s_vec.append(otu)
            s_names.append(taxon_names[otu])
            if lineages:
                taxonomies.append(lineages[otu])

    h_vec = []
    h_names = []

    for sample in range(data.shape[1]):
        if data[numpy.ix_(s_vec), sample].sum() >= 1:
            h_vec.append(sample)
            h_names.append(sample_names[sample])

    # slice data
    data = data[numpy.ix_(s_vec, h_vec)]

    return h_names, s_names, data, taxonomies


def write_results(results_dict, acc_dict, output_dir, potu, test, host_tree):
    # print results_dict
    # print acc_dict

    results_file = open(join(output_dir,('%s_%s_results.txt' % (potu, test))), 'w')

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

    # filter cOTU table by samples present in host_tree/dm

    cotu_table_filtered = cotu_table.filter(host_dist[0], axis = 'sample', inplace = False)
    

    # filter cOTU table again to get rid of absent cOTUs

    filtered_cotu_table.filter(lambda val, id_, metadata: 1 <= val.sum(), axis='observation', inplace = True)
    
    # Filter the host_dists to match the newly trimmed subtree
    # Note: this is requiring the modified filter_dist method which
    # returns a native dm tuple rather than a string.

    host_dist_filtered = filter_samples_from_distance_matrix(
        host_dist, sample_names, negate=True)

    

    return cotu_table_filtered, host_dist_filtered


def test_cospeciation(potu_table_fp, subcluster_dir, host_tree_fp, mapping_fp, mapping_category, output_dir, significance_level, test, permutations, taxonomy_fp, force):

    # Convert inputs to absolute paths
    output_dir = os.path.abspath(output_dir)
    host_tree_fp = os.path.abspath(host_tree_fp)
    mapping_fp = os.path.abspath(mapping_fp)
    potu_table_fp = os.path.abspath(potu_table_fp)
    subcluster_dir = os.path.abspath(subcluster_dir)

    # Check Host Tree
    try:
        with open(host_tree_fp) as f:
            pass

    except IOError as e:
        print 'Host Data could not be opened! Are you sure it is located at ' + host_tree_fp + '  ?'
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
    host_tree, host_dist = make_dists_and_tree(sample_names, host_tree_fp)

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
            cotu_tree_fp = os.path.join(potu,"rep_set.tre")
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

            # Load up and filter cOTU sequences
            aligned_otu_seqs = LoadSeqs(
                os.path.join(potu,'seqs_rep_set_aligned.fasta'), moltype=DNA, label_to_name=lambda x: x.split()[0])
            cotu_seqs_filtered = aligned_otu_seqs.takeSeqs(cotu_names_filtered)

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