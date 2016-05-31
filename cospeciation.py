
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

from qiime.util import load_qiime_config, parse_command_line_parameters, get_options_lookup, make_option, write_biom_table
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_distmat, make_envs_dict, fields_to_dict
from qiime.filter import filter_samples_from_otu_table, filter_samples_from_distance_matrix
from qiime.workflow.upstream import run_pick_de_novo_otus
from qiime.stats import benjamini_hochberg_step_down, bonferroni_correction, fdr_correction
from qiime.workflow.util import call_commands_serially, no_status_updates
from qiime.group import collapse_samples

from biom import load_table

from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent.phylo import distance, nj
from cogent.evolve.models import HKY85
from cogent.evolve.pairwise_distance import TN93Pair
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent import LoadTree, LoadSeqs, DNA
from cogent.util.dict2d import Dict2D, largest

#from skbio.stats.evolve import hommola_cospeciation

from skbio import DistanceMatrix
from scipy.stats import pearsonr
import numpy as np

def hommola_cospeciation(host_dist, par_dist, interaction, permutations=999):
    """Perform Hommola et al (2009) host/parasite cospeciation test.
    This test for host/parasite cospeciation is as described in [1]_. This test
    is a modification of a Mantel test, expanded to accept the case where
    multiple hosts map to a single parasite (and vice versa).
    For a basic Mantel test, the distance matrices being compared must have the
    same number of values. To determine the significance of the correlations
    between distances in the two matrices, the correlation coefficient of those
    distances is calculated and compared to the correlation coefficients
    calculated from a set of matrices in which rows and columns have been
    permuted.
    In this test, rather than comparing host-host to parasite-parasite
    distances directly (requiring one host per parasite), the distances are
    compared for each interaction edge between host and parasite. Thus, a host
    interacting with two different parasites will be represented in two
    different edges, with the host-host distance for the comparison between
    those edges equal to zero, and the parasite-parasite distance equal to the
    distance between those two parasites. Like in the Mantel test, significance
    of the interaction is assessed by permutation, in this case permutation of
    the host-symbiont interaction links.
    Note that the null hypothesis being tested here is that the hosts and
    parasites have evolved independently of one another. The alternative to
    this is a somewhat weaker case than what is often implied with the term
    'cospeciation,' which is that each incidence of host speciation is
    recapitulated in an incidence of symbiont speciation (strict
    co-cladogenesis). Although there may be many factors that could contribute
    to non-independence of host and symbiont phylogenies, this loss of
    explanatory specificity comes with increased robustness to phylogenetic
    uncertainty. Thus, this test may be especially useful for cases where host
    and/or symbiont phylogenies are poorly resolved, or when simple correlation
    between host and symbiont evolution is of more interest than strict
    co-cladogenesis.
    This test requires pairwise distance matrices for hosts and symbionts, as
    well as an interaction matrix specifying links between hosts (in columns)
    and symbionts (in rows). This interaction matrix should have the same
    number of columns as the host distance matrix, and the same number of rows
    as the symbiont distance matrix. Interactions between hosts and symbionts
    should be indicated by values of ``1`` or ``True``, with non-interactions
    indicated by values of ``0`` or ``False``.
    Parameters
    ----------
    host_dist : 2-D array_like or DistanceMatrix
        Symmetric matrix of m x m pairwise distances between hosts.
    par_dist : 2-D array_like or DistanceMatrix
        Symmetric matrix of n x n pairwise distances between parasites.
    interaction : 2-D array_like, bool
        n x m binary matrix of parasite x host interactions. Order of hosts
        (columns) should be identical to order of hosts in `host_dist`, as
        should order of parasites (rows) be identical to order of parasites in
        `par_dist`.
    permutations : int, optional
        Number of permutations used to compute p-value. Must be greater than or
        equal to zero. If zero, statistical significance calculations will be
        skipped and the p-value will be ``np.nan``.
    Returns
    -------
    corr_coeff : float
        Pearson correlation coefficient of host : parasite association.
    p_value : float
        Significance of host : parasite association computed using
        `permutations` and a one-sided (greater) alternative hypothesis.
    perm_stats : 1-D numpy.ndarray, float
        Correlation coefficients observed using permuted host : parasite
        interactions. Length will be equal to the number of permutations used
        to compute p-value (see `permutations` parameter above).
    See Also
    --------
    skbio.stats.distance.mantel
    scipy.stats.pearsonr
    Notes
    -----
    It is assumed that the ordering of parasites in `par_dist` and hosts in
    `host_dist` are identical to their ordering in the rows and columns,
    respectively, of the interaction matrix.
    This code is loosely based on the original R code from [1]_.
    References
    ----------
    .. [1] Hommola K, Smith JE, Qiu Y, Gilks WR (2009) A Permutation Test of
       Host-Parasite Cospeciation. Molecular Biology and Evolution, 26,
       1457-1468.
    Examples
    --------
    >>> from skbio.stats.evolve import hommola_cospeciation
    Create arrays for host distances, parasite distances, and their
    interactions (data taken from example in [1]_):
    >>> hdist = [[0,3,8,8,9], [3,0,7,7,8], [8,7,0,6,7], [8,7,6,0,3],
    ...          [9,8,7,3,0]]
    >>> pdist = [[0,5,8,8,8], [5,0,7,7,7], [8,7,0,4,4], [8,7,4,0,2],
    ...          [8,7,4,2,0]]
    >>> interaction = [[1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0],
    ...                [0,0,0,1,1]]
    Run the cospeciation test with 99 permutations. Note that the correlation
    coefficient for the observed values counts against the final reported
    p-value:
    >>> corr_coeff, p_value, perm_stats = hommola_cospeciation(
    ...     hdist, pdist, interaction, permutations=99)
    >>> corr_coeff
    0.83170965463247903
    In this case, the host distances have a fairly strong positive correlation
    with the symbiont distances. However, this may also reflect structure
    inherent in the phylogeny, and is not itself indicative of significance.
    >>> p_value <= 0.05
    True
    After permuting host : parasite interactions, we find that the observed
    correlation is indeed greater than we would expect by chance.
    """
    host_dist = DistanceMatrix(host_dist)
    par_dist = DistanceMatrix(par_dist)
    interaction = np.asarray(interaction, dtype=bool)

    num_hosts = host_dist.shape[0]
    num_pars = par_dist.shape[0]

    if num_hosts < 3 or num_pars < 3:
        raise ValueError("Distance matrices must be a minimum of 3x3 in size.")
    if num_hosts != interaction.shape[1]:
        raise ValueError("Number of interaction matrix columns must match "
                         "number of hosts in `host_dist`.")
    if num_pars != interaction.shape[0]:
        raise ValueError("Number of interaction matrix rows must match "
                         "number of parasites in `par_dist`.")
    if permutations < 0:
        raise ValueError("Number of permutations must be greater than or "
                         "equal to zero.")
    if interaction.sum() < 3:
        raise ValueError("Must have at least 3 host-parasite interactions in "
                         "`interaction`.")

    # shortcut to eliminate nested for-loops specifying pairwise interaction
    # partners as randomizeable indices
    pars, hosts = np.nonzero(interaction)
    pars_k_labels, pars_t_labels = _gen_lists(pars)
    hosts_k_labels, hosts_t_labels = _gen_lists(hosts)

    # get a vector of pairwise distances for each interaction edge
    x = _get_dist(hosts_k_labels, hosts_t_labels, host_dist.data,
                  np.arange(num_hosts))
    y = _get_dist(pars_k_labels, pars_t_labels, par_dist.data,
                  np.arange(num_pars))

    # calculate the observed correlation coefficient for these hosts/symbionts
    corr_coeff = pearsonr(x, y)[0]

    # now do permutatitons. initialize index lists of the appropriate size
    mp = np.arange(num_pars)
    mh = np.arange(num_hosts)

    # initialize list of shuffled correlation vals
    perm_stats = np.empty(permutations)

    if permutations == 0 or np.isnan(corr_coeff):
        p_value = np.nan
        perm_stats.fill(np.nan)
    else:
        for i in range(permutations):
            # generate a shuffled list of indexes for each permutation. this
            # effectively randomizes which host is associated with which
            # symbiont, but maintains the distribution of genetic distances
            np.random.shuffle(mp)
            np.random.shuffle(mh)

            # get pairwise distances in shuffled order
            y_p = _get_dist(pars_k_labels, pars_t_labels, par_dist.data, mp)
            x_p = _get_dist(hosts_k_labels, hosts_t_labels, host_dist.data, mh)

            # calculate shuffled correlation coefficient
            perm_stats[i] = pearsonr(x_p, y_p)[0]

        p_value = ((perm_stats >= corr_coeff).sum() + 1) / (permutations + 1)

    return corr_coeff, p_value, perm_stats


def hommola_cospeciation_host(host_dist, par_dist, interaction, permutations=999):
    """Perform Hommola et al (2009) host/parasite cospeciation test.

    Performs a modification of the Hommola et al cospeciation test in which
    only the host distance matrix is permuted, leaving all relationships
    among host and symbionts the same. 

    """
    host_dist = DistanceMatrix(host_dist)
    par_dist = DistanceMatrix(par_dist)
    interaction = np.asarray(interaction, dtype=bool)

    num_hosts = host_dist.shape[0]
    num_pars = par_dist.shape[0]

    if num_hosts < 3 or num_pars < 3:
        raise ValueError("Distance matrices must be a minimum of 3x3 in size.")
    if num_hosts != interaction.shape[1]:
        raise ValueError("Number of interaction matrix columns must match "
                         "number of hosts in `host_dist`.")
    if num_pars != interaction.shape[0]:
        raise ValueError("Number of interaction matrix rows must match "
                         "number of parasites in `par_dist`.")
    if permutations < 0:
        raise ValueError("Number of permutations must be greater than or "
                         "equal to zero.")
    if interaction.sum() < 3:
        raise ValueError("Must have at least 3 host-parasite interactions in "
                         "`interaction`.")

    # shortcut to eliminate nested for-loops specifying pairwise interaction
    # partners as randomizeable indices
    pars, hosts = np.nonzero(interaction)
    pars_k_labels, pars_t_labels = _gen_lists(pars)
    hosts_k_labels, hosts_t_labels = _gen_lists(hosts)

    # get a vector of pairwise distances for each interaction edge
    x = _get_dist(hosts_k_labels, hosts_t_labels, host_dist.data,
                  np.arange(num_hosts))
    y = _get_dist(pars_k_labels, pars_t_labels, par_dist.data,
                  np.arange(num_pars))

    # calculate the observed correlation coefficient for these hosts/symbionts
    corr_coeff = pearsonr(x, y)[0]

    # initialize list of shuffled correlation vals
    perm_stats = np.empty(permutations)

    if permutations == 0 or np.isnan(corr_coeff):
        p_value = np.nan
        perm_stats.fill(np.nan)
    else:
        for i in range(permutations):
            # generate a shuffled host distance matrix
            host_dist_perm = host_dist.permute()
      
            # get pairwise distances in shuffled order
            x_p = _get_dist(hosts_k_labels, hosts_t_labels, host_dist_perm.data,
                          np.arange(num_hosts))
            y = _get_dist(pars_k_labels, pars_t_labels, par_dist.data,
                          np.arange(num_pars))

            # calculate shuffled correlation coefficient
            perm_stats[i] = pearsonr(x_p, y)[0]

        p_value = ((perm_stats >= corr_coeff).sum() + 1) / (permutations + 1)

    return corr_coeff, p_value, perm_stats


def _get_dist(k_labels, t_labels, dists, index):
    """Subset a distance matrix using a set of (randomizable) index labels.
    Parameters
    ----------
    k_labels : numpy.array
        index labels specifying row-wise member of pairwise interaction
    t_labels : numpy.array
        index labels specifying column-wise member of pairwise interaction
    dists : numpy.array
        pairwise distance matrix
    index : numpy.array of int
        permutable indices for changing order in pairwise distance matrix
    Returns
    -------
    vec : list of float
        List of distances associated with host:parasite edges.
    """
    return dists[index[k_labels], index[t_labels]]


def _gen_lists(labels):
    """Generate matched lists of row and column index labels.
    Shortcut function for generating matched lists of row and col index
    labels for the set of pairwise comparisons specified by the list of those
    indices recovered using ``np.nonzero(interaction)``.
    Reproduces values of iterated indices from the nested for-loops contained
    in ``get_dist`` function in original code from [1]_.
    Parameters
    ----------
    labels : numpy.array
        array containing the indices of nonzero elements in one dimension of an
        interaction matrix
    Returns
    -------
    k_labels : numpy.array
        index labels specifying row-wise member of pairwise interaction
    t_labels : numpy.array
        index labels specifying column-wise member of pairwise interaction
    References
    ----------
    .. [1] Hommola K, Smith JE, Qiu Y, Gilks WR (2009) A Permutation Test of
       Host-Parasite Cospeciation. Molecular Biology and Evolution, 26,
       1457-1468.
    """
    i_array, j_array = np.transpose(np.tri(len(labels)-1)).nonzero()
    j_array = j_array + 1
    return labels[i_array], labels[j_array]


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

def cache_tipnames(tree):
    """
    Function to traverse a tree and store dependent tip names for each node.
    Replaces the getTipNames() method in cogent, should be replaced if/when this
    functionality is native to skbio. 
    """

    for n in tree.postorder(include_self=True):
        if n.isTip():
            n._tip_names = [n.Name]
        else:
            n._tip_names = reduce(add, [c._tip_names for c in n.Children])

def recursive_hommola(aligned_otu_seqs, host_subtree, host_dm, otu_tree, otu_table, 
                permutations=10000, perm_type='hommola', recurse=False):
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

    # convert OTU table to presence/absence interaction matrix
    presence_absence = otu_table.pa(inplace=False)
    interaction = numpy.asarray(list(presence_absence.iter_data(axis='observation')))
    
    # calculate pairise distances between OTUs
    dist_calc = TN93Pair(DNA, alignment=aligned_otu_seqs)
    dist_calc.run()

    otu_dists = dist_calc.getPairwiseDistances()

    otu_dm = cogent_dist_to_qiime_dist(otu_dists)

    # Sort dms into same order as sample and taxon name lists
    host_dm = sort_dm_by_sample(host_dm, sample_names)
    otu_dm = sort_dm_by_sample(otu_dm, taxon_names)

    # initialize our output lists
    s_nodes = []
    h_nodes = []
    p_vals = []
    s_tips = []
    h_tips = []
    r_vals = []
    r_distro_vals = []

    # store the child tip names for each node on the tree
    cache_tipnames(otu_tree)

    # iterate over the tree of child OTUs
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

            # calculate permutation p value for hommola test for this node
            if perm_type == 'hommola':
                r, p, r_distro = hommola_cospeciation(host_dm_sub[1], otu_dm_sub[1],
                                                          interaction_sub, permutations)
            elif perm_type == 'host':
                r, p, r_distro = hommola_cospeciation_host(host_dm_sub[1], otu_dm_sub[1],
                                                          interaction_sub, permutations)

            # append to results list
            p_vals.append(p)
            r_vals.append(r)
            r_distro_vals.append(r_distro)

        # If only testing top-level node, break out of tree traverse.
        if not recurse:
            break

    results_list = [p_vals, s_tips, h_tips, s_nodes, h_nodes, r_vals]

    results_header = ['p_vals', 's_tips', 'h_tips', 's_nodes', 'h_nodes', 'r_vals']

    return (results_list, results_header)

def make_dists_and_tree(sample_names, host_fp, host_input_type):
    """
    This routine reads in your host information (tree, alignment, or distance 
    matrix) and converts it to a distance matrix and a tree. These are subsetted
    to just the samples passed to the routine. Both the 
    distance matrix and host subtree are passed back to the main routine for 
    testing.
    """
    with open(host_fp, 'r') as host_f:
        host_str = host_f.read()

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

    dm_names_dict = {x:i for i, x in enumerate(dm[0])}

    name_slice = [dm_names_dict[name] for name in sample_names]

    sorted_dm = dm[1][numpy.ix_(name_slice, name_slice)]

    return (sample_names, sorted_dm)


def filter_dms(otu_dm, host_dm, interaction, otu_subset):
    """This filters a host dm, symbiont dm, and interaction matrix by a set of
    symbionts (otus) defined by otu_subset, and returns the sliced values.
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

def add_corrections_to_results_dict(results_dict, results_header):
    # Takes a results dictionary and adds multiple test corrections.

    pvals = []
    potus = []

    for potu_name in results_dict:
        potus = reduce(add,[potus,[potu_name]*len(results_dict[potu_name][0])])
        pvals = reduce(add,[pvals,results_dict[potu_name][results_header.index('p_vals')]])

    b_h_fdr_p_vals = benjamini_hochberg_step_down(pvals)
    bonferroni_p_vals = bonferroni_correction(pvals)
    fdr_p_vals = fdr_correction(pvals)

    for potu_name in results_dict:
        potu_b_h_fdr_p_vals = [b_h_fdr_p_vals[i] for i, x 
                               in enumerate(potus) if x == potu_name]
        results_dict[potu_name].append(potu_b_h_fdr_p_vals)

        potu_fdr_p_vals = [fdr_p_vals[i] for i, x 
                           in enumerate(potus) if x == potu_name]
        results_dict[potu_name].append(potu_fdr_p_vals)
        
        potu_bonferroni_p_vals = [bonferroni_p_vals[i] for i, x 
                                  in enumerate(potus) if x == potu_name]
        results_dict[potu_name].append(potu_bonferroni_p_vals)


    results_header = reduce(add,[results_header,['B&H_FDR_pvals','FDR_pvals','Bonferroni_pvals']])
        
    return(results_dict, results_header)


def add_h_span_to_results_dict(results_dict, results_header, host_tree):
    # Takes a results dictionary and adds multiple test corrections.

    for potu_name in results_dict:
        h_span_list =[]
        for i in range(len(results_dict[potu_name][results_header.index('h_nodes')])):
            h_span = str(len(host_tree.lowestCommonAncestor(results_dict[potu_name][results_header.index('h_nodes')][i].getTipNames()).getTipNames()))
            h_span_list.append(h_span)
        results_dict[potu_name].append(h_span_list)

    results_header.append('h_span')
        
    return(results_dict, results_header)


def write_cospeciation_results(results_dict, results_header, significance_level, output_dir, host_tree, otu_to_taxonomy, test):
    # Takes a results dictionary, which is a dictionary keyed by parent OTU 
    # names, with values composed of a list of lists. Each internal list 
    # corresponds to an output of the test, with values ordered according to
    # the order in which each node of the pOTU was tested. In the case of non-
    # recursive tests, these lists will only have a single value. 
    #
    # Also takes a results_header, which is a list of header names associated 
    # with these values, and an output directory path.
    # 
    # Does:
    # - Performs significance correction (FDR, B&H, Bonferroni)
    # - Calculates number of significant nodes under each correction
    #
    # Writes:
    # - a results file for each pOTU, indicating each node and its values
    # - a results summary file, indicating for each pOTU the number of nodes
    #   tested and the number significant under each correction factor.
    # - a file of significant nodes under each multiple test correction
    #
    results_dict, results_header = add_corrections_to_results_dict(results_dict, results_header)
    results_dict, results_header = add_h_span_to_results_dict(results_dict, results_header, host_tree)
   
    sig_nodes, fdr_sig_nodes, bh_fdr_sig_nodes, bonferroni_sig_nodes = get_sig_nodes(results_dict, results_header, significance_level)

    write_per_otu_results_file(results_dict, results_header, output_dir, test)
    write_summary_file(results_dict, results_header, output_dir, otu_to_taxonomy, sig_nodes, fdr_sig_nodes, bh_fdr_sig_nodes, bonferroni_sig_nodes)
    write_sig_nodes_files(results_dict, results_header, output_dir, otu_to_taxonomy, sig_nodes, fdr_sig_nodes, bh_fdr_sig_nodes, bonferroni_sig_nodes)

def get_sig_nodes(results_dict, results_header, significance_level):
    sig_nodes = []
    fdr_sig_nodes = []
    bh_fdr_sig_nodes = []
    bonferroni_sig_nodes = []
    
    potu_names = list(results_dict.keys())
    for potu in potu_names:
        num_nodes = len(results_dict[potu][0])
        for i in range(num_nodes):
            if results_dict[potu][results_header.index('p_vals')][i] < significance_level:
                sig_nodes.append((potu,i))
            if results_dict[potu][results_header.index('FDR_pvals')][i] < significance_level:
                fdr_sig_nodes.append((potu,i))
            if results_dict[potu][results_header.index('Bonferroni_pvals')][i] < significance_level:
                bonferroni_sig_nodes.append((potu,i))
            if results_dict[potu][results_header.index('B&H_FDR_pvals')][i] < significance_level:
                bh_fdr_sig_nodes.append((potu,i))
    return sig_nodes, fdr_sig_nodes, bh_fdr_sig_nodes, bonferroni_sig_nodes

def collapse_and_write_otu_table(otu_table_fp, mapping_fp, collapse_fields, collapse_mode):

    otu_table = load_table(otu_table_fp)

    sample_ids_to_keep = otu_table.ids()

    with open(mapping_fp, 'U') as mapping_f:
        sample_id_f_ids = set([l.strip().split()[0] for l in mapping_f if not
                               l.startswith('#')])

        sample_ids_to_keep = set(sample_ids_to_keep) & sample_id_f_ids

    filtered_otu_table = filter_samples_from_otu_table(
        otu_table, sample_ids_to_keep, 0, np.inf,
        negate_ids_to_keep=False)

    collapsed_metadata, collapsed_table = \
        collapse_samples(filtered_otu_table,
                         open(mapping_fp, 'U'),
                         collapse_fields,
                         collapse_mode)

    output_biom_fp = '_'.join([os.path.splitext(otu_table_fp)[0]] + 
                                collapse_fields) + os.path.splitext(otu_table_fp)[1]

    #print collapsed_table

    write_biom_table(collapsed_table, output_biom_fp, write_hdf5=False)

    return output_biom_fp

def write_per_otu_results_file(results_dict, results_header, output_dir, test):
    #Write per-OTU results files:
    potu_names = list(results_dict.keys())
    for potu in potu_names:
        results_file = open(os.path.join(output_dir,('%s_%s_results.txt' % (potu, test))), 'w')
        for header in results_header:
            results_file.write(header + "\t")
        results_file.write("\n")
        for i in range(len(results_dict[potu][0])):
            results_file.write("\t".join(str(x[i]) for x in results_dict[potu]))
            results_file.write("\n")
        results_file.close()

def write_summary_file(results_dict, results_header, output_dir, otu_to_taxonomy, sig_nodes, fdr_sig_nodes, bh_fdr_sig_nodes, bonferroni_sig_nodes):
    #Write results summary file
    potu_names = list(results_dict.keys())
    summary_file = open(os.path.join(output_dir, "cospeciation_results_summary.txt"),'w')
    summary_file.write("pOTU\tUncorrected sig nodes\tFDR sig nodes\tB&H FDR sig nodes\tBonferroni sig nodes\tNodes tested\tTaxonomy\n")
    for potu in potu_names:
        num_nodes = len(results_dict[potu][0])
        outline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(potu, 
            len([item for item in sig_nodes if item[0] == potu]), 
            len([item for item in fdr_sig_nodes if item[0] == potu]), 
            len([item for item in bh_fdr_sig_nodes if item[0] == potu]), 
            len([item for item in bonferroni_sig_nodes if item[0] == potu]), 
            num_nodes, otu_to_taxonomy[potu])
        summary_file.write(outline)
    summary_file.close()

def write_sig_nodes_files(results_dict, results_header, output_dir, otu_to_taxonomy, sig_nodes, fdr_sig_nodes, bh_fdr_sig_nodes, bonferroni_sig_nodes):
    #Write lists of significant nodes for each standard of significance
    for shortname, name, var in [('uncorrected', 'p_vals', sig_nodes), ('FDR', 
        'FDR_pvals', fdr_sig_nodes), ('bonferroni', 'Bonferroni_pvals', 
        bonferroni_sig_nodes), ('bh_FDR', 'B&H_FDR_pvals', bh_fdr_sig_nodes)]:
        
        sig_nodes_file = open(os.path.join(output_dir, "%s_sig_nodes.txt" % shortname),'w')
        sig_nodes_file.write("pOTU\ttaxonomy\t")
        sig_nodes_file.write("\t".join(results_header) + "\n")

        for sig_node in var:
            sig_nodes_file.write("{0}\t{1}\t".format(sig_node[0],otu_to_taxonomy[sig_node[0]]))
            sig_nodes_file.write("\t".join(str(x[sig_node[1]]) for x in results_dict[sig_node[0]]))
            sig_nodes_file.write("\n")
        sig_nodes_file.close()

def reconcile_hosts_symbionts(cotu_table, host_dist, min_val=1):

    shared_hosts = set(cotu_table.ids(axis='sample')).intersection(host_dist[0])

    # filter cOTU table by samples present in host_tree/dm

    cotu_table_filtered = cotu_table.filter(shared_hosts, axis='sample', inplace=False)
    

    # filter cOTU table again to get rid of absent cOTUs. BIOM should do this.

    cotu_table_filtered.filter(lambda val, id_, metadata: min_val <= val.sum(), axis='observation', inplace=True)
    
    # Filter the host_dists to match the newly trimmed subtree
    # Note: this is requiring the modified filter_dist method which
    # returns a native dm tuple rather than a string.

    host_dist_filtered = filter_samples_from_distance_matrix(
        host_dist, cotu_table_filtered.ids(axis = 'sample'), negate=True)

    

    return cotu_table_filtered, host_dist_filtered
