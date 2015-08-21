#!/usr/bin/env python
# file test_filter_otus_by_sample.py


# augmented & updated 2013-07-09 and 2014-05-01 by Aaron Behr
__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["Jesse Stombaugh"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"


import os
import sys
import re
from StringIO import StringIO
import numpy
import numpy as np
import numpy.testing as npt

from skbio.stats.distance import mantel

from random import shuffle
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_distmat, make_envs_dict
from qiime.filter import filter_samples_from_otu_table, filter_samples_from_distance_matrix

from biom import load_table, Table

from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from cogent.phylo import distance, nj
from cogent.evolve.models import HKY85
from cogent.evolve.pairwise_distance import TN93Pair
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent import LoadTree, LoadSeqs, DNA
from cogent.util.dict2d import Dict2D, largest
# For the HommolaTests, need to get this from cogent in
# tests/test_maths/test_stats (actual file test_test.py)
from cogent.util.unit_test import TestCase, main
from cogent.util.array import array
from cospeciation import *
from cospeciation import _get_dist, _gen_lists

import unittest

class HommolaCospeciationTests(unittest.TestCase):
    def setUp(self):
        # Test matrices, as presented in original paper by Hommola et al.
        self.hdist = np.array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
            8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
        self.pdist = np.array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
            8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
        self.interact = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
            0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]])

        # Reduced-size host matrix for testing asymmetric interaction matrix
        self.hdist_4x4 = np.array([[0, 3, 8, 8], [3, 0, 7, 7], [8, 7, 0, 6],
                                  [8, 7, 6, 0]])
        self.interact_5x4 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                                      [0, 0, 0, 1], [0, 0, 0, 1]])

        # One to one interaction matrix for comparing against Mantel output
        self.interact_1to1 = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
            0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])

        # interaction matrix yielding non-significant results.
        # this matrix was picked because it will generate an r value that's
        # less than a standard deviation away from the mean of the normal
        # distribution of r vals
        self.interact_ns = np.array(
            [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 0],
             [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])

        # minimal size matrices for sanity checks of inputs
        self.h_dist_3x3 = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        self.h_dist_2x2 = np.array([[0, 3], [3, 0]])
        self.p_dist_3x3 = np.array([[0, 3, 2], [3, 0, 1], [2, 1, 0]])
        self.interact_3x3 = np.array([[0, 1, 1], [1, 0, 1], [0, 0, 1]])
        self.interact_3x2 = np.array([[0, 1], [1, 0], [1, 1]])
        self.interact_2x3 = np.array([[0, 1, 1], [1, 0, 1]])
        self.interact_zero = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    def test_hommola_cospeciation_sig(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact, 9)
        exp_p = .1
        exp_r = 0.83170965463247915
        exp_perm_stats = np.array([-0.14928122, 0.26299538, -0.21125858,
                                   0.24143838, 0.61557855, -0.24258293,
                                   0.09885203, 0.02858, 0.42742399])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_cospeciation_asymmetric(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist_4x4, self.pdist, self.interact_5x4, 9)
        exp_p = 0.2
        exp_r = 0.85732140997411233
        exp_perm_stats = np.array([-0.315244162496, -0.039405520312,
                                   0.093429386594, -0.387835875941,
                                   0.183711730709,  0.056057631956,
                                   0.945732487487,  0.056057631956,
                                   -0.020412414523])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_cospeciation_host_sig(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation_host(
            self.hdist, self.pdist, self.interact, 9)
        exp_p = .1
        exp_r = 0.83170965463247915
        exp_perm_stats = np.array([ 0.06451683,  0.31988833,  0.00816571,
                                0.07187642,  0.07187642, 0.46916955,
                                -0.0678993 ,  0.57579899,  0.31895164])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats, rtol=1e-4)

    def test_hommola_cospeciation_host_asymmetric(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation_host(
            self.hdist_4x4, self.pdist, self.interact_5x4, 9)
        exp_p = 0.6
        exp_r = 0.85732140997411233
        exp_perm_stats = np.array([0.945732,  0.428661,  0.945732,  0.945732,
                                   0.945732,  0.35465 , 0.428661,  0.168173,
                                   0.95298 ])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats, rtol=1e-4)

    def test_hommola_cospeciation_no_sig(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact_ns, 9)
        exp_p = .6
        exp_r = -0.013679391379114569
        exp_perm_stats = np.array([-0.22216543, -0.14836061, -0.04434843,
                                   0.1478281, -0.29105645, 0.56395839,
                                   0.47304992, 0.79125657, 0.06804138])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)
        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_vs_mantel(self):
        # we don't compare p-values because the two methods use different
        # permutation strategies
        r_mantel, p_mantel, _ = mantel(
            self.hdist, self.pdist, method='pearson', permutations=0,
            alternative='greater'
        )
        r_hommola, p_hommola, _ = hommola_cospeciation(
            self.hdist, self.pdist, self.interact_1to1, permutations=0
        )

        self.assertAlmostEqual(r_hommola, r_mantel)
        npt.assert_equal(p_hommola, p_mantel)

    def test_zero_permutations(self):
        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact, 0)

        exp_p = np.nan
        exp_r = 0.83170965463247915
        exp_perm_stats = np.array([])

        npt.assert_equal(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)
        npt.assert_equal(obs_perm_stats, exp_perm_stats)

    def test_get_dist(self):
        labels = np.array([0, 1, 1, 2, 3])
        k_labels, t_labels = _gen_lists(labels)
        dists = np.array([[0, 2, 6, 3], [2, 0, 5, 4], [6, 5, 0, 7],
                          [3, 4, 7, 0]])
        index = np.array([2, 3, 1, 0])

        expected_vec = np.array([7, 7, 5, 6, 0, 4, 3, 4, 3, 2])
        actual_vec = _get_dist(k_labels, t_labels, dists, index)

        npt.assert_allclose(actual_vec, expected_vec)

    def test_gen_lists(self):
        exp_pars_k_labels = np.array([0, 0, 0, 0, 0, 1, 1, 1,
                                      1, 2, 2, 2, 3, 3, 4])
        exp_pars_t_labels = np.array([1, 2, 3, 4, 4, 2, 3, 4,
                                      4, 3, 4, 4, 4, 4, 4])
        exp_host_k_labels = np.array([0, 0, 0, 0, 0, 1, 1, 1,
                                      1, 2, 2, 2, 3, 3, 3])
        exp_host_t_labels = np.array([1, 2, 3, 3, 4, 2, 3, 3,
                                      4, 3, 3, 4, 3, 4, 4])

        pars, hosts = np.nonzero(self.interact)

        obs_pars_k_labels, obs_pars_t_labels = _gen_lists(pars)
        obs_hosts_k_labels, obs_hosts_t_labels = _gen_lists(hosts)

        npt.assert_allclose(exp_pars_k_labels, obs_pars_k_labels)
        npt.assert_allclose(exp_pars_t_labels, obs_pars_t_labels)
        npt.assert_allclose(exp_host_k_labels, obs_hosts_k_labels)
        npt.assert_allclose(exp_host_t_labels, obs_hosts_t_labels)

    def test_dm_too_small(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_2x2, self.p_dist_3x3,
                                 self.interact_3x3)

    def test_host_interaction_not_equal(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_2x3)

    def test_par_interaction_not_equal(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_3x2)

    def test_interaction_too_few(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_zero)

    def test_permutations_too_few(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_3x3, -1)


class HommolaTests(TestCase):
    # be careful with redefining the 'setUp' method, because it's
    # already defined in the superclass (TestsHelper), and you don't
    # want to fully override it here
    def test_cogent_dist_to_qiime_dist(self):

        input_dict = {('a','b'): 4, ('a','c'): 5, ('a','d'): 6,
            ('b','a'): 4, ('b','c'): 7, ('b','d'): 8,
            ('c','a'): 5, ('c','b'): 7, ('c','d'): 9,
            ('d','a'): 6, ('d','b'): 8,('d','c'): 9}
        
        # RUN METHOD TO BE TESTED
        actual_output = cogent_dist_to_qiime_dist(input_dict)

        # Generate expected output
        matrix_order = ['a','b','c','d']
        expected_output = []
        for x in matrix_order:
            row = []
            for y in matrix_order:
                if x != y:
                    # use input tuple dist dict to populate expected output matrix
                    row.append(input_dict[x,y])
                else:
                    # input tuple dist dict doesn't store the null self-distances
                    row.append(0.)  
            expected_output.append(row)
        expected_output = (matrix_order, array(expected_output))

        self.assertEqual(actual_output, expected_output)

    def test_add_corrections_to_results(self):
        h_nodes = [LoadTree(treestring="(((((SHNT:0.0824904475935,(SHOG:0.0687626692246,SHOB:0.0687626692247):0.0137277783688):0.00410114147696,(SHNP:0.0718459904766,(SHNO:0.0482357525131,SHOF:0.0482357525131):0.0236102379636):0.0147455985938):0.00675384740989,SHNY:0.0933454364804):0.00997699844309,SHNI:0.103322434923):0.0100297779958,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="(((((SHNT:0.0824904475935,(SHOG:0.0687626692246,SHOB:0.0687626692247):0.0137277783688):0.00410114147696,(SHNO:0.0482357525131,SHOF:0.0482357525131):0.0383558365574):0.00675384740989,SHNY:0.0933454364804):0.00997699844309,SHNI:0.103322434923):0.0100297779958,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="((((SHOG:0.0687626692246,SHOB:0.0687626692247):0.0178289198458,SHNO:0.0865915890705):0.00675384740989,SHNY:0.0933454364804):0.0200067764389,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="((SHOG:0.0933454364803,SHNY:0.0933454364804):0.0200067764389,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="(SHNY:0.113352212919,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="((((SHNT:0.0824904475935,SHOB:0.0824904475935):0.00410114147696,SHOF:0.0865915890705):0.016730845853,SHNI:0.103322434923):0.0100297779958,SHNU:0.113352212919);")]

        p_vals = [0.295704295704,0.271728271728,0.444555444555,0.652347652348,0.784215784216,0.355644355644]

        s_nodes = [LoadTree(treestring="(1:0.038895,(((6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014,(2:0.02835,3:0.00239)0.855:0.01023)1.000:0.00018,((5:0.00014,7:0.00014)0.924:0.00499,(0:0.01261,9:0.0025)0.550:0.00015)1.000.2:0.00483):0.000135);"),LoadTree(treestring="(((6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014,(2:0.02835,3:0.00239)0.855:0.01023)1.000:0.00018,((5:0.00014,7:0.00014)0.924:0.00499,(0:0.01261,9:0.0025)0.550:0.00015)1.000.2:0.00483):0.000135;"),LoadTree(treestring="((6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014,(2:0.02835,3:0.00239)0.855:0.01023)1.000:0.00018;"),LoadTree(treestring="(6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014;"),LoadTree(treestring="(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689;"),LoadTree(treestring="((5:0.00014,7:0.00014)0.924:0.00499,(0:0.01261,9:0.0025)0.550:0.00015)1.000.2:0.00483;")]

        s_tips = [11,10,6,4,3,4]

        h_tips = [10,9,6,4,3,5]

        r_vals = [0.0174371319933,0.0083461905486,-0.0258459108602,-0.271865301768,-0.411161344273,-0.115766076854]

        h_span = [15,15,15,15,15,15]

        h_nodes2 = [LoadTree(treestring="(((((SHNT:0.0824904475935,(SHOG:0.0687626692246,SHOB:0.0687626692247):0.0137277783688):0.00410114147696,(SHNP:0.0718459904766,(SHNO:0.0482357525131,SHOF:0.0482357525131):0.0236102379636):0.0147455985938):0.00675384740989,SHNY:0.0933454364804):0.00997699844309,SHNI:0.103322434923):0.0100297779958,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="(((((SHNT:0.0824904475935,(SHOG:0.0687626692246,SHOB:0.0687626692247):0.0137277783688):0.00410114147696,(SHNO:0.0482357525131,SHOF:0.0482357525131):0.0383558365574):0.00675384740989,SHNY:0.0933454364804):0.00997699844309,SHNI:0.103322434923):0.0100297779958,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="((((SHOG:0.0687626692246,SHOB:0.0687626692247):0.0178289198458,SHNO:0.0865915890705):0.00675384740989,SHNY:0.0933454364804):0.0200067764389,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="((SHOG:0.0933454364803,SHNY:0.0933454364804):0.0200067764389,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="(SHNY:0.113352212919,(SHNU:0.0987231206604,SHNW:0.0987231206604):0.0146290922588);"),LoadTree(treestring="((((SHNT:0.0824904475935,SHOB:0.0824904475935):0.00410114147696,SHOF:0.0865915890705):0.016730845853,SHNI:0.103322434923):0.0100297779958,SHNU:0.113352212919);")]

        p_vals2 = [0.0001,0.002,0.01,0.05,0.1,0.353]

        s_nodes2 = [LoadTree(treestring="(1:0.038895,(((6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014,(2:0.02835,3:0.00239)0.855:0.01023)1.000:0.00018,((5:0.00014,7:0.00014)0.924:0.00499,(0:0.01261,9:0.0025)0.550:0.00015)1.000.2:0.00483):0.000135);"),LoadTree(treestring="(((6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014,(2:0.02835,3:0.00239)0.855:0.01023)1.000:0.00018,((5:0.00014,7:0.00014)0.924:0.00499,(0:0.01261,9:0.0025)0.550:0.00015)1.000.2:0.00483):0.000135;"),LoadTree(treestring="((6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014,(2:0.02835,3:0.00239)0.855:0.01023)1.000:0.00018;"),LoadTree(treestring="(6:0.0025,(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689)0.890:0.00014;"),LoadTree(treestring="(8:0.00317,(10:0.00943,4:0.00318)0.679:0.00759)0.873:0.00689;"),LoadTree(treestring="((5:0.00014,7:0.00014)0.924:0.00499,(0:0.01261,9:0.0025)0.550:0.00015)1.000.2:0.00483;")]

        s_tips2 = [11,10,6,4,3,4]

        h_tips2 = [10,9,6,4,3,5]

        r_vals2 = [0.0174371319933,0.0083461905486,-0.0258459108602,-0.271865301768,-0.411161344273,-0.115766076854]

        h_span2 = [15,15,15,15,15,15]


        results_list = [p_vals, s_tips, h_tips, s_nodes, h_nodes, r_vals]
        results_list2 = [p_vals2, s_tips2, h_tips2, s_nodes2, h_nodes2, r_vals2]

        results_dict = {'potu1': results_list, 'potu2': results_list2}

        results_header = ['p_vals', 's_tips', 'h_tips', 's_nodes', 'h_nodes', 'r_vals']

        exp_b_h_fdr_p_vals = [0.47419247419200006,   0.47419247419200006,   0.53346653346600004,   0.71165198437963639,   0.78421578421599991,   0.47419247419200006]
        exp_fdr_p_vals = [0.50692164977828569,   0.54345654345600003,   0.53346653346600004,   0.71165198437963639,   0.78421578421599991,   0.47419247419200006]
        exp_bonferroni_p_vals = [3.5484515484479999,   3.260739260736,   5.3346653346600004,   7.828171828176,   9.4105894105919994,   4.2677322677280003]

        exp_b_h_fdr_p_vals_2 = [0.0012000000000000001,   0.012,   0.040000000000000001,   0.15000000000000002,   0.24000000000000005,   0.47419247419200006]
        exp_fdr_p_vals_2 = [0.0012000000000000001,   0.012,   0.040000000000000001,   0.15000000000000002,   0.24000000000000005,   0.52949999999999997]
        exp_bonferroni_p_vals_2 = [0.0012000000000000001,   0.024,   0.12,   0.60000000000000009,   1.2000000000000002,   4.2359999999999998]

        exp_results_list = [p_vals, s_tips, h_tips, s_nodes, h_nodes, r_vals, exp_b_h_fdr_p_vals, exp_fdr_p_vals, exp_bonferroni_p_vals]
        exp_results_list2 = [p_vals2, s_tips2, h_tips2, s_nodes2, h_nodes2, r_vals2, exp_b_h_fdr_p_vals_2, exp_fdr_p_vals_2, exp_bonferroni_p_vals_2]

        exp_results_dict = {'potu1': exp_results_list, 'potu2': exp_results_list2}

        exp_results_header = ['p_vals', 's_tips', 'h_tips', 's_nodes', 'h_nodes', 'r_vals','B&H_FDR_pvals','FDR_pvals','Bonferroni_pvals']

        obs_results_dict, obs_results_header = add_corrections_to_results_dict(results_dict, results_header)

        self.assertItemsEqual(exp_results_header,obs_results_header)


    def test_filter_dms(self):

        # Test that filtering with all otu IDs returns original DMs
        otu_subset = otu_dm[0]

        expected_otu_dm = otu_dm
        expected_host_dm = host_dm
        expected_interaction = interaction

        sliced_otu_dm, sliced_host_dm, sliced_interaction = filter_dms(otu_dm, host_dm, interaction, otu_subset)

        self.assertEqual(sliced_otu_dm,expected_otu_dm)
        self.assertEqual(sliced_host_dm,expected_host_dm)
        self.assertEqual(sliced_interaction,expected_interaction)

        # Test that filtering with a subset returns the correct subset

        otu_subset = ['1','3','5']

        expected_interaction = array([[1, 1, 0],
                                    [0, 1, 1],
                                    [0, 1, 0]])

        expected_otu_dm = (['1','3','5'], array([[0., 0.02030519, 0.01781229],
                                [0.02030519, 0., 0.0075794],
                                [0.01781229, 0.0075794, 0.]]))

        expected_host_dm = (['SHNO', 'SHNP', 'SHNT'], array([[0., 0.14369198, 0.17318318],
                            [0.14369198, 0., 0.17318318],
                            [0.17318318, 0.17318318, 0.]]))
        
        sliced_otu_dm, sliced_host_dm, sliced_interaction = filter_dms(otu_dm, host_dm, interaction, otu_subset)

        self.assertEqual(sliced_otu_dm,expected_otu_dm)
        self.assertEqual(sliced_host_dm,expected_host_dm)
        self.assertEqual(sliced_interaction,expected_interaction)
    
    def test_sort_dm_by_sample(self):
        unsorted_dm = (['VXNO', 'VXNP', 'VXNT'], array([[0., 0.704, 0.157],
                            [0.704, 0., 0.339],
                            [0.157, 0.339, 0.]]))

        sample_order = ['VXNP', 'VXNO', 'VXNT']

        expected_dm = (['VXNP', 'VXNO', 'VXNT'], array([[0., 0.704, 0.339],
                            [0.704, 0., 0.157],
                            [0.339, 0.157, 0.]]))
        
        dm = sort_dm_by_sample(unsorted_dm, sample_order)
        
        self.assertEqual(dm, expected_dm)
	
    def test_recursive_hommola(self):
        exp_h_nodes = LoadTree(treestring="((SHNT:0.0865915890705,(SHNP:0.0718459904766,SHNO:0.0718459904767):0.0147455985938):0.0267606238488,SHNW:0.113352212919);")
        exp_s_nodes = LoadTree(treestring="((1:0.00015,2:0.02813)0.894:0.00235,(0:0.02491,(3:0.00499,(4:0.00503,5:0.0025)0.927:0.00014)0.655:0.00787):0.00557);")
        exp_s_tips = 6
        exp_h_tips = 4
        exp_r_vals = 0.49739280928

        aligned_otu_seqs = LoadSeqs(data=test_seqs)

        host_subtree = LoadTree(treestring="((SHNT:0.0865915890705,(SHNP:0.0718459904766,SHNO:0.0718459904767):0.0147455985938):0.0267606238488,SHNW:0.113352212919);")
        otu_tree = LoadTree(treestring="(3:0.00499,(0:0.02491,(1:0.00015,2:0.02813)0.894:0.00792)0.655:0.00787,(4:0.00503,5:0.0025)0.927:0.00014);")
        otu_table = Table.from_tsv(StringIO(otu_table_lines), None, None, lambda x : x)

        results_list, results_header = recursive_hommola(aligned_otu_seqs,host_subtree,host_dm,otu_tree,
                        otu_table,permutations=1000,recurse=False)
        
        self.assertEqual(results_list[results_header.index('s_tips')][0],exp_s_tips)
        self.assertEqual(results_list[results_header.index('h_tips')][0],exp_h_tips)
        self.assertEqual(str(results_list[results_header.index('s_nodes')][0].rootAtMidpoint()),str(exp_s_nodes.rootAtMidpoint()))
        self.assertEqual(str(results_list[results_header.index('h_nodes')][0].rootAtMidpoint()),str(exp_h_nodes.rootAtMidpoint()))
        self.assertAlmostEqual(results_list[results_header.index('r_vals')][0],exp_r_vals)
        
    def test_distmat_to_tree(self):
        pass

otu_table_lines = "#OTU ID\tSHNO\tSHNP\tSHNT\tSHNW\n0\t1\t0\t0\t0\n1\t35\t1\t0\t0\n2\t0\t0\t0\t1\n3\t0\t1\t1\t0\n4\t0\t1\t0\t0\n5\t0\t4\t0\t0"
test_seqs = """>2
--TTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGGGAAAGTAGCTTGCTACTGACCTTAGCGGCGGACGGGTGAGTAATACTTAGGAATCTACCTATTAATGGGGGACAACGTTTCGAAAGGGACGCTAATACCGCATACGCCCTACGGGGGAAAGCAGGGGATCTTCGGACCTTGCGTTAATAGATGAGCCTAAGCCGGATTAGCTAGTTGGTGGGGTAAAGGCCTACCAAGGCGACGATCTGTAGCGGGTTTGAGAGGATGATCCGCCACACT-GGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAAGC
>1
GATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGGGAAGGTAGCTTGCTACTGGACCTAGCGGCGGACGGGTGAGTAATACTTAGGAATCTGCCTATTAGTGGGGGACAACGTTCCGAAAGGAGCGCTAATACCGCATACGCCCTACGGGGGAAAGCAGGGGATCTTCGGACCTTGCGCTAATAGATGAGCCTAAGTCGGATTAGCTAGTTGGTGGGGTAAAGGCCTACCAAGGCGACGATCTGTAGCGGGTTTGAGAGGATGATCCGCCACACT-GGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAA--
>0
GATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGGAAGGGTAGCTTGCTACCTGACCTAGCGGCGGACGGGTGAGTAATGCTTAGGAATCTGCCTATTAGTGGGGGACAACATTCCGAAAGGAATGCTAATACCGCATACGCCCTACGGGGGAAAGCAGGGGATCTTCGGACCTTGCGCTAATAGATGAGCCTAAGTCAGATTAGCTAGTTGGTGGGGTAAAGGCCTACCAAGGCGACGATCTGTAGCGGGTCTGAGAGGATGATCCGCCACACT-GGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTTTGGTTGTAAA--
>3
GATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGGAATTGTAGCTTGCTACATGACCTAGCGGCGGACGGGTGAGTAATACTTAGGAATCTGCCTATTAGTGGGGGACAACGTTCCGAAAGGAGCGCTAATACCGCATACGCCCTACGGGGGAAAGCAGGGGATCTTCGGACCTTGCGCTAATAGATGAGCCTAAGTCGGATTAGCTAGTTGGTAGGGTAAAGGCCTACCAAGGCGACGATCTGTGGCGGGTTTGAGAGGATGATCCGCCACACT-GGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAA--
>5
--TTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGGAATTGTAGCTTGCTACATGACCTAGCGGCGGACGGGTGAGTAATACTTAGGAATCTGCCTATTAGTGTGGGACAACGTTCCGAAAGGAGCGCTAATACCGCATACGCCCTACGGGGGAAAGCAGGGGATCTTCGGACCTTGCGCTAATAGATGAGCCTAAGTCGGATTAGCTAGTTGGTAGGGTAAAGGCCTACCAAGGCGACGATCTGTAGCGGGTTTGAGAGGATGATCCGCCACACT-GGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAAAGC
>4
GATTGAACGCTGGCGGCAGGCTTAACACATGCAAGTCGAGCGGGGAATTGTAGCTTGCTACATGACCTAGCGGCGGACGGGTGAGTAATACTTAGGAATCTGCCTATTAGTGGGGGACAACGTTCCGAAAGGAGCGCTAATACCGCATACGCCCTACGGGGGAAAGCAGGGGATCTTCGGACCTTGCGCTAATAGATGAGCCTAAGTCGGATTAGCTAGTTGGTAGGGTAAAGGCCTACCAAGGCGACGATCTGTAGCGGGTTTGAGAGGATGATCCGCCACACTGGGGGGTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTATGGTTGTAA---
"""

host_dm = (['SHNO', 'SHNP', 'SHNT', 'SHNW'], array([[ 0., 0.14369198, 0.17318318, 0.22670443],
       [ 0.14369198, 0., 0.17318318, 0.22670443],
       [ 0.17318318, 0.17318318, 0., 0.22670443],
       [ 0.22670443, 0.22670443, 0.22670443, 0.]]))

otu_dm = (['0', '1', '2', '3', '4', '5'], array([[ 0.        ,  0.03076205,  0.05527554,  0.03074287,  0.0360373 ,
         0.03345466],
       [ 0.03076205,  0.        ,  0.02837217,  0.02030519,  0.02035744,
         0.01781229],
       [ 0.05527554,  0.02837217,  0.        ,  0.04687171,  0.04699685,
         0.0439181 ],
       [ 0.03074287,  0.02030519,  0.04687171,  0.        ,  0.0101108 ,
         0.0075794 ],
       [ 0.0360373 ,  0.02035744,  0.04699685,  0.0101108 ,  0.        ,
         0.00759866],
       [ 0.03345466,  0.01781229,  0.0439181 ,  0.0075794 ,  0.00759866,
         0.        ]]))

interaction = array([[1, 0, 0, 0], [1, 1, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0], [0, 1, 0, 0], [0, 1, 0, 0]])

# run tests if called from command line
if __name__ == "__main__":
    main()
