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
from test_test import *
''' NOT SURE WHY IT ONLY WORKS THIS WAY '''
#from test_cospeciation import *
from test_cospeciation import (hommola_cospeciation_test,
                               get_dist, cogent_dist_to_qiime_dist,recursive_hommola)





class HommolaTests(TestsHelper):
    # be careful with redefining the 'setUp' method, because it's
    # already defined in the superclass (TestsHelper), and you don't
    # want to fully override it here

    def test_hommola_cospeciation_test(self):
        hdist = array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
                      8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
        pdist = array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
                      8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
        matrix = array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
                       0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]])

        # this matrix was picked because it will generate an r value that's less than
        # a standard deviation away from the mean of the normal distribution of
        # r vals
        randomized_matrix = array(
            [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])

        self.assertCorrectPValue(0.0, 0.05, hommola_cospeciation_test,
                                 (hdist, pdist, matrix, 1000))

        self.assertCorrectPValue(0.2, 0.8, hommola_cospeciation_test,
                                 (hdist, pdist, randomized_matrix, 1000))

    def test_get_dist(self):
        labels = [0, 1, 1, 2, 3]
        dists = array([[0, 2, 6, 3], [2, 0, 5, 4], [6, 5, 0, 7], [3, 4, 7, 0]])
        index = [2, 3, 1, 0]

        expected_vec = [7, 7, 5, 6, 0, 4, 3, 4, 3, 2]
        actual_vec = get_dist(labels, dists, index)

        self.assertEqual(actual_vec, expected_vec)

    def test_recursive_hommola(self):
        
        aligned_otu_seqs = LoadSeqs(data=test_seqs)
        dm = (['SHNO', 'SHNP', 'SHNT', 'SHNW'], array([[ 0., 0.14369198, 0.17318318, 0.22670443],
       [ 0.14369198, 0., 0.17318318, 0.22670443],
       [ 0.17318318, 0.17318318, 0., 0.22670443],
       [ 0.22670443, 0.22670443, 0.22670443, 0.]]))
        host_subtree = LoadTree(treestring="((SHNT:0.0865915890705,(SHNP:0.0718459904766,SHNO:0.0718459904767):0.0147455985938):0.0267606238488,SHNW:0.113352212919);")
        otu_tree = LoadTree(treestring="(3:0.00499,(0:0.02491,(1:0.00015,2:0.02813)0.894:0.00792)0.655:0.00787,(4:0.00503,5:0.0025)0.927:0.00014);")
        host_dm = dm
        otu_table = Table.from_tsv(StringIO(otu_table_lines), None, None, lambda x : x)

        results_dict, acc_dict = recursive_hommola(aligned_otu_seqs,host_subtree,host_dm,otu_tree,
                        otu_table,permutations=1000,recurse=False)
        
        self.assertEqual(results_dict, expected_results)
        
        self.assertEqual(acc_dict,expected_acc)
        
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

# run tests if called from command line
if __name__ == "__main__":
    main()
