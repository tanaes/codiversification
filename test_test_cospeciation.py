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
from cospeciation import *



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

    def test_recursive_hommola(self):
        exp_h_nodes = LoadTree(treestring="((SHNT:0.0865915890705,(SHNP:0.0718459904766,SHNO:0.0718459904767):0.0147455985938):0.0267606238488,SHNW:0.113352212919);")
        exp_s_nodes = LoadTree(treestring="((1:0.00015,2:0.02813)0.894:0.00235,(0:0.02491,(3:0.00499,(4:0.00503,5:0.0025)0.927:0.00014)0.655:0.00787):0.00557);")
        exp_s_tips = [6]
        exp_h_tips = [4]
        exp_r_vals = 0.49739280928

        aligned_otu_seqs = LoadSeqs(data=test_seqs)

        host_subtree = LoadTree(treestring="((SHNT:0.0865915890705,(SHNP:0.0718459904766,SHNO:0.0718459904767):0.0147455985938):0.0267606238488,SHNW:0.113352212919);")
        otu_tree = LoadTree(treestring="(3:0.00499,(0:0.02491,(1:0.00015,2:0.02813)0.894:0.00792)0.655:0.00787,(4:0.00503,5:0.0025)0.927:0.00014);")
        otu_table = Table.from_tsv(StringIO(otu_table_lines), None, None, lambda x : x)

        results_dict, acc_dict = recursive_hommola(aligned_otu_seqs,host_subtree,host_dm,otu_tree,
                        otu_table,permutations=1000,recurse=False)
        self.assertEqual(results_dict['s_tips'],exp_s_tips)
        self.assertEqual(results_dict['h_tips'],exp_h_tips)
        self.assertEqual(str(results_dict['s_nodes'][0].rootAtMidpoint()),str(exp_s_nodes.rootAtMidpoint()))
        self.assertEqual(str(results_dict['h_nodes'][0].rootAtMidpoint()),str(exp_h_nodes.rootAtMidpoint()))
        self.assertAlmostEqual(acc_dict['r_vals'][0],exp_r_vals)

        
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
