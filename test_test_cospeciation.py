#!/usr/bin/env python
#file test_filter_otus_by_sample.py


# augmented & updated 2013-07-09 and 2014-05-01 by Aaron Behr
__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from numpy import array
from StringIO import StringIO
from os.path import exists
from cogent.util.unit_test import TestCase, main
from os import remove
from cogent import LoadSeqs
import shutil
from qiime.filter_otus_by_sample import (filter_otus,filter_aln_by_otus,\
                                 process_extract_samples)

# For the HommolaTests, need to get this from cogent in 
#tests/test_maths/test_stats (actual file test_test.py)
from test_test import TestsHelper

''' NOT SURE WHY IT ONLY WORKS THIS WAY '''
#from test_cospeciation import *
from test_cospeciation import (hommola_cospeciation_test,\
                         get_dist, cogent_dist_to_qiime_dist)


class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.aln=[]
        self.aln.append(('SampleA','AAAAAAAAAAAAAAA'))
        self.aln.append(('SampleB','CCCCCCC'))
        self.aln.append(('SampleC','GGGGGGGGGGGGGG'))
        
        self.otus={'0':['SampleC'],'1':['SampleC','SampleA'],'2':['SampleB',\
                                                                  'SampleA']}

        self.prefs={}
        self.prefs={'0':'SampleC','1':'SampleB'}

    def test_take_(self):
        """process_extract_samples: parses the cmd line and determines which
        samples should be removed"""

        self.sample_to_extract='SampleA,SampleB'
        exp1={'0':'SampleA','1':'SampleB'}

        obs1=process_extract_samples(self.sample_to_extract)

        self.assertEqual(obs1,exp1)

    def test_filter_otus(self):
        """filter_otus: determines which sequences should be removed and
        generates a new otus list"""

        exp1=[('1',['SampleA']),('2',['SampleA'])]
        obs1=filter_otus(self.otus,self.prefs)
            
        self.assertEqual(obs1,exp1)

        
    def test_filter_aln_by_otus(self):
        """filter_aln_by_otus: determines which sequences to keep and which
        sequences to remove"""

        self.sample_to_extract='SampleA,SampleB'
        exp1=[]
        exp1.append(('SampleA','AAAAAAAAAAAAAAA'))
        exp2=[]
        exp2.append(('SampleB','CCCCCCC'))
        exp2.append(('SampleC','GGGGGGGGGGGGGG'))
        aln=LoadSeqs(data=self.aln,aligned=False)

        obs1,obs2=filter_aln_by_otus(aln,self.prefs)

        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)


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



class HommolaTests(TestsHelper):
    # be careful with redefining the 'setUp' method, because it's
    # already defined in the superclass (TestsHelper), and you don't
    # want to fully override it here

    def test_hommola_cospecitation_test(self):
        hdist = array([[0,3,8,8,9],[3,0,7,7,8],[8,7,0,6,7],[8,7,6,0,3],[9,8,7,3,0]])
        pdist = array([[0,5,8,8,8],[5,0,7,7,7],[8,7,0,4,4],[8,7,4,0,2],[8,7,4,2,0]])
        matrix = array([[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,1,1]])
        
        # this matrix was picked because it will generate an r value that's less than 
        # a standard deviation away from the mean of the normal distribution of r vals 
        randomized_matrix = array([[0,0,0,1,0],[0,0,0,0,1],[1,0,0,0,0],[1,0,0,0,0],[0,0,0,0,1]])

        self.assertCorrectPValue(0.0, 0.05, hommola_cospeciation_test, 
                                                (hdist, pdist, matrix, 1000))
        
        self.assertCorrectPValue(0.2, 0.8, hommola_cospeciation_test, 
                                                (hdist, pdist, randomized_matrix, 1000))


    def test_get_dist(self):
        labels = [0,1,1,2,3]
        dists = array([[0,2,6,3],[2,0,5,4],[6,5,0,7],[3,4,7,0]])
        index = [2,3,1,0]

        expected_vec = [7,7,5,6,0,4,3,4,3,2]
        actual_vec = get_dist(labels,dists,index)

        self.assertEqual(actual_vec, expected_vec)


    def test_distmat_to_tree(self):
        pass
        






#run tests if called from command line
if __name__ == "__main__":
    main()
