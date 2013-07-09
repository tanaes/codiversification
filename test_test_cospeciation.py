#!/usr/bin/env python
#file test_filter_otus_by_sample.py

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

# the following import added 2013-07-09 by Aaron Behr
# for method test_cogent_dist_to_qiime_dist 
from test_cospeciation import *


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
        """Tests cogent_dist_to_qiime_dist() in dist_convert.py.
        
        1)  An actual output Qiime distance matrix is generated,
            using the input data.
        2)  An expected output matrix is generated. The header row
            from the output matrix is used to populate the output
            matrix in proper order.
        3)  The actual and expected numpy array objects are compared.

        Note that the actual output qiime distance matrix in its entirety
        is not analyzed. This would be redundant since its header row is
        used to generate the expected output.
        """

        input_dict = {('a','b'): 4, ('a','c'): 5, ('a','d'): 6,
            ('b','a'): 4, ('b','c'): 7, ('b','d'): 8,
            ('c','a'): 5, ('c','b'): 7, ('c','d'): 9,
            ('d','a'): 6, ('d','b'): 8,('d','c'): 9}
        
        # Run actual method
        actual_output = cogent_dist_to_qiime_dist(input_dict)

        # Generate expected output
        matrix_order = actual_output[0] # get matrix order from actual output
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
        expected_output = numpy.array(expected_output) 

        self.assertEqual(actual_output[1], expected_output)

        

#run tests if called from command line
if __name__ == "__main__":
    main()
