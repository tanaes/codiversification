#!/usr/bin/python

'''
Takes in a (txt) OTU table generated from 
reference OTU picking (database comparison), and makes
a copy of it with those OTUs that weren't found in the
removed (i.e. the ones whose names start with 'None').
'''
import sys
import os

def remove_unmatched_otus(otu_table_fp):
	if os.path.splitext(otu_table_fp)[1] != '.txt':
		print 'OTU table must be a text file. Use convert_biom.py.'

	output_fp = os.path.splitext(otu_table_fp)[0] + '_unmatched_otus_removed.txt'
	output = open(output_fp, 'w')

	for line in open(otu_table_fp, 'r'):
		if line[:4] != 'None':
			output.write(line)
	output.close()
	

if __name__ == '__main__':
	remove_unmatched_otus(sys.argv[1])