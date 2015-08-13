#!/usr/bin/env python

from __future__ import division

__author__ = "Nate Bresnick"
__copyright__ = "Copyright 2015, The QIIME Project"
__credits__ = ["Jon Sanders", "Nate Bresnick", "Aaron Behr"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Development"

import os
import sys
from qiime.util import make_option
from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup, write_biom_table
from qiime.parse import parse_qiime_parameters, parse_taxonomy
import os.path
import os
import csv

options_lookup = get_options_lookup()
script_info = {}
script_info['brief_description'] = """
A script to visualize the results of the test cospeciation script"""
script_info['script_description'] = """
.
"""
script_info['script_usage'] = []
script_info['output_description'] = """
"""
script_info['required_options'] = [
    make_option('-i', '--results_dir', type="existing_dirpath",
                help='the input folder, which contains results directories. IMPORTANT: this should be the grandparent folder of any hommola_results folders [REQUIRED]'),
    make_option('-f', '--folder_name', type="string",
                help='the name of the results folder inside of the directories in results_dir i.e. hommola_test_HostSpecies'),
    options_lookup["output_dir"]
    ]
script_info['optional_options'] = [
    make_option('--force', action='store_true',
                dest='force', help='Force overwrite of existing output directory'
                ' (note: existing files in output_dir will not be removed)'
                ' [default: %default]')
    ]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    results_dir = opts.results_dir
    output_dir = opts.output_dir    
    force = opts.force
    
    output_dir = os.path.abspath(output_dir)
    results_dir = os.path.abspath(results_dir)
    folder_name = opts.folder_name
    
    try:
        os.makedirs(output_dir)
    except OSError:
        if not force:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            raise OSError("Output directory already exists. Please choose "
                "a different directory, or force overwrite with -f.")
    cluster_widths = sorted([name for name in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir,name))])
    
    with open(os.path.join(output_dir, 'summary.txt'), 'w') as outfile:
        fieldnames = ["pOTU Width", "uncorrected_sig_nodes", "FDR_sig_nodes", "bh_FDR_sig_nodes", "bonferroni_sig_nodes"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for width in cluster_widths:
            rdict = {results_file: len(read_results(os.path.join(results_dir, width, folder_name, results_file + ".txt"))) for results_file in fieldnames[1:]}
            rdict['pOTU Width'] = width
            writer.writerow(rdict)
        
    
def read_results(fname):
    with open(fname, mode='r') as infile:
        reader = csv.DictReader(infile, delimiter="\t", )
        results_dicts = [row for row in reader]
        return results_dicts

if __name__ == "__main__":
    main()
