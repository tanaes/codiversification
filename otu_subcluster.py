#!/usr/bin/env python
# File created on 09 Aug 2012

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jon Sanders", "Nate Bresnick"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Development"

# importing modules
from sys import exit, stderr, stdout

import os
import sys

from qiime.util import make_option
from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup
from qiime.parse import parse_qiime_parameters, parse_taxonomy
import os.path
import os
from qiime.workflow.upstream import run_pick_de_novo_otus
from qiime.parse import parse_qiime_parameters, parse_taxonomy, parse_distmat, make_envs_dict, fields_to_dict
from qiime.workflow.util import call_commands_serially, no_status_updates
from cogent import LoadTree, LoadSeqs, DNA
from biom import load_table

options_lookup = get_options_lookup()

# otu_subcluster_v0.py
script_info = {}
script_info[
    'brief_description'] = """Subcluster sequences within an OTU table"""
script_info['script_description'] = """This workflow script takes the constituent 
sequences of OTUs from an OTU map ('parent OTUs') and re-clusters them at a higher
percent ID (into 'child OTUs') using typical QIIME workflow commands, outputting
new OTU tables, maps, and rep_set files (including an alignment and tree) for 
each OTU.

Changes to the default QIIME workflow behavior must be specified by a parameters 
file. If, for example, you want to subcluster OTU sequences at 99% ID, you would
 modify the pick_otus:similarity line to 0.99."""
script_info['script_usage'] = []


script_info['script_usage'].append(("", """Subcluster all sequences: """, """%prog 
-i $PWD/otu_map.txt -o $PWD/subclustered_otus -f $PWD/seqs.fna -p $PWD/subcluster_params.txt 
-t $PWD/Example_input/otu_table.txt"""))
script_info['output_description'] = """Output is standard QIIME otu tables, maps, 
and rep_set fasta files, one per original (parent) OTU. Files are all output to 
a new directory, specified by -o or --output_dir. These files follow the naming 
convention pOTU#_[fasta basename]_[filetype]. For example, the sequences which 
constitute parent OTU #90 will end up in a file called '90_seqs.fasta.' Clusters
 for these child OTUs will be defined in 90_seqs_otus.txt, and summarized in an 
 OTU table in 90_seqs_otu_table.txt."""

script_info['required_options'] = [

    options_lookup["otu_map_as_primary_input"],
    options_lookup["output_dir"],
    options_lookup["input_fasta"],
    make_option('-p', '--parameter_fp', type="existing_filepath",
                help='path to parameters file [REQUIRED]'),
    make_option('-b', '--biom_table_fp', type="existing_filepath",
                help='path to OTU table file [REQUIRED]')
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

    output_dir = opts.output_dir
    otu_map_fp = opts.otu_map_fp
    otu_table_fp = opts.biom_table_fp
    fasta_fp = opts.input_fasta_fp
    parameter_fp = opts.parameter_fp
    force = opts.force

    qiime_config = load_qiime_config()

    # Verify that parameters file exists, if it is specified
    with open(parameter_fp) as parameter_f:
        params = parse_qiime_parameters(parameter_f)

    try:
        os.makedirs(output_dir)
    except OSError:
        if force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            raise OSError("Output directory already exists. Please choose "
                "a different directory, or force overwrite with --force")

    # these are hardcoded from options selection in pick_otus_through_otu...py
    command_handler = call_commands_serially
    status_update_callback = no_status_updates
    parallel = False

    # get parent OTU map and load it into dict otu_to_seqid
    with open(otu_map_fp) as otu_map_f:
        otu_to_seqid = fields_to_dict(otu_map_f)

    # get the seqs.fna for all the sequences in the whole set
    fasta_collection = LoadSeqs(fasta_fp, moltype=DNA, aligned=False,
                                label_to_name=lambda x: x.split()[0])

    # get the pOTUs in the filtered otu table
    potu_table = load_table(otu_table_fp)

    # for each of the OTUs in this filtered parent table,
    for pOTU in potu_table.ids(axis='observation'):

        potu_dir = os.path.join(output_dir, str(pOTU))

        try:
            os.makedirs(potu_dir)
        except OSError:
            pass

        seqs_in_otu = fasta_collection.takeSeqs(otu_to_seqid[pOTU])
        output_fna_fp = os.path.join(potu_dir, 'seqs.fasta')
        seqs_in_otu.writeToFile(output_fna_fp)

        # pre process seqs from pOTUs
        run_pick_de_novo_otus(
            output_fna_fp,
            potu_dir,
            command_handler=command_handler,
            params=params,
            qiime_config=qiime_config,
            parallel=parallel,
            status_update_callback=status_update_callback,
            run_assign_tax=False)

if __name__ == "__main__":
    main()
