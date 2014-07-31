#!/n/sw/python-2.7.1/bin/python
# File created on 09 Aug 2012

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2012"
__credits__ = ["Rob Knight", "Justin Kuczynski",
               "Jesse Stombaugh", "Jon Sanders"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Experimental"

# Nate's test edit
# Jon's test addition

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
-i otu_map.txt -o subclustered_otus -f seqs.fna -p subcluster_params.txt 
-t Example_input/otu_table.txt"""))
script_info['output_description'] = """Output is standard QIIME otu tables, maps, 
and rep_set fasta files, one per original (parent) OTU. Files are all output to 
a new directory, specified by -o or --output_dir. These files follow the naming 
convention pOTU#_[fasta basename]_[filetype]. For example, the sequences which 
constitute parent OTU #90 will end up in a file called '90_seqs.fasta.' Clusters
 for these child OTUs will be defined in 90_seqs_otus.txt, and summarized in an 
 OTU table in 90_seqs_otu_table.txt."""

script_info['required_options'] = [
    make_option('-i', '--otu_map_fp',
                help='the input OTU map file [REQUIRED]'),
    make_option('-o', '--output_dir',
                help='path to output directory [REQUIRED]'),
    make_option('-f', '--fasta_fp',
                help='path to sequence fasta [REQUIRED]'),
    make_option('-p', '--parameter_fp',
                help='path to parameters file [REQUIRED]'),
    make_option('-t', '--otu_table_fp',
                help='path to OTU table file [REQUIRED]')
]
script_info['optional_options'] = [

    make_option('--force', action='store_true',
                dest='force', help='Force overwrite of existing output directory' +
                ' (note: existing files in output_dir will not be removed)' +
                ' [default: %default]'), options_lookup['jobs_to_start_workflow']
]

script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_dir = opts.output_dir
    otu_map_fp = opts.otu_map_fp
    otu_table_fp = opts.otu_table_fp
    fasta_fp = opts.fasta_fp
    parameter_fp = opts.parameter_fp
    force = opts.force

    qiime_config = load_qiime_config()

    # Check that specified input files do, in fact, exist
    try:
        with open(otu_map_fp) as f:
            pass
    except IOError as e:
        print 'OTU Map could not be opened! Are you sure it is located at %s?' % otu_map_fp
        exit(1)

    # Check OTU Table
    try:
        with open(otu_table_fp) as f:
            pass
    except IOError as e:
        print 'OTU Table could not be opened! Are you sure it is located at %s?' % otu_table_fp
        exit(1)

    # Check Sequences FASTA
    try:
        with open(fasta_fp) as f:
            pass
    except IOError as e:
        print 'FASTA Sequences could not be opened! Are you sure it is located at %s?' % fasta_fp
        exit(1)

    # Verify that parameters file exists, if it is specified
    if parameter_fp:
        try:
            parameter_f = open(parameter_fp)
        except IOError:
            raise IOError,\
                "Can't open parameters file (%s). Does it exist? Do you have read access?"\
                % parameter_fp
        params = parse_qiime_parameters(parameter_f)
    else:
        params = parse_qiime_parameters([])
        # empty list returns empty defaultdict for now

    try:
        os.makedirs(output_dir)
    except OSError:
        if force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose " +\
                "a different directory, or force overwrite with --force"
            exit(1)

    # these are hardcoded from options selection in pick_otus_through_otu...py
    command_handler = call_commands_serially
    status_update_callback = no_status_updates
    parallel = False

    # get parent OTU map and load it into dict otu_to_seqid
    otu_to_seqid = fields_to_dict(open(otu_map_fp, 'U'))

    print "Loading seqs..."

    # get the seqs.fna for all the sequences in the whole set
    fasta_collection = LoadSeqs(fasta_fp, moltype=DNA, aligned=False,
                                label_to_name=lambda x: x.split()[0])


    # get the pOTUs in the filtered otu table
    potu_table = load_table(otu_table_fp)

    # for each of the OTUs in this filtered parent table,
    for pOTU in potu_table.ids(axis = 'observation'):

        print "Reclustering pOTU# %s..." % str(pOTU)

        potu_dir = os.path.join(output_dir,str(pOTU))

        try:
            os.makedirs(potu_dir)
        except OSError:
            if force:
                pass
            else:
                # Since the analysis can take quite a while, I put this check
                # in to help users avoid overwriting previous output.
                print "Output directory already exists. Please choose " +\
                    "a different directory, or force overwrite with --force"
                exit(1)


        seqs_in_otu = fasta_collection.takeSeqs(otu_to_seqid[pOTU])
        output_fna_fp = os.path.join(potu_dir,'seqs.fasta')
        seqs_in_otu.writeToFile(output_fna_fp)

        # pre process seqs from pOTUs

        #try:
        run_pick_de_novo_otus(
            output_fna_fp,
            potu_dir,
            command_handler=command_handler,
            params=params,
            qiime_config=qiime_config,
            parallel=parallel,
            status_update_callback=status_update_callback,
            run_assign_tax=False)
        #except Exception, err:
        #    sys.stderr.write('ERROR: %s\n' % str(err))
            # return 1


if __name__ == "__main__":
    main()
