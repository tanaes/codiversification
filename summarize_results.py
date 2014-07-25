#!/n/sw/python-2.7.1/bin/python
# File created on 3 October 2011.
from __future__ import division

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2011, Jon Sanders"
__credits__ = ["Jon Sanders"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Experimental"

from qiime.util import make_option
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 get_options_lookup
from qiime.parse import parse_qiime_parameters, parse_taxonomy
import os.path
import os




#test_cospeciation.py
qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""
Corrects cospeciation test results for multiple tests, and generates annotated trees."""
script_info['script_description']="""
This script reads in the output folder from test_cospeciation.py, corrects for 
multiple tests using Bonferroni & False Discovery Rate (fdr) correction, and 
returns a corrected results file for each previously-significant parent OTU. In 
addition, it returns a summary file (corrected_results_summary.txt) with 
taxonomy info and # of significant nodes for each parent OTU, plus files with a 
list of all significant nodes for each type of multiple test correction."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""
The following command reads in a test_cospeciation output folder, applies 
multiple test correction to each node's significance value, optionally looks up
taxonomy information for each parent OTU, and returns a summary file indicating 
the number of significant child nodes within each parent OTU for each type of 
multiple test correction at the given threshold, and a list of all significant
nodes for each multiple test correction.

""","""
python summarize_results.py -i otu_dir -r results_dir -o out_dir 
"""))
script_info['output_description']="""
This script outputs a new folder containing a new summary overview text file 
indicating the number of significant nodes in each pOTU at the given threshold, 
the number of nodes significant under bonferroni and FDR correction, the total 
number of nodes in that pOTU, and the assigned taxonomy. In addition, this file 
includes a column indicating 'max_fdr_h_ratio,' which reports the maximum ratio 
of h_tips to h_span for a node significant under FDR correction in that pOTU. 

In addition, a new nodes file is returned for each pOTU with bonferroni and FDR-
corrected p-vals for each node, as well as files containing all significant 
nodes under each threshold of correction (none, FDR, and bonferroni).
"""
script_info['required_options']=[\
 make_option('-i','--cOTUs_dir',\
            help='the cOTUs directory used for the cospecitaion analysis [REQUIRED]'),\
 make_option('-r','--results_dir',\
            help='the directory(ies) of output from test_cospeciation.py. [REQUIRED]  '+\
            '(If multiple directories are given, results will be combined)'),\
 make_option('-o','--output_dir',\
            help='path to output directory [REQUIRED]'),\
]

script_info['optional_options']=[\
 make_option('-t','--taxonomy_fp', \
            help='parent OTU taxonomy filepath '),\
 make_option('-s','--significance_level',type='float',\
            help='Desired level of significance for permutation test '+\
            '[default: %default]', default=0.05),\
 make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
 options_lookup['jobs_to_start_workflow'] \
]

script_info['version'] = __version__
def read_results_keys(results_fp):
    #This reads in the per-pOTU results file and sticks it in dicts in a list
    results_file = open(results_fp,'Ur')
    
    results_list = []
    results_keys = results_file.next().strip().split('\t')
    
    return results_keys

def read_results_file(results_fp):
    #This reads in the per-pOTU results file and sticks it in dicts in a list
    results_file = open(results_fp,'Ur')
    
    results_list = []
    results_keys = results_file.next().strip().split('\t')
    
    pOTU = os.path.basename(results_fp).split('_')[0]
    
    for line in results_file:
        entry_dict = dict(zip(results_keys,line.strip().split('\t')))
        entry_dict['pOTU'] = pOTU
        results_list.append(entry_dict)
    
    results_file.close()
    
    return results_list

def print_corrected_results_files(results_list,results_keys,p_dict):
    
    #initialize pOTU counter to 0
    last_pOTU = results_list[0]['pOTU']
    outline = ''
    
    headers = results_keys
    
    headerline = ''
    for key in headers:
        headerline += key + '\t'
    
    headerline += 'fdr_p\tbonf_p\t\n'
    
    #iterate through each node
    for node in range(len(results_list)):
        thisline = ''
        
        #for each value, 
        for key in headers:
            thisline += results_list[node][key] + '\t'
        
        #Append the fdr_p and bonf_p to the results line
        thisline += str(p_dict[node][2]) + '\t' + str(p_dict[node][3]) + '\t\n'
        this_pOTU = results_list[node]['pOTU']
        
        #if it's a new pOTU, write the previous one.
        if this_pOTU != last_pOTU:
            outfile = open(last_pOTU + '_corrected_results.txt','w')
            outfile.write(headerline)
            outfile.write(outline)
            outfile.close()
            
            outline = thisline
            last_pOTU = this_pOTU
        else:
            outline += thisline
            last_pOTU = this_pOTU
    
    #write the final pOTU file.
    outfile = open(last_pOTU + '_corrected_results.txt','w')
    outfile.write(headerline)
    outfile.write(outline)
    outfile.close()
    
    return True


def print_sig_lists(results_list,results_keys,p_dict,significance_level=0.05,otu_to_taxonomy={}):
    
    alpha = significance_level
    
    #initialize empty outstrings
    p_list = ''
    fdr_p_list = ''
    bonf_p_list = ''
    pOTU_summary = ''
    
    p_list_headers = ''
    for key in results_keys:
        p_list_headers += key + '\t'
    
    pOTU_dict = {}
    
    #initialize headers for pOTU summary table and lists of nodes
    if otu_to_taxonomy:
        
        pOTU_summary = 'pOTU\tsig_nodes\tfdr_sig_nodes\tbonf_sig_nodes\ttotal_nodes\tmax_h_ratio\ttaxonomy\t\n'
        
        p_list = 'pOTU\ttaxonomy\t' + p_list_headers + '\n'
        bonf_p_list = 'pOTU\tbonf_p\ttaxonomy\t' + p_list_headers + '\n'
        fdr_p_list = 'pOTU\tfdr_p\ttaxonomy\t' + p_list_headers + '\n'
    else:
        pOTU_summary = 'pOTU\tsig_nodes\tfdr_sig_nodes\tbonf_sig_nodes\ttotal_nodes\tmax_fdr_h_ratio\t\n'
        
        p_list = 'pOTU\t' + p_list_headers + '\n'
        bonf_p_list = 'pOTU\tbonf_p\t' + p_list_headers + '\n'
        fdr_p_list = 'pOTU\tfdr_p\t' + p_list_headers + '\n'
    
    #now iterate over each node in the results_list
    for node in range(len(results_list)):
        
        node_data = ''
        
        pOTU = results_list[node]['pOTU']
        
        #convert data for each node to a tab-delimited string
        #add taxonomy info if appropriate
        if otu_to_taxonomy:
            node_data += otu_to_taxonomy[pOTU] + '\t'
            
        for key in results_keys:
            node_data += results_list[node][key] + '\t'
        
        h_ratio = float(results_list[node]['h_tips'])/float(results_list[node]['h_span'])
        
        #If it's a new pOTU, initialize the dict with a five-element list. 
        #elements are: [0] number of nodes, [1] uncorrected signficant nodes,
        #[2] fdr_significant nodes, [3] bonf_significant nodes, and [4] the max
        #h_ratio for that pOTU
        if pOTU not in pOTU_dict:
            pOTU_dict[pOTU] = [0,0,0,0,0]
        
        pOTU_dict[pOTU][0] += 1
        
        if p_dict[node][1] < alpha:
            pOTU_dict[pOTU][1] += 1
            p_list += str(pOTU) + '\t' + node_data + '\n'
            
        if p_dict[node][2] < alpha:
            pOTU_dict[pOTU][2] += 1
            fdr_p_list += str(pOTU) + '\t' + str(p_dict[node][2]) + '\t' + node_data + '\n'
            
        if p_dict[node][3] < alpha:
            pOTU_dict[pOTU][3] += 1
            bonf_p_list += str(pOTU) + '\t' + str(p_dict[node][3]) + '\t' + node_data + '\n'
            
        if h_ratio > pOTU_dict[pOTU][4] and p_dict[node][2] < alpha:
            pOTU_dict[pOTU][4] = h_ratio
    
    #make an integer sorted list of pOTU keys
    pOTUs = sorted([int(key) for key in pOTU_dict.keys()])
    
    for pOTU in pOTUs:
        pOTU_summary += str(pOTU) + '\t'
        pOTU_summary += str(pOTU_dict[str(pOTU)][1]) + '\t' 
        pOTU_summary += str(pOTU_dict[str(pOTU)][2]) + '\t' 
        pOTU_summary += str(pOTU_dict[str(pOTU)][3]) + '\t' 
        pOTU_summary += str(pOTU_dict[str(pOTU)][0]) + '\t' 
        pOTU_summary += str(pOTU_dict[str(pOTU)][4])
        
        if otu_to_taxonomy:
            pOTU_summary += '\t' + otu_to_taxonomy[str(pOTU)] + '\n'
        else:
            pOTU_summary += '\n'
    
    p_list_file = open('uncorrected_sig_nodes.txt','w')
    p_list_file.write(p_list)
    p_list_file.close()
    
    fdr_p_list_file = open('fdr_sig_nodes.txt','w')
    fdr_p_list_file.write(fdr_p_list)
    fdr_p_list_file.close()
    
    bonf_p_list_file = open('bonferroni_sig_nodes.txt','w')
    bonf_p_list_file.write(bonf_p_list)
    bonf_p_list_file.close()
    
    summary_file = open('corrected_results_summary.txt','w')
    summary_file.write(pOTU_summary)
    summary_file.close()
    
    return

def de_NA(p_dict):
    for key in p_dict.keys():
        for i in range(len(p_dict[key])):
            if p_dict[key][1] == 0 and p_dict[key][i] == 'NA':
                p_dict[key][i] = 0
    return p_dict

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    cOTUs_dir = opts.cOTUs_dir
    results_dir = opts.results_dir
    output_dir = opts.output_dir
    significance_level = float(opts.significance_level)
    taxonomy_fp = opts.taxonomy_fp
    force = opts.force
    
    from qiime.otu_category_significance import add_fdr_correction_to_results, \
                                                add_bonferroni_to_results, \
                                                fdr_correction
    
    #test input and output dirs
    
    if opts.taxonomy_fp:
        try:
            taxonomy_file = open(opts.taxonomy_fp,'Ur')
            otu_to_taxonomy = parse_taxonomy(open(taxonomy_fp,'Ur'))
            #
        except IOError:
            raise IOError,\
             "Can't open taxonomy file (%s). Does it exist? Do you have read access?"\
             % opts.taxonomy_fp
    else:
        otu_to_taxonomy=None
    
    if not os.path.isdir(cOTUs_dir):
        print "cOTUs_directory not a directory. Please try again."
        exit(1)
    
    cOTUs_dir = os.path.abspath(cOTUs_dir)
    
    if not os.path.isdir(results_dir):
        print "results_directory not a directory. Please try again."
        exit(1)
    
    results_dir = os.path.abspath(results_dir)
    
    try:
        os.makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)
    
    output_dir = os.path.abspath(output_dir)
    
    #get results dict
    #   results dict is 2D, with key a sequential per-node UID and fields for:
    #   pOTU, uncorrected pval,  taxonomy, plus other results values
    
    results_list = []
    results_keys = []
    os.chdir(results_dir)
    for file in os.listdir('.'):
        if file.endswith("results.txt"):
            
            results_list += read_results_file(file)
            
            if results_keys == []:
                results_keys = read_results_keys(file)
        
        #do FDR correction
        #   now results dict has FDR and bonfo vals
        
    p_dict = {}
    
    #we're making a dict here for all the nodes that have been tested, with the
    #key corresponding to the position in the results_list array of dicts. 
    for node in range(len(results_list)):
        pval = float(results_list[node]['p_vals'])
        p_dict[node] = [pval,pval]
    
    add_fdr_correction_to_results(p_dict)
    add_bonferroni_to_results(p_dict)
    
    #a previous iteration of the permutation test allowed 0.0 p_vals, which gave
    #'NA' results for highly significant nodes. Retained for legacy purposes. 
    p_dict = de_NA(p_dict)
    
    os.chdir(output_dir)
    
    print_corrected_results_files(results_list,results_keys,p_dict)
    
    print_sig_lists(results_list,results_keys,p_dict,significance_level,otu_to_taxonomy)
    

if __name__ == "__main__":
    main()
