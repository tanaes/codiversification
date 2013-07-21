#!/n/sw/Python-2.6.2/bin/python
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
import os
import numpy
from qiime.util import load_qiime_config, parse_command_line_parameters,\
 get_options_lookup, parse_otu_table
from qiime.parse import parse_qiime_parameters, parse_taxonomy


#test_cospeciation.py
qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""
Makes iTOL data annotation files for trees."""
script_info['script_description']="""
This script reads in an OTU table and an OTU tree, and outputs a number of 
data annotation options for upload and display via iTOL."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""
The following steps are performed by the command below:

1. Read in OTU table

2. Read OTU tree to determine represented terminals

3. Generate a heatmap datafile for iTOL based on the OTU table
""","""
make_itol_annotations.py -i filtered_otu_table.txt -t input_tree.tre -o iTOL_files """))
script_info['output_description']="""
We'll see!"""
script_info['required_options']=[\
 make_option('-i','--otu_table_fp',\
            help='the input OTU table file [REQUIRED]'),\
 make_option('-t','--tree_fp',\
            help='an OTU tree [REQUIRED]'),\
 make_option('-o','--output_dir',\
            help='path to output directory [REQUIRED]'),\
]
script_info['optional_options']=[\
 make_option('-m','--heatmap',\
            help='make iTOL heatmap datafile'+\
            '[default: %default]',\
            default=True),\
 make_option('-d','--data_fp',\
            help='data filepath for barchart/piechart'+\
            '[default: %default]',\
            default=None),\
 make_option('--data_otu_header',\
            help='column header for OTU name in data file'+\
            '[default: %default]',\
            default='pOTU'),\
 make_option('-b','--barchart_headers',\
            help='datafile header(s) for barcharts, comma separated'+\
            '[default: %default]',\
            default=None),\
 make_option('-g','--piechart_headers',\
            help='datafile header(s) for pie charts, comma separated'+\
            '[default: %default]',\
            default=None),\
 make_option('--barchart_colors',\
            help='hex colors for barcharts, comma separated'+\
            '[default: %default]',\
            default=None),\
 make_option('--piechart_colors',\
            help='hex colors for pie charts, comma separated'+\
            '[default: %default]',\
            default=None),\
 make_option('--bar_denominator',\
            help='use given data header as denominator for others in barchart'+\
            '[default: %default]',\
            default=False),\
 make_option('--pie_denominator',\
            help='use given data header as denominator for previous in piechart'+\
            '[default: %default]',\
            default=False),\
 make_option('--suppress_empty_pies',\
            help='do not output piecharts with no numerator'+\
            '[default: %default]',\
            default=True),\
 make_option('--color_by_taxonomy',\
            help='make iTOL color datafile for taxonomy'+\
            '[default: %default]',\
            default=None),\
 make_option('--taxonomy_colors_fp',\
            help='color mapping for taxonomy colors'+\
            '[default: %default]',\
            default=None),\
 make_option('--taxonomy_colors_level',\
            help='level to color taxonomy by'+\
            '[default: %efault]',\
            default=4),\
 make_option('--force',action='store_true',\
            dest='force',help='Force overwrite of existing output directory'+\
            ' (note: existing files in output_dir will not be removed)'+\
            ' [default: %default]'),\
 make_option('-p','--parameter_fp',
            help='path to the parameter file, which specifies changes'+\
            ' to the default behavior. '+\
            'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters.'+\
            ' [if omitted, default values will be used]',\
            default=None),
]

script_info['version'] = __version__




def get_terminals(tree_fp):
    from cogent.core.tree import PhyloNode
    from cogent import LoadTree
    
    tree = LoadTree(tree_fp)
    terminals = tree.getTipNames()
    
    return terminals



def iTOL_piechart(results_dict,terminals,piechart_headers,pie_denominator=False,suppress_empty_pies=True,piechart_colors=None):
    #make pie chart file
    
    if not piechart_colors:
        piechart_colors = easycolor(len(piechart_headers))
    
    colors_line = 'COLORS'
    for color in piechart_colors:
        colors_line += ',' + color 
    
    headers_line = 'LABELS'
    for header in piechart_headers:
        headers_line += ',' + header
    
    outstring = headers_line + '\n' + colors_line + '\n'
    
    values_dict = {}
    totals_dict = {}
    
    for otu in set(results_dict.keys()) & set(terminals):
        values_dict[otu] = {}
        if pie_denominator:
            totals_dict[otu] = float(results_dict[otu][piechart_headers[-1]])
            for header in piechart_headers[0:-1]:
                values_dict[otu][header] = float(results_dict[otu][header])
            values_dict[otu][piechart_headers[-1]] = totals_dict[otu] - sum([float(results_dict[otu][x]) for x in piechart_headers[0:-1]])
        else:
            totals_dict[otu] = sum([float(results_dict[otu][x]) for x in piechart_headers])
            for header in piechart_headers:
                values_dict[otu][header] = float(results_dict[otu][header])
        
    for otu in values_dict.keys():
        if (values_dict[otu][piechart_headers[-1]] == totals_dict[otu]) and suppress_empty_pies:
            continue
        else:
            outstring += str(otu) + ',R' + str(int(((float(totals_dict[otu])/(max(totals_dict.values())*3.14))**(0.5))*100)) + '|1.0'
            for header in piechart_headers:
                outstring += ',' + str(int(values_dict[otu][header]))
            outstring += '\n'
    
    outfile = open('iTOL_piechart.txt','w')
    outfile.write(outstring)
    outfile.close()
    
    return outstring

def iTOL_barchart():
    """
    #make bar chart file
            
    header = 'LABELS,sig,total\n' + 'COLORS,#ff0000,#dddddd\n'
    notax = header
    tax = header
    
    for otu in set(results_dict.keys()) & set(potu_tree.getTipNames()):
        if int(results_dict[otu]['fdr_sig_nodes']) >= 0:
            notax += otu + ',' +  results_dict[otu]['fdr_sig_nodes'] + ',' + str(int(results_dict[otu]['total_nodes']) - int(results_dict[otu]['fdr_sig_nodes'])) + '\n'
    """
    return True

def easycolor(numcolors):
    basecolors = [[46,87,140],\
                [93,50,73],\
                [232,160,60],\
                [188,44,48],\
                [111,60,121],\
                [125,127,127]]
    
    mycolors=[]
    iteration = 0
    while iteration < numcolors:
        for shade in range(6):
            
            for hue in range(6):
                
                if iteration >= numcolors:
                    break
                
                mycolors.append([int((x + (255-x)*((float(shade))/6))) for x in basecolors[hue]])
                iteration += 1
    outcolors = []
    for mycolor in mycolors:
        outcolors.append('#' + str(hex(mycolor[0])[2:]) + str(hex(mycolor[1])[2:]) + str(hex(mycolor[2])[2:]))
        
    return outcolors


def color_file_to_dict(taxonomy_colors_fp):
    color_file = open(taxonomy_colors_fp,'r')
    color_dict = {}
    for line in color_file:
        
        color_dict[line.strip().split('\t')[0]] = line.strip().split('\t')[1]
    color_file.close()
    return color_dict

def iTOL_taxcolors(taxon_names,lineages,terminals,level=4,taxonomy_colors_fp=None):
    outstring = 'NODE_ID\tTYPE\tCOLOR\tLABEL\n'
    
    color_dict = {}
    
    if taxonomy_colors_fp:
        color_dict = color_file_to_dict(taxonomy_colors_fp)
    
    taxon_colors = easycolor(len(taxon_names))
    taxon_colors.reverse()
    
    for i in range(len(taxon_names)):
        
        otu = taxon_names[i]
        
        if otu in terminals:
            if len(lineages[i]) < (level + 1):
                taxonomy = lineages[i][-1]
            else:
                taxonomy = lineages[i][level]
            if taxonomy not in color_dict.keys():
                color_dict[taxonomy] = taxon_colors.pop()
            outstring += str(otu) + '\trange\t' + color_dict[taxonomy] + '\t' + taxonomy + '\n'
    
    outfile = open('iTOL_taxcolors.txt','w')
    outfile.write(outstring)
    outfile.close()
    
    return outstring


def iTOL_heatmap(sample_names, taxon_names, data, lineages, terminals):
    
    header = 'LABELS'
    
    for sample in sample_names:
        header += ',' + sample
    
    persample_string = header + '\n'
    perotu_string = header + '\n'
    
    per_sample_sums = data.sum(0)
    
    per_sample_norm = []
    
    data = data.astype(float)
    
    #normalize overall reads to 10000 per category
    for j in range(len(sample_names)):
        per_sample_norm.append(float(10000)/float(per_sample_sums[j]))
    
    for j in range(len(per_sample_norm)):
        data[:,j] = data[:,j]*per_sample_norm[j]
    
    data = numpy.ceil(data)
    
    perotu_data = data.copy()
    
    #normalize to 255 reads per OTU
    for i in range(len(taxon_names)):
        perotu_data[i,:] = data[i,:]/data[i,:].sum()*255
    
    perotu_data = numpy.ceil(perotu_data)
    
    for i in range(len(taxon_names)):
        otu = taxon_names[i]
        if otu in terminals:
            sampleline = otu
            otuline = otu
            for j in range(len(sample_names)):
                sampleline += ',' + str(int(data[i,j]))
                otuline += ',' + str(int(perotu_data[i,j]))
            persample_string += sampleline + '\n'
            perotu_string += otuline + '\n'
            #print sampleline
            #print otuline
    heatmap_persample_file = open('iTOL_heatmap_per_sample.txt','w')
    heatmap_perotu_file = open('iTOL_heatmap_per_otu.txt','w')
    heatmap_persample_file.write(persample_string)
    heatmap_perotu_file.write(perotu_string)
    heatmap_persample_file.close()
    heatmap_perotu_file.close()
    
    return True

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    
    
    otu_table_fp = opts.otu_table_fp
    tree_fp = opts.tree_fp
    output_dir = opts.output_dir
    data_fp = opts.data_fp
    data_otu_header = opts.data_otu_header
    heatmap = opts.heatmap
    
    pie_denominator = opts.pie_denominator
    bar_denominator = opts.bar_denominator
    suppress_empty_pies = opts.suppress_empty_pies
    
    color_by_taxonomy = opts.color_by_taxonomy
    taxonomy_colors_fp = opts.taxonomy_colors_fp
    taxonomy_colors_level = int(opts.taxonomy_colors_level)
    force = opts.force
    
    barchart_headers = opts.barchart_headers
    piechart_headers = opts.piechart_headers
    barchart_colors = opts.barchart_colors
    piechart_colors = opts.piechart_colors
    
    if opts.barchart_headers:
        print opts.barchart_headers
        barchart_headers = opts.barchart_headers.split(',')
    if opts.piechart_headers:
        print opts.piechart_headers
        piechart_headers = opts.piechart_headers.split(',')
    if opts.barchart_colors:
        print opts.barchart_colors
        barchart_colors = opts.barchart_colors.split(',')
    if opts.piechart_colors:
        print opts.piechart_colors
        piechart_colors = opts.piechart_colors.split(',')
        print piechart_colors
    
    
    
    try:
        os.makedirs(output_dir)
    except OSError:
        if force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)
        
    otu_table_fp=os.path.abspath(otu_table_fp)
    tree_fp=os.path.abspath(tree_fp)
    output_dir=os.path.abspath(output_dir)
    os.chdir(output_dir)
    
    
    if data_fp:
        results_file = open(data_fp,'r')
        results_dict = {}
        headers = results_file.next().split('\t')
        for line in results_file:
            #print line
            values = line.split('\t')
            line_dict = {}
            for i in range(len(values)):
                line_dict[headers[i]] = values[i]
            results_dict[line_dict[data_otu_header]] = line_dict
        
        results_file.close()
    
    
    
    if os.path.isfile(otu_table_fp):
        #run test on all otu tables in directory
        otu_file = open(otu_table_fp, 'Ur')
        
        #parse OTU table 
        try:
            sample_names, taxon_names, data, lineages = parse_otu_table(otu_file)
        except:
            print otu_file + " not an OTU table?"
        
        #get tree taxa
        try:
            terminals = get_terminals(tree_fp)
        except:
            print tree_fp + " not a good tree?"
        
        
        if heatmap:
            iTOL_heatmap(sample_names, taxon_names, data, lineages, terminals)
        
        if piechart_headers:
            iTOL_piechart(results_dict,terminals,piechart_headers,pie_denominator,suppress_empty_pies,piechart_colors)
        
        if color_by_taxonomy:
            iTOL_taxcolors(taxon_names,lineages,terminals,taxonomy_colors_level,taxonomy_colors_fp)
        
        #print_results
        
    else:
        print "Can't find OTU table!"
    




if __name__ == "__main__":
    main()



