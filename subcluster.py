#!/n/sw/python-2.7.1/bin/python
# File created on 09 Aug 2012

__author__ = "Jon Sanders"
__copyright__ = "Copyright 2012"
__credits__ = ["Rob Knight", "Justin Kuczynski","Jesse Stombaugh","Jon Sanders"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Jon Sanders"
__email__ = "jonsan@gmail.com"
__status__ = "Experimental"


#importing modules 
from sys import exit, stderr, stdout
from qiime.parse import fields_to_dict, parse_qiime_parameters
from os import makedirs
from os.path import split, splitext, join, dirname, abspath
from qiime.util import (make_option,
						qiime_system_call,
						create_dir,
						get_qiime_scripts_dir,
						load_qiime_config, 
						parse_command_line_parameters,
						parse_otu_table,
						get_options_lookup)
from qiime.workflow import (run_qiime_data_preparation, print_commands,
	call_commands_serially, print_to_stdout, no_status_updates,
	validate_and_set_jobs_to_start, WorkflowLogger, generate_log_fp, get_params_str)

from cogent import LoadSeqs, DNA

qiime_config = load_qiime_config()

def sub_otu_prep(input_fp, 
							   output_dir, 
							   command_handler,
							   params, 
							   qiime_config,
							   parallel=False,
							   logger=None,
							   status_update_callback=print_to_stdout):
	""" Run the data preparation steps of Qiime for OTU subclustering
	
		The steps performed by this function are:
		  1) Pick OTUs;
		  2) Pick a representative set;
		  3) Align the representative set;
		  4) Filter the alignment prior to tree building - remove positions
			 which are all gaps, and specified as 0 in the lanemask
		  5) Build a phylogenetic tree;
		  6) Build an OTU table.
	
	"""
	
	# Prepare some variables for the later steps
	input_dir, input_filename = split(input_fp)
	input_basename, input_ext = splitext(input_filename)
	create_dir(output_dir)
	commands = []
	python_exe_fp = qiime_config['python_exe_fp']
	script_dir = get_qiime_scripts_dir()
	if logger == None:
		logger = WorkflowLogger(generate_log_fp(output_dir),
								params=params,
								qiime_config=qiime_config)
		close_logger_on_success = True
	else:
		close_logger_on_success = False
	
	# Prep the OTU picking command
	try:
		otu_picking_method = params['pick_otus']['otu_picking_method']
	except KeyError:
		otu_picking_method = 'uclust'
		
	#pick_otu_dir = '%s/%s_picked_otus' % (output_dir, otu_picking_method)
	
	#Don't want lots of folders, dump all cOTU files into output_dir
	pick_otu_dir = output_dir
	
	#Name cOTU files based on pOTU filename
	otu_fp = '%s/%s_otus.txt' % (pick_otu_dir,input_basename)


	# Get OTU picking parameters
	
	try:
		params_str = get_params_str(params['pick_otus'])
	except KeyError:
		params_str = ''
	# Build the OTU picking command
	pick_otus_cmd = '%s %s/pick_otus.py -i %s -o %s %s' %\
	 (python_exe_fp, script_dir, input_fp, pick_otu_dir, params_str)

	commands.append([('Pick OTUs', pick_otus_cmd)])
	
	# Prep the representative set picking command
	#rep_set_dir = '%s/rep_set/' % output_dir
	#again, dump all rep sets into the output directory
	rep_set_dir = output_dir
	
	rep_set_fp = '%s/%s_rep_set.fasta' % (rep_set_dir,input_basename)
	rep_set_log_fp = '%s/%s_rep_set.log' % (rep_set_dir,input_basename)
	
	try:
		params_str = get_params_str(params['pick_rep_set'])
	except KeyError:
		params_str = ''
	# Build the representative set picking command
	pick_rep_set_cmd = '%s %s/pick_rep_set.py -i %s -f %s -l %s -o %s %s' %\
	 (python_exe_fp, script_dir, otu_fp, input_fp, rep_set_log_fp,\
	  rep_set_fp, params_str)
	commands.append([('Pick representative set', pick_rep_set_cmd)])
	
	
	
	#Don't want to re-assign taxonomy at this point, but leaving in just in case
	"""
	# Prep the taxonomy assignment command
	try:
		assignment_method = params['assign_taxonomy']['assignment_method']
	except KeyError:
		assignment_method = 'rdp'
	assign_taxonomy_dir = '%s/%s_assigned_taxonomy' %\
	 (output_dir,assignment_method)
	taxonomy_fp = '%s/%s_rep_set_tax_assignments.txt' % \
	 (assign_taxonomy_dir,input_basename)
	if parallel and (assignment_method == 'rdp' or assignment_method == 'blast'):
		# Grab the parallel-specific parameters
		try:
			params_str = get_params_str(params['parallel'])
		except KeyError:
			params_str = ''
		
		# Grab the OTU picker parameters
		try:
			# Want to find a cleaner strategy for this: the parallel script
			# is method-specific, so doesn't take a --assignment_method
			# option. This works for now though.
			d = params['assign_taxonomy'].copy()
			del d['assignment_method']
			params_str += ' %s' % get_params_str(d)
		except KeyError:
			pass
			
		# Build the parallel taxonomy assignment command
		assign_taxonomy_cmd = \
		 '%s %s/parallel_assign_taxonomy_%s.py -i %s -o %s -T %s' %\
		 (python_exe_fp, script_dir, assignment_method, rep_set_fp,\
		  assign_taxonomy_dir, params_str)
	else:
		try:
			params_str = get_params_str(params['assign_taxonomy'])
		except KeyError:
			params_str = ''
		# Build the taxonomy assignment command
		assign_taxonomy_cmd = '%s %s/assign_taxonomy.py -o %s -i %s %s' %\
		 (python_exe_fp, script_dir, assign_taxonomy_dir,\
		  rep_set_fp, params_str)
	
	commands.append([('Assign taxonomy',assign_taxonomy_cmd)])
	"""
	
	# Prep the OTU table building command
	otu_table_fp = '%s/%s_otu_table.txt' % (output_dir, input_basename)
	try:
		params_str = get_params_str(params['make_otu_table'])
	except KeyError:
		params_str = ''
	# Build the OTU table building command
	make_otu_table_cmd = '%s %s/make_otu_table.py -i %s -o %s %s' %\
	 (python_exe_fp, script_dir, otu_fp, otu_table_fp, params_str)
	
	commands.append([('Make OTU table', make_otu_table_cmd)])
	
	
	#Leaving in all alignment options; set in params file
	
	# Prep the pynast alignment command
	try:
		alignment_method = params['align_seqs']['alignment_method']
	except KeyError:
		alignment_method = 'pynast'
	pynast_dir = output_dir
	aln_fp = '%s/%s_rep_set_aligned.fasta' % (pynast_dir,input_basename)
	
	#Currently hardcoding serial alignment, so this is not relevant.
	if parallel and alignment_method == 'pynast':
		# Grab the parallel-specific parameters
		try:
			params_str = get_params_str(params['parallel'])
		except KeyError:
			params_str = ''
		
		# Grab the OTU picker parameters
		try:
			# Want to find a cleaner strategy for this: the parallel script
			# is method-specific, so doesn't take a --alignment_method
			# option. This works for now though.
			d = params['align_seqs'].copy()
			del d['alignment_method']
			params_str += ' %s' % get_params_str(d)
		except KeyError:
			pass
			
		# Build the parallel pynast alignment command
		align_seqs_cmd = '%s %s/parallel_align_seqs_pynast.py -i %s -o %s -T %s' %\
		 (python_exe_fp, script_dir, rep_set_fp, pynast_dir, params_str)
	else:
		try:
			params_str = get_params_str(params['align_seqs'])
		except KeyError:
			params_str = ''
		# Build the pynast alignment command
		align_seqs_cmd = '%s %s/align_seqs.py -i %s -o %s %s' %\
		 (python_exe_fp, script_dir, rep_set_fp, pynast_dir, params_str)
	commands.append([('Align sequences', align_seqs_cmd)])
	
	if alignment_method == 'pynast':
		# Prep the alignment filtering command (only applicable when aligned
		# with pynast)
		filtered_aln_fp = '%s/%s_rep_set_aligned_pfiltered.fasta' %\
		 (pynast_dir,input_basename)
		try:
			params_str = get_params_str(params['filter_alignment'])
		except KeyError:
			params_str = ''
		# Build the alignment filtering command
		filter_alignment_cmd = '%s %s/filter_alignment.py -o %s -i %s %s' %\
		 (python_exe_fp, script_dir, pynast_dir, aln_fp, params_str)
		commands.append([('Filter alignment', filter_alignment_cmd)])
	else: 
		filtered_aln_fp = aln_fp
	
	# Prep the tree building command
	tree_fp = '%s/%s_rep_set.tre' % (output_dir, input_basename)
	try:
		params_str = get_params_str(params['make_phylogeny'])
	except KeyError:
		params_str = ''
	# Build the tree building command
	make_phylogeny_cmd = '%s %s/make_phylogeny.py -i %s -o %s %s' %\
	 (python_exe_fp, script_dir, filtered_aln_fp, tree_fp,\
	 params_str)
	commands.append([('Build phylogenetic tree', make_phylogeny_cmd)])
	
	
	# Call the command handler on the list of commands. BINGO!
	command_handler(commands,
					status_update_callback,
					logger=logger,
					close_logger_on_success=close_logger_on_success)
	
	return abspath(otu_table_fp)


def otu_subcluster(output_dir,otu_map_fp,otu_table_fp,fasta_fp,parameter_fp,force):

	
	#Check that specified input files do, in fact, exist
	try:
		with open(otu_map_fp) as f:
			pass
	except IOError as e:
   		print 'OTU Map could not be opened! Are you sure it is located at ' + otu_map_fp + '  ?'
		exit(1)
	
	#Check OTU Table
	try:
   		with open(otu_table_fp) as f:
			pass
	except IOError as e:
   		print 'OTU Table could not be opened! Are you sure it is located at ' + otu_table_fp + '  ?'
		exit(1)
	
	#Check Sequences FASTA
	try:
   		with open(fasta_fp) as f:
			pass
	except IOError as e:
   		print 'FASTA Sequences could not be opened! Are you sure it is located at ' + fasta_fp + '  ?'
		exit(1)
	
	#Verify that parameters file exists, if it is specified
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
		makedirs(output_dir)
	except OSError:
		if force:
			pass
		else:
			# Since the analysis can take quite a while, I put this check
			# in to help users avoid overwriting previous output.
			print "Output directory already exists. Please choose "+\
			 "a different directory, or force overwrite with --force"
			exit(1)
	
	
	#these are hardcoded from options selection in pick_otus_through_otu...py
	command_handler = call_commands_serially
	status_update_callback = no_status_updates
	parallel = False
	
	#get parent OTU map and load it into dict otu_to_seqid
	otu_to_seqid = fields_to_dict(open(otu_map_fp, 'U'))
	
	print "Loading seqs..."
	
	#get the seqs.fna for all the sequences in the whole set
	fasta_collection = LoadSeqs(fasta_fp, moltype=DNA, aligned=False,
	 label_to_name=lambda x: x.split()[0])
		
	#testing code:
	#print "OTU dict keys:"
	
	#print otu_to_seqid.keys()[0:3]
	
	#pause = raw_input('Pause!')

	print "Parsing OTU table..."	
		
	#get the OTUs in the filtered otu table
	sample_names, passed_otus, data, lineages = parse_otu_table(open(otu_table_fp, 'Ur'))
	
	#for each of the OTUs in this filtered parent table,
	for OTU in passed_otus:
		
		print "Reclustering OTU# " + str(OTU) + "..."
		
		seqs_in_otu = fasta_collection.takeSeqs(otu_to_seqid[OTU])
		output_fna_fp = output_dir + '/' + str(OTU) + '_seqs.fasta'
		seqs_in_otu.writeToFile(output_fna_fp)
		
		
		#pre process seqs from OTUs
		
		try:
			sub_otu_prep(
			 output_fna_fp, 
			 output_dir,
			 command_handler=command_handler,
			 params=params,
			 qiime_config=qiime_config,
			 parallel=parallel,\
			 status_update_callback=status_update_callback)
		except Exception, err:
			stderr.write('ERROR: %s\n' % str(err))
        	#return 1
	