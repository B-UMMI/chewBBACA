#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module 

Code documentation
------------------
"""


import os

import pandas as pd

try:
	from utils import (constants as ct,
					   gene_prediction as gp,
					   file_operations as fo,
					   fasta_operations as fao,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 gene_prediction as gp,
								 file_operations as fo,
								 fasta_operations as fao,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def write_coordinates_file(coordinates_file, output_file):
	"""Write genome CDS coordinates to a TSV file.

	Parameters
	----------
	coordinates_file : str
		Path to the pickle file that contains data about
		the CDSs coordinates.
	output_file : str
		Path to the output TSV file.
	"""
	data = fo.pickle_loader(coordinates_file)
	lines = [coords for h, coords in data[0].items()]
	lines = im.flatten_list(lines)
	lines = ['\t'.join(line) for line in lines]
	fo.write_lines(lines, output_file)


def predict_genes(fasta_files, ptf_path, translation_table,
				  prodigal_mode, cpu_cores, output_directory,
				  output_formats):
	"""Execute Prodigal to predict coding sequences from Fasta files.

	Parameters
	----------
	fasta_files : list
		List of paths to FASTA files with genomic
		sequences.
	ptf_path : str
		Path to the Prodigal training file. Should
		be NoneType if a training file is not provided.
	translation_table : int
		Genetic code used to predict and translate
		coding sequences.
	prodigal_mode : str
		Prodigal execution mode.
	cpu_cores : int
		Number of processes that will run Prodigal in
		parallel.
	output_directory : str
		Path to the directory where output files
		with Prodigal's results will be stored in.

	Returns
	-------
	failed_info : list
		List that contains a list with the stderr for the
		cases that Prodigal failed to predict genes for
		and the path to the file with information about
		failed cases. Returns NoneType if gene prediction
		succeeded for all inputs.
	"""
	if ptf_path is not None:
		# Read training file to create GeneFinder object
		training_data = gp.read_training_file(ptf_path)
		# Create GeneFinder object based on training data
		gene_finder = gp.create_gene_finder(training_data, True, True, False)
	elif ptf_path is None and prodigal_mode == 'meta':
		# Create GeneFinder object to run in meta mode
		gene_finder = gp.create_gene_finder(None, True, True, True)
	else:
		gene_finder = None

	common_args = [output_directory, gene_finder, translation_table, output_formats]

	# Divide into equal number of sublists for maximum progress resolution
	pyrodigal_inputs = im.divide_list_into_n_chunks(list(fasta_files.items()),
													len(fasta_files))

	# Add common arguments to all sublists
	pyrodigal_inputs = im.multiprocessing_inputs(pyrodigal_inputs,
												 common_args,
												 gp.predict_genome_genes)

	# Run Pyrodigal to predict genes
	# Need to use ThreadPool. Pyrodigal might hang when using Pool
	pyrodigal_results = mo.map_async_parallelizer(pyrodigal_inputs,
												  mo.function_helper,
												  cpu_cores,
												  show_progress=True,
												  pool_type='threadpool')

	# Get number of inputs for which gene prediction failed
	# Inputs with 0 CDSs and inputs with error messages
	failed = {line[0][0]: line[1]
			  for line in pyrodigal_results
			  if line[1] == 0
			  or isinstance(line[1], str) is True}

	# Get number of CDSs predicted per valid input
	cds_counts = {line[0][1]: line[1]
				  for line in pyrodigal_results
				  if isinstance(line[1], int) is True}

	# Get total number of CDSs predicted
	total_cds = sum([line[1]
					 for line in pyrodigal_results
					 if isinstance(line[1], int) is True])

	# Get paths to FASTA files with the extracted CDSs
	cds_fastas = [line[-1][0] for line in pyrodigal_results if line[-1][0] is not None]
	# Get paths to files with the coordinates of the CDSs extracted for each input
	cds_hashes = {line[0][1]: line[-1][-1] for line in pyrodigal_results if line[-1][-1] is not None}

	# Merge dictionaries with info about CDSs close to contig tips
	close_to_tip = [line[2] for line in pyrodigal_results if len(line[2]) > 0]
	close_to_tip = im.merge_dictionaries(close_to_tip)

	return [failed, total_cds, cds_fastas, cds_hashes, cds_counts, close_to_tip]


def main(input_files, output_directory, training_file, translation_table, prodigal_mode,
		 training_reference, just_training, output_formats, minimum_confidence, cpu_cores):
	# Read file with paths to input files
	if isinstance(input_files, str):
		input_files = fo.read_lines(input_files, strip=True)
	# Passed list of file paths
	elif isinstance(input_files, list):
		pass
	# Map full paths to unique identifier (prefix before first '.')
	full_to_basename = im.mapping_function(input_files, fo.file_basename, [False])
	full_to_unique = {k: fo.split_joiner(v, [0], '.')
						for k, v in full_to_basename.items()}

	# Create directory to store files with Pyrodigal results
	pyrodigal_path = fo.join_paths(output_directory, ['CDS_files'])
	fo.create_directory(pyrodigal_path)

	# Gene prediction step
	print(f'Predicting CDSs for {len(input_files)} inputs...')
	pyrodigal_results = predict_genes(full_to_unique, training_file, translation_table,
									  prodigal_mode, cpu_cores, pyrodigal_path, output_formats)

	# Dictionary with info about inputs for which gene prediction failed
	# Total number of CDSs identified in the inputs
	# Paths to FASTA files with the extracted CDSs
	# Paths to files with the coordinates of the CDSs extracted for each input
	# Total number of CDSs identified per input
	# Dictionary with info about the CDSs closer to contig tips per input
	failed, total_extracted, cds_fastas, cds_coordinates, cds_counts, close_to_tip = pyrodigal_results
	print()

	if len(failed) > 0:
		print(f'Failed to predict CDSs for {len(failed)} inputs.')
		print('Make sure that Pyrodigal runs in meta mode (--pm meta) '
				'if any input file has less than 100kbp.')
	if len(cds_fastas) == 0:
		sys.exit(f'{ct.CANNOT_PREDICT}')

	# Create TSV file with CDS coordinates
	print(f'Creating file with the coordinates of CDSs identified in inputs ({ct.CDS_COORDINATES_BASENAME})...')
	files = []
	for gid, file in cds_coordinates.items():
		tsv_file = fo.join_paths(os.path.dirname(file), [f'{gid}_coordinates.tsv'])
		write_coordinates_file(file, tsv_file)
		files.append(tsv_file)

	# Concatenate all TSV files with CDS coordinates
	merged_coordinates = fo.join_paths(output_directory, [ct.CDS_COORDINATES_BASENAME])
	fo.concatenate_files(files, merged_coordinates, header=ct.CDS_TABLE_HEADER+'\n')
	fo.remove_files(files)
	print(f'Extracted a total of {total_extracted} CDSs from {len(input_files)-len(failed)} inputs.')

	if len(failed) > 0:
		# Write Prodigal stderr for inputs that failed gene prediction
		failed_lines = [f'{k}\t{v}' for k, v in failed.items()]
		failed_outfile = fo.join_paths(os.path.dirname(output_directory),
									   ['gene_prediction_failures.tsv'])
		fo.write_lines(failed_lines, failed_outfile)

	return failed, total_extracted, cds_fastas, [merged_coordinates, cds_coordinates], cds_counts, close_to_tip
