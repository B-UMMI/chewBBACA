#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module predicts genes from a set of input genome assemblies.

Code documentation
------------------
"""


import os
import sys
import csv
import math

import pandas as pd

try:
	from utils import (constants as ct,
					   gene_prediction as gp,
					   file_operations as fo,
					   fasta_operations as fao,
					   assembly_statistics as ast,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 gene_prediction as gp,
								 file_operations as fo,
								 fasta_operations as fao,
								 assembly_statistics as ast,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def write_coordinates_file(coordinates_file, output_file):
	"""Write the coordinates of CDSs identified in a genome assembly to a TSV file.

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
	# Sort lines by genome and protein ID (converting to int ensures natural sort)
	lines.sort(key=lambda x: (x[2], int(x[6])))
	lines = ['\t'.join(line) for line in lines]
	fo.write_lines(lines, output_file)


def predict_genes(fasta_files, ptf_path, translation_table, pyrodigal_mode, cpu_cores,
				  output_directory, output_formats):
	"""Predict genes from a set of input genome assemblies.

	Parameters
	----------
	fasta_files : list
		List of paths to FASTA files with genomic
		sequences.
	ptf_path : str
		Path to the Pyrodigal training file. Should
		be NoneType if a training file is not provided.
	translation_table : int
		Genetic code used to predict and translate
		coding sequences.
	pyrodigal_mode : str
		Pyrodigal execution mode.
	cpu_cores : int
		Number of processes that will run Pyrodigal in
		parallel.
	output_directory : str
		Path to the directory where output files
		with Pyrodigal's results will be stored in.

	Returns
	-------
	failed_info : list
		List that contains a list with the stderr for the
		cases that Pyrodigal failed to predict genes for
		and the path to the file with information about
		failed cases. Returns NoneType if gene prediction
		succeeded for all inputs.
	"""
	if ptf_path is not None:
		# Read training file to create GeneFinder object
		training_data = gp.read_training_file(ptf_path)
		# Create GeneFinder object based on training data
		gene_finder = gp.create_gene_finder(training_data, True, True, False)
	elif ptf_path is None and pyrodigal_mode == 'meta':
		# Create GeneFinder object to run in meta mode
		gene_finder = gp.create_gene_finder(None, True, True, True)
	else:
		gene_finder = None

	pyrodigal_inputs = im.divide_list_into_n_chunks(fasta_files, len(fasta_files))
	# Add path to TSV file to store CDS coordinate data per input subset
	for i in range(len(pyrodigal_inputs)):
		pyrodigal_inputs[i].append(fo.join_paths(output_directory, [f'gene_coordinates_{i}']))

	# Create subfolders to store each output format
	for output_format in output_formats:
		format_outdir = fo.join_paths(output_directory, [output_format])
		fo.create_directory(format_outdir)

	# Add common arguments to all sublists
	common_args = [output_directory, gene_finder, translation_table, output_formats]
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

	# Merge TSV files with CDS coordinates
	coordinate_files = [i[1] for i in pyrodigal_inputs]
	merged_coordinates = fo.join_paths(output_directory, ['gene_coordinates.tsv'])
	fo.concatenate_files(coordinate_files, merged_coordinates, header=ct.GENE_TABLE_HEADER)
	# Delete separate coordinate files
	fo.remove_files(coordinate_files)

	# Merge dictionaries with info about CDSs close to contig tips
	close_to_tip = [line[2] for line in pyrodigal_results if len(line[2]) > 0]
	close_to_tip = im.merge_dictionaries(close_to_tip)

	return [failed, total_cds, cds_fastas, merged_coordinates, cds_counts, close_to_tip]


def exclude_records(input_file, input_id, records_ids, output_directory, delete_original=False):
	"""
	"""
	records = fao.sequence_generator(input_file)
	selected_records = []
	excluded_records = []
	for rec in records:
		if rec.id not in records_ids:
			selected_records.append(rec)
		else:
			excluded_records.append(rec)

	selected_outfile = fo.join_paths(output_directory, [f'{input_id}.fasta'])
	fao.write_records(selected_records, selected_outfile, mode='w')

	if delete_original:
		# Delete FASTA file with all CDSs to release disk space
		fo.remove_files([input_file])

	return excluded_records


def main(input_files, output_directory, pyrodigal_training_file, translation_table, pyrodigal_mode,
		 pyrodigal_output_formats, pyrodigal_minimum_confidence, cpu_cores):
	# Read file with paths to input files
	if isinstance(input_files, str):
		input_files = fo.read_lines(input_files, strip=True)
	# Passed list of file paths
	elif isinstance(input_files, list):
		pass

	# Sort file paths
	input_files = im.sort_iterable(input_files, sort_key=str.lower)
	print('Number of inputs: {0}'.format(len(input_files)))

	# Compute assembly statistics
	print('Computing assembly statistics...')
	assembly_statistics = ast.compute_assembly_stats(input_files, cpu_cores)
	print()

	# Map input file paths to file basename without extension and MD5 file hash
	input_file_ids = [(file, fo.file_basename(file, False)) for file in input_files]

	# Gene prediction step
	print(f'Predicting CDSs for {len(input_files)} inputs...')
	pyrodigal_results = predict_genes(input_file_ids, pyrodigal_training_file, translation_table,
									  pyrodigal_mode, cpu_cores, output_directory, pyrodigal_output_formats)

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

	print(f'Predicted a total of {total_extracted} CDSs from {len(input_files)-len(failed)} inputs.')
	print(f'CDS coordinate data saved to {cds_coordinates}')

	# Add the CDS count to the assembly statistics
	for gid in assembly_statistics:
		# Add 0 if gene prediction failed or did not identify any genes
		assembly_statistics[gid].append(cds_counts.get(gid, 0))

	if len(failed) > 0:
		# Write Pyrodigal stderr for inputs that failed gene prediction
		failed_lines = [f'{k}\t{v}' for k, v in failed.items()]
		failed_outfile = fo.join_paths(os.path.dirname(output_directory), ['gene_prediction_failures.tsv'])
		fo.write_lines(failed_lines, failed_outfile)

	# Filter FASTA files based on confidence threshold
	grouped = {}
	if 'genes' in pyrodigal_output_formats and pyrodigal_minimum_confidence:
		print(f'Identifying the CDSs with a confidence value < {pyrodigal_minimum_confidence}...')
		# Create folder to store FASTA files containing CDSs above confidence threshold
		filtered_genes_outdir = fo.join_paths(output_directory, ['filtered_genes'])
		fo.create_directory(filtered_genes_outdir)

		# Read coordinate file to identify CDSs below confidence threshold
		with open(cds_coordinates, 'r') as infile:
			reader = csv.reader(infile, delimiter='\t')
			header = next(reader, None)
			# Get lines where the confidence value is below the threshold
			below = [line for line in reader if float(line[-2]) < pyrodigal_minimum_confidence]

		print(f'Identified {len(below)} CDSs with a confidence value < {pyrodigal_minimum_confidence}.')
		# Save lines matching excluded CDSs to separate TSV file
		below_outlines = ['\t'.join(line) for line in [header]+below]
		below_lines_outfile = fo.join_paths(output_directory, [ct.GENE_COORDINATES_EXCLUDED_BASENAME])
		fo.write_lines(below_outlines, below_lines_outfile)
		print(f'Wrote CDS coordinates for the CDSs excluded based on the confidence threshold to {below_lines_outfile}')

		print(f'Filtering FASTA records to create FASTA files containing only the {total_extracted-len(below)} CDSs with a confidence value >= {pyrodigal_minimum_confidence}...')
		# Group CDS IDs based on genome of origin
		for line in below:
			grouped.setdefault(line[0].rsplit('_', 1)[0], []).append(line[0])

		filtering_inputs = []
		for gid in grouped:
			filtering_inputs.append([fo.join_paths(output_directory, ['genes', f'{gid}.fasta']), gid, grouped[gid], filtered_genes_outdir, True])

		# Add common arguments to all sublists
		filtering_inputs = im.multiprocessing_inputs(filtering_inputs, [], exclude_records)

		# Filter records
		filtering_results = mo.map_async_parallelizer(filtering_inputs,
													  mo.function_helper,
													  cpu_cores,
													  show_progress=True,
													  pool_type='threadpool')
		print()

		# Save excluded records from all inputs to the same FASTA file
		excluded_records = im.flatten_list(filtering_results)
		excluded_outfile = fo.join_paths(output_directory, ['excluded_genes.fasta'])
		fao.write_records(excluded_records, excluded_outfile, mode='w')

		print(f'FASTA files with filtered records saved to {filtered_genes_outdir}')
		print(f'FASTA file containing all excluded records saved to {excluded_outfile}')

		# Delete folder that contained the FASTA files with all CDSs
		fo.delete_directory(fo.join_paths(output_directory, ['genes']))

		# Update data returned by the function to account for the CDSs excluded based on the confidence threshold
		cds_fastas = fo.listdir_fullpath(filtered_genes_outdir)
		for gid in grouped:
			cds_counts[gid] = cds_counts[gid] - len(grouped[gid])

	# Add the number of excluded CDSs to the assembly statistics
	for gid in assembly_statistics:
		assembly_statistics[gid].append(len(grouped.get(gid, [])))

	# Save assembly statistics
	assembly_statistics_df = pd.DataFrame.from_dict(assembly_statistics, orient='index').reset_index()
	assembly_statistics_df.columns = ct.ASSEMBLY_STATS_COLUMNS
	assembly_statistics_outfile = fo.join_paths(output_directory, ['assembly_stats.tsv'])
	assembly_statistics_df.to_csv(assembly_statistics_outfile, sep='\t', index=False)
	print(f'Assembly statistics saved to {assembly_statistics_outfile}')

	print(f'Gene prediction results available in {output_directory}')

	return failed, total_extracted, cds_fastas, cds_coordinates, assembly_statistics_outfile, cds_counts, close_to_tip
