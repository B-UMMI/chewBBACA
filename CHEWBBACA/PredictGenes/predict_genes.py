#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module predicts genes from a set of input genome assemblies.

Code documentation
------------------
"""


import sys
import csv

import pandas as pd

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   fasta_operations as fao,
					   assembly_statistics as ast,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo,
					   pyrodigal_gene_prediction as pgp,
					   augustus_gene_prediction as agp)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 file_operations as fo,
								 fasta_operations as fao,
								 assembly_statistics as ast,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo,
								 pyrodigal_gene_prediction as pgp,
								 augustus_gene_prediction as agp)


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


def main(input_files, output_directory, gene_predictor, gene_predictor_parameters, cpu_cores):
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
	if gene_predictor == "pyrodigal":
		print(f'Using Pyrodigal to predict genes for {len(input_files)} inputs...')
		results = pgp.predict_genes(input_file_ids, output_directory, gene_predictor_parameters, cpu_cores)
	elif gene_predictor == "augustus":
		print(f'Using AUGUSTUS to predict genes for {len(input_files)} inputs...')
		results = agp.predict_genes(input_file_ids, output_directory, gene_predictor_parameters, cpu_cores)

	# Dictionary with info about inputs for which gene prediction failed
	# Total number of CDSs identified in the inputs
	# Paths to FASTA files with the extracted CDSs
	# Paths to files with the coordinates of the CDSs extracted for each input
	# Total number of CDSs identified per input
	# Dictionary with info about the CDSs closer to contig tips per input
	failed, total_extracted, cds_fastas, cds_coordinates, cds_counts, close_to_tip = results
	print()

	if len(failed) > 0:
		print(f'Failed to predict CDSs for {len(failed)} inputs.')
		print('Make sure that Pyrodigal runs in meta mode (--pm meta) if any input file has less than 100kbp.')
		# Write Pyrodigal stderr for inputs that failed gene prediction
		failed_lines = [f'{k}\t{v}' for k, v in failed.items()]
		failed_outfile = fo.join_paths(output_directory, ['gene_prediction_failures.tsv'])
		fo.write_lines(failed_lines, failed_outfile)

	if len(cds_fastas) == 0:
		sys.exit(f'{ct.CANNOT_PREDICT}')

	print(f'Predicted a total of {total_extracted} CDSs from {len(input_files)-len(failed)} inputs.')
	print(f'CDS coordinate data saved to {cds_coordinates}')

	# Add the CDS count to the assembly statistics
	for gid in assembly_statistics:
		# Add 0 if gene prediction failed or did not identify any genes
		assembly_statistics[gid].append(cds_counts.get(gid, 0))

	# Filter FASTA files based on confidence threshold
	grouped = {}
	if gene_predictor == "pyrodigal":
		confidence_threshold = gene_predictor_parameters["confidence_threshold"]
		if 'genes' in gene_predictor_parameters["output_formats"] and confidence_threshold:
			print(f'Identifying the CDSs with a confidence value < {confidence_threshold}...')
			# Create folder to store FASTA files containing CDSs above confidence threshold
			filtered_genes_outdir = fo.join_paths(output_directory, ['filtered_genes'])
			fo.create_directory(filtered_genes_outdir)

			# Read coordinate file to identify CDSs below confidence threshold
			with open(cds_coordinates, 'r') as infile:
				reader = csv.reader(infile, delimiter='\t')
				header = next(reader, None)
				# Get lines where the confidence value is below the threshold
				below = [line for line in reader if float(line[-2]) < confidence_threshold]

			print(f'Identified {len(below)} CDSs with a confidence value < {confidence_threshold}.')
			# Save lines matching excluded CDSs to separate TSV file
			below_outlines = ['\t'.join(line) for line in [header]+below]
			below_lines_outfile = fo.join_paths(output_directory, [ct.GENE_COORDINATES_EXCLUDED_BASENAME])
			fo.write_lines(below_outlines, below_lines_outfile)
			print(f'Wrote CDS coordinates for the CDSs excluded based on the confidence threshold to {below_lines_outfile}')

			print(f'Filtering FASTA records to create FASTA files containing only the {total_extracted-len(below)} CDSs with a confidence value >= {confidence_threshold}...')
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
