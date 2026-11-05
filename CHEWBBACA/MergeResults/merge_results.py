#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module merges results files created by chewBBACA. It can merge files 
with results for the same set of loci or identify the set of shared loci to 
create merged files containing the results only for the set of shared loci.

Code documentation
------------------
"""


from functools import reduce

import pandas as pd

try:
	import global_constants as gct
	from SubsetResults import subset_results
	from utils import file_operations as fo
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct
	from CHEWBBACA.SubsetResults import subset_results
	from CHEWBBACA.utils import file_operations as fo


def main(config):
	"""Merge results files created by chewBBACA.

	Parameters
	----------
	input_files : list
		List with paths to folders containing results files created by chewBBACA.
	output_diretory : str
		Path to the output directory.
	common : bool
		If the process should join results data only for the set of shared
		loci when the results do not share the same loci sets.
	"""
	print(f'Received paths to {len(config.results_directories)} folders:\n{"\n".join(config.results_directories)}')
	# Work with the set of common loci if --common was provided
	if config.common[0] and len(config.common[0]) != config.common[2]:
		print(f'The input data share results for {len(config.common[1])} loci. Will subset the results in each input directory before merging.')
		# Save list of loci to pass to the SubsetResults module
		loci_list_outfile = fo.join_paths(config.output_directory, ['loci.txt'])
		fo.write_lines(config.common[1], loci_list_outfile)
		subsetted_inputs = []
		for d in config.results_directories:
			print(f'Calling the SubsetResults module to subset the results in {d}...')
			subset_outdir = fo.join_paths(config.output_directory, [fo.file_basename(d)])
			subset_results.main(d, subset_outdir, loci_list_outfile, None, None, None)
			subsetted_inputs.append(subset_outdir)
		input_files = subsetted_inputs
		# Delete the temporary loci list file
		fo.remove_files([loci_list_outfile])
	else:
		input_files = config.results_directories

	# Need to list the files in each directory and group them based on filename
	file_list = {ftype: [] for ftype in gct.SUPPORTED_RESULTS_FILES}
	for d in input_files:
		# List all files in the folder
		files = fo.listdir(d, full_path=True)
		for f in files:
			fbasename = fo.file_basename(f)
			if fbasename in file_list:
				file_list[fbasename].append(f)

	# Do not merge if a file type was not in all input folders
	# or if due to some unknown reason there are more files than expected
	to_exclude = []
	for ftype, paths in file_list.items():
		if len(paths) == 0:
			print(f'Did not find matching files for the {ftype} file type.')
			to_exclude.append(ftype)
		elif len(paths) != len(input_files):
			# These file types are expected to have the same number of files as input folders.
			# If absent, merging would just create incomplete results.
			if ftype in ['results_alleles.tsv', 'results_contigsInfo.tsv',
						 'presence_absence.tsv', 'loci_summary_stats.tsv',
						 'results_statistics.tsv', 'cds_coordinates.tsv']:
				print(f'Will not merge files ending in {ftype} because their '
					  f'number ({len(paths)}) is different than expected ({len(input_files)}).')
				to_exclude.append(ftype)
			# These file types may not always be present depending on the options used to run chewBBACA and merging them when some are missing is not critical.
			elif ftype in ['paralogous_counts.tsv', 'paralogous_loci.tsv',
				  		   'missing_classes.tsv', 'missing_classes.fasta',
						   'novel_alleles.fasta', 'invalid_cds.txt',
						   'unclassified_sequences.fasta']:
				print(f'Found {len(paths)} files for the {ftype} file type.')
		else:
			print(f'Found {len(paths)} files for the {ftype} file type.')

	# Remove from the dictionary the file types without matching files or an unexpected number of files
	for ftype in to_exclude:
		del(file_list[ftype])

	print(f'Files in the input directories share results for {len(config.common[1])} loci.')
	print(f'Merging results for each file type.')

	# Iterate over each file type to merge all
	for ftype, fpaths in file_list.items():
		# Define name for the output file with the merged results
		output_file = fo.join_paths(config.output_directory, [ftype])
		print(f'Merging {len(fpaths)} files for the {ftype} file type...')
		# Merge TSV files with profiles
		if ftype in ['results_alleles.tsv', 'results_contigsInfo.tsv', 'presence_absence.tsv']:
			# Merge all files without duplicating the header
			# Add the FILE column to read it as first column
			# Pass the full or common set of loci
			dataframes = [pd.read_csv(p, delimiter='\t', dtype=str, usecols=['FILE']+config.common[1]) for p in fpaths]
			merged_df = pd.concat(dataframes, ignore_index=True)
			# Save merged files
			merged_df.to_csv(output_file, index=False, sep='\t')
		elif ftype in ['loci_summary_stats.tsv', 'paralogous_counts.tsv']:
			# Read and use FILE column as index to filter based on loci sets if necessary
			dataframes = [pd.read_csv(p, delimiter='\t', index_col=0) for p in fpaths]
			# Sum all values to get merged dataframe
			merged_df = reduce(lambda a, b: a.add(b, fill_value=0), dataframes)
			# Reset the index to set Index as the first column
			# Need to rename the index/first column to set the original name
			merged_df = merged_df.reset_index().rename(columns={'index': 'Locus'})
			# Save merged files
			merged_df.to_csv(output_file, index=False, sep='\t')
		elif ftype in ['results_statistics.tsv', 'missing_classes.tsv', 'paralogous_loci.tsv', 'cds_coordinates.tsv']:
			# Get header and write it to output file
			header = fo.read_lines(fpaths[0], num_lines=1)
			fo.write_lines(header, output_file)
			for p in fpaths:
				inlines = fo.read_lines(p)[1:]
				fo.write_lines(inlines, output_file, joiner='\n', write_mode='a')
		elif ftype in ['missing_classes.fasta', 'novel_alleles.fasta', 'invalid_cds.txt', 'unclassified_sequences.fasta']:
			fo.concatenate_files(fpaths, output_file)

		print(f'Saved merged results to {output_file}')

	# Delete the temporary files created to subset the results if --common was provided
	if config.common[0] and len(config.common[0]) != config.common[2]:
		for d in input_files:
			fo.delete_directory(d)

	print(f'Done. Merged results available in {config.output_directory}')
