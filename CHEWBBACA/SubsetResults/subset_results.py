#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module subsets the data in files created by chewBBACA based on a list of loci and/or sample identifiers.

Code documentation
------------------
"""


import re
import sys

import pandas as pd

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   fasta_operations as fao)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
					   			 file_operations as fo,
								 fasta_operations as fao)


def main(input_directory, output_directory, loci_list, sample_list, inverse_loci, inverse_samples):
	"""Subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers

	Parameters
	----------
	input_directory : str
		Path to the directory containing the files to be subsetted.
	output_directory : str
		Path to the output directory.
	loci_list : str
		Path to a file with a list of loci to select, one identifier per line (without file extension).
	sample_list : str
		Path to a file with a list of samples to select, one identifier per line (without file extension).
	inverse_loci : bool
		If provided, the process will select the loci that are not in the input loci list.
	inverse_samples : bool
		If provided, the process will select the samples that are not in the input sample list.
	"""
	# Do not proceed if user provided neither a list of loci nor of samples
	if not loci_list and not sample_list:
		sys.exit(ct.SUBSETRESULTS_MISSING_LISTS_EXCEPTION)

	print(f'Listing files in {input_directory}...')

	# Need to list the files in each directory and group them based on filename
	file_list = {file_type: None for file_type in ct.ALLELECALL_OUTFILES}
	input_files = fo.listdir_fullpath(input_directory)
	for f in input_files:
		fbasename = fo.file_basename(f)
		if fbasename in file_list:
			file_list[fbasename] = f

	# Determine if there any file types missing
	missing = []
	for ftype, fpath in file_list.items():
		if fpath is None:
			print(f'Did not find matching files for the {ftype} file type. Output '
		 		  'directory will not include that file.')
			missing.append(ftype)
		else:
			print(f'Found the {ftype} file type to subset.')

	# Remove from the dictionary the file types without matching files
	for ftype in missing:
		del(file_list[ftype])

	if len(file_list) == 0:
		sys.exit(ct.SUBSETRESULTS_MISSING_FILES_EXCEPTION)
	elif 'results_alleles.tsv' not in file_list:
		sys.exit(ct.SUBSETRESULTS_MISSING_PROFILES_EXCEPTION)

	# Get list of loci and samples from the allelic profiles
	print(f'Getting the complete list of loci from {file_list['results_alleles.tsv']}...')
	loci = fo.read_lines(file_list['results_alleles.tsv'], strip=True, num_lines=1)
	loci = loci[0].split('\t')[1:]
	print(f'Input data includes results for {len(loci)} loci.')

	print(f'Getting the complete list of samples from {file_list['results_alleles.tsv']}...')
	samples = pd.read_csv(file_list['results_alleles.tsv'], usecols=['FILE'], sep='\t', dtype=str)
	samples = samples['FILE'].tolist()
	print(f'Input data includes results for {len(samples)} samples.')

	if loci_list:
		# Read loci list
		print(f'Getting the list of loci IDs to subset the results from {loci_list}...')
		# File containing loci list is a TSV file
		if loci_list.endswith('.tsv'):
			loci_subset = pd.read_csv(loci_list, sep='\t', dtype=str, index_col=0).index.tolist()
		# TXT file without header
		elif loci_list.endswith('.txt'):
			loci_subset = [line[0] for line in fo.read_tabular(loci_list)]

		print(f'Got {len(loci_subset)} loci IDs from {loci_list}')

		# Determine if any of the loci IDs in the list are not present in input results
		loci_to_keep = set(loci_subset).intersection(set(loci))
		if len(loci_to_keep) < len(loci_subset):
			missing_loci = set(loci_subset).difference(set(loci))
			print(f'{len(missing_loci)} loci IDs in the provided list are not '
		 		  f'present in the input results: {', '.join(missing_loci)}')
			sys.exit(ct.SUBSETRESULTS_ABSENT_LOCI_EXCEPTION)

		if inverse_loci is True:
			print('Provided `--inverse-loci`. Inverting loci list...')
			loci_to_keep = set(loci).difference(set(loci_subset))
	else:
		loci_to_keep = loci

	# Sort list of loci
	loci_to_keep = sorted(loci_to_keep)
	print(f'Will subset the results based on a list of {len(loci_to_keep)} loci IDs '
	   	  f'(from a total of {len(loci)} loci).')

	if sample_list:
		# Read samples list
		print(f'Getting the complete list of samples from {sample_list}...')
		# File containing loci list is a TSV file
		if sample_list.endswith('.tsv'):
			sample_subset = pd.read_csv(sample_list, sep='\t', dtype=str, index_col=0).index.tolist()
		# TXT file without header
		elif sample_list.endswith('.txt'):
			sample_subset = [line[0] for line in fo.read_tabular(sample_list)]

		print(f'Got {len(sample_subset)} sample IDs from {sample_list}')

		# Determine if any of the provided sample IDs in the list are not present in the input dataset
		samples_to_keep = set(sample_subset).intersection(set(samples))
		if len(samples_to_keep) < len(sample_subset):
			missing_samples = set(sample_subset).difference(set(samples))
			print(f'{len(missing_samples)} sample IDs in the provided list are not '
				  f'present in the input results: {', '.join(missing_samples)}')
			sys.exit(ct.SUBSETRESULTS_ABSENT_SAMPLES_EXCEPTION)

		if inverse_samples is True:
			print('Provided `--inverse-samples`. Inverting list of samples...')
			# It is easier to use a list of sample indexes to skip to subset the TSV files with profiles
			# So, we use that for that type of file and the actual loci subset for other file types
			# Need to +1 since Pandas uses 1-index
			samples_to_skip = [i+1 for i, sample in enumerate(samples) if sample in sample_subset]
		else:
			samples_to_skip = [i+1 for i, sample in enumerate(samples) if sample not in sample_subset]
	else:
		samples_to_skip = []
		samples_to_keep = samples

	# Sort list of samples
	samples_to_keep = sorted(samples_to_keep)
	print(f'Will subset the results based on a list of {len(samples_to_keep)} sample IDs '
		  f'(from a total of {len(samples)} samples).')

	# Create regex objects to filter based on loci and sample lists
	loci_finder = re.compile('|'.join(loci_to_keep))
	sample_finder = re.compile('|'.join(samples_to_keep))

	# Create output directory
	outdir_nonexists = fo.create_directory(output_directory)
	if outdir_nonexists is False:
		print('Output directory already exists. Will save subsetted data files to it anyway.')

	# Iterate over each file type to merge all
	for ftype, fpath in file_list.items():
		print(f'Subsetting {fpath}...')
		# Define name for the output file with the merged results
		output_file = fo.join_paths(output_directory, [ftype])
		# Subset TSV files with profiles
		if ftype in ['results_alleles.tsv', 'results_contigsInfo.tsv', 'presence_absence.tsv']:
			df = pd.read_csv(fpath, usecols=['FILE']+loci_to_keep, skiprows=samples_to_skip, sep='\t', dtype=str)
			# Save dataframe to file after filtering loci and samples
			df.to_csv(output_file, sep='\t', index=False)
		elif ftype in ['paralogous_counts.tsv']:
			# Get the allelic profiles for the loci and sample subsets
			profiles = pd.read_csv(file_list['results_alleles.tsv'], delimiter='\t', usecols=loci_to_keep, skiprows=samples_to_skip, dtype=str)
			paralogous_counts = (profiles=='NIPH').sum(axis=0) + (profiles=='NIPHEM').sum(axis=0)
			# Only keep entries for loci with a count > 0
			counts_df = pd.DataFrame(data={'Locus': loci_to_keep, 'Count': paralogous_counts})
			counts_df = counts_df[counts_df['Count']>0]
			counts_df.to_csv(output_file, index=False, sep='\t')
		elif ftype in ['missing_classes.fasta']:
			records = fao.sequence_generator(fpath)
			records = {rec.id: str(rec.seq) for rec in records}
			# Filter based on loci and sample subsets
			filtered_recids = [recid for recid in records if loci_finder.search(recid) and sample_finder.search(recid)]
			filtered_records = {recid: records[recid] for recid in filtered_recids}
			out_records = [f'>{recid}\n{seq}' for recid, seq in filtered_records.items()]
			fo.write_lines(out_records, output_file, joiner='\n')
		elif ftype in ['unclassified_sequences.fasta']:
			records = fao.sequence_generator(fpath)
			records = {rec.id: str(rec.seq) for rec in records}
			# Filter based on sample subset
			filtered_recids = [recid for recid in records if sample_finder.search(recid)]
			filtered_records = {recid: records[recid] for recid in filtered_recids}
			out_records = [f'>{recid}\n{seq}' for recid, seq in filtered_records.items()]
			fo.write_lines(out_records, output_file, joiner='\n')
		elif ftype in ['novel_alleles.fasta']:
			records = fao.sequence_generator(fpath)
			records = {rec.id: str(rec.seq) for rec in records}
			# Filter based on loci subset
			filtered_recids = [recid for recid in records if loci_finder.search(recid)]
			filtered_records = {recid: records[recid] for recid in filtered_recids}
			out_records = [f'>{recid}\n{seq}' for recid, seq in filtered_records.items()]
			fo.write_lines(out_records, output_file, joiner='\n')
		elif ftype in ['missing_classes.tsv', 'paralogous_loci.tsv']:
			# Get header and write it to output file
			header = fo.read_lines(fpath, num_lines=1)
			fo.write_lines(header, output_file)
			inlines = fo.read_lines(fpath)[1:]
			# Identify lines that include loci and samples in the subsets
			inlines = [line for line in inlines if loci_finder.search(line) and sample_finder.search(line)]
			fo.write_lines(inlines, output_file, joiner='\n', write_mode='a')
		elif ftype in ['cds_coordinates.tsv']:
			# Get header and write it to output file
			header = fo.read_lines(fpath, num_lines=1)
			fo.write_lines(header, output_file)
			inlines = fo.read_lines(fpath)[1:]
			# Filter based on samples subset
			if sample_list:
				inlines = [line for line in inlines if sample_finder.search(line)]
			fo.write_lines(inlines, output_file, joiner='\n', write_mode='a')
		elif ftype in ['invalid_cds.txt']:
			inlines = fo.read_lines(fpath)
			# Identify lines for samples in the subset
			if sample_list:
				inlines = [line for line in inlines if sample_finder.search(line)]
			fo.write_lines(inlines, output_file, joiner='\n')
		elif ftype in ['loci_summary_stats.tsv']:
			# Get the lines matching the loci subset
			statistics = pd.read_csv(fpath, delimiter='\t', index_col=0, dtype=str)
			# Filter rows to only get index entries matching the loci_to_keep
			statistics = statistics.filter(items=loci_to_keep, axis=0)
			# Reset index to add it as the first column
			statistics = statistics.reset_index().rename(columns={'index': 'Locus'})

			# Get the allelic profiles for the loci and sample subsets
			profiles = pd.read_csv(file_list['results_alleles.tsv'], delimiter='\t', usecols=loci_to_keep, skiprows=samples_to_skip, dtype=str)
			counts = {}

			# Count the number of INFs and EXC
			counts['EXC'] = (~profiles.isin(ct.ALLELECALL_CLASSIFICATIONS[2:])).sum(axis=0)
			counts['INF'] = profiles.apply(lambda col: col.str.startswith('INF')).sum(axis=0)
			counts['EXC'] = counts['EXC'] - counts['INF']
			for l in ct.ALLELECALL_CLASSIFICATIONS[2:]:
				counts[l] = (profiles==l).sum(axis=0)

			counts_df = pd.concat([s.reset_index(drop=True) for k, s in counts.items()], axis=1)
			counts_df.columns = ct.ALLELECALL_CLASSIFICATIONS

			# This sets the FILE column as Index
			counts_df = pd.concat([statistics['Locus'], counts_df], axis=1)
			# Compute the values for the Total CDSs Classified column
			# Cannot compute this value when subsetting a dataset because it is
			# not possible to know the number of total CDSs classified per NIPH/NIPHEM
			counts_df['Total CDSs Classified'] = counts_df[ct.ALLELECALL_CLASSIFICATIONS[:-1]].sum(axis=1) + counts_df[['NIPH', 'NIPHEM']].sum(axis=1)
			# Add '~' to indicate the values only approximate the real values
			#counts_df['Total_CDSs'] = '~' + counts_df['Total_CDSs'].astype(str)

			# Save merged files
			counts_df.to_csv(output_file, index=False, sep='\t')
		elif ftype in ['results_statistics.tsv']:
			# Get the allelic profiles per dataframe
			statistics = pd.read_csv(fpath, delimiter='\t', skiprows=samples_to_skip, dtype=str)
			profiles = pd.read_csv(file_list['results_alleles.tsv'], delimiter='\t', usecols=loci_to_keep, skiprows=samples_to_skip, dtype=str)

			counts = {}
			# Count the number of INFs and EXC
			counts['EXC'] = (~profiles.isin(ct.ALLELECALL_CLASSIFICATIONS[2:])).sum(axis=1)
			counts['INF'] = profiles.apply(lambda col: col.str.startswith('INF')).sum(axis=1)
			counts['EXC'] = counts['EXC'] - counts['INF']
			for l in ct.ALLELECALL_CLASSIFICATIONS[2:]:
				counts[l] = (profiles==l).sum(axis=1)

			counts_df = pd.concat([s.reset_index(drop=True) for k, s in counts.items()], axis=1)
			counts_df.columns = ct.ALLELECALL_CLASSIFICATIONS

			# This sets the FILE column as Index
			counts_df = pd.concat([statistics['FILE'], counts_df], axis=1)

			counts_df['Invalid CDSs'] = statistics['Invalid CDSs']
			# Cannot compute this value when subsetting a dataset because it is
			# not possible to know the number of CDSs classified per NIPH/NIPHEM
			counts_df['Total CDSs Classified'] = counts_df[ct.ALLELECALL_CLASSIFICATIONS[:-1]].sum(axis=1) + counts_df[['NIPH', 'NIPHEM']].sum(axis=1)
			# Add '~' to indicate the values only approximate the real values
			#counts_df['Classified_CDSs'] = '~' + counts_df['Classified_CDSs'].astype(str)
			counts_df['Total_CDSs'] = statistics['Total_CDSs']

			# Save merged files
			counts_df.to_csv(output_file, index=False, sep='\t')

		print(f'Saved subsetted {ftype} to {output_file}')

	print(f'Done. Subsetted results available in {output_directory}')
