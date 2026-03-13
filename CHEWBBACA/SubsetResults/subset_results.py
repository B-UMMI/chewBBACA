#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module subsets results generated with chewBBACA based on loci and sample lists.

Code documentation
------------------
"""


import sys
import csv

import pandas as pd

try:
	from utils import file_operations as fo
except ModuleNotFoundError:
	from CHEWBBACA.utils import file_operations as fo


def main(input_file, output_file, loci_list, sample_list, inverse_loci, inverse_samples):
	"""Subset results generated with chewBBACA.

	Parameters
	----------
	input_file : str
		Path to the input TSV file. The first column in the file must have sample
		identifiers and the remaining columns must have locus identifiers. Currently,
		the module can subset results from the `results_alleles.tsv` and `results_contigsInfo.tsv`
		files generated with the AlleleCall module and the `presence_absence.tsv` file generated 
		with the ExtractCgMLST module.
	output_file : str
		Path to the output file used to save the subsetted results.
	loci_list : str
		Path to a file with a list of loci to select, one identifier per line (without file extension).
	sample_list : str
		Path to a file with a list of samples to select, one identifier per line (without file extension).
	inverse_loci : bool
		If provided, the process will select the loci that are not in the input loci list.
	inverse_samples : bool
		If provided, the process will select the samples that are not in the input samples list.
	"""
	# Do not proceed if user provided neither a list of loci nor of samples
	if not loci_list and not sample_list:
		sys.exit('Did not provide a list of loci or samples. '
				 'Please provide at least one of those lists to subset results.')

	# Get list of loci in allele calling results
	print(f'Getting list of loci from {input_file}...')
	loci = fo.read_lines(input_file, strip=True, num_lines=1)
	loci = loci[0].split('\t')[1:]
	print('Total loci: {0}'.format(len(loci)))

	# Get list of samples in allele calling results
	print(f'Getting list of samples from {input_file}...')
	samples = pd.read_csv(input_file, usecols=['FILE'], sep='\t', dtype=str)
	samples = samples['FILE'].tolist()
	print('Total samples: {0}'.format(len(samples)))

	if loci_list:
		# Read loci list
		print(f'Reading list of loci ({loci_list}) to subset results...')
		with open(loci_list, 'r') as infile:
			loci_list_data = list(csv.reader(infile, delimiter='\t'))
			# Get only values in the first column to support TSV files with multiple columns
			loci_list_ids = [line[0] for line in loci_list_data]

		# Determine if any of the provided loci IDs in the list are not present in the input file
		common_loci = set(loci_list_ids).intersection(set(loci))
		if len(common_loci) < len(loci_list_ids):
			missing_loci = set(loci_list_ids).difference(set(loci))
			print(f'Warning: {len(missing_loci)} loci identifiers in the provided list are not present in the input file: {missing_loci}')
			print('Please check the provided list of loci identifiers and make sure all loci identifiers match those in the input file.')
			sys.exit(1)

		if inverse_loci is True:
			print('Provided --inverse-loci. Inverting loci list...')
			loci_to_keep = [locus for locus in loci if locus not in loci_list_ids]
		else:
			loci_to_keep = [locus for locus in loci if locus in loci_list_ids]
	else:
		loci_to_keep = loci

	print(f'Loci to keep in subsetted results: {len(loci_to_keep)}')

	if sample_list:
		# Read samples list
		print(f'Reading list of samples ({sample_list}) to subset results...')
		with open(sample_list, 'r') as infile:
			sample_list_data = list(csv.reader(infile, delimiter='\t'))
			# Get only values in the first column to support TSV files with multiple columns
			sample_list_ids = [line[0] for line in sample_list_data]

		# Determine if any of the provided sample IDs in the list are not present in the input file
		common_samples = set(sample_list_ids).intersection(set(samples))
		if len(common_samples) < len(sample_list_ids):
			missing_samples = set(sample_list_ids).difference(set(samples))
			print(f'Warning: {len(missing_samples)} sample identifiers in the provided list are not present in the input file: {missing_samples}')
			print('Please check the provided list of sample identifiers and make sure all sample identifiers match those in the input file.')
			sys.exit(1)

		if inverse_samples is True:
			print('Provided --inverse-samples. Inverting list of samples...')
			samples_to_skip = [i+1 for i, sample in enumerate(samples) if sample in sample_list_ids]
		else:
			samples_to_skip = [i+1 for i, sample in enumerate(samples) if sample not in sample_list_ids]
	else:
		samples_to_skip = []

	print(f'Samples to keep in subsetted results: {len(samples)-len(samples_to_skip)}')

	# Include first column with sample ids
	print(f'Subsetting results and saving to file ({output_file})...')
	loci_to_keep = ['FILE'] + loci_to_keep
	df = pd.read_csv(input_file, usecols=loci_to_keep, skiprows=samples_to_skip, sep='\t', dtype=str)
	# Save dataframe to file after filtering loci and samples
	df.to_csv(output_file, header=True, sep='\t', index=False)
	print('Done!')
