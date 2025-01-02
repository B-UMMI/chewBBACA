#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module subsets allele calling results based on loci and sample lists.

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


def main(input_file, output_file, loci_list, samples_list, inverse_loci, inverse_samples):
	"""Subset allele calling results.

	Parameters
	----------
	input_file : str
		Path to a TSV file with allelic profiles.
	output_file : str
		Path to the output file.
	loci_list : str
		Path to a file with a list of loci to select, one identifier per line.
	loci_list : str
		Path to a file with a list of samples to select, one identifier per line.
	inverse_loci : bool
		If provided, the process will select the loci that are not in the input loci list.
	inverse_samples : bool
		If provided, the process will select the samples that are not in the input samples list.
	"""
	# Do not proceed if user provided neither a list of loci nor of samples
	if not loci_list and not samples_list:
		sys.exit('Did not provide a list of loci or of samples. Please provide at least one of those lists.')

	# Get list of loci in allele calling results
	print('Getting list of loci from input file...')
	loci = fo.read_lines(input_file, strip=True, num_lines=1)
	loci = loci[0].split('\t')[1:]
	print('Total loci: {0}'.format(len(loci)))

	# Get list of samples in allele calling results
	print('Getting list of samples from input file...')
	samples = pd.read_csv(input_file, usecols=['FILE'], sep='\t', dtype=str)
	samples = samples['FILE'].tolist()
	print('Total samples: {0}'.format(len(samples)))

	if loci_list:
		# Read loci list
		print('Reading list of loci to keep...')
		with open(loci_list, 'r') as infile:
			loci_list_data = list(csv.reader(infile, delimiter='\t'))
			# Get only values in the first column to support TSV files with multiple columns
			loci_list_ids = [line[0] for line in loci_list_data]

		if inverse_loci is True:
			print('Provided --inverse-loci. Inverting loci list...')
			loci_to_keep = [locus for locus in loci if locus not in loci_list_ids]
		else:
			loci_to_keep = [locus for locus in loci if locus in loci_list_ids]
	else:
		loci_to_keep = loci

	print('Loci to keep: {0}'.format(len(loci_to_keep)))

	if samples_list:
		# Read samples list
		print('Reading list of samples to keep...')
		with open(samples_list, 'r') as infile:
			samples_list_data = list(csv.reader(infile, delimiter='\t'))
			# Get only values in the first column to support TSV files with multiple columns
			samples_list_ids = [line[0] for line in samples_list_data]

		if inverse_samples is True:
			print('Provided --inverse-samples. Inverting samples list...')
			samples_to_skip = [i+1 for i, sample in enumerate(samples) if sample in samples_list_ids]
		else:
			samples_to_skip = [i+1 for i, sample in enumerate(samples) if sample not in samples_list_ids]
	else:
		samples_to_skip = []

	print('Samples to keep: {0}'.format(len(samples)-len(samples_to_skip)))

	# Include first column with sample ids
	loci_to_keep = ['FILE'] + loci_to_keep
	df = pd.read_csv(input_file, usecols=loci_to_keep, skiprows=samples_to_skip, sep='\t', dtype=str)
	
	# Save dataframe to file after filtering loci and samples
	df.to_csv(output_file, header=True, sep='\t', index=False)
