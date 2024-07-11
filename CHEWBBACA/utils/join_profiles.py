#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module joins allele calling results from different
runs. It can concatenate files with allelic profiles
for the same set of loci or create a new file with the
allelic profiles for the loci shared by all input files.

Code documentation
------------------
"""


import os
import sys

import pandas as pd

try:
	from utils import file_operations as fo
except ModuleNotFoundError:
	from CHEWBBACA.utils import file_operations as fo


def concatenate_profiles(files, loci_list, output_file):
	"""Concatenate allele calling results for a common set of loci.

	Parameters
	----------
	files : list
		List with the paths to the TSV files with
		allele calling results.
	loci_list : list
		List with the identifiers of common loci.
	output_file : str
		Path to the output file.

	Returns
	-------
	total_profiles : int
		Number of profiles written to the output file.
	"""
	total_profiles = 0
	for f in files:
		# create dataframe with columns for loci in list
		profiles = pd.read_csv(f, sep='\t',
							   usecols=loci_list, dtype=str)
		# reorder columns to prevent concatenating results
		# with same loci but with different loci/column order
		profiles = profiles[loci_list]
		total_profiles += len(profiles)
		# save dataframe to file
		profiles.to_csv(output_file, mode='a',
						header=not os.path.exists(output_file),
						sep='\t', index=False)

	return total_profiles


def main(profiles, output_file, common):
	"""Join files with allelic profiles.

	Parameters
	----------
	profiles : list
		List with paths to TSV files with allelic profiles.
	output_file : str
		Path to the output file.
	common : bool
		If the process should join profile data only for shared loci
		when the profiles do not share the same loci sets.
	"""
	if len(profiles) == 1:
		sys.exit('Provided a single file. Nothing to do.')

	headers = []
	for file in profiles:
		header = fo.read_lines(file, strip=True, num_lines=1)
		headers.append(header[0].split('\t'))

	if common is False:
		# check if headers are equal
		if all([set(headers[0]) == set(h) for h in headers[1:]]) is True:
			print('Profiles have {0} loci.'.format(len(headers[0])-1))
			total_profiles = concatenate_profiles(profiles,
												  headers[0],
												  output_file)
		else:
			sys.exit('Files have different sets of loci. Please provide '
					 'files with the results for the same set of loci or '
					 'provide the "--common" parameter to create a file '
					 'with the results for the set of common loci.')
	else:
		# determine set of common loci
		common_loci = headers[0]
		for h in headers[1:]:
			# determine common loci ordered based on headers in first file
			latest_common = [locus for locus in common_loci if locus in h]
			# update set of common loci so far
			common_loci = latest_common

		if len(latest_common) <= 1:
			sys.exit('Profiles do not have loci in common.')

		print('Profiles have {0} loci in common.'
			  ''.format(len(latest_common)-1))
		total_profiles = concatenate_profiles(profiles,
											  latest_common,
											  output_file)

	print('Joined {0} files with a total of {1} '
		  'profiles.'.format(len(profiles), total_profiles))
