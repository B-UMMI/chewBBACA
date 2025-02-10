#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

Determines the pairwise allelic differences based on a TSV file
with allelic profiles determined by the AlleleCall module to
create a distance matrix. The 'INF-' prefix is removed and ASM,
ALM, NIPH, NIPHEM, PLOT3, PLOT5, LNF and LOTSC classifications
are substituted by '0' before computing the pairwise distances.

Code documentation
------------------
"""


import os
import csv
import math

import random

import numpy as np
import pandas as pd

try:
	from utils import (
		constants as ct,
		file_operations as fo,
		iterables_manipulation as im,
		multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (
		constants as ct,
		file_operations as fo,
		iterables_manipulation as im,
		multiprocessing_operations as mo)


def compute_hamming(current_row, permutation_rows, similarity):
	""" Compute pairwise Hamming distances.
	"""
	# multiply 1D-array per whole matrix
	# all non-shared loci will be converted to 0
	# values different than 0 correspond to shared loci
	multiplied = current_row * permutation_rows
	# subtraction will lead to values different than 0 for loci that have different alleles
	# multiplying ensures that we only keep results for shared loci and not for
	# loci that are not shared and that had value different than 0 from subtraction
	allelic_distances = np.count_nonzero(multiplied * (current_row - permutation_rows), axis=-1)
	# Get number of shared alleles if similarity is True
	if similarity:
		shared_alleles = current_row.shape[1] - allelic_distances
		shared_alleles = shared_alleles.astype('int32')
		return shared_alleles
	else:
		allelic_distances = allelic_distances.astype('int32')
		return allelic_distances


def compute_jaccard(current_row, permutation_rows, total_loci, similarity):
	""" Compute pairwise Jaccard distances.
	"""
	# Determine shared 0's
	shared_zeros = ((current_row==0) & (permutation_rows==0)).sum(1)
	# Determine shared values, including 0's
	shared_values = np.count_nonzero(current_row==permutation_rows, axis=-1)
	# Determined non-zero shared values
	shared_alleles = shared_values - shared_zeros

	# Compute Jaccard similarity
	jaccard_similarity = shared_alleles / (total_loci-shared_zeros)
	if not similarity:
		# Compute Jaccard distance
		jaccard_distances = 1 - jaccard_similarity
		# Change values to dtype=float16
		# Reduces output string size and is precise enough
		jaccard_distances = jaccard_distances.astype('float16')
		return jaccard_distances
	else:
		jaccard_similarity = jaccard_similarity.astype('float16')
		return jaccard_similarity


def compute_different_loci(current_row, permutation_rows, similarity):
	""" Compute number of loci not shared.
	"""
	not_shared = ((current_row==0) & (permutation_rows!=0)).sum(1)
	if similarity:
		shared_loci = current_row.shape[1] - not_shared
		shared_loci = shared_loci.astype('int32')
		return shared_loci
	else:
		not_shared = not_shared.astype('int32')
		return not_shared


def compute_distances(indexes, np_profiles, sample_ids, tmp_directory, method, similarity):
	""" Compute pairwise distances.

	Parameters
	----------
	indexes : list
		List with the line index of the allelic profiles
		that will be processed.
	np_profiles : ndarray
		Numpy array with dtype=int32 values for allelic profiles.
	sample_ids : list
		List with sample identifiers.
	tmp_directory : str
		Path to temporary directory where pickle files with
		results will be stored.
	method : str
		Type of distance to compute.
	similarity : bool
		Compute similarity instead of distance.

	Returns
	-------
	output_files : list
		List with the paths to all pickle files that were created
		to store results.
	"""
	# Get total number of loci
	total_loci = np_profiles.shape[1]

	# Multiply one row per cycle to avoid memory overflow
	output_files = {}
	for i in indexes:
		current_genome = sample_ids[i]
		# Get one row to perform pairwise comparisons against whole matrix
		current_row = np_profiles[i:i+1, :]
		permutation_rows = np_profiles[i:, :]

		if method == 'hamming':
			distances = compute_hamming(current_row, permutation_rows, similarity)
		elif method == 'jaccard':
			distances = compute_jaccard(current_row, permutation_rows, total_loci, similarity)
		elif method == 'loci':
			distances = compute_different_loci(current_row, permutation_rows, similarity)

		# Save computed distances for current genome
		output_file = os.path.join(tmp_directory, current_genome)
		fo.pickle_dumper(distances, output_file)
		output_files[current_genome] = output_file

	return output_files


def write_matrix(pickled_results, genome_ids, output_directory, col_ids, output_type):
	"""Write upper triangular distance matrix.

	Parameters
	----------
	pickled_results : dict
		Dictionary with sample identifiers as keys
		and paths to binary files with pickled results
		as values.
	genome_ids : list
		List with sample identifiers.
	output_file : str
		Path to the output file to which the distance
		matrix will be saved.
	col_ids: list
		List with sample identifiers to add as headers.

	Returns
	-------
	True if there are no errors.
	"""
	upper_triangular = os.path.join(output_directory, 'distances_upper.tsv')
	ad_lines = [col_ids]
	limit = 300
	# Create based on genome order in input matrix
	for g in genome_ids:
		current_file = pickled_results[g]
		# Load data
		data = fo.pickle_loader(current_file)
		allele_diffs = list(data)
		allele_diffs = list(map(str, allele_diffs))
		# Add padding before distance values
		# This creates a right/upper triangular matrix
		padding = [''] * (len(genome_ids)-len(allele_diffs))
		ad_line = [g] + padding + allele_diffs
		ad_lines.append(ad_line)
		if len(ad_lines) >= limit or g == genome_ids[-1]:
			ad_text = [im.join_list(l, '\t') for l in ad_lines]
			fo.write_lines(ad_text, upper_triangular, write_mode='a')
			ad_lines = []

	if output_type != 'upper_triangular':
		print('Transposing upper triangular matrix...')
		lower_triangular = transpose_matrix(upper_triangular, output_directory)
		if output_type == 'symmetric':
			print('Creating symmetric matrix...')
			output_file = merge_triangular_matrices(upper_triangular, lower_triangular, output_directory, len(col_ids))
			os.remove(lower_triangular)
		else:
			output_file = lower_triangular
		os.remove(upper_triangular)
	else:
		output_file = upper_triangular

	return output_file


def transpose_matrix(input_file, output_directory):
	"""Transpose lines in a TSV file.

	Parameters
	----------
	input_file : str
		Path to the input TSV file.
	output_directory : str
		Path to the directory to which intermediate files
		with transposed lines will be written.

	Returns
	-------
	output_transpose : str
		Path to the file with the transposed matrix.
		This file is created by concatenating all
		files saved into `output_directory`.
	"""
	file_id = 1
	transpose_files = []
	input_basename = os.path.basename(input_file)
	with open(input_file, 'r') as infile:
		# get columns names
		columns = [e.strip() for e in (infile.__next__()).split('\t')]
		# divide into smaller sets to avoid loading huge files
		num_col_sets = math.ceil(len(columns)/500)
		col_sets = im.divide_list_into_n_chunks(columns, num_col_sets)
		# use Pandas to read columns sets and save transpose
		for c in col_sets:
			# dtype=str or Pandas converts values into floats
			df = pd.read_csv(input_file, usecols=c, delimiter='\t', dtype=str)
			output_basename = input_basename.replace('.tsv', '_{0}.tsv'.format(file_id))
			output_file = os.path.join(output_directory, output_basename)
			# transpose columns
			df = df.T
			# do not save header that contains row indexes
			df.to_csv(output_file, sep='\t', header=False)
			transpose_files.append(output_file)
			file_id += 1

	# concatenate all files with transposed lines
	output_transpose = input_file.replace('.tsv', '_transpose.tsv')
	fo.concatenate_files(transpose_files, output_transpose)
	# Delete intermediate files
	for file in transpose_files:
		os.remove(file)

	return output_transpose


def merge_triangular_matrices(upper_matrix, lower_matrix, output_directory, matrix_size):
	"""Merge two triangular matrices to create a symmetric matrix.

	Parameters
	----------
	upper_matrix : str
		Path to the TSV file that contains the upper
		triangular matrix.
	lower_matrix : str
		Path to the TSV file that contains the lower
		triangular matrix.
	output_file : str
		Path to the output file to which the symmetric
		matrix will be saved.
	matrix_size : int
		Total number of lines in the triangular matrix.

	Returns
	-------
	None.
	"""
	output_file = os.path.join(output_directory, 'distances_symmetric.tsv')
	with open(upper_matrix, 'r') as upper_handle, open(lower_matrix, 'r') as lower_handle:
		upper_reader = csv.reader(upper_handle, delimiter='\t')
		lower_reader = csv.reader(lower_handle, delimiter='\t')
		merged_lines = []
		for i in range(matrix_size):
			upper_line = upper_reader.__next__()
			lower_line = lower_reader.__next__()
			merged_line = lower_line[0:i] + upper_line[i:]
			merged_lines.append(merged_line)
			if len(merged_lines) >= 200 or i == (matrix_size-1):
				ad_text = [im.join_list(l, '\t') for l in merged_lines]
				fo.write_lines(ad_text, output_file, write_mode='a')
				merged_lines = []

	return output_file


def write_table(pickled_results, genome_ids, output_directory):
	"""Write TSV file with distance values.

	Parameters
	----------
	pickled_results : dict
		Dictionary with sample identifiers as keys
		and paths to binary files with pickled results
		as values.
	genome_ids : list
		List with sample identifiers.
	output_file : str
		Path to the output file to which the distance
		matrix will be saved.

	Returns
	-------
	True if there are no errors.
	"""
	output_file = os.path.join(output_directory, 'distances.tsv')
	ad_lines = []
	limit = 100000
	# Create based on genome order in input matrix
	for g in genome_ids:
		current_file = pickled_results[g]
		# Load data
		data = fo.pickle_loader(current_file)
		allele_diffs = list(data)
		allele_diffs = list(map(str, allele_diffs))
		# Determine sample index
		sample_index = len(genome_ids)-len(allele_diffs)
		sample_lines = [[g, genome_ids[i], allele_diffs[i-sample_index]] for i in range(sample_index, len(genome_ids))]
		ad_lines.extend(sample_lines)

		if len(ad_lines) >= limit or g == genome_ids[-1]:
			ad_text = [im.join_list(l, '\t') for l in ad_lines]
			fo.write_lines(ad_text, output_file, write_mode='a')
			ad_lines = []

	return output_file


def main(input_file, output_directory, method, output_format, no_mask, similarity, cpu_cores):
	"""Compute a distance matrix based on allelic profiles.

	Parameters
	----------
	input_file : str
		Path to a TSV file with allelic profiles determined by
		the AlleleCall module.
	output_directory : str
		Path to the output directory.
	method : str
		Method to compute distances.
	output_format : str
		Output format for the output file.
	no_mask : bool
		True if the input matrix values are masked, False otherwise.
		The process will mask the matrix values it this value is False.
	similarity : bool
		Compute similarity instead of distance.
	cpu_cores : int
		Number of CPU cores used to compute distances.
	"""
	# Create output directory
	fo.create_directory(output_directory)

	# Determine input basename
	input_basename = fo.file_basename(input_file, False)

	# Import matrix
	print('Reading input file...')
	profiles = pd.read_csv(input_file, sep='\t', dtype=str, index_col=0)
	# Get sample identifiers
	sample_ids = profiles.index.tolist()
	total_samples = len(sample_ids)
	print(f'Total samples: {total_samples}')
	total_loci = profiles.shape[1] - 1
	print(f'Total loci: {total_loci}')

	# Mask matrix values
	if no_mask is False:
		output_masked = os.path.join(output_directory, ct.MASKED_PROFILES_BASENAME)
		# Mask special classifications and remove 'INF-' prefixes
		print('Masking profiles...')
		masked_profiles = profiles.apply(im.replace_chars)
		# Save masked profiles
		masked_profiles.to_csv(output_masked, sep='\t')
		print(f'Masked {total_samples} profiles.')
	else:
		masked_profiles = profiles

	# Create temp directory to store pairwise distances per genome
	tmp_directory = os.path.join(output_directory, 'pairwise_distances')
	fo.create_directory(tmp_directory)

	# Drop column with sample identifiers
	masked_profiles = masked_profiles.drop(columns=[masked_profiles.columns[0]])
	# Convert to numpy array with int32 values
	np_profiles = masked_profiles.to_numpy(dtype='int32')

	# Divide input into 20 lists for 5% progress resolution
	rows_indexes = [i for i in range(np_profiles.shape[0])]
	random.shuffle(rows_indexes)
	parallel_inputs = im.divide_list_into_n_chunks(rows_indexes, 20)

	# Create common arguments for parallel processing
	common_args = [[l, np_profiles, sample_ids, tmp_directory, method, similarity, compute_distances] for l in parallel_inputs]
	print(f'Divided inputs to process into {len(parallel_inputs)} tasks...')

	# Increasing cpu cores can greatly increase memory usage
	print('Computing pairwise distances...')
	results = mo.map_async_parallelizer(common_args,
										mo.function_helper,
										cpu_cores,
										show_progress=True)
	print()

	merged = im.merge_dictionaries(results)

	if output_format != 'table':
		print('Creating output matrix...')
		# Import arrays per genome and save to matrix file
		col_ids = ['FILE'] + sample_ids
		output_file = write_matrix(merged, sample_ids, output_directory, col_ids, output_format)
	elif output_format == 'table':
		print('Creating output table...')
		# Wrtie TSV with one pairwise distance per line
		output_file = write_table(merged, sample_ids, output_directory)

	print('Results available in {0}'.format(output_directory))

	# Delete folder with intermediate pickles
	fo.delete_directory(tmp_directory)

	return output_file
