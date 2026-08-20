#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

Compute pairwise distances based on allele calling results.

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
	from utils import (constants as ct,
					   file_operations as fo,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
							  	 file_operations as fo,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def compute_hamming(current_row, permutation_rows, total_loci, similarity):
	"""Compute pairwise Hamming distances or similarities.
	
	Parameters
	----------
	current_row : ndarray
		Numpy array with dtype=int32 values for the allelic profile of the current genome.
	permutation_rows : ndarray
		Numpy array with dtype=int32 values for the allelic profiles of the genomes that will be compared against the current genome.
	total_loci : int
		Total number of loci in the allelic profiles.
	similarity : bool
		Compute similarity instead of distance.

	Returns
	-------
	ndarray
		Numpy array with dtype=int32 values for the pairwise distances or similarities.
	"""
	# Get number of shared alleles if similarity is True
	# Determine shared 0's
	# Determine each position the current row and permutation rows are both 0 and sum the number of shared 0's per row
	shared_zeros = ((current_row==0) & (permutation_rows==0)).sum(1)
	# Determine shared values, including 0's
	shared_values = np.count_nonzero(current_row==permutation_rows, axis=-1)
	# Determine non-zero shared values
	shared_alleles = shared_values - shared_zeros
	if similarity:
		shared_alleles = shared_alleles.astype('int32')
		return shared_alleles
	else:
		allelic_distances = total_loci - (shared_alleles+shared_zeros)
		allelic_distances = allelic_distances.astype('int32')
		return allelic_distances


def compute_jaccard(current_row, permutation_rows, total_loci, similarity):
	"""Compute pairwise Jaccard distances or similarities.
		
	Parameters
	----------
	current_row : ndarray
		Numpy array with dtype=int32 values for the allelic profile of the current genome.
	permutation_rows : ndarray
		Numpy array with dtype=int32 values for the allelic profiles of the genomes that will be compared against the current genome.
	total_loci : int
		Total number of loci in the allelic profiles.
	similarity : bool
		Compute similarity instead of distance.

	Returns
	-------
	ndarray
		Numpy array with dtype=int32 values for the pairwise distances or similarities.
	"""
	# Determine shared 0's
	shared_zeros = ((current_row==0) & (permutation_rows==0)).sum(1)
	# Determine shared values, including 0's
	shared_values = np.count_nonzero(current_row==permutation_rows, axis=-1)
	# Determine non-zero shared values
	shared_alleles = shared_values - shared_zeros

	# Compute Jaccard similarity
	# Subtract the number of shared zeros from the total number of loci to avoid biasing the similarity value with loci that are absent in both samples
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


def compute_different_loci(current_row, permutation_rows, total_loci, similarity):
	"""Compute the number of shared or not shared loci.
			
	Parameters
	----------
	current_row : ndarray
		Numpy array with dtype=int32 values for the allelic profile of the current genome.
	permutation_rows : ndarray
		Numpy array with dtype=int32 values for the allelic profiles of the genomes that will be compared against the current genome.
	total_loci : int
		Total number of loci in the allelic profiles.
	similarity : bool
		Compute the number of shared loci instead of the number of not shared loci.

	Returns
	-------
	ndarray
		Numpy array with dtype=int32 values for the number of shared or not shared loci.
	"""
	shared_loci = ((current_row!=0) & (permutation_rows!=0)).sum(1)
	if not similarity:
		shared_zeros = ((current_row==0) & (permutation_rows==0)).sum(1)
		not_shared = total_loci - (shared_loci+shared_zeros)
		not_shared = not_shared.astype('int32')
		return not_shared
	else:
		shared_loci = shared_loci.astype('int32')
		return shared_loci


def compute_distances(indexes, np_profiles, sample_ids, tmp_directory, method, similarity):
	""" Compute pairwise distances.

	Parameters
	----------
	indexes : list
		List with the line indexes of the allelic profiles to select from the input matrix and process.
	np_profiles : ndarray
		Numpy array with dtype=int32 values for the allelic profiles.
	sample_ids : list
		List with the sample identifiers.
	tmp_directory : str
		Path to a temporary directory where pickle files with intermediate results will be stored.
	method : str
		Type of distance to compute.
	similarity : bool
		Compute similarity values instead of distance values.

	Returns
	-------
	output_files : list
		List with the paths to all pickle files that were created to store intermediate results.
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
			distances = compute_hamming(current_row, permutation_rows, total_loci, similarity)
		elif method == 'jaccard':
			distances = compute_jaccard(current_row, permutation_rows, total_loci, similarity)
		elif method == 'loci':
			distances = compute_different_loci(current_row, permutation_rows, total_loci, similarity)

		# Save computed distances for current genome
		output_file = os.path.join(tmp_directory, current_genome)
		fo.pickle_dumper(distances, output_file)
		output_files[current_genome] = output_file

	return output_files


def write_matrix(pickled_results, genome_ids, output_directory, col_ids, output_type):
	"""Write matrix with the computed distances.

	Parameters
	----------
	pickled_results : dict
		Dictionary with sample identifiers as keys and paths to binary files with pickled results
		as values.
	genome_ids : list
		List with the sample identifiers.
	output_file : str
		Path to the output file to which the distance matrix will be saved.
	col_ids: list
		List with sample identifiers to add as headers.
	output_type: str
		Type of output matrix to create (upper_triangular, lower_triangular, symmetric).

	Returns
	-------
	output_file : str
		Path to the file with the distance matrix.
	"""
	# Default output is upper triangular matrix
	upper_triangular = os.path.join(output_directory, 'distances_upper.tsv')
	ad_lines = [col_ids]
	# Limit the number of lines to write at once to avoid memory overflow
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

	# User specified different output type
	if output_type != 'upper_triangular':
		# Transpose to get lower triangular matrix
		print('Transposing upper triangular matrix...')
		lower_triangular_outfile = fo.join_paths(output_directory, ['distances_lower.tsv'])
		lower_triangular = transpose_matrix(upper_triangular, lower_triangular_outfile)
		# Merge upper and lower triangular matrices to create a symmetric matrix
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


def transpose_matrix(input_file, output_file):
	"""Transpose lines in a TSV file.

	Parameters
	----------
	input_file : str
		Path to the input TSV file.
	output_file : str
		Path to the output file.

	Returns
	-------
	output_file : str
		Path to the file with the transposed matrix.
	"""
	# Divide into smaller sets (ncol=500) to avoid memory overflow when transposing huge files
	file_id = 1
	transpose_files = []
	with open(input_file, 'r') as infile:
		# Get columns names from first line
		columns = [e.strip() for e in (infile.__next__()).split('\t')]
		# Define smaller sets of columns to read and transpose
		num_col_sets = math.ceil(len(columns)/500)
		col_sets = im.divide_list_into_n_chunks(columns, num_col_sets)
		# Use Pandas to read columns sets and save transpose
		for c in col_sets:
			# dtype=str or Pandas converts values into floats
			df = pd.read_csv(input_file, usecols=c, delimiter='\t', dtype=str)
			# Transpose columns
			df = df.T
			# Do not save header that contains row indexes
			intermediate_file = fo.join_paths(os.path.dirname(output_file), [f'transpose{file_id}.tsv'])
			df.to_csv(intermediate_file, sep='\t', header=False)
			transpose_files.append(intermediate_file)
			file_id += 1

	# Concatenate all files with transposed lines
	fo.concatenate_files(transpose_files, output_file)
	# Delete intermediate files
	for file in transpose_files:
		os.remove(file)

	return output_file


def merge_triangular_matrices(upper_matrix, lower_matrix, output_directory, matrix_size):
	"""Merge a upper and lower triangular matrices to create a symmetric matrix.

	Parameters
	----------
	upper_matrix : str
		Path to the TSV file that contains the upper triangular matrix.
	lower_matrix : str
		Path to the TSV file that contains the lower triangular matrix.
	output_file : str
		Path to the output file to which the symmetric matrix will be saved.
	matrix_size : int
		Total number of lines in the triangular matrix.

	Returns
	-------
	output_file : str
		Path to the file with the symmetric matrix.
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
	"""Write a TSV file with distance values.

	Parameters
	----------
	pickled_results : dict
		Dictionary with sample identifiers as keys and paths to binary files with pickled results
		as values.
	genome_ids : list
		List with the sample identifiers.
	output_file : str
		Path to the output file to which the distances will be saved.

	Returns
	-------
	output_file : str
		Path to the output file with the distance values.
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
	"""Compute pairwise distances based on allele calling results.

	Parameters
	----------
	input_file : str
		Path to a TSV file containing allelic profiles determined by the AlleleCall module.
	output_directory : str
		Path to the output directory where the process will store intermediate and final results.
	method : str
		Distance method used to compute the distance matrix. The module supports the hamming, 
		jaccard and loci (number of loci not shared) methods.
	output_format : str
		Output format for the distance matrix (upper_triangular, lower_triangular, symmetric, table).
	no_mask : bool
		Do not mask missing data when computing the distance matrix. This option is useful when 
		the input profiles are already masked.
	similarity : bool
		Compute similarity values instead of distance values.
	cpu_cores : int
		Number of CPU cores used to compute distances.

	Returns
	-------
	output_file : str
		Path to the output file containing the distance or similarity values.
	"""
	# Create output directory
	fo.create_directory(output_directory)

	print(f'Distance method: {method}')
	print(f'Output format: {output_format}')
	# Import profiles
	print('Reading input file...')
	profiles = pd.read_csv(input_file, sep='\t', dtype=str, index_col=0)
	# Get sample identifiers
	sample_ids = profiles.index.tolist()
	total_samples = len(sample_ids)
	print(f'Total samples: {total_samples}')
	total_loci = profiles.shape[1]
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

	if similarity:
		print('Computing pairwise similarities...')
		tmp_directory_basename = 'pairwise_similarities'
	else:
		print('Computing pairwise distances...')
		tmp_directory_basename = 'pairwise_distances'

	# Create temp directory to store pairwise distances per genome
	tmp_directory = os.path.join(output_directory, tmp_directory_basename)
	fo.create_directory(tmp_directory)

	# Convert to numpy array with int32 values to reduce memory usage and make it easier to compute distances
	np_profiles = masked_profiles.to_numpy(dtype='int32')

	# Divide input into 20 lists for 5% progress resolution
	rows_indexes = [i for i in range(np_profiles.shape[0])]
	random.shuffle(rows_indexes)
	parallel_inputs = im.divide_list_into_n_chunks(rows_indexes, 20)

	# Create common arguments for parallel processing
	common_args = [[l, np_profiles, sample_ids, tmp_directory, method, similarity, compute_distances] for l in parallel_inputs]

	# Increasing cpu cores can greatly increase memory usage
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
		# Write TSV with one pairwise distance per line
		output_file = write_table(merged, sample_ids, output_directory)

	print(f'Results available in {output_directory}')

	# Delete folder with intermediate pickles
	fo.delete_directory(tmp_directory)

	return output_file
