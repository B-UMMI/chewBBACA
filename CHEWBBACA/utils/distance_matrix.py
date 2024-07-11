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


def tsv_to_nparray(input_file, array_dtype='int32'):
	"""Read matrix of allelic profiles and convert to Numpy array.

	File lines are read as a generator to exclude the sample identifier
	and avoid loading huge files into memory.

	Parameters
	----------
	input_file : str
		Path to the TSV file with the matrix with
		allelic profiles.
	array_dtype : str
		Array data type.

	Returns
	-------
	np_array : ndarray
		Numpy array with the numeric values for
		all allelic profiles.
	"""
	# import matrix without column and row identifiers
	with open(input_file, 'r') as infile:
		lines = ('\t'.join(line.split('\t')[1:])
				 for line in infile)
		# dtype=float32 should be faster than integer dtypes
		# but runs faster with dtype=int32 in test setup
		# dtype=int32 supports max integer of 2,147,483,647
		# should be safe even when arrays are multiplied
		np_array = np.genfromtxt(fname=lines, delimiter='\t',
								 dtype=array_dtype, skip_header=1)

		# need to reshape array to 2D if it only has info for one sample
		try:
			# if it does not have columns it needs to be reshaped
			np_array.shape[1]
		except Exception:
			np_array = np.array([np_array])

	return np_array


def compute_distances(indexes, np_matrix, genome_ids, tmp_directory):
	"""Compute pairwise allelic differences and number of shared loci.

	Parameters
	----------
	indexes : list
		List with the line index of the allelic profiles
		that will be processed.
	np_matrix : ndarray
		Numpy array with dtype=int32 values for allelic profiles.
	genome_ids : list
		List with sample identifiers.
	tmp_directory : str
		Path to temporary directory where pickle files with
		results will be stored.

	Returns
	-------
	output_files : list
		List with the paths to all pickle files that were created
		to store results.
	"""
	# multiply one row per cycle to avoid memory overflow
	# read only part of the matrix for huge files and process in chunks?
	output_files = {}
	for i in indexes:
		current_genome = genome_ids[i]
		# get one row to perform pairwise comparisons against whole matrix
		current_row = np_matrix[i:i+1, :]
		# do not multiply against rows that were multiplied against
		# matrix's rows in previous iterations
		# combinations instead of permutations
		permutation_rows = np_matrix[i:, :]

		# multiply 1D-array per whole matrix
		# all non-shared loci will be converted to 0
		# values different than 0 correspond to shared loci
		multiplied = current_row * permutation_rows
		# count number of shared loci, non-zero values
		# pairwise_shared_loci = np.count_nonzero(multiplied, axis=-1)

		# subtraction will lead to values different than 0 for loci that have different alleles
		# multiplying ensures that we only keep results for shared loci and not for
		# loci that are not shared and that had value different than 0 from subtraction
		pairwise_allelic_differences = np.count_nonzero(multiplied * (current_row - permutation_rows), axis=-1)

		output_file = os.path.join(tmp_directory, current_genome)
		# fo.pickle_dumper([pairwise_shared_loci, pairwise_allelic_differences], output_file)
		fo.pickle_dumper([pairwise_allelic_differences], output_file)
		output_files[current_genome] = output_file

	return output_files


def get_sample_ids(input_file, delimiter='\t'):
	r"""Extract the sample identifiers from a matrix with allelic profiles.

	Parameters
	----------
	input_file : str
		Path to the input file that contains a matrix
		with allelic profiles.
	delimiter : str, optional
		Field delimiter. The default is '\t'.

	Returns
	-------
	sample_ids : list
		List with the sample identifiers.
	"""
	with open(input_file, 'r') as infile:
		reader = csv.reader(infile, delimiter=delimiter)
		sample_ids = [line[0] for line in reader][1:]

	return sample_ids


# def write_matrices(pickled_results, genome_ids, output_pairwise,
#                    output_p, col_ids):
def write_matrices(pickled_results, genome_ids, output_pairwise, col_ids):
	"""Write above diagonal matrices with allelic differences and shared loci.

	Parameters
	----------
	pickled_results : dict
		Dictionary with sample identifiers as keys
		and paths to binary files with pickled results
		as values.
	genome_ids : list
		List with sample identifiers.
	output_pairwise : str
		Path to the output file to which the matrix
		with pairwise allelic differences will be saved.
	output_p : str
		Path to the output file to which the matrix
		with pairwise shared loci will be saved.
	col_ids: list
		List with sample identifiers to add as headers.

	Returns
	-------
	None.
	"""
	# sl_lines = [col_ids]
	ad_lines = [col_ids]
	limit = 300
	for g in genome_ids:
		current_file = pickled_results[g]
		# load data
		data = fo.pickle_loader(current_file)

		# shared_loci = list(data[0])
		# shared_loci = list(map(str, shared_loci))
		allele_diffs = list(data[0])
		allele_diffs = list(map(str, allele_diffs))

		padding = [''] * (len(genome_ids)-len(allele_diffs))

		# sl_line = [g] + padding + shared_loci
		# sl_lines.append(sl_line)
		ad_line = [g] + padding + allele_diffs
		ad_lines.append(ad_line)

		# if len(sl_lines) >= limit or g == genome_ids[-1]:
		if len(ad_lines) >= limit or g == genome_ids[-1]:
			ad_lines = [im.join_list(line, '\t') for line in ad_lines]
			fo.write_lines(ad_lines, output_pairwise, joiner='\n', write_mode='a')
			ad_lines = []
			# write_lines(sl_lines, output_p, mode='a')
			# sl_lines = []

	return True


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

	return output_transpose


def merge_triangular_matrices(upper_matrix, lower_matrix, output_file, matrix_size):
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
	with open(upper_matrix, 'r') as upper_handle, open(lower_matrix, 'r') as lower_handle:
		upper_reader = csv.reader(upper_handle, delimiter='\t')
		lower_reader = csv.reader(lower_handle, delimiter='\t')

		merged_lines = []
		for i in range(matrix_size):
			upper_line = upper_reader.__next__()
			lower_line = lower_reader.__next__()
			merged_line = [e
						   if e != ''
						   else lower_line[i]
						   for i, e in enumerate(upper_line)]
			merged_lines.append(merged_line)

			if len(merged_lines) >= 200 or i == (matrix_size-1):
				merged_lines = [im.join_list(line, '\t') for line in merged_lines]
				fo.write_lines(merged_lines, output_file, joiner='\n', write_mode='a')
				merged_lines = []


def symmetrify_matrix(input_matrix, matrix_size, tmp_directory):
	"""Symmetrify a triangular matrix.

	Parameters
	----------
	input_matrix : str
		Path to TSV file that contains the triangular matrix.
	matrix_size : int
		Total number of lines in input file.
	tmp_directory : str
		Path to the output temporary directory.

	Returns
	-------
	symmetric_output : str
		Path to the output file that contains the symmetric
		matrix.
	"""
	output_transpose = transpose_matrix(input_matrix, tmp_directory)

	# merge upper and lower diagonal matrices into symmetric matrix
	symmetric_output = input_matrix.replace('.tsv', '_symmetric.tsv')

	merge_triangular_matrices(input_matrix, output_transpose,
							  symmetric_output, matrix_size)

	# delete files with triangular matrices
	os.remove(input_matrix)
	os.remove(output_transpose)

	return symmetric_output


def main(input_matrix, output_directory, cpu_cores, symmetric, masked):
	"""Compute a distance matrix based on allelic profiles.

	Parameters
	----------
	input_matrix : str
		Path to a TSV file with allelic profiles determined by
		the AlleleCall module.
	output_directory : str
		Path to the output directory.
	cpu_cores : int
		Number of CPU cores used to compute distances.
	symmetric : bool
		Determine a symmetric pairwise distance matrix, instead
		of a triangular matrix.
	masked : bool
		False if the input matrix values are masked, True otherwise.
		The process will mask the matrix values it this value is False.

	Returns
	-------
	output_pairwise : str
		Path to the TSV file that contains the distance matrix.
	"""
	# Create output directory if it does not exist
	if os.path.isdir(output_directory) is False:
		os.mkdir(output_directory)

	# determine input basename
	input_basename = os.path.basename(input_matrix)
	# remove extension that is after last '.'
	input_basename = '.'.join(input_basename.split('.')[0:-1])

	output_masked = os.path.join(output_directory, ct.MASKED_PROFILES_BASENAME)
	if masked is False:
		print('Masking profile matrix...', end='')
		profiles_matrix = pd.read_csv(input_matrix,
									  header=0, index_col=0,
									  sep='\t', low_memory=False)
		masked_profiles = profiles_matrix.apply(im.replace_chars)
		masked_profiles.to_csv(output_masked, sep='\t')
		print('masked matrix available at {0}'.format(output_masked))
	else:
		output_masked = input_matrix

	# create temp directory to store pairwise distances per genome
	tmp_directory = os.path.join(output_directory, 'temp', 'pairwise_distances')
	if os.path.isdir(tmp_directory) is False:
		os.mkdir(tmp_directory)

	# get sample identifiers
	genome_ids = get_sample_ids(input_matrix, delimiter='\t')
	total_genomes = len(genome_ids)

	np_matrix = tsv_to_nparray(output_masked)

	rows_indexes = [i for i in range(len(np_matrix))]
	random.shuffle(rows_indexes)
	# divide inputs into 20 lists for 5% progress resolution
	parallel_inputs = im.divide_list_into_n_chunks(rows_indexes, 20)

	common_args = [[l, np_matrix, genome_ids, tmp_directory, compute_distances]
				   for l in parallel_inputs]

	# increasing cpu cores can greatly increase memory usage
	print('Computing pairwise distances...')
	results = mo.map_async_parallelizer(common_args,
										mo.function_helper,
										cpu_cores,
										show_progress=True)

	merged = im.merge_dictionaries(results)

	print('\nCreating distance matrix...', end='')
	# create files with headers
	col_ids = ['FILE'] + genome_ids
	output_pairwise = os.path.join(output_directory, ct.DISTANCE_MATRIX_BASENAME)
	# output_p = os.path.join(output_directory,
	#                         '{0}_shared_loci.tsv'.format(input_basename))

	# Import arrays per genome and save to matrix file
	# results = write_matrices(merged, genome_ids, output_pairwise, output_p, col_ids)
	results = write_matrices(merged, genome_ids, output_pairwise, col_ids)

	if symmetric is True:
		# add 1 to include header
		output_pairwise = symmetrify_matrix(output_pairwise,
											len(genome_ids)+1,
											tmp_directory)
		# output_p = symmetrify_matrix(output_p,
		#                              len(genome_ids)+1,
		#                              tmp_directory)

	print('done.')

	return [output_pairwise]
