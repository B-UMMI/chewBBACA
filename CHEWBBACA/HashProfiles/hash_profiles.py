#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module hashes allelic profiles.

Code documentation
------------------
"""


import sys
import zlib
import hashlib

import pandas as pd

try:
	from utils import (file_operations as fo,
					   fasta_operations as fao,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (file_operations as fo,
								 fasta_operations as fao,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def hash_column(column, locus_file, hashing_function):
	"""Substitute allele identifiers by allele sequence hashes.

	Parameters
	----------
	column : pandas.core.series.Series
		Column with allele identifiers to substitute by the
		hash values computed from each allele sequence.
	locus_file : str
		Path to the FASTA file that contains the locus alleles.
	hashing_function : func
		Hashing function used to hash each allele.

	Returns
	-------
	hashed_column : pandas.core.series.Series
		Column where each allele identifier was substituted by
		the hash computed from the allele sequence.
	"""
	# Read locus Fasta file
	# Replace '*' added to allele identifiers if schema was downloaded from Chewie-NS
	locus_alleles = {(rec.id).split('_')[-1].replace('*', ''): str(rec.seq)
					 for rec in fao.sequence_generator(locus_file)}

	hashed_alleles = {}
	for seqid, seq in locus_alleles.items():
		# Hash function does not accept string object, encode to get bytes object
		hashed_seq = hashing_function(seq.encode())
		# Bitwise operation to convert crc32 and adler32 hashes to unsigned
		# integer and ensure the computed value is the same for Python 2 & 3
		if isinstance(hashed_seq, int):
			hashed_seq &= 0xffffffff
		else:
			hashed_seq = hashed_seq.hexdigest()
		hashed_alleles[seqid] = hashed_seq

	# Faster than replace or map with update to avoid adding NaN
	hashed_column = column.apply(lambda x: hashed_alleles.get(x, x))

	return hashed_column


def hash_profiles(input_file, loci_ids, loci_files, hashing_function,
				  nrows, skiprows, output_directory):
	"""Hash a set of allelic profiles read from a TSV file.

	Parameters
	----------
	input_file : str
		Path to the TSV file that contains the allelic profiles.
	loci_ids : list
		List with the loci identifiers.
	loci_files : dict
		Dictionary with locus identifiers as keys and path to loci
		FASTA files as values.
	hashing_function : func
		Hashing function used to hash each allele.
	nrows : int
		Number of rows/allelic profiles to read from the input
		file.
	skiprows : range
		Range of rows to skip.
	output_directory : str
		Path to the output directory.

	Returns
	-------
	output_file : str
		Path to the output file with the hashed profiles.
	"""
	current_rows = pd.read_csv(input_file, delimiter='\t', dtype=str,
							   skiprows=skiprows, nrows=nrows, index_col=0)

	# Remove all 'INF-' prefixes, missing data and '*' from identifiers
	# Replace values by '-'
	current_rows = current_rows.apply(im.replace_chars, args=('-'))

	hashed_profiles = []
	for locus in loci_ids:
		locus_column = current_rows[locus]
		hashed_column = hash_column(locus_column, loci_files[locus],
									hashing_function)
		hashed_profiles.append(hashed_column)

	hashed_df = pd.concat(hashed_profiles, axis=1)
	start = skiprows.stop
	stop = skiprows.stop+len(hashed_df)-1
	input_basename = fo.file_basename(input_file, False)
	output_file = fo.join_paths(output_directory,
								['{0}_{1}-{2}_hashed.tsv'.format(input_basename, start, stop)])
	hashed_df.to_csv(output_file, sep='\t', index=True, header=False)

	return output_file


def main(input_file, schema_directory, output_directory, hash_type, cpu_cores, nrows):
	"""Hash allele identifiers in a matrix of allelic profiles.

	Parameters
	----------
	input_file : str
		Path to the TSV file with allelic profiles.
	schema_directory : str
		Path to the schema directory.
	output_directory : str
		Path to the output directory.
	hash_type : str
		Hashing algorithm to use. The module supports
		all hashing algorithms implemented in the hashlib
		and zlib Python libraries.
	cpu_cores : int
		Number of CPU cores used by the process.
	nrows : int
		Divide input file into subsets to process more efficiently.
	"""
	# Get hash function
	hashing_function = getattr(hashlib, hash_type, None)
	if hashing_function is None:
		hashing_function = getattr(zlib, hash_type, None)

	if hashing_function is None:
		sys.exit('{0} hash function is not available in the hashlib '
				 '(https://docs.python.org/3/library/hashlib.html) and '
				 'zlib (https://docs.python.org/3/library/zlib.html) modules.'.format(hash_type))
	else:
		print(f'Hash function: {hash_type}')

	# Get loci identifiers
	with open(input_file, 'r') as infile:
		header = infile.readline()
		loci_ids = header.split()[1:]
	print(f'Total loci: {len(loci_ids)}')

	# Get input/sample identifiers
	sample_ids = pd.read_csv(input_file, delimiter='\t',
							 dtype=str, usecols=['FILE'])
	sample_ids = sample_ids['FILE'].tolist()
	print(f'Total samples: {len(sample_ids)}')

	print(f'Dividing input file into subsets of {nrows} rows to process more efficiently...')
	loci_files = {}
	for locus in loci_ids:
		locus_file = fo.join_paths(schema_directory, [locus])
		# Add .fasta extension if file headers did not include it
		if locus_file.endswith('.fasta') is False:
			locus_file += '.fasta'
		loci_files[locus] = locus_file

	# Create multiprocessing inputs
	multi_inputs = []
	# Divide and process by row chunks
	for i in range(0, len(sample_ids), nrows):
		multi_inputs.append([input_file, loci_ids, loci_files,
							 hashing_function, nrows, range(1, i+1),
							 output_directory, hash_profiles])

	print('Hashing profiles...')
	hashed_files = mo.map_async_parallelizer(multi_inputs, mo.function_helper, cpu_cores)

	print('Creating output file...')
	# Create file with header
	header_basename = fo.file_basename(input_file, False) + '_header.tsv'
	header_file = fo.join_paths(output_directory, [header_basename])
	fo.write_to_file(header, header_file, 'w', '')
	# Concatenate all files
	output_basename = fo.file_basename(input_file) + '_hashed.tsv'
	output_file = fo.join_paths(output_directory, [output_basename])
	fo.concatenate_files([header_file]+hashed_files, output_file)

	# Delete intermediate dataframes
	print('Removing intermediate files...')
	fo.remove_files([header_file]+hashed_files)
	print(f'Hashed profiles saved to {output_file}')
