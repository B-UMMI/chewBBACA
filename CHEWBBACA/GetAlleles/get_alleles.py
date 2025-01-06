#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module gets the alleles identified in a dataset and saves them to FASTA files.

Code documentation
------------------
"""


import pandas as pd

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   fasta_operations as fao,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 file_operations as fo,
								 fasta_operations as fao,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def get_locus_alleles(locus_id, allele_ids, locus_file, output_directory):
	"""Get alleles from a locus FASTA file.

	Parameters
	----------
	locus_id : str
		The locus identifier.
	allele_ids : list
		The list of allele IDs.
	locus_file : str
		The path to the locus FASTA file.
	output_directory : str
		The path to the output directory.

	Returns
	-------
	output_file : str
		The path to the output FASTA file.
	"""
	# Import sequences
	locus_alleles = fao.import_sequences(locus_file)
	# Only keep allele ID
	locus_alleles = {(k.split('_')[-1]).replace('*', ''): v for k, v in locus_alleles.items()}
	dataset_alleles = [(t, locus_alleles.get(t[1])) for t in allele_ids]
	# Create records
	records = []
	for a in dataset_alleles:
		if a[0][0] is not None:
			header = f'>{a[0][0]}_{locus_id}_{a[0][1]}'
		else:
			header = f'>{locus_id}_{a[0][1]}'
		sequence = a[1]
		record = '\n'.join([header, sequence])
		records.append(record)

	# Save to file
	output_file = fo.join_paths(output_directory, [f'{locus_id}.fasta'])
	fo.write_lines(records, output_file, joiner='\n', write_mode='w')

	return output_file


def main(input_file, schema_directory, output_directory, cpu_cores, distinct, translate, translation_table):
	# Read input file
	print('Reading input file with allelic profiles...')
	profiles = pd.read_csv(input_file, delimiter='\t', dtype=str, index_col=0)
	nsamples, nloci = profiles.shape
	print(f'Total loci: {nloci}')
	print(f'Total samples: {nsamples}')
	# Remove all 'INF-' prefixes, missing data and '*' from identifiers
	# Replace values by '0'
	print('Masking profiles...')
	masked_profiles = profiles.apply(im.replace_chars)

	# Create folder to store FASTA files
	alleles_directory = fo.join_paths(output_directory, ['alleles'])
	fo.create_directory(alleles_directory)

	# Get alleles per locus
	loci = masked_profiles.columns.tolist()
	# Create inputs for multiprocessing
	loci_files = [fo.join_paths(schema_directory, [f'{locus}.fasta']) for locus in loci]
	inputs = []
	loci_stats = {locus: [] for locus in loci}
	for i, locus_id in enumerate(loci):
		# Get alleles identifiers
		locus_column = masked_profiles[locus_id]
		# Remove lines with 0
		locus_column = locus_column[locus_column != '0']
		# Store sample and allele ID
		if not distinct:
			allele_ids = list(zip(locus_column.index, locus_column))
		# Only store allele ID
		else:
			allele_ids = pd.unique(locus_column).tolist()
			# Sort values based
			allele_ids = sorted(allele_ids, key=lambda x: int(x))
			allele_ids = [(None, i) for i in allele_ids]
		inputs.append((locus_id, allele_ids, loci_files[i], alleles_directory, get_locus_alleles))
		# Get stats
		loci_stats[locus_id].extend([fao.count_sequences(loci_files[i]), len(locus_column), len(pd.unique(locus_column))])

	# Save loci statistics
	loci_stats_file = fo.join_paths(output_directory, ['summary_stats.tsv'])
	loci_stats_header = ct.GETALLELES_LOCI_STATS_HEADER
	loci_stats_lines = [[k]+list(map(str, v)) for k, v in loci_stats.items()]
	loci_stats_lines = [loci_stats_header] + ['\t'.join(l) for l in loci_stats_lines]
	fo.write_lines(loci_stats_lines, loci_stats_file)

	# Get and save alleles to FASTA files
	if not distinct:
		print('Creating FASTA files with the alleles identified in the dataset...')
	else:
		print('Creating FASTA files with the distinct alleles identified in the dataset...')
	fasta_files = mo.map_async_parallelizer(inputs, mo.function_helper, cpu_cores)

	if translate:
		# Create folder to store FASTA files
		translated_alleles_directory = fo.join_paths(output_directory, ['translated_alleles'])
		fo.create_directory(translated_alleles_directory)
		print('Translating alleles...')
		inputs = [[file, translated_alleles_directory, translation_table, fao.translate_fasta] for file in fasta_files]
		protein_files = mo.map_async_parallelizer(inputs, mo.function_helper, cpu_cores)

	print(f'Output files available in {output_directory}')
