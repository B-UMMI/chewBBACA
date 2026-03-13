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
	# Only keep the allele ID, do not keep '*' added to new allele IDs in schemas downloaded from Chewie-NS
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
		# Alleles must be in the schema!
		try:
			record = '\n'.join([header, sequence])
		except:
			print(f'Error creating record for {locus_id} with allele {a[0][1]}')
		records.append(record)

	# Save to file
	output_file = fo.join_paths(output_directory, [f'{locus_id}.fasta'])
	fo.write_lines(records, output_file, joiner='\n', write_mode='w')

	return output_file


def main(input_file, schema_directory, genes_list, output_directory, cpu_cores, distinct, translate, translation_table):
	# Read input file
	print('Importing allelic profiles...')
	if not genes_list:
		profiles = pd.read_csv(input_file, delimiter='\t', dtype=str, index_col=0)
	# User provided a list of loci to create FASTA files for
	else:
		# Read list of loci from file
		print(f'User provided list of loci: {genes_list}')
		loci_list = fo.read_lines(genes_list, strip=True)
		profiles = pd.read_csv(input_file, delimiter='\t', dtype=str, index_col=0, usecols=['FILE'] + loci_list)

	nsamples, nloci = profiles.shape
	print(f'Total loci: {nloci}')
	print(f'Total samples: {nsamples}')
	# Remove all 'INF-' prefixes, missing data and '*' from identifiers
	# Replace values by '0'
	print('Masking profiles...')
	masked_profiles = profiles.apply(im.replace_chars)

	# Create folder to store FASTA files
	alleles_directory = fo.join_paths(output_directory, ['fastas'])
	fo.create_directory(alleles_directory)

	# Get alleles per locus
	loci = masked_profiles.columns.tolist()
	# Create inputs for multiprocessing
	loci_files = [fo.join_paths(schema_directory, [f'{locus}.fasta']) for locus in loci]
	inputs = []
	loci_stats = {locus: [] for locus in loci}
	absent_loci = 0
	print('Determining the list of alleles identified for each locus...')
	for i, locus_id in enumerate(loci):
		# Get alleles identifiers
		locus_column = masked_profiles[locus_id]
		# Remove lines with 0
		locus_column = locus_column[locus_column != '0']
		# Store stats
		loci_stats[locus_id].extend([fao.count_sequences(loci_files[i]), len(locus_column), len(pd.unique(locus_column))])
		# Only try to fetch if the locus is in at least one sample
		if len(locus_column) == 0:
			print(f'No alleles for locus {locus_id} in the dataset. Skipping...')
			absent_loci += 1
			continue
		else:
			# Store sample and allele ID
			if not distinct:
				allele_ids = list(zip(locus_column.index, locus_column))
			# Only store allele ID
			else:
				allele_ids = pd.unique(locus_column).tolist()
				# Sort values based
				allele_ids = sorted(allele_ids, key=lambda x: int(x))
				allele_ids = [(None, i) for i in allele_ids]
			# Append input data for multiprocessing
			inputs.append((locus_id, allele_ids, loci_files[i], alleles_directory, get_locus_alleles))

	print(f'Number of loci that were not identified in the dataset: {absent_loci}')

	output_files = []

	# Save loci statistics
	loci_stats_file = fo.join_paths(output_directory, ['summary_stats.tsv'])
	loci_stats_header = ct.GETALLELES_LOCI_STATS_HEADER
	loci_stats_lines = [[k]+list(map(str, v)) for k, v in loci_stats.items()]
	loci_stats_lines = [loci_stats_header] + ['\t'.join(l) for l in loci_stats_lines]
	fo.write_lines(loci_stats_lines, loci_stats_file)
	output_files.append(loci_stats_file)

	# Get and save alleles to FASTA files
	if not distinct:
		print('Creating FASTA files with the alleles identified in the dataset...')
	else:
		print('Creating FASTA files with the distinct alleles identified in the dataset...')
	fasta_files = mo.map_async_parallelizer(inputs, mo.function_helper, cpu_cores)
	fasta_files = [file for file in fasta_files if file is not None]
	output_files.append(fasta_files)

	if translate:
		# Create folder to store FASTA files
		translated_alleles_directory = fo.join_paths(output_directory, ['translated_fastas'])
		fo.create_directory(translated_alleles_directory)
		print('Translating alleles...')
		inputs = [[file, translated_alleles_directory, translation_table, False, fao.translate_fasta] for file in fasta_files]
		protein_files = mo.map_async_parallelizer(inputs, mo.function_helper, cpu_cores)
		# Only keep the translated files
		protein_files = [f[1] for f in protein_files]
		output_files.append(protein_files)

	print(f'Output files available in {output_directory}')
	return output_files
