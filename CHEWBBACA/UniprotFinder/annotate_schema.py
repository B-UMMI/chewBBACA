#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module retrieves annotations for the loci in a schema.
The process can retrieve annotations through UniProt's SPARQL
endpoint to find exact matches. If users provide a taxon/taxa
name/s, the process will also search for reference proteomes
for the specified taxon/taxa and use BLASTp to align local
sequences against reference sequences to assign annotations
based on the BSR value computed for each alignment.

Code documentation
------------------
"""


import os
import sys
from http.client import HTTPResponse

try:
	from utils import (
		constants as ct,
		blast_wrapper as bw,
		core_functions as cf,
		file_operations as fo,
		uniprot_requests as ur,
		fasta_operations as fao,
		parameters_validation as pv,
		iterables_manipulation as im,
		multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (
		constants as ct,
		blast_wrapper as bw,
		core_functions as cf,
		file_operations as fo,
		uniprot_requests as ur,
		fasta_operations as fao,
		parameters_validation as pv,
		iterables_manipulation as im,
		multiprocessing_operations as mo)


def extract_annotations(blastout_files, indexed_proteome, self_scores,
						blast_score_ratio, proteome_matches):
	"""Extract the annotation terms from high-scoring proteome matches.

	Parameters
	----------
	blastout_files : list
		List with the paths to the TSV files withBLASTp results
		in tabular format. One per locus.
	indexed_proteome : Bio.File._IndexedSeqFileDict
		Fasta file index created with BioPython.
	self_scores : dict
		Dictionary with the identifiers of schema representatives
		as keys and the self-alignment raw socre as value.
	blast_score_ratio : float
		BLAST Score Ratio value. Hits with a BSR value
		>= than this value will be considered as high
		scoring hits that can be included in the final
		table according to the maximum number of matches
		to report.
	proteome_matches : int
		Maximum number of proteome matches to report.

	Returns
	-------
	proteome_results : dict
		Dictionary with loci identifiers as keys and a list
		with information about loci retrieved from the most
		similar records in UniProt's reference proteomes.
	"""
	proteome_results = {}
	for file in blastout_files:
		locus_id = fo.file_basename(file).split('_short')[0]
		results = fo.read_tabular(file)
		if len(results) > 0:
			# compute BSR values
			for r in results:
				r.append(round(float(r[6])/float(self_scores[r[0]][1]), 2))

			# sort based on decreasing BSR
			sorted_results = sorted(results, key=lambda x: x[-1], reverse=True)
			# get results equal or above BSR
			high_bsr_results = [r for r in sorted_results
								if r[-1] >= blast_score_ratio]

			for res in high_bsr_results[0:proteome_matches]:
				# get record and extract relevant info
				hit = indexed_proteome[res[4]]
				hit_dict = vars(hit)
				hit_terms = ur.extract_proteome_terms(hit_dict)
				proteome_results.setdefault(locus_id, []).append(hit_terms+[res[-1]])

	return proteome_results


def proteome_annotations(schema_directory, temp_directory, taxa,
						 blast_score_ratio, cpu_cores, proteome_matches,
						 blast_path, translation_table):
	"""Get annotations based on matches against UniProt's reference proteomes.

	Determine loci annotations based on alignment against UniProt's
	reference proteomes.

	Parameters
	----------
	schema_directory : str
		Path to the schema's directory.
	temp_directory : str
		Path to the temporary directory where intermediate
		files will be written to.
	taxa : list
		List of taxa scientific names. The process will
		search for reference proteomes whose "Species Name"
		field contain any of the provided taxa names.
	blast_score_ratio : float
		BLAST Score Ratio value. Hits with a BSR value
		>= than this value will be considered as high
		scoring hits that can be included in the final
		table according to the maximum number of matches
		to report.
	cpu_cores : int
		Number of threads used to run BLASTp.
	proteome_matches : int
		Maximum number of proteome matches to report.
	blast_path : str
		Path to BLAST executables.
	translation_table : int
		Genetic code used to translate DNA sequences.

	Returns
	-------
	proteome_results : dict
		Dictionary with loci identifiers as keys and a list
		with information about loci retrieved from the most
		similar records in UniProt's reference proteomes.
	"""
	# Get paths to files with representative sequences
	short_directory = fo.join_paths(schema_directory, ['short'])
	reps_paths = [fo.join_paths(short_directory, [file])
				  for file in os.listdir(short_directory)
				  if file.endswith('.fasta') is True]

	print('Translating representative sequences...')
	# Translate representatives for all loci
	translated_reps = fo.join_paths(temp_directory, ['translated_reps'])
	fo.create_directory(translated_reps)

	reps_protein_files = mo.parallelize_function(fao.translate_fasta,
												 reps_paths,
												 [translated_reps, translation_table],
												 cpu_cores, True)
	reps_protein_files = [r[1] for r in reps_protein_files]

	print('\nDownloading list of reference proteomes...', end='')
	remote_readme = fo.join_paths(ct.UNIPROT_PROTEOMES_FTP, ['README'])
	local_readme = fo.join_paths(temp_directory,
								 ['reference_proteomes_readme.txt'])

	# get README file with list of reference proteomes
	res = fo.download_file(remote_readme, local_readme)
	print('done.')

	# get lines with proteomes info for species of interest
	readme_lines = fo.read_lines(local_readme, strip=False)

	selected_proteomes = im.contained_terms(readme_lines, taxa)
	selected_proteomes = [line.strip('\n') for line in selected_proteomes]
	selected_proteomes = [line.split('\t') for line in selected_proteomes]
	print('Found {0} reference proteomes for '
		  '{1}.'.format(len(selected_proteomes), taxa))
	proteome_results = {}
	if len(selected_proteomes) > 0:
		# create directory to store proteomes
		proteomes_directory = fo.join_paths(temp_directory, ['proteomes'])
		fo.create_directory(proteomes_directory)

		proteomes_files = ur.get_proteomes(selected_proteomes,
										   proteomes_directory)

		# uncompress files and concatenate into single FASTA
		uncompressed_proteomes = [fo.unzip_file(file) for file in proteomes_files]
		proteomes_concat = fo.join_paths(proteomes_directory,
										 ['full_proteome.fasta'])
		proteomes_concat = fo.concatenate_files(uncompressed_proteomes,
												proteomes_concat)

		# get self-scores
		# concatenate protein files
		reps_concat = fo.concatenate_files(reps_protein_files,
										   fo.join_paths(temp_directory,
														 ['reps_concat.fasta']))

		print('\nDetermining self-score of representatives...', end='')
		blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
		makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)
		blast_version = pv.get_blast_version(blast_path)
		blastdb_aliastool_path = fo.join_paths(blast_path, [ct.BLASTDB_ALIASTOOL_ALIAS])
		self_scores = cf.determine_self_scores(reps_concat, temp_directory,
			makeblastdb_path, blastp_path, 'prot', cpu_cores, blastdb_aliastool_path)
		print('done.')

		# create BLASTdb with proteome sequences
		proteome_blastdb = fo.join_paths(proteomes_directory,
										 ['proteomes_db'])
		db_std = bw.make_blast_db(makeblastdb_path, proteomes_concat, proteome_blastdb, 'prot')

		# BLASTp to determine annotations
		blast_inputs = [[blastp_path, proteome_blastdb, file, file+'_blastout.tsv',
						 1, 1, None, None, 10, None, bw.run_blast]
						for file in reps_protein_files]

		print('\nBLASTing representatives against proteomes...')
		blast_results = mo.map_async_parallelizer(blast_inputs,
												  mo.function_helper,
												  cpu_cores,
												  show_progress=True)

		blastout_files = [fo.join_paths(translated_reps, [file])
						  for file in os.listdir(translated_reps)
						  if 'blastout' in file]

		# index proteome file
		indexed_proteome = fao.index_fasta(proteomes_concat)

		# process results for each BLASTp
		proteome_results = extract_annotations(blastout_files,
											   indexed_proteome,
											   self_scores,
											   blast_score_ratio,
											   proteome_matches)

	return proteome_results


def sparql_annotations(loci_files, translation_table, cpu_cores):
	"""Retrieve annotations from UniProt's SPARQL endpoint.

	Parameters
	----------
	loci_files : list
		List with the paths to the loci FASTA files.

	Returns
	-------
	annotations : list
		List with sublists. Each sublist contains
		the path to the FASTA file of a locus, the
		product name found for that locus and the
		URL to the page of the record that matched the
		locus.
	"""
	# Create inputs to multiprocessing
	uniprot_args = [[gene, translation_table, ur.get_annotation]
					for gene in loci_files]

	# This works with all alleles in the loci to maximize
	# chance of finding annotations
	workers = cpu_cores if cpu_cores <= ct.UNIPROT_SPARQL_THREADS else ct.UNIPROT_SPARQL_THREADS
	annotations = mo.map_async_parallelizer(uniprot_args,
											mo.function_helper,
											workers,
											show_progress=True)

	return annotations


def create_annotations_table(annotations, output_directory, header,
							 schema_name):
	"""Create output table with loci information.

	Parameters
	----------
	annotations : dcit
		Dictionary with loci identifiers as keys and
		lists with information about loci as values (each
		list contains the information extracted from the
		"cds_info.tsv" table, if it was passed to the process,
		and the product and URL link for the match found
		through UniProt's SPARQL endpoint).
	output_directory : str
		Path to the output directory where the table
		will be written to.
	header : list
		File header (first line with column names).
	schema_name : str
		Name of the schema.
	loci_info : bool
		True if the user passed the "cds_info.tsv" table
		to the process, false otherwise.

	Returns
	-------
	output_table : str
		Path to the table with loci information.
	"""
	annotation_lines = [header]
	for locus, data in annotations.items():
		locus_annotations = [locus]
		if isinstance(data[-1], list) is False:
			locus_annotations.extend(data)
		else:
			locus_annotations.extend(data[:-1])
			proteome_data = list(zip(*data[-1]))
			proteome_data = [';'.join(list(map(str, d))) for d in proteome_data]
			proteome_data = ['' if set(d) == {';'} else d for d in proteome_data]
			locus_annotations.extend(proteome_data)

		annotation_lines.append(locus_annotations)

	annotation_lines = ['\t'.join(line) for line in annotation_lines]
	output_basename = '{0}_annotations.tsv'.format(schema_name)
	output_file = fo.join_paths(output_directory, [output_basename])
	fo.write_lines(annotation_lines, output_file)

	return output_file


def main(schema_directory, output_directory, genes_list, protein_table,
		 blast_score_ratio, cpu_cores, taxa, proteome_matches, no_sparql,
		 no_cleanup, blast_path):
	"""Annotate loci in a schema.

	Parameters
	----------
	schema_directory
		Path to the schema directory.
	output_directory
		Path to the output directory where the process will store
		intermediate files and save the results.
	genes_list
		Path to a file that contains a list of schema loci to
		annotate.
	protein_table
		Path to the 'cds_coordinates.tsv' file created by the
		'CreateSchema' process.
	blast_score_ratio
		BLAST Score Ratio value. This value is only used to evaluate
		matches against reference proteomes when a taxon/taxa name/s
		are provided.
	cpu_cores
		Number of CPU cores used by the process.
	taxa
		List of scientific names for a set of taxa. The process will
		download reference proteomes from UniProt and align schema
		translated alleles against the proteomes to find annotations
		for the loci.
	proteome_matches
		Maximum number of proteome matches per locus to report.
	no_sparql
		Do not search for annotations through UniProt's SPARQL
		endpoint.
	no_cleanup
		Do not keep intermediate files.
	blast_path
		Path to the directory that contains the BLAST executables.
	"""
	# Create output directory
	created = fo.create_directory(output_directory)
	if created is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)

	# Create temp directory
	temp_directory = fo.join_paths(output_directory, ['temp'])
	fo.create_directory(temp_directory)

	# Validate input files
	# User provided a list of genes to call
	loci_list = fo.join_paths(temp_directory, [ct.LOCI_LIST])
	if genes_list is not False:
		loci_list = pv.validate_loci_list(genes_list, loci_list, schema_directory)
	# Working with the whole schema
	else:
		loci_list, total_loci = pv.check_input_type(schema_directory, loci_list)

	loci_paths = fo.read_lines(loci_list)
	loci_basenames = [fo.file_basename(locus, False) for locus in loci_paths]

	schema_directory = os.path.dirname(loci_paths[0])
	schema_basename = fo.file_basename(schema_directory)
	config_file = fo.join_paths(schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
	if os.path.isfile(config_file) is True:
		config = fo.pickle_loader(config_file)
		translation_table = config.get('translation_table', [ct.GENETIC_CODES_DEFAULT])[0]
	else:
		translation_table = ct.GENETIC_CODES_DEFAULT

	print(f'Schema: {schema_directory}')
	print(f'Number of loci: {len(loci_paths)}')
	print(f'Translation table: {translation_table}')

	# Find annotations based on reference proteomes for species
	proteome_results = {}
	if taxa is not None:
		proteome_results = proteome_annotations(schema_directory,
												temp_directory,
												taxa,
												blast_score_ratio,
												cpu_cores,
												proteome_matches,
												blast_path,
												translation_table)
	else:
		print('No taxa names provided. Will not annotate based on '
			  'UniProt\'s reference proteomes.')

	sparql_annotated = {}
	sparql_failed = {}
	if no_sparql is False:
		# Check if SPARQL endpoint is up
		available = ur.website_availability(ct.UNIPROT_SPARQL)
		if type(available) != HTTPResponse or available.code != 200:
			print(f'Could not connect to {ct.UNIPROT_SPARQL}')
			print(available)
			print('Cannot retrieve annotations through UniProt\'s SPARQL endpoint.')
		else:
			# Search for annotations through the SPARQL endpoint
			print('\nQuerying UniProt\'s SPARQL endpoint...')

			# Get annotations through UniProt SPARQL endpoint
			results = sparql_annotations(loci_paths, translation_table, cpu_cores)
			for i, r in enumerate(results):
				if fo.file_basename(r[0], False) in loci_basenames:
					if r[1] != '' or r[2] != '':
						sparql_annotated[loci_basenames[i]] = r[1:3]
					else:
						if len(r[-1]) > 0:
							messages = list(set([exc.msg for exc in r[-1]]))
							sparql_failed[loci_basenames[i]] = '\n'.join(messages)
				else:
					sparql_failed[loci_basenames[i]] = r[1]

			print('\nFound annotations for {0}/{1} loci.'.format(len(sparql_annotated), len(loci_paths)))
	else:
		print('\nProvided "--no-sparql" argument. Skipped step to '
			  'search for annotations through UniProt\'s SPARQL '
			  'endpoint.')

	if len(proteome_results) == 0 and len(sparql_annotated) == 0:
		exists = fo.delete_directory(output_directory)
		sys.exit('Could not retrieve annotations for any loci.')
	else:
		header = ['Locus']
		loci_info = {}
		if protein_table is not None:
			# read "cds_coordinates.tsv" file created by CreateSchema
			table_lines = fo.read_tabular(protein_table)
			header += table_lines[0]
			for l in table_lines[1:]:
				# Create locus ID based on genome ID and CDS ID in file
				locus_id = l[0].replace('_', '-')
				locus_id = locus_id + '-protein{0}'.format(l[-2])
				loci_info[locus_id] = l

		if no_sparql is False:
			header += ['Uniprot_Name', 'UniProt_URL']

		if taxa is not None:
			header.extend(['Proteome_ID', 'Proteome_Product',
						   'Proteome_Gene_Name', 'Proteome_Species',
						   'Proteome_BSR'])

		annotations = {}
		for locus in loci_basenames:
			annotations[locus] = []
			if protein_table is not None:
				annotations[locus] += loci_info.get(locus, ['']*6)

			if no_sparql is False:
				annotations[locus] += sparql_annotated.get(locus, ['']*2)

			if taxa is not None:
				annotations[locus].append(proteome_results.get(locus, [['']*5]))

		output_file = create_annotations_table(annotations, output_directory,
											   header, schema_basename)

		print(f'\nOutput file with loci annotations available at {output_file}')

		# Write file with information about cases that failed
		if len(sparql_failed) > 0:
			failed_lines = [f'Locus: {k}\n{str(v)}' for k, v in sparql_failed.items()]
			failed_outfile = fo.join_paths(output_directory, ['failed.txt'])
			failed_text = '\n'.join(failed_lines)
			fo.write_to_file(failed_text, failed_outfile, 'w', '\n')
			print(f'Output file with information about exceptions available at {failed_outfile}')

		if no_cleanup is False:
			exists = fo.delete_directory(temp_directory)
