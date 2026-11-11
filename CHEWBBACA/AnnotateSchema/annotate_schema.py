#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module retrieves annotations for the loci in a schema. The process 
can retrieve annotations through UniProt's Proteins API endpoint to find 
exact matches. If users provide a taxon/taxa name/s, the process will 
also search for and download reference proteomes for the specified taxon/taxa 
and use BLASTp to align loci representative alleles against the proteomes 
to assign annotations based on the BSR value computed for each alignment.

Code documentation
------------------
"""


import os

from loguru import logger

import CHEWBBACA.global_constants as gct
import CHEWBBACA.AnnotateSchema.constants as ct
from CHEWBBACA.utils import (
	blast_wrapper as bw,
	core_functions as cf,
	file_operations as fo,
	uniprot_requests as ur,
	fasta_operations as fao,
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
		locus_id = fo.file_basename(file).split("_short")[0]
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
				proteome_results.setdefault(locus_id, []).append(hit_terms+[str(res[-1])])

	return proteome_results


def proteome_annotations(schema_directory, loci, temp_directory, taxa, blast_score_ratio,
                         cpu_cores, proteome_matches, blast_path, translation_table):
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
	short_directory = fo.join_paths(schema_directory, ["short"])
	# Chewie schemas
	if fo.exists(short_directory):
		reps_paths = [fo.join_paths(short_directory, [file+"_short.fasta"]) for file in loci[0]]
	# External schemas
	# Get first allele sequence for each locus
	else:
		# Create temporary folder to store FASTA files with the first/representative allele for all loci
		temp_short = fo.join_paths(temp_directory, ["temp_short"])
		fo.create_directory(temp_short)
		# Get the first allele from all FASTA files
		reps_paths = []
		for i, file in enumerate(loci[1]):
			records = fao.sequence_generator(file)
			first_record = records.__next__()
			rep_file = fo.join_paths(temp_short, [loci[0][i]+"_short.fasta"])
			fao.write_records([first_record], rep_file, mode='w')
			reps_paths.append(rep_file)

	logger.info("Translating representative sequences...")
	# Translate representatives for all loci
	translated_reps = fo.join_paths(temp_directory, ["translated_reps"])
	fo.create_directory(translated_reps)

	reps_protein_files = mo.parallelize_function(fao.translate_fasta, reps_paths, [translated_reps, translation_table, False], cpu_cores)
	reps_protein_files = [r[1] for r in reps_protein_files]

	logger.info("Downloading list of reference proteomes...")
	remote_readme = fo.join_paths(ct.UNIPROT_PROTEOMES_FTP, ["README"])
	local_readme = fo.join_paths(temp_directory, ["reference_proteomes_readme.txt"])

	# Get README file with list of reference proteomes
	res = fo.download_file(remote_readme, local_readme)

	# Get lines with proteome info for species of interest
	readme_lines = fo.read_lines(local_readme, strip=False)

	selected_proteomes = im.contained_terms(readme_lines, taxa)
	selected_proteomes = [line.strip("\n") for line in selected_proteomes]
	selected_proteomes = [line.split("\t") for line in selected_proteomes]
	logger.info(f"Found {len(selected_proteomes)} reference proteomes for {taxa}.")
	proteome_results = {}
	if len(selected_proteomes) > 0:
		# Create directory to store proteomes
		proteomes_directory = fo.join_paths(temp_directory, ["proteomes"])
		fo.create_directory(proteomes_directory)

		proteomes_files = ur.get_proteomes(selected_proteomes,
										   proteomes_directory)

		# Uncompress files and concatenate into single FASTA
		uncompressed_proteomes = [fo.unzip_file(file) for file in proteomes_files]
		proteomes_concat = fo.join_paths(proteomes_directory, ["full_proteome.fasta"])
		proteomes_concat = fo.concatenate_files(uncompressed_proteomes, proteomes_concat)

		# Get self-scores
		# Concatenate protein files
		reps_concat = fo.concatenate_files(reps_protein_files, fo.join_paths(temp_directory, ["reps_concat.fasta"]))

		logger.info("Determining self-score of representatives...")
		blastp_path = os.path.join(blast_path, gct.BLASTP_ALIAS)
		makeblastdb_path = os.path.join(blast_path, gct.MAKEBLASTDB_ALIAS)
		blast_version = bw.get_blast_version(blast_path)
		blastdb_aliastool_path = fo.join_paths(blast_path, [gct.BLASTDB_ALIASTOOL_ALIAS])
		self_scores = cf.determine_self_scores(reps_concat, temp_directory,
			makeblastdb_path, blastp_path, "prot", cpu_cores, blastdb_aliastool_path)

		# Create BLASTdb with proteome sequences
		proteome_blastdb = fo.join_paths(proteomes_directory, ["proteomes_db"])
		db_std = bw.make_blast_db(makeblastdb_path, proteomes_concat, proteome_blastdb, "prot")

		# BLASTp to determine annotations
		blast_inputs = [[blastp_path, proteome_blastdb, file, file+"_blastout.tsv",
						 1, 1, None, None, 10, None, bw.run_blast]
						for file in reps_protein_files]

		logger.info("BLASTing representatives against proteomes...")
		blast_results = mo.map_async_parallelizer(blast_inputs, mo.function_helper, cpu_cores)

		blastout_files = [fo.join_paths(translated_reps, [file])
						  for file in os.listdir(translated_reps)
						  if "blastout" in file]

		# index proteome file
		indexed_proteome = fao.index_fasta(proteomes_concat)

		# process results for each BLASTp
		proteome_results = extract_annotations(blastout_files,
											   indexed_proteome,
											   self_scores,
											   blast_score_ratio,
											   proteome_matches)

	return proteome_results


def create_annotations_table(annotations, output_directory, header, schema_name):
	"""Create output table with loci information.

	Parameters
	----------
	annotations : dcit
		Dictionary with loci identifiers as keys and lists containing loci 
		annotation terms as values (each list can contain the information 
		extracted from the TSV file with CDS coordinates, the annotation 
		terms identified based on BLASTp alignment against UniProt's 
		reference proteomes, and the annotation terms retrieved through 
		UniProt's Proteins API).
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
			proteome_data = [";".join(list(map(str, d))) for d in proteome_data]
			proteome_data = ["" if set(d) == {";"} else d for d in proteome_data]
			locus_annotations.extend(proteome_data)

		annotation_lines.append(locus_annotations)

	annotation_lines = ["\t".join(line) for line in annotation_lines]
	output_basename = "{0}_annotations.tsv".format(schema_name)
	output_file = fo.join_paths(output_directory, [output_basename])
	fo.write_lines(annotation_lines, output_file)

	return output_file


def main(config):
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
	no_cleanup
		Do not keep intermediate files.
	blast_path
		Path to the directory that contains the BLAST executables.
	"""
	# Create temp directory
	temp_directory = fo.join_paths(config.output_directory, ["temp"])
	fo.create_directory(temp_directory)

	logger.info(f"Schema directory: {config.schema_directory}")
	logger.info(f"Number of loci: {len(config.loci_list[0])}")
	logger.info(f"BLAST Score Ratio: {config.blast_score_ratio}")
	logger.info(f"Translation table: {config.translation_table}")

	# Search for annotations based on the reference proteomes for the specified taxa
	proteome_results = {}
	if config.proteome_taxa is not None:
		logger.info(f"Searching for annotations based on UniProt's reference proteomes for {config.proteome_taxa}...")
		proteome_results = proteome_annotations(config.schema_directory,
		                                        config.loci_list,
												temp_directory,
												config.proteome_taxa,
												config.blast_score_ratio,
												config.cpu_cores,
												config.proteome_matches,
												config.blast_path,
												config.translation_table)
		logger.info(f"Retrieved annotations for {len(proteome_results)}/{len(config.loci_list[0])} loci.")
	else:
		logger.info("Will not annotate based on UniProt\'s reference proteomes because no taxa names were provided to '--ptax, --proteome-taxa'.")

	rest_results = {}
	if config.rest_api:
		logger.info("Querying UniProt\'s Proteins API...")
		# Get annotations through UniProt Proteins API
		# Create inputs for multiprocessing
		results = mo.parallelize_function(ur.get_annotation, config.loci_list[1], [config.translation_table, config.rest_api_matches], ct.UNIPROT_REST_API_THREADS)
		for r in results:
			if len(r[1]) > 0:
				locus_id = fo.file_basename(r[0], False)
				rest_results[locus_id] = r[1]

		logger.info(f"Retrieved annotations for {len(rest_results)}/{len(config.loci_list[0])} loci.")
	else:
		logger.info("Did not provide the '--rapi, --rest-api' option. Will not query UniProt's Proteins API.")

	if len(proteome_results) == 0 and len(rest_results) == 0:
		# Logged as warning because it is kinda suspicious
		logger.warning("Could not retrieve annotations for any loci in the schema.")
	else:
		# Get column headers
		annotations_header = [ct.ANNOTATE_SCHEMA_FIRST_COLUMN]
		# Read TSV file with CDS coordinates
		loci_coordinates = {}
		if config.protein_table is not None:
			# Use iterator to avoid loading huge files into memory
			table_lines = fo.read_tabular_lazy(config.protein_table)
			# Add column headers from file
			annotations_header += table_lines.__next__()
			for l in table_lines:
				# Create locus ID based on genome ID and protein identifier
				genome_id = l[1].replace("_", "-")
				locus_id = f"{genome_id}-protein{l[5]}"
				if locus_id in config.loci_list[0]:
					loci_coordinates[locus_id] = l

		# Add remaining column headers
		if config.proteome_taxa is not None:
			annotations_header += ct.ANNOTATE_SCHEMA_PROTEOME_COLUMNS

		if config.rest_api:
			annotations_header += ct.ANNOTATE_SCHEMA_RESTAPI_COLUMNS

		# Get annotation terms for all loci
		annotations_lines = {}
		for locus in config.loci_list[0]:
			annotations_lines[locus] = []
			if config.protein_table:
				annotations_lines[locus].extend(loci_coordinates.get(locus, [""]*len(gct.GENE_TABLE_HEADERS[1:])))
			if config.proteome_taxa:
				# Get proteome annotations for locus and combine data from multiple matches
				proteome_data = proteome_results.get(locus, None)
				if proteome_data:
					combined_data = list(zip(*proteome_data))
					combined_data = [";".join(d) if set(d) != {""} else "" for d in combined_data]
					annotations_lines[locus].extend(combined_data)
				else:
					annotations_lines[locus].extend([""]*len(ct.ANNOTATE_SCHEMA_PROTEOME_COLUMNS))
			if config.rest_api:
				# Get REST API annotations for locus and combine data from multiple matches
				rest_api_data = rest_results.get(locus, None)
				if rest_api_data:
					combined_data = list(zip(*rest_api_data))
					combined_data = [";".join(d) if set(d) != {""} else "" for d in combined_data]
					annotations_lines[locus].extend(combined_data)
				else:
					annotations_lines[locus].extend([""]*len(ct.ANNOTATE_SCHEMA_RESTAPI_COLUMNS))

		schema_basename = fo.file_basename(config.schema_directory)
		output_file = create_annotations_table(annotations_lines, config.output_directory, annotations_header, schema_basename)

		logger.info(f"Output file with loci annotations available at {output_file}")

	if config.no_cleanup is False:
		fo.delete_directory(temp_directory)
