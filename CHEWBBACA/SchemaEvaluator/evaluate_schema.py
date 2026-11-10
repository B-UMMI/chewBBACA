#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module generates an interactive HTML report for a schema with
statistics about the allele size per locus and an analysis to
identify alleles that are not complete coding sequences or that
deviate from the allele mode size and minimum length. There is
also the option to include loci annotations provided in a TSV
file and to create a detailed HTML report for each locus that
can include a protein MSA and a NJ Tree.

Code documentation
------------------
"""


import os
import math
import json
import shutil
import statistics
from collections import Counter

import CHEWBBACA.global_constants as gct
from CHEWBBACA.ComputeMSA import compute_msa
import CHEWBBACA.SchemaEvaluator.constants as ct
from CHEWBBACA.ComputeMSA.validation_models import ComputeMSAValidator
from CHEWBBACA.utils import (
	file_operations as fo,
	fasta_operations as fao,
	sequence_manipulation as sm,
	iterables_manipulation as im,
	multiprocessing_operations as mo)


def compute_quartiles(values):
	"""Compute the q1, median and q3.

	Parameters
	----------
	values : list
		List with integers and/or floats.

	Returns
	-------
	A list with the q1, median and q3 computed
	based on the input list.
	"""
	median = round(statistics.median(values))
	# q1 and q3
	if len(values) > 1:
		half = int(len(values) // 2)
		q1 = statistics.median(values[:half])
		q3 = statistics.median(values[-half:])
	else:
		q1 = values[0]
		q3 = values[0]

	return [q1, median, q3]


def compute_mean(values, decimal_place=None):
	"""Compute the mean.

	Parameters
	----------
	values : list
		List with integers and/or floats.
	decimal_place : int
		Rounding precision.

	Returns
	-------
	mean : int or float
		The mean computed based on the input list.
	"""
	mean = round(sum(values) / len(values), decimal_place)

	return mean


def compute_sd(values):
	"""Compute the standard deviation.

	Parameters
	----------
	values : list
		List with integers and/or floats.

	Returns
	-------
	sd : float
		The standard deviation computed based on the input list.
	"""
	# standard deviation
	if len(values) > 1:
		sd = statistics.stdev(values)
	else:
		sd = 0.0

	return sd


def outside_threshold(values, threshold, limit):
	"""Determine if values are outside a threshold.

	Parameters
	----------
	values : list
		List with integers and/or floats.
	threshold : int or float
		Threshold value.
	limit : str
		Check if the values are above ('top') or below
		('bot') the threshold.

	Returns
	-------
	outside_indices : list
		List with the indices of the values that were
		outside the threshold.
	"""
	outside_indices = []
	for i, v in enumerate(values):
		if limit == 'top':
			if v > threshold:
				outside_indices.append((i, v))
		elif limit == 'bot':
			if v < threshold:
				outside_indices.append((i, v))

	return outside_indices


def count_translation_exception_categories(exceptions):
	"""Count the number of ocurrences for each translation exception.

	Parameters
	----------
	exceptions : list
		List with the strings corresponding to the exceptions captured
		during sequence translation.

	Returns
	-------
	inframe_stop : int
		Count for in-frame stop codons.
	no_start : int
		Count for sequences that do not start with a start codon.
	no_stop : int
		Count for sequences that do not end with a stop codon.
	incomplete : int
		Count for sequences with a length value that is not a
		multiple of 3.
	ambiguos : int
		Count for sequences that include ambiguous bases.
	"""
	inframe_stop = sum([1 for exc in exceptions
						if ct.TRANSLATION_EXCEPTIONS[0] in exc.split(',')[0]])
	no_start = sum([1 for exc in exceptions
					if ct.TRANSLATION_EXCEPTIONS[1] in exc.split(',')[0]])
	no_stop = sum([1 for exc in exceptions
				   if ct.TRANSLATION_EXCEPTIONS[2] in exc.split(',')[0]])
	incomplete = sum([1 for exc in exceptions
					  if ct.TRANSLATION_EXCEPTIONS[3] in exc.split(',')[0]])
	ambiguous = sum([1 for exc in exceptions
					 if ct.TRANSLATION_EXCEPTIONS[4] in exc.split(',')[0]])

	return [inframe_stop, no_start, no_stop, incomplete, ambiguous]


def reformat_translation_exceptions(exceptions, allele_lengths):
	"""Reformat translation exceptions.

	Parameters
	----------
	exceptions : list
		List with one sublist for each exception. Each sublist
		includes an allele identifier and the exception string.
	allele_lengths : dict
		Dictionary mapping the allele identifier to the allele
		sequence length.

	Returns
	-------
	formatted_exceptions : dict
		Dictionary mapping the allele identifiers to a list
		with the allele identifier, the exception category of
		the first exception that was captured and the list of
		all exceptions captured for each allele.
	"""
	formatted_exceptions = {}
	for exception in exceptions:
		allele_id = exception[0]
		sequence_length = allele_lengths[int(allele_id)]
		raw_exception = exception[1]
		if raw_exception == 'sequence length is not a multiple of 3':
			exception_category = 'Incomplete ORF'
			exception_description = ('Sequence length is not a multiple '
									 f'of 3 ({sequence_length}bp);')
		if raw_exception == 'ambiguous or invalid characters':
			exception_category = 'Ambiguous Bases'
			exception_description = ('Sequence contains ambiguous bases;')
		if 'codon' in raw_exception:
			# split exception string
			split_exception = raw_exception.split(',')
			clean_exceptions = []
			for e in split_exception:
				estr = e.replace('(', ':')
				estr = estr.replace(')', ';')
				estr = estr.replace('.', '')
				clean_exceptions.append(estr)
			if 'Extra' in clean_exceptions[0]:
				exception_category = 'In-frame Stop Codon'
			else:
				exception_category = 'Missing Start/Stop Codon'
			exception_description = ' '.join(clean_exceptions)

		formatted_exceptions.setdefault(allele_id,
										[allele_id,
										 exception_category]).append(exception_description)

	return formatted_exceptions


def get_alleleID(allele_seqid):
	"""
	"""
	if '_*' not in allele_seqid:
		allele_id = int(allele_seqid.split('_')[-1])
	else:
		allele_id = int(allele_seqid.split('_*')[-1])

	return allele_id


def compute_locus_statistics(locus, translation_table, minimum_length,
							 size_threshold, translation_dir, distinct_dir):
	"""Compute sequence length statistics for a locus.

	Parameters
	----------
	locus : str
		Path to the locus FASTA file.
	translation_table : int
		Translation table to use for sequence translation.
	minimum_length : int
		Minimum sequence length value used to exclude sequences.
	size_threshold : float
		Sequence size variation threshold value used to compute bot
		and top sequence length thresholds.
	translation_dir : str
		Path to the directory where the FASTA files with the
		translated alleles will be saved to.

	Returns
	-------
	results : list
		A list with the length statistics and the sequence translation
		exceptions captured for each sequence in the FASTA file.
	"""
	locus_id = fo.file_basename(locus, False)
	sequence_lengths = fao.sequence_lengths(locus)
	ns_alleles = [get_alleleID(k) for k in sequence_lengths if '*' in k]
	ns_alleles = sorted(ns_alleles)
	# sort based on sequence length
	# Determine if any allele identifiers include "*"
	# "*" is added to novel alleles in schemas downloaded from Chewie-NS
	allele_lengths = {}
	for x in sorted(sequence_lengths.items(), key=lambda item: item[1]):
		allele_id = get_alleleID(x[0])
		allele_lengths[allele_id] = x[1]

	lengths = list(allele_lengths.values())
	allele_ids = list(allele_lengths.keys())

	# Determine missing allele ids
	missing_ids = im.find_missing(allele_ids)

	# Get total number of alleles
	nr_alleles = len(lengths)

	# Summary statistics
	max_length = max(lengths)
	min_length = min(lengths)
	size_range = f'{min_length}-{max_length}'
	mean_length = compute_mean(lengths, 0)
	locus_sd = compute_sd(lengths)
	q1, median, q3 = compute_quartiles(lengths)
	mode_length = sm.determine_mode(lengths)[0]

	# Get index of alleles above size threshold
	top_threshold = math.floor(mode_length*(1+size_threshold))
	above_threshold = outside_threshold(lengths, top_threshold, 'top')
	# Get ids based on indices
	above_threshold = [(allele_ids[i[0]], i[1]) for i in above_threshold]
	# Get index of alleles below size threshold
	bot_threshold = math.ceil(mode_length*(1-size_threshold))
	below_threshold = outside_threshold(lengths, bot_threshold, 'bot')
	# Get ids based on indices
	below_threshold = [(allele_ids[i[0]], i[1]) for i in below_threshold]

	# Translate alleles and capture translation exceptions
	_, protein_file, _, exceptions = fao.translate_fasta(locus,
														 translation_dir,
														 translation_table,
														 False)
	# If some sequence headers include "*", reformat FASTA file to remove "*"
	# Reformatting protein files uses less disk space than reformatting DNA
	# Do not attempt to reformat if no sequences could be translated and protei_file is a NoneType
	if len(ns_alleles) > 0 and len(protein_file) < nr_alleles:
		records = fao.sequence_generator(protein_file)
		records = [[(rec.id).replace('*', ''), str(rec.seq)] for rec in records]
		output_lines = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, records)
		fo.remove_files([protein_file])
		fo.write_lines(output_lines, protein_file)

	# Determine distinct proteins
	if protein_file is not None:
		translated_alleles = fao.import_sequences(protein_file)
		distinct_proteins = sm.determine_duplicated_seqs(translated_alleles)
		distinct_records = [[v[0], k]
							for k, v in distinct_proteins.items()]
		distinct_lines = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, distinct_records)
		distinct_file = fo.join_paths(distinct_dir, [f"{locus_id}_distinct.fasta"])
		fo.write_lines(distinct_lines, distinct_file)
		distinct_ids = [[get_alleleID(i) for i in v]
						for k, v in distinct_proteins.items()]
		distinct_ids = [[v[0], v] for v in distinct_ids]
		# Sort in order of decreasing length
		distinct_ids = sorted(distinct_ids, key=lambda x: len(x[1]), reverse=True)
		distinct_ids = sorted(distinct_ids, key=lambda x: x[0])
	else:
		distinct_file = None
		distinct_ids = []

	exceptions = {str(get_alleleID(exc[0])): exc[1] for exc in exceptions}
	exceptions_values = list(exceptions.values())
	exceptions_lines = [[k, v] for k, v in exceptions.items()]
	# Count number of ocurrences for each translation exception
	exception_counts = count_translation_exception_categories(exceptions_values)

	# Determine list of valid allele ids
	valid_ids = [str(i) for i in allele_ids if str(i) not in exceptions]

	# Do not count shorter alleles as invalid alleles
	# Only count as invalid alleles that cannot be translated
	invalidCDS = sum(exception_counts)
	validCDS = nr_alleles - invalidCDS
	validated_proportion = round(validCDS / nr_alleles, 3)

	# Determine sequences shorter than minimum length
	short_ids = [(allele_ids[i], v)
				 for i, v in enumerate(lengths)
				 if v < minimum_length]

	# Format exception strings for Exceptions Table Component in loci reports
	formatted_exceptions = reformat_translation_exceptions(exceptions_lines,
														   allele_lengths)

	# Add info about short alleles
	short_category = f"Alleles < {minimum_length}bp"
	short_description = "Sequence shorter than {0}bp ({1}bp);"
	for i in short_ids:
		if str(i[0]) in formatted_exceptions:
			formatted_exceptions[str(i[0])].append(short_description.format(minimum_length, i[1]))
		else:
			formatted_exceptions[str(i[0])] = [str(i[0]),
											   short_category,
											   short_description.format(minimum_length, i[1])]

	# Add info about alleles above threshold
	above_category = "Alleles above size threshold"
	above_description = "Sequence above size threshold ({0}bp)"
	for i in above_threshold:
		if str(i[0]) in formatted_exceptions:
			formatted_exceptions[str(i[0])].append(above_description.format(i[1]))
		else:
			formatted_exceptions[str(i[0])] = [str(i[0]),
											   above_category,
											   above_description.format(i[1])]

	# Add info about alleles below threshold
	below_category = "Alleles below size threshold"
	below_description = "Sequence below size threshold ({0}bp)"
	for i in below_threshold:
		if str(i[0]) in formatted_exceptions:
			formatted_exceptions[str(i[0])].append(below_description.format(i[1]))
		else:
			formatted_exceptions[str(i[0])] = [str(i[0]),
											   below_category,
											   below_description.format(i[1])]

	# Join exceptions per allele
	formatted_exceptions = [[v[0], v[1], ' '.join(v[2:])]
							for k, v in formatted_exceptions.items()]

	# Unpack exception counts
	inframe_stop, no_start, no_stop, incomplete, ambiguous = exception_counts

	results = [locus_id, nr_alleles, max_length, min_length,
			   size_range, median, mean_length, mode_length,
			   locus_sd, q1, q3, lengths, allele_ids,
			   [locus_id, nr_alleles, validCDS, invalidCDS,
				validated_proportion, len(distinct_ids), incomplete, ambiguous,
				no_start+no_stop, inframe_stop, len(short_ids),
				len(below_threshold), len(above_threshold), len(missing_ids)],
			   bot_threshold, top_threshold, formatted_exceptions,
			   protein_file, distinct_file, missing_ids, valid_ids,
			   ns_alleles, distinct_ids]

	return results


def locus_report(locus_file, locus_data, annotation_columns,
				 annotation_values, msa_file, html_dir, minimum_length,
				 light, add_sequences):
	"""Create the HTML report for a locus.

	Parameters
	----------
	locus_file : str
		Path to the locus FASTA file.
	locus_data : list
		List with the locus data returned by `compute_locus_statistics`.
	annotation_columns : list
		List with the column headers for the annotations table.
	annotation_values : list
		List with the locus annotations.
	html_dir : str
		Path to the directory where the loci HTML reports will
		be saved to.
	minimum_length : int
		Minimum sequence length.
	light : bool
		True to compute and add the MSA and NJ to the report.
		False otherwise.
	add_sequences : bool
		True to add the allele DNA sequences and the translated
		allele sequences to the report.

	Returns
	-------
	locus_html_file : str
		Path to the locus report HTML file.
	"""
	locus = locus_data[0]
	allele_lengths = locus_data[11]
	allele_ids = locus_data[12]

	# determine counts per distinct allele size
	counts = list(Counter(allele_lengths).items())
	sorted_counts = sorted(counts, key=lambda x: x[0])
	counts_data = list(zip(*sorted_counts))

	locus_rows = [locus,
				  locus_data[1],
				  locus_data[13][2],
				  locus_data[13][3],
				  locus_data[13][4],
				  locus_data[13][5],
				  locus_data[13][6],
				  locus_data[13][7],
				  locus_data[13][8],
				  locus_data[13][9],
				  locus_data[13][10],
				  locus_data[4],
				  locus_data[5],
				  locus_data[7],
				  locus_data[13][11],
				  locus_data[13][12],
				  locus_data[13][13]]

	dna_sequences = {"sequences": []}
	protein_sequences = {"sequences": []}
	# Include DNA and Protein sequences if --add-sequences was provided
	if add_sequences is True:
		if locus_data[17] is not None:
			protein_records = fao.sequence_generator(locus_data[17])
			for record in protein_records:
				allele_id = get_alleleID(record.id)
				protein_sequences["sequences"].append({"name": allele_id,
													"sequence": str(record.seq)})
		dna_records = fao.sequence_generator(locus_file)
		for record in dna_records:
			allele_id = get_alleleID(record.id)
			dna_sequences["sequences"].append({"name": allele_id,
											   "sequence": str(record.seq)})

	# Get data for MSA and NJ Tree if --light flag was not provided
	phylo_data = {"phylo_data": []}
	msa_data = {"sequences": []}
	if light is False:
		if locus_data[13][2] > 1 and locus_data[13][5] > 1:
			# Get MSA data
			alignment_text = fo.read_file(msa_file)
			# Remove locus name from allele identifiers in MSA
			alignment_text = alignment_text.replace(f'>{locus}_', '>')
			msa_data['sequences'] = alignment_text

			# Get Tree data
			# Get the phylocanvas data
			tree_file = msa_file+'.tree'
			phylo_text = fo.read_file(tree_file)
			# Remove locus name from allele identifiers in Newick tree
			phylo_text = phylo_text.replace(f'{locus}_', '')
			# Remove node numbers from Newick tree
			# Start by substituting greatest value to avoid substituting
			# smaller values contained in greater values
			for i in range(locus_data[1], 0, -1):
				phylo_text = phylo_text.replace(f'{i}_', '')

			phylo_text = phylo_text.replace('\n', '')
			phylo_data['phylo_data'] = phylo_text

	locus_columns = ct.LOCI_REPORTS_SUMMARY_TABLE_HEADER.format(minimum_length,
											locus_data[14],
											locus_data[15])
	locus_columns = locus_columns.split('\t')

	# Gather data in dictionary to store in HTML
	locus_html_data = {"summaryData": [{"columns": locus_columns},
									   {"rows": [locus_rows]}],
					   "annotations": [{"columns": annotation_columns},
									   {"rows": [annotation_values]}],
					   "lengths": allele_lengths,
					   "ids": allele_ids,
					   "counts": [list(counts_data[0]), list(counts_data[1])],
					   "phylo": phylo_data,
					   "validIDs": locus_data[20],
					   "msa": msa_data,
					   "dna": dna_sequences,
					   "protein": protein_sequences,
					   "botThreshold": locus_data[14],
					   "topThreshold": locus_data[15],
					   "invalidAlleles": [{"columns": ct.INVALID_ALLELES_TABLE_HEADER},
										  {"rows": locus_data[16]}],
					   "nsAlleles": locus_data[21],
					   "distinctAlleles": [{"columns": ct.DISTINCT_ALLELES_TABLE_HEADER},
										   {"rows": locus_data[22]}],
					   }

	# Add data to HTML string and save HTML string into file
	locus_html = ct.LOCUS_REPORT_HTML
	locus_html = locus_html.format(json.dumps(locus_html_data))

	locus_html_file = fo.join_paths(html_dir, [f'{locus}.html'])
	fo.write_to_file(locus_html, locus_html_file, 'w', '\n')

	return locus_html_file


def main(config):
	"""Evaluate a schema and create a HTML report.

	Parameters
	----------
	schema_directory : str
		Path to the schema directory.
	output_directory : str
		Path to the output directory where the report files will be created.
	genes_list : str
		Path to a file with the list of loci to evaluate.
	annotations : str
		Path to a TSV file created by the UniprotFinder module.
	translation_table : int
		Genetic code used to translate alleles.
	size_threshold : float
		Allele size variation threshold. Used to determine if the size
		of an allele is within the interval of the locus size mode +/-
		the size threshold.
	minimum_length : int
		Minimum sequence length accepted for an allele to be included
		in the schema.
	cpu_cores : int
		Number of CPU cores used to run the process.
	loci_reports : bool
		Create individual reports for the loci.
	light : bool
		Skips MSA and NJ tree computations for loci reports.
	add_sequences : bool
		Adds Code Editor components with the allele DNA and
		protein sequences.
	"""
	# Create directory to store intermediate files
	temp_directory = fo.join_paths(config.output_directory, ["temp"])
	fo.create_directory(temp_directory)

	# Sort list of schema loci
	schema_files = sorted(config.loci_list[1])
	# Map basename to full path
	loci_basenames = {fo.file_basename(file, False): file for file in config.loci_list[1]}

	# Check if the schema was created with chewBBACA
	if config.chewie_schema is not False:
		print("The schema was created with chewBBACA v{0}.".format(config.chewie_schema))
		# Message displayed in the Schema Creation Alert
		creation_message = (f"Created with chewBBACA v{config.chewie_schema}\n"
							f"BLAST Score Ratio: {', '.join(map(str, config.blast_score_ratio_values))}\n"
							f"Minimum length: {', '.join(map(str, config.minimum_length_values))}\n"
							f"Size threshold: {', '.join(map(str, config.size_threshold_values))}\n"
							f"Translation table: {', '.join(map(str, config.translation_table_values))}\n"
							f"Training file: {config.training_file}")
	else:
		# Message displayed in the Schema Creation Alert
		creation_message = ("Could not find a config file. Evaluated as external schema.")

	# Message displayed in the Schema Evaluation Alert
	evaluation_message = (f"Schema evaluated with minimum length of {config.minimum_length}, "
	                      f"size threshold of {config.size_threshold}, and translation table "
						  f"{config.translation_table}.")

	# Create temporary directory to store translated alleles and distinct proteins
	translation_dir = fo.join_paths(temp_directory, ["translated_loci"])
	distinct_dir = fo.join_paths(translation_dir, ["distinct_proteins"])
	fo.create_directory(distinct_dir)

	# Calculate the summary statistics and other information about each locus.
	inputs = im.divide_list_into_n_chunks(schema_files, len(schema_files))

	common_args = [config.translation_table, config.minimum_length, config.size_threshold, translation_dir, distinct_dir]

	# Add common arguments to all sublists
	inputs = im.multiprocessing_inputs(inputs, common_args, compute_locus_statistics)

	# Compute loci statistics
	print("Computing loci statistics...")
	results = mo.map_async_parallelizer(inputs, mo.function_helper, config.cpu_cores, show_progress=True)

	# Group values for the same statistic
	data = list(zip(*results))

	# Read loci annotations from TSV file
	annotation_values = [[], []]
	if config.annotations_file is not None:
		annotation_lines = fo.read_tabular(config.annotations_file)
		# Add file columns
		annotation_values[0].extend(annotation_lines[0])
		# Add sublists with lines
		# Only keep lines for loci in the schema
		loci_annotations = [line for line in annotation_lines[1:] if line[0] in loci_basenames]
		annotation_values[1].extend(loci_annotations)

	print(f"\nProvided annotations for {len(annotation_values[1])} loci in the schema.")

	# Columns for the Summary Data Table in the Schema Report
	summary_columns = ct.SCHEMA_SUMMARY_TABLE_HEADER.format(config.minimum_length)
	summary_columns = summary_columns.split('\t')

	# Get total number of alleles with a length value not multiple of 3
	notMultiple_sum = sum([subdata[6] for subdata in data[13]])
	# Get total number of alleles with an in-frame stop codon
	stopC_sum = sum([subdata[9] for subdata in data[13]])
	# Get total number of alleles with no start or stop codon
	notStart_sum = sum([subdata[8] for subdata in data[13]])
	# Get total number of alleles shorter than the minimum length value
	shorter_sum = sum([subdata[10] for subdata in data[13]])
	# Get totla number of alleles with ambiguous bases
	ambiguous_sum = sum([subdata[7] for subdata in data[13]])
	# Get total number of alleles below or above the sequence length thresholds
	below_sum = sum([subdata[11] for subdata in data[13]])
	above_sum = sum([subdata[12] for subdata in data[13]])
	# Get total number of valid and invalid alleles
	invalid_sum = sum([subdata[3] for subdata in data[13]])
	valid_sum = sum(data[1]) - invalid_sum
	# Create list with row values for Summary Data Table in the Schema Report
	summary_rows = [len(data[0]),
					sum(data[1]),
					valid_sum,
					invalid_sum,
					notMultiple_sum,
					ambiguous_sum,
					notStart_sum,
					stopC_sum,
					shorter_sum,
					below_sum,
					above_sum]

	# Columns for the Allele Analysis Table component
	analysis_columns = ct.LOCI_ANALYSIS_TABLE_HEADER.format(config.minimum_length)
	analysis_columns = analysis_columns.split("\t")
	analysis_rows = list(data[13])

	# Message displayed if schema was downloaded from Chewie-NS
	ns_schema = []
	if config.ns_config:
		total_ns = sum([len(subdata) for subdata in data[21]])
		ns_schema.append(f"This schema was downloaded from Chewie-NS and contains {total_ns} "
		                 "alleles added since the last synchronization with the remote schema.")

	# Data in the Schema Report HTML
	schema_data = {"summaryData": [{"columns": summary_columns},
								   {"rows": [summary_rows]}],
				   "loci": list(data[0]),
				   "total_alleles": list(data[1]),
				   "max": list(data[2]),
				   "min": list(data[3]),
				   "median": list(data[5]),
				   "mode": list(data[7]),
				   "q1": list(data[9]),
				   "q3": list(data[10]),
				   "annotations": [{"columns": annotation_values[0]},
								   {"rows": annotation_values[1]}],
				   "analysis": [{"columns": analysis_columns},
								{"rows": analysis_rows}],
				   "lociReports": 1 if config.loci_reports else 0,
				   "evaluationConfig": evaluation_message,
				   "creationConfig": creation_message,
				   "nsSchema": ns_schema}

	# Write Schema Report HTML file
	schema_html = ct.SCHEMA_REPORT_HTML
	schema_html = schema_html.format(json.dumps(schema_data))
	schema_html_file = fo.join_paths(config.output_directory, [ct.SCHEMA_REPORT_BASENAME])
	fo.write_to_file(schema_html, schema_html_file, 'w', '\n')

	# Compute data for loci reports
	if config.loci_reports is True:
		# Create mapping between locus ID and annotations
		if config.annotations_file is not None:
			annotations_dict = {a[0]: a for a in annotation_values[1]}
		else:
			annotations_dict = {}

		# Create directory to store loci HTML reports
		html_dir = fo.join_paths(config.output_directory, ['loci_reports'])
		fo.create_directory(html_dir)

		# Previous data includes values per locus
		loci_data = results

		# Call the ComputeMSA module to compute the loci MSAs
		protein_files_dir = distinct_dir
		# Check if there are any protein files to compute MSAs for before calling the ComputeMSA module
		if len(fo.listdir(protein_files_dir)) > 0:
			# Need to pass complete MAFFT default parameters defined in constants.py to create tree output
			print("Calling the ComputeMSA module to compute the loci MSAs...")
			# Create config object to pass to ComputeMSA
			compute_msa_outdir = fo.join_paths(temp_directory, ["MSAs"])
			compute_msa_config = ComputeMSAValidator(output_directory=compute_msa_outdir,
			                                         input_path=protein_files_dir,
													 schema_directory=None,
													 dna_msa=False,
													 output_variable=False,
													 translation_table=config.translation_table,
													 cpu_cores=config.cpu_cores,
													 only_loci_msas=False,
													 gaps="exclude",
													 ambiguous="exclude",
													 custom_mafft_options=' '.join(gct.MAFFT_DEFAULT_PARAMETERS),
													 protein_input=True,
													 no_cleanup=True)

			compute_msa.main(compute_msa_config)

			msa_files_dir = fo.join_paths(compute_msa_outdir, ["MSAs", "protein"])
			msa_files = fo.listdir(msa_files_dir, full_path=True)
			msa_files_dict = {fo.file_basename(f, False).replace("_distinct", ''): f for f in msa_files}

			# Copy tree files to diretory with MSA files
			tree_files_dir = fo.join_paths(temp_directory, ["MSAs", "input_files"])
			tree_files = fo.listdir(tree_files_dir, full_path=True)
			tree_files = fo.filter_by_extension(tree_files, [".tree"])[0]
			for tf in tree_files:
				shutil.copy(tf, msa_files_dir)

			# Delete temporary directory with input files
			fo.delete_directory(tree_files_dir)
		else:
			msa_files_dict = {}

		inputs = [[loci_basenames[d[0]],
			 	   d,
				   annotation_values[0],
				   annotations_dict.get(d[0], []),
				   msa_files_dict.get(d[0], None)] for d in loci_data]

		common_args = [html_dir, config.minimum_length, config.light, config.add_sequences]

		# Add common arguments to all sublists
		inputs = im.multiprocessing_inputs(inputs, common_args, locus_report)

		# Create loci reports
		# Compute MSA and NJ tree if --light is not True
		# Add DNA and Protein sequences if --add-sequences is provided
		print("Creating loci reports...")
		loci_htmls = mo.map_async_parallelizer(inputs,
											   mo.function_helper,
											   config.cpu_cores,
											   show_progress=True)

	# Copy the JS bundle files to the respective directories
	# JS bundle used by schema report
	# When chewBBACA is installed
	script_path = os.path.dirname(os.path.abspath(__file__))
	# For development
	if script_path.endswith("CHEWBBACA") is False:
		script_path = os.path.dirname(script_path)

	shutil.copy(fo.join_paths(script_path, [ct.SCHEMA_EVALUATOR_SCHEMA_BUNDLE]), config.output_directory)

	if config.loci_reports is True:
		shutil.copy(fo.join_paths(script_path, [ct.SCHEMA_EVALUATOR_LOCI_BUNDLE]), html_dir)

	# Delete all temporary files
	fo.delete_directory(temp_directory)

	print(f'\nResults available in {os.path.abspath(config.output_directory)}.')
