#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to gene prediction with Pyrodigal.

Code documentation
------------------
"""


import math

import pyrodigal

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


def create_gene_finder(training_data, closed, mask, meta):
	"""Create a Pyrodigal GeneFinder object.

	Parameters
	----------
	training_data : pyrodigal.TrainingInfo
		A training info instance used to predict genes in single
		mode.
	closed : bool
		True to prevent prediction of partial genes at edges of
		sequences, False otherwise.
	meta: bool
		True to run Prodigal in `meta` mode (uses pre-trained
		profiles).

	Returns
	-------
	gene_finder : pyrodigal.GeneFinder
		A GeneFinder object configured based on provided arguments.
	"""
	gene_finder = pyrodigal.GeneFinder(training_info=training_data,
									   closed=closed,
									   mask=mask,
									   meta=meta)

	return gene_finder


def train_gene_finder(gene_finder, sequences, translation_table):
	"""Train a Pyrodigal GeneFinder object based on a set of sequences.

	Parameters
	----------
	gene_finder : pyrodigal.GeneFinder
		A GeneFinder object.
	sequences : list
		Sequences used to train the GeneFinder (list
		of bytes objects).
	translation_table : int
		Translation table to use.

	Return
	------
	gene_finder : pyrodigal.GeneFinder
		A GeneFinder object configured based on provided arguments.
	"""
	training_info = gene_finder.train(*sequences, translation_table=translation_table)

	return training_info


def create_training_file(input_file, output_directory, translation_table):
	"""
	"""
	records = fao.sequence_generator(input_file)
	records = {rec.id: bytes(rec.seq) for rec in records}
	gene_finder = create_gene_finder(None, True, True, False)
	training_info = train_gene_finder(gene_finder, records.values(), translation_table)
	training_file = fo.join_paths(output_directory, [fo.file_basename(input_file, False)+'.trn'])
	fo.pickle_dumper(training_info, training_file)

	return training_file


def read_training_file(training_file):
	"""Load training info for Pyrodigal from Prodigal training file.

	Parameters
	----------
	training_file : str
		Path to Prodigal training file.

	Returns
	-------
	training_data : pyrodigal.TrainingInfo
		The deserialized training info.
	"""
	######### Need to improve training file type detection
	# Pyrodigal training file must be read like this if created with pickle
	try:
		training_data = fo.pickle_loader(training_file)
	except Exception as e:
		# Prodigal training file must be read like this
		with open(training_file, 'rb') as infile:
			training_data = pyrodigal.TrainingInfo.load(infile)

	return training_data


def get_gene_info(contig_id, genome_id, protid, genes):
	"""Get genes information from a pyrodigal.Genes object.

	Parameters
	----------
	contig_id : str
		The unique identifier of the sequence/contig.
	genome_id : str
		The unique identifier of the genome/file.
	protid : int
		The integer identifier to attriute to the first gene.
	genes : pyrodigal.Genes
		The list of genes predicted by Prodigal.

	Returns
	-------
	gene_info : list
		List with one sublist per gene predicted. Each sublist
		includes the sequence SHA256 hash, the DNA sequence, the
		genome identifier, the contig identifier, the start position
		in the sequence, the end position, the integer identifier and
		the strand the gene was identified in.
	protid : int
		The integer identifier to attribute to the first gene
		in the next sequence/contig.
	"""
	gene_info = []
	for gene in genes:
		sequence = gene.sequence()
		confidence = round(gene.confidence(), 2)
		sequence_hash = im.hash_sequence(sequence)
		# Store CDS IDs used by Pyrodigal
		pyrodigal_id = f'{genome_id}_{protid}'
		gene_info.append([sequence_hash, sequence, pyrodigal_id,
						  genome_id, contig_id, str(gene.begin), str(gene.end),
						  str(protid), str(gene.strand), str(confidence)])
		protid += 1

	return gene_info, protid


def predict_genome_genes(input_file, coordinates_outfile, output_directory, gene_finder, translation_table, output_formats):
	"""Predict genes for sequences in a FASTA file.

	Parameters
	----------
	input_file : str
		Path to the FASTA file.
	output_directory : str
		Path to the output_directory to store files with
		the results.
	gene_finder : pyrodigal.GeneFinder
		A GeneFinder object.
	translation_table : int
		Translation table used to configure the GeneFinder
		(None type if the GeneFinder does not need to be
		configured).

	Returns
	-------
	input_file : str
		Path to the input FASTA file.
	total_genome : int
		Total number of genes predicted.
	fasta_outfile : str
		Path to the output FASTA file that contains the
		predited gene sequences.
	coordinates_outfile : str
		Path to the output pickle file that contains the gene
		coordinates and contig size data.
	"""
	# Get genome unique identifier
	genome_basename = input_file[1]
	records = fao.sequence_generator(input_file[0])
	records = {rec.id: bytes(rec.seq) for rec in records}
	contig_sizes = {recid: len(sequence) for recid, sequence in records.items()}

	# Train based on input sequences
	# Only train if there is no GeneFinder object
	if gene_finder is not None:
		current_gene_finder = gene_finder
	else:
		current_gene_finder = create_gene_finder(None, True, True, False)
		training_info = train_gene_finder(current_gene_finder, records.values(), translation_table)

	# Predict genes for all input contigs
	contig_genes = {}
	for recid, sequence in records.items():
		genes = current_gene_finder.find_genes(sequence)
		contig_genes[recid] = genes

	# Extract data from Gene objects
	protid = 1
	gene_info = []
	# Store data about first and last CDSs in each sequence to speedup PLOT classification
	close_to_tip = {genome_basename: {}}
	for recid, genes in contig_genes.items():
		# Get coordinate data for CDSs
		data = get_gene_info(recid, genome_basename, protid, genes)
		gene_info.extend(data[0])
		# Get data for first and last CDS in the contig
		if len(data[0]) > 0:
			first_cds = data[0][0]
			close_to_tip[genome_basename].setdefault(first_cds[0], []).append((contig_sizes[first_cds[4]], int(first_cds[5]), int(first_cds[6]), first_cds[-1]))
			if first_cds != data[0][-1]:
				last_cds = data[0][-1]
				close_to_tip[genome_basename].setdefault(last_cds[0], []).append((contig_sizes[last_cds[4]], int(last_cds[5]), int(last_cds[6]), last_cds[-1]))
		# Reset protid based on the number of CDSs predicted for the sequence
		protid = data[1]

	# Get total number of CDSs predicted
	total_genome = len(gene_info)

	# Save data if Pyrodigal was able to predict genes
	output_files = [None, None, None, None, None]
	if total_genome > 0:
		if 'genes' in output_formats:
			# Define path to output FASTA
			fasta_outfile = fo.join_paths(output_directory, ['genes', f'{genome_basename}.fasta'])
			with open(fasta_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					# Use math.inf to write sequences in a single line
					genes.write_genes(outfile, sequence_id=genome_basename, width=math.inf)
			output_files[0] = (fasta_outfile)

		if 'translations' in output_formats:
			translations_outfile = fo.join_paths(output_directory, ['translations', f'{genome_basename}_translated.fasta'])
			with open(translations_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_translations(outfile, sequence_id=genome_basename, width=math.inf, include_stop=False)
			output_files[1] = translations_outfile

		if 'gff' in output_formats:
			gff_outfile = fo.join_paths(output_directory, ['gff', f'{genome_basename}.gff'])
			with open(gff_outfile, 'w') as outfile:
				i = 0
				for recid, genes in contig_genes.items():
					genes.write_gff(outfile, sequence_id=genome_basename, header=(i==0), include_translation_table=True)
					i += 1
			output_files[2] = gff_outfile

		if 'genbank' in output_formats:
			gbk_outfile = fo.join_paths(output_directory, ['genbank', f'{genome_basename}.gbk'])
			with open(gbk_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_genbank(outfile, sequence_id=genome_basename)
			output_files[3] = gbk_outfile

		if 'scores' in output_formats:
			scores_outfile = fo.join_paths(output_directory, ['scores', f'{genome_basename}.scores'])
			with open(scores_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_scores(outfile, sequence_id=genome_basename)
			output_files[4] = scores_outfile

		# Append gene coordinates to process TSV file
		coordinate_data = []
		for gene in gene_info:
			coordinate_data.append(gene[2:]+[gene[0]])
		coordinate_outlines = ['\t'.join(l) for l in coordinate_data]
		fo.write_lines(coordinate_outlines, coordinates_outfile, write_mode='a')

	return [input_file, total_genome, close_to_tip, output_files]


def write_coordinates_file(coordinates_file, output_file):
	"""Write the coordinates of CDSs identified in a genome assembly to a TSV file.

	Parameters
	----------
	coordinates_file : str
		Path to the pickle file that contains data about
		the CDSs coordinates.
	output_file : str
		Path to the output TSV file.
	"""
	data = fo.pickle_loader(coordinates_file)
	lines = [coords for h, coords in data[0].items()]
	lines = im.flatten_list(lines)
	# Sort lines by genome and protein ID (converting to int ensures natural sort)
	lines.sort(key=lambda x: (x[2], int(x[6])))
	lines = ['\t'.join(line) for line in lines]
	fo.write_lines(lines, output_file)


def predict_genes(fasta_files, output_directory, gene_prediction_parameters, cpu_cores):
	"""Predict genes from a set of input genome assemblies.

	Parameters
	----------
	fasta_files : list
		List of paths to FASTA files with genomic
		sequences.
	ptf_path : str
		Path to the Pyrodigal training file. Should
		be NoneType if a training file is not provided.
	translation_table : int
		Genetic code used to predict and translate
		coding sequences.
	pyrodigal_mode : str
		Pyrodigal execution mode.
	cpu_cores : int
		Number of processes that will run Pyrodigal in
		parallel.
	output_directory : str
		Path to the directory where output files
		with Pyrodigal's results will be stored in.

	Returns
	-------
	failed_info : list
		List that contains a list with the stderr for the
		cases that Pyrodigal failed to predict genes for
		and the path to the file with information about
		failed cases. Returns NoneType if gene prediction
		succeeded for all inputs.
	"""
	if gene_prediction_parameters["training_file"] is not None:
		# Read training file to create GeneFinder object
		training_data = read_training_file(gene_prediction_parameters["training_file"])
		# Create GeneFinder object based on training data
		gene_finder = create_gene_finder(training_data, True, True, False)
	elif gene_prediction_parameters["training_file"] is None and gene_prediction_parameters["mode"] == 'meta':
		# Create GeneFinder object to run in meta mode
		gene_finder = create_gene_finder(None, True, True, True)
	else:
		gene_finder = None

	pyrodigal_inputs = im.divide_list_into_n_chunks(fasta_files, len(fasta_files))
	# Add path to TSV file to store CDS coordinate data per input subset
	for i in range(len(pyrodigal_inputs)):
		pyrodigal_inputs[i].append(fo.join_paths(output_directory, [f'gene_coordinates_{i}']))

	# Create subfolders to store each output format
	for output_format in gene_prediction_parameters["output_formats"]:
		format_outdir = fo.join_paths(output_directory, [output_format])
		fo.create_directory(format_outdir)

	# Add common arguments to all sublists
	common_args = [output_directory, gene_finder, gene_prediction_parameters["translation_table"], gene_prediction_parameters["output_formats"]]
	pyrodigal_inputs = im.multiprocessing_inputs(pyrodigal_inputs, common_args, predict_genome_genes)

	# Run Pyrodigal to predict genes
	# Need to use ThreadPool. Pyrodigal might hang when using Pool
	pyrodigal_results = mo.map_async_parallelizer(pyrodigal_inputs,
												  mo.function_helper,
												  cpu_cores,
												  show_progress=True,
												  pool_type='threadpool')

	# Get number of inputs for which gene prediction failed
	# Inputs with 0 CDSs and inputs with error messages
	failed = {line[0][0]: line[1]
			  for line in pyrodigal_results
			  if line[1] == 0
			  or isinstance(line[1], str) is True}

	# Get number of CDSs predicted per valid input
	cds_counts = {line[0][1]: line[1]
				  for line in pyrodigal_results
				  if isinstance(line[1], int) is True}

	# Get total number of CDSs predicted
	total_cds = sum(cds_counts.values())

	# Get paths to FASTA files with the extracted CDSs
	cds_fastas = [line[-1][0] for line in pyrodigal_results if line[-1][0] is not None]

	# Merge TSV files with CDS coordinates
	coordinate_files = [i[1] for i in pyrodigal_inputs]
	merged_coordinates = fo.join_paths(output_directory, ['gene_coordinates.tsv'])
	fo.concatenate_files(coordinate_files, merged_coordinates, header=ct.GENE_TABLE_HEADER)
	# Delete separate coordinate files
	fo.remove_files(coordinate_files)

	# Merge dictionaries with info about CDSs close to contig tips
	close_to_tip = [line[2] for line in pyrodigal_results if len(line[2]) > 0]
	close_to_tip = im.merge_dictionaries(close_to_tip)

	return [failed, total_cds, cds_fastas, merged_coordinates, cds_counts, close_to_tip]
