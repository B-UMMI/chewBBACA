#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to gene prediction
with Pyrodigal.

Code documentation
------------------
"""


import math

import Bio.SeqIO
import pyrodigal

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   fasta_operations as fao,
					   iterables_manipulation as im)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 file_operations as fo,
								 fasta_operations as fao,
								 iterables_manipulation as im)


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
	######### Need to improce training file type detection
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
		# Store CDS ID used by chewBBACA
		cds_id = f'{genome_id}-protein{protid}'
		gene_info.append([sequence_hash, sequence, cds_id, genome_id, contig_id,
						  str(gene.begin), str(gene.end), str(protid),
						  str(gene.strand), str(confidence)])
		protid += 1

	return gene_info, protid


def write_coordinates_pickle(gene_info, contig_sizes, output_file):
	"""Write gene coordinates to a pickle file.

	Parameters
	----------
	gene_info : list
		List with the data for the genes returned by `get_gene_info`.
	contig_sizes : dict
		Dictionary with contig/sequence identifiers as keys and
		contig/sequence size as values.
	output_file : str
	Path to the output file.
	"""
	gene_coordinates = {}
	for gene in gene_info:
		gene_coordinates.setdefault(gene[0], []).append(gene[2:])
	fo.pickle_dumper([gene_coordinates, contig_sizes], output_file)


def predict_genome_genes(input_file, output_directory, gene_finder,
						 translation_table, output_formats):
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
	contig_sizes = {recid: len(sequence)
					for recid, sequence in records.items()}

	# Train based on input sequences
	# Only train if there is no GeneFinder object
	if gene_finder is not None:
		current_gene_finder = gene_finder
	else:
		current_gene_finder = create_gene_finder(None, True, True, False)
		training_info = train_gene_finder(current_gene_finder,
												records.values(),
												translation_table)

	# Predict genes for all input contigs
	contig_genes = {}
	for recid, sequence in records.items():
		genes = current_gene_finder.find_genes(sequence)
		contig_genes[recid] = genes

	# Extract data from Gene objects
	protid = 1
	gene_info = []
	# Store data about first and last CDS in each sequence to speedup PLOT classification
	close_to_tip = {genome_basename: {}}
	for recid, genes in contig_genes.items():
		data = get_gene_info(recid, genome_basename, protid, genes)
		gene_info.extend(data[0])
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
	output_files = [None, None, None, None, None, None]
	if total_genome > 0:
		if 'genes' in output_formats:
			fasta_outfile = fo.join_paths(output_directory, [f'{genome_basename}.fasta'])
			with open(fasta_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_genes(outfile, sequence_id=genome_basename, width=math.inf)
			output_files[0] = fasta_outfile

		if 'translations' in output_formats:
			translations_outfile = fo.join_paths(output_directory, [f'{genome_basename}.translations'])
			with open(translations_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_translations(outfile, sequence_id=genome_basename, width=math.inf, include_stop=False)
			output_files[1] = translations_outfile

		if 'gff' in output_formats:
			gff_outfile = fo.join_paths(output_directory, [f'{genome_basename}.gff'])
			with open(gff_outfile, 'w') as outfile:
				i = 0
				for recid, genes in contig_genes.items():
					genes.write_gff(outfile, sequence_id=genome_basename, header=(i==0), include_translation_table=True)
					i += 1
			output_files[2] = gff_outfile

		if 'genbank' in output_formats:
			gbk_outfile = fo.join_paths(output_directory, [f'{genome_basename}.gbk'])
			with open(gbk_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_genbank(outfile, sequence_id=genome_basename)
			output_files[3] = gbk_outfile

		if 'scores' in output_formats:
			scores_outfile = fo.join_paths(output_directory, [f'{genome_basename}.scores'])
			with open(scores_outfile, 'w') as outfile:
				for recid, genes in contig_genes.items():
					genes.write_scores(outfile, sequence_id=genome_basename)
			output_files[4] = scores_outfile

		# Save gene coordinates and contig sizes to pickle
		coordinates_outfile = fo.join_paths(output_directory,
											[f'{genome_basename}_coordinates'])
		write_coordinates_pickle(gene_info, contig_sizes, coordinates_outfile)
		output_files[5] = coordinates_outfile

	return [input_file, total_genome, close_to_tip, output_files]
