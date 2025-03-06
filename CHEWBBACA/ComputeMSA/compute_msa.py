#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------



Code documentation
------------------
"""


import os
import sys

import pandas as pd

try:
	from utils import (
		constants as ct,
		mafft_wrapper as mw,
		file_operations as fo,
		fasta_operations as fao,
		iterables_manipulation as im,
		multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (
		constants as ct,
		mafft_wrapper as mw,
		file_operations as fo,
		fasta_operations as fao,
		iterables_manipulation as im,
		multiprocessing_operations as mo)


def profile_column_to_fasta(locus, input_file, schema_directory, output_directory):
	"""Create a FASTA file with the locus alleles identified in a dataset.

	Parameters
	----------
	locus : tuple
		A tuple with the locus identifier and column index.
	input_file : str
		Path to the TSV file with the allelic profiles.
	schema_directory : str
		Path to the schema directory.
	output_directory : str
		Path to the output directory.

	Returns
	-------
	fasta_file : str
		Path to the FASTA file that contains the allele sequences
		identified for each sample.
	"""
	# Read locus column
	df = pd.read_csv(input_file, usecols=[0, locus[1]],
					 delimiter='\t', dtype=str)
	locus_column = df[locus[0]].tolist()
	sample_ids = df[df.columns[0]].tolist()

	# Read locus schema FASTA file
	locus_fasta = fo.join_paths(schema_directory, [f'{locus[0]}.fasta'])
	alleles = fao.import_sequences(locus_fasta)
	# Map allele int IDs to sequence
	alleles = {k.split('_')[-1]: v for k, v in alleles.items()}

	failed = []
	sequences = []
	seq_hashes = set()
	for i, allele in enumerate(locus_column):
		# Remove INF- prefixes to be able to get alleles
		clean_class = allele.split('INF-')[-1]
		if clean_class in alleles:
			current_sample = sample_ids[i]
			sequence = alleles[clean_class]
			seqid = f'{locus[0]}_{clean_class}_{current_sample}'
			record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [seqid, sequence])
			sequences.append(record)
			seq_hashes.add(im.hash_sequence(sequence))
		# FASTA file does not contain the allele or it is a special classification
		# Using the --no-inferred option can cause this
		else:
			failed.append(clean_class)

	# Only create FASTA file if locus was identified in at least one sample
	if len(sequences) > 0:
		fasta_file = fo.join_paths(output_directory, [f'{locus[0]}.fasta'])
		fo.write_lines(sequences, fasta_file)
		return [locus, fasta_file, failed, len(sequences), len(seq_hashes)]
	else:
		return [locus, failed, len(sequences), len(seq_hashes)]


def concatenate_loci_alignments(sample, loci, sample_profile, fasta_index, output_directory):
	"""Concatenate the aligned sequences for a sample.

	Parameters
	----------
	sample : str
		Sample identifier.
	loci : list
		Loci identifiers.
	fasta_index : Bio.File._IndexedSeqFileDict
		Indexed FASTA file to get sequences from.
	output_directory : str
		Path to the output directory.

	Returns
	-------
	alignment_outfile : str
		Path to the FASTA file with the concatenated aligned sequences.
	"""
	alignment = ''
	for locus in loci:
		# Sequence headers include locus, allele IDs, and sample IDs joined by '_'
		# Get allele ID
		allele_id = sample_profile[locus].tolist()[0]
		allele_id = allele_id.replace('INF-', '')
		# Get aligned sequence from index
		try:
			seqid = f'{locus}_{allele_id}_{sample}'
			alignment += str(fasta_index[seqid].seq)
		except Exception as e:
			seqid = f'{locus}_0_{sample}'
			alignment += str(fasta_index[seqid].seq)
	# Save alignment for sample
	alignment_outfile = fo.join_paths(output_directory,
									  [f'{sample}_cgMLST_alignment.fasta'])
	alignment_record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
											[sample, alignment])
	fo.write_lines([alignment_record], alignment_outfile)

	return alignment_outfile


def add_gaps(input_file, locus_id, gap_char, sample_ids, output_directory):
	"""Add gap sequences to a MSA for samples that do not have an allele.

	Parameters
	----------
	input_file : str
		Path to the FASTA file with the aligned sequences.
	locus_id : str
		Locus identifier.
	gap_char : str
		Character to use as gap.
	sample_ids : list
		Sample identifiers.
	output_directory : str
		Path to the output directory.

	Returns
	-------
	gapped_fasta : str
		Path to the FASTA file containing the updated MSA.
	"""
	records = fao.import_sequences(input_file)
	# Get list of samples where locus was identified
	dataset_sids = {k.split('_')[-1]: k for k in records}
	# Get length of alignment to create gapped sequence
	msa_len = len(records[list(records.keys())[0]])
	gapped_seq = gap_char * msa_len
	gapped_records = {}
	for sid in sample_ids:
		# Sample contains locus
		if sid in dataset_sids:
			gapped_records[dataset_sids[sid]] = records[dataset_sids[sid]]
		# Sample does not contain locus
		# Add gapped sequence
		else:
			gapped_seq_sid = f'{locus_id}_0_{sid}'
			gapped_records[gapped_seq_sid] = gapped_seq

	# Save Fasta file with gapped sequences
	gapped_fasta = fo.join_paths(output_directory, [f'{locus_id}_gapped.fasta'])
	outrecords = [f'>{k}\n{v}' for k, v in gapped_records.items()]
	fo.write_lines(outrecords, gapped_fasta)
	
	return gapped_fasta


def convert_msa_to_dna(input_file, schema_file, locus_id, gap_char, output_directory):
	"""
	"""
	gapped_records = fao.import_sequences(input_file)
	schema_records = fao.import_sequences(schema_file)
	dna_records = {}
	for seqid, sequence in gapped_records.items():
		# Get allele ID
		allele_id = '_'.join(seqid.split('_')[:2])
		# Check if it matches any record in the schema
		if allele_id in schema_records:
			# Get allele sequence
			allele = schema_records[allele_id]
			# Iterate over gapped protein sequence to create gapped DNA
			dna_index = 0
			gapped_dna = ''
			for i, char in enumerate(sequence):
				# Add codon if it is not a gap
				if char != gap_char:
					gapped_dna += allele[dna_index:dna_index+3]
					dna_index += 3
				# Add '---' if it is a gap
				elif char == gap_char:
					gapped_dna += gap_char * 3
			dna_records[seqid] = gapped_dna
		else:
			dna_records[seqid] = sequence * 3

	# Save DNA MSA
	dna_fasta = fo.join_paths(output_directory, [f'{locus_id}_gapped_dna.fasta'])
	dna_msa_recs = [f'>{k}\n{v}' for k, v in dna_records.items()]
	fo.write_lines(dna_msa_recs, dna_fasta)

	return dna_fasta


input_file = '/home/rmamede/test_chewie/features/ComputeMSA/results_alleles.tsv'
schema_directory = '/home/rmamede/Chewie_Schemas/spyogenes_wgMLST'
output_directory = '/home/rmamede/test_chewie/features/ComputeMSA/test'
dna_msa = True
output_variable = False
gap_char = '-'
translation_table = 11
cpu_cores = 12
def main(input_file, schema_directory, output_directory, dna_msa, output_variable, gap_char, translation_table, cpu_cores):
	# Create output directory
	fo.create_directory(output_directory)
	# Get sample and loci IDs
	sample_ids = fo.extract_column(input_file, delimiter='\t', column_index=0)
	loci_ids = fo.read_lines(input_file, strip=True, num_lines=1)[0].split('\t')[1:]
	# Create FASTA files with the alleles identified in the dataset
	# Create temporary directory to store FASTA files
	print('Creating FASTA files with the alleles identified in the samples per locus...')
	fasta_dir = fo.join_paths(output_directory, ['fastas'])
	fo.create_directory(fasta_dir)
	# Divide into groups and process in parallel
	# Get loci indexes to read locus column with Pandas
	loci_indexes = [(l, loci_ids.index(l)+1) for l in loci_ids]
	inputs = im.divide_list_into_n_chunks(loci_indexes, len(loci_indexes))
	common_args = [input_file, schema_directory, fasta_dir]
	# Add common arguments to all sublists
	inputs = im.multiprocessing_inputs(inputs, common_args, profile_column_to_fasta)
	results = mo.map_async_parallelizer(inputs,
										mo.function_helper,
										cpu_cores,
										show_progress=True)

	# Save the number of alleles and distinct alleles identified per locus
	stats_lines = [f'{r[0][0]}\t{r[-2]}\t{r[-1]}' for r in results]
	stats_lines = [ct.COMPUTEMSA_STATS_HEADER] + stats_lines
	stats_file = fo.join_paths(output_directory, [ct.COMPUTEMSA_STATS_FILE])
	fo.write_lines(stats_lines, stats_file)

	# Get FASTA files for loci identified in at least one sample
	loci_to_msa = [r[1] for r in results if r[-1] > 0]
	if len(loci_to_msa) == 0:
		sys.exit(ct.COMPUTEMSA_NO_ALLELES)

	# Translate FASTA files
	print('Translating FASTA files...')
	translation_inputs = im.divide_list_into_n_chunks(loci_to_msa, len(loci_to_msa))
	common_args = [fasta_dir, translation_table]
	# Add common arguments to all sublists
	inputs = im.multiprocessing_inputs(translation_inputs, common_args, fao.translate_fasta)
	translation_results = mo.map_async_parallelizer(inputs,
												    mo.function_helper,
													cpu_cores,
													show_progress=True)

	##### Should compute MSA only for loci with more than one distinct allele
	protein_files = [r[1] for r in translation_results]

	# Run MAFFT to compute MSA
	print('\nRunning MAFFT to compute the MSA for each locus...')
	mafft_outdir = fo.join_paths(output_directory, ['mafft_results'])
	fo.create_directory(mafft_outdir)
	mafft_outfiles = [os.path.basename(file) for file in protein_files]
	mafft_outfiles = [file.replace('.fasta', '_aligned.fasta') for file in mafft_outfiles]
	mafft_outfiles = [fo.join_paths(mafft_outdir, [file]) for file in mafft_outfiles]
	mafft_inputs = [[file, mafft_outfiles[i]] for i, file in enumerate(protein_files)]
	common_args = []
	# Add common arguments to all sublists
	inputs = im.multiprocessing_inputs(mafft_inputs, common_args, mw.call_mafft)
	mafft_results = mo.map_async_parallelizer(inputs,
											  mo.function_helper,
											  cpu_cores,
											  show_progress=True)

	# Identify cases where MAFFT failed
	mafft_failed = [r[0] for r in mafft_results if r[1] is False]
	if len(mafft_failed) > 0:
		print(f'\nCould not determine MSA for {len(mafft_failed)} loci.')

	# Get files created by MAFFT
	mafft_success = [r[0] for r in mafft_results if r[1] is True]

	# Add gap sequences when sample did not have an allele
	gapped_inputs = []
	for file in mafft_success:
		locus_id = fo.file_basename(file, False).split('_protein')[0]
		gapped_inputs.append([file, locus_id])
		
	common_args = [gap_char, sample_ids, mafft_outdir]
	gapped_inputs = im.multiprocessing_inputs(gapped_inputs, common_args, add_gaps)
	gapped_results = mo.map_async_parallelizer(gapped_inputs,
											   mo.function_helper,
											   cpu_cores,
											   show_progress=True)
	
	# Delete original file with ungapped MSA
	fo.remove_files(mafft_success)

	# Convert protein MSAs to DNA MSAs
	if dna_msa:
		dna_inputs = []
		for file in gapped_results:
			locus_id = fo.file_basename(file).split('_gapped')[0]
			schema_file = fo.join_paths(schema_directory, [f'{locus_id}.fasta'])
			dna_inputs.append([file, schema_file, locus_id])
			
		common_args = [gap_char, mafft_outdir]
		dna_inputs = im.multiprocessing_inputs(dna_inputs, common_args, convert_msa_to_dna)
		dna_results = mo.map_async_parallelizer(dna_inputs,
										        mo.function_helper,
												cpu_cores,
												show_progress=True)

	# Identify variable positions to get SNP MSA
	if output_variable:
		pass

	# Create folder to store sample MSAs
	sample_msa_folder = fo.join_paths(output_directory, ['sample_alignments'])
	fo.create_directory(sample_msa_folder)

	# Create the full protein MSA
	print('\nCreating file with the full protein MSA...')
	# Concatenate all alignment files and index with BioPython
	protein_concat = fo.join_paths(sample_msa_folder, [ct.COMPUTEMSA_PROTEIN_CONCAT])
	fo.concatenate_files(gapped_results, protein_concat)
	# Index file
	indexed_protein_concat = fao.index_fasta(protein_concat)
	sample_alignment_files = []
	# Get loci IDs
	loci_ids = [fo.file_basename(file, False).split('_gapped')[0] for file in gapped_results]
	for sid in sample_ids:
		# Get sample profile
		sample_profile = pd.read_csv(input_file,
								     skiprows=range(1,sample_ids.index(sid)+1), # Only get header and sample profile
									 nrows=1,
									 delimiter='\t',
									 header=0,
									 dtype=str)
		alignment_file = concatenate_loci_alignments(sid,
													 loci_ids,
													 sample_profile,
													 indexed_protein_concat,
													 sample_msa_folder)
		sample_alignment_files.append(alignment_file)

	# Concatenate sample protein alignments
	full_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_PROTEIN_MSA])
	fo.concatenate_files(sample_alignment_files, full_alignment)

	if dna_msa:
		# Create the full DNA MSA
		print('Creating file with the full DNA MSA...')
		# Concatenate all alignment files and index with BioPython
		dna_concat = fo.join_paths(sample_msa_folder, [ct.COMPUTEMSA_DNA_CONCAT])
		fo.concatenate_files(dna_results, dna_concat)
		# Index file
		indexed_dna_concat = fao.index_fasta(dna_concat)
		sample_alignment_files = []
		# Get loci IDs
		loci_ids = [fo.file_basename(file, False).split('_gapped')[0] for file in dna_results]
		for sid in sample_ids:
			# Get sample profile
			sample_profile = pd.read_csv(input_file,
										 skiprows=range(1,sample_ids.index(sid)+1),
										 nrows=1,
										 delimiter='\t',
										 header=0,
										 dtype=str)
			alignment_file = concatenate_loci_alignments(sid,
														 loci_ids,
														 sample_profile,
														 indexed_dna_concat,
														 sample_msa_folder)
			sample_alignment_files.append(alignment_file)

		# Concatenate all cgMLST alignmnet records
		full_alignment = fo.join_paths(output_directory, [ct.COMPUTEMSA_DNA_MSA])
		fo.concatenate_files(sample_alignment_files, full_alignment)
