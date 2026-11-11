#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module computes a multiple sequence alignment (MSA)
based on allele calling results. It uses the GetAlleles module
to extract allele sequences from a given schema for the samples
in the dataset. It then runs MAFFT to compute the MSA for each locus.
The final output is a concatenated MSA for all samples in the dataset.
The module can create both protein and DNA MSAs, as well as MSAs
containing only variable positions (SNPs). It also handles cases
where samples do not have alleles identified for certain loci by
adding gap sequences in the MSA. Additionally, it accepts the path to
a directory containing FASTA files as input, in which case it
computes the MSA directly from the sequences in the FASTA files.

Code documentation
------------------
"""


import os
import sys

import pandas as pd
from loguru import logger

import CHEWBBACA.global_constants as gct
import CHEWBBACA.ComputeMSA.constants as ct
from CHEWBBACA.GetAlleles import get_alleles
from CHEWBBACA.GetAlleles.validation_models import GetAllelesValidator
from CHEWBBACA.utils import (
	mafft_wrapper as mw,
	file_operations as fo,
	fasta_operations as fao,
	iterables_manipulation as im,
	multiprocessing_operations as mo)


def concatenate_loci_alignments(sample, loci, sample_profile, fasta_index, output_directory):
	"""Concatenate the aligned sequences for a sample.

	Parameters
	----------
	sample : str
		Sample identifier for which to create the concatenated alignment.
	loci : list
		Identifiers of the loci to include in the concatenated alignment.
	fasta_index : Bio.File._IndexedSeqFileDict
		Indexed FASTA file to get sequences from.
	output_directory : str
		Path to the output directory where the FASTA file containing the
		sample MSA will be saved.

	Returns
	-------
	alignment_outfile : str
		Path to the FASTA containing the sample MSA.
	"""
	alignment = ''
	for locus in loci:
		# Sequence headers include sample, locus, and allele IDs joined by '_'
		# Get allele ID
		allele_id = sample_profile[locus].tolist()[0]
		# Remove '*' from allele ID
		allele_id = allele_id.replace('*', '')
		# Get aligned sequence from index
		try:
			seqid = f'{sample}_{locus}_{allele_id}'
			alignment += str(fasta_index[seqid].seq)
		except Exception as e:
			seqid = f'{sample}_{locus}_0'
			alignment += str(fasta_index[seqid].seq)
	# Save FASTA file with sample alignment
	alignment_outfile = fo.join_paths(output_directory,
									  [f'{sample}.fasta'])
	alignment_record = fao.fasta_str_record(gct.FASTA_RECORD_TEMPLATE,
											[sample, alignment])
	fo.write_lines([alignment_record], alignment_outfile)

	return alignment_outfile


def add_gaps(input_file, locus_id, output_file, gap_char, sample_ids):
	"""Add gap sequences to a MSA for samples that do not have a locus.

	Parameters
	----------
	input_file : str
		Path to the FASTA file with the aligned sequences.
	locus_id : str
		Locus identifier.
	gap_char : str
		Character used to fill gaps.
	sample_ids : list
		Sample identifiers.

	Returns
	-------
	output_file : str
		Path to the FASTA file containing the MSA updated with gap sequences.
	"""
	records = fao.import_sequences(input_file)
	# Get list of samples where locus was identified
	dataset_sids = {k.split(f'_{locus_id}')[0]: k for k in records}
	# Get length of alignment to create gapped sequence
	# Just get length of the first record
	msa_len = len(records[list(records.keys())[0]])
	gap_seq = gap_char * msa_len
	gap_records = {}
	for sid in sample_ids:
		# Sample contains locus
		if sid in dataset_sids:
			gap_records[dataset_sids[sid]] = records[dataset_sids[sid]]
		# Sample does not contain locus
		# Add gap sequence
		else:
			gap_seq_sid = f'{sid}_{locus_id}_0'
			gap_records[gap_seq_sid] = gap_seq

	# Save Fasta file with gap sequences
	outrecords = [f'>{k}\n{v}' for k, v in gap_records.items()]
	fo.write_lines(outrecords, output_file)

	return output_file


def convert_msa_to_dna(input_file, dna_file, output_file, gap_char):
	"""Convert a protein MSA to a DNA MSA.

	Parameters
	----------
	input_file : str
		Path to the FASTA file with the protein MSA.
	dna_file : str
		Path to the FASTA file with the DNA sequences.
	output_file : str
		Path to the FASTA file where the DNA MSA will be saved.
	gap_char : str
		Character used to fill gaps.

	Returns
	-------
	output_file : str
		Path to the FASTA file containing the DNA MSA.
	"""
	# Import aligned protein sequences
	protein_records = fao.import_sequences(input_file)
	# Import DNA sequences (unaligned)
	dna_sequences = fao.import_sequences(dna_file)
	# Create DNA MSA based on protein MSA
	# by substituting each amino acid by its codon
	dna_records = {}
	for seqid, sequence in protein_records.items():
		# Check if it matches any record in the schema
		if seqid in dna_sequences:
			# Get allele sequence
			allele = dna_sequences[seqid]
			# Iterate over gapped protein sequence to create gapped DNA
			dna_index = 0
			gap_dna = ''
			for i, char in enumerate(sequence):
				# Add codon if it is not a gap
				if char != gap_char:
					gap_dna += allele[dna_index:dna_index+3]
					dna_index += 3
				# Add '---' if it is a gap
				elif char == gap_char:
					gap_dna += gap_char * 3
			dna_records[seqid] = gap_dna
		# If no matching record found, create a gap sequence
		else:
			dna_records[seqid] = sequence * 3

	# Save DNA MSA
	dna_msa_recs = [f'>{k}\n{v}' for k, v in dna_records.items()]
	fo.write_lines(dna_msa_recs, output_file)

	return output_file


def determine_variable_positions(msa_file, output_file, gap_char, ambiguous_chars, gaps, ambiguous):
	"""Determine the variable positions in a MSA.

	Parameters
	----------
	msa_file : str
		Path to the FASTA file with the complete MSA.
	output_file : str
		Path to the FASTA file where the MSA with only variable positions will be saved.
	gap_char : str
		Character used to fill gaps.
	ambiguous_chars : list
		List of characters considered ambiguous.
	gaps : str
		How to handle gaps. Options: 'ignore', 'exclude'. If 'ignore', positions
		with gaps are considered but gaps are not counted as a distinct character.
		If 'exclude', positions with gaps are not considered.
	ambiguous : str
		How to handle ambiguous characters. Options: 'ignore', 'exclude'. If 'ignore',
		positions with ambiguous characters are considered but ambiguous characters
		are not counted as distinct characters. If 'exclude', positions with ambiguous
		characters are not considered.

	Returns
	-------
	Returns a list with two elements:
	- bool : True if variable positions were found, False otherwise.
	- output_file : str
		Path to the FASTA file containing the MSA with only variable positions,
		if there are any. Returns the input msa_file if no variable positions were found.
	"""
	# Read MSA
	msa_records = fao.import_sequences(msa_file)
	msa_sequences = list(msa_records.values())
	# Zip to pair chars in the same position
	zipped = list(zip(*msa_sequences))
	variable = []
	for i, pos in enumerate(zipped):
		# Get unique characters in position
		distinct = set(pos)

		# Ignore positions if any sequence has a gap or an ambiguous char
		if gaps == 'ignore' and gap_char in distinct:
			distinct.remove('-')
		if ambiguous == 'ignore' and any([c in distinct for c in ambiguous_chars]):
			for c in ambiguous_chars:
				if c in distinct:
					distinct.remove(c)

		# Exclude positions if any sequence has a gap or an ambiguous char
		# Simply continue to next position
		if gaps == 'exclude' and gap_char in distinct:
			continue
		if ambiguous == 'exclude' and any([c in distinct for c in ambiguous_chars]):
			continue

		# If there are more than 1 unique character, it is a variable position
		if len(distinct) > 1:
			variable.append([i, pos])

	# Save variable positions
	# Create FASTA records with only variable positions
	# If no variable positions, do not create file
	if len(variable) > 0:
		variable_pos = list(zip(*[i[1] for i in variable]))
		variable_pos = [''.join(i) for i in variable_pos]
		variable_records = [[allele_id, variable_pos[i]] for i, allele_id in enumerate(msa_records.keys())]
		variable_records = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, variable_records)
		fo.write_lines(variable_records, output_file)
		return [True, output_file]
	else:
		return [False, msa_file]


def create_sample_msa(sample_id, sample_index, fasta_file, output_directory, loci_ids, profiles):
	"""Concatenate loci alignments for a given sample.

	Parameters
	----------
	sample_id : str
		Sample identifier for which to create the concatenated alignment.
	sample_index : list
		Sample column index in the profiles file.
	fasta_file : str
		Path to a FASTA file containing the aligned sequences for all loci.
	output_directory : str
		Path to the output directory where the FASTA file containing the
		sample MSA will be saved.
	loci_ids : list
		Identifiers of the loci to include in the concatenated alignment.
	profiles : str
		Path to the TSV file containing the allelic profiles.

	Returns
	-------
	alignment_file : str
		Path to the FASTA containing the sample MSA.
	"""
	# This tries to create an index for the FASTA file
	# but ends up loading an existing index created with SeqIO.index_db
	fasta_index = fao.index_fasta(fasta_file, True)
	# Get sample profile
	sample_profile = pd.read_csv(profiles,
								 skiprows=sample_index, # Only get header and sample profile
								 nrows=1,
								 delimiter='\t',
								 header=0,
								 dtype=str)

	# Remove all 'INF-' prefixes and '*' from identifiers
	# Replace special classifications by '0'
	masked_profile = sample_profile.apply(im.replace_chars)
	alignment_file = concatenate_loci_alignments(sample_id,
												 loci_ids,
												 masked_profile,
												 fasta_index,
												 output_directory)

	return alignment_file


def main(config):
	# Create output directory
	fo.create_directory(config.output_directory)

	if os.path.isfile(config.input_path):
		logger.info("Input is a TSV file with allelic profiles.")
		if config.schema_directory is None:
			sys.exit(ct.COMPUTEMSA_NO_SCHEMA)
		# Get sample IDs
		sample_ids = fo.extract_column(config.input_path, delimiter="\t", column_index=0)
		# Create FASTA files with the alleles identified in the dataset
		# Call GetAlleles module
		# Removes the '*' in allele IDs
		logger.info("Calling the GetAlleles module to create FASTA files with the allele sequences for the samples in the dataset...")
		# Create config object pass to GetAlleles
		get_alleles_outdir = fo.join_paths(config.output_directory, ["alleles"])
		get_alleles_config = GetAllelesValidator(output_directory=get_alleles_outdir,
		                                         allelic_profiles=config.input_path,
												 schema_directory=config.schema_directory,
												 cpu_cores=config.cpu_cores,
												 translate_alleles=True,
												 translation_table=config.translation_table,
												 distinct=False)
		_, dna_files, protein_files = get_alleles.main(get_alleles_config)
	elif os.path.isdir(config.input_path):
		logger.info("Input is a directory with FASTA files.")
		if config.schema_directory is not None:
			logger.info("Warning: Schema directory provided will be ignored since input is a directory with FASTA files.")
		# Input is a directory with FASTA files
		# Copy input files to output directory
		logger.info("Copying FASTA files to temp directory...")
		files_to_copy = fo.listdir(config.input_path, full_path=True)
		files_to_copy = fo.filter_by_extension(files_to_copy, gct.FASTA_EXTENSIONS)[0]
		copied_dir = fo.join_paths(config.output_directory, ["input_files"])
		fo.create_directory(copied_dir)
		copied_files = []
		for file in files_to_copy:
			fo.copy_file(file, copied_dir)
			copied_files.append(fo.join_paths(copied_dir, [os.path.basename(file)]))

		# Define list of files with DNA and protein sequences
		# If input files contain DNA sequences, translate to protein
		if config.protein_input is False:
			dna_files = copied_files
			translated_outdir = fo.join_paths(config.output_directory, ["translated"])
			fo.create_directory(translated_outdir)
			protein_files = [fao.translate_fasta(file, translated_outdir, config.translation_table, False) for file in copied_files]
			# Do not store result if files included no alleles that could be translated
			protein_files = [result[1] for result in protein_files if result[1] is not None]
		# If input files contain protein sequences
		else:
			dna_files = []
			protein_files = copied_files

	# Exit if no sequences were found
	if len(dna_files) == 0 and len(protein_files) == 0:
		sys.exit(ct.NO_MSAS_CREATED)

	# Run MAFFT to compute protein MSAs
	logger.info("Running MAFFT to compute the MSA for each input file...")
	mafft_outdir = fo.join_paths(config.output_directory, ["MSAs"])
	mafft_protein_outdir = fo.join_paths(mafft_outdir, ["protein"])
	fo.create_directory(mafft_protein_outdir)
	mafft_outfiles = [os.path.basename(file) for file in protein_files]
	mafft_outfiles = [fo.join_paths(mafft_protein_outdir, [file]) for file in mafft_outfiles]
	mafft_inputs = [[file, mafft_outfiles[i]] for i, file in enumerate(protein_files)]
	common_args = [config.custom_mafft_options.split()] if config.custom_mafft_options is not None else [gct.MAFFT_DEFAULT_PARAMETERS[:-1]]
	# Add common arguments to all sublists
	inputs = im.multiprocessing_inputs(mafft_inputs, common_args, mw.call_mafft)
	mafft_results = mo.map_async_parallelizer(inputs, mo.function_helper, config.cpu_cores)

	# Identify cases where MAFFT failed
	mafft_failed = [r[0] for r in mafft_results if r[1] is False]
	if len(mafft_failed) > 0:
		logger.info(f"Could not determine MSA for {len(mafft_failed)} loci.")

	# Get files created by MAFFT
	mafft_success = [r[0] for r in mafft_results if r[1] is True]
	if len(mafft_success) == 0:
		sys.exit(ct.NO_MSAS_CREATED)

	# Add gap sequences when samples did not have an allele identified
	# Only do this if input is a TSV file with allelic profiles
	# If it is a folder with FASTA files there is no information about samples missing loci
	if not os.path.isdir(config.input_path):
		logger.info("Adding gap sequences for samples missing loci...")
		# This will create new files with MSAs containing gap sequences (samples missing loci will show as gap sequences)
		gapped_inputs = []
		for file in mafft_success:
			locus_id = fo.file_basename(file, False).split("_protein")[0]
			gapped_inputs.append([file, locus_id, file])
		common_args = [ct.GAP_CHAR, sample_ids]
		gapped_inputs = im.multiprocessing_inputs(gapped_inputs, common_args, add_gaps)
		gapped_results = mo.map_async_parallelizer(gapped_inputs, mo.function_helper, config.cpu_cores)
	# If the input was a directory with FASTA files, there is no need to add gap sequences
	else:
		gapped_results = mafft_success

	# Convert protein MSAs to DNA MSAs
	if config.dna_msa:
		dna_inputs = []
		if config.protein_input is False:
			logger.info("Converting protein MSAs to DNA MSAs...")
			aligned_dna_outdir = fo.join_paths(mafft_outdir, ["dna"])
			fo.create_directory(aligned_dna_outdir)
			for i, file in enumerate(gapped_results):
				locus_id = fo.file_basename(file, False).split("_protein")[0]
				outfile_dna = fo.join_paths(aligned_dna_outdir, [f"{locus_id}.fasta"])
				dna_inputs.append([file, dna_files[i], outfile_dna])

			# Convert to DNA MSA
			common_args = [ct.GAP_CHAR]
			dna_inputs = im.multiprocessing_inputs(dna_inputs, common_args, convert_msa_to_dna)
			dna_results = mo.map_async_parallelizer(dna_inputs, mo.function_helper, config.cpu_cores)
		else:
			logger.info("DNA MSAs were requested, but input files contain protein sequences. Skipping DNA MSA creation.")

	# Delete directories containing the FASTA files with DNA and protein sequences
	if not config.no_cleanup:
		if config.protein_input is False:
			fo.delete_directory(os.path.dirname(dna_files[0]))
		fo.delete_directory(os.path.dirname(protein_files[0]))

	# Identify variable positions to get SNP MSA
	if config.output_variable:
		logger.info("Determining variable positions in the protein MSAs...")
		# Create folder to store variable position MSAs
		variable_protein_outfolder = fo.join_paths(mafft_outdir, ['variable_protein'])
		fo.create_directory(variable_protein_outfolder)
		# Define paths to output files
		variable_protein_outfiles = [fo.join_paths(variable_protein_outfolder, [fo.file_basename(file)]) for file in gapped_results]
		common_args = [ct.GAP_CHAR, gct.PROTEIN_AMBIGUOUS_CHARS, config.gaps, config.ambiguous]
		variable_inputs = im.multiprocessing_inputs([list(i) for i in zip(gapped_results, variable_protein_outfiles)], common_args, determine_variable_positions)
		variable_results = mo.map_async_parallelizer(variable_inputs, mo.function_helper, config.cpu_cores)

		# Collect variable and non-variable protein MSA files
		variable_protein = [r[1] for r in variable_results if r[0] is True]
		non_variable_protein = [r[1] for r in variable_results if r[0] is False]

		logger.info(f"{len(variable_protein)} loci include variable positions at the protein level.")

		# Determine variable positions for DNA MSAs
		if config.dna_msa:
			if config.protein_input is False:
				logger.info("Determining variable positions in the DNA MSAs...")
				# Create folder to store variable position MSAs
				variable_dna_outfolder = fo.join_paths(mafft_outdir, ["variable_dna"])
				fo.create_directory(variable_dna_outfolder)
				variable_dna_outfiles = [fo.join_paths(variable_dna_outfolder, [fo.file_basename(file)]) for file in dna_results]
				common_args = [ct.GAP_CHAR, gct.PROTEIN_AMBIGUOUS_CHARS, config.gaps, config.ambiguous]
				variable_inputs = im.multiprocessing_inputs([list(i) for i in zip(dna_results, variable_dna_outfiles)], common_args, determine_variable_positions)
				variable_results = mo.map_async_parallelizer(variable_inputs, mo.function_helper, config.cpu_cores)

				# Collect variable and non-variable DNA MSA files
				variable_dna = [r[1] for r in variable_results if r[0] is True]
				non_variable_dna = [r[1] for r in variable_results if r[0] is False]
				
				logger.info(f"{len(variable_dna)} loci include variable positions at the DNA level.")
			else:
				logger.info("DNA MSAs were requested, but input files contain protein sequences. Skipping variable position determination for DNA MSAs.")

	# User only wants the loci MSAs or the input is a folder with FASTA files
	if config.only_loci_msas or os.path.isdir(config.input_path):
		logger.info(f"MSAs for each input file are available in {mafft_outdir}")
		return mafft_outdir

	# Create folder to store sample MSAs
	sample_msas_outdir = fo.join_paths(config.output_directory, ["sample_MSAs"])
	fo.create_directory(sample_msas_outdir)

	# Create the full protein MSA
	logger.info("Creating file with the full protein MSA...")
	sample_protein_msas_outdir = fo.join_paths(sample_msas_outdir, ["protein"])
	fo.create_directory(sample_protein_msas_outdir)
	# Get loci IDs for which a gapped protein MSA was created
	gapped_results_loci_ids = [fo.file_basename(file).split("_protein")[0] for file in gapped_results]
	# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
	loci_msa_concat = fo.join_paths(sample_protein_msas_outdir, ["concat.fasta"])
	# Concatenate all alignment files and index with BioPython
	fo.concatenate_files(gapped_results, loci_msa_concat)

	# Index file with SeqIO.index_db to store record information as a file on disk
	# This allows to reload the index with multiprocessing
	# SeqIO.index creates the index in memory which cannot be shared between processes
	loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

	# Pass the FASTA file used to create the index to each process and the function reloads the index for each process
	common_args = [loci_msa_concat, sample_protein_msas_outdir, gapped_results_loci_ids, config.input_path]
	inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
	inputs = im.multiprocessing_inputs(inputs, common_args, create_sample_msa)
	sample_protein_MSA_outfiles = mo.map_async_parallelizer(inputs, mo.function_helper, config.cpu_cores)

	# Concatenate sample protein alignments
	full_protein_alignment = fo.join_paths(config.output_directory, [ct.COMPUTEMSA_PROTEIN_MSA])
	fo.concatenate_files(sample_protein_MSA_outfiles, full_protein_alignment)

	# Get length of the full alignment
	protein_msa_length = len(fo.read_lines(full_protein_alignment, strip=True, num_lines=2)[1])
	logger.info(f"Protein MSA length: {protein_msa_length}")

	# Delete file indexed to get individual sequences
	fo.remove_files([loci_msa_concat])

	if not config.no_cleanup:
		fo.delete_directory(sample_protein_msas_outdir)

	if config.dna_msa:
		sample_dna_msas_outdir = fo.join_paths(sample_msas_outdir, ["dna"])
		fo.create_directory(sample_dna_msas_outdir)

		# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
		loci_msa_concat = fo.join_paths(sample_dna_msas_outdir, ["concat.fasta"])
		# Concatenate all alignment files and index with BioPython
		fo.concatenate_files(dna_results, loci_msa_concat)

		# Index file with SeqIO.index_db to store record information as a file on disk
		# This allows to reload the index with multiprocessing
		# SeqIO.index creates the index in memory which cannot be shared between processes
		loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

		# Create the full DNA MSA
		logger.info("Creating file with the full DNA MSA...")
		common_args = [loci_msa_concat, sample_dna_msas_outdir, gapped_results_loci_ids, config.input_path]
		inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
		inputs = im.multiprocessing_inputs(inputs, common_args, create_sample_msa)
		sample_dna_MSA_outfiles = mo.map_async_parallelizer(inputs, mo.function_helper, config.cpu_cores)

		# Concatenate sample dna alignments
		full_dna_alignment = fo.join_paths(config.output_directory, [ct.COMPUTEMSA_DNA_MSA])
		fo.concatenate_files(sample_dna_MSA_outfiles, full_dna_alignment)

		# Get length of the full alignment
		dna_msa_length = len(fo.read_lines(full_dna_alignment, strip=True, num_lines=2)[1])
		logger.info(f"DNA MSA length: {dna_msa_length}")

		# Delete folder with intermediate sample MSAs
		if not config.no_cleanup:
			fo.delete_directory(sample_dna_msas_outdir)

	# Create full MSAs with only variable positions if requested
	if config.output_variable and len(variable_protein) > 0:
		sample_protein_variable_msas_folder = fo.join_paths(sample_msas_outdir, ["variable_protein"])
		fo.create_directory(sample_protein_variable_msas_folder)

		# Get loci IDs for which a gapped protein MSA was created and had variable positions
		variable_results_loci_ids = [fo.file_basename(file).split("_protein")[0] for file in variable_protein]

		logger.info(f"{len(variable_results_loci_ids)} loci include variable positions.")

		# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
		loci_msa_concat = fo.join_paths(sample_protein_variable_msas_folder, ["concat.fasta"])
		# Concatenate all alignment files and index with BioPython
		fo.concatenate_files(variable_protein, loci_msa_concat)

		# Index file with SeqIO.index_db to store record information as a file on disk
		# This allows to reload the index with multiprocessing
		# SeqIO.index creates the index in memory which cannot be shared between processes
		loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

		# Create the full protein MSA for the variable positions
		logger.info("Creating file with the full protein MSA for the variable positions...")
		common_args = [loci_msa_concat, sample_protein_variable_msas_folder, variable_results_loci_ids, config.input_path]
		inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
		inputs = im.multiprocessing_inputs(inputs, common_args, create_sample_msa)
		sample_protein_variable_MSA_outfiles = mo.map_async_parallelizer(inputs, mo.function_helper, config.cpu_cores)

		# Concatenate sample protein alignments
		full_variable_protein_alignment = fo.join_paths(config.output_directory, [ct.COMPUTEMSA_PROTEIN_MSA_VARIABLE])
		fo.concatenate_files(sample_protein_variable_MSA_outfiles, full_variable_protein_alignment)

		# Get length of the full alignment
		protein_variable_msa_length = len(fo.read_lines(full_variable_protein_alignment, strip=True, num_lines=2)[1])
		logger.info(f"Protein variable MSA length: {protein_variable_msa_length}")

		# Delete folder with intermediate sample protein variable MSAs
		if not config.no_cleanup:
			fo.delete_directory(sample_protein_variable_msas_folder)

	if config.output_variable and config.dna_msa and len(variable_dna) > 0:
		sample_dna_variable_msas_folder = fo.join_paths(sample_msas_outdir, ["variable_dna"])
		fo.create_directory(sample_dna_variable_msas_folder)

		# Get loci IDs for which a gapped DNA MSA was created and had variable positions
		variable_results_loci_ids = [fo.file_basename(file, False) for file in variable_dna]

		# Concatenate all loci MSAs and index with BioPython to retrieve sequences easily
		loci_msa_concat = fo.join_paths(sample_dna_variable_msas_folder, ["concat.fasta"])
		# Concatenate all alignment files and index with BioPython
		fo.concatenate_files(variable_dna, loci_msa_concat)

		# Index file with SeqIO.index_db to store record information as a file on disk
		# This allows to reload the index with multiprocessing
		# SeqIO.index creates the index in memory which cannot be shared between processes
		loci_msa_concat_index = fao.index_fasta(loci_msa_concat, True)

		# Create the full DNA MSA for the variable positions
		logger.info("Creating file with the full DNA MSA for the variable positions...")
		common_args = [loci_msa_concat, sample_dna_variable_msas_folder, variable_results_loci_ids, config.input_path]
		inputs = [[sid, range(1,sample_ids.index(sid)+1)] for sid in sample_ids]
		inputs = im.multiprocessing_inputs(inputs, common_args, create_sample_msa)
		sample_dna_variable_MSA_outfiles = mo.map_async_parallelizer(inputs, mo.function_helper, config.cpu_cores)

		# Concatenate sample dna alignments
		full_variable_dna_alignment = fo.join_paths(config.output_directory, [ct.COMPUTEMSA_DNA_MSA_VARIABLE])
		fo.concatenate_files(sample_dna_variable_MSA_outfiles, full_variable_dna_alignment)

		# Get length of the full alignment
		dna_variable_msa_length = len(fo.read_lines(full_variable_dna_alignment, strip=True, num_lines=2)[1])
		logger.info(f"DNA variable MSA length: {dna_variable_msa_length}")

		# Delete folder with intermediate sample DNA variable MSAs
		if not config.no_cleanup:
			fo.delete_directory(sample_dna_variable_msas_folder)

	# Delete folders with gapped MSAs for full and only variable positions
	if not config.no_cleanup:
		fo.delete_directory(get_alleles_outdir)
		fo.delete_directory(mafft_outdir)
		fo.delete_directory(sample_msas_outdir)

	logger.info(f"Results are available in {config.output_directory}")
