#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module performs allele calling to determine the allelic profiles
for a set of bacterial strains.

Code documentation
------------------
"""


import os
import re
import csv
import sys
import math
from collections import Counter

try:
	from utils import (constants as ct,
					   blast_wrapper as bw,
					   profile_hasher as ph,
					   core_functions as cf,
					   file_operations as fo,
					   fasta_operations as fao,
					   process_datetime as pdt,
					   sequence_manipulation as sm,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 blast_wrapper as bw,
								 profile_hasher as ph,
								 core_functions as cf,
								 file_operations as fo,
								 fasta_operations as fao,
								 process_datetime as pdt,
								 sequence_manipulation as sm,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def compute_loci_modes(loci_files, output_file):
	"""Determine the allele size mode for a set of loci.

	Parameters
	----------
	loci_files : list
		List with the full paths to loci FASTA files.
	output_file : str
		Path to the output file created to store the allele size mode
		values (created with the Pickle module).

	Returns
	-------
	loci_modes : str
		Path to the output file with the allele size mode values (a
		dictionary with loci identifiers as keys and the allele size
		mode and the list of allele sizes as values).
	"""
	loci_modes = {}
	for file in loci_files:
		locus_id = fo.file_basename(file, False)
		allele_sizes = list(fao.sequence_lengths(file).values())
		# Select first if there are several values with same frequency
		loci_modes[locus_id] = [sm.determine_mode(allele_sizes)[0], allele_sizes]

	fo.pickle_dumper(loci_modes, output_file)

	return loci_modes


def create_hash_table(fasta_files, output_file, translation_table):
	"""Create a hash table with allele hashes mapped to allele identifiers.

	Parameters
	----------
	fasta_files : list
		List with paths to loci FASTA files.
	output_file : str
		Path to the output file.
	translation_table : int
		Genetic code used to translate the alleles if the hashes should be
		determined based on the aminoacid sequences. Will use the DNA
		sequences if this argument is not provided.

	Returns
	-------
	output_file : str
		Path to the output file (a dictionary with allele hashes as keys
		and the locus integer index and allele identifier as values is
		saved with the Pickle module).
	"""
	hashtable = {}
	for file in fasta_files:
		locus_index = file[1]
		alleles = [(rec.id, str(rec.seq))
				   for rec in fao.sequence_generator(file[0])]
		if translation_table is not None:
			alleles = [(rec[0], str(sm.translate_sequence(rec[1], translation_table)))
					   for rec in alleles]

		for record in alleles:
			# Needs to be string because of '*' added to schemas from Chewie-NS
			allele_id = record[0].split('_')[-1]
			allele_hash = im.hash_sequence(record[1])
			hashtable.setdefault(allele_hash, []).append((locus_index, allele_id))

	fo.pickle_dumper(hashtable, output_file)

	return output_file


def precompute_hash_tables(output_directory, loci_files, translation_table,
						   cpu_cores, max_sequences=ct.HASH_TABLE_MAXIMUM_ALLELES):
	"""Create allele and translated allele hash tables.

	Parameters
	----------
	output_directory : str
		Path to the output directory.
	loci_files : list
		List with paths to loci FASTA files.
	translation_table : int
		Genetic code used to translate the alleles if the hashes should be
		determined based on the aminoacid sequences. Will use the DNA
		sequences if this argument is not provided.
	cpu_cores : int
		Number of CPU cores used to compute the hash tables.
	max_sequences : int
		Maximum number of sequences per hash table.

	Returns
	-------
	hash_tables : list
		List with the paths for the files that contain the hash tables (
		files are created in the output diretory with the Pickle module).
	"""
	input_groups = []
	current_group = []
	total_alleles = 0
	for file in loci_files:
		# Count number of sequences per file
		num_alleles = fao.count_sequences(file)
		# Increment total sequences and reset if it reaches limit
		total_alleles += num_alleles
		current_group.append((file, loci_files.index(file)))
		if total_alleles >= max_sequences or file == loci_files[-1]:
			input_groups.append(current_group)
			current_group = []
			total_alleles = 0

	hash_tables = []
	# Create inputs to parallelize hash table creation
	# Hash tables for DNA sequences
	output_files = [fo.join_paths(output_directory, [f"{'DNAtable'}{i+1}"])
					for i, g in enumerate(input_groups)]

	inputs = [[g, output_files[i], None, create_hash_table]
			  for i, g in enumerate(input_groups)]

	dna_tables = mo.map_async_parallelizer(inputs,
										   mo.function_helper,
										   cpu_cores,
										   show_progress=False)
	hash_tables.append(dna_tables)

	# Hash tables for translated alleles
	output_files = [fo.join_paths(output_directory, [f"{'PROTEINtable'}{i+1}"])
					for i, g in enumerate(input_groups)]

	inputs = [[g, output_files[i], translation_table, create_hash_table]
			  for i, g in enumerate(input_groups)]

	protein_tables = mo.map_async_parallelizer(inputs,
											   mo.function_helper,
											   cpu_cores,
											   show_progress=False)
	hash_tables.append(protein_tables)

	return hash_tables


def update_hash_tables(loci_files, loci_to_call, translation_table,
					   pre_computed_dir):
	"""Update pre-computed hash tables.

	Parameters
	----------
	loci_files : dict
		Dictionary with paths to schema loci FASTA files as keys and
		a list with the paths to the temporary FASTA files that contain
		the newly inferred alleles for the loci as values.
	loci_to_call : dict
		Dictionary with paths to schema loci FASTA files as keys and
		loci integer identifiers as values.
	translation_table : int
		Genetic code used to translate the alleles.
	pre_computed_dir : str
		Path to the directory that contains the files with the pre-computed
		hash tables.

	Returns
	-------
	Total number of new allele hashes added to the pre-computed hash tables.
	"""
	# Create hash tables with data for new alleles
	novel_dna_hashes = {}
	novel_protein_hashes = {}
	for key, value in loci_files.items():
		# Get locus integer identifier
		locus_index = loci_to_call[key]
		# Import new alleles
		records = fao.sequence_generator(value[0])
		for rec in records:
			# Compute allele SHA256 hash
			allele_id = (rec.id).split('_')[-1]
			sequence = str(rec.seq)
			seq_hash = im.hash_sequence(sequence)
			prot_hash = im.hash_sequence(str(sm.translate_sequence(sequence, translation_table)))
			# Add to hash tables that will be used to update pre-computed data
			novel_dna_hashes.setdefault(seq_hash, []).append((locus_index, allele_id))
			novel_protein_hashes.setdefault(seq_hash, []).append(prot_hash)

	# Update pre-computed hash tables
	# List files with pre-computed hash tables and select latest
	dna_tables = fo.listdir_fullpath(pre_computed_dir, 'DNAtable')
	prot_tables = fo.listdir_fullpath(pre_computed_dir, 'PROTEINtable')
	latest_dna_table = sorted(dna_tables,
							  key=lambda x: int(x.split('table')[-1]))[-1]
	latest_prot_table = sorted(prot_tables,
							   key=lambda x: int(x.split('table')[-1]))[-1]
	current_dna_table = fo.pickle_loader(latest_dna_table)
	current_prot_table = fo.pickle_loader(latest_prot_table)
	# Update hash tables
	for key, value in novel_dna_hashes.items():
		# Add entry for each new allele
		current_dna_table.setdefault(key, []).extend(value)
		current_prot_table.setdefault(novel_protein_hashes[key][0], []).extend(value)
		# Check if it reached the maximum number of alleles per table file
		if len(current_dna_table) >= ct.HASH_TABLE_MAXIMUM_ALLELES:
			# Save current hash tables and create new files
			fo.pickle_dumper(current_dna_table, latest_dna_table)
			fo.pickle_dumper(current_prot_table, latest_prot_table)
			new_index = int(latest_dna_table.split('table')[-1]) + 1
			latest_dna_table = fo.join_paths(pre_computed_dir,
											 ['DNAtable{0}'.format(new_index)])
			current_dna_table = {}
			latest_prot_table = fo.join_paths(pre_computed_dir,
											  ['PROTEINtable{0}'.format(new_index)])
			current_prot_table = {}

	if len(current_dna_table) > 0:
		fo.pickle_dumper(current_dna_table, latest_dna_table)
		fo.pickle_dumper(current_prot_table, latest_prot_table)

	return len(novel_dna_hashes)


def update_classification(genome_id, locus_results, match_info):
	"""Update locus classification for an input.

	Parameters
	----------
	genome_id : int
		Integer identifier attributed to the input.
	locus_results : dict
		Dictionary with the matches found for the locus
		in the inputs.
	match_info : list
		List with information about the match found for
		the locus.

	Returns
	-------
	locus_results : dict
		Updated results.
	"""
	# Add data about match
	locus_results.setdefault(genome_id, [match_info[3]]).append(match_info)

	# Get all classifications
	classes_list = [c[3] for c in locus_results[genome_id][1:]]
	# Evaluate classification for genomes with multiple matches
	if len(classes_list) > 1:
		classes_counts = Counter(classes_list)
		# Multiple matches, single class
		if len(classes_counts) == 1:
			if 'EXC' in classes_counts:
				locus_results[genome_id][0] = 'NIPHEM'
			# Multiple INF, ASM, ALM, etc classes are classified as NIPH
			else:
				locus_results[genome_id][0] = 'NIPH'
		# Multiple matches and classes
		elif len(classes_counts) > 1:
			# Inputs that include both EXC and INF are classified as NIPH
			if 'EXC' in classes_counts and 'INF' in classes_counts:
				locus_results[genome_id][0] = 'NIPH'
			# Any class plus PLOT3, PLOT5 or LOTSC are classified as NIPH
			elif any([c in ['PLOT3', 'PLOT5', 'LOTSC'] for c in classes_counts]) is True:
				locus_results[genome_id][0] = 'NIPH'
			# EXC or INF with ASM/ALM
			elif 'EXC' in classes_counts or 'INF' in classes_counts:
				match_count = classes_counts.get('EXC', classes_counts['INF'])
				# Single EXC or INF classified as EXC or INF even if there are ASM/ALM
				if match_count == 1:
					locus_results[genome_id][0] = 'EXC' if 'EXC' in classes_counts else 'INF'
				# Multiple EXC or INF classified as NIPH
				else:
					locus_results[genome_id][0] = 'NIPH'
			# Multiple ASM and ALM are classified as NIPH
			else:
				locus_results[genome_id][0] = 'NIPH'

	return locus_results


def count_global_classifications(classification_files, classification_labels):
	"""Determine counts for each classification type.

	Parameters
	----------
	classification_files : list
		List of paths to pickled files that contain the
		classifications for a set of loci.
	classification_labels : list
		List with the possible class labels assigned by
		chewBBACA.

	Returns
	-------
	classification_counts : dict
		Dicitonary with classification types as keys
		and the total number of inputs classified per
		type as values.
	total_cds : int
		The total number of coding sequences that
		were classified.
	"""
	total_cds = 0
	classification_counts = Counter()
	for file in classification_files:
		locus_results = fo.pickle_loader(file)
		total_cds += sum([len([r for r in c if isinstance(r, tuple)])
						  for g, c in locus_results.items()])
		locus_classifications = [c[0] for g, c in locus_results.items()]
		locus_counts = Counter(locus_classifications)
		classification_counts += locus_counts

	# Add classification types that might be missing
	classification_counts.update(Counter({k: 0 for k in classification_labels[:-1]
										  if k not in classification_counts}))

	return [classification_counts, total_cds]


def dna_exact_matches(table_file, cds_hashtable, loci_files,
					  classification_files, input_ids):
	"""Find exact matches at DNA level.

	Parameters
	----------
	table_file : str
		Path to the file with the pre-computed SHA-256 hashes for the
		alleles in the schema.
	cds_hashtable : dict
		Dictionary with SHA-256 hashes for distinct DNA sequences
		extracted from the inputs and lists of genome integer
		identifiers enconded with the polyline algorithm as values.
	loci_files : dict
		Dictionary with the loci integer identifiers as keys and the
		paths to the loci FASTA files as values.
	classification_files : dict
		Dictionary with the paths to the loci FASTA files as keys and
		the paths to the files with classification results as values.
	input_ids : dict
		Dictionary with input integer identifiers as keys
		and input sequence identifiers as values.

	Returns
	-------
	matched_seqids : list
		Sequence identifiers of the distinct CDSs that matched
		alleles in the schema.
	total_matches : int
		Total number of exact matches.
	distinct_seqids_total : int
		Total number of distinct CDSs that matched alleles in
		the schema
	"""
	hashes_table = fo.pickle_loader(table_file)
	matches = set(hashes_table.keys()).intersection(set(cds_hashtable.keys()))
	matches_data = {m: hashes_table[m] for m in matches}
	# Sort based on locus integer identifier
	sorted_matches = sorted(list(matches_data.items()), key=lambda x: x[1][0][0])
	# Group matches per locus
	loci_matches = {}
	for match in sorted_matches:
		# Only get information about the first locus that contains alleles
		# Does not get additional info if allele is present in several loci
		loci_matches.setdefault(match[1][0][0], []).append((match[0], match[1][0][1]))

	loci_matches = {loci_files[k]: v for k, v in loci_matches.items() if k in loci_files}

	total_matches = 0
	matched_seqids = []
	distinct_seqids_total = 0
	for g, results in loci_matches.items():
		locus_file = classification_files[g]
		locus_classifications = fo.pickle_loader(locus_file)

		locus_seqids = []
		locus_total_matches = 0
		for r in results:
			# Decode list of inputs that contain allele
			matched_inputs = im.polyline_decoding(cds_hashtable[r[0]])
			# Get seqid chosen as representative during sequence deduplication
			representative_seqid = '{0}-protein{1}'.format(input_ids[matched_inputs[1]], matched_inputs[0])
			match_data = (r[1], representative_seqid, r[0], 'EXC', 1.0)
			# Classify as exact matches
			# Skip first value in list, it is the protein identifier
			for gid in matched_inputs[1:]:
				locus_classifications = update_classification(gid,
															  locus_classifications,
															  match_data)

			locus_total_matches += len(matched_inputs[1:])
			# Store representative id for the sequences
			locus_seqids.append(representative_seqid)

		# Save updated classifications
		fo.pickle_dumper(locus_classifications, locus_file)
		# Extend list of matched seqids
		matched_seqids.extend(locus_seqids)
		# Increment number of distinct CDSs
		distinct_seqids_total += len(locus_seqids)
		# Increment total number of CDSs
		total_matches += locus_total_matches

	return [matched_seqids, total_matches, distinct_seqids_total]


def protein_exact_matches(table_file, translated_hashtable, loci_files,
						  classification_files, input_ids, distinct_cds_index,
						  cds_hashtable, previous_hashes):
	"""Find exact matches at protein level.

	Parameters
	----------
	table_file : str
		Path to the file with the pre-computed SHA-256 hashes for the
		alleles in the schema.
	translated_hashtable : dict
		Dictionary with SHA-256 hashes for distinct protein sequences
		extracted from the inputs and lists of sequence identifiers
		enconded with the polyline algorithm as values.
	loci_files : dict
		Dictionary with the loci integer identifiers as keys and the
		paths to the loci FASTA files as values.
	classification_files : dict
		Dictionary with the paths to the loci FASTA files as keys and
		the paths to the files with classification results as values.
	input_ids : dict
		Dictionary with input integer identifiers as keys
		and input sequence identifiers as values.
	distinct_cds_index : Bio.File._IndexedSeqFileDict
		Biopython index for the Fasta file with distinct CDSs.
	cds_hashtable : dict
		Dictionary with SHA-256 hashes for distinct DNA sequences
		extracted from the inputs and lists of genome integer
		identifiers enconded with the polyline algorithm as values.
	previous_hashes : set
		Set of protein hashes that were matched in the previous
		function calls.

	Returns
	-------
	matched_seqids : list
		Sequence identifiers of the distinct CDSs that matched.
	distinct_seqids_total : int
		Total number of distinct CDSs that matched.
	total_matches : int
		Total number of CDSs that matched.
	distinct_protids_total : int
		Total number of distinct protids that matched.
	inferred_lengths : dict
		Dictionary with loci identifiers as keys and lists with
		the size of newly inferred alleles as values. Used to
		update the loci mode values.
	previous_hashes : set
		Set of protein hashes that were matched in the previous
		function calls and in the current call.
	"""
	hashes_table = fo.pickle_loader(table_file)
	common = set(hashes_table.keys()).intersection(set(translated_hashtable.keys()))
	common_data = {c: hashes_table[c] for c in common}
	sorted_common = sorted(list(common_data.items()), key=lambda x: x[1][0][0])
	# Group matches by locus
	loci_matches = {}
	for r in sorted_common:
		# Only gets information about the first locus and first match
		# Does not get additional info if translated allele matches multiple
		# alleles in a single locus or multiple alleles in different loci
		loci_matches.setdefault(r[1][0][0], []).append((r[0], r[1][0][1]))

	loci_matches = {loci_files[k]: v for k, v in loci_matches.items() if k in loci_files}

	# Total number of matches
	total_matches = 0
	matched_seqids = []
	# Total number of distinct seqids matched
	distinct_seqids_total = 0
	# Total number of distinct protids matched
	distinct_protids_total = 0
	inferred_lengths = {}
	for g, results in loci_matches.items():
		locus_file = classification_files[g]
		locus_classifications = fo.pickle_loader(locus_file)
		# Total locus matches
		locus_total_matches = 0
		# List of distinct seqids matched
		locus_distinct_seqids = []
		# Total number of distinct protids matched
		locus_distinct_protids_total = 0
		# Set of CDS hashes that matched
		locus_matched_cds_hashes = set()
		# Store length of new alleles to update mode value
		locus_inferred_lengths = []
		for r in results:
			# Do not proceed if the protein hash was in another pre-computed table
			if r[0] in translated_hashtable and r[0] not in previous_hashes:
				# Get protids for distinct DNA CDSs
				matched_protids = im.polyline_decoding(translated_hashtable[r[0]])
				matched_protids = ['{0}-protein{1}'.format(input_ids[matched_protids[i+1]], matched_protids[i])
								   for i in range(0, len(matched_protids), 2)]
				locus_distinct_seqids.extend(matched_protids)
				locus_distinct_protids_total += 1
				# For each distinct CDS that codes for the protein
				for m in matched_protids:
					cds = str(distinct_cds_index.get(m).seq)
					cds_hash = im.hash_sequence(cds)
					# Get IDs of genomes that contain the CDS
					matched_inputs = im.polyline_decoding(cds_hashtable[cds_hash])
					locus_total_matches += len(matched_inputs)-1
					# For each genome ID that contains the CDS
					for gid in matched_inputs[1:]:
						# First time seeing CDS
						if cds_hash not in locus_matched_cds_hashes:
							current_class = 'INF'
							locus_matched_cds_hashes.add(cds_hash)
							locus_inferred_lengths.append(len(cds))
						else:
							current_class = 'EXC'
						# For protein exact matches, the seqid of the translated allele,
						# the seqid of the protein chosen as representative during sequence deduplication,
						match_data = (r[1], m, cds_hash, current_class, 1.0)
						locus_classifications = update_classification(gid, locus_classifications,
																	  match_data)

		# Save updated results for locus
		fo.pickle_dumper(locus_classifications, locus_file)
		matched_seqids.extend(locus_distinct_seqids)
		distinct_seqids_total += len(locus_distinct_seqids)
		total_matches += locus_total_matches
		distinct_protids_total += locus_distinct_protids_total
		# Update locus mode
		if len(locus_inferred_lengths) > 0:
			locus_id = fo.file_basename(g, False)
			inferred_lengths[locus_id] = locus_inferred_lengths

	# Save hashes to avoid classifying a protein hash that is in multiple PROTEINtables
	previous_hashes = previous_hashes.union(common)

	return [matched_seqids, distinct_seqids_total, total_matches,
			distinct_protids_total, inferred_lengths, previous_hashes]


def contig_position_classification(representative_length, representative_leftmost_pos,
								   representative_rightmost_pos, contig_length,
								   contig_leftmost_pos, contig_rightmost_pos, strand):
	"""Determine classification based on the alignment position on the contig.

	Parameters
	----------
	representative_length : int
		Length of the representative allele that matched a
		CDS identified in an input.
	representative_leftmost_pos : int
		Representative sequence leftmost aligned position.
	representative_rightmost_pos : int
		Representative sequence rightmost aligned position.
	contig_length : int
		Length of the contig that contains the coding sequence
		that matched with the representative allele.
	contig_leftmost_pos : int
		Contig leftmost aligned position.
	contig_rightmost_pos : int
		Contig rightmost aligned position.

	Returns
	-------
	'LOTSC' if the contig is smaller than the matched representative
	allele, 'PLOT5' or 'PLOT3' if the matched allele unaligned part
	exceeds one of the contig ends, None otherwise.
	"""
	# Check if contig is smaller than the representative allele (LOTSC)
	if contig_length < representative_length:
		return 'LOTSC'

	# Check if it is a PLOT
	# Determine rightmost aligned position in contig
	contig_rightmost_rest = contig_length - contig_rightmost_pos
	# Determine leftmost aligned position in contig
	contig_leftmost_rest = contig_leftmost_pos
	# Determine number of rightmost bases in the target that did not align
	representative_rightmost_rest = representative_length - representative_rightmost_pos
	# Determine number of leftmost bases in the target that did not align
	representative_leftmost_rest = representative_leftmost_pos

	# Determine if the unaligned region of the representative allele exceeds
	# one of the contig ends
	if strand == '1':
		if contig_leftmost_rest < representative_leftmost_rest:
			return 'PLOT5'
		elif contig_rightmost_rest < representative_rightmost_rest:
			return 'PLOT3'
	# CDS is in the reverse strand
	else:
		if contig_rightmost_rest < representative_leftmost_rest:
			return 'PLOT5'
		elif contig_leftmost_rest < representative_rightmost_rest:
			return 'PLOT3'


def allele_size_classification(sequence_length, locus_mode, size_threshold):
	"""Determine classification based on sequence size varaition.

	Parameters
	----------
	sequence_length : int
		Length of the CDS.
	locus_mode : int
		Locus allele size mode.
	size_threshold : float
		Sequence size variation threshold.

	Returns
	-------
	'ASM' if sequence size value is below computed sequence
	size interval, 'ALM' if it is above and None if it is
	contained in the interval.
	"""
	if size_threshold is not None:
		if sequence_length < (locus_mode[0]-(locus_mode[0])*size_threshold):
			return 'ASM'
		elif sequence_length > (locus_mode[0]+(locus_mode[0])*size_threshold):
			return 'ALM'


def write_loci_summary(classification_files, output_directory, total_inputs,
					   classification_labels, loci_finder):
	"""Write a TSV file with classification counts per locus.

	Parameters
	----------
	classification_files : dict
		Dictionary with the paths to loci FASTA files as keys
		and paths to loci classification files as values.
	output_directory : str
		Path to the output directory where the TSV file will
		be created.
	total_inputs : int
		Total number of inputs.
	classification_labels : list
		List with the possible class labels assigned by
		chewBBACA.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.


	Returns
	-------
	output_file : str
		Path to the output file.
	"""
	loci_stats = [ct.LOCI_STATS_HEADER]
	for k, v in classification_files.items():
		locus_id = loci_finder.search(k).group()
		locus_results = fo.pickle_loader(v)

		# Count locus classifications
		current_counts = count_global_classifications([v], classification_labels)
		counts_list = [locus_id]
		for c in classification_labels[:-1]:
			counts_list.append(str(current_counts[0][c]))
		# Add LNF count
		counts_list.append(str(total_inputs-len(locus_results)))
		# Add total number of classified CDSs
		counts_list.append(str(current_counts[1]))
		locus_line = im.join_list(counts_list, '\t')
		loci_stats.append(locus_line)

	output_file = fo.join_paths(output_directory, [ct.LOCI_STATS_BASENAME])
	fo.write_lines(loci_stats, output_file)

	return output_file


def write_logfile(start_time, end_time, total_inputs,
				  total_loci, cpu_cores, blast_score_ratio,
				  output_directory):
	"""Write the log file.

	Parameters
	----------
	start_time : datetime.datetime
		Datetime object with the date and hour
		determined when the process started running.
	end_time : datetime.datetime
		Datetime object with the date and hour
		determined when the process concluded.
	total_inputs : int
		Number of inputs passed to the process.
	total_loci : int
		Number of schema loci used for allele calling.
	cpu_cores : int
		Number of CPU cores/threads used by the
		process.
	blast_score_ratio : float
		BLAST Score Ratio value used by the
		process.
	output_directory : str
		Path to the output directory where the
		log file will be created.

	Returns
	-------
	log_outfile : str
		Path to the log file.
	"""
	start_time_str = pdt.datetime_str(start_time,
									  date_format='%H:%M:%S-%d/%m/%Y')

	end_time_str = pdt.datetime_str(end_time,
									date_format='%H:%M:%S-%d/%m/%Y')

	log_outfile = fo.join_paths(output_directory, [ct.LOGFILE_BASENAME])
	logfile_text = ct.LOGFILE_TEMPLATE.format(start_time_str, end_time_str,
											  total_inputs, total_loci,
											  cpu_cores, blast_score_ratio)

	fo.write_to_file(logfile_text, log_outfile, 'w', '')

	return log_outfile


def write_results_alleles(classification_files, input_identifiers,
						  output_directory, missing_class, loci_finder):
	"""Write a TSV file with the allelic profiles for the input samples.

	Parameters
	----------
	classification_files : list
		List with the paths to loci classification files.
	input_identifiers : list
		Sorted list that contains input string identifiers.
	output_directory : str
		Path to the output directory.
	missing_class : str
		'LNF' if execution mode is 4, 'PLNF' otherwise.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	output_file : str
		Path to the output file.
	"""
	# Limit the number of values to store in memory
	values_limit = ct.RESULTS_MAXVALS
	# Define intermediate file path
	intermediate_file = fo.join_paths(output_directory,
									  ['inter_results_alleles.tsv'])
	# Add first column with input identifiers
	columns = [['FILE'] + input_identifiers]
	for i, file in enumerate(classification_files):
		# Get locus identifier to add as column header
		locus_id = loci_finder.search(file).group()
		locus_results = fo.pickle_loader(file)
		locus_column = [locus_id]
		for ri in range(1, len(input_identifiers)+1):
			# Determine if locus was found in each input
			if ri in locus_results:
				current_result = locus_results[ri]
				# Exact or inferred, append assigned allele id
				if current_result[0] in ['EXC', 'INF']:
					locus_column.append(current_result[-1])
				# Ambiguous classes (PLOT, ASM, ALM, ...)
				else:
					locus_column.append(current_result[0])
			# Locus was not identified in the input
			else:
				locus_column.append(missing_class)

		columns.append(locus_column)

		if (len(columns)*len(input_identifiers)) >= values_limit or (i+1) == len(classification_files):
			inter_lines = [im.join_list(c, '\t') for c in columns]
			fo.write_lines(inter_lines, intermediate_file, write_mode='a')
			columns = []

	# Transpose intermediate file
	transposed_file = fo.transpose_matrix(intermediate_file, output_directory)
	# Change basename of transposed file
	output_file = fo.join_paths(output_directory, [ct.RESULTS_ALLELES_BASENAME])
	fo.move_file(transposed_file, output_file)

	# Delete intermediate files
	fo.remove_files([intermediate_file])

	return output_file


def write_results_statistics(classification_files, input_identifiers,
							 cds_counts, output_directory, classification_labels,
							 repeated_counts, invalid_data):
	"""Write a TSV file with classification counts per input.

	Parameters
	----------
	classification_files : dict
		Dictionary with the paths to loci FASTA files as keys
		and paths to loci classification files as values.
	input_identifiers : dict
		Dictionary with input integer identifiers as keys
		and input string identifiers as values.
	cds_counts : dict
		Dictionary with input integer identifiers as keys
		and the number of CDSs identified in each input as values.
	output_directory : str
		Path to the output directory where the TSV file will
		be created.
	classification_labels : list
		List with the possible class labels assigned by
		chewBBACA.
	repeated_counts : dict
		Dictionary with input identifiers as keys and the total
		number of times a CDS matches multiple loci -1 as values.
	invalid_data : dict
		Dictionary with input identifiers as keys and the total
		number of invalid CDSs as values.

	Returns
	-------
	Path to the output file.
	"""
	# Store total number of classified CDSs per input
	classified = {i: 0 for i in input_identifiers}
	# Initialize classification counts per input
	class_counts = {i: {c: 0 for c in classification_labels}
					for i in input_identifiers}
	for file in classification_files.values():
		locus_results = fo.pickle_loader(file)

		for i in class_counts:
			if i in locus_results:
				class_counts[i][locus_results[i][0]] += 1
				# Increment classified CDSs count
				classified[i] += len([c for c in locus_results[i][1:] if isinstance(c, tuple) is True])
			else:
				class_counts[i][classification_labels[-1]] += 1

	for i in class_counts:
		class_counts[i]['Classified_CDSs'] = classified[i]
		class_counts[i]['Total_CDSs'] = cds_counts[i]

	# Substitute integer identifiers by string identifiers
	class_counts = {input_identifiers[i]: v for i, v in class_counts.items()}
	# Subtract number of times a CDS is repeated and add invalid CDSs per input
	for i in class_counts:
		class_counts[i]['Classified_CDSs'] -= repeated_counts[i]
		class_counts[i]['Invalid CDSs'] = 0 if invalid_data is None else invalid_data[1][i]
	# Initialize with header line
	header_line = ['FILE'] + classification_labels + ['Invalid CDSs', 'Classified_CDSs', 'Total_CDSs']
	lines = [header_line]
	for i, v in class_counts.items():
		input_line = [i] + [str(v[c]) for c in header_line[1:]]
		lines.append(input_line)

	outlines = ['\t'.join(line) for line in lines]
	output_file = fo.join_paths(output_directory, [ct.RESULTS_STATISTICS_BASENAME])
	fo.write_lines(outlines, output_file)

	return output_file


def write_results_contigs(classification_files, input_identifiers,
						  output_directory, cds_coordinates_files,
						  classification_labels, loci_finder):
	"""Write a TSV file with the CDS coordinates for each input.

	Writes a TSV file with coding sequence coordinates (contig
	identifier, start and stop positions, and coding strand) for
	EXC and INF classifications or with the classification type
	if it is not EXC or INF.

	Parameters
	----------
	classification_files : list
		List with the paths to loci classification files.
	input_identifiers : dict
		Dictionary with input integer identifiers as keys
		and input string identifiers as values.
	output_directory : str
		Path to the output directory where the TSV file will
		be created.
	cds_coordinates_files : dict
		Dictionary with input string identifiers as keys
		and paths to pickled files with coding sequence
		coordinates as values.
	classification_labels : list
		List with the class labels attributed by chewBBACA.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	output_file : str
		Path to the output file that contains the sequence
		coordinates.
	repeated_info : dict
		Dictionary with hashes for the CDSs that matched multiple loci
		as keys and a list with information about the genome of origin
		and matched loci as values.
	repeated_counts : dict
		Dictionary with input identifiers as keys and the total
		number of times a CDS matches multiple loci -1 as values.
	"""
	# Do not include EXC, INF and LNF/PLNF classes
	invalid_classes = classification_labels[2:-1]
	intermediate_file = fo.join_paths(output_directory,
									  ['inter_results_contigsInfo.tsv'])
	columns = [['FILE'] + list(input_identifiers.values())]
	# Limit the number of values to store in memory
	values_limit = ct.RESULTS_MAXVALS
	# Get hash if coordinates are available, seqid otherwise
	id_index = 2 if cds_coordinates_files is not None else 1
	for i, file in enumerate(classification_files):
		locus_id = loci_finder.search(file).group()
		locus_results = fo.pickle_loader(file)
		column = [locus_id]
		for gid in input_identifiers:
			# Get sequence hash for exact and inferred
			if gid in locus_results and locus_results[gid][0] not in invalid_classes:
				column.append(locus_results[gid][1][id_index])
			# Get classification for other cases or LNF/PLNF for no classification
			else:
				column.append(locus_results.get(gid, [classification_labels[-1]])[0])

		columns.append(column)

		if (len(columns)*len(input_identifiers)) >= values_limit or (i+1) == len(classification_files):
			inter_lines = [im.join_list(c, '\t') for c in columns]
			fo.write_lines(inter_lines, intermediate_file, write_mode='a')
			columns = []

	# Transpose intermediate file
	transposed_file = fo.transpose_matrix(intermediate_file, output_directory)

	# Include LNF class to exclude from counting and coordinate retrieval
	invalid_classes.append(classification_labels[-1])

	# Get all hashes that match more than one locus per input
	# Alleles that match more than one locus are only identified if they are
	# classified as EXC or INF for both loci
	# If they are classified as EXC/INF for one locus and NIPH for the other,
	# they will not be identified based on the results for that genome
	# We need to get all alleles that match multiple loci before deciding if
	# they should be added to the schema.
	repeated_hashes, repeated_counts = fo.count_repeated_matrix(transposed_file, invalid_classes)

	# Define path to output file
	output_file = fo.join_paths(output_directory, [ct.RESULTS_COORDINATES_BASENAME])
	# Create dictionary to store info about repeated CDSs
	repeated_info = {}
	with open(transposed_file, 'r') as infile:
		csv_reader = csv.reader(infile, delimiter='\t')
		header = csv_reader.__next__()
		output_lines = [header]
		# Convert hashes to CDS coordinates
		for i, l in enumerate(csv_reader):
			genome_id = l[0]
			coordinates = {}
			if cds_coordinates_files is not None:
				# Open file with loci coordinates
				coordinates = fo.pickle_loader(cds_coordinates_files[genome_id])[0]
				# Start position is 0-based, stop position is upper-bound exclusive
				# Convert to PAMA CDSs that matched multiple loci
			cds_coordinates = []
			for j, c in enumerate(l[1:]):
				current_coordinates = coordinates.get(c, [c])[0]
				# Contig identifier, start and stop positions and strand
				# 1 for sense, -1 for reverse
				coordinates_str = (c if current_coordinates in invalid_classes or isinstance(current_coordinates, list) is False
								   else '{0}&{1}-{2}&{3}'.format(*current_coordinates[1:4], current_coordinates[5]))
				if c not in repeated_hashes:
					cds_coordinates.append(coordinates_str)
				else:
					cds_coordinates.append(classification_labels[-2])
					repeated_info.setdefault(c, []).append([genome_id, header[j+1], coordinates_str])

			output_lines.append([genome_id]+cds_coordinates)

			if (len(output_lines)*len(classification_files)) >= values_limit or (i+1) == len(input_identifiers):
				output_lines = ['\t'.join(l) for l in output_lines]
				fo.write_lines(output_lines, output_file, write_mode='a')
				output_lines = []

	# Delete intermediate files
	fo.remove_files([intermediate_file, transposed_file])

	return [output_file, repeated_info, repeated_counts]


def create_unclassified_fasta(fasta_file, prot_file, unclassified_protids,
							  protein_hashtable, output_directory, inv_map):
	"""Write a FASTA file with the CDSs that were not classified.

	Parameters
	----------
	fasta_file : str
		Path to FASTA file that contains the distinct coding
		sequences identified in the inputs.
	prot_file : str
		Path to FASTA file that contains the distinct translated
		coding sequences identified in the inputs.
	unclassified_protids : list
		List with the sequence identifiers of the representative
		sequences that were not classified.
	protein_hashtable : dict
		Dictionary with SHA-256 hashes for distinct DNA
		sequences extracted from the inputs and lists of
		genome integer identifiers enconded with the
		polyline algorithm as values.
	output_directory : str
		Path to the output directory where the file will be
		created.
	inv_map : dict
		Dictionary with input integer identifiers as keys
		and input string identifiers as values.

	Returns
	-------
	output_file : str
		Path to the output FASTA file.
	"""
	# Mode > 1, there is a FASTA file for CDSs and another for translated CDSs
	# We need to get the distinct seqids from the list of unclassified protids
	if fasta_file != prot_file:
		unclassified_seqids = []
		prot_distinct_index = fao.index_fasta(prot_file)
		for protid in unclassified_protids:
			prot_seq = str(prot_distinct_index[protid].seq)
			# Compute SHA-256 hash
			prot_hash = im.hash_sequence(prot_seq)
			# Get all seqids for the distinct CDSs that code for the protein
			seqids = im.polyline_decoding(protein_hashtable[prot_hash])
			# Create CDS identifiers to fetch from FASTA file
			seqids = ['{0}-protein{1}'.format(inv_map[seqids[i+1]], seqids[i])
					  for i in range(0, len(seqids), 2)]
			unclassified_seqids.extend(seqids)
	# Mode 1 or no valid CDSs after translation
	# There is no FASTA file with translated CDSs
	# We can fetch sequences directly
	else:
		unclassified_seqids = unclassified_protids

	output_file = fo.join_paths(output_directory, [ct.UNCLASSIFIED_BASENAME])
	dna_index = fao.index_fasta(fasta_file)
	# Create FASTA file with unclassified CDSs
	fao.get_sequences_by_id(dna_index, unclassified_seqids, output_file)

	return output_file


def assign_allele_ids(locus_files, ns, repeated, output_directory, loci_finder):
	"""Assign allele identifiers to CDSs classified as EXC or INF.

	Parameters
	----------
	locus_files : list
		List with the path to the locus FASTA file and the path to
		the locus classification file.
	ns : bool
		If the schema was downloaded from Chewie-NS. If True,
		adds '*' before new allele integer identifiers.
	repeated : list
		List of allele SHA256 hashes for the alleles that matched
		multiple loci. Those alleles are not assigned identifers
		and the classification changes to PAMA (PAralogous MAtch).

	Returns
	-------
	novel_alleles : dict
		Dictionary with paths to loci FASTA files as keys and
		lists with SHA-256 hashes and allele integer identifiers
		for each novel allele.
	"""
	# Dictionary to store new allele identifiers
	novel_alleles = []
	total_novel = 0
	# Import allele calling results
	locus_id = loci_finder.search(locus_files[0]).group()
	locus_results = fo.pickle_loader(locus_files[1])
	# Sort by input order
	sorted_results = sorted(locus_results.items(), key=lambda x: x[0])
	# Only keep INF and EXC classifications
	sorted_results = [r for r in sorted_results if r[1][0] in ['EXC', 'INF']]
	# Sort to get INF classifications first
	sorted_results = sorted(sorted_results,
							key=lambda x: x[1][0] == 'INF',
							reverse=True)
	# Add * to new alleles added to schemas downloaded from Chewie-NS
	alleleid_prefix = '*' if ns is True else ''
	# Assign allele identifiers if there are novel alleles
	if len(sorted_results) > 0:
		# Import locus records in the schema
		records = fao.import_sequences(locus_files[0])
		# Map allele SHA-256 hash to allele integer ID
		matched_alleles = {im.hash_sequence(v): k.split('_')[-1]
						   for k, v in records.items()}
		# Get greatest allele integer identifier
		# Remove * from alleles not yet synchronized to Chewie-NS 
		max_alleleid = max([int(rec.replace('*', '').split('_')[-1])
							for rec in records])
		for k in sorted_results:
			genome_id = k[0]
			current_results = k[1]
			# Get match that was EXC or INF
			current_match = [c for c in current_results[1:] if c[3] in ['EXC', 'INF']][0]
			cds_hash = current_match[2]
			# Do not add to schema CDSs that matched several loci
			if cds_hash not in repeated:
				if cds_hash in matched_alleles:
					locus_results[genome_id].append(matched_alleles[cds_hash])
				else:
					max_alleleid += 1
					current_alleleid_str = f'{alleleid_prefix}{max_alleleid}'
					locus_results[genome_id].append('INF-{0}'.format(current_alleleid_str))
					matched_alleles[cds_hash] = current_alleleid_str
					novel_alleles.append([cds_hash, current_alleleid_str])
					total_novel += 1
					# EXC to INF to enable accurate count of INF classifications
					# Some INF classifications might be converted to NIPH based on similar
					# matches on the same genome. Matches to the new INF/NIPH will be
					# classified as EXC and need to be added as new alleles and converted to INF
					if current_results[0] == 'EXC':
						locus_results[genome_id][0] = 'INF'
			# Classify as PAMA when a CDS matches multiple loci
			else:
				locus_results[genome_id][0] = ct.ALLELECALL_CLASSIFICATIONS[9]

		# Save updated results
		fo.pickle_dumper(locus_results, locus_files[1])

	# Save novel alleles data
	if len(novel_alleles) > 0:
		novel_outfile = fo.join_paths(output_directory, [f'{locus_id}'])
		fo.pickle_dumper(novel_alleles, novel_outfile)
		return [locus_files[0], novel_outfile, total_novel]


def create_novel_fastas(inferred_alleles, inferred_representatives,
						sequences_file, output_directory, loci_finder):
	"""Create FASTA files with the novel alleles for each locus.

	Parameters
	----------
	inferred_alleles : dict
		Dictionary with paths to loci FASTA files as keys and
		lists with SHA-256 hashes, allele integer identifiers and
		sequence identifiers for each novel allele.
	inferred_representatives : dict
		Dictionary with loci identifiers as keys and lists with
		sequence identifiers, SHA-256 hashes and allele integer
		identifiers for each novel representative allele.
	sequences_file : str
		Path to FASTA file that contains the distinct coding
		sequences identified in the inputs.
	output_directory : str
		Path to the directory where the FASTA files that contain
		the novel alleles will be saved to.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	total_inferred : int
		Total number of inferred alleles added to the schema.
	total_representatives : int
		Total number of representative alleles added to the
		schema.
	updated_novel : dict
		Dictionary with loci identifiers as keys and paths to
		the FASTA files that contain the novel alleles as values.
	"""
	# Create index for Fasta file with distinct CDSs
	sequence_index = fao.index_fasta(sequences_file)

	# Store paths to FASTA files with new alleles
	updated_novel = {}
	# Count number of novel and representative alleles added to schema
	total_inferred = 0
	total_representatives = 0
	for locus_novel in inferred_alleles:
		locus_id = loci_finder.search(locus_novel[1]).group()
		current_novel = fo.pickle_loader(locus_novel[1])

		updated_novel[locus_novel[0]] = []
		# Get novel alleles through indexed Fasta file
		novel_alleles = ['>{0}_{1}\n{2}'.format(locus_id, a[1],
												str(sequence_index.get(a[2]).seq))
						 for a in current_novel]
		# Create Fasta file with novel alleles
		novel_file = fo.join_paths(output_directory, ['{0}.fasta'.format(locus_id)])
		fo.write_lines(novel_alleles, novel_file)
		updated_novel[locus_novel[0]].append(novel_file)
		total_inferred += len(novel_alleles)

		# Add representatives
		novel_representatives = inferred_representatives.get(locus_id, None)
		if novel_representatives is not None:
			reps_sequences = ['>{0}_{1}\n{2}'.format(locus_id, a[2], str(sequence_index.get(a[0]).seq))
							  for a in novel_representatives]
			# Create Fasta file with novel representative alleles
			novel_rep_file = fo.join_paths(output_directory, ['short', '{0}_short.fasta'.format(locus_id)])
			fo.write_lines(reps_sequences, novel_rep_file)
			updated_novel[locus_novel[0]].append(novel_rep_file)
			total_representatives += len(reps_sequences)

	return [total_inferred, total_representatives, updated_novel]


def add_inferred_alleles(inferred_alleles, loci_finder):
	"""Add inferred alleles to a schema.

	Prameters
	---------
	inferred_alleles . dict
		Dictionary with paths to loci FASTA files in the schema as
		keys and paths to the FASTA files that contain the novel
		alleles as values.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	alleles_added : list
		List that contains the number of novel alleles added for each
		locus (same order as `inferred_alleles` keys).
	"""
	alleles_added = {}
	for locus, files in inferred_alleles.items():
		locus_id = loci_finder.search(locus).group()
		# Read novel alleles
		novel_alleles = fo.read_lines(files[0])
		# Append novel alleles to locus FASTA file in the schema
		fo.write_lines(novel_alleles, locus, write_mode='a')
		alleles_added.setdefault(locus_id, []).append(len(novel_alleles)//2)

		# Add representative alleles to FASTA files in the 'short' directory
		if len(files) > 1:
			novel_representatives = fo.read_lines(files[1])
			locus_short_path = fo.join_paths(os.path.dirname(locus),
											 ['short', locus_id+'_short.fasta'])
			fo.write_lines(novel_representatives, locus_short_path, write_mode='a')
			alleles_added[locus_id].append(len(novel_representatives)//2)
		else:
			alleles_added[locus_id].append(0)

	return alleles_added


def select_highest_scores(blast_outfile):
	"""Select the highest-scoring match for each BLAST target.

	Parameters
	----------
	blast_outfile : str
		Path to the TSV file created by BLAST.

	Returns
	-------
	best_matches : list
		List with the highest-scoring match/line for each
		distinct target.
	"""
	blast_results = fo.read_tabular(blast_outfile)
	# Sort results based on decreasing raw score
	blast_results = im.sort_iterable(blast_results,
									 lambda x: int(x[5]),
									 reverse=True)

	# Select matches with highest score for each target
	best_matches = {}
	for result in blast_results:
		# Only get the best raw score for each target
		if result[4] not in best_matches:
			best_matches[result[4]] = result

	return best_matches


def process_blast_results(blast_results, bsr_threshold, query_scores):
	"""Process BLAST results to get the data relevant for classification.

	Parameters
	----------
	blast_results : list
		List with one sublist per BLAST match (must have one
		sublist per target with the highest-scoring match).
	bsr_threshold : float
		BLAST Score Ratio (BSR) value to select matches. Matches
		will be kept if the computed BSR is equal or greater than
		this value.
	query_scores :  dict
		Dictionary with loci representative sequence identifiers
		as values and a tuple with sequence length and the raw
		score for the self-alignment as value.

	Returns
	-------
	match_info : dict
		Dictionary with the distinct target sequence identifiers
		as keys and a tuple with the BSR value, target sequence
		length, query sequence length and query sequence identifier
		for the highest-scoring match for each target as values.
	"""
	match_info = {}
	for r in blast_results:
		query_id = r[0]
		target_id = r[4]
		raw_score = float(r[6])
		# Compute BSR
		# Might fail if it is not possible to get the query self-score
		try:
			bsr = cf.compute_bsr(raw_score, query_scores[query_id][1])
		except Exception as e:
			print('Could not get the self-score for the representative '
				  f'allele {query_id}', e)
			continue
		# Only keep matches above BSR threshold
		if bsr >= bsr_threshold:
			# BLAST has 1-based positions
			qstart = (int(r[1])-1)*3  # Subtract 3 to exclude start position
			qend = (int(r[2])*3)+3  # Add 3 to count stop codon
			target_length = (int(r[5])*3)+3
			query_length = query_scores[query_id][0]

			match_info[target_id] = (bsr, qstart, qend,
									 target_length, query_length, query_id)

	return match_info


def expand_matches(match_info, pfasta_index, dfasta_index, dhashtable,
				   phashtable, inv_map, close_to_tip):
	"""Expand distinct matches to create matches for all matched inputs.

	Parameters
	----------
	match_info : dict
		Dictionary with the distinct target sequence identifiers
		as keys and a tuple with the BSR value, target sequence
		length, query sequence length and query sequence identifier
		for the highest-scoring match for each target as values.
	pfasta_index : Bio.File._IndexedSeqFileDict
		Fasta file index created with BioPython. Index for the
		distinct protein sequences.
	dfasta_index : Bio.File._IndexedSeqFileDict
		Fasta file index created with BioPython. Index for the
		distinct DNA coding sequences.
	dhashtable : dict
		Dictionary with SHA-256 hashes for distinct DNA
		sequences extracted from the inputs and lists of
		genome integer identifiers enconded with the
		polyline algorithm as values.
	phashtable : dict
		Dictionary with SHA-256 hashes for distinct protein
		sequences extracted from the inputs and lists of
		sequence identifiers enconded with the polyline
		algorithm as values.
	inv_map : dict
		Dictionary with input integer identifiers as keys
		and input string identifiers as values.
	close_to_tip : dict
		Dictionary with input identifiers as keys and dictionaries with
		the SHA-256 hashes of the CDSs that are close to contig tips as
		keys and tuples with the CDS coordinates as values.

	Returns
	-------
	input_matches : dict
		Dictionary with input integer identifiers as keys
		and tuples with information about matches identified
		in the inputs as values.
	"""
	input_matches = {}
	# Expand each distinct protid match into the distinct CDS matches
	for target_id in match_info:
		# Get distinct protein sequence and expand into distinct CDS seqids
		target_protein = str(pfasta_index.get(target_id).seq)
		target_phash = im.hash_sequence(target_protein)
		target_integers = im.polyline_decoding(phashtable[target_phash])
		target_seqids = ['{0}-protein{1}'.format(inv_map[target_integers[i+1]], target_integers[i])
						 for i in range(0, len(target_integers), 2)]
		# Associate match data to each genome that contains one of the distinct CDSs
		for seqid in target_seqids:
			target_cds = str(dfasta_index.get(seqid).seq)
			target_dhash = im.hash_sequence(target_cds)
			# Get ids for all genomes with same CDS as representative
			target_inputs = im.polyline_decoding(dhashtable[target_dhash])[1:]
			for i in target_inputs:
				# There will be no data about CDSs close to contig tips if input FASTA contain CDSs
				genome_coordinates = close_to_tip.get(inv_map[i], {}).get(target_dhash, [()])[0]
				input_matches.setdefault(i, []).append((target_id, target_phash,
														target_dhash, *match_info[target_id],
														genome_coordinates))

	return input_matches


def identify_paralogous(repeated, output_directory):
	"""Identifiy paralogous loci based on CDSs that matched multiple loci.

	Parameters
	----------
	repeated : dict
		Dictionary with hashes for the CDSs that matched multiple loci
		as keys and a list with information about the genome of origin
		and matched loci as values.
	output_directory : str
		Path to the output directory where the file with
		the list of paralogus loci will be created.

	Returns
	-------
	The total number of paralogous loci detected.
	The path to the file that contains the list of paralogous loci.
	The path to the file that contains paralogous matches per input.
	"""
	paralogous_data = {}
	paralogous_counts = {}
	for value in repeated.values():
		for p in value:
			# CDS identifier as key
			paralogous_data.setdefault(p[2], {})
			# Genome identifier and locus identifier
			paralogous_data[p[2]].setdefault(p[0], []).append(p[1])
			paralogous_counts.setdefault(p[1], []).append(1)

	# Write paralogous loci counts
	paralogous_counts_lines = [ct.PARALOGOUS_COUNTS_HEADER]
	paralogous_counts_lines.extend(['{0}\t{1}'.format(k, sum(v))
									for k, v in paralogous_counts.items()])
	paralogous_counts_outfile = fo.join_paths(output_directory,
											  [ct.PARALOGOUS_COUNTS_BASENAME])
	fo.write_lines(paralogous_counts_lines, paralogous_counts_outfile)

	# Write groups of paralogous loci per input
	paralogous_lines = [ct.PARALOGOUS_LIST_HEADER]
	for cds, g in paralogous_data.items():
		paralogous_lines.extend(['{0}\t{1}\t{2}'.format(k, '|'.join(v), cds)
								 for k, v in g.items()])

	paralogous_loci_outfile = fo.join_paths(output_directory,
											[ct.PARALOGOUS_LOCI_BASENAME])
	fo.write_lines(paralogous_lines, paralogous_loci_outfile)

	return [len(paralogous_counts), paralogous_counts_outfile, paralogous_loci_outfile]


def classify_inexact_matches(locus, genomes_matches, inv_map,
							 locus_results_file, locus_mode, temp_directory,
							 size_threshold, blast_score_ratio, output_directory,
							 cds_input):
	"""Classify inexact matches to a locus.

	Parameters
	----------
	locus : str
		Locus identifier.
	genomes_matches : str
		Path to file with data about matches found in the inputs.
	inv_map : dict
		Dictionary with input integer identifiers as keys and
		input string identifiers as values.
	locus_results_file : str
		Path to file with classification results for the locus.
	locus_mode : list
		List where wthe first element is the locus allele size mode
		and the second element is a list with the length values for
		all alleles.
	temp_directory : str
		Path to the directory where temporary files will be stored.
	size_threshold : float or None
		Sequence size variation threshold.
	blast_score_ratio : float
		BLAST Score Ratio value.
	output_directory : str
		Path to the directory where the files with information for
		each locus will be saved to.
	cds_input : bool
		True if the input files contained CDSs, False if they contained
		contigs.

	Returns
	-------
	locus_info_file : str
		Path to pickle file with the data about the classification
		of inexact matches (contains dictionary with locus identifier
		as key and a list with the path to the pickle file with the
		locus classifications, locus allele size mode, sequence
		identifiers of the distinct sequences that were classified
		and a list with data about representative candidates as value).
	"""
	# Import locus classification data to update
	locus_results = fo.pickle_loader(locus_results_file)

	# Import locus match data to classify
	genomes_matches = fo.pickle_loader(genomes_matches)

	# Initialize lists to store hashes of CDSs that have been classified
	seen_dna = {}
	seen_prot = set()
	# Initialize list to store sequence identifiers that have been classified
	# and that should be excluded from next steps
	excluded = []
	representative_candidates = []
	for genome, matches in genomes_matches.items():
		for match in matches:
			# Get sequence identifier for the distinct protid that matched
			target_seqid = match[0]
			# Get allele identifier for the schema or novel representative
			rep_alleleid = match[8]
			# The representative identifier can be the integer ID associated to
			# the allele if it is in the schema or the full CDS identifier if the
			# representative allele was selected in a previous iteration and is
			# not in the schema
			# Representative is in the schema and has allele identifier
			try:
				# Split to get allele identifier
				# Need to replace '*' for novel alleles added to schemas from Chewie-NS
				int(rep_alleleid.replace('*', '').split('_')[-1])
				rep_alleleid = rep_alleleid.split('_')[-1]
			# Representative allele was selected in a previous iteration and is not in the schema
			except Exception as e:
				pass

			# Get hash of the CDS DNA sequence
			target_dna_hash = match[2]
			# Get hash of the translated CDS sequence
			target_prot_hash = match[1]
			# Get the BSR value
			bsr = match[3]

			# CDS DNA sequence was identified in one of the previous inputs
			# This might change classification to NIPH if the input
			# already had a classification for the current locus
			if target_dna_hash in seen_dna:
				locus_results = update_classification(genome, locus_results,
													  (seen_dna[target_dna_hash], target_seqid,
													   target_dna_hash, 'EXC', 1.0))
				continue

			# Translated CDS matches other translated CDS that was classified
			if target_prot_hash in seen_prot:
				locus_results = update_classification(genome, locus_results,
													  (rep_alleleid, target_seqid,
													   target_dna_hash, 'INF', 1.0))
				# Add DNA hash to classify the next match as EXC
				seen_dna[target_dna_hash] = target_seqid
				continue

			if cds_input is False and len(match[-1]) > 0:
				# Get representative length
				representative_length = match[7]
				# Get target left and right positions that aligned
				representative_leftmost_pos = match[4]
				representative_rightmost_pos = match[5]
				# Cannot be a PLOT3/PLOT5/LOTSC if representative fully aligns
				if representative_leftmost_pos > 0 or representative_rightmost_pos < representative_length:
					# Determine if it is PLOT3, PLOT5 or LOTSC
					relative_pos = contig_position_classification(representative_length,
																  representative_leftmost_pos,
																  representative_rightmost_pos,
																  *match[-1])

					if relative_pos is not None:
						locus_results = update_classification(genome, locus_results,
																(rep_alleleid, target_seqid,
																target_dna_hash, relative_pos, bsr))
						# Exclude CDSs classified as PLOT3/5
						excluded.append(target_seqid)
						continue

			target_dna_len = match[6]
			# Check if ASM or ALM
			relative_size = allele_size_classification(target_dna_len, locus_mode, size_threshold)
			if relative_size is not None:
				locus_results = update_classification(genome, locus_results,
													  (rep_alleleid, target_seqid,
													   target_dna_hash, relative_size, bsr))
				# Exclude CDSs classified as ASM/ALM
				excluded.append(target_seqid)
				continue

			# Add INF
			# This will turn into NIPH if there are multiple hits for the same input
			locus_results = update_classification(genome, locus_results,
												  (rep_alleleid, target_seqid,
												   target_dna_hash, 'INF', bsr))

			seen_dna[target_dna_hash] = target_seqid
			excluded.append(target_seqid)
			seen_prot.add(target_prot_hash)

		# Update locus mode value if classification for genome is INF
		if genome in locus_results and locus_results[genome][0] == 'INF':
			# Append length of inferred allele to list with allele sizes
			locus_mode[1].append(target_dna_len)
			# Compute mode
			locus_mode[0] = sm.determine_mode(locus_mode[1])[0]
			# Only add as representative candidate if classification is not NIPH
			inf_bsr = locus_results[genome][1][4]
			if inf_bsr >= blast_score_ratio and inf_bsr < blast_score_ratio+0.1:
				representative_candidates.append((genome, target_seqid,
												  match[8], target_dna_hash))

	# Save updated results
	fo.pickle_dumper(locus_results, locus_results_file)

	# Save info about updated mode, excluded ids and representative candidates
	locus_info = {locus: [locus_results_file, locus_mode,
						  excluded, representative_candidates]}
	locus_info_file = fo.join_paths(output_directory, ['{0}_classification_info'.format(locus)])
	fo.pickle_dumper(locus_info, locus_info_file)

	return locus_info_file


def create_missing_fasta(class_files, fasta_file, input_map, dna_hashtable,
						 output_directory, classification_labels, cds_input,
						 loci_finder):
	"""Create Fasta file with sequences for missing data classes.

	Parameters
	----------
	class_files : dict
		Dictionary with paths to loci files as keys and paths to
		pickled files with classification results as values.
	fasta_file : str
		Path to Fasta file with the distinct CDS extracted from
		the input genomes.
	input_map : dict
		Dictionary with the mapping between the input integer
		identifiers and input string identifiers.
	dna_hashtable : dict
		Dictionary with hashes of the distinct CDS extracted from
		input genomes as keys and lists containing the integer
		identifiers for te inputs that contained the CDS encoded
		with the polyline algorithm.
	output_directory : str
		Path to the output directory where the Fasta file will
		be saved to.
	classification_labels : list
		List with the class labels attributed by chewBBACA.
	cds_input : bool
		False if there are files with CDS coordinates, True otherwise.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	output_fasta_file : str
		Path to the output FASTA file.
	output_tsv_file : str
		Path to the output TSV file.
	"""
	invalid_cases = classification_labels[2:-1]

	# Index FASTA file to get CDSs
	dna_index = fao.index_fasta(fasta_file)

	# Define path to output FASTA file
	output_fasta_file = fo.join_paths(output_directory, [ct.MISSING_FASTA_BASENAME])

	# Define path to output TSV file
	output_tsv_file = fo.join_paths(output_directory, [ct.MISSING_TSV_BASENAME])
	# Write header line
	fo.write_lines([ct.MISSING_HEADER], output_tsv_file, write_mode='a')

	# Add integer index to FASTA header and first TSV column
	index = 1
	# Get hash if coordinates are available, seqid otherwise
	id_index = 2 if cds_input is False else 1
	for locus, file in class_files.items():
		locus_lines = []
		locus_records = []
		locus_id = loci_finder.search(locus).group()
		locus_classifications = fo.pickle_loader(file)
		# Get data for genomes that do not have EXC or INF classifications
		# it will not get invalid classes if a genome is classified as EXC|INF
		for gid, v in locus_classifications.items():
			genome_id = input_map[gid]
			if v[0] in invalid_cases:
				matches = [[e[id_index], e[3]] for e in v[1:]]
				for match in matches:
					# Get hash or seqid depending on input type
					match_id = match[0]
					match_classification = match[1]
					# Get seqid of the representative CDS
					if id_index == 2:
						match_protid, match_gid = im.polyline_decoding(dna_hashtable[match_id])[:2]
						match_id = f'{input_map[match_gid]}-protein{match_protid}'

					sequence_header = (f'{index}|{genome_id}|{locus_id}&{v[0]}|{match_id}&{match_classification}')
					sequence = str(dna_index[match_id].seq)
					record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [sequence_header, sequence])
					locus_records.append(record)
					locus_lines.append(f'{index}\t{genome_id}\t{locus_id}\t{v[0]}\t{match_id}\t{match_classification}')
					index += 1

		if len(locus_records) > 0:
			fo.write_lines(locus_records, output_fasta_file, write_mode='a')
			fo.write_lines(locus_lines, output_tsv_file, write_mode='a')

	return [output_fasta_file, output_tsv_file]


def select_representatives(representative_candidates, locus, fasta_file,
						   iteration, output_directory, blastp_path,
						   blast_db, blast_score_ratio, threads,
						   blastdb_aliastool_path):
	"""Select new representative alleles for a locus.

	Parameters
	----------
	representative_candidates : dict
		Dictionary with sequence identifiers as keys and sequence
		hashes as values.
	locus : str
		Locus identifier.
	fasta_file : path
		Path to Fasta file that contains the translated sequences
		of the representative candidates.
	iteration : int
		Iteration number to add to generated files.
	output_directory : str
		Path to the output directory.
	blastp_path : str
		Path to the BLASTp executable.
	blast_db : str
		Path to the BLAST database.
	blast_score_ratio : float
		BLAST Score Ratio value.
	threads : int
		Number of threads passed to BLAST.
	blastdb_aliastool_path : str or None
		Path to the `blastdb_aliastool` executable used to convert
		the list of seqids to binary format.

	Returns
	-------
	locus : str
		Locus identifier.
	selected : list
		List that contains one tuple per selected candidate (tuples
		contain the sequence identifier and the sequence hash for
		each new representative).
	"""
	# Create file with candidate seqids
	ids_file = fo.join_paths(output_directory,
							 ['{0}_candidates_ids_{1}.fasta'.format(locus, iteration)])
	fo.write_lines(list(representative_candidates.keys()), ids_file)
	# Convert to binary format if BLAST>=2.10
	binary_file = f'{ids_file}.bin'
	blastp_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
											ids_file,
											binary_file)
	ids_file = binary_file

	# BLASTp to compare all candidates
	blast_output = fo.join_paths(output_directory,
								 ['{0}_candidates_{1}_blastout.tsv'.format(locus, iteration)])
	# Define -max_target_seqs to reduce execution time
	blastp_std = bw.run_blast(blastp_path, blast_db, fasta_file,
							  blast_output, threads=threads,
							  ids_file=ids_file, max_targets=100)

	blast_results = fo.read_tabular(blast_output)
	# Get self-score for all candidates
	candidates_self_scores = {line[0]: ((int(line[3])*3)+3, float(line[6]))
							  for line in blast_results if line[0] == line[4]}
	# Select results between different candidates
	blast_results = [line for line in blast_results if line[0] != line[4]]

	# Compute BSR
	for line in blast_results:
		line.append(cf.compute_bsr(float(line[6]), candidates_self_scores[line[0]][1]))
	# Sort by sequence length to process longest candidates first
	blast_results = sorted(blast_results, key=lambda x: int(x[3]))

	excluded_candidates = set()
	for r in blast_results:
		if r[7] >= blast_score_ratio+0.1 and r[0] not in excluded_candidates:
			excluded_candidates.add(r[4])

	selected_candidates = set(representative_candidates.keys()) - excluded_candidates

	selected = [(seqid, representative_candidates[seqid]) for seqid in selected_candidates]

	return [locus, selected]


def count_invalid(input_ids, invalid_seqids, cds_index, distinct_htable):
	"""Count the number of invalid CDSs per input.

	Parameters
	----------
	input_ids : list
		Dicionaries with input unique identifiers mapped to
		input integer identifiers and vice versa.
	invalid_seqids : list
		List of seqids for the invalid CDSs.
	cds_index : Bio.File._IndexedSeqFileDict
		Biopython index for the Fasta file that contains the
		CDSs.
	distinct_htable : dict
		Dictionary with SHA-256 hashes for distinct DNA sequences
		extracted from the inputs and lists of genome integer
		identifiers enconded with the polyline algorithm as values.
	"""
	# Count number of invalid CDSs per input
	invalid_counts = {k: 0 for k in input_ids[0]}
	for seqid in invalid_seqids:
		current_seq = str(cds_index.get(seqid).seq)
		seq_hash = im.hash_sequence(current_seq)
		# Get list of inputs that contain CDS
		genome_list = im.polyline_decoding(distinct_htable[seq_hash])[1:]
		# Increment count for all inputs that contain the CDS
		for gid in genome_list:
			invalid_counts[input_ids[1][gid]] += 1

	return invalid_counts


def merge_blast_results(blast_outfiles, output_directory, loci_finder):
	"""Concatenate BLAST output files based on locus identifier.

	Parameters
	----------
	blast_outfiles : list
		List with paths to BLAST output files. Paths must contain
		the locus identifier.
	output_directory : str
		Path to the output directory.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	concatenated_files : list
		List with paths to the concatenated files.
	"""
	# Group BLAST results by locus
	loci_files = {}
	for file in blast_outfiles:
		# Filenames include allele identifier
		locus_id = loci_finder.search(file).group()
		loci_files.setdefault(locus_id, []).append(file)

	# Concatenate files with results for the same locus
	concatenated_files = []
	for locus, files in loci_files.items():
		outfile = fo.join_paths(output_directory, [f'{locus}_concatenated_blastout.tsv'])
		fo.concatenate_files(files, outfile)
		concatenated_files.append(outfile)
		
	return concatenated_files


def allele_calling(fasta_files, schema_directory, temp_directory,
				   loci_modes, loci_files, config, pre_computed_dir,
				   loci_finder):
	"""Perform allele calling for a set of inputs.

	Parameters
	----------
	fasta_files : list
		List with the full paths to the input FASTA files, one per
		strain.
	schema_directory : str
		Path to the schema directory.
	temp_directory : str
		Path to the temporary directory in which the intermediate files
		are created.
	loci_modes : dict
		Dictionary with the allele size mode for each locus.
	loci_files : dict
		Dictionary with the paths to the loci FASTA files as keys and
		the loci integer identifiers as values.
	config : dict
		Dictionary with the config values that will be used to perform
		allele calling.
	pre_computed_dir : str
		Path to the the directory that contains the files with the
		pre-computed hash tables.
	loci_finder : re.Pattern
		Regular expression object to search for loci identifiers
		in paths and filenames.

	Returns
	-------
	template_dict : dict
		Dictionary used to store data and paths to files with results
		(paths to files with the classification results, dictionary with
		the mapping between input identifiers and unique integer identifiers,
		paths to files with the CDS coordinates for each input, path to a
		FASTA file with all distinct DNA sequences, path to a FASTA file
		with all distinct protein sequences, a hash table with the mapping
		between the SHA256 hash for each distinct CDS extracted from the
		input genomes and the list of input genomes that contain each CDS,
		a hash table with the mapping between the SHA256 hash for each
		translated distinct CDS and the distinct DNA CDSs that code for each
		protein, list with the paths to input files for which it was not
		possible to predict CDSs with Pyrodigal, path to a file with the list
		of invalid CDSs, list with the sequence identifiers of the distinct
		CDSs that were not classified, dictionary with the mapping between the
		sequence identifiers of the representative alleles and their
		self-alignment BLASTp raw score, dictionary with information about
		new representatives for each locus).
	"""
	# Get dictionary template to store variables to return
	template_dict = ct.ALLELECALL_DICT

	# Map input file paths to unique identifier (prefix before first '.' in basename)
	full_to_basename = im.mapping_function(fasta_files, fo.file_basename, [False])
	full_to_unique = {k: fo.split_joiner(v, [0], '.') for k, v in full_to_basename.items()}

	# Create directory to store files with Pyrodigal results
	pyrodigal_path = fo.join_paths(temp_directory, ['1_cds_prediction'])
	fo.create_directory(pyrodigal_path)

	# Inputs are genome assemblies
	if config['CDS input'] is False:
		# Run Pyrodigal to determine CDSs for all input genomes
		print(f'\n {ct.CDS_PREDICTION} ')
		print('='*(len(ct.CDS_PREDICTION)+2))

		# Gene prediction step
		print(f'Predicting CDSs for {len(fasta_files)} inputs...')
		pyrodigal_results = cf.predict_genes(full_to_unique,
											 config['Prodigal training file'],
											 config['Translation table'],
											 config['Prodigal mode'],
											 config['CPU cores'],
											 pyrodigal_path)

		# Dictionary with info about inputs for which gene prediction failed
		# Total number of CDSs identified in the inputs
		# Paths to FASTA files with the extracted CDSs
		# Paths to files with the coordinates of the CDSs extracted for each input
		# Total number of CDSs identified per input
		# Dictionary with info about the CDSs closer to contig tips per input
		failed, total_extracted, cds_fastas, cds_coordinates, cds_counts, close_to_tip = pyrodigal_results

		if len(failed) > 0:
			print(f'\nFailed to predict CDSs for {len(failed)} inputs.')
			print('Make sure that Pyrodigal runs in meta mode (--pm meta) '
				  'if any input file has less than 100kbp.')
		if len(cds_fastas) == 0:
			sys.exit(f'\n{ct.CANNOT_PREDICT}')

		print(f'\nExtracted a total of {total_extracted} CDSs from {len(fasta_files)} inputs.')
	# Inputs are Fasta files with the predicted CDSs
	else:
		# Rename the CDSs in each file based on the input unique identifiers
		print(f'\nRenaming CDSs for {len(full_to_unique)} input files...')

		renaming_inputs = []
		cds_fastas = []
		for k, v in full_to_unique.items():
			output_file = fo.join_paths(pyrodigal_path, [f'{v}.fasta'])
			cds_prefix = f'{v}-protein'
			renaming_inputs.append([k, output_file, 1, 50000,
									cds_prefix, False, fao.integer_headers])
			cds_fastas.append(output_file)

		# Rename CDSs in files
		renaming_results = mo.map_async_parallelizer(renaming_inputs,
													 mo.function_helper,
													 config['CPU cores'],
													 show_progress=False)

		# No inputs failed gene prediction
		failed = []
		# Cannot get CDS coordinates if skipping gene prediction
		cds_coordinates = None
		close_to_tip = {}
		cds_counts = {r[0]: r[1] for r in renaming_results}
		cds_counts = {full_to_unique[k]: v for k, v in cds_counts.items()}
		total_cdss = sum([r[1] for r in renaming_results])
		print(f'Input files contain a total of {total_cdss} coding sequences.')

	if len(failed) > 0:
		# Exclude inputs that failed gene prediction
		full_to_unique = im.prune_dictionary(full_to_unique, failed.keys())
		# Write Prodigal stderr for inputs that failed gene prediction
		failed_lines = [f'{k}\t{v}' for k, v in failed.items()]
		failed_outfile = fo.join_paths(os.path.dirname(temp_directory),
									   ['gene_prediction_failures.tsv'])
		fo.write_lines(failed_lines, failed_outfile)

	template_dict['cds_coordinates'] = cds_coordinates

	# Map input identifiers to integers
	# Use the mapped integers to refer to each input
	# This reduces memory usage compared to using string identifiers
	unique_to_int = im.integer_mapping(full_to_unique.values())
	int_to_unique = im.invert_dictionary(unique_to_int)
	template_dict['int_to_unique'] = int_to_unique

	# Change to unique integer identifiers
	cds_counts = {unique_to_int[k]: v for k, v in cds_counts.items()}
	template_dict['cds_counts'] = cds_counts

	# Concatenate subgroups of FASTA files before deduplication
	num_chunks = 20 if config['CPU cores'] <= 20 else config['CPU cores']
	concatenation_inputs = im.divide_list_into_n_chunks(cds_fastas, num_chunks)
	file_index = 1
	cds_files = []
	for group in concatenation_inputs:
		output_file = fo.join_paths(pyrodigal_path,
									['cds_{0}.fasta'.format(file_index)])
		fo.concatenate_files(group, output_file)
		cds_files.append(output_file)
		file_index += 1

	# Create directory to store files from pre-process steps
	preprocess_dir = fo.join_paths(temp_directory, ['2_cds_preprocess'])
	fo.create_directory(preprocess_dir)

	# DNA sequences deduplication step
	# keep hash of unique sequences and a list with the integer
	# identifiers of genomes that have those sequences
	# lists of integers are encoded with polyline algorithm
	print(f'\n {ct.CDS_DEDUPLICATION} ')
	print('='*(len(ct.CDS_DEDUPLICATION)+2))
	# Create directory to store files from DNA deduplication
	dna_dedup_dir = fo.join_paths(preprocess_dir, ['cds_deduplication'])
	fo.create_directory(dna_dedup_dir)
	print('Identifying distinct CDSs...')
	dna_dedup_results = cf.exclude_duplicates(cds_files, dna_dedup_dir,
											  config['CPU cores'],
											  [unique_to_int, int_to_unique])

	dna_distinct_htable, distinct_seqids, distinct_file, repeated = dna_dedup_results
	print(f'Identified {len(dna_distinct_htable)} distinct CDSs.')
	template_dict['dna_fasta'] = distinct_file
	template_dict['dna_hashtable'] = dna_distinct_htable

	# Delete concatenated FASTA files
	fo.remove_files(cds_files)

	# Get mapping between locus file path and locus identifier
	loci_basenames = im.mapping_function(loci_files, fo.file_basename, [False])

	# Create directory to store files with classification results
	classification_dir = fo.join_paths(temp_directory, ['classification_files'])
	fo.create_directory(classification_dir)
	# Create files with empty dict to store results per locus
	empty_results = {}
	classification_files = {}
	for file in loci_files:
		file_path = fo.join_paths(classification_dir, [f'{loci_basenames[file]}_results'])
		fo.pickle_dumper(empty_results, file_path)
		classification_files[file] = file_path

	# CDS exact matching
	print(f'\n {ct.CDS_EXACT} ')
	print('='*(len(ct.CDS_EXACT)+2))

	matched_seqids = []
	total_matches = 0
	distinct_seqids_count = 0
	print('Searching for CDS exact matches...')
	# List files in pre-computed directory
	dna_tables = fo.listdir_fullpath(pre_computed_dir, 'DNAtable')
	for file in dna_tables:
		dna_matches = dna_exact_matches(file,
										dna_distinct_htable,
										im.invert_dictionary(loci_files),
										classification_files,
										int_to_unique)
		# Might have duplicate IDs
		matched_seqids.extend(dna_matches[0])
		total_matches += dna_matches[1]
		distinct_seqids_count += dna_matches[2]

	print(f'Found {total_matches} exact matches ({distinct_seqids_count} distinct schema alleles).')

	# Save list of seqids that matched
	dna_exact_outdir = fo.join_paths(preprocess_dir, ['cds_exact_matching'])
	fo.create_directory(dna_exact_outdir)
	dna_exact_outfile = fo.join_paths(dna_exact_outdir, ['cds_exact_matches.txt'])
	fo.write_lines(matched_seqids, dna_exact_outfile)

	# Get seqids for unclassified sequences
	# Exclude seqids that matched from distinct seqids
	unclassified_seqids = im.filter_list(distinct_seqids, matched_seqids)

	print(f'Unclassified CDSs: {len(unclassified_seqids)}')

	# User only wants to determine exact matches or all sequences were classified
	if config['Mode'] == 1 or len(unclassified_seqids) == 0:
		template_dict['classification_files'] = classification_files
		template_dict['protein_fasta'] = distinct_file
		template_dict['unclassified_ids'] = unclassified_seqids
		template_dict['representatives'] = {}

		return template_dict

	# Index FASTA file with distinct DNA sequences
	dna_index = fao.index_fasta(distinct_file)

	# Translate CDSs
	print(f'\n {ct.CDS_TRANSLATION} ')
	print('='*(len(ct.CDS_TRANSLATION)+2))

	# Create directory to store translation results
	cds_translation_dir = fo.join_paths(preprocess_dir, ['cds_translation'])
	fo.create_directory(cds_translation_dir)
	print(f'Translating {len(unclassified_seqids)} CDSs...')
	# This step excludes small sequences
	ts_results = cf.translate_sequences(unclassified_seqids, distinct_file,
										cds_translation_dir,
										config['Translation table'],
										config['Minimum sequence length'],
										config['CPU cores'])

	protein_file, ut_seqids, ut_lines = ts_results
	print(f'\n{len(ut_seqids)} CDSs could not be translated.')

	# Create file with list of invalid CDSs
	invalid_file = fo.join_paths(temp_directory, [ct.INVALID_CDS_BASENAME])
	invalid_alleles = im.join_list(im.sort_iterable(ut_lines), '\n')
	fo.write_to_file(invalid_alleles, invalid_file, 'w', '\n')
	# Count the number of invalid CDSs per input
	invalid_counts = count_invalid([unique_to_int, int_to_unique], ut_seqids,
								   dna_index, dna_distinct_htable)
	template_dict['invalid_alleles'] = [invalid_file, invalid_counts]
	print(f'Unclassified CDSs: {len(unclassified_seqids)-len(ut_seqids)}')

	# All sequences were excluded during translation
	if len(unclassified_seqids)-len(ut_seqids) == 0:
		template_dict['classification_files'] = classification_files
		template_dict['protein_fasta'] = distinct_file
		template_dict['unclassified_ids'] = unclassified_seqids
		template_dict['representatives'] = {}

		return template_dict

	# Protein deduplication step
	print(f'\n {ct.PROTEIN_DEDUPLICATION} ')
	print('='*(len(ct.PROTEIN_DEDUPLICATION)+2))

	# Create directory to store files from protein deduplication
	protein_dedup_dir = fo.join_paths(preprocess_dir, ['protein_deduplication'])
	fo.create_directory(protein_dedup_dir)
	print('Identifying distinct proteins...')
	ds_results = cf.exclude_duplicates([protein_file], protein_dedup_dir, 1,
									   [unique_to_int, int_to_unique], True)

	distinct_pseqids, representative_pseqids, representative_pfasta, _ = ds_results

	print(f'Identified {len(distinct_pseqids)} distinct proteins.')
	template_dict['protein_hashtable'] = distinct_pseqids

	# Identify exact matches at protein level
	# Exact matches at protein level are novel alleles
	print(f'\n {ct.PROTEIN_EXACT} ')
	print('='*(len(ct.PROTEIN_EXACT)+2))
	total_matches = 0
	distinct_seqids_count = 0
	distinct_protids = 0
	matched_seqids = []
	previous_hashes = set()
	print('Searching for Protein exact matches...')
	# List files in pre-computed dir
	protein_tables = fo.listdir_fullpath(pre_computed_dir, 'PROTEINtable')
	for file in protein_tables:
		protein_matches = protein_exact_matches(file,
												distinct_pseqids,
												im.invert_dictionary(loci_files),
												classification_files,
												int_to_unique,
												dna_index,
												dna_distinct_htable,
												previous_hashes)

		matched_seqids.extend(protein_matches[0])
		distinct_seqids_count += protein_matches[1]
		total_matches += protein_matches[2]
		distinct_protids += protein_matches[3]
		# Update locus mode
		if len(protein_matches[4]) > 0:
			for k, v in protein_matches[4].items():
				loci_modes[k][1].extend(v)
				loci_modes[k][0] = sm.determine_mode(loci_modes[k][1])[0]
		# Save protein hashes that were already processed
		# Different hash tables may contain the same hash because
		# different CDSs may encode the same protein
		previous_hashes = protein_matches[5]

	print(f'Found {distinct_protids} exact matches ({distinct_seqids_count} distinct CDSs, {total_matches} total CDSs).')

	# Save list of seqids that matched
	protein_exact_outdir = fo.join_paths(preprocess_dir, ['protein_exact_matching'])
	fo.create_directory(protein_exact_outdir)
	protein_exact_outfile = fo.join_paths(protein_exact_outdir, ['protein_exact_matches.txt'])
	fo.write_lines(matched_seqids, protein_exact_outfile)

	# Create protein file index
	protein_index = fao.index_fasta(representative_pfasta)
	# Create Fasta file without the protein that matched
	unique_pfasta = fo.join_paths(protein_exact_outdir, ['distinct_proteins.fasta'])
	unclassified_seqids = im.filter_list(representative_pseqids, matched_seqids)
	total_selected = fao.get_sequences_by_id(protein_index, unclassified_seqids, unique_pfasta)

	print(f'Unclassified proteins: {total_selected}')

	# User only wanted DNA and Protein exact matches or all sequences were classified
	if config['Mode'] == 2 or len(unclassified_seqids) == 0:
		template_dict['classification_files'] = classification_files
		template_dict['protein_fasta'] = unique_pfasta
		template_dict['unclassified_ids'] = unclassified_seqids
		template_dict['representatives'] = {}

		return template_dict

	# Translate schema representatives
	print(f'\n {ct.PROTEIN_CLUSTERING} ')
	print('='*(len(ct.PROTEIN_CLUSTERING)+2))
	print('Translating schema representative alleles...')
	rep_dir = fo.join_paths(schema_directory, ['short'])
	rep_full_list = fo.listdir_fullpath(rep_dir, '.fasta')
	reps_protein_dir = fo.join_paths(temp_directory, ['3_translated_representatives'])
	fo.create_directory(reps_protein_dir)
	protein_files = mo.parallelize_function(fao.translate_fasta, rep_full_list,
											[reps_protein_dir, config['Translation table']],
											config['CPU cores'], False)

	# Translated representative FASTA files mapped to loci basenames
	repprot_to_locibase = {file[1]: fo.file_basename(file[0], False).replace('_short', '') for file in protein_files}
	# Translated representative FASTA files for loci being called
	repprot_fastas = [k for k, v in repprot_to_locibase.items() if v in loci_basenames.values()]

	# Create directory to store clustering data
	clustering_dir = fo.join_paths(temp_directory, ['4_clustering'])
	fo.create_directory(clustering_dir)

	# Define BLASTp, makeblastdb and blastdb_aliastool paths
	blastp_path = fo.join_paths(config['BLAST path'], [ct.BLASTP_ALIAS])
	makeblastdb_path = fo.join_paths(config['BLAST path'], [ct.MAKEBLASTDB_ALIAS])
	blastdb_aliastool_path = fo.join_paths(config['BLAST path'], [ct.BLASTDB_ALIASTOOL_ALIAS])

	# Concatenate representative FASTA files
	concat_reps = fo.join_paths(reps_protein_dir, ['loci_to_call_translated_representatives.fasta'])
	fo.concatenate_files(repprot_fastas, concat_reps)

	# Determine self-score for representatives if file is missing
	self_score_file = fo.join_paths(schema_directory, ['short', 'self_scores'])
	if os.path.isfile(self_score_file) is False:
		print('Determining BLASTp self-score for each representative...')
		# Determine for all loci, not just target loci
		if len(repprot_to_locibase) > len(repprot_fastas):
			concat_full_reps = fo.join_paths(reps_protein_dir, ['schema_loci_translated_representatives.fasta'])
			fo.concatenate_files(list(repprot_to_locibase.keys()), concat_full_reps)
		else:
			concat_full_reps = concat_reps
		self_score_dir = fo.join_paths(reps_protein_dir, ['self_scores'])
		fo.create_directory(self_score_dir)
		self_scores = cf.determine_self_scores(concat_full_reps, self_score_dir,
											   makeblastdb_path, blastp_path,
											   'prot',
											   config['CPU cores'],
											   blastdb_aliastool_path)
		fo.pickle_dumper(self_scores, self_score_file)
	else:
		self_scores = fo.pickle_loader(self_score_file)
	print(f'Representative BLASTp self-scores stored in {self_score_file}')

	# Create Kmer index for representatives
	print('Creating minimizer index for representative alleles...')
	representatives = im.kmer_index(concat_reps, 5)
	print(f'Created index with {len(representatives)} distinct minimizers for {len(loci_files)} loci.')

	# Import unclassified proteins
	proteins = fao.import_sequences(unique_pfasta)
	# Cluster proteins into representative clusters
	print('Clustering proteins...')
	# Define input group size based on number of available CPU cores
	group_size = math.ceil(len(proteins)/config['CPU cores'])
	cs_results = cf.cluster_sequences(proteins,
									  config['Word size'],
									  config['Window size'],
									  config['Clustering similarity'],
									  representatives,
									  False,
									  1,
									  30,
									  clustering_dir,
									  config['CPU cores'],
									  group_size,
									  False)

	# Exclude singletons (>0 because clusters do not include representative)
	clusters = {k: v for k, v in cs_results.items() if len(v) > 0}
	# Determine number of clustered sequences
	distinct_clustered = set()
	for k, v in clusters.items():
		distinct_clustered = distinct_clustered.union(set([s[0] for s in v]))
	total_clustered = len(distinct_clustered)
	print(f'\nClustered {total_clustered} proteins into {len(clusters)} clusters.')
	print(f'{len(proteins)-total_clustered} proteins were not added to any cluster.')

	# Create Fasta file and index for unclassified proteins and schema representatives
	all_prots = fo.join_paths(clustering_dir, ['distinct_proteins.fasta'])
	fo.concatenate_files([unique_pfasta, concat_reps], all_prots)
	prot_index = fao.index_fasta(all_prots)

	# BLASTp if there are clusters with n>1
	excluded = []
	if len(clusters) > 0:
		# BLAST representatives against clustered sequences
		print('Aligning cluster representatives against clustered proteins...')
		blast_results, blast_results_dir = cf.blast_clusters(clusters, all_prots,
															 clustering_dir, blastp_path,
															 makeblastdb_path,
															 config['CPU cores'],
															 blastdb_aliastool_path
															 True)

		blast_files = im.flatten_list(blast_results)

		# Concatenate files based on locus identifier included in file paths
		blast_merged_dir = fo.join_paths(blast_results_dir, ['concatenated'])
		fo.create_directory(blast_merged_dir)
		concatenated_files = merge_blast_results(blast_files, blast_merged_dir, loci_finder)

		# Select best hit per target, filter based on BSR, expand matches
		# and get relevant data for classification
		loci_results = {}
		blast_matches_dir = fo.join_paths(clustering_dir, ['match_data'])
		fo.create_directory(blast_matches_dir)
		for file in concatenated_files:
			locus_id = fo.file_basename(file).split('_concatenated')[0]
			# Get best match per target
			best_matches = select_highest_scores(file)
			best_matches = list(best_matches.values())
			# Exclude results in the BSR+0.1 threshold
			# to process representative candidates in later stage
			match_info = process_blast_results(best_matches, config['BLAST Score Ratio']+0.1, self_scores)
			# Expand distinct protein matches to all inputs
			# that contain a CDS that encodes the protein
			locus_results = expand_matches(match_info, prot_index, dna_index,
										   dna_distinct_htable, distinct_pseqids, int_to_unique,
										   close_to_tip)

			if len(locus_results) > 0:
				# Save results to file
				locus_file = fo.join_paths(blast_matches_dir, ['{0}_matches'.format(locus_id)])
				fo.pickle_dumper(locus_results, locus_file)
				loci_results[locus_id] = locus_file

		# Classify matches
		if len(loci_results) > 0:
			print('\nClassifying high-scoring matches...')
			classification_inputs = []
			blast_clusters_results_dir = fo.join_paths(clustering_dir, ['classification_data'])
			fo.create_directory(blast_clusters_results_dir)
			for locus, file in loci_results.items():
				# Get locus length mode
				locus_mode = loci_modes[locus]
				# Import file with locus classifications
				locus_results_file = fo.join_paths(classification_dir, [f'{locus}_results'])
				# Create input lists
				classification_inputs.append([locus, file,
											  int_to_unique,
											  locus_results_file, locus_mode,
											  temp_directory,
											  config['Size threshold'],
											  config['BLAST Score Ratio'],
											  blast_clusters_results_dir,
											  config['CDS input'],
											  classify_inexact_matches])

			# Classify, parallelized per locus
			class_results = mo.map_async_parallelizer(classification_inputs,
													  mo.function_helper,
													  config['CPU cores'],
													  show_progress=True)

			# Update loci mode values and exclude sequences that were classified
			for r in class_results:
				current_results = fo.pickle_loader(r)
				for locus, v in current_results.items():
					loci_modes[locus] = v[1]
					excluded.extend(v[2])

			# May have repeated elements due to same CDS matching different loci
			excluded = set(excluded)

		print(f'\nClassified {len(excluded)} distinct proteins.')

	# Get seqids of remaining unclassified sequences
	unclassified_seqids = list(set(unclassified_seqids)-set(excluded))
	print(f'Unclassified proteins: {len(unclassified_seqids)}')

	# User only wanted exact matches and clustering or all sequences were classified
	if config['Mode'] == 3 or len(unclassified_seqids) == 0:
		template_dict['classification_files'] = classification_files
		template_dict['protein_fasta'] = all_prots
		template_dict['unclassified_ids'] = unclassified_seqids
		template_dict['representatives'] = {}

		return template_dict

	print(f'\n {ct.REPRESENTATIVE_DETERMINATION} ')
	print('='*(len(ct.REPRESENTATIVE_DETERMINATION)+2))
	print('Aligning representative alleles against unclassified proteins...')

	# Create directory to store data for each iteration
	iterative_rep_dir = fo.join_paths(temp_directory, ['5_representative_determination'])
	fo.create_directory(iterative_rep_dir)

	remaining_seqs_file = fo.join_paths(iterative_rep_dir, ['unclassified_proteins.fasta'])
	# Create Fasta with unclassified sequences
	fao.get_sequences_by_id(prot_index, unclassified_seqids,
							remaining_seqs_file, limit=50000)

	# Create BLAST DB
	blast_db_dir = fo.join_paths(iterative_rep_dir, ['BLASTp_db'])
	fo.create_directory(blast_db_dir)
	blast_db = fo.join_paths(blast_db_dir, ['unclassified_proteins'])
	db_std = bw.make_blast_db(makeblastdb_path, remaining_seqs_file, blast_db, 'prot')

	# Map representative allele header to locus ID
	rep_recs = fao.sequence_generator(concat_reps)
	rep_map = {}
	for rec in rep_recs:
		rep_map[rec.id] = '_'.join((rec.id).split('_')[:-1])

	# BLAST schema representatives against remaining unclassified CDSs
	new_reps = {}
	iteration = 1
	exausted = False
	# Keep iterating while there are sequences being classified
	rep_iter_header = '{:^11} {:^9} {:^14} {:^12} {:^10} {:^14}'.format('Iteration', 'Loci', 'High-Scoring',
								'Classified', 'Selected', 'Unclassified')
	print('='*len(rep_iter_header))
	print(rep_iter_header)
	print('='*len(rep_iter_header))
	while exausted is False:
		print('\r', f'{iteration}\t...', end='')
		# Create directory for current iteration
		iteration_directory = fo.join_paths(iterative_rep_dir, ['iteration_{0}'.format(iteration)])
		fo.create_directory(iteration_directory)
		# Create text file with unclassified seqids
		remaining_seqids_file = fo.join_paths(iteration_directory, ['unclassified_seqids_{0}.txt'.format(iteration)])
		fo.write_lines(unclassified_seqids, remaining_seqids_file)
		binary_file = f'{remaining_seqids_file}.bin'
		blastp_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
												remaining_seqids_file,
												binary_file)
		remaining_seqids_file = binary_file

		# BLAST representatives against unclassified sequences
		# Iterative process until no more sequences are classified
		print('\r', '{:^11} {:^9} {:^14}'.format(iteration, len(repprot_fastas), '...'), end='')

		# Concatenate to create groups of 100 loci
		iteration_blastin_dir = fo.join_paths(iteration_directory, ['BLASTp_infiles'])
		fo.create_directory(iteration_blastin_dir)
		concat_repfiles = []
		for i in range(0, len(repprot_fastas), 100):
			concat_file = fo.join_paths(iteration_blastin_dir, ['concat_reps{0}-{1}.fasta'.format(i+1, i+len(repprot_fastas[i:i+100]))])
			fo.concatenate_files(repprot_fastas[i:i+100], concat_file)
			concat_repfiles.append(concat_file)

		# Create BLASTp inputs
		output_files = []
		blast_inputs = []
		# Create directory to store BLASTp results
		iteration_blast_dir = fo.join_paths(iteration_directory, ['BLASTp_outfiles'])
		fo.create_directory(iteration_blast_dir)
		for file in concat_repfiles:
			concat_rep_basename = fo.file_basename(file, False)
			outfile = fo.join_paths(iteration_blast_dir,
									[concat_rep_basename+'_blastout_iter{0}.tsv'.format(iteration)])
			output_files.append(outfile)

			# If max_targets is set to None, BLAST defaults to 500
			blast_inputs.append([blastp_path, blast_db, file,
								 outfile, 1, 1,
								 remaining_seqids_file, 'blastp', 500,
								 None, bw.run_blast])

		# BLAST representatives against unclassified sequences
		blastp_results = mo.map_async_parallelizer(blast_inputs,
												   mo.function_helper,
												   config['CPU cores'],
												   show_progress=False)

		# Get BLASTp results per locus
		blast_merged_dir = fo.join_paths(iteration_blast_dir, ['concatenated'])
		fo.create_directory(blast_merged_dir)

		loci_output_files = []
		for f in output_files:
			concat_results = fo.read_tabular(f)
			loci_separate_results = {}
			# Get locus based on representative seqid
			for line in concat_results:
				loci_separate_results.setdefault(rep_map[line[0]], []).append(line)
			# Save results per locus
			for k, v in loci_separate_results.items():
				loci_separate_outfile = fo.join_paths(blast_merged_dir, ['{0}_concatenated_blastout_iter{1}.tsv'.format(k, iteration)])
				lines = ['\t'.join(l) for l in v]
				fo.write_lines(lines, loci_separate_outfile)
				loci_output_files.append(loci_separate_outfile)

		loci_results = {}
		# Create directory to store files with matches
		iteration_matches_dir = fo.join_paths(iteration_directory, ['match_data'])
		fo.create_directory(iteration_matches_dir)
		for file in loci_output_files:
			locus_id = loci_finder.search(file).group()
			best_matches = select_highest_scores(file)
			best_matches = list(best_matches.values())
			match_info = process_blast_results(best_matches, config['BLAST Score Ratio'], self_scores)
			locus_results = expand_matches(match_info, prot_index, dna_index,
										   dna_distinct_htable, distinct_pseqids, int_to_unique,
										   close_to_tip)

			if len(locus_results) > 0:
				locus_file = fo.join_paths(iteration_matches_dir, ['{0}_matches'.format(locus_id)])
				fo.pickle_dumper(locus_results, locus_file)
				loci_results[locus_id] = locus_file

		print('\r', '{:^11} {:^9} {:^14} {:^12}'.format(iteration, len(repprot_fastas), len(loci_results), '...'), end='')

		if len(loci_results) == 0:
			exausted = True
			print('\r', '{:^11} {:^9} {:^14} {:^12} {:^10} {:^14}'.format(iteration, len(repprot_fastas), len(loci_results), 0, 0, len(unclassified_seqids)))
			continue

		# Process results per genome and per locus
		classification_inputs = []
		blast_iteration_results_dir = fo.join_paths(iteration_directory, ['classification_data'])
		fo.create_directory(blast_iteration_results_dir)
		for locus, file in loci_results.items():
			# Get locus length mode
			locus_mode = loci_modes[locus]
			# Import file with locus classifications
			locus_results_file = fo.join_paths(classification_dir, [locus+'_results'])

			classification_inputs.append([locus, file,
										  int_to_unique,
										  locus_results_file, locus_mode,
										  temp_directory,
										  config['Size threshold'],
										  config['BLAST Score Ratio'],
										  blast_iteration_results_dir,
										  config['CDS input'],
										  classify_inexact_matches])

		class_results = mo.map_async_parallelizer(classification_inputs,
												  mo.function_helper,
												  config['CPU cores'],
												  show_progress=False)

		# Need to identify representative candidates that match several
		# loci and remove them from the analysis
		excluded = []
		representative_candidates = {}
		for r in class_results:
			current_results = fo.pickle_loader(r)
			for locus, v in current_results.items():
				# Update allele size mode value
				loci_modes[locus] = v[1]
				# Exclude classified seqids
				# This includes seqids for representative candidates
				excluded.extend(v[2])
				# Add representative candidates
				if len(v[3]) > 0:
					representative_candidates[locus] = v[3]

		# May have repeated elements due to same CDS matching different loci
		excluded = set(excluded)

		print('\r', '{:^11} {:^9} {:^14} {:^12} {:^10}'.format(iteration, len(repprot_fastas), len(loci_results), len(excluded), '...'), end='')

		# Exclude sequences that were excluded
		unclassified_seqids = set(unclassified_seqids) - excluded

		# Create directory to store new representatives
		new_reps_directory = fo.join_paths(iteration_directory, ['representative_candidates'])
		fo.create_directory(new_reps_directory)

		# Select new representatives for next iteration
		representatives = {}
		representative_inputs = []
		total_representatives = 0
		if len(representative_candidates) > 0:
			candidates_dir = fo.join_paths(new_reps_directory, ['candidates'])
			fo.create_directory(candidates_dir)
			selection_dir = fo.join_paths(new_reps_directory, ['selection'])
			fo.create_directory(selection_dir)
			blast_selection_dir = fo.join_paths(selection_dir, ['BLASTp_results'])
			fo.create_directory(blast_selection_dir)
			for k, v in representative_candidates.items():
				current_candidates = {e[1]: e[3] for e in v}
				fasta_file = fo.join_paths(candidates_dir,
										   ['{0}_candidates.fasta'.format(k)])
				# Create file with sequences
				fao.get_sequences_by_id(prot_index, list(current_candidates.keys()), fasta_file)
				# If multiple candidates, compare and select
				if len(v) > 1:
					representative_inputs.append([current_candidates, k, fasta_file,
												  iteration, blast_selection_dir, blastp_path,
												  blast_db, config['BLAST Score Ratio'], 1,
												  blastdb_aliastool_path, select_representatives])
				# Single candidate
				else:
					representatives[k] = [(v[0][1], v[0][3])]

			# Select new representatives for loci with multiple candidates
			selected_candidates = mo.map_async_parallelizer(representative_inputs,
															mo.function_helper,
															config['CPU cores'],
															show_progress=False)

			for c in selected_candidates:
				representatives[c[0]] = c[1]

			for k, v in representatives.items():
				new_reps.setdefault(k, []).extend(v)

			total_representatives = sum([len(v) for k, v in representatives.items()])

		print('\r', '{:^11} {:^9} {:^14} {:^12} {:^10} {:^14}'.format(iteration, len(repprot_fastas), len(loci_results), len(excluded), total_representatives, '...'), end='')
		print('\r', '{:^11} {:^9} {:^14} {:^12} {:^10} {:^14}'.format(iteration, len(repprot_fastas), len(loci_results), len(excluded), total_representatives, len(unclassified_seqids)))

		# Stop iterating if there are no new representatives
		if len(representatives) == 0:
			exausted = True
		# Prepare representative data for next iteration
		else:
			selected_dir = fo.join_paths(new_reps_directory, ['selected'])
			fo.create_directory(selected_dir)
			# Create files with representative sequences
			rep_map = {}
			reps_ids = []
			repprot_fastas = []
			for k, v in representatives.items():
				# Get new representative for locus
				current_new_reps = [e[0] for e in v]
				reps_ids.extend(current_new_reps)
				for c in current_new_reps:
					rep_map[c] = k

				# Need to add 'short' or locus id will not be split
				rep_file = fo.join_paths(selected_dir,
										 ['{0}_short_reps_iter.fasta'.format(k)])
				fao.get_sequences_by_id(prot_index, current_new_reps, rep_file)
				repprot_fastas.append(rep_file)

			# Concatenate reps
			concat_repy = fo.join_paths(new_reps_directory, ['concat_reps.fasta'])
			fao.get_sequences_by_id(prot_index, set(reps_ids), concat_repy, limit=50000)
			# Determine self-score for new reps
			candidates_blast_dir = fo.join_paths(new_reps_directory, ['representatives_self_score'])
			fo.create_directory(candidates_blast_dir)
			new_self_scores = cf.determine_self_scores(concat_repy, candidates_blast_dir,
													   makeblastdb_path, blastp_path,
													   'prot', config['CPU cores'],
													   blastdb_aliastool_path)

			# This includes self-score for candidates that are not added
			# (e.g. classification changes due to multiple matches)
			self_scores = {**self_scores, **new_self_scores}

		iteration += 1

	print('='*len(rep_iter_header))

	template_dict['classification_files'] = classification_files
	template_dict['protein_fasta'] = all_prots
	template_dict['unclassified_ids'] = unclassified_seqids
	template_dict['self_scores'] = self_scores
	template_dict['representatives'] = new_reps

	return template_dict


def main(input_file, loci_list, schema_directory, output_directory,
		 no_inferred, output_unclassified, output_missing, output_novel,
		 no_cleanup, hash_profiles, ns, config):

	start_time = pdt.get_datetime()

	print(f' {ct.CONFIG_VALUES} ')
	print('='*(len(ct.CONFIG_VALUES)+2))

	if config['Prodigal mode'] == 'meta' and config['Prodigal training file'] is not None:
		print('Prodigal mode is set to "meta". Will not use '
			  f'the training file in {config["Prodigal training file"]}')
		config['Prodigal training file'] = None

	# Print config parameters
	for k, v in config.items():
		print('{0}: {1}'.format(k, v))

	# Read list of paths to input FASTA files
	input_files = fo.read_lines(input_file, strip=True)
	# Sort file paths
	input_files = im.sort_iterable(input_files, sort_key=str.lower)
	print('Number of inputs: {0}'.format(len(input_files)))

	# Get list of loci in the schema
	schema_loci = fo.pickle_loader(fo.join_paths(schema_directory,
												 [ct.GENE_LIST_BASENAME]))
	schema_loci_paths = [fo.join_paths(schema_directory, [file])
							for file in schema_loci]

	# Read list of loci to call
	loci_to_call = fo.read_lines(loci_list)
	loci_to_call = {file: schema_loci_paths.index(file)
					for file in loci_to_call}
	print('Number of loci: {0}'.format(len(loci_to_call)))

	# Create regex compiler to find longest locus identifier in paths/strings
	loci_ids = [fo.file_basename(file, False) for file in loci_to_call]
	# Sort to find longest first
	loci_ids = sorted(loci_ids, key=lambda x: len(x), reverse=True)
	# Create regex object to search for loci identifiers in paths/strings
	loci_finder = re.compile('|'.join(loci_ids))

	# Create directory to store intermediate files
	temp_directory = fo.join_paths(output_directory, ['temp'])
	fo.create_directory(temp_directory)
	print(f'Intermediate files will be stored in {temp_directory}')

	print(f'\n {ct.PRECOMPUTED_DATA} ')
	print('='*(len(ct.PRECOMPUTED_DATA)+2))

	# Read or compute locus allele size mode
	loci_modes_file = fo.join_paths(schema_directory, ['loci_modes'])
	if os.path.isfile(loci_modes_file) is True:
		loci_modes = fo.pickle_loader(loci_modes_file)
	else:
		print('Determining allele size mode for all loci...')
		# Compute for all loci, not just for the subset of loci to call
		loci_modes = compute_loci_modes(schema_loci_paths, loci_modes_file)
	print(f'Loci allele size mode values stored in {loci_modes_file}')

	# Check if schema contains folder with pre-computed hash tables
	pre_computed_dir = fo.join_paths(schema_directory, ['pre_computed'])
	if os.path.isdir(pre_computed_dir) is False:
		print('Could not find pre-computed hash tables used for exact matching.')
		print('Creating hash tables...')
		# Create hash tables for DNA and Protein exact matching
		# This avoids translating all the schema alleles in each run
		fo.create_directory(pre_computed_dir)
		precompute_hash_tables(pre_computed_dir,
							   schema_loci_paths,
							   config['Translation table'],
							   config['CPU cores'])
	print(f'Hash tables stored in {pre_computed_dir}')

	# Perform allele calling
	results = allele_calling(input_files, schema_directory,
							 temp_directory, loci_modes.copy(),
							 loci_to_call, config, pre_computed_dir,
							 loci_finder)

	# Assign allele identifiers, add alleles to schema and create output files
	print(f'\n {ct.WRAPPING_UP} ')
	print('='*(len(ct.WRAPPING_UP)+2))

	# Adjust missing locus classification based on mode
	classification_labels = ct.ALLELECALL_CLASSIFICATIONS
	# Only mode 4 performs a final exhaustive search
	# Other modes do not and LNF classifications are considered Probable LNFs (PLNF)
	if config['Mode'] != 4:
		classification_labels[-1] = ct.PROBABLE_LNF
		print(f"Running in mode {config['Mode']}. Renaming Locus Not Found (LNF) "
			  "class to Probable LNF (PLNF).")

	# Sort to get output order similar to chewBBACA v2
	results['classification_files'] = dict(sorted(results['classification_files'].items()))

	print(f'Creating file with genome coordinates profiles ({ct.RESULTS_COORDINATES_BASENAME})...')
	results_contigs = write_results_contigs(list(results['classification_files'].values()),
											results['int_to_unique'],
											output_directory,
											results['cds_coordinates'],
											classification_labels,
											loci_finder)
	outfile, repeated_info, repeated_counts = results_contigs

	# Identify paralogous loci
	print('Identifying paralogous loci and creating files with the list of paralogous '
		  f'loci ({ct.PARALOGOUS_COUNTS_BASENAME} & {ct.PARALOGOUS_LOCI_BASENAME})...')
	paralogous_info = identify_paralogous(repeated_info, output_directory)
	print(f'Identified {paralogous_info[0]} paralogous loci.')

	# Assign allele identifiers to inferred alleles
	print('Assigning allele identifiers to inferred alleles...')
	# Create directory to store data for novel alleles
	novel_directory = fo.join_paths(temp_directory, ['novel_alleles'])
	novel_data_directory = fo.join_paths(novel_directory, ['data'])
	fo.create_directory(novel_data_directory)
	assignment_inputs = list(results['classification_files'].items())
	repeated_hashes = set(repeated_info.keys())
	assignment_inputs = [[g, ns, repeated_hashes, novel_data_directory, loci_finder, assign_allele_ids]
						 for g in assignment_inputs]

	novel_alleles = mo.map_async_parallelizer(assignment_inputs,
											  mo.function_helper,
											  config['CPU cores'],
											  show_progress=False)
	# Only keep data for loci that have novel alleles
	novel_alleles = [r for r in novel_alleles if r is not None]
	novel_alleles_count = sum([locus_novel[2] for locus_novel in novel_alleles])
	print(f'Assigned identifiers to {novel_alleles_count} new alleles for {len(novel_alleles)} loci.')

	updated_files = {}
	if config['Mode'] != 1:
		print('Getting original sequence identifiers for new alleles...')
		for locus_novel in novel_alleles:
			current_novel = fo.pickle_loader(locus_novel[1])
			for l in current_novel:
				# Get seqids that match hashes
				rep_seqid = im.polyline_decoding(results['dna_hashtable'][l[0]])[0:2]
				rep_seqid = '{0}-protein{1}'.format(results['int_to_unique'][rep_seqid[1]], rep_seqid[0])
				l.append(rep_seqid)
			fo.pickle_dumper(current_novel, locus_novel[1])

		reps_info = {}
		if config['Mode'] == 4:
			print('Getting data for new representative alleles...')
			# Get info for new representative alleles that must be added to files in the short directory
			for locus_novel in novel_alleles:
				locus_id = loci_finder.search(locus_novel[1]).group()
				current_novel = fo.pickle_loader(locus_novel[1])
				current_results = results['representatives'].get(locus_id, None)
				if current_results is not None:
					for e in current_results:
						allele_id = [line[1] for line in current_novel if line[0] == e[1]]
						# We might have representatives that were converted to NIPH but still appear in the list
						if len(allele_id) > 0:
							reps_info.setdefault(locus_id, []).append(list(e)+allele_id)

			if no_inferred is False:
				self_score_file = fo.join_paths(schema_directory, ['short', 'self_scores'])
				print(f'Adding the BLASTp self-score for the new representatives to {self_score_file}')
				# Update self_scores
				reps_to_del = set()
				for k, v in reps_info.items():
					for r in v:
						new_id = k+'_'+r[-1]
						results['self_scores'][new_id] = results['self_scores'][r[0]]
						# Delete old entries
						# Does not delete entries from representative candidates that were converted to NIPH
						if r[0] not in reps_to_del:
							reps_to_del.add(r[0])

				for r in reps_to_del:
					del results['self_scores'][r]

				# Save updated self-scores
				fo.pickle_dumper(results['self_scores'], self_score_file)

		if len(novel_alleles) > 0:
			print('Creating FASTA files with the new alleles...')
			# Create Fasta files with novel alleles
			novel_fastas_directory = fo.join_paths(novel_directory, ['novel_fastas'])
			novel_rep_directory = fo.join_paths(novel_fastas_directory, ['short'])
			fo.create_directory(novel_rep_directory)
			novel_data = create_novel_fastas(novel_alleles, reps_info, results['dna_fasta'], novel_fastas_directory, loci_finder)
			total_inferred, total_representatives, updated_novel = novel_data
			updated_files = updated_novel
			if no_inferred is False:
				print('Adding new alleles to schema...')
				# Add inferred alleles to schema
				alleles_added = add_inferred_alleles(updated_novel, loci_finder)
				# Recompute mode for loci with novel alleles
				print(f'Updating allele size mode values stored in {loci_modes_file}')
				for locus_novel in novel_alleles:
					alleles_sizes = list(fao.sequence_lengths(locus_novel[0]).values())
					# Select first value in list if there are several values with same frequency
					loci_modes[fo.file_basename(locus_novel[0], False)] = [sm.determine_mode(alleles_sizes)[0], alleles_sizes]
				fo.pickle_dumper(loci_modes, loci_modes_file)
				# Add novel alleles hashes to pre-computed hash tables
				print(f'Updating pre-computed hash tables in {pre_computed_dir}')
				total_hashes = update_hash_tables(updated_novel, loci_to_call,
								   config['Translation table'], pre_computed_dir)

	# Create file with allelic profiles
	print(f'Creating file with the allelic profiles ({ct.RESULTS_ALLELES_BASENAME})...')
	profiles_table = write_results_alleles(list(results['classification_files'].values()),
										   list(results['int_to_unique'].values()),
										   output_directory,
										   classification_labels[-1],
										   loci_finder)

	# Create file with class counts per input file
	print(f'Creating file with class counts per input ({ct.RESULTS_STATISTICS_BASENAME})...')
	input_stats_file = write_results_statistics(results['classification_files'],
												results['int_to_unique'],
												results['cds_counts'],
												output_directory,
												classification_labels,
												repeated_counts,
												results['invalid_alleles'])

	# Create file with class counts per locus called
	print(f'Creating file with class counts per locus ({ct.LOCI_STATS_BASENAME})...')
	loci_stats_file = write_loci_summary(results['classification_files'],
										 output_directory,
										 len(input_files),
										 classification_labels,
										 loci_finder)

	# Create FASTA file with unclassified CDSs
	if output_unclassified is True:
		print(f'Creating file with unclassified CDSs ({ct.UNCLASSIFIED_BASENAME})...')
		unclassified_file = create_unclassified_fasta(results['dna_fasta'],
													  results['protein_fasta'],
													  results['unclassified_ids'],
													  results['protein_hashtable'],
													  output_directory,
													  results['int_to_unique'])

	# Create FASTA and TSV files with data about the CDSs classified as non-EXC/non-INF
	if output_missing is True:
		print('Creating files with data about the CDSs classified as non-EXC/non-INF '
			  f'({ct.MISSING_FASTA_BASENAME} & {ct.MISSING_TSV_BASENAME})...')
		missing_outfiles = create_missing_fasta(results['classification_files'],
												results['dna_fasta'],
												results['int_to_unique'],
												results['dna_hashtable'],
												output_directory,
												classification_labels,
												config['CDS input'],
												loci_finder)

	# Create FASTA file with inferred alleles
	if len(novel_alleles) > 0 and output_novel is True:
		print(f'Creating file with new alleles ({ct.NOVEL_BASENAME})...')
		novel_fasta = fo.join_paths(output_directory, [ct.NOVEL_BASENAME])
		novel_files = [v[0] for v in updated_files.values()]
		# Concatenate all FASTA files with inferred alleles
		fo.concatenate_files(novel_files, novel_fasta)

	# Create TSV file with hashed profiles
	if hash_profiles is not None:
		print(f'Creating file with {hash_profiles} hashed profiles...')
		hashed_profiles_file = ph.main(profiles_table, schema_directory, output_directory,
									   hash_profiles, config['CPU cores'], 100, updated_files,
									   no_inferred)

	# Create TSV file with CDS coordinates
	# Will not be created if input files contain set of CDS instead of contigs
	if config['CDS input'] is False:
		print(f'Creating file with the coordinates of CDSs identified in inputs ({ct.CDS_COORDINATES_BASENAME})...')
		files = []
		for gid, file in results['cds_coordinates'].items():
			tsv_file = fo.join_paths(os.path.dirname(file), [f'{gid}_coordinates.tsv'])
			cf.write_coordinates_file(file, tsv_file)
			files.append(tsv_file)
		# Concatenate all TSV files with CDS coordinates
		cds_coordinates = fo.join_paths(output_directory,
										[ct.CDS_COORDINATES_BASENAME])
		fo.concatenate_files(files, cds_coordinates,
							 header=ct.CDS_TABLE_HEADER)

	# Move file with list of excluded CDS
	# File is not created if we only search for exact matches
	if config['Mode'] != 1:
		print(f'Creating file with invalid CDSs ({ct.INVALID_CDS_BASENAME})...')
		fo.move_file(results['invalid_alleles'][0], output_directory)

	# Count total for each classification type
	print('Counting number of classified CDSs...')
	global_counts, total_cds = count_global_classifications(results['classification_files'].values(),
													 classification_labels)

	# Subtract number of times a CDSs is repeated
	print(f'Classified a total of {total_cds-sum(repeated_counts.values())} CDSs.')
	class_counts_header = [f'{c:^8}' for c in ct.ALLELECALL_CLASSIFICATIONS[:-1]]
	class_counts_header = ' '.join(class_counts_header)
	class_counts_values = [f'{global_counts.get(c, 0):^8}' for c in ct.ALLELECALL_CLASSIFICATIONS[:-1]]
	class_counts_values = ' '.join(class_counts_values)
	print('='*len(class_counts_header))
	print(class_counts_header)
	print('='*len(class_counts_header))
	print(class_counts_values)
	print('='*len(class_counts_header))

	if no_inferred is False and config['Mode'] != 1 and len(novel_alleles) > 0:
		print('Added {0} new alleles to the schema.'.format(total_inferred))
		print('Added {0} new representative alleles to the schema.'.format(total_representatives))
	elif no_inferred is True:
		print('User passed "--no-inferred". No alleles added to the schema.')
	elif len(novel_alleles) == 0:
		print('No newly inferred alleles to add to schema.')

	# Remove temporary files
	if no_cleanup is False:
		print('Removing temporary directory with intermediate files...')
		fo.delete_directory(temp_directory)

	end_time = pdt.get_datetime()

	# Create basic log file
	print(f'Creating log file ({ct.LOGFILE_BASENAME})...')
	logfile_path = write_logfile(start_time,
								 end_time,
								 len(results['int_to_unique']),
								 len(results['classification_files']),
								 config['CPU cores'],
								 config['BLAST Score Ratio'],
								 output_directory)

	print(f'\nResults available in {output_directory}')
