#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module enables the adaptation of external schemas so that the loci and
alleles present in those schemas can be used with chewBBACA. During the
process, alleles that do not correspond to a complete CDS or that cannot be
translated are discarded from the final schema. One or more alleles of each
gene/locus will be chosen as representatives and included in the 'short'
directory of the adapted schema.

Code documentation
------------------
"""

import os
import shutil
import itertools
import multiprocessing

from loguru import logger

import CHEWBBACA.global_constants as gct
import CHEWBBACA.AdaptSchema.constants as ct
from CHEWBBACA.utils import (blast_wrapper as bw,
								core_functions as cf,
								file_operations as fo,
								fasta_operations as fao,
								sequence_manipulation as sm,
								iterables_manipulation as im,
								multiprocessing_operations as mo)


def bsr_categorizer(blast_results, representatives,
					representatives_scores, min_bsr, max_bsr):
	"""Identify BLAST hits below and above the BSR min and max thresholds.

	Parameters
	----------
	blast_results : list of list
		A list with sublists, each sublist contains information
		about a BLAST hit.
	representatives : list
		List with sequence identifiers of representative
		sequences.
	representatives_scores : dict
		Dictionary with self BLAST raw score for every
		representative.
	min_bsr : float
		Minimum BSR value accepted to consider a sequence as
		a possible new representative.
	max_bsr : float
		Maximum BSR value accepted to consider a sequence as
		a possible new representative.

	Returns
	-------
	List with the following elements:
		high_bsr : list
			list with all sequence identifiers of subject
			sequences that had hits with a BSR higher than the
			maximum defined threshold.
		low_bsr : list
			list with all sequence identifiers of subject
			sequences that had hits with a BSR lower than the
			minimum defined threshold.
	"""
	high_bsr = []
	hotspot_bsr = []
	low_bsr = []

	high_reps = {}
	hot_reps = {}
	low_reps = {}

	filtered_results = [res for res in blast_results
						if res[0] != res[4] and res[4] not in representatives]
	bsr_values = [float(res[-1])/representatives_scores[res[0]]
				  for res in filtered_results]

	high_bsr = [res[4] for ind, res in enumerate(filtered_results)
				if bsr_values[ind] >= max_bsr]
	low_bsr = [res[4] for ind, res in enumerate(filtered_results)
			   if bsr_values[ind] < min_bsr]
	hotspot_bsr = [res[4] for ind, res in enumerate(filtered_results)
				   if bsr_values[ind] >= min_bsr and bsr_values[ind] < max_bsr]

	for ind, res in enumerate(filtered_results):
		if bsr_values[ind] >= min_bsr:
			high_reps.setdefault(res[0], []).append(res[4])
		if bsr_values[ind] < min_bsr:
			low_reps.setdefault(res[0], []).append(res[4])
		if bsr_values[ind] >= min_bsr and bsr_values[ind] < max_bsr:
			hot_reps.setdefault(res[0], []).append(res[4])

	# Identify representatives that only led to low BSR
	low_reps = list(set(low_reps) - set(high_reps))

	return [high_bsr, low_bsr, hotspot_bsr, high_reps, low_reps, hot_reps]


def select_candidate(candidates, proteins, seqids,
					 representatives, final_representatives):
	"""Select a new representative sequence.

	Parameters
	----------
	candidates : list
		List with the sequence identifiers of all candidates.
	proteins : dict
		A dictionary with protein identifiers as keys and
		protein sequences as values.
	seqids : list
		A list with the sequence identifiers that still have
		no representative (representatives identifiers are
		included because they have to be BLASTed in order to
		determine their self score).
	representatives : list
		The sequence identifiers of all representatives.

	Returns
	-------
	representatives : list
		The set of all representatives, including the new
		representative that was chosen by the function.
	"""
	# With more than one sequence as candidate, select longest
	if len(candidates) > 1:
		# Determine length of all candidates
		candidates_len = [(seqid, len(proteins[seqid]))
						  for seqid in candidates]

		# Order representative candidates by length descending order
		candidates_len = sorted(candidates_len, key=lambda x: x[1],
								reverse=True)

		# Longest allele is the new representative
		representatives.append(candidates_len[0][0])
		final_representatives.append(candidates_len[0][0])
	# If there is only one candidate, keep that
	elif len(candidates) == 1:
		representatives.append(candidates[0])
		final_representatives.append(candidates[0])
	# If no hit qualifies and there are still sequences
	# without representative
	elif len(candidates) == 0 and \
			len(seqids) > len(representatives):
		# Determine length of remaining sequences
		# (representatives not included)
		candidates_len = [(seqid, len(proteins[seqid]))
						  for seqid in seqids
						  if seqid not in representatives]
		# Sort by descending length
		candidates_len = sorted(candidates_len, key=lambda x: x[1],
								reverse=True)
		# Longest of remaining sequences is new representative
		representatives.append(candidates_len[0][0])
		final_representatives.append(candidates_len[0][0])

	return [representatives, final_representatives]


def adapt_loci(loci, schema_path, schema_short_path, bsr, min_len,
			   table_id, size_threshold, blastp_path, makeblastdb_path,
			   blastdb_aliastool_path):
	"""Adapts a set of loci from an external schema.

	Adapts an external schema for usage with chewBBACA. Removes invalid
	alleles and selects representative alleles to include in the "short"
	directory.

	Parameters
	----------
	genes_list : list
		A list with the following elements:

		- List with paths to the files to be processed.
		- Path to the schema directory.
		- Path to the "short" directory.
		- BLAST Score Ratio value.
		- Minimum sequence length value.
		- Genetic code.
		- Sequence size variation threshold.

	Returns
	-------
	invalid_alleles : list
		List with the identifiers of the alleles that were
		determined to be invalid.
	invalid_loci : list
		List with the identifiers of the loci that had no
		valid alleles.
	summary_stats : list of list
		List with one sublist per processed locus. Each
		sublist has four elements:

		- The identifier of the locus.
		- The number of alleles in the external file.
		- The number of alleles that were a valid CDS.
		- The number of representatives determined determined
		  by the process.

	The function writes the schema files .
	"""
	summary_stats = []
	invalid_loci = []
	invalid_alleles = []
	for locus in loci:
		representatives = []
		final_representatives = []
		rep_self_scores = {}

		# Get locus identifier (does not include extension)
		locus_id = fo.file_basename(locus, file_extension=False)

		# Create paths to gene files in new schema
		locus_file = fo.join_paths(schema_path, [f'{locus_id}.fasta'])
		locus_short_file = fo.join_paths(schema_short_path,
										 [f'{locus_id}_short.fasta'])

		# Create temp directory for current gene
		locus_temp_dir = fo.join_paths(schema_path, [f'{locus_id}_temp'])
		# Create temp directory for the current gene
		fo.create_directory(locus_temp_dir)

		# Translate sequences
		dna_file, protein_file, _, invalid = fao.translate_fasta(locus, locus_temp_dir, table_id, True)
		# Get seqids of alleles that could not be translated
		excluded = set([s[0] for s in invalid])

		# FASTA file included valid alleles
		if protein_file:
			# Import DNA sequences only after translating to avoid importing sequences that cannot be translated
			dna_seqs = fao.import_sequences(dna_file)
			total_sequences = len(dna_seqs) + len(invalid)
			prot_seqs = fao.import_sequences(protein_file)
			# Get the IDs of alleles that could not be translated

			if size_threshold is not None and len(prot_seqs) > 0:
				# Remove alleles based on length mode and size threshold
				modes, alm, asm, allele_sizes = sm.mode_filter(dna_seqs, size_threshold)
				excluded.add(set(asm + alm))

				modes_concat = ':'.join(map(str, modes))
				st_percentage = int(size_threshold*100)
				invalid += [[s, gct.ALM_MSG.format(st_percentage, allele_sizes[s], modes_concat)] for s in alm]
				invalid += [[s, gct.ASM_MSG.format(st_percentage, allele_sizes[s], modes_concat)] for s in asm]

			# Remove excluded alleles
			dna_seqs = {seqid: seq
						for seqid, seq in dna_seqs.items()
						if seqid not in excluded}
			prot_seqs = {seqid: seq
						for seqid, seq in prot_seqs.items()
						if seqid not in excluded}
		# If no sequences could be translated, remove temp directory, add locus to invalid loci and continue to next locus
		else:
			dna_seqs = {}
			prot_seqs = {}
			total_sequences = len(invalid)

		# Prepare list of invalid alleles to write
		# Add locus identifier as prefix to allele IDs if it is not included in the original ID
		invalid = [[i[0] if locus_id in i[0] else f'{locus_id}_{i[0]}', i[1]] for i in invalid]
		invalid_alleles.extend(invalid)

		# Continue to next locus if there are no valid CDSs for current locus
		if len(prot_seqs) == 0:
			shutil.rmtree(locus_temp_dir)
			invalid_loci.append(locus_id)
			summary_stats.append([locus_id, str(total_sequences), '0', '0'])
			continue

		if len(dna_seqs) > 1:
			# Identify DNA sequences that code for same protein
			equal_prots = sm.determine_duplicated_seqs(prot_seqs)

			# Get only one identifier per protein
			ids_to_blast = [protids[0] for protein, protids in equal_prots.items()]

			# Get longest sequence as first representative
			longest = sm.determine_longest(ids_to_blast, prot_seqs)
			representatives.append(longest)
			final_representatives.append(longest)

			# Create FASTA file with distinct protein sequences
			protein_file = fo.join_paths(locus_temp_dir,
										 [f'{locus_id}_protein.fasta'])
			protein_data = [[i, prot_seqs[i]] for i in ids_to_blast]
			protein_lines = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, protein_data)
			fo.write_lines(protein_lines, protein_file)

			# Shorten sequence IDs to avoid issues with long identifiers when creating BLAST DBs
			lcl_proteins = fo.join_paths(locus_temp_dir, [protein_file.replace('.fasta', '_LCL.fasta')])
			# Use 'lcl|SEQ' as prefix to avoid issues where makeblastdb modifies the IDs
			_ = fao.integer_headers(protein_file, lcl_proteins, start=1, limit=50000, prefix=gct.BLASTDB_LCL_PREFIX, id_map=True)
			# Create BLASTdb with all distinct proteins
			blastp_db = os.path.join(locus_temp_dir, locus_id)
			db_std = bw.make_blast_db(makeblastdb_path, lcl_proteins, blastp_db, 'prot')
			# Delete FASTA file with LCL IDs
			fo.remove_files([lcl_proteins])

			# Create a second FASTA file with the sequence ID format returned by BLAST ('SEQ' instead of 'lcl|SEQ')
			seq_proteins = fo.join_paths(locus_temp_dir, [protein_file.replace('.fasta', '_SEQ.fasta')])
			id_mapping = fao.integer_headers(protein_file, seq_proteins, start=1, limit=50000, prefix=gct.BLASTDB_SEQ_PREFIX, id_map=True)
			# Return inverse mapping to convert original IDs to renamed IDs for BLASTp
			inverse_id_mapping = im.invert_dictionary(id_mapping)

			# Determine appropriate blastp task (proteins < 30aa need blastp-short)
			blastp_task = bw.determine_blast_task(equal_prots)

			# Cycle BLAST representatives against non-representatives until
			# all non-representatives have a representative
			while len(set(ids_to_blast) - set(representatives)) != 0:
				# Create FASTA file with representative alleles
				rep_file = fo.join_paths(locus_temp_dir,
										 [f'{locus_id}_rep_protein.fasta'])
				rep_protein_data = [[r, prot_seqs[r]] for r in representatives]
				rep_protein_lines = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, rep_protein_data)
				fo.write_lines(rep_protein_lines, rep_file)

				# Compute self-score for representative alleles
				for seqid in representatives:
					if seqid not in rep_self_scores:
						record = fao.fasta_str_record(gct.FASTA_RECORD_TEMPLATE, [seqid, prot_seqs[seqid]])
						current_rep_file = fo.join_paths(locus_temp_dir, [f'{seqid}_solo.fasta'])
						fo.write_lines([record], current_rep_file)
						# Create file with representative seqid to only compare against self
						id_file = fo.join_paths(locus_temp_dir, [f'{seqid}_ids.txt'])
						fo.write_lines([inverse_id_mapping[seqid]], id_file)
						binary_file = f'{id_file}.bin'
						blast_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
															id_file,
															binary_file)
						id_file = binary_file

						rep_blastout = fo.join_paths(locus_temp_dir, [f'{seqid}_blastout.tsv'])
						# Cannot get self-alignemnt for some sequences if composition-based stats is enabled
						blast_std = bw.run_blast(blastp_path, blastp_db, current_rep_file,
												rep_blastout, 1, 1,
												id_file, 'blastp', None, 0)
						rep_results = fo.read_tabular(rep_blastout)
						# Convert short ID back to original ID
						rep_results = [[r[0]]+r[1:4]+[id_mapping[r[4]]]+r[5:] for r in rep_results]
						if len(rep_results) > 0:
							rep_self_scores[rep_results[0][0]] = float(rep_results[0][6])
						else:
							logger.error("Could not determine the self-alignment raw score for {rep_results[0][0]}")

				# Create file with seqids to BLAST against
				ids_str = im.join_list([str(inverse_id_mapping[i]) for i in ids_to_blast if i not in representatives], '\n')
				ids_file = fo.join_paths(locus_temp_dir, [f'{locus_id}_ids.txt'])
				fo.write_to_file(ids_str, ids_file, 'w', '')
				binary_file = f'{ids_file}.bin'
				blast_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
														ids_file,
														binary_file)
				ids_file = binary_file

				# BLAST representatives against non-represented
				blast_output = fo.join_paths(locus_temp_dir,
											 [f'{locus_id}_blast_out.tsv'])
				# Set 'max_target_seqs' to huge number because BLAST only
				# returns 500 hits by default
				blast_std = bw.run_blast(blastp_path, blastp_db, rep_file,
										 blast_output, 1, 1, ids_file,
										 blastp_task, 100000)

				# Import BLAST results
				blast_results = fo.read_tabular(blast_output)
				# Convert short ID back to original ID
				blast_results = [[r[0]]+r[1:4]+[id_mapping[r[4]]]+r[5:] for r in blast_results]

				# Divide results into high, low and hot BSR values
				hitting_high, hitting_low, hotspots, high_reps, low_reps, hot_reps = \
					bsr_categorizer(blast_results, representatives,
									rep_self_scores, bsr, bsr+0.1)

				excluded_reps = []
				# Remove high BSR hits that have representative
				hitting_high = set(hitting_high)
				ids_to_blast = [i for i in ids_to_blast if i not in hitting_high]

				# Remove representatives that led to high BSR with subjects that were removed
				prunned_high_reps = {k: [r for r in v if r in ids_to_blast] for k, v in high_reps.items()}
				reps_to_remove = [k for k, v in prunned_high_reps.items() if len(v) == 0]

				excluded_reps.extend(reps_to_remove)

				# Determine smallest set of representatives that allow to get all cycle candidates
				excluded = []
				hotspot_reps = set(im.flatten_list(list(hot_reps.values())))
				for rep, hits in hot_reps.items():
					common = hotspot_reps.intersection(set(hits))
					if len(common) > 0:
						hotspot_reps = hotspot_reps - common
					else:
						excluded.append(rep)

				excluded_reps.extend(excluded)

				# Remove representatives that only led to low BSR
				excluded_reps.extend(low_reps)

				representatives = [rep for rep in representatives if rep not in excluded_reps]
				ids_to_blast = [i for i in ids_to_blast if i not in excluded_reps]

				# Determine next representative from candidates
				rep_candidates = list(set(hotspots) - hitting_high)
				# Sort to guarantee reproducible results with same datasets
				rep_candidates = sorted(rep_candidates, key=lambda x: x.split('_')[-1])
				representatives, final_representatives = select_candidate(rep_candidates,
																		  prot_seqs,
																		  ids_to_blast,
																		  representatives,
																		  final_representatives)

				# Remove files created for current gene iteration
				os.remove(rep_file)
				os.remove(blast_output)
				os.remove(ids_file)
		else:
			final_representatives = list(prot_seqs.keys())

		# Write schema file with all alleles
		# Add locus identifier as prefix to sequence IDs
		# Some external schemas, such as schemas from cgMLST.org, do not include the locus identifier in the sequence headers
		locus_data = [[k, v] if locus_id in k else [f'{locus_id}_{k}', v] for k, v in dna_seqs.items()]
		locus_lines = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, locus_data)
		fo.write_lines(locus_lines, locus_file)

		# Get total number of valid sequences
		valid_sequences = len(locus_lines)

		# Write schema file with representatives
		# Do not forget to add locus identifier as prefix to sequence IDs if it is not included in the original ID
		locus_rep_data = [[r, dna_seqs[r]] if locus_id in r else [f'{locus_id}_{r}', dna_seqs[r]] for r in final_representatives]
		locus_rep_lines = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE,
										  locus_rep_data)
		fo.write_lines(locus_rep_lines, locus_short_file)

		# Determine the number of representatives for current locus
		representatives_number = len(locus_rep_lines)

		summary_stats.append([locus_id,
							  str(total_sequences),
							  str(valid_sequences),
							  str(representatives_number)])

		shutil.rmtree(locus_temp_dir)

	return [invalid_alleles, invalid_loci, summary_stats]


def main(config):
	"""
	Adapt a schema to be used with chewBBACA.

	Parameters
	----------
	input_files : str
		Path to a TXT file with the list of schema loci to adapt.
	output_directories :  list
		Path to the output directories to create (the main schema
		directory and the 'short' directory to store representative
		alleles).
	cpu_cores : int
		Number of CPU cores that will be used to run the process.
	blast_score_ratio : float
		The BLAST Score Ratio value that will be used to evaluate
		allele similarity and select representative alleles.
	minimum_length : int
		Minimum sequence length value stored in the schema config file.
		The schema adaptation process will only discard sequences smaller
		that this value if the `--size-filter` parameter is provided.
	translation_table : int
		Genetic code used to translate alleles.
	size_threshold : float
		Allele size variation threshold value stored in the schema
		config file. The schema adaptation process will only discard
		alleles below or above the locus size threshold if the
		´--size-filter´ parameter is provided.
	blast_path : str
		Path to the directory that contains the BLAST executables.
	"""
	schema_path = fo.join_paths(config.output_directory, ["adapted_schema"])
	# Create the "short" directory to store the FASTA files with representative alleles
	schema_short_path = fo.join_paths(schema_path, ["short"])
	fo.create_directory(schema_short_path)

	logger.info(f"Number of loci to adapt: {len(config.loci_list[0])}")
	logger.info("Determining the total number of alleles and allele mean length per gene...")

	# Count number of sequences and mean length per locus
	loci_info = []
	loci_pools = multiprocessing.Pool(processes=config.cpu_cores)
	gp = loci_pools.map_async(fao.fasta_stats, config.loci_list[1], callback=loci_info.extend)
	gp.wait()

	# Split files according to number of sequences and sequence mean length
	# in each file to pass even groups of sequences to all cores
	# divide into 100 input sets to get 1% progress resolution
	even_loci_groups = mo.distribute_loci(loci_info, 100, "seqcount")
	# with few inputs, some sublists might be empty
	even_loci_groups = [i for i in even_loci_groups if len(i) > 0]
	# Add common arguments
	blastp_path = os.path.join(config.blast_path, gct.BLASTP_ALIAS)
	makeblastdb_path = os.path.join(config.blast_path, gct.MAKEBLASTDB_ALIAS)
	blastdb_aliastool_path = fo.join_paths(config.blast_path, [gct.BLASTDB_ALIASTOOL_ALIAS])
	even_loci_groups = [[i, schema_path, schema_short_path,
						 config.blast_score_ratio, config.adaptation_minimum_length,
						 config.translation_table, config.adaptation_size_threshold,
						 blastp_path, makeblastdb_path, blastdb_aliastool_path,
						 adapt_loci] for i in even_loci_groups]

	logger.info(f"Adapting {len(config.loci_list[0])} loci...")
	invalid_data = mo.map_async_parallelizer(even_loci_groups, mo.function_helper, config.cpu_cores)

	# Write files with list of invalid alleles and invalid loci
	schema_basename = fo.file_basename(schema_path.rstrip('/'))
	parent_directory = os.path.dirname(schema_path)

	# Write file with alleles that were determined to be invalid
	invalid_alleles = [sub[0] for sub in invalid_data]
	invalid_alleles = list(itertools.chain.from_iterable(invalid_alleles))
	invalid_alleles_file = os.path.join(parent_directory, f"{schema_basename}_invalid_alleles.txt")
	invalid_alleles = [f'{allele[0]}: {allele[1]}' for allele in invalid_alleles]
	fo.write_lines(invalid_alleles, invalid_alleles_file)

	# Write file with identifiers of loci that had no valid alleles
	invalid_loci = [sub[1] for sub in invalid_data]
	invalid_loci = list(itertools.chain.from_iterable(invalid_loci))
	invalid_loci_file = os.path.join(parent_directory,
									  "{0}_{1}".format(schema_basename, "invalid_loci.txt"))
	fo.write_lines(invalid_loci, invalid_loci_file)

	stats_lines = [sub[2] for sub in invalid_data]
	stats_lines = list(itertools.chain.from_iterable(stats_lines))
	stats_lines = ['\t'.join(line) for line in stats_lines]
	stats_loci_file = '{0}/{1}_{2}'.format(parent_directory, schema_basename, "summary_stats.tsv")
	stats_lines = [im.join_list(ct.ADAPTSCHEMA_STATS_HEADER, "\t")] + stats_lines
	fo.write_lines(stats_lines, stats_loci_file)

	logger.info(f"Number of invalid loci: {len(invalid_loci)}")
	logger.info(f"Wrote list of invalid loci to {invalid_loci_file}")

	logger.info(f"Number of invalid alleles: {len(invalid_alleles)}")
	logger.info(f"Wrote list of invalid alleles to {invalid_alleles_file}")
	logger.info(f"Wrote file with summary statistics to {stats_loci_file}")

	logger.info("Successfully adapted {0}/{1} loci present in the "
		  		"input schema.".format(len(config.loci_list[0])-len(invalid_loci), len(config.loci_list[0])))

	# Copy Pyrodigal Training File (PTF) to schema directory
	if config.gene_predictor == "pyrodigal" and config.validated_gene_prediction_options.pyrodigal_training_file is not None:
		fo.copy_file(config.validated_gene_prediction_options.pyrodigal_training_file, schema_path)
		logger.info(f"Copied Pyrodigal training file to {schema_path}")

	# Write schema config file
	schema_config = cf.write_schema_config(config, gct.__version__, schema_path)
	logger.info(f"Wrote schema config values to {schema_config[1]}")

	# Create the file with the list of genes/loci
	cf.write_gene_list(schema_path)
	logger.info(f"Wrote list of loci to {os.path.join(schema_path, gct.GENE_LIST_BASENAME)}")
