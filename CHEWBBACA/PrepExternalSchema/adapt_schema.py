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

try:
	from utils import (constants as ct,
					   blast_wrapper as bw,
					   file_operations as fo,
					   fasta_operations as fao,
					   sequence_manipulation as sm,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 blast_wrapper as bw,
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

		# Dictionaries mapping gene identifiers to DNA and Protein sequences
		locus_seqs, prot_seqs, locus_invalid, seqids_map, total_sequences = \
			sm.get_seqs_dicts(locus, locus_id, table_id, min_len, size_threshold)
		invalid_alleles.extend(locus_invalid)

		# Continue to next locus if there are no valid CDSs for current locus
		if len(prot_seqs) == 0:
			shutil.rmtree(locus_temp_dir)
			invalid_loci.append(locus_id)
			summary_stats.append([locus_id, str(total_sequences), '0', '0'])
			continue

		if len(locus_seqs) > 1:
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
			protein_lines = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, protein_data)
			fo.write_lines(protein_lines, protein_file)

			# Create blastdb with all distinct proteins
			blastp_db = os.path.join(locus_temp_dir, locus_id)
			db_std = bw.make_blast_db(makeblastdb_path, protein_file, blastp_db, 'prot')

			# Determine appropriate blastp task (proteins < 30aa need blastp-short)
			blastp_task = bw.determine_blast_task(equal_prots)

			# Cycle BLAST representatives against non-representatives until
			# all non-representatives have a representative
			while len(set(ids_to_blast) - set(representatives)) != 0:
				# create FASTA file with representative sequences
				rep_file = fo.join_paths(locus_temp_dir,
										 [f'{locus_id}_rep_protein.fasta'])
				rep_protein_data = [[r, prot_seqs[r]] for r in representatives]
				rep_protein_lines = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, rep_protein_data)
				fo.write_lines(rep_protein_lines, rep_file)

				# Compute self-score for representative alleles
				for seqid in representatives:
					if seqid not in rep_self_scores:
						record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [seqid, prot_seqs[seqid]])
						current_rep_file = fo.join_paths(locus_temp_dir, [f'{seqid}_solo.fasta'])
						fo.write_lines([record], current_rep_file)
						# Create file with representative seqid to only compare against self
						id_file = fo.join_paths(locus_temp_dir, [f'{seqid}_ids.txt'])
						fo.write_lines([seqid], id_file)
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
						if len(rep_results) > 0:
							rep_self_scores[rep_results[0][0]] = float(rep_results[0][6])
						else:
							print('Could not determine the self-alignment raw '
								f'score for {rep_results[0][0]}')

				# Create file with seqids to BLAST against
				ids_str = im.join_list([str(i) for i in ids_to_blast if i not in representatives], '\n')
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
				rep_candidates = sorted(rep_candidates, key=lambda x: int(x))
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
		locus_data = [[k, v] for k, v in locus_seqs.items()]
		locus_lines = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, locus_data)
		fo.write_lines(locus_lines, locus_file)

		# Get total number of valid sequences
		valid_sequences = len(locus_lines)

		# Write schema file with representatives
		final_representatives = [seqids_map[rep]
								 for rep in final_representatives]
		locus_rep_data = [[r, locus_seqs[r]]
						  for r in final_representatives]
		locus_rep_lines = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE,
										  locus_rep_data)
		fo.write_lines(locus_rep_lines, locus_short_file)

		# get number of representatives
		representatives_number = len(locus_rep_lines)

		summary_stats.append([locus_id,
							  str(total_sequences),
							  str(valid_sequences),
							  str(representatives_number)])

		shutil.rmtree(locus_temp_dir)

	return [invalid_alleles, invalid_loci, summary_stats]


def main(input_files, output_directories, cpu_cores, blast_score_ratio,
		 minimum_length, translation_table, size_threshold, blast_path):
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
	schema_path, schema_short_path = output_directories

	# Import list of loci to adapt
	genes_list = fo.read_lines(input_files, strip=True)

	print('Number of loci to adapt: {0}\n'.format(len(genes_list)))
	print('Determining the total number of alleles and '
		  'allele mean length per gene...\n'.format())

	# Count number of sequences and mean length per locus
	loci_info = []
	loci_pools = multiprocessing.Pool(processes=cpu_cores)
	gp = loci_pools.map_async(fao.fasta_stats, genes_list,
							  callback=loci_info.extend)
	gp.wait()

	# split files according to number of sequences and sequence mean length
	# in each file to pass even groups of sequences to all cores
	# divide into 100 input sets to get 1% progress resolution
	even_loci_groups = mo.distribute_loci(loci_info, 100, 'seqcount')
	# with few inputs, some sublists might be empty
	even_loci_groups = [i for i in even_loci_groups if len(i) > 0]
	# Add common arguments
	blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
	makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)
	blastdb_aliastool_path = fo.join_paths(blast_path, [ct.BLASTDB_ALIASTOOL_ALIAS])
	even_loci_groups = [[i, schema_path, schema_short_path,
						 blast_score_ratio, minimum_length,
						 translation_table, size_threshold,
						 blastp_path, makeblastdb_path, blastdb_aliastool_path,
						 adapt_loci] for i in even_loci_groups]

	print('Adapting {0} loci...\n'.format(len(genes_list)))
	invalid_data = mo.map_async_parallelizer(even_loci_groups,
											 mo.function_helper,
											 cpu_cores,
											 show_progress=True)

	# Write files with list of invalid alleles and invalid loci
	schema_basename = fo.file_basename(schema_path.rstrip('/'))
	parent_directory = os.path.dirname(schema_path)

	# write file with alleles that were determined to be invalid
	invalid_alleles = [sub[0] for sub in invalid_data]
	invalid_alleles = list(itertools.chain.from_iterable(invalid_alleles))
	invalid_alleles_file = os.path.join(parent_directory,
										f'{schema_basename}_invalid_alleles.txt')
	invalid_alleles = [f'{allele[0]}: {allele[1]}' for allele in invalid_alleles]
	fo.write_lines(invalid_alleles, invalid_alleles_file)

	# Write file with identifiers of loci that had no valid alleles
	invalid_loci = [sub[1] for sub in invalid_data]
	invalid_loci = list(itertools.chain.from_iterable(invalid_loci))
	invalid_loci_file = os.path.join(parent_directory,
									  '{0}_{1}'.format(schema_basename, 'invalid_loci.txt'))
	fo.write_lines(invalid_loci, invalid_loci_file)

	stats_lines = [sub[2] for sub in invalid_data]
	stats_lines = list(itertools.chain.from_iterable(stats_lines))
	stats_lines = ['\t'.join(line) for line in stats_lines]
	stats_loci_file = '{0}/{1}_{2}'.format(parent_directory,
										   schema_basename,
										   'summary_stats.tsv')
	stats_lines = [ct.PREPEXTERNAL_SUMMARY_STATS_HEADER] + stats_lines
	fo.write_lines(stats_lines, stats_loci_file)

	print(f'\n\nNumber of invalid loci: {len(invalid_loci)}')
	print(f'Wrote list of invalid loci to {invalid_loci_file}')

	print(f'Number of invalid alleles: {len(invalid_alleles)}')
	print(f'Wrote list of invalid alleles to {invalid_alleles_file}')
	print(f'Wrote file with summary statistics to {stats_loci_file}')

	print('\nSuccessfully adapted {0}/{1} loci present in the '
		  'input schema.'.format(len(genes_list)-len(invalid_loci),
								 len(genes_list)))
