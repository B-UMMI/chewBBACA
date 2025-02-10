#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the creation of a whole genome multi locus sequence
typing (wgMLST) schema seed. The process selects one representative allele
per distinct locus identified in the input files. The schema seed corresponds
to a wgMLST schema with one FASTA file per distinct locus, each FASTA file
containing the representative allele selected by the process.

Code documentation
------------------
"""


import os
import sys
import math

try:
	from PredictCDSs import predict_cdss
	from utils import (constants as ct,
					   blast_wrapper as bw,
					   core_functions as cf,
					   file_operations as fo,
					   fasta_operations as fao,
					   sequence_manipulation as sm,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.PredictCDSs import predict_cdss
	from CHEWBBACA.utils import (constants as ct,
								 blast_wrapper as bw,
								 core_functions as cf,
								 file_operations as fo,
								 fasta_operations as fao,
								 sequence_manipulation as sm,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def create_schema_structure(schema_seed_fasta, output_directory, schema_name):
	"""Create the schema directory structure.

	Creates the schema seed directory with one FASTA file per
	distinct locus and the `short` directory with the FASTA files
	used to save the representative sequences.

	Parameters
	----------
	schema_seed_fasta : str
		Path to the FASTA file that contains the sequences that
		constitute the schema seed. Each FASTA record in the file
		is a representative sequence chosen for a locus.
	output_directory : str
		Path to the main output directory of the process.
	schema_name : str
		Name for the schema's directory.

	Returns
	-------
	schema_files : list
		List with the paths to the FASTA files in the schema seed.
	"""
	schema_dir = fo.join_paths(output_directory, [schema_name])
	fo.create_directory(schema_dir)

	# add allele identifier to all sequences
	schema_records = {im.replace_multiple_characters(rec.id, ct.CHAR_REPLACEMENTS): str(rec.seq)
					  for rec in fao.sequence_generator(schema_seed_fasta)}

	loci_basenames = {k: k+'.fasta' for k in schema_records}
	loci_paths = {k: fo.join_paths(schema_dir, [v])
				  for k, v in loci_basenames.items()}

	for k, v in schema_records.items():
		current_representative = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [k+'_1', v])
		fo.write_to_file(current_representative, loci_paths[k], 'w', '\n')

	# create 'short' directory
	fo.create_short(loci_paths.values(), schema_dir)

	return loci_paths


def create_schema_seed(fasta_files, output_directory, schema_name, ptf_path,
					   blast_score_ratio, minimum_length, translation_table,
					   size_threshold, word_size, window_size, clustering_sim,
					   representative_filter, intra_filter, cpu_cores, blast_path,
					   prodigal_mode, cds_input):
	"""Create a schema seed based on a set of input FASTA files."""
	# Map full paths to unique identifier (prefix before first '.')
	full_to_basename = im.mapping_function(fasta_files, fo.file_basename, [False])
	full_to_unique = {k: fo.split_joiner(v, [0], '.')
						for k, v in full_to_basename.items()}

	# Create directory to store temporary files
	temp_directory = fo.join_paths(output_directory, ['temp'])
	fo.create_directory(temp_directory)

	# Create directory to store files with Pyrodigal results
	pyrodigal_path = fo.join_paths(temp_directory, ['1_cds_prediction'])
	fo.create_directory(pyrodigal_path)
	if cds_input is False:
		# Run Pyrodigal to predict genes for all input genomes
		print(f'\n {ct.CDS_PREDICTION} ')
		print('='*(len(ct.CDS_PREDICTION)+2))

		# Gene prediction step
		pyrodigal_results = predict_cdss.main(fasta_files, pyrodigal_path, ptf_path,
											  translation_table, prodigal_mode, None,
											  None, ['genes'], None, cpu_cores)

		# Dictionary with info about inputs for which gene prediction failed
		# Total number of CDSs identified in the inputs
		# Paths to FASTA files with the extracted CDSs
		# Paths to files with the coordinates of the CDSs extracted for each input
		# Total number of CDSs identified per input
		# Dictionary with info about the CDSs closer to contig tips per input
		failed, total_extracted, cds_fastas, cds_coordinates, cds_counts, _ = pyrodigal_results
		if len(failed) > 0:
			print(f'\nFailed to predict genes for {len(failed)} inputs')
			print('Make sure that Prodigal runs in meta mode (--pm meta) '
				  'if any input file has less than 100kbp.')
		if len(cds_fastas) == 0:
			sys.exit(f'\n{ct.CANNOT_PREDICT}')

		renaming_inputs = []
		for i, file in enumerate(cds_fastas):
			basename = fo.file_basename(file, False)
			parent_dir = os.path.dirname(os.path.dirname(file))
			output_file = fo.join_paths(parent_dir, [f'{basename}.fasta'])
			cds_prefix = f'{basename}-protein'
			renaming_inputs.append([file, output_file, 1, 50000,
									cds_prefix, False, True, fao.integer_headers])
			cds_fastas[i] = output_file

		# Rename CDSs in files
		renaming_results = mo.map_async_parallelizer(renaming_inputs,
													 mo.function_helper,
													 cpu_cores,
													 show_progress=False)

		print(f'\nExtracted a total of {total_extracted} CDSs from {len(fasta_files)} inputs.')
	# Inputs are Fasta files with the predicted CDSs
	else:
		# Rename the CDSs in each file based on the input unique identifiers
		print(f'\nRenaming CDSs for {len(full_to_unique)} input files...')
		# Rename CDSs in each file to ensure IDs are compatible with chewBBACA
		renaming_inputs = []
		cds_fastas = []
		for k, v in full_to_unique.items():
			output_file = fo.join_paths(pyrodigal_path, [f'{v}.fasta'])
			cds_prefix = f'{v}-protein'
			renaming_inputs.append([k, output_file, 1, 50000,
									cds_prefix, False, True, fao.integer_headers])
			cds_fastas.append(output_file)

		# Rename CDSs in files
		renaming_results = mo.map_async_parallelizer(renaming_inputs,
														mo.function_helper,
														cpu_cores,
														show_progress=False)

		# No inputs failed gene prediction
		failed = []
		# Cannot get CDS coordinates if skipping gene prediction
		cds_coordinates = None
		total_cdss = sum([r[1] for r in renaming_results])
		print(f'Input files contain a total of {total_cdss} coding sequences.')

	if len(failed) > 0:
		# Exclude inputs that failed gene prediction
		full_to_unique = im.prune_dictionary(full_to_unique, failed.keys())
		# Write Prodigal stderr for inputs that failed gene prediction
		failed_lines = [f'{k}\t{v}' for k, v in failed.items()]
		failed_outfile = fo.join_paths(output_directory,
									   ['gene_prediction_failures.tsv'])
		fo.write_lines(failed_lines, failed_outfile)

	# Map input identifiers to integers
	# Use the mapped integers to refer to each input
	# This reduces memory usage compared to using string identifiers
	unique_to_int = im.integer_mapping(full_to_unique.values())
	int_to_unique = im.invert_dictionary(unique_to_int)

	# Concatenate subgroups of FASTA files before deduplication
	num_chunks = 20 if cpu_cores <= 20 else cpu_cores
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
	# create directory to store files from DNA deduplication
	dna_dedup_dir = fo.join_paths(preprocess_dir, ['cds_deduplication'])
	fo.create_directory(dna_dedup_dir)
	print('Identifying distinct CDSs...')
	dna_dedup_results = cf.exclude_duplicates(cds_files, dna_dedup_dir, cpu_cores,
									   [unique_to_int, int_to_unique],
									   False, True)

	_, distinct_seqids, distinct_file, repeated = dna_dedup_results
	print(f'Identified {len(distinct_seqids)} distinct CDSs.')

	# Delete concatenated FASTA files
	fo.remove_files(cds_files)

	# Index FASTA file with distinct DNA sequences
	dna_index = fao.index_fasta(distinct_file)

	# Translate CDSs
	print(f'\n {ct.CDS_TRANSLATION} ')
	print('='*(len(ct.CDS_TRANSLATION)+2))

	# Create directory to store translation results
	cds_translation_dir = fo.join_paths(preprocess_dir, ['cds_translation'])
	fo.create_directory(cds_translation_dir)
	print(f'Translating {len(distinct_seqids)} CDS...')
	# This step excludes small sequences
	ts_results = cf.translate_sequences(distinct_seqids, distinct_file,
										cds_translation_dir,
										translation_table,
										minimum_length,
										cpu_cores)

	protein_file, ut_seqids, ut_lines = ts_results
	print(f'\n{len(ut_seqids)} CDSs could not be translated.')

	# Create file with list of invalid CDSs
	invalid_file = fo.join_paths(output_directory, [ct.INVALID_CDS_BASENAME])
	invalid_alleles = im.join_list(im.sort_iterable(ut_lines), '\n')
	fo.write_to_file(invalid_alleles, invalid_file, 'w', '\n')

	# Protein deduplication
	print(f'\n {ct.PROTEIN_DEDUPLICATION} ')
	print('='*(len(ct.PROTEIN_DEDUPLICATION)+2))
	# Create directory to store files from protein deduplication
	protein_dedup_dir = fo.join_paths(preprocess_dir, ['protein_deduplication'])
	fo.create_directory(protein_dedup_dir)
	print('Identifying distinct proteins...')
	ds_results = cf.exclude_duplicates([protein_file], protein_dedup_dir, 1,
									   [unique_to_int, int_to_unique],
									   True, True)

	_, representative_pseqids, representative_pfasta, _ = ds_results

	print(f'Identified {len(representative_pseqids)} distinct proteins.')

	# Sort remaining seqids
	schema_seqids = im.sort_iterable(representative_pseqids, sort_key=lambda x: x.lower())
	print(f'Kept {len(schema_seqids)} sequences after filtering the initial sequences.')

	# Protein clustering
	print(f'\n {ct.PROTEIN_CLUSTERING} ')
	print('='*(len(ct.PROTEIN_CLUSTERING)+2))

	# Create directory to store clustering data
	clustering_dir = fo.join_paths(temp_directory, ['3_clustering'])
	fo.create_directory(clustering_dir)

	# Read protein sequences
	proteins = fao.import_sequences(representative_pfasta)

	print('Clustering proteins...')
	# Do not divide based on the number of available cores
	# It can lead to different results depending on the number of cores
	# Divide into fixed number of groups
	group_size = math.ceil(len(proteins)/ct.CREATESCHEMA_CLUSTERING_NGROUPS)
	cs_results = cf.cluster_sequences(proteins, word_size, window_size,
									  clustering_sim, None, True,
									  1, 1, clustering_dir, cpu_cores,
									  group_size, False)
	print(f'\nClustered {len(proteins)} proteins into {len(cs_results)} clusters.')

	# Exclude based on high similarity to cluster representatives
	rep_filter_dir = fo.join_paths(clustering_dir, ['representative_filter'])
	fo.create_directory(rep_filter_dir)
	print('Removing proteins highly similar to the cluster representative...')
	cp_results = cf.cluster_representative_filter(cs_results,
												  representative_filter,
												  rep_filter_dir)
	clusters, excluded_seqids, singletons, clustered_sequences = cp_results
	print(f'Removed {len(excluded_seqids)} sequences.')
	print(f'Identified {len(singletons)} singletons.')
	print(f'Remaining sequences after representative and singleton pruning: {clustered_sequences}')

	# Remove excluded seqids
	schema_seqids = list(set(schema_seqids) - excluded_seqids)

	# Exclude based on high similarity to other clustered sequences
	intra_filter_dir = fo.join_paths(clustering_dir, ['intracluster_filter'])
	fo.create_directory(intra_filter_dir)
	print('Removing sequences highly similar to other clustered sequences...')
	cip_results = cf.cluster_intra_filter(clusters, proteins,
										  word_size, intra_filter,
										  intra_filter_dir)
	clusters, intra_excluded = cip_results
	print(f'Removed {len(intra_excluded)} sequences.')

	# Remove excluded seqids - we get set of sequences from clusters
	# plus singletons
	schema_seqids = list(set(schema_seqids) - set(intra_excluded))

	# Define BLASTp and makeblastdb paths
	blastp_path = fo.join_paths(blast_path, [ct.BLASTP_ALIAS])
	makeblastdb_path = fo.join_paths(blast_path, [ct.MAKEBLASTDB_ALIAS])
	blastdb_aliastool_path = fo.join_paths(blast_path, [ct.BLASTDB_ALIASTOOL_ALIAS])

	# All-vs-all BLASTp per cluster
	if len(clusters) > 0:
		print(f'Clusters to BLAST: {len(clusters)}')
		print('Performing all-vs-all BLASTp per cluster...')
		blast_results, blast_results_dir = cf.blast_clusters(clusters, representative_pfasta,
													clustering_dir, blastp_path,
													makeblastdb_path, cpu_cores,
													blastdb_aliastool_path)

		blast_files = im.flatten_list(blast_results)

		# Compute and exclude based on BSR
		print('\nRemoving sequences based on high BSR...')
		bsr_excluded = [sm.apply_bsr(fo.read_tabular(file),
									 dna_index,
									 blast_score_ratio)
						for file in blast_files]

		# Merge BSR results
		bsr_excluded = set(im.flatten_list(bsr_excluded))
		schema_seqids = list(set(schema_seqids) - bsr_excluded)
		print(f'Removed {len(bsr_excluded)} sequences.')

		# Write list of excluded to file
		blast_excluded_outfile = fo.join_paths(blast_results_dir, ['bsr_excluded.txt'])
		fo.write_lines(bsr_excluded, blast_excluded_outfile)

	# Perform final BLAST to identify similar sequences that do not share many/any kmers
	print(f'\n {ct.FINAL_BLASTp} ')
	print('='*(len(ct.FINAL_BLASTp)+2))
	print(f'{len(schema_seqids)} sequences to compare in final BLASTp.'.format(len(schema_seqids)))

	# Sort seqids before final BLASTp to ensure consistent results
	schema_seqids = im.sort_iterable(schema_seqids, sort_key=lambda x: x.lower())

	# Create directory for final BLASTp
	final_blast_dir = fo.join_paths(temp_directory, ['4_final_blast'])
	fo.create_directory(final_blast_dir)

	# Create FASTA file with remaining sequences
	quasi_schema_file = os.path.join(final_blast_dir, 'remaining_sequences.fasta')
	fao.get_sequences_by_id(proteins, schema_seqids, quasi_schema_file)

	# Create BLASTp database
	blast_db = fo.join_paths(final_blast_dir, ['remaining_sequences'])
	db_std = bw.make_blast_db(makeblastdb_path, quasi_schema_file, blast_db, 'prot')

	# Divide FASTA file into groups of 100 sequences to reduce
	# execution time for large sequence sets
	split_dir = fo.join_paths(final_blast_dir, ['cds_subsets'])
	fo.create_directory(split_dir)
	splitted_fastas = fao.split_seqcount(quasi_schema_file, split_dir, 100)

	# Create directory to store results from final BLASTp
	final_blastp_dir = fo.join_paths(final_blast_dir, ['BLAST_results'])
	fo.create_directory(final_blastp_dir)
	blast_outputs = ['{0}/{1}_blast_out.tsv'.format(final_blastp_dir,
													fo.file_basename(i[0], False))
					 for i in splitted_fastas]

	# Add common arguments to all sublists
	blast_inputs = [[blastp_path, blast_db, file[0],
					 blast_outputs[i], 1, 1, bw.run_blast]
					for i, file in enumerate(splitted_fastas)]

	print('Performing final BLASTp...')
	blast_results = mo.map_async_parallelizer(blast_inputs,
											  mo.function_helper,
											  cpu_cores,
											  show_progress=True)

	# Concatenate files with BLASTp results
	blast_output = fo.join_paths(final_blast_dir, ['blast_out_concat.tsv'])
	blast_output = fo.concatenate_files(blast_outputs, blast_output)
	final_excluded = sm.apply_bsr(fo.read_tabular(blast_output),
								  dna_index,
								  blast_score_ratio)
	schema_seqids = list(set(schema_seqids) - set(final_excluded))
	print(f'\nRemoved {len(final_excluded)} sequences highly similar to other sequences.')

	# Create file with the schema representative sequences
	loci_representatives = os.path.join(final_blast_dir, 'loci_representatives.fasta')
	fao.get_sequences_by_id(dna_index, schema_seqids, loci_representatives)

	# Create schema directory and FASTA files
	print(f'Creating schema seed in {output_directory}')
	schema_files = create_schema_structure(loci_representatives, output_directory,
										   schema_name)

	# Copy file with CDS coordinates to output directory
	# Will not be created if input files contain predicted CDSs
	if cds_input is False:
		fo.copy_file(cds_coordinates[0], output_directory)

		# Create TSV with allele ID to CDS ID
		clines = []
		for seqid in schema_seqids:
			allele_id = f'{im. replace_multiple_characters(seqid, ct.CHAR_REPLACEMENTS)}_1'
			clines.append(f'{allele_id}\t{seqid}\tY')

		coutfile = fo.join_paths(output_directory, ['selected_ids.tsv'])
		clines = [ct.ALLELE_TO_CDS_HEADER] + clines
		fo.write_lines(clines, coutfile)

	return [schema_files, temp_directory]


def main(input_files, output_directory, schema_name, ptf_path,
		 blast_score_ratio, minimum_length, translation_table,
		 size_threshold, word_size, window_size, clustering_sim,
		 representative_filter, intra_filter, cpu_cores, blast_path,
		 cds_input, prodigal_mode, no_cleanup):
	"""Create a wgMLST schema seed.

	Parameters
	----------
	input_files : str
		Path to the directory that contains the input FASTA files.
		Alternatively, a single file with a list of paths to FASTA
		files, one per line.
	output_directory : str
		Output directory where the process will store intermediate
		files and create the schema seed.
	schema_name : str
		Name given to the folder that will store the schema seed files.
	ptf_path : str
		Path to the Prodigal training file.
	blast_score_ratio : float
		BLAST Score Ratio value.
	minimum_length : int
		Minimum sequence length. Coding sequences shorter than this
		value are excluded.
	translation_table : int
		Genetic code used to predict genes and to translate coding
		sequences.
	size_threshold : float
		CDS size variation threshold. Added to the schema's config
		file and used to identify alleles with a length value that
		deviates from the locus length mode during the allele calling
		process.
	word_size : int
		K-mer size used during minimizer clustering.
	window_size : int
		Number of consecutive k-mers included in each window to
		determine a minimizer.
	clustering_sim :float
		Minimum decimal proportion of shared distinct minimizers for
		a sequence to be added to a cluster.
	representative_filter : float
		Clustered sequences are excluded if they share this proportion
		of distinct minimizers with the cluster representative.
	intra_filter : float
		Clustered sequences are excluded if they share this proportion
		of distinct minimizers with another clustered sequence of equal
		or greater length.
	cpu_cores : int
		Number of CPU cores used to run the process.
	blast_path : str
		Path to the BLAST executables.
	cds_input : bool
		If provided, input is a single or several FASTA files with
		coding sequences (skips gene prediction and CDS extraction).
	prodigal_mode : str
		Prodigal running mode ("single" or "meta").
	no_cleanup : bool
		If provided, intermediate files generated during process
		execution are not removed at the end.
	"""
	# Read file with paths to input files
	input_files = fo.read_lines(input_files, strip=True)

	# Sort paths to FASTA files
	input_files = im.sort_iterable(input_files, sort_key=lambda x: x.lower())

	results = create_schema_seed(input_files, output_directory, schema_name,
								 ptf_path, blast_score_ratio, minimum_length,
								 translation_table, size_threshold, word_size,
								 window_size, clustering_sim, representative_filter,
								 intra_filter, cpu_cores, blast_path,
								 prodigal_mode, cds_input)

	# Remove temporary files
	if no_cleanup is False:
		exists = fo.delete_directory(results[1])

	return len(results[0])
