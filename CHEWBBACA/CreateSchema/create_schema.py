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
import math

from loguru import logger

from CHEWBBACA import global_constants as gct
from CHEWBBACA import global_validation as gval
from CHEWBBACA.CreateSchema import constants as ct
from CHEWBBACA.PredictGenes import predict_genes
from CHEWBBACA.utils import (blast_wrapper as bw,
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
	# Create folder to store schema files
	schema_directory = fo.join_paths(output_directory, [schema_name])
	fo.create_directory(schema_directory)

	# Import schema representative alleles
	schema_records = fao.import_sequences(schema_seed_fasta)
	# Convert the CDS identifiers into allele identifiers (add "-protein")
	# Sequence IDs may include multiple "-", so always split only based on the last one
	schema_records = {f"{seqid.rsplit('_', 1)[0]}-protein{seqid.rsplit('_', 1)[1]}": seq for seqid, seq in schema_records.items()}
	# Replace special chars
	schema_records = {im.replace_multiple_characters(seqid, gct.CHAR_REPLACEMENTS): seq for seqid, seq in schema_records.items()}

	# Define basenames and paths for loci FASTA files (the file basename is used as locus ID)
	loci_paths = {seqid: fo.join_paths(schema_directory, [seqid+'.fasta']) for seqid in schema_records}
	# Create loci FASTA files
	for seqid, seq in schema_records.items():
		# Representative is the first allele, so it is necessary to add "_1" to the sequence ID
		current_representative = fao.fasta_str_record(gct.FASTA_RECORD_TEMPLATE, [seqid+'_1', seq])
		fo.write_to_file(current_representative, loci_paths[seqid])

	# Create 'short' directory with the FASTA files that contain only the representative alleles
	fo.create_short(loci_paths.values(), schema_directory)

	return loci_paths


def create_schema_seed(config):
	"""Create a schema seed based on a set of input FASTA files."""
	# Map input file paths to file basename without extension and MD5 file hash
	input_file_ids = [(file, fo.file_basename(file, False)) for file in config.input_path]

	# Create directory to store temporary files
	temp_directory = fo.join_paths(config.output_directory, ["temp"])
	fo.create_directory(temp_directory)

	# Create directory to store files with Pyrodigal results
	if config.cds_input is False:
		# Run Pyrodigal to predict genes for all input genomes
		# Gene prediction step
		logger.info(f"Calling the PredictGenes module to predict genes for {len(config.input_path)} inputs...")
		pyrodigal_results = predict_genes.main(config)

		# Dictionary with info about inputs for which gene prediction failed
		# Total number of CDSs identified in the inputs
		# Paths to FASTA files with the extracted CDSs
		# Paths to files with the coordinates of the CDSs extracted for each input
		# Total number of CDSs identified per input
		# Dictionary with info about the CDSs closer to contig tips per input
		failed, _, cds_fastas, _, _, cds_counts, _ = pyrodigal_results
		if len(failed) > 0:
			logger.warning(f"Failed to predict genes for {len(failed)} inputs.")
			logger.warning("Make sure that Prodigal runs in meta mode if any input file has less than 100kbp.")
		if len(cds_fastas) == 0:
			logger.warning(f"{gct.PREDICTGENES_NO_PREDICTIONS}")
			raise gval.InputsNoPredictions
	# Inputs are FASTA files with the predicted CDSs

	##### With the new sequence deduplication approach, maybe it is not necessary to rename the sequences?
	##### Or maybe it is still necessary to ensure the consistency of locus naming
	else:
		if not config.no_cds_renaming:
			# Rename the input CDSs based on the unique identifier for each input file
			logger.info(f"Renaming CDSs for {len(input_file_ids)} input files...")
			# Create a temporary directory to store renamed CDSs
			renamed_dir = fo.join_paths(config.output_directory, ["renamed_genes"])
			fo.create_directory(renamed_dir)
			renaming_inputs = []
			cds_fastas = []
			for file in input_file_ids:
				output_file = fo.join_paths(renamed_dir, [f"{file[1]}.fasta"])
				# Use the basename of each input file as unique identifier
				cds_prefix = f"{file[1]}_"
				renaming_inputs.append([file[0], output_file, 1, 50000,
										cds_prefix, False, False, fao.integer_headers])
				cds_fastas.append(output_file)

			# Rename CDSs in files
			renaming_results = mo.map_async_parallelizer(renaming_inputs,
														 mo.function_helper,
														 config.cpu_cores)
		else:
			logger.info(f"User provided `--no-cds-renaming`. Will not rename CDS identifiers in input files. "
		 		  		"The process may not complete successfully if the CDS identifiers do not conform to the "
				  		"format used by chewBBACA.")

			cds_fastas = [file[0] for file in input_file_ids]

			# Try to get the TSV file with the CDS coordinates if the CDS sets were created by the PredictGenes module
			input_files_parent_dir = input_file_ids[0][0].rsplit('/', 2)[0]
			coordinates_file = fo.join_paths(input_files_parent_dir, [gct.GENE_COORDINATES_BASENAME])
			if os.path.isfile(coordinates_file):
				logger.info(f"Found a TSV file with CDS coordinate data in {coordinates_file}. Copying it to the output directory.")
				# Copy file with CDS coordinates to output directory
				destination = fo.join_paths(os.path.dirname(temp_directory), [gct.GENE_COORDINATES_BASENAME])
				fo.copy_file(coordinates_file, destination)

			# Try to get the TSV file with assembly statistics if the CDS sets were created by the PredictGenes module 
			assembly_stats_file = fo.join_paths(input_files_parent_dir, [gct.ASSEMBLY_STATS_BASENAME])
			if os.path.isfile(assembly_stats_file):
				logger.info(f"Found a TSV file with assembly statistics in {assembly_stats_file}. Copying it to the output directory.")
				destination = fo.join_paths(os.path.dirname(temp_directory), [gct.ASSEMBLY_STATS_BASENAME])
				fo.copy_file(assembly_stats_file, destination)

		# No inputs failed gene prediction
		failed = []
		# Need to count the number of CDSs in each input file to get total
		cds_counts = {fo.file_basename(f, False): fao.count_sequences(f) for f in cds_fastas}

	total_cdss = sum(cds_counts.values())
	logger.info(f"Total number of CDSs identified from all input files: {total_cdss}")

	if len(failed) > 0:
		# Exclude inputs that failed gene prediction
		input_file_ids = [file for file in input_file_ids if file[0] not in failed]
		# Write Prodigal stderr for inputs that failed gene prediction
		failed_lines = [f"{k}\t{v}" for k, v in failed.items()]
		failed_outfile = fo.join_paths(config.output_directory, [gct.GENE_PREDICTION_FAILURES_BASENAME])
		fo.write_lines(failed_lines, failed_outfile)

	# Create directory to store files from pre-process steps
	preprocess_dir = fo.join_paths(temp_directory, ["1_cds_preprocess"])
	fo.create_directory(preprocess_dir)

#############################

	distinct_file = fo.join_paths(preprocess_dir, ["distinct.fasta"])
	distinct_seqids = sm.deduplicate_predicted_genes(cds_fastas, distinct_file, config.cpu_cores)
	#fo.pickle_dumper(distinct_seqids, fo.join_paths(preprocess_dir, ["distinct"]))

	# Only need the seqids
	distinct_seqids = [v[0] for k, v in distinct_seqids.items()]

##############################

	logger.info(f"Identified {len(distinct_seqids)} distinct CDSs and excluded {total_cdss-len(distinct_seqids)} duplicated CDSs.")

	# Index FASTA file with distinct DNA sequences
	dna_index = fao.index_fasta(distinct_file)

	# Create directory to store translation results
	cds_translation_dir = fo.join_paths(preprocess_dir, ["cds_translation"])
	fo.create_directory(cds_translation_dir)
	logger.info(f"Translating {len(distinct_seqids)} CDS...")
	# This step excludes small sequences
	ts_results = cf.translate_sequences(distinct_seqids, distinct_file,
										cds_translation_dir,
										config.translation_table,
										config.minimum_length,
										config.cpu_cores)

	protein_file, ut_seqids, ut_lines = ts_results
	logger.info(f"{len(ut_seqids)} CDSs could not be translated.")

	# Create file with list of invalid distinct CDSs
	invalid_file = fo.join_paths(config.output_directory, [gct.INVALID_CDS_BASENAME])
	invalid_alleles = im.join_list(im.sort_iterable(ut_lines), "\n")
	fo.write_to_file(invalid_alleles, invalid_file, "w", "\n")
##### Should count the number of invalid CDSs per input like in AlleleCall

	# Protein deduplication
	# Create directory to store files from protein deduplication
	protein_dedup_dir = fo.join_paths(preprocess_dir, ["protein_deduplication"])
	fo.create_directory(protein_dedup_dir)
	logger.info("Identifying distinct proteins...")

### Should scan the FASTA file with protein sequences to divide it in a number of intervals equal to the number of CPU cores to process in parallel

#################### This works for schema creation, but what we need is not a BitMap with the list of genome IDs but a list with the seqids for all distinct DNA sequences encoding each protein
#################### Instead of STR seqids I can probably store a list of BLAKE2b hashes to redice memory usage a little

	representative_pfasta = fo.join_paths(protein_dedup_dir, ["distinct_proteins.fasta"])
	representative_pseqids = sm.deduplicate_predicted_genes([protein_file], representative_pfasta, config.cpu_cores)

	# Only need the seqids
	representative_pseqids = [v[0] for k, v in representative_pseqids.items()]

####################

	logger.info(f"Identified {len(representative_pseqids)} distinct proteins and excluded {len(distinct_seqids)-len(representative_pseqids)} duplicated proteins.")

	# Sort remaining seqids
	schema_seqids = im.sort_iterable(representative_pseqids, sort_key=lambda x: x.lower())
	logger.info(f"Kept {len(schema_seqids)} sequences after filtering the initial sequences.")

	# Protein clustering
	# Create directory to store clustering data
	clustering_dir = fo.join_paths(temp_directory, ["2_clustering"])
	fo.create_directory(clustering_dir)

	# Read protein sequences
	proteins = fao.import_sequences(representative_pfasta)

	logger.info("Clustering proteins...")
	# Do not divide based on the number of available cores
	# It can lead to different results depending on the number of cores
	# Divide into fixed number of groups
	group_size = math.ceil(len(proteins)/ct.CREATESCHEMA_CLUSTERING_NGROUPS)
	cs_results = cf.cluster_sequences(proteins,
	                                  config.validated_clustering_options.word_size,
	                                  config.validated_clustering_options.window_size,
									  config.validated_clustering_options.clustering_similarity,
									  None,
									  True,
									  1,
									  1,
									  clustering_dir,
									  config.cpu_cores,
									  group_size,
									  False)
	logger.info(f"Clustered {len(proteins)} proteins into {len(cs_results)} clusters.")

	# Exclude based on high similarity to cluster representatives
	rep_filter_dir = fo.join_paths(clustering_dir, ["representative_filter"])
	fo.create_directory(rep_filter_dir)
	logger.info("Removing proteins highly similar to the cluster representative...")
	cp_results = cf.cluster_representative_filter(cs_results, config.validated_clustering_options.representative_filter, rep_filter_dir)
	clusters, excluded_seqids, singletons, clustered_sequences = cp_results
	logger.info(f"Removed {len(excluded_seqids)} sequences.")
	logger.info(f"Identified {len(singletons)} singletons.")
	logger.info(f"Remaining sequences after representative and singleton pruning: {clustered_sequences}")

	# Remove excluded seqids
	schema_seqids = list(set(schema_seqids) - excluded_seqids)

	# Exclude based on high similarity to other clustered sequences
	intra_filter_dir = fo.join_paths(clustering_dir, ["intracluster_filter"])
	fo.create_directory(intra_filter_dir)
	logger.info("Removing sequences highly similar to other clustered sequences...")
	cip_results = cf.cluster_intra_filter(clusters, proteins,
										  config.validated_clustering_options.word_size,
										  config.validated_clustering_options.intra_filter,
										  intra_filter_dir)
	clusters, intra_excluded = cip_results
	logger.info(f"Removed {len(intra_excluded)} sequences.")

	# Remove excluded seqids - we get set of sequences from clusters
	# plus singletons
	schema_seqids = list(set(schema_seqids) - set(intra_excluded))

	# Define BLASTp and makeblastdb paths
	blastp_path = fo.join_paths(config.blast_path, [gct.BLASTP_ALIAS])
	makeblastdb_path = fo.join_paths(config.blast_path, [gct.MAKEBLASTDB_ALIAS])
	blastdb_aliastool_path = fo.join_paths(config.blast_path, [gct.BLASTDB_ALIASTOOL_ALIAS])

	# Shorten sequence IDs to avoid issues with long identifiers when creating BLAST DBs
	lcl_distinct_prots = fo.join_paths(clustering_dir, ["distinct_proteins_lcl.fasta"])
	# Return mapping between new short IDs and original IDs
	_ = fao.integer_headers(representative_pfasta, lcl_distinct_prots, start=1, limit=50000, prefix=gct.BLASTDB_LCL_PREFIX, id_map=True)
	# Create BLAST DB
	blast_db_dir = fo.join_paths(clustering_dir, ["BLASTp_db"])
	fo.create_directory(blast_db_dir)
	blast_db = fo.join_paths(blast_db_dir, ["distinct_proteins"])
	db_std = bw.make_blast_db(makeblastdb_path, lcl_distinct_prots, blast_db, "prot")
	# Delete FASTA file with LCL IDs
	fo.remove_files([lcl_distinct_prots])

	# Create a second FASTA file with the sequence ID format returned by BLAST ('SEQ' instead of 'lcl|SEQ')
	seq_distinct_prots = fo.join_paths(clustering_dir, ["distinct_proteins_SEQ.fasta"])
	id_mapping = fao.integer_headers(representative_pfasta, seq_distinct_prots, start=1, limit=50000, prefix=gct.BLASTDB_SEQ_PREFIX, id_map=True)
	# Return inverse mapping to convert original IDs to renamed IDs for BLASTp
	inverse_id_mapping = im.invert_dictionary(id_mapping)

	# All-vs-all BLASTp per cluster
	if len(clusters) > 0:
		logger.info(f"Clusters to BLAST: {len(clusters)}")
		logger.info("Performing all-vs-all BLASTp per cluster...")
		blast_results, blast_results_dir = cf.blast_clusters(clusters,
													seq_distinct_prots,
													inverse_id_mapping,
													clustering_dir, blast_db,
													blastp_path, config.cpu_cores,
													blastdb_aliastool_path)

		blast_files = im.flatten_list(blast_results)

		# Compute and exclude based on BSR
		logger.info("Removing sequences based on high BSR...")
		bsr_excluded = []
		for file in blast_files:
			blast_lines = fo.read_tabular(file)
			# Convert renamed IDs back to original IDs
			blast_lines = [[id_mapping[r[0]]] + r[1:4] + [id_mapping[r[4]]] + r[5:] for r in blast_lines]
			bsr_excluded.append(sm.apply_bsr(blast_lines,
										 dna_index,
										 config.blast_score_ratio))

		# Merge BSR results
		bsr_excluded = set(im.flatten_list(bsr_excluded))
		schema_seqids = list(set(schema_seqids) - bsr_excluded)
		logger.info(f"Removed {len(bsr_excluded)} sequences.")

		# Write list of excluded to file
		blast_excluded_outfile = fo.join_paths(blast_results_dir, ["bsr_excluded.txt"])
		fo.write_lines(bsr_excluded, blast_excluded_outfile)

	# Perform final BLAST to identify similar sequences that do not share many/any kmers
	logger.info(f"{len(schema_seqids)} sequences to compare in final BLASTp.")

	# Sort seqids before final BLASTp to ensure consistent results
	schema_seqids = im.sort_iterable(schema_seqids, sort_key=lambda x: x.lower())

	# Create directory for final BLASTp
	final_blast_dir = fo.join_paths(temp_directory, ["3_final_blast"])
	fo.create_directory(final_blast_dir)

	# Create FASTA file with remaining sequences
	quasi_schema_file = os.path.join(final_blast_dir, "remaining_sequences.fasta")
	fao.get_sequences_by_id(proteins, schema_seqids, quasi_schema_file)

	# Shorten sequence IDs to avoid issues with long identifiers when creating BLAST DBs
	lcl_quasi_schema = fo.join_paths(final_blast_dir, ["remaining_sequences_LCL.fasta"])
	_ = fao.integer_headers(quasi_schema_file, lcl_quasi_schema, start=1, limit=50000, prefix=gct.BLASTDB_LCL_PREFIX, id_map=True)
	# Create BLASTp database
	blast_db = fo.join_paths(final_blast_dir, ["remaining_sequences"])
	db_std = bw.make_blast_db(makeblastdb_path, lcl_quasi_schema, blast_db, "prot")
	# Delete FASTA file with LCL IDs
	fo.remove_files([lcl_quasi_schema])

	# Create a second FASTA file with the sequence ID format returned by BLAST ('SEQ' instead of 'lcl|SEQ')
	seq_quasi_schema = fo.join_paths(final_blast_dir, ["remaining_sequences_SEQ.fasta"])
	id_mapping = fao.integer_headers(quasi_schema_file, seq_quasi_schema, start=1, limit=50000, prefix=gct.BLASTDB_SEQ_PREFIX, id_map=True)
	# Return inverse mapping to convert original IDs to renamed IDs for BLASTp
	inverse_id_mapping = im.invert_dictionary(id_mapping)

	# Divide FASTA file into groups of 100 sequences to reduce
	# execution time for large sequence sets
	split_dir = fo.join_paths(final_blast_dir, ["cds_subsets"])
	fo.create_directory(split_dir)
	splitted_fastas = fao.split_seqcount(seq_quasi_schema, split_dir, 100)

	# Create directory to store results from final BLASTp
	final_blastp_dir = fo.join_paths(final_blast_dir, ["BLAST_results"])
	fo.create_directory(final_blastp_dir)
	blast_outputs = [f"{final_blastp_dir}/{fo.file_basename(i[0], False)}_blast_out.tsv" for i in splitted_fastas]

	# Add common arguments to all sublists
	blast_inputs = [[blastp_path, blast_db, file[0],
					 blast_outputs[i], 1, 1, bw.run_blast]
					for i, file in enumerate(splitted_fastas)]

	logger.info("Performing final BLASTp...")
	blast_results = mo.map_async_parallelizer(blast_inputs,
											  mo.function_helper,
											  config.cpu_cores)

	# Concatenate files with BLASTp results
	blast_output = fo.join_paths(final_blast_dir, ["blast_out_concat.tsv"])
	blast_output = fo.concatenate_files(blast_outputs, blast_output)
	blast_lines = fo.read_tabular(blast_output)
	# Convert renamed IDs back to original IDs
	blast_lines = [[id_mapping[r[0]]] + r[1:4] + [id_mapping[r[4]]] + r[5:] for r in blast_lines]
	final_excluded = sm.apply_bsr(blast_lines,
								  dna_index,
								  config.blast_score_ratio)
	schema_seqids = list(set(schema_seqids) - set(final_excluded))
	logger.info(f"Removed {len(final_excluded)} sequences highly similar to other sequences.")

	# Create file with the schema representative sequences
	loci_representatives = os.path.join(final_blast_dir, "loci_representatives.fasta")
	fao.get_sequences_by_id(dna_index, schema_seqids, loci_representatives)

	# Create schema directory and FASTA files
	logger.info(f"Creating schema seed with {len(schema_seqids)} loci in {config.output_directory}")
	schema_files = create_schema_structure(loci_representatives, config.output_directory, config.schema_name)

	return [schema_files, temp_directory]


def main(config):
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
	clustering_similarity :float
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
	pyrodigal_mode : str
		Pyrodigal running mode ("single" or "meta").
	no_cleanup : bool
		If provided, intermediate files generated during process
		execution are not removed at the end.
	"""
	# Sort paths to FASTA files
	config.input_path = im.sort_iterable(config.input_path, sort_key=lambda x: x.lower())
	results = create_schema_seed(config)

	# Copy Pyrodigal Training File (PTF) to schema directory
	schema_dir = os.path.join(config.output_directory, config.schema_name)
	if config.gene_predictor == "pyrodigal" and config.validated_gene_prediction_options.pyrodigal_training_file is not None:
		fo.copy_file(config.validated_gene_prediction_options.pyrodigal_training_file, schema_dir)
		logger.info(f"Copied Pyrodigal training file to {schema_dir}")

	# Write schema config file
	schema_config = cf.write_schema_config(config, gct.__version__, schema_dir)
	logger.info(f"Wrote schema config values to {schema_config[1]}")

	# Create the file with the list of genes/loci
	cf.write_gene_list(schema_dir)
	logger.info(f"Wrote list of loci to {os.path.join(schema_dir, gct.GENE_LIST_BASENAME)}")

	# Remove temporary files
	if config.no_cleanup is False:
		exists = fo.delete_directory(results[1])

	return len(results[0])
