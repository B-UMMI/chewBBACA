#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This is the main script of the chewBBACA suite.
It parses the options and arguments provided through the command
line and calls the specified module.
"""


import os
import sys
import shutil
import argparse

try:
	from __init__ import __version__
	from PredictGenes import predict_genes
	from CreateSchema import create_schema
	from AlleleCall import allele_call
	from SchemaEvaluator import evaluate_schema
	from AlleleCallEvaluator import evaluate_calls
	from PrepExternalSchema import adapt_schema
	from UniprotFinder import annotate_schema
	from ExtractCgMLST import determine_cgmlst
	from SubsetResults import subset_results
	from HashProfiles import hash_profiles
	from GetAlleles import get_alleles
	from ComputeDistances import compute_distances
	from ComputeMSA import compute_msa
	from MergeResults import merge_results
	from utils import (gene_prediction as gp,
					   # profiles_sqlitedb as ps,
					   process_datetime as pdt,
					   constants as ct,
					   parameters_validation as pv,
					   file_operations as fo)

	from utils.parameters_validation import ModifiedHelpFormatter

	from CHEWBBACA_NS import (download_schema, upload_schema,
							  synchronize_schema, stats_requests)
except ModuleNotFoundError:
	from CHEWBBACA import __version__
	from CHEWBBACA.PredictGenes import predict_genes
	from CHEWBBACA.CreateSchema import create_schema
	from CHEWBBACA.AlleleCall import allele_call
	from CHEWBBACA.SchemaEvaluator import evaluate_schema
	from CHEWBBACA.AlleleCallEvaluator import evaluate_calls
	from CHEWBBACA.PrepExternalSchema import adapt_schema
	from CHEWBBACA.UniprotFinder import annotate_schema
	from CHEWBBACA.ExtractCgMLST import determine_cgmlst
	from CHEWBBACA.SubsetResults import subset_results
	from CHEWBBACA.HashProfiles import hash_profiles
	from CHEWBBACA.GetAlleles import get_alleles
	from CHEWBBACA.ComputeDistances import compute_distances
	from CHEWBBACA.ComputeMSA import compute_msa
	from CHEWBBACA.MergeResults import merge_results
	from CHEWBBACA.utils import (gene_prediction as gp,
								 # profiles_sqlitedb as ps,
								 process_datetime as pdt,
								 constants as ct,
								 parameters_validation as pv,
								 file_operations as fo)

	from CHEWBBACA.utils.parameters_validation import ModifiedHelpFormatter

	from CHEWBBACA.CHEWBBACA_NS import (download_schema, upload_schema,
										synchronize_schema, stats_requests)


@pdt.process_timer
def run_predict_genes():
	"""Run the PredictGenes module to predict genes from a set of input genome assemblies."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py PredictGenes --input-files <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='PredictGenes',
									 description='Predict genes from a set of input genome assemblies.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/PredictGenes.html')

	parser.add_argument('PredictGenes', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-files', type=str,
						required=False, dest='input_files',
						help='Path to the directory that contains the input '
							 'files or to a file with a list of full paths to '
							 'the input files, one per line. Input files can be '
							 'in FASTA or GenBank format.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory where the process will store the files with the predicted CDSs.')

	parser.add_argument('--ptf', '--pyrodigal-training-file', type=str,
						required=False, dest='pyrodigal_training_file',
						help='Path to the Pyrodigal training file used to predict '
							 'CDSs. The translation table used to create this file '
							 'overrides any value passed to `--t`, `--translation-table`.')

	parser.add_argument('--t', '--translation-table', type=pv.translation_table_type,
						required=False, default=ct.GENETIC_CODES_DEFAULT, dest='translation_table',
						help='Genetic code used for gene prediction. This value is ignored if a valid '
							 'training file is passed to `--ptf`, `--training-file`.')

	parser.add_argument('--pm', '--pyrodigal-mode', required=False,
						choices=['single', 'meta'],
						default='single', dest='pyrodigal_mode',
						help='Pyrodigal running mode ("single" for finished genomes, reasonable quality '
							 'draft genomes and big viruses. "meta" for metagenomes, low quality draft '
							 'genomes, small viruses, and small plasmids).')

	parser.add_argument('--tr', '--training-reference', type=str,
						required=False, dest='training_reference',
						help='Path to a reference genome in FASTA format used to create a Pyrodigal training file to predict CDSs.')

	parser.add_argument('--jt', '--just-training', action='store_true',
						required=False, dest='just_training',
						help='Create a training file based on the reference genome and exit.')

	parser.add_argument('--pof', '--pyrodigal-output-formats', nargs='+', type=str,
						required=False, default=['genes'], choices=['genes', 'translations', 'gff', 'genbank', 'scores'],
						dest='pyrodigal_output_formats',
						help='Output file formats created by Pyrodigal. Users can select a single or multiple options from '
							 '`genes` (CDSs in FASTA format), `translations` (translated CDSs in FASTA format), '
							 '`gff` (GFF file format), `genbank` (GenBank file format), and `scores` (TSV file '
							 'with the scores for all predicted CDSs). Default is `genes`.')

	parser.add_argument('--pmc', '--pyrodigal-minimum-confidence', type=float,
						required=False, dest='pyrodigal_minimum_confidence',
						help='Minimum confidence value for CDSs predicted with Pyrodigal. Predicted CDSs with '
							 'a confidence score lower than this value are excluded.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores that will be used to run the process (chewie resets to a lower '
							 'value if it is equal to or exceeds the total number of available CPU cores).')

	args = parser.parse_args()
	del args.PredictGenes

	# Create output directory
	created = fo.create_directory(args.output_directory)
	if created is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)
	print(f'Output directory: {args.output_directory}')

	if args.training_reference:
		print(f'Creating training file based on {args.training_reference}...')
		training_file = gp.create_training_file(args.training_reference, args.output_directory, args.translation_table)
		print(f'Training file saved to {training_file}')
		if args.just_training:
			sys.exit(ct.JUST_TRAINING)
		else:
			args.pyrodigal_training_file = training_file

	# Delete training_reference and just_training arguments
	# These arguments are no passed to the main script of the PredictGenes module
	del args.training_reference
	del args.just_training

	# Check if user passed PTF
	if args.pyrodigal_training_file:
		# Check if PTF exists
		if not os.path.isfile(args.pyrodigal_training_file):
			sys.exit(ct.INVALID_PTF_PATH)
		else:
			# Get translation table used to create training file
			ptf_table = gp.read_training_file(args.pyrodigal_training_file).translation_table
			args.translation_table = ptf_table
	else:
		if not args.translation_table:
			args.translation_table = ct.GENETIC_CODES_DEFAULT
			print(f'Did not provide training file and translation table. Using default translation table ({ct.GENETIC_CODES_DEFAULT})')

	print(f'Pyrodigal training file: {args.pyrodigal_training_file}')
	print(f'Translation table: {args.translation_table}')
	print(f'Pyrodigal mode: {args.pyrodigal_mode}')
	if args.pyrodigal_mode == 'meta' and args.pyrodigal_training_file is not None:
		print('Pyrodigal mode is set to "meta". Will not use provided training file.')
		args.pyrodigal_training_file = None

	print(f'CPU cores: {args.cpu_cores}')

	genome_list = fo.join_paths(args.output_directory, [ct.GENOME_LIST])
	args.input_files, total_inputs = pv.check_input_type(args.input_files, genome_list)

	# Detect if some inputs share the same unique prefix
	repeated_prefixes = pv.check_unique_prefixes(genome_list)
	# Detect if filenames include blank spaces
	blank_spaces = pv.check_blanks(genome_list)

	# Predict CDSs
	predict_genes.main(**vars(args))

	# Delete temporary file with paths to input genomes
	fo.remove_files([genome_list])


@pdt.process_timer
def run_create_schema():
	"""Run the CreateSchema module to create a schema seed."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py CreateSchema --input-files <path> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='CreateSchema',
									 description='Create a schema seed.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='It is strongly advised to provide a training file to '
											'create a schema. Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/CreateSchema.html')

	parser.add_argument('CreateSchema', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-files', type=str,
						required=True, dest='input_files',
						help='Path to the directory that contains the input '
							 'FASTA files or to a file with a list of full '
							 'paths to FASTA files, one per line.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Output directory where the process will store '
							 'intermediate files and create the schema\'s '
							 'directory.')

	parser.add_argument('--n', '--schema-name', type=str,
						required=False, default='schema_seed',
						dest='schema_name',
						help='Name given to the schema folder.')

	parser.add_argument('--ptf', '--training-file', type=str,
						required=False, dest='ptf_path',
						help='Path to the Pyrodigal training file used by Pyrodigal '
							 'to predict genes. The translation table used to create '
							 'this file overrides any value passed to `--t`, '
							 '`--translation-table`. This file is copied '
							 'to the schema folder to be used for allele calling.')

	parser.add_argument('--bsr', '--blast-score-ratio', type=pv.bsr_type,
						required=False, default=ct.DEFAULT_BSR,
						dest='blast_score_ratio',
						help='BLAST Score Ratio (BSR) value. The BSR is computed '
							 'for each BLASTp alignment and aligned sequences with '
							 'a BSR >= than the defined value are considered to be '
							 'alleles of the same gene.')

	parser.add_argument('--l', '--minimum-length', type=pv.minimum_sequence_length_type,
						required=False, default=201, dest='minimum_length',
						help='Minimum sequence length value. Predicted coding '
							 'sequences (CDSs) shorter than this value are excluded.')

	parser.add_argument('--t', '--translation-table', type=pv.translation_table_type,
						required=False, dest='translation_table',
						help='Genetic code used to predict genes and'
							 ' to translate coding DNA sequences (CDSs). '
							 'This value is ignored if a valid training file '
							 'is passed to `--ptf`, `--training-file`.')

	parser.add_argument('--st', '--size-threshold', type=pv.size_threshold_type,
						required=False, default=0.2, dest='size_threshold',
						help='Coding sequence (CDS) size variation threshold. '
							 'Added to the schema\'s config file to identify '
							 'alleles with a size that deviates from the locus '
							 'length mode during the allele calling process.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores that will be '
							 'used to run the process (chewie '
							 'resets to a lower value if it is equal to '
							 'or exceeds the total number of available '
							 'CPU cores).')

	parser.add_argument('--b', '--blast-path', type=pv.check_blast,
						required=False, default='', dest='blast_path',
						help='Path to the directory that contains the '
							 'BLAST executables.')

	parser.add_argument('--pm', '--pyrodigal-mode', required=False,
						choices=['single', 'meta'],
						default='single', dest='pyrodigal_mode',
						help='Pyrodigal running mode ("single" for '
							 'finished genomes, reasonable quality '
							 'draft genomes and big viruses. "meta" '
							 'for metagenomes, low quality draft '
							 'genomes, small viruses, and small '
						 	 'plasmids).')

	parser.add_argument('--pmc', '--pyrodigal-minimum-confidence', type=float,
						required=False, dest='pyrodigal_minimum_confidence',
						help='Minimum confidence value for CDSs predicted with Pyrodigal. Predicted CDSs with '
							 'a confidence score lower than this value are excluded.')

	parser.add_argument('--cds', '--cds-input', required=False,
						action='store_true', dest='cds_input',
						help='If provided, chewBBACA skips the gene '
							 'prediction step and assumes the input FASTA '
							 'files contain coding sequences.')

	parser.add_argument('--no-cds-renaming', required=False, action='store_true',
						dest='no_cds_renaming',
						help='Do not rename the sequence/CDS identifiers when using the `--cds` '
							 'option. Provide this parameter when the input FASTA files containing '
							 'CDSs were generated by the PredictGenes module or if you are sure that '
							 'the CDS identifiers conform to the format used by chewBBACA (the input '
							 'file basename and an integer joined by "_").')

	parser.add_argument('--no-cleanup', required=False, action='store_true',
						dest='no_cleanup',
						help='If provided, intermediate files generated '
							 'during process execution are not deleted at '
							 'the end.')

	args = parser.parse_args()
	del args.CreateSchema

	# Check if user passed PTF
	if args.ptf_path:
		# Check if PTF exists
		if not os.path.isfile(args.ptf_path):
			sys.exit(ct.INVALID_PTF_PATH)
		else:
			# Get translation table used to create training file
			ptf_table = gp.read_training_file(args.ptf_path).translation_table
			args.translation_table = ptf_table
			print('Provided training file. Using translation table used to create training file.')
	else:
		if not args.translation_table:
			args.translation_table = ct.GENETIC_CODES_DEFAULT
			print(f'Did not provide training file and translation table. Using default translation table ({ct.GENETIC_CODES_DEFAULT})')

	print(f'Pyrodigal training file: {args.ptf_path}')
	print(f'Pyrodigal mode: {args.pyrodigal_mode}')
	if args.pyrodigal_mode == 'meta' and args.ptf_path is not None:
		print('Pyrodigal mode is set to "meta". Will add training file to '
			  'the schema, but will not use it for gene prediction during '
			  'schema creation.')
		args.ptf_path = None

	# Check if translation table is supported
	pv.translation_table_type([args.translation_table])
	print(f'Translation table: {args.translation_table}')

	print(f'CPU cores: {args.cpu_cores}')
	print(f'BLAST Score Ratio: {args.blast_score_ratio}')
	print(f'Minimum sequence length: {args.minimum_length}')
	print(f'Size threshold: {args.size_threshold}')

	# Create output directory
	created = fo.create_directory(args.output_directory)
	if created is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)
	print(f'Output directory: {args.output_directory}')

	genome_list = fo.join_paths(args.output_directory, [ct.GENOME_LIST])
	args.input_files, total_inputs = pv.check_input_type(args.input_files, genome_list)
	# Detect if any input files share the same basename
	repeated_prefixes = pv.check_unique_prefixes(genome_list)
	# Detect if filenames include blank spaces
	blank_spaces = pv.check_blanks(genome_list)
	# Since v3.5.1, it is not mandatory for basenames to be shorter than 30 chars
	# This means that the loci IDs defined by the CreateSchema module can be longer
	# than the ones defined by previous versions if users provide input files with
	# long basenames

	print(f'Number of inputs: {total_inputs}')

	# Add clustering parameters
	args.word_size = ct.WORD_SIZE_DEFAULT
	args.window_size = ct.WINDOW_SIZE_DEFAULT
	args.clustering_sim = ct.CLUSTERING_SIMILARITY_DEFAULT
	args.representative_filter = ct.REPRESENTATIVE_FILTER_DEFAULT
	args.intra_filter = ct.INTRA_CLUSTER_DEFAULT

	# Run the CreateSchema process
	nloci = create_schema.main(**vars(args))
	print(f'Created schema seed with {nloci} loci.')

	schema_dir = os.path.join(args.output_directory, args.schema_name)
	# Copy Pyrodigal Training File (PTF) to schema directory
	ptf_hash = None
	if args.ptf_path is not None:
		shutil.copy(args.ptf_path, schema_dir)
		# Determine PTF checksum
		ptf_hash = fo.hash_file(args.ptf_path, 'blake2b')
		print(f'Copied Pyrodigal training file to {schema_dir}')

	# Write schema config file
	args.minimum_length = ct.MSL_MIN
	args.ptf_path = ptf_hash
	schema_config = pv.write_schema_config(vars(args), __version__, schema_dir)
	print(f'Wrote schema config values to {schema_config[1]}')

	# Create the file with the list of genes/loci
	pv.write_gene_list(schema_dir)
	print(f'Wrote list of loci to {os.path.join(schema_dir, ct.GENE_LIST_BASENAME)}')

	# Remove temporary file with paths to input genomes
	fo.remove_files([genome_list])


@pdt.process_timer
def run_allele_call():
	"""Run the AlleleCall module to perform allele calling."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py AlleleCall --input-files <path> --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='AlleleCall',
									 description='Determine the allelic profiles of a set of genomes.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='It is strongly advised to perform allele calling '
											'with the default schema parameters to ensure '
											'more consistent results. Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/AlleleCall.html')

	parser.add_argument('AlleleCall', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-files', nargs='?', type=str,
						required=True, dest='input_files',
						help='Path to the directory that contains the input '
							 'FASTA files or to a file with a list of full '
							 'paths to FASTA files, one per line.')

	parser.add_argument('-g', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the schema directory. The schema '
							 'directory contains the loci FASTA files and '
							 'a folder named "short" that contains the '
							 'FASTA files with the loci representative '
							 'alleles.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Output directory where the process will store '
							 'intermediate files and allele calling results '
							 '(will create a subdirectory named "results_<TIMESTAMP>" '
							 'if the path passed by the user already exists).')

	parser.add_argument('--ptf', '--training-file', type=str,
						required=False, dest='ptf_path',
						help='Path to the Pyrodigal training file used by Pyrodigal '
							 'to predict genes. Default is to use the training file '
							 'included in the schema\'s directory. The translation '
							 'table used to create this file overrides any value '
							 'passed to `--t`, `--translation-table`.')

	parser.add_argument('--gl', '--genes-list', type=str,
						required=False, default=False, dest='genes_list',
						help='Path to a file with the list of genes/loci to '
							 'perform allele calling. The file must include '
							 'the full paths to the loci FASTA files or the loci '
							 'IDs, one per line. The process will perform allele '
							 'calling only for the subset of genes provided in '
							 'the file.')

	parser.add_argument('--bsr', '--blast-score-ratio', type=pv.bsr_type,
						required=False, dest='blast_score_ratio',
						help='BLAST Score Ratio (BSR) value. The BSR is computed '
							 'for each BLASTp alignment and aligned sequences with '
							 'a BSR >= than the defined value are considered to be '
							 'alleles of the same gene.')

	parser.add_argument('--l', '--minimum-length', type=pv.minimum_sequence_length_type,
						required=False, dest='minimum_length',
						help='Minimum sequence length value. Predicted coding '
							 'sequences (CDSs) shorter than this value are excluded.')

	parser.add_argument('--t', '--translation-table', type=pv.translation_table_type,
						required=False, dest='translation_table',
						help='Genetic code used to predict genes and'
							 ' to translate coding DNA sequences (CDSs). '
							 'This value will be ignored if a training file is used.')

	parser.add_argument('--st', '--size-threshold', type=pv.size_threshold_type,
						required=False, dest='size_threshold',
						help='Coding sequence (CDS) size variation threshold. '
							 'At the default value of 0.2, CDSs with a size that '
							 'deviates +-20 percent from the locus length mode '
							 'are classified as ASM/ALM.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores that will be '
							 'used to run the process (chewie '
							 'resets to a lower value if it is equal to '
							 'or exceeds the total number of available '
							 'CPU cores).')

	parser.add_argument('--b', '--blast-path', type=pv.check_blast,
						required=False, default='', dest='blast_path',
						help='Path to the directory that contains the '
							 'BLAST executables.')

	parser.add_argument('--pm', '--pyrodigal-mode', type=str,
						required=False, choices=['single', 'meta'],
						default='single', dest='pyrodigal_mode',
						help='Pyrodigal running mode ("single" for '
							 'finished genomes, reasonable quality '
							 'draft genomes and big viruses. "meta" '
							 'for metagenomes, low quality draft '
							 'genomes, small viruses, and small '
							 'plasmids).')

	parser.add_argument('--pmc', '--pyrodigal-minimum-confidence', type=float,
						required=False, dest='pyrodigal_minimum_confidence',
						help='Minimum confidence value for CDSs predicted with Pyrodigal. Predicted CDSs with '
							 'a confidence score lower than this value are excluded.')

	parser.add_argument('--cds', '--cds-input', action='store_true',
						required=False, dest='cds_input',
						help='If provided, chewBBACA skips the gene '
							 'prediction step and assumes the input FASTA '
							 'files contain coding sequences (one FASTA '
							 'file per strain).')

	parser.add_argument('--no-inferred', required=False,
						action='store_true', dest='no_inferred',
						help='If provided, the process will not add '
							 'the sequences of inferred alleles (INF) to the '
							 'schema. Allelic profiles will still include '
							 'the allele identifiers attributed to the '
							 'inferred alleles. Use this parameter if the '
							 'schema is being accessed by multiple '
							 'processes/users simultaneously.')

	parser.add_argument('--output-unclassified', required=False,
						action='store_true', dest='output_unclassified',
						help='Create a Fasta file with the coding sequences '
							 '(CDSs) that were not classified.')

	parser.add_argument('--output-missing', required=False,
						action='store_true', dest='output_missing',
						help='Create a Fasta file with coding sequences (CDSs) '
							 'classified as NIPH, NIPHEM, ASM, ALM, PLOT3, '
							 'PLOT5 and LOTSC.')

	parser.add_argument('--output-novel', required=False,
						action='store_true', dest='output_novel',
						help='Create a Fasta file with the novel alleles '
							 'inferred during allele calling. The '
							 'sequence headers include the locus and allele '
							 'identifiers attributed by chewBBACA based '
							 'on the allele calling results.')

	parser.add_argument('--output-masked', required=False,
						action='store_true', dest='output_masked',
						help='Create a TSV file with the masked allelic profiles. '
							 'The masking process removes the `INF-` prefix from '
							 'inferred alleles and substitutes all special classes '
							 '(NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5, LOTSC, PAMA) '
							 'with `0`.')

	parser.add_argument('--no-cleanup', required=False,
						action='store_true', dest='no_cleanup',
						help='If provided, intermediate files generated '
							 'during process execution are not removed at '
							 'the end.')

	# parser.add_argument('--db', '--store-profiles', required=False,
	#                     action='store_true', dest='store_profiles',
	#                     help='If the profiles in the output matrix '
	#                          'should be stored in the local SQLite '
	#                          'database.')

	parser.add_argument('--no-cds-renaming', required=False, action='store_true',
						dest='no_cds_renaming',
						help='Do not rename the sequence/CDS identifiers when using the `--cds` '
							 'option. Provide this parameter when the input FASTA files containing '
							 'CDSs were generated by the PredictGenes module or if you are sure that '
							 'the CDS identifiers conform to the format used by chewBBACA (the input '
							 'file basename and an integer joined by "_").')

	parser.add_argument('--force-continue', required=False,
						action='store_true', dest='force_continue',
						help='If provided, chewie will not warn users and ask '
							 'for permission to continue if any of the provided '
							 'argument values does not match the values in the '
							 'config file.')

	parser.add_argument('--mode', type=int, required=False,
						choices=[1, 2, 3, 4], default=4,
						help='Execution mode (1: only exact matches at DNA '
							 'level; 2: exact matches at DNA and Protein '
							 'level; 3: exact matches and minimizer-based '
							 'clustering to find similar alleles based on '
							 'BSR+0.1; 4: run the full process to find '
							 'exact matches and similar matches based on '
							 'BSR value, including the determination of new '
							 'representative alleles to add to the schema).')

	args = parser.parse_args()

	# Check if input schema path exists
	if not os.path.exists(args.schema_directory):
		sys.exit(ct.SCHEMA_PATH_MISSING)

	# Verify schema files
	schema_files = os.listdir(args.schema_directory)
	# Exit if there is no 'short' directory or if there are no FASTA files
	if 'short' not in schema_files or len(fo.filter_by_extension(schema_files, ['.fasta'])[0]) == 0:
		sys.exit(ct.SCHEMA_INVALID_PATH)
	# Check if 'short' directory includes files terminating in 'bsr.txt'
	schema_short_path = fo.join_paths(args.schema_directory, ['short'])
	schema_short_files = os.listdir(schema_short_path)
	if any([file.endswith('bsr.txt') for file in schema_short_files]):
		sys.exit(ct.ADAPT_LEGACY_SCHEMA)

	config_file = os.path.join(args.schema_directory, ct.SCHEMA_CONFIG_BASENAME)
	# Legacy schemas do not have config file
	# Tell users to adapt with PrepExternalSchema module
	if os.path.isfile(config_file) is False:
		sys.exit(ct.ADAPT_LEGACY_SCHEMA)
	else:
		schema_params = fo.pickle_loader(config_file)
		# Chek if user provided different values
		run_params = pv.solve_conflicting_arguments(schema_params, args.ptf_path,
													args.blast_score_ratio, args.translation_table,
													args.minimum_length, args.size_threshold,
													args.force_continue, config_file, args.schema_directory)
		args.ptf_path = run_params['ptf_path']
		args.blast_score_ratio = run_params['bsr']
		args.translation_table = run_params['translation_table']
		args.minimum_length = run_params['minimum_locus_length']
		args.size_threshold = run_params['size_threshold']

	# Create output directory
	created = fo.create_directory(args.output_directory)
	# Output directory exists
	# Create a subdirectory to store intermediate files and results
	if created is False:
		current_time = pdt.get_datetime()
		current_time_str = pdt.datetime_str(current_time,
											date_format='%Y%m%dT%H%M%S')
		results_dir = fo.join_paths(args.output_directory,
									['results_{0}'.format(current_time_str)])
		created = fo.create_directory(results_dir)
		args.output_directory = results_dir
		print(f'Output directory exists. Will store results in {results_dir}\n')

	loci_list = fo.join_paths(args.output_directory, [ct.LOCI_LIST])
	# User provided a list of genes to call
	if args.genes_list is not False:
		loci_list = pv.validate_loci_list(args.genes_list, loci_list,
											args.schema_directory)
	# Working with the whole schema
	else:
		loci_list, total_loci = pv.check_input_type(args.schema_directory, loci_list)

	genome_list = fo.join_paths(args.output_directory, [ct.GENOME_LIST])
	genome_list, total_inputs = pv.check_input_type(args.input_files, genome_list)
	# Detect if any input files share the same basename
	repeated_prefixes = pv.check_unique_prefixes(genome_list)
	# Detect if filenames include blank spaces
	blank_spaces = pv.check_blanks(genome_list)
	# Since v3.5.1, it is not mandatory for basenames to be shorter than 30 chars
	# This means that the CDS IDs defined by the AlleleCall module can be longer
	# than the ones defined by previous versions if users provide input files with
	# long basenames. This also means that the output files will include the full
	# basename of the input files

	# Determine if schema was downloaded from Chewie-NS
	ns_config = fo.join_paths(args.schema_directory, ['.ns_config'])
	args.ns = os.path.isfile(ns_config)

	# Add clustering arguments
	args.word_size = ct.WORD_SIZE_DEFAULT
	args.window_size = ct.WINDOW_SIZE_DEFAULT
	args.clustering_sim = ct.CLUSTERING_SIMILARITY_DEFAULT

	# Single dictionary with most arguments
	config = {'Minimum sequence length': args.minimum_length,
			  'Size threshold': args.size_threshold,
			  'Translation table': args.translation_table,
			  'BLAST Score Ratio': args.blast_score_ratio,
			  'Word size': args.word_size,
			  'Window size': args.window_size,
			  'Clustering similarity': args.clustering_sim,
			  'Pyrodigal training file': args.ptf_path,
			  'CPU cores': args.cpu_cores,
			  'BLAST path': args.blast_path,
			  'CDS input': args.cds_input,
			  'Pyrodigal mode': args.pyrodigal_mode,
			  'Pyrodigal minimum confidence': args.pyrodigal_minimum_confidence,
			  'Mode': args.mode}

	allele_call.main(genome_list, loci_list, args.schema_directory,
					 args.output_directory, args.no_inferred,
					 args.output_unclassified, args.output_missing,
					 args.output_novel, args.output_masked,
					 args.no_cleanup, args.no_cds_renaming, args.ns, config)

	# if args.store_profiles is True:
	#     updated = ps.store_allelecall_results(args.output_directory, args.schema_directory)

	# Remove temporary files with paths to input genomes and schema files
	fo.remove_files([loci_list])
	fo.remove_files([genome_list])


@pdt.process_timer
def run_evaluate_schema():
	"""Run the SchemaEvaluator module to evaluate a typing schema."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py SchemaEvaluator --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='SchemaEvaluator',
									 description='Build an interactive report for schema evaluation.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='The module can evaluate schemas created with chewBBACA or other external MLST platforms. Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/SchemaEvaluator.html')

	parser.add_argument('SchemaEvaluator', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-g', '--schema-directory', type=str, required=True,
						dest='schema_directory',
						help='Path to the schema\'s directory.')

	parser.add_argument('-o', '--output-directory', type=str, required=True,
						dest='output_directory',
						help='Path to the output directory where the report '
							 'HTML files will be created.')

	parser.add_argument('--gl', '--genes-list', type=str,
						required=False, default=False, dest='genes_list',
						help='Path to a file with the list of loci '
							 'in the schema that the process should '
							 'analyse (one per line, full paths or loci IDs).')

	parser.add_argument('-a', '--annotations', type=str, required=False,
						dest='annotations',
						help='Path to the TSV file created by the '
							 'UniprotFinder module. The annotation data '
							 'is included in a table component.')

	parser.add_argument('--ta', '--translation-table',
						type=pv.translation_table_type, required=False,
						dest='translation_table',
						help='Genetic code used to translate coding '
							 'sequences (CDSs).')

	parser.add_argument('--st', '--size-threshold',
						type=pv.size_threshold_type,
						required=False, dest='size_threshold',
						help='Coding sequence (CDS) size variation threshold. '
							 'The module identifies the alleles with size '
							 'that deviates from the locus length mode +- '
							 'the size threshold.')

	parser.add_argument('--ml', '--minimum-length',
						type=pv.minimum_sequence_length_type,
						required=False, dest='minimum_length',
						help='Minimum sequence length value. The module '
							 'identifies alleles shorter than this value.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	parser.add_argument('--loci-reports', required=False,
						action='store_true', dest='loci_reports',
						help='Create a detailed report page for each locus. '
							 'The locus report includes components with '
							 'relevant data and analysis results, such as '
							 'allele diversity charts, a MSA for the alignment '
							 'of the distinct translated alleles and a tree '
							 'drawn with Phylocanvas based on the MAFFT guide tree.')

	parser.add_argument('--light', action='store_true', required=False,
						dest='light',
						help='Skips MSA computation with MAFFT and does not add '
							 'the Phylogenetic Tree and MSA components to the loci reports.')

	parser.add_argument('--add-sequences', required=False,
						action='store_true', dest='add_sequences',
						help='Adds Code Editor components with the DNA and Protein '
							 'sequences to the loci reports. The Code Editor '
							 'is in readonly mode (allows to search for '
							 'and copy text).')

	args = parser.parse_args()
	del args.SchemaEvaluator

	# Check if input schema path exists
	if not os.path.exists(args.schema_directory):
		sys.exit(ct.SCHEMA_PATH_MISSING)

	# Create output directory
	created = fo.create_directory(args.output_directory)
	if created is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)

	loci_list = fo.join_paths(args.output_directory, [ct.LOCI_LIST])
	# User provided a loci subset to analyse
	if args.genes_list is not False:
		loci_list = pv.validate_loci_list(args.genes_list, loci_list,
										  args.schema_directory)
	# Working with the whole schema
	else:
		loci_list, total_loci = pv.check_input_type(args.schema_directory, loci_list)

	args.genes_list = loci_list

	evaluate_schema.main(**vars(args))

	# Delete file with list of loci that were evaluated
	fo.remove_files([loci_list])


@pdt.process_timer
def run_evaluate_calls():
	"""Run the AlleleCallEvaluator module to evaluate allele calling results."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py AlleleCallEvaluator --input-files <dir> --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='AlleleCallEvaluator',
									 description='Build an interactive report for allele calling results evaluation.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/AlleleCallEvaluator.html')

	parser.add_argument('AlleleCallEvaluator', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-files', type=str, required=True,
						dest='input_files',
						help='Path to the directory that contains the allele '
							 'calling results generated by the AlleleCall module.')

	parser.add_argument('-g', '--schema-directory', type=str, required=True,
						dest='schema_directory',
						help='Path to the schema\'s directory.')

	parser.add_argument('-o', '--output-directory', type=str, required=True,
						dest='output_directory',
						help='Path to the output directory where the module will '
							 'store intermediate files and create the report HTML files.')

	parser.add_argument('-a', '--annotations', type=str, required=False,
						dest='annotations',
						help='Path to the TSV file created by the '
							 'UniprotFinder module.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	parser.add_argument('--light', action='store_true', required=False,
						dest='light',
						help='Do not compute the presence-absence matrix, '
							 'the distance matrix and the Neighbor-Joining tree.')

	parser.add_argument('--no-pa', action='store_true', required=False,
						dest='no_pa',
						help='Do not compute the presence-absence matrix.')

	parser.add_argument('--no-dm', action='store_true', required=False,
						dest='no_dm',
						help='Do not compute the distance matrix.')

	parser.add_argument('--no-tree', action='store_true', required=False,
						dest='no_tree',
						help='Do not compute the Neighbor-Joining tree.')

	parser.add_argument('--cg-alignment', action='store_true', required=False,
						dest='cg_alignment',
						help='Compute the MSA of the core genome loci, even '
							 'if `--no-tree` is provided.')

	args = parser.parse_args()
	del args.AlleleCallEvaluator

	# Check if path to input files exists
	if not os.path.exists(args.input_files):
		sys.exit(ct.MISSING_INPUT_ARG)

	# Create output directory
	created = fo.create_directory(args.output_directory)
	if created is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)

	evaluate_calls.main(**vars(args))


@pdt.process_timer
def run_determine_cgmlst():
	"""Run the ExtractCgMLST module to determine the set of core loci based on allele calling results."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py ExtractCgMLST --input-file <file> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='ExtractCgMLST',
									 description='Determine the set of core loci based on allele calling results.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/ExtractCgMLST.html')

	parser.add_argument('ExtractCgMLST', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-file', type=str,
						required=True, dest='input_file',
						help='Path to the TSV file that contains the allelic '
							 'profiles determined by the AlleleCall module.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the directory where the process '
							 'will store the output files.')

	parser.add_argument('--t', '--threshold', nargs='+', type=float,
						required=False, default=ct.CGMLST_THRESHOLDS,
						dest='threshold',
						help='Loci/genes that constitute the core genome '
							 'must be in a proportion of genomes that is '
							 'at least equal to this value. Provide multiple '
							 'values to compute the core genome for multiple '
							 'threshold values.')

	parser.add_argument('--s', '--step', type=int, required=False,
						default=1, dest='step',
						help='The allele calling results are processed '
							 'iteratively to evaluate the impact of '
							 'adding subsets of the results in computing '
							 'the core genome. The step value controls '
							 'the number of allelic profiles added in each '
							 'iteration until all profiles are included.')

	parser.add_argument('--ca', '--compute-accessory', action='store_true', required=False,
					    dest='compute_accessory',
						help='Determine the set of accessory loci.'
							 'The accessory genome corresponds to all the '
							 'loci not included in the core genome. The '
							 'accessory genome is determined for each core '
							 'genome threshold.')

	parser.add_argument('--ra', '--rarefaction-analysis', action='store_true', required=False,
					    dest='rarefaction_analysis',
						help='Perform rarefaction analysis to evaluate the stability of the core genome '
							 'and the openness of the pangenome. For the rarefaction analysis, samples '
							 'are randomly selected and the number of core loci (for the core genome) or '
							 'total loci (for the pangenome) is computed for each subset from 1 to the '
							 'total number of samples. The average number of core loci (for the core genome) '
							 'or total loci (for the pangenome) is computed for each sample size subset '
							 'using the values of all permutations to plot the rarefaction curve. The '
							 'rarefaction analysis is performed for each threshold. A power law model '
							 'is fitted to the rarefaction curve to estimate the stability of the core genome '
							 'and the openness of the pangenome.')

	parser.add_argument('--pn', '--permutation-number', type=int, required=False, default=100,
					    dest='permutation_number',
						help='Number of permutations for the rarefaction analysis. The rarefaction '
							 'analysis is repeated a number of times equal to the value provided to this '
							 'parameter.')

	parser.add_argument('--ps', '--permutation-samples', type=int, required=False,
					    dest='permutation_samples',
						help='Number of samples randomly selected for each permutation. All '
							 'samples will be used if this value is not provided.')

	parser.add_argument('--el', '--exclude-loci', type=str,
						required=False, default=False, dest='exclude_loci',
						help='Path to a file with a list of loci identifiers to '
							 'exclude from the analysis (one locus identifier '
							 'per line).')

	parser.add_argument('--eg', '--exclude-genomes', type=str,
						required=False, default=False, dest='exclude_genomes',
						help='Path to a file with a list of genome identifiers to '
							 'exclude from the analysis (one genome identifier '
							 'per line).')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Maximum number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a lower '
							 'value if it is equal to or exceeds the total number '
							 'of available CPU cores/threads).')

	args = parser.parse_args()
	del args.ExtractCgMLST

	determine_cgmlst.main(**vars(args))


@pdt.process_timer
def run_subset_results():
	"""Run the SubsetResults module to subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py SubsetResults --input-file <file> --loci-list <file> --samples-list <file> --output-file <file> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='SubsetResults',
									 description='Subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/SubsetResults.html')

	parser.add_argument('SubsetResults', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-directory', type=str,
						required=True, dest='input_directory',
						help='Path to the directory containing the files to be subsetted.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory.')

	parser.add_argument('-l', '--loci-list', type=str,
						required=False, default=False, dest='loci_list',
						help='Path to a TXT/TSV file containing a list of loci to select, one locus identifier per '
							 'line. If the file contains multiple columns, the loci identifiers must be in the first column.')

	parser.add_argument('-s', '--sample-list', type=str,
						required=False, default=False, dest='sample_list',
						help='Path to a TXT/TSV file containing a list of samples to select, one sample identifier '
							 'per line. If the file contains multiple columns, the sample identifiers must be in the first column.')

	parser.add_argument('--inverse-loci', action='store_true',
						required=False, dest='inverse_loci',
						help='If provided, the process will select the loci that are not in the input loci list.')

	parser.add_argument('--inverse-samples', action='store_true',
						required=False, dest='inverse_samples',
						help='If provided, the process will select the samples that are not in the input samples list.')

	args = parser.parse_args()
	del args.SubsetResults

	subset_results.main(**vars(args))


@pdt.process_timer
def run_merge_results():
	"""Run the MergeResults module to merge results files created by chewBBACA."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py MergeResults --input-directories <dir> <dir> ... --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='MergeResults',
									 description='Merge results files created by chewBBACA.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/MergeResults.html')

	parser.add_argument('MergeResults', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-files', nargs='+', type=str,
						required=True, dest='input_files',
						help='Paths to the directories containing the results files '
							 'created by chewBBACA. The results must have been '
							 'determined with the same schema and share all the loci '
							 'or a subset of the loci if using the --common parameter.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory.')

	parser.add_argument('--common', action='store_true',
						required=False, dest='common',
						help='Merge the results based on the subset of '
							 'loci shared between all inputs.')

	args = parser.parse_args()
	del args.MergeResults

	merge_results.main(**vars(args))


@pdt.process_timer
def run_hash_profiles():
	"""Run the HashProfiles module to hash allelic profiles."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py HashProfiles --input-file <file> --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='HashProfiles',
									 description='Hash allelic profiles.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/HashProfiles.html')

	parser.add_argument('HashProfiles', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-file', type=str,
						required=True, dest='input_file',
						help='Path to the TSV file that contains the allelic profiles determined by the AlleleCall module.')

	parser.add_argument('-g', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the schema\'s directory to get the allele sequences and compute the hashes.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory.')

	parser.add_argument('--hash-type', type=str, required=False,
						default='crc32', dest='hash_type',
						help='Hashing algorithm used to hash the profiles. The hashing '
							 'algorithms implemented in the hashlib and zlib Python '
							 'libraries are supported.')

	parser.add_argument('--nrows', type=int, required=False,
						default=100, dest='nrows',
						help='Divide the input file into chunks of this many rows '
							 'to process larger files more efficiently.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a '
							 'lower value if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	args = parser.parse_args()
	del args.HashProfiles

	hash_profiles.main(**vars(args))


@pdt.process_timer
def run_get_alleles():
	"""Run the GetAlleles module to create FASTA files containing the alleles identified by the AlleleCall module."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py GetAlleles --input-file <file> --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='GetAlleles',
									 description='Create FASTA files containing the alleles identified by the AlleleCall module.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/GetAlleles.html')

	parser.add_argument('GetAlleles', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-file', type=str,
						required=True, dest='input_file',
						help='Path to the TSV file containing the allelic profiles.')

	parser.add_argument('-g', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the schema directory.')

	parser.add_argument('--gl', '--genes-list', type=str,
						required=False, dest='genes_list',
						help='Path to a file with the list of genes/loci to '
							 'create FASTA files for. The file must include '
							 'the identifiers of the loci, one per line, '
							 'without the .fasta extension.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a '
							 'lower value if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	parser.add_argument('--distinct', action='store_true',
						required=False, dest='distinct',
						help='Only get distinct alleles.')

	parser.add_argument('--translate', action='store_true',
						required=False, dest='translate',
						help='Create FASTA files with the translated alleles.')

	parser.add_argument('--ta', '--translation-table',
						type=pv.translation_table_type, required=False,
						dest='translation_table',
						help='Genetic code used to translate coding DNA'
							 'sequences (CDSs). If no value is specified, the '
							 'process tries to get the value stored in the '
							 'schema config file. If the schema does not '
							 'include a config file, the process uses the '
							 'default translation table (11).')

	args = parser.parse_args()
	del args.GetAlleles

	if args.translate and not args.translation_table:
		# Check if schema includes config file
		config_file = os.path.join(args.schema_directory, ct.SCHEMA_CONFIG_BASENAME)
		if os.path.isfile(config_file):
			schema_params = fo.pickle_loader(config_file)
			args.translation_table = schema_params['translation_table'][0]
		else:
			# Use default translation table
			args.translation_table = ct.GENETIC_CODES_DEFAULT

	# Create output directory
	created = fo.create_directory(args.output_directory)
	if created is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)

	get_alleles.main(**vars(args))


@pdt.process_timer
def run_adapt_schema():
	"""Run the PrepExternalSchema module to adapt a typing schema."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py PrepExternalSchema --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='PrepExternalSchema',
									 description='Adapt an external schema to be used with chewBBACA.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/PrepExternalSchema.html')

	parser.add_argument('PrepExternalSchema', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-g', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the directory of the schema to adapt. '
							 'The schema must contain one FASTA file per gene/locus.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory where the adapted schema will '
							 'be created.')

	parser.add_argument('--gl', '--genes-list', type=str,
						required=False, default=False, dest='genes_list',
						help='Path to a file with the list of loci '
							 'in the schema that the process should '
							 'adapt (one per line, full paths or loci IDs).')

	parser.add_argument('--ptf', '--training-file', type=str,
						required=False, dest='ptf_path',
						help='Path to the Pyrodigal training file that '
							 'will be included in the directory of the '
							 'adapted schema. The translation table used to create '
							 'this file overrides any value passed to `--t`, '
							 '`--translation-table`.')

	parser.add_argument('--bsr', '--blast-score-ratio', type=pv.bsr_type,
						required=False, default=ct.DEFAULT_BSR,
						dest='blast_score_ratio',
						help='BLAST Score Ratio (BSR) value. The process '
							 'selects representative alleles for each locus '
							 'based on this value. Representative alleles '
							 'are selected until all alleles in a locus align '
							 'against one of the representatives with a BSR '
							 '>= than the specified value.')

	parser.add_argument('--l', '--minimum-length',
						type=pv.minimum_sequence_length_type, required=False,
						default=ct.MSL_MIN, dest='minimum_length',
						help='Minimum sequence length value stored in the '
							 'schema config file. The schema adaptation '
							 'process will only discard sequences smaller '
							 'than this value if the --size-filter parameter '
							 'is provided.')

	parser.add_argument('--t', '--translation-table',
						type=pv.translation_table_type, required=False,
						dest='translation_table',
						help='Genetic code used for allele translation. This '
							 'value is ignored if a valid training file '
							 'is passed to `--ptf`, `--training-file`.')

	parser.add_argument('--st', '--size-threshold', type=pv.size_threshold_type,
						required=False, default=ct.SIZE_THRESHOLD_DEFAULT,
						dest='size_threshold',
						help='Allele size variation threshold value stored in '
							 'the schema config file. The schema adaptation '
							 'process will only discard alleles with a size '
							 'that deviates from the locus length mode +- the '
							 'size theshold value if the --size-filter parameter '
							 'is provided.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	parser.add_argument('--b', '--blast-path', type=pv.check_blast,
						required=False, default='', dest='blast_path',
						help='Path to the directory that contains the '
							 'BLAST executables.')

	parser.add_argument('--size-filter', action='store_true',
						required=False, dest='size_filter',
						help='Apply the minimum length and size threshold'
							 ' values to filter out alleles during schema '
							 'adaptation.')

	args = parser.parse_args()
	del args.PrepExternalSchema

	# Check if user passed PTF
	if args.ptf_path:
		# Check if PTF exists
		if not os.path.isfile(args.ptf_path):
			sys.exit(ct.INVALID_PTF_PATH)
		else:
			# Get translation table used to create training file
			ptf_table = gp.read_training_file(args.ptf_path).translation_table
			args.translation_table = ptf_table
	else:
		if not args.translation_table:
			args.translation_table = ct.GENETIC_CODES_DEFAULT

	# Define output paths
	schema_path = os.path.abspath(args.output_directory)
	schema_short_path = fo.join_paths(schema_path, ['short'])
	output_dirs = [schema_path, schema_short_path]

	# Create output directories
	schema_path_exists = fo.create_directory(schema_path)
	if schema_path_exists is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)
	fo.create_directory(schema_short_path)

	# User provided a list of loci to adapt
	loci_list = fo.join_paths(schema_path, [ct.LOCI_LIST])
	if args.genes_list is not False:
		loci_list = pv.validate_loci_list(args.genes_list, loci_list,
										  args.schema_directory)
	# Working with the whole schema
	else:
		loci_list, total_loci = pv.check_input_type(args.schema_directory, loci_list)

	# Detect if any input files share the same basename
	repeated_prefixes = pv.check_unique_prefixes(loci_list)
	# Detect if filenames include blank spaces
	blank_spaces = pv.check_blanks(loci_list)
	# Since v3.5.1, there is no limitation for the length of the loci basenames
	# This means that chewBBACA can adapt schemas with longer loci IDs

	print(f'Number of cores: {args.cpu_cores}')
	print(f'BLAST Score Ratio: {args.blast_score_ratio}')
	print(f'Translation table: {args.translation_table}')

	# Only apply minimum length and size threshold during schema
	# adaptation if --size-filter parameter is True
	if args.size_filter:
		adaptation_st = args.size_threshold
		adaptation_ml = args.minimum_length
	else:
		adaptation_st = None
		adaptation_ml = 0

	print(f'Using a minimum length value of {adaptation_ml} for schema '
		  f'adaptation and {args.minimum_length} to store in the schema '
		  'config file.')
	print(f'Using a size threshold value of {adaptation_st} for schema '
		  f'adaptation and {args.size_threshold} to store in the schema '
		  'config file.')

	adapt_schema.main(loci_list, output_dirs,
					  args.cpu_cores, args.blast_score_ratio,
					  adaptation_ml, args.translation_table,
					  adaptation_st, args.blast_path)

	# Copy training file to schema directory
	ptf_hash = None
	if args.ptf_path is not None:
		shutil.copy(args.ptf_path, schema_path)
		# Determine PTF checksum
		ptf_hash = fo.hash_file(args.ptf_path, 'blake2b')
		print('Copied Pyrodigal training file to schema directory.')

	# Write schema config file
	args.ptf_path = ptf_hash
	args.word_size = ct.WORD_SIZE_DEFAULT
	args.window_size = ct.WINDOW_SIZE_DEFAULT
	args.clustering_sim = ct.CLUSTERING_SIMILARITY_DEFAULT
	args.representative_filter = ct.REPRESENTATIVE_FILTER_DEFAULT
	args.intra_filter = ct.INTRA_CLUSTER_DEFAULT
	schema_config = pv.write_schema_config(vars(args), __version__, schema_path)

	# Create hidden file with list of loci
	genes_list_file = pv.write_gene_list(schema_path)

	# Delete file with list of loci to adapt
	os.remove(loci_list)


@pdt.process_timer
def run_annotate_schema():
	"""Run the UniprotFinder module to annotate loci in a schema."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py UniprotFinder --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='UniprotFinder',
									 description='Retrieve annotations for loci in a schema.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/UniprotFinder.html')

	parser.add_argument('UniprotFinder', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-g', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the schema\'s directory.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory where the process will '
							 'store intermediate files and save the final '
							 'TSV file with the loci annotations.')

	parser.add_argument('--gl', '--genes-list', type=str,
						required=False, default=False, dest='genes_list',
						help='Path to a file with the list of loci '
							 'in the schema that the process should '
							 'find annotations for (one per line, full '
							 'paths or loci IDs).')

	parser.add_argument('--t', '--protein-table', type=str,
						required=False, dest='protein_table',
						help='Path to the TSV file with coding sequence (CDS) '
							 'coordinate data, "cds_coordinates.tsv", created by '
							 'the CreateSchema process.')

	parser.add_argument('--bsr', type=float, required=False,
						dest='blast_score_ratio',
						default=0.6,
						help='BLAST Score Ratio value. The BSR is only '
							 'used when taxa names are provided to the --taxa '
							 'parameter and local sequences are aligned against '
							 'reference proteomes downloaded from UniProt. Annotations '
							 'are selected based on a BSR >= than the specified value.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	parser.add_argument('--taxa', nargs='+', type=str,
						required=False, dest='taxa',
						help='List of scientific names for a set of taxa. The '
							 'process will download reference proteomes from UniProt '
							 'associated to taxa names that contain any of the '
							 'provided terms. The schema representative alleles are '
							 'aligned against the reference proteomes to assign '
							 'annotations based on high-BSR matches.')

	parser.add_argument('--pm', type=int, required=False,
						default=1, dest='proteome_matches',
						help='Maximum number of proteome matches to report.')

	parser.add_argument('--no-sparql', action='store_true',
						required=False, dest='no_sparql',
						help='Do not search for annotations through '
							 'the UniProt SPARQL endpoint.')

	parser.add_argument('--no-cleanup', action='store_true',
						required=False, dest='no_cleanup',
						help='If provided, intermediate files generated '
							 'during process execution are not removed '
							 'at the end.')

	parser.add_argument('--b', '--blast-path', type=pv.check_blast,
						required=False, default='', dest='blast_path',
						help='Path to the directory that contains the '
							 'BLAST executables.')

	args = parser.parse_args()
	del args.UniprotFinder

	annotate_schema.main(**vars(args))


@pdt.process_timer
def run_compute_distances():
	"""Run the ComputeDistances module to compute pairwise distances based on allele calling results."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py ComputeDistances --input-file <file> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='ComputeDistances',
									 description='Compute pairwise distances based on allele calling results.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/ComputeDistances.html')

	parser.add_argument('ComputeDistances', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-file', type=str,
						required=True, dest='input_file',
						help='Path to a TSV file containing allelic profiles determined by the AlleleCall module.')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory where the process will '
							 'store intermediate and final results.')

	parser.add_argument('--m', '--method', type=str, required=False,
						default='hamming', choices=['hamming', 'jaccard', 'loci', 'core'],
						dest='method',
						help='Distance method used to compute the distance matrix. The module supports '
							 'the hamming, jaccard, loci (number of loci not shared), and core (number of '
							 'different alleles for core loci) methods.')

	parser.add_argument('--outfmt', '--output-format', type=str, required=False,
						choices=['upper_triangular', 'lower_triangular', 'symmetric', 'table'],
						default='upper_triangular',
						dest='output_format',
						help='Output format for the distance matrix (upper_triangular, lower_triangular, symmetric, table).')

	parser.add_argument('--no-mask', action='store_true', required=False,
						dest='no_mask',
						help='Do not mask missing data when computing the distance matrix. '
							 'This option is useful when the input profiles are already masked.')

	parser.add_argument('--similarity', action='store_true', required=False,
						dest='similarity',
						help='Compute similarity values instead of distance values.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a '
							 'lower value if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	args = parser.parse_args()
	del args.ComputeDistances

	compute_distances.main(**vars(args))


@pdt.process_timer
def run_compute_msa():
	"""Run the ComputeMSA module to compute a Multiple Sequence Alignment based on allele calling results."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py ComputeMSA --input-file <file> --schema-directory <dir> --output-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='ComputeMSA',
									 description='Compute a Multiple Sequence Alignment based on allele calling results.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/ComputeMSA.html')

	parser.add_argument('ComputeMSA', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--input-path', type=str,
						required=True, dest='input_path',
						help='Path to a TSV file containing allelic profiles or '
							 'to a folder containing FASTA files. If a TSV file '
							 'containing allelic profiles is provided, it is necessary '
							 'to provide the path to the schema to the `--schema-directory` '
							 'parameter. The module will create a FASTA file with the alleles '
							 'identified in the samples for each schema locus and compute a MSA. '
							 'The loci MSAs are joined to create the complete MSA based on the '
							 'allele calling results. If a path to a folder is provided, the '
							 'module computes a MSA for each FASTA file in the folder, but will '
							 'not attempt to join the MSAs as it does not have the sample '
							 'information (in this case, it is not necessary to pass the schema '
							 'path).')

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest='output_directory',
						help='Path to the output directory where the process will '
							 'store intermediate and final results.')

	parser.add_argument('-g', '--schema-directory', type=str,
						required=False, dest='schema_directory',
						help='Path to the schema\'s directory. This parameter is '
							 'only required if the input is a TSV file with allelic '
							 'profiles.')

	parser.add_argument('--dna-msa', action='store_true',
						required=False, dest='dna_msa',
						help='Converts the protein MSA back to DNA to create an additional '
							 'output file with the DNA MSA.')

	parser.add_argument('--output-variable', action='store_true',
					 	required=False, dest='output_variable',
						help='Output a reduced MSA including only the variable '
							 'positions. If the `--dna-msa` parameter is provided, the '
							 'process will output a reduced MSA for both the protein '
							 'and DNA MSAs.')

	parser.add_argument('--t', '--translation-table', type=int,
						required=False, default=11, dest='translation_table',
						help='Genetic code used for sequence translation.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a '
							 'lower value if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	parser.add_argument('--only-loci-msas', action='store_true',
						required=False, dest='only_loci_msas',
						help='Do not compute the full MSA when the input file is a '
							 'TSV file containing allelic profiles (this is already '
							 'the default when the input is a path to a folder with '
							 'FASTA files).')

	parser.add_argument('--gaps', type=str, choices=['ignore', 'exclude'],
						required=False, default='exclude', dest='gaps',
						help='How to treat gaps when determining the reduced MSA '
							 'for the variable positions. The default value, '
							 '"exclude", removes variable positions if any of the '
							 'aligned sequences contain a gap. The "ignore" option '
							 'allows to consider variable positions that include '
							 'gaps in some sequences as long as other sequences '
							 'include variable non-gap characters. The character '
							 'used to represent gaps is "-".')

	parser.add_argument('--ambiguous', type=str, choices=['ignore', 'exclude'],
						required=False, default='exclude', dest='ambiguous',
						help='How to treat ambiguous amino acids or nucleotides '
							 'when determining the reduced MSA for the variable '
							 'positions. The default value, "exclude", removes '
							 'variable positions if any of the aligned sequences '
							 'contain an ambiguous amino acid or nucleotide. The '
							 '"ignore" option allows to consider variable positions '
							 'that include ambiguous amino acids or nucleotides in '
							 'some sequences as long as other sequences include '
							 'variable non-ambiguous characters. The characters '
							 'interpreted as ambiguous amino acids are [B, Z, X, J].'
							 ' The characters interpreted as ambiguous nucleotides are'
                             ' [R, Y, S, W, K, M, B, D, H, V, N]')

	parser.add_argument('--custom-mafft-params', type=str,
						required=False, dest='custom_mafft_params',
						help='Custom parameters to pass to MAFFT when computing '
							 'the loci MSAs. The value must be a single string with all '
							 'parameters enclosed in quotes (e.g. "--retree 1 --maxiterate 0").')

	parser.add_argument('--protein-input', action='store_true',
						required=False, dest='protein_input',
						help='Input files contain protein sequences. This option is only valid '
							 'for cases when users provide a path to a directory containing FASTA files.')

	parser.add_argument('--no-cleanup', action='store_true',
						required=False, dest='no_cleanup',
						help='Keep intermediate files with locus/file MSAs and sample MSAs '
							 'if input is a TSV file containing allelic profiles.')

	args = parser.parse_args()
	del args.ComputeMSA

	compute_msa.main(**vars(args))


@pdt.process_timer
def run_download_schema():
	"""Run the DownloadSchema module to download a schema from Chewie-NS."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py DownloadSchema --species-id <id> --schema-id <id> --download-folder <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='DownloadSchema',
									 description='Download a schema from Chewie-NS.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/DownloadSchema.html')

	parser.add_argument('DownloadSchema', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-sp', '--species-id', type=str,
						required=True, dest='species_id',
						help='The integer identifier or name of the species '
							 'that the schema is associated to in Chewie-NS.')

	parser.add_argument('-sc', '--schema-id', type=str,
						required=True, dest='schema_id',
						help='The URI, integer identifier or name of '
							 'the schema to download from Chewie-NS.')

	parser.add_argument('-o', '--download-folder', type=str,
						required=True, dest='download_folder',
						help='Output folder to which the schema will '
							 'be saved.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads). '
							 'This value is only used if it is '
							 'necessary to construct the schema locally.')

	parser.add_argument('--ns', '--nomenclature-server', type=pv.validate_ns_url,
						required=False, default='main', dest='nomenclature_server',
						help='The base URL for the Chewie-NS instance. '
							 'The default value, "main", will establish a '
							 'connection to "https://chewbbaca.online/", '
							 '"tutorial" to "https://tutorial.chewbbaca.online/" '
							 'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
							 'Users may also provide the IP address to other '
							 'Chewie-NS instances.')

	parser.add_argument('--b', '--blast-path', type=pv.check_blast,
						required=False, default='', dest='blast_path',
						help='Path to the directory that contains the '
							 'BLAST executables.')

	parser.add_argument('--d', '--date', type=str,
						required=False, default=None, dest='date',
						help='Download schema with state from specified date. '
							 'Must be in the format "Y-m-dTH:M:S".')

	parser.add_argument('--latest', action='store_true',
						required=False, dest='latest',
						help='If the compressed version that is available is '
							 'not the latest, downloads all loci FASTA files '
							 'and constructs schema locally.')

	args = parser.parse_args()
	del args.DownloadSchema

	download_schema.main(**vars(args))


@pdt.process_timer
def run_upload_schema():
	"""Run the LoadSchema module to upload a schema to Chewie-NS."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py LoadSchema --schema-directory <dir> --species-id <id> --schema-name <name> --loci-prefix <prefix> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='LoadSchema',
									 description='Upload a schema to Chewie-NS.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/LoadSchema.html')

	parser.add_argument('LoadSchema', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-i', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the directory of the schema to upload.')

	parser.add_argument('-sp', '--species-id', type=str,
						required=True, dest='species_id',
						help='The integer identifier or name of the species '
							 'that the schema will be associated to in '
							 'Chewie-NS.')

	parser.add_argument('-sn', '--schema-name', type=str,
						required=True, dest='schema_name',
						help='A brief and meaningful name that '
							 'should help understand the type and content '
							 'of the schema.')

	parser.add_argument('-lp', '--loci-prefix', type=str,
						required=True, dest='loci_prefix',
						help='Prefix included in the name of each locus of '
							 'the schema.')

	parser.add_argument('--df', '--description-file', type=str,
						required=False, dest='description_file', default=None,
						help='Path to a text file with a description '
							 'about the schema. Markdown syntax is supported '
							 'in order to offer greater customizability of '
							 'the rendered description in the Frontend. '
							 'Will default to the schema\'s name if the user '
							 'does not provide a valid path for a file.')

	parser.add_argument('--a', '--annotations', type=str,
						required=False, dest='annotations', default=None,
						help='Path to a TSV file with loci annotations. The first column has '
							 'loci identifiers (w/o .fasta extension), the second has UniProt '
							 'protein names, the third has UniProt gene names, the fourth has '
							 'UniProt URIs, the fifth has user annotations, and the sixth has '
							 'custom annotations.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, dest='cpu_cores', default=1,
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads). '
							 'This value is used to accelerate the '
							 'quality control step that checks all alleles '
							 'in the schema.')

	parser.add_argument('--ns', '--nomenclature-server', type=pv.validate_ns_url,
						required=False, default='main', dest='nomenclature_server',
						help='The base URL for the Chewie-NS instance. '
							 'The default value, "main", will establish a '
							 'connection to "https://chewbbaca.online/", '
							 '"tutorial" to "https://tutorial.chewbbaca.online/" '
							 'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
							 'Users may also provide the IP address to other '
							 'Chewie-NS instances.')

	parser.add_argument('--continue_up', required=False, action='store_true',
						dest='continue_up',
						help='Check if the schema upload was interrupted and '
							 'attempt to continue upload.')

	args = parser.parse_args()
	del args.LoadSchema

	upload_schema.main(**vars(args))


@pdt.process_timer
def run_synchronize_schema():
	"""Run the SyncSchema module to synchronize a local schema with the remote version in Chewie-NS."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py SyncSchema --schema-directory <dir> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='SyncSchema',
									 description='Synchronize a schema with its remote version in Chewie-NS.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/SyncSchema.html')

	parser.add_argument('SyncSchema', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-sc', '--schema-directory', type=str,
						required=True, dest='schema_directory',
						help='Path to the directory with the schema to be '
							 'synced.')

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process '
							 '(chewie resets to a lower value '
							 'if it is equal to or exceeds the total '
							 'number of available CPU cores/threads). '
							 'This value is only used if the process '
							 'retrieves novel alleles from the remote '
							 'schema and needs to redetermine the set '
							 'of representative alleles for the local '
							 'schema.')

	parser.add_argument('--ns', '--nomenclature-server', type=pv.validate_ns_url,
						required=False, default=None, dest='nomenclature_server',
						help='The base URL for the Chewie-NS instance. '
							 'The default option will get the base URL from the '
							 'schema\'s URI. It is also possible to specify other '
							 'options that are available in chewBBACA\'s configs, '
							 'such as: "main" will establish a connection to '
							 '"https://chewbbaca.online/", "tutorial" to '
							 '"https://tutorial.chewbbaca.online/" and "local" '
							 'to "http://127.0.0.1:5000/NS/api/" (localhost). '
							 'Users may also provide the IP address to other '
							 'Chewie-NS instances.')

	parser.add_argument('--b', '--blast-path', type=pv.check_blast,
						required=False, default='', dest='blast_path',
						help='Path to the directory that contains the '
							 'BLAST executables.')

	parser.add_argument('--submit', required=False,
						action='store_true', dest='submit',
						help='If the process should identify new alleles '
							 'in the local schema and send them to the '
							 'Chewie-NS instance. (only authorized users can submit '
							 'new alleles).')

	# parser.add_argument('--update-profiles', required=False,
	#                     action='store_true', dest='update_profiles',
	#                     help='If the process should update local profiles '
	#                          'stored in the SQLite database.')

	args = parser.parse_args()
	del args.SyncSchema

	synchronize_schema.main(**vars(args))


@pdt.process_timer
def run_stats_requests():
	"""Run the NSStats module to get information about schemas in Chewie-NS."""

	def msg(name=None):
		usage_msg = 'chewBBACA.py NSStats --mode <mode> [options]'

		return usage_msg

	parser = argparse.ArgumentParser(prog='NSStats',
									 description='Retrieve basic information about the species and schemas in Chewie-NS.',
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog='Module documentation available at '
											'https://chewbbaca.readthedocs.io/en/latest/user/modules/NSStats.html')

	parser.add_argument('NSStats', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument('-m', '--mode', type=str,
						required=True, dest='mode',
						choices=['species', 'schemas'],
						help='The process can retrieve the list of species '
							 '("species" option) in Chewie-NS or the '
							 'list of schemas for a species '
							 '("schemas" option).')

	parser.add_argument('--sp', '--species-id', type=str,
						required=False, dest='species_id', default=None,
						help='The integer identifier of a '
							 'species in Chewie-NS.')

	parser.add_argument('--sc', '--schema-id', type=str,
						required=False, dest='schema_id', default=None,
						help='The integer identifier of a schema in '
							 'Chewie-NS.')

	parser.add_argument('--ns', '--nomenclature-server', type=pv.validate_ns_url,
						required=False, default='main', dest='nomenclature_server',
						help='The base URL for the Chewie-NS instance. '
							 'The default value, "main", will establish a '
							 'connection to "https://chewbbaca.online/", '
							 '"tutorial" to "https://tutorial.chewbbaca.online/" '
							 'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
							 'Users may also provide the IP address to other '
							 'Chewie-NS instances.')

	args = parser.parse_args()
	del args.NSStats

	stats_requests.main(**vars(args))


def main():

	functions_info = {'PredictGenes': ['Predict genes from a set of input genome assemblies.',
									  run_predict_genes],
					  'CreateSchema': ['Create a gene-by-gene schema based on '
									   'a set of genome assemblies or coding sequences.',
									   run_create_schema],
					  'AlleleCall': ['Determine the allelic profiles of a set of '
									 'bacterial genomes based on a schema.',
									 run_allele_call],
					  'SchemaEvaluator': ['Build an interactive report for schema evaluation.',
										  run_evaluate_schema],
					  'AlleleCallEvaluator': ['Build an interactive report for allele calling results evaluation.',
											  run_evaluate_calls],
					  'ExtractCgMLST': ['Determines the set of '
										'loci that constitute the '
										'core genome based on loci '
										'presence thresholds.',
										run_determine_cgmlst],
					  'SubsetResults': ['Subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers.',
										run_subset_results],
					  'PrepExternalSchema': ['Adapt an external schema to be used with chewBBACA.',
											 run_adapt_schema],
					  'MergeResults': ['Merge results files created by chewBBACA.',
									   run_merge_results],
					  'HashProfiles': ['Hash allelic profiles.',
					  					run_hash_profiles],
					  'GetAlleles': ['Create FASTA files containing the alleles identified by the AlleleCall module.',
					  				 run_get_alleles],
					  'UniprotFinder': ['Retrieve annotations for loci in a schema.',
										run_annotate_schema],
					  'ComputeDistances': ['Compute pairwise distances based on allele calling results.',
										  run_compute_distances],
					  'ComputeMSA': ['Compute a Multiple Sequence Alignment based on allele calling results.',
									 run_compute_msa],
					  'DownloadSchema': ['Download a schema from Chewie-NS.',
										 run_download_schema],
					  'LoadSchema': ['Upload a schema to Chewie-NS.',
									 run_upload_schema],
					  'SyncSchema': ['Synchronize a schema with its remote version '
									 'in Chewie-NS.',
									 run_synchronize_schema],
					  'NSStats': ['Retrieve basic information about the species '
								  'and schemas in Chewie-NS.',
								  run_stats_requests]}

	print(f'chewBBACA version: {__version__}')
	version_triggers = ['-v', '--v', '-version', '--version']
	if len(sys.argv) > 1 and sys.argv[1] in version_triggers:
		# Exit after printing version
		sys.exit(0)

	print(f'Authors: {ct.AUTHORS}')
	print(f'Github: {ct.REPOSITORY}')
	print(f'Documentation: {ct.DOCUMENTATION}')
	print(f'Contacts: {ct.CONTACTS}\n')

	# Display help message if selected process is not valid
	help_triggers = ['-h', '--h', '-help', '--help']
	if len(sys.argv) == 1 or sys.argv[1] not in functions_info or sys.argv[1] in help_triggers:
		exit_code = 0
		# Detect if user passed module name that does not exist
		if len(sys.argv) > 1 and sys.argv[1] not in help_triggers:
			print(f'No module named {sys.argv[1]}.\n')
			exit_code = 1
		print('USAGE: chewBBACA.py [module] -h, --help\n')
		print('Select one of the following modules:')
		for f in functions_info:
			print('{0}: {1}'.format(f, functions_info[f][0]))
		sys.exit(exit_code)

	# Check python version
	python_version = pv.validate_python_version()

	# Trigger module help message if no arguments are provided
	if len(sys.argv) == 2 and sys.argv[1] in functions_info:
		sys.argv.append('-h')

	process = sys.argv[1]
	functions_info[process][1]()


if __name__ == "__main__":

	main()
