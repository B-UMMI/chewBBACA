#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This is the main script of the chewBBACA suite. It parses the options and 
arguments provided through the command line and calls the specified module.
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
	from CHEWBBACA_NS import (download_schema, upload_schema,
							  synchronize_schema, stats_requests)
	from utils import (process_datetime as pdt,
					   constants as ct,
					   parameters_validation as pv,
					   file_operations as fo,
					   iterables_manipulation as im,
					   pyrodigal_gene_prediction as pgp)
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
	from CHEWBBACA.CHEWBBACA_NS import (download_schema, upload_schema,
										synchronize_schema, stats_requests)
	from CHEWBBACA.utils import (process_datetime as pdt,
								 constants as ct,
								 parameters_validation as pv,
								 file_operations as fo,
								 iterables_manipulation as im,
								 pyrodigal_gene_prediction as pgp)


@pdt.process_timer
def run_predict_genes():
	"""Run the PredictGenes module to predict genes from a set of input genome assemblies."""

	def msg(name=None):
		usage_msg = "chewBBACA.py PredictGenes --input-files <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="PredictGenes",
									 description="Predict genes from a set of input genome assemblies.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {ct.PredictGenesDocs}")

	parser.add_argument('PredictGenes', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-files", type=ct.ARGUMENT_TYPES[ct.INPUT_FILES_ARGNAME],
						required=False, dest=ct.INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the input files or to a file "
							 "with a list of full paths to the input files, one per line. Input "
							 "files must be in FASTA.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store the "
							 "files with the predicted CDSs.")

	parser.add_argument("-gp", "--gene-predictor", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTOR_ARGNAME],
					 	required=False, dest=ct.GENE_PREDICTOR_ARGNAME,
						help="Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")

	parser.add_argument("--gpa", "--gene-prediction-arguments", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTION_STR_ARGNAME],
					    nargs="+", required=False, dest=ct.GENE_PREDICTION_STR_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When providing "
							 "genome assemblies in FASTA format, the list of arguments for each parameter "
							 "used to configure the gene prediction can be passed as the long format of "
							 "the parameter name followed by the argument value (e.g., pyrodigal-training"
							 "-file=/path/to/file).")

	parser.add_argument("--t", "--translation-table", type=ct.ARGUMENT_TYPES[ct.GENETIC_CODE_ARGNAME],
						required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used for gene prediction. This value is ignored if a valid "
							 "training file is passed to `--ptf`, `--training-file`.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores that will be used to run the process (chewie resets "
							 "to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores).")

	args = parser.parse_args()
	# Use Pydantic model to validate argument values
	args = pv.PredictGenesValidator(**vars(args))

	# Exit if user only requested to create training file
	if args.validated_gene_prediction_arguments.just_training:
		sys.exit(ct.JUST_TRAINING)

	sys.exit()

	# Predict CDSs
	predict_genes.main(args.input_files, args.output_directory, args.gene_predictor, gene_predictor_parameters, args.cpu_cores)


@pdt.process_timer
def run_create_schema():
	"""Run the CreateSchema module to create a schema seed."""

	def msg(name=None):
		usage_msg = "chewBBACA.py CreateSchema --input-files <path> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="CreateSchema",
									 description="Create a schema seed.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="It is strongly advised to provide a training file to create a schema. "
											f"Module documentation available at {ct.CreateSchemaDocs}")

	parser.add_argument('CreateSchema', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-files", type=ct.ARGUMENT_TYPES[ct.INPUT_FILES_ARGNAME],
						required=True, dest=ct.INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the input FASTA files or to a file"
							 " with a list of full paths to FASTA files, one per line.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Output directory where the process will store intermediate files and "
							 "create the schema's directory.")

	parser.add_argument("--n", "--schema-name", type=ct.ARGUMENT_TYPES[ct.SCHEMA_NAME_ARGNAME],
						required=False, dest=ct.SCHEMA_NAME_ARGNAME,
						help="Name given to the schema folder.")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ARGUMENT_TYPES[ct.BLAST_SCORE_RATIO_ARGNAME],
						required=False, dest=ct.BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp "
							 "alignment and aligned sequences with a BSR >= than the defined value "
							 "are considered to be alleles of the same gene.")

	parser.add_argument("--l", "--minimum-length", type=ct.ARGUMENT_TYPES[ct.MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value. Predicted coding sequences (CDSs) "
							 "shorter than this value are excluded.")

	parser.add_argument("--t", "--translation-table", type=ct.ARGUMENT_TYPES[ct.GENETIC_CODE_ARGNAME],
						required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used to predict genes and to translate coding DNA "
							 "sequences (CDSs). This value is ignored if a valid training file "
							 "is passed to `--gpa`.")

	parser.add_argument("--st", "--size-threshold", type=ct.ARGUMENT_TYPES[ct.SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.SIZE_THRESHOLD_ARGNAME,
						help="Coding sequence (CDS) size variation threshold. Added to the "
							 "schema's config file to identify alleles with a size that deviates "
							 "from the locus length mode during the allele calling process.")

	parser.add_argument("-gp", "--gene-predictor", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTOR_ARGNAME],
					 	required=False, dest=ct.GENE_PREDICTOR_ARGNAME,
						help="Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")

	parser.add_argument("--gpa", "--gene-prediction-arguments", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTION_STR_ARGNAME],
					    nargs="+", required=False, dest=ct.GENE_PREDICTION_STR_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When "
							 "providing genome assemblies in FASTA format, the list of arguments "
							 "for each parameter used to configure the gene prediction can be "
							 "passed as the long format of the parameter name followed by the "
							 "argument value (e.g., pyrodigal-training-file=/path/to/file).")

	parser.add_argument("--cp", "--clustering-parameters", type=ct.ARGUMENT_TYPES[ct.CLUSTERING_STR_ARGNAME],
					 	nargs="+", required=False, dest=ct.CLUSTERING_STR_ARGNAME,
						help="List of arguments passed to configure the clustering of predicted "
							 "coding sequences (CDSs). The list of arguments for each parameter "
							 "used to configure the clustering can be passed as the long format "
							 "of the parameter name followed by the argument value (e.g., word-size=).")

	parser.add_argument("--b", "--blast-path", type=ct.ARGUMENT_TYPES[ct.BLAST_PATH_ARGNAME],
						required=False, dest=ct.BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--cds", "--cds-input", action="store_true",
					 	required=False, dest=ct.CDS_INPUT_ARGNAME,
						help="If provided, chewBBACA skips the gene prediction step and "
							 "assumes the input FASTA files contain coding sequences.")

	parser.add_argument("--no-cds-renaming", action="store_true",
					 	required=False, dest=ct.NO_CDS_RENAMING_ARGNAME,
						help="Do not rename the sequence/CDS identifiers when using the `--cds` "
							 "option. Provide this parameter when the input FASTA files containing "
							 "CDSs were generated by the PredictGenes module or if you are sure that "
							 "the CDS identifiers conform to the format used by chewBBACA (the input "
							 "file basename and an integer joined by "_").")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores).")

	parser.add_argument("--no-cleanup", action='store_true',
					 	required=False, dest=ct.NO_CLEANUP_ARGNAME,
						help="If provided, intermediate files generated during process execution "
							 "are not deleted at the end.")

	args = parser.parse_args()
	# Use Pydantic model to validate argument values
	args = pv.CreateSchemaValidator(**vars(args))
	print(args)

	sys.exit()

	# Run the CreateSchema process
	nloci = create_schema.main(**vars(args))
	print(f'Created schema seed with {nloci} loci.')


@pdt.process_timer
def run_allele_call():
	"""Run the AlleleCall module to perform allele calling."""

	def msg(name=None):
		usage_msg = "chewBBACA.py AlleleCall --input-files <path> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="AlleleCall",
									 description="Determine the allelic profiles of a set of genomes.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="It is strongly advised to perform allele calling with the schema's "
									 		"parameters to ensure the consistency of the results. "
											f"Module documentation available at {ct.AlleleCallDocs}")

	parser.add_argument("AlleleCall", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-files", type=ct.ARGUMENT_TYPES[ct.INPUT_FILES_ARGNAME],
						required=True, dest=ct.INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the input FASTA files or to "
							 "a file with a list of full paths to FASTA files, one per line.")

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema directory. The schema directory contains the "
							 "loci FASTA files and a folder named `short` that contains the FASTA "
							 "files with the loci representative alleles.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Output directory where the process will store intermediate files "
							 "and allele calling results (will create a subdirectory named "
							 "`results_<TIMESTAMP>` if the path passed by the user already exists).")

	parser.add_argument("--gl", "--genes-list", type=ct.ARGUMENT_TYPES[ct.LOCI_LIST_ARGNAME],
						required=False, dest=ct.LOCI_LIST_ARGNAME,
						help="Path to a file with the list of genes/loci to perform allele "
							 "calling. The file must include the full paths to the loci FASTA "
							 "files or the loci IDs, one per line. The process will perform "
							 "allele calling only for the subset of genes provided in the file.")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ARGUMENT_TYPES[ct.BLAST_SCORE_RATIO_ARGNAME],
						required=False, dest=ct.BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp "
							 "alignment and aligned sequences with a BSR >= than the defined value "
							 "are considered to be alleles of the same gene.")

	parser.add_argument("--l", "--minimum-length", type=ct.ARGUMENT_TYPES[ct.MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value. Predicted coding sequences (CDSs) "
							 "shorter than this value are excluded.")

	parser.add_argument("--t", "--translation-table", type=ct.ARGUMENT_TYPES[ct.GENETIC_CODE_ARGNAME],
						required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used to predict genes and to translate coding DNA "
							 "sequences (CDSs). This value will be ignored if a training file "
							 "is used.")

	parser.add_argument("--st", "--size-threshold", type=ct.ARGUMENT_TYPES[ct.SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.SIZE_THRESHOLD_ARGNAME,
						help="Coding sequence (CDS) size variation threshold. At the default "
							 "value of 0.2, CDSs with a size that deviates +-20 percent from "
							 "the locus length mode are classified as ASM/ALM.")

	parser.add_argument("-gp", "--gene-predictor", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTOR_ARGNAME],
					 	required=False, dest=ct.GENE_PREDICTION_STR_ARGNAME,
						help="Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")	

	parser.add_argument("--gpa", "--gene-prediction-arguments", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTION_STR_ARGNAME],
					 	nargs="+", required=False, dest=ct.GENE_PREDICTION_STR_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When "
							 "providing genome assemblies in FASTA format, the list of arguments "
							 "for each parameter used to configure the gene prediction can be "
							 "passed as the long format of the parameter name followed by the "
							 "argument value (e.g., pyrodigal-training-file=/path/to/file).")

	parser.add_argument("--cp", "--clustering-parameters", type=ct.ARGUMENT_TYPES[ct.CLUSTERING_STR_ARGNAME],
					 	nargs="+", required=False, dest=ct.CLUSTERING_STR_ARGNAME,
						help="List of arguments passed to configure the clustering of predicted "
							 "coding sequences (CDSs). The list of arguments for each parameter "
							 "used to configure the clustering can be passed as the long format "
							 "of the parameter name followed by the argument value (e.g., word-size=).")

	parser.add_argument("--b", "--blast-path", type=ct.ARGUMENT_TYPES[ct.BLAST_PATH_ARGNAME],
						required=False, dest=ct.BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--cds", "--cds-input", action='store_true',
						required=False, dest=ct.CDS_INPUT_ARGNAME,
						help="If provided, chewBBACA skips the gene prediction step and "
							 "assumes the input FASTA files contain coding sequences "
							 "(one FASTA file per strain).")

	parser.add_argument("--no-inferred", action="store_true",
						required=False, dest=ct.NO_INFERRED_ARGNAME,
						help="If provided, the process will not add the sequences of "
							 "inferred alleles (INF) to the schema. Allelic profiles "
							 "will still include the allele identifiers attributed to "
							 "the inferred alleles. Use this parameter if the schema "
							 "is being accessed by multiple processes/users simultaneously.")

	parser.add_argument("--output-unclassified", action="store_true",
						required=False, dest=ct.OUTPUT_UNCLASSIFIED_ARGNAME,
						help="Create a Fasta file with the coding sequences (CDSs) that "
							 "were not classified.")

	parser.add_argument("--output-missing", action="store_true",
					 	required=False, dest=ct.OUTPUT_MISSING_ARGNAME,
						help="Create a Fasta file with coding sequences (CDSs) classified "
							 "as NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5 and LOTSC.")

	parser.add_argument("--output-novel", action='store_true',
						required=False, dest=ct.OUTPUT_NOVEL_ARGNAME,
						help="Create a Fasta file with the novel alleles inferred during "
							 "allele calling. The sequence headers include the locus and "
							 "allele identifiers attributed by chewBBACA based on the "
							 "allele calling results.")

	parser.add_argument("--output-masked", action="store_true",
						required=False, dest=ct.OUTPUT_MASKED_ARGNAME,
						help="Create a TSV file with the masked allelic profiles. The "
							 "masking process removes the `INF-` prefix from inferred "
							 "alleles and substitutes all special classes (NIPH, NIPHEM, "
							 "ASM, ALM, PLOT3, PLOT5, LOTSC, PAMA) with `0`.")

	parser.add_argument('--no-cds-renaming', action='store_true',
					 	required=False, dest=ct.NO_CDS_RENAMING_ARGNAME,
						help="Do not rename the sequence/CDS identifiers when using "
							 "the `--cds` option. Provide this parameter when the input "
							 "FASTA files containing CDSs were generated by the PredictGenes "
							 "module or if you are sure that the CDS identifiers conform "
							 "to the format used by chewBBACA (the input file basename "
							 "and an integer joined by `_`).")

	parser.add_argument("--force-continue", action='store_true',
						required=False, dest=ct.FORCE_CONTINUE_ARGNAME,
						help="If provided, chewie will not warn users and ask for "
							 "permission to continue if any of the provided argument "
							 "values does not match the values in the config file.")

	parser.add_argument("--mode", type=int,
					 	required=False, dest=ct.ALLELECALL_MODE_ARGNAME,
						help="Execution mode (1: only exact matches at DNA level; 2: "
							 "exact matches at DNA and Protein level; 3: exact matches "
							 "and minimizer-based clustering to find similar alleles "
							 "based on BSR+0.1; 4: run the full process to find exact "
							 "matches and similar matches based on BSR value, including "
							 "the determination of new representative alleles to add to "
							 "the schema).")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores).")

	parser.add_argument("--no-cleanup", action='store_true',
						required=False, dest=ct.NO_CLEANUP_ARGNAME,
						help="If provided, intermediate files generated during process execution "
							 "are not removed at the end.")

	args = parser.parse_args()
	# Use Pydantic model to validate argument values
	args = pv.AlleleCallValidator(**vars(args))
	print(args)

	sys.exit()

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


@pdt.process_timer
def run_evaluate_schema():
	"""Run the SchemaEvaluator module to evaluate a typing schema."""

	def msg(name=None):
		usage_msg = "chewBBACA.py SchemaEvaluator --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="SchemaEvaluator",
									 description="Build an interactive report for schema evaluation.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="The module can evaluate schemas created with chewBBACA or from external platforms. "
									 		"Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/SchemaEvaluator.html")

	parser.add_argument("SchemaEvaluator", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the report HTML files will be created.")

	parser.add_argument("--ll", "--loci-list", type=ct.ARGUMENT_TYPES[ct.LOCI_LIST_ARGNAME],
						required=False, dest=ct.LOCI_LIST_ARGNAME,
						help="Path to a file with the list of loci in the schema that the process "
							 "should analyse (one per line, full paths or loci IDs).")

	parser.add_argument("-a", "--annotations", type=ct.ARGUMENT_TYPES[ct.ANNOTATIONS_ARGNAME],
					 	required=False, dest=ct.ANNOTATIONS_ARGNAME,
						help="Path to the TSV file created by the UniprotFinder module. The "
							 "annotation data is included in a table component.")

	parser.add_argument("--ta", "--translation-table", type=ct.ARGUMENT_TYPES[ct.GENETIC_CODE_ARGNAME],
						required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used to translate coding sequences (CDSs).")

	parser.add_argument("--st", "--size-threshold", type=ct.ARGUMENT_TYPES[ct.SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.SIZE_THRESHOLD_ARGNAME,
						help="Coding sequence (CDS) size variation threshold. The module identifies "
							 "the alleles with size that deviates from the locus length mode +- the "
							 "size threshold.")

	parser.add_argument("--ml", "--minimum-length", type=ct.ARGUMENT_TYPES[ct.MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value. The module identifies alleles shorter "
							 "than this value.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--loci-reports", action="store_true",
					 	required=False, dest=ct.LOCI_REPORTS_ARGNAME,
						help="Create a detailed report page for each locus. The locus report "
							 "includes components with relevant data and analysis results, such "
							 "as allele diversity charts, a MSA for the alignment of the distinct "
							 "translated alleles and a tree drawn with Phylocanvas based on the "
							 "MAFFT guide tree.")

	parser.add_argument("--light", action="store_true",
					 	required=False, dest=ct.LIGTH_ARGNAME,
						help="Skips MSA computation with MAFFT and does not add the Phylogenetic "
							 "Tree and MSA components to the loci reports.")

	parser.add_argument("--add-sequences", action="store_true",
					 	required=False, dest=ct.ADD_SEQUENCES_ARGNAME,
						help="Adds Code Editor components with the DNA and Protein sequences to "
							 "the loci reports. The Code Editor is in readonly mode (allows to "
							 "search for and copy text).")

	args = parser.parse_args()

	# Use Pydantic model to validate argument values
	args = pv.SchemaEvaluatorValidator(**vars(args))
	print(args)

	sys.exit()


	evaluate_schema.main(**vars(args))


@pdt.process_timer
def run_evaluate_calls():
	"""Run the AlleleCallEvaluator module to evaluate allele calling results."""

	def msg(name=None):
		usage_msg = "chewBBACA.py AlleleCallEvaluator --input-files <dir> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="AlleleCallEvaluator",
									 description="Build an interactive report for allele calling results evaluation.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/AlleleCallEvaluator.html")

	parser.add_argument("AlleleCallEvaluator", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-files", type=ct.ARGUMENT_TYPES[ct.INPUT_FILES_ARGNAME],
					 	required=True, dest=ct.INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the allele calling results "
							 "generated by the AlleleCall module.")

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the module will store intermediate "
							 "files and create the report HTML files.")

	parser.add_argument("-a", "--annotations", type=ct.ARGUMENT_TYPES[ct.ANNOTATIONS_ARGNAME],
					 	required=False, dest=ct.ANNOTATIONS_ARGNAME,
						help="Path to the TSV file created by the UniprotFinder module.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--light", action="store_true",
					 	required=False, dest=ct.LIGTH_ARGNAME,
						help="Do not compute the presence-absence matrix, the distance matrix "
							 "and the Neighbor-Joining tree.")

	parser.add_argument("--no-pa", action="store_true",
					 	required=False, dest=ct.NO_PA_ARGNAME,
						help="Do not compute the presence-absence matrix.")

	parser.add_argument("--no-dm", action='store_true',
					 	required=False, dest=ct.NO_DM_ARGNAME,
						help="Do not compute the distance matrix.")

	parser.add_argument("--no-tree", action="store_true",
					 	required=False, dest=ct.NO_TREE_ARGNAME,
						help="Do not compute the Neighbor-Joining tree.")

	parser.add_argument("--cg-alignment", action="store_true",
					 	required=False, dest=ct.CG_ALIGNMENT_ARGNAME,
						help="Compute the MSA of the core genome loci, even if `--no-tree` "
							 "is provided.")

	args = parser.parse_args()

	# Use Pydantic model to validate argument values
	args = pv.AlleleCallEvaluatorValidator(**vars(args))
	print(args)

	sys.exit()

	evaluate_calls.main(**vars(args))


@pdt.process_timer
def run_determine_cgmlst():
	"""Run the ExtractCgMLST module to determine the set of core loci based on allele calling results."""

	def msg(name=None):
		usage_msg = "chewBBACA.py ExtractCgMLST --input-file <file> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="ExtractCgMLST",
									 description="Determine the set of core loci based on allele calling results.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/ExtractCgMLST.html")

	parser.add_argument("ExtractCgMLST", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-files", type=ct.ARGUMENT_TYPES[ct.RESULTS_FILES_ARGNAME],
						required=True, dest=ct.RESULTS_FILES_ARGNAME,
						help="Path to the TSV file that contains the allelic profiles determined by "
							 "the AlleleCall module.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the directory where the process will store the output files.")

	parser.add_argument("--t", "--threshold", type=float, nargs="+",
						required=False, dest=ct.THRESHOLD_ARGNAME,
						help="Loci/genes that constitute the core genome must be in a proportion of "
							 "genomes that is at least equal to this value. Provide multiple values "
							 "to compute the core genome for multiple threshold values.")

	parser.add_argument("--s", "--step", type=ct.ARGUMENT_TYPES[ct.STEP_ARGNAME],
					 	required=False, dest=ct.STEP_ARGNAME,
						help="The allele calling results are processed iteratively to evaluate the "
							 "impact of adding subsets of the results in computing the core genome. "
							 "The step value controls the number of allelic profiles added in each "
							 "iteration until all profiles are included.")

	parser.add_argument("--ca", "--compute-accessory", action="store_true",
					 	required=False, dest=ct.COMPUTE_ACCESSORY_ARGNAME,
						help="Determine the set of accessory loci. The accessory genome corresponds "
							 "to all the loci not included in the core genome. The accessory genome "
							 "is determined for each core genome threshold.")

	parser.add_argument("--ra", "--rarefaction-analysis", action="store_true",
					 	required=False, dest=ct.RAREFACTION_ANALYSIS_ARGNAME,
						help="Perform rarefaction analysis to evaluate the stability of the core genome "
							 "are randomly selected and the number of core loci (for the core genome) or "
							 "total loci (for the pangenome) is computed for each subset from 1 to the "
							 "total number of samples. The average number of core loci (for the core genome) "
							 "or total loci (for the pangenome) is computed for each sample size subset "
							 "using the values of all permutations to plot the rarefaction curve. The "
							 "rarefaction analysis is performed for each threshold. A power law model "
							 "is fitted to the rarefaction curve to estimate the stability of the core genome "
							 "and the openness of the pangenome.")

	parser.add_argument("--pn", "--permutation-number", type=ct.ARGUMENT_TYPES[ct.PERMUTATION_NUMBER_ARGNAME],
					 	required=False, dest=ct.PERMUTATION_NUMBER_ARGNAME,
						help="Number of permutations for the rarefaction analysis. The rarefaction analysis "
							 "is repeated a number of times equal to the value provided to this parameter.")

	parser.add_argument("--ps", "--permutation-samples", type=ct.ARGUMENT_TYPES[ct.PERMUTATION_SAMPLES_ARGNAME],
					 	required=False, dest=ct.PERMUTATION_SAMPLES_ARGNAME,
						help="Number of samples randomly selected for each permutation. All samples will be "
							 "used if this value is not provided.")

	parser.add_argument("--el", "--exclude-loci", type=ct.ARGUMENT_TYPES[ct.EXCLUDE_LOCI_ARGNAME],
						required=False, dest=ct.EXCLUDE_LOCI_ARGNAME,
						help="Path to a file with a list of loci identifiers to exclude from the analysis "
							 "(one locus identifier per line).")

	parser.add_argument("--eg", "--exclude-genomes", type=ct.ARGUMENT_TYPES[ct.EXCLUDE_GENOMES_ARGNAME],
						required=False, dest=ct.EXCLUDE_GENOMES_ARGNAME,
						help="Path to a file with a list of genome identifiers to exclude from the analysis "
							 "(one genome identifier per line).")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Maximum number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads).")

	args = parser.parse_args()
	args = pv.ExtractCgMLSTValidator(**vars(args))
	print(args)

	sys.exit(0)

	determine_cgmlst.main(**vars(args))


@pdt.process_timer
def run_subset_results():
	"""Run the SubsetResults module to subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers."""

	def msg(name=None):
		usage_msg = "chewBBACA.py SubsetResults --input-file <file> --loci-list <file> --samples-list <file> --output-file <file> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="SubsetResults",
									 description="Subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/SubsetResults.html")

	parser.add_argument("SubsetResults", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-files", type=ct.ARGUMENT_TYPES[ct.RESULTS_FILES_ARGNAME],
						required=True, dest=ct.RESULTS_FILES_ARGNAME,
						help="Path to the directory containing the files to be subsetted.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("-l", "--loci-list", type=ct.ARGUMENT_TYPES[ct.LOCI_LIST_ARGNAME],
						required=False, dest=ct.LOCI_LIST_ARGNAME,
						help="Path to a TXT/TSV file containing a list of loci to select, one locus identifier "
							 "per line. If the file contains multiple columns, the loci identifiers must be in "
							 "the first column.")

	parser.add_argument("-s", "--sample-list", type=ct.ARGUMENT_TYPES[ct.SAMPLE_LIST_ARGNAME],
						required=False, dest=ct.SAMPLE_LIST_ARGNAME,
						help="Path to a TXT/TSV file containing a list of samples to select, one sample identifier "
							 "per line. If the file contains multiple columns, the sample identifiers must be in "
							 "the first column.")

	parser.add_argument("--inverse-loci", action="store_true",
						required=False, dest=ct.INVERSE_LOCI_ARGNAME,
						help="If provided, the process will select the loci that are not in the input loci list.")

	parser.add_argument("--inverse-samples", action="store_true",
						required=False, dest=ct.INVERSE_SAMPLES_ARGNAME,
						help="If provided, the process will select the samples that are not in the input samples list.")

	args = parser.parse_args()
	args = pv.SubsetResultsValidator(**vars(args))
	print(args)

	sys.exit(0)

	subset_results.main(**vars(args))


@pdt.process_timer
def run_merge_results():
	"""Run the MergeResults module to merge results files created by chewBBACA."""

	def msg(name=None):
		usage_msg = "chewBBACA.py MergeResults --input-directories <dir> <dir> ... --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="MergeResults",
									 description="Merge results files created by chewBBACA.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/MergeResults.html")

	parser.add_argument("MergeResults", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-files", nargs="+", type=ct.ARGUMENT_TYPES[ct.RESULTS_FILES_ARGNAME],
						required=True, dest=ct.RESULTS_FILES_ARGNAME,
						help="Paths to the directories containing the results files created by chewBBACA. The "
							 "results must have been determined with the same schema and share all the loci "
							 "or a subset of the loci if using the --common parameter.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--common", action="store_true",
						required=False, dest=ct.COMMON_ARGNAME,
						help="Merge the results based on the subset of loci shared between all inputs.")

	args = parser.parse_args()
	args = pv.MergeResults(**vars(args))
	print(args)

	sys.exit(0)

	merge_results.main(**vars(args))


@pdt.process_timer
def run_hash_profiles():
	"""Run the HashProfiles module to hash allelic profiles."""

	def msg(name=None):
		usage_msg = "chewBBACA.py HashProfiles --input-file <file> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="HashProfiles",
									 description="Hash allelic profiles.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/HashProfiles.html")

	parser.add_argument("HashProfiles", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-a", "--allelic-profiles", type=ct.ARGUMENT_TYPES[ct.ALLELIC_PROFILES_ARGNAME],
						required=True, dest=ct.ALLELIC_PROFILES_ARGNAME,
						help="Path to the TSV file that contains the allelic profiles determined by "
							 "the AlleleCall module.")

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory to get the allele sequences and compute the hashes.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--hash-type", type=ct.ARGUMENT_TYPES[ct.HASH_TYPE_ARGNAME],
					 	required=False, dest=ct.HASH_TYPE_ARGNAME,
						help="Hashing algorithm used to hash the profiles. The hashing algorithms implemented "
							 "in the hashlib and zlib Python libraries are supported.")

	parser.add_argument("--nrows", type=ct.ARGUMENT_TYPES[ct.NROWS_ARGNAME],
					 	required=False, dest=ct.NROWS_ARGNAME,
						help="Divide the input file into chunks of this many rows to process larger files"
						 	 " more efficiently.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

	args = parser.parse_args()
	args = pv.HashProfilesValidator(**vars(args))
	print(args)

	sys.exit(0)

	hash_profiles.main(**vars(args))


@pdt.process_timer
def run_get_alleles():
	"""Run the GetAlleles module to create FASTA files containing the alleles identified by the AlleleCall module."""

	def msg(name=None):
		usage_msg = "chewBBACA.py GetAlleles --input-file <file> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="GetAlleles",
									 description="Create FASTA files containing the alleles identified by the AlleleCall module.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/GetAlleles.html")

	parser.add_argument("GetAlleles", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-a", "--allelic-profiles", type=ct.ARGUMENT_TYPES[ct.ALLELIC_PROFILES_ARGNAME],
						required=True, dest=ct.ALLELIC_PROFILES_ARGNAME,
						help="Path to the TSV file containing the allelic profiles.")

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema directory.")

	parser.add_argument("-l", "--loci-list", type=ct.ARGUMENT_TYPES[ct.LOCI_LIST_ARGNAME],
						required=False, dest=ct.LOCI_LIST_ARGNAME,
						help="Path to a file with the list of genes/loci to create FASTA files for. The "
							 "file must include the identifiers of the loci, one per line, without the "
							 ".fasta extension.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

	parser.add_argument("--distinct", action="store_true",
						required=False, dest=ct.DISTINCT_ARGNAME,
						help="Only get distinct alleles.")

	parser.add_argument("--translate", action="store_true",
						required=False, dest=ct.TRANSLATE_ARGNAME,
						help="Create FASTA files with the translated alleles.")

	parser.add_argument("--ta", "--translation-table", type=ct.GENETIC_CODE_ARGNAME,
					 	required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used to translate coding DNA sequences (CDSs). If no value"
							 " is specified, the process tries to get the value stored in the schema "
							 "config file. If the schema does not include a config file, the process "
							 "uses the default translation table (11).")

	args = parser.parse_args()
	args = pv.GetAllelesValidator(**vars(args))
	print(args)

	sys.exit(0)

	get_alleles.main(**vars(args))


@pdt.process_timer
def run_adapt_schema():
	"""Run the PrepExternalSchema module to adapt a typing schema."""

	def msg(name=None):
		usage_msg = "chewBBACA.py PrepExternalSchema --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="PrepExternalSchema",
									 description="Adapt an external schema to be used with chewBBACA.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/PrepExternalSchema.html")

	parser.add_argument("PrepExternalSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the directory of the schema to adapt. The schema must contain one "
							 "FASTA file per gene/locus.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the adapted schema will be created.")

	parser.add_argument("--l", "--loci-list", type=ct.ARGUMENT_TYPES[ct.LOCI_LIST_ARGNAME],
						required=False, dest=ct.LOCI_LIST_ARGNAME,
						help="Path to a file with the list of loci in the schema that the process should "
							 "adapt (one per line, full paths or loci IDs).")

	parser.add_argument("-gp", "--gene-predictor", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTOR_ARGNAME],
							required=False, dest=ct.GENE_PREDICTOR_ARGNAME,
							help="Specify which gene prediction software to use. Default is Pyrodigal "
								 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
								 "for prokaryotic and eukaryotic genomes.")
	
	parser.add_argument("--gpa", "--gene-prediction-arguments", type=ct.ARGUMENT_TYPES[ct.GENE_PREDICTION_STR_ARGNAME],
						nargs="+", required=False, dest=ct.GENE_PREDICTION_STR_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When providing "
							 "genome assemblies in FASTA format, the list of arguments for each parameter "
							 "used to configure the gene prediction can be passed as the long format of "
							 "the parameter name followed by the argument value (e.g., pyrodigal-training"
							 "-file=/path/to/file).")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ARGUMENT_TYPES[ct.BLAST_SCORE_RATIO_ARGNAME],
						required=False, dest=ct.BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio (BSR) value. The process selects representative alleles "
							 "for each locus based on this value. Representative alleles are selected "
							 "until all alleles in a locus align against one of the representatives "
							 "with a BSR >= than the specified value.")

	parser.add_argument("--l", "--minimum-length", type=ct.ARGUMENT_TYPES[ct.MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value stored in the schema config file. The "
							 "schema adaptation process will only discard sequences smaller than "
							 "this value if the --size-filter parameter is provided.")

	parser.add_argument("--t", "--translation-table", type=ct.ARGUMENT_TYPES[ct.GENETIC_CODE_ARGNAME],
					 	required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used for allele translation. This value is ignored if "
							 "a valid training file is passed to `--ptf`, `--training-file`.")

	parser.add_argument("--st", "--size-threshold", type=ct.ARGUMENT_TYPES[ct.SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.SIZE_THRESHOLD_ARGNAME,
						help="Allele size variation threshold value stored in the schema config "
							 "file. The schema adaptation process will only discard alleles with "
							 "a size that deviates from the locus length mode +- the size theshold "
							 "value if the --size-filter parameter is provided.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--b", "--blast-path", type=ct.ARGUMENT_TYPES[ct.BLAST_PATH_ARGNAME],
						required=False, dest=ct.BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--size-filter", action="store_true",
						required=False, dest=ct.SIZE_FILTER_ARGNAME,
						help="Apply the minimum length and size threshold values to filter out "
							 "alleles during schema adaptation.")

	args = parser.parse_args()
	args = pv.PrepExternalSchemaValidator(**vars(args))
	print(args)

	sys.exit(0)

	# Define output paths
	schema_path = os.path.abspath(args.output_directory)
	schema_short_path = fo.join_paths(schema_path, ['short'])
	output_dirs = [schema_path, schema_short_path]

	# Create output directories
	schema_path_exists = fo.create_directory(schema_path)
	if schema_path_exists is False:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)
	fo.create_directory(schema_short_path)

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


@pdt.process_timer
def run_annotate_schema():
	"""Run the UniprotFinder module to annotate loci in a schema."""

	def msg(name=None):
		usage_msg = "chewBBACA.py UniprotFinder --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="UniprotFinder",
									 description="Retrieve annotations for loci in a schema.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/UniprotFinder.html")

	parser.add_argument("UniprotFinder", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store intermediate "
							 "files and save the final TSV file with the loci annotations.")

	parser.add_argument("--l", "--loci-list", type=ct.ARGUMENT_TYPES[ct.LOCI_LIST_ARGNAME],
						required=False, dest=ct.LOCI_LIST_ARGNAME,
						help="Path to a file with the list of loci in the schema that the process "
							 "should find annotations for (one per line, full paths or loci IDs).")

	parser.add_argument("--t", "--protein-table", type=ct.ARGUMENT_TYPES[ct.PROTEIN_TABLE_ARGNAME],
						required=False, dest=ct.PROTEIN_TABLE_ARGNAME,
						help="Path to the TSV file with coding sequence (CDS) coordinate data, "
							 "`cds_coordinates.tsv`, created by the CreateSchema process.")

	parser.add_argument("--bsr", type=ct.ARGUMENT_TYPES[ct.BLAST_SCORE_RATIO_ARGNAME],
					 	required=False, dest=ct.BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio value. The BSR is only used when taxa names are "
							 "provided to the --taxa parameter and local sequences are aligned "
							 "against reference proteomes downloaded from UniProt. Annotations "
							 "are selected based on a BSR >= than the specified value.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--taxa", nargs="+", type=ct.ARGUMENT_TYPES[ct.TAXA_ARGNAME],
						required=False, dest=ct.TAXA_ARGNAME,
						help="List of scientific names for a set of taxa. The process will download "
							 "reference proteomes from UniProt associated to taxa names that contain "
							 "any of the provided terms. The schema representative alleles are aligned "
							 "against the reference proteomes to assign annotations based on high-BSR "
							 "matches.")

	parser.add_argument("--pm", type=ct.ARGUMENT_TYPES[ct.PROTEOME_MATCHES_ARGNAME],
					 	required=False, dest=ct.PROTEOME_MATCHES_ARGNAME,
						help="Maximum number of proteome matches to report.")

	parser.add_argument("--no-sparql", action="store_true",
						required=False, dest=ct.NO_SPARQL_ARGNAME,
						help="Do not search for annotations through the UniProt SPARQL endpoint.")

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, dest=ct.NO_CLEANUP_ARGNAME,
						help="If provided, intermediate files generated during process execution "
							 "are not removed at the end.")

	parser.add_argument("--b", "--blast-path", type=ct.ARGUMENT_TYPES[ct.BLAST_PATH_ARGNAME],
						required=False, default='', dest=ct.BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	args = parser.parse_args()
	args = pv.UniprotFinderValidator(**vars(args))
	print(args)

	sys.exit(0)

	annotate_schema.main(**vars(args))


@pdt.process_timer
def run_compute_distances():
	"""Run the ComputeDistances module to compute pairwise distances based on allele calling results."""

	def msg(name=None):
		usage_msg = "chewBBACA.py ComputeDistances --input-file <file> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="ComputeDistances",
									 description="Compute pairwise distances based on allele calling results.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/ComputeDistances.html")

	parser.add_argument("ComputeDistances", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-a", "--allelic-profiles", type=ct.ARGUMENT_TYPES[ct.ALLELIC_PROFILES_ARGNAME],
						required=True, dest=ct.ALLELIC_PROFILES_ARGNAME,
						help="Path to a TSV file containing allelic profiles determined by the AlleleCall "
							 "module.")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store intermediate "
							 "and final results.")

	parser.add_argument("--m", "--method", type=ct.ARGUMENT_TYPES[ct.METHOD_ARGNAME],
					 	required=False, dest=ct.METHOD_ARGNAME,
						help="Distance method used to compute the distance matrix. The module supports "
							 "the hamming, jaccard, loci (number of loci not shared), and core (number of "
							 "different alleles for core loci) methods.")

	parser.add_argument('--outfmt', '--output-format', type=ct.ARGUMENT_TYPES[ct.OUTPUT_FORMAT_ARGNAME],
					 	required=False, dest=ct.OUTPUT_FORMAT_ARGNAME,
						help="Output format for the distance matrix (upper_triangular, lower_triangular, "
							 "symmetric, table).")

	parser.add_argument("--no-mask", action="store_true",
					 	required=False, dest=ct.NO_MASK_ARGNAME,
						help="Do not mask missing data when computing the distance matrix. This option "
							 "is useful when the input profiles are already masked.")

	parser.add_argument("--similarity", action="store_true",
					 	required=False, dest=ct.SIMILARITY_ARGNAME,
						help="Compute similarity values instead of distance values.")

	parser.add_argument('--cpu', '--cpu-cores', type=pv.verify_cpu_usage,
						required=False, default=1, dest='cpu_cores',
						help='Number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a '
							 'lower value if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	args = parser.parse_args()
	args = pv.ComputeDistancesValidator(**vars(args))
	print(args)
	sys.exit(0)

	compute_distances.main(**vars(args))


@pdt.process_timer
def run_compute_msa():
	"""Run the ComputeMSA module to compute a Multiple Sequence Alignment based on allele calling results."""

	def msg(name=None):
		usage_msg = "chewBBACA.py ComputeMSA --input-file <file> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="ComputeMSA",
									 description="Compute a Multiple Sequence Alignment based on allele calling results.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/ComputeMSA.html")

	parser.add_argument("ComputeMSA", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-path", type=ct.ARGUMENT_TYPES[ct.INPUT_PATH_ARGNAME],
						required=True, dest=ct.INPUT_PATH_ARGNAME,
						help="Path to a TSV file containing allelic profiles or to a folder containing "
							 "FASTA files. If a TSV file containing allelic profiles is provided, it is "
							 "necessary to provide the path to the schema to the `--schema-directory` "
							 "parameter. The module will create a FASTA file with the alleles identified "
							 "in the samples for each schema locus and compute a MSA. The loci MSAs are "
							 "joined to create the complete MSA based on the allele calling results. If "
							 "a path to a folder is provided, the module computes a MSA for each FASTA "
							 "file in the folder, but will not attempt to join the MSAs as it does not "
							 "have the sample information (in this case, it is not necessary to pass the "
							 "schema path).")

	parser.add_argument("-o", "--output-directory", type=ct.ARGUMENT_TYPES[ct.OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store intermediate and "
							 "final results.")

	parser.add_argument("-g", "--schema-directory", type=ct.ARGUMENT_TYPES[ct.SCHEMA_DIRECTORY_ARGNAME],
						required=False, dest=ct.SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema\'s directory. This parameter is only required if the "
							 "input is a TSV file with allelic profiles.")

	parser.add_argument("--dna-msa", action="store_true",
						required=False, dest=ct.DNA_MSA_ARGNAME,
						help="Converts the protein MSA back to DNA to create an additional output file "
							 "with the DNA MSA.")

	parser.add_argument("--output-variable", action="store_true",
					 	required=False, dest=ct.OUTPUT_VARIABLE_ARGNAME,
						help="Output a reduced MSA including only the variable positions. If the "
							 "`--dna-msa` parameter is provided, the process will output a reduced MSA "
							 "for both the protein and DNA MSAs.")

	parser.add_argument("--t", "--translation-table", type=ct.ARGUMENT_TYPES[ct.GENETIC_CODE_ARGNAME],
						required=False, dest=ct.GENETIC_CODE_ARGNAME,
						help="Genetic code used for sequence translation.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CPU_CORES_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

	parser.add_argument("--only-loci-msas", action="store_true",
						required=False, dest=ct.ONLY_LOCI_MSAS_ARGNAME,
						help="Do not compute the full MSA when the input file is a TSV file containing "
							 "allelic profiles (this is already the default when the input is a path to "
							 "a folder with FASTA files).")

	parser.add_argument("--gaps", type=ct.ARGUMENT_TYPES[ct.GAPS_ARGNAME],
						required=False, dest=ct.GAPS_ARGNAME,
						help="How to treat gaps when determining the reduced MSA for the variable "
							 "positions. The default value, `exclude`, removes variable positions if "
							 "any of the aligned sequences contain a gap. The `ignore` option allows "
							 "to consider variable positions that include gaps in some sequences as "
							 "long as other sequences include variable non-gap characters. The character "
							 "used to represent gaps is `-`.")

	parser.add_argument("--ambiguous", type=ct.ARGUMENT_TYPES[ct.AMBIGUOUS_ARGNAME],
						required=False, dest=ct.DEFAULT_AMBIGUOUS,
						help="How to treat ambiguous amino acids or nucleotides when determining the "
							 "reduced MSA for the variable positions. The default value, `exclude`, "
							 "removes variable positions if any of the aligned sequences contain an "
							 "ambiguous amino acid or nucleotide. The `ignore` option allows to consider "
							 "variable positions that include ambiguous amino acids or nucleotides in "
							 "some sequences as long as other sequences include variable non-ambiguous "
							 "characters. The characters interpreted as ambiguous amino acids are "
							 "[B, Z, X, J]. The characters interpreted as ambiguous nucleotides are"
                             " [R, Y, S, W, K, M, B, D, H, V, N].")

	parser.add_argument("--custom-mafft-params", type=ct.ARGUMENT_TYPES[ct.CUSTOM_MAFFT_PARAMETERS_ARGNAME],
						required=False, dest=ct.CUSTOM_MAFFT_PARAMETERS_ARGNAME,
						help="Custom parameters to pass to MAFFT when computing the loci MSAs. The "
							 "value must be a single string with all parameters enclosed in quotes "
							 "(e.g. `--retree 1 --maxiterate 0`).")

	parser.add_argument("--protein-input", action="store_true",
						required=False, dest=ct.PROTEIN_INPUT_ARGNAME,
						help="Input files contain protein sequences. This option is only valid for "
							 "cases when users provide a path to a directory containing FASTA files.")

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, dest=ct.NO_CLEANUP_ARGNAME,
						help="Keep intermediate files with locus/file MSAs and sample MSAs if input "
							 "is a TSV file containing allelic profiles.")

	args = parser.parse_args()
	args = pv.ComputeMSAValidator(**vars(args))
	print(args)

	sys.exit(0)

	compute_msa.main(**vars(args))


@pdt.process_timer
def run_download_schema():
	"""Run the DownloadSchema module to download a schema from Chewie-NS."""

	def msg(name=None):
		usage_msg = "chewBBACA.py DownloadSchema --species-id <id> --schema-id <id> --download-folder <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="DownloadSchema",
									 description="Download a schema from Chewie-NS.",
									 usage=msg(),
									 formatter_class=pv.ModifiedHelpFormatter,
									 epilog="Module documentation available at https://chewbbaca.readthedocs.io/en/latest/user/modules/DownloadSchema.html")

	parser.add_argument("DownloadSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-sp", "--species-id", type=ct.ARGUMENT_TYPES[ct.SPECIES_ID_ARGNAME],
						required=True, dest=ct.SPECIES_ID_ARGNAME,
						help="The integer identifier or name of the species that the schema is "
							 "associated to in Chewie-NS.")

	parser.add_argument("-sc", "--schema-id", type=ct.ARGUMENT_TYPES[ct.SCHEMA_ID_ARGNAME],
						required=True, dest=ct.SCHEMA_ID_ARGNAME,
						help="The URI, integer identifier or name of the schema to download from "
							 "Chewie-NS.")

	parser.add_argument("-o", "--download-folder", type=ct.ARGUMENT_TYPES[ct.DOWNLOAD_FOLDER_ARGNAME],
						required=True, dest=ct.DOWNLOAD_FOLDER_ARGNAME,
						help="Output folder to which the schema will be saved.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.ARGUMENT_TYPES[ct.CDS_INPUT_ARGNAME],
						required=False, dest=ct.CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads). This value is only used if it is "
							 "necessary to construct the schema locally.")

	parser.add_argument("--ns", "--nomenclature-server", type=ct.ARGUMENT_TYPES[ct.NOMENCLATURE_SERVER_ARGNAME],
						required=False, dest=ct.NOMENCLATURE_SERVER_ARGNAME,
						help="The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https://"
							 "tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")

	parser.add_argument("--b", "--blast-path", type=ct.ARGUMENT_TYPES[ct.BLAST_PATH_ARGNAME],
						required=False, dest=ct.BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--d", "--date", type=ct.ARGUMENT_TYPES[ct.DATE_ARGNAME],
						required=False, dest=ct.DATE_ARGNAME,
						help="Download schema with state from specified date. Must be in the format "
							 "`Y-m-dTH:M:S`.")

	parser.add_argument("--latest", action="store_true",
						required=False, dest=ct.LATEST_ARGNAME,
						help="If the compressed version that is available is not the latest, downloads "
							 "all loci FASTA files and constructs schema locally.")

	args = parser.parse_args()
	args = pv.DownloadSchemaValidator(**vars(args))
	print(args)

	sys.exit(0)

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
									 formatter_class=pv.ModifiedHelpFormatter,
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
									 formatter_class=pv.ModifiedHelpFormatter,
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
									 formatter_class=pv.ModifiedHelpFormatter,
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
					  'CreateSchema': ['Create a gene-by-gene schema based on a set of genome assemblies or coding sequences.',
									   run_create_schema],
					  'AlleleCall': ['Determine the allelic profiles of a set of bacterial genomes based on a schema.',
									 run_allele_call],
					  'SchemaEvaluator': ['Build an interactive report for schema evaluation.',
										  run_evaluate_schema],
					  'AlleleCallEvaluator': ['Build an interactive report for allele calling results evaluation.',
											  run_evaluate_calls],
					  'ExtractCgMLST': ['Determines the set of loci that constitute the core genome based on loci presence thresholds.',
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
					  'SyncSchema': ['Synchronize a schema with its remote version in Chewie-NS.',
									 run_synchronize_schema],
					  'NSStats': ['Retrieve basic information about the species and schemas in Chewie-NS.',
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
	python_uptodate, python_version = pv.validate_python_version(ct.MIN_PYTHON)
	if not python_uptodate:
		sys.exit(ct.PYTHON_VERSION.format(python_version, ct.MIN_PYTHON))

	# Trigger module help message if no arguments are provided
	if len(sys.argv) == 2 and sys.argv[1] in functions_info:
		sys.argv.append('-h')

	process = sys.argv[1]
	functions_info[process][1]()


if __name__ == "__main__":

	main()
