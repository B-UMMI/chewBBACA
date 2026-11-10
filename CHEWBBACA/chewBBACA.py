#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This is the main script of the chewBBACA suite. It parses the options and arguments 
provided through the command line to call the specified module.
"""


import os
import sys
import shlex
import shutil
import argparse

from loguru import logger
from pydantic import ValidationError

from CHEWBBACA import global_constants as gct
from CHEWBBACA import global_validation as gval


# Initialize logging
# Start with console logging in case Pydantic validation fails
logger.remove() # CLear default handlers before adding custom one
logger.add(sys.stderr, level="INFO", format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}")


# Function to initialize file logging
def setup_file_logging(output_dir):
	"""
	"""
	logfile_path = os.path.join(output_dir, "chewie.log")
	logger.add(
		logfile_path,
		level="DEBUG",
		enqueue=True,
		format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {name}:{function}:{line} - {message}"
	)

	logger.info(f"File logging initialized at: {logfile_path}")

	return logfile_path


def print_initial_header():
	"""
	"""
	print("="*len(f"Authors: {gct.AUTHORS}"))
	print(f"Authors: {gct.AUTHORS}")
	print(f"Github: {gct.REPOSITORY}")
	print(f"Documentation: {gct.DOCUMENTATION}")
	print(f"Contacts: {gct.CONTACTS}")
	print("="*len(f"Authors: {gct.AUTHORS}"))


def log_cli_command(process_name):
	"""
	"""
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[process_name]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")


def log_validated_args(validated_args):
	""""""
	logger.info("Validated arguments:")
	for parameter, argument in validated_args.items():
		logger.info(f"{parameter}: {argument}")


class ModifiedHelpFormatter(argparse.HelpFormatter):

	# prog is the name of the program 'ex: chewBBACA.py'
	def __init__(self, prog, indent_increment=1, max_help_position=56, width=None):
		if width is None:
			# Try to determine terminal width to optimize width of displayed text
			# Fall back to 80 if not running in an interactive environment such as when using Pipes and CI/CD
			width = shutil.get_terminal_size(fallback=(80, 24)).columns

		super().__init__(prog, indent_increment, max_help_position, width)

	def add_usage(self, usage, actions, groups, prefix=None):
		if prefix is None:
			# Change standard "usage:" prefix to uppercase
			prefix = "USAGE: "

		super().add_usage(usage, actions, groups, prefix)

	# Override split lines method
	def _split_lines(self, text, width):
		lines = super()._split_lines(text, width) + ['']
		return lines

	def _format_action_invocation(self, action):
		if not action.option_strings:
			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			return metavar
		else:
			parts = []
			parts.extend(action.option_strings)
			parts_text = ', '.join(parts)
			return f'{parts_text}'


def run_predict_genes():
	"""Predict genes for genome assemblies."""

	# For development without installing run "python -m CHEWBBACA.chewBBACA PredictGenes -h" from the top directory
	from CHEWBBACA.PredictGenes import constants as ct
	from CHEWBBACA.PredictGenes import validation_models as vm
	from CHEWBBACA.PredictGenes import predict_genes

	parser = argparse.ArgumentParser(prog=gct.PREDICTGENES_MODULE_NAME,
									 description=gct.PREDICTGENES_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --input-path <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.PREDICTGENES_MODULE_NAME)}")

	parser.add_argument("-i", "--input-path", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_INPUT_PATH_ARGNAME],
						required=False, default=vm.PredictGenesValidator.model_fields[ct.PREDICTGENES_INPUT_PATH_ARGNAME].default,
						dest=ct.PREDICTGENES_INPUT_PATH_ARGNAME,
						help=gval.get_field_help(vm.PredictGenesValidator, ct.PREDICTGENES_INPUT_PATH_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.PredictGenesValidator.model_fields[ct.PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.PredictGenesValidator, ct.PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--gp", "--gene-predictor", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_GENE_PREDICTOR_ARGNAME],
					 	required=False, default=vm.PredictGenesValidator.model_fields[ct.PREDICTGENES_GENE_PREDICTOR_ARGNAME].default,
						dest=ct.PREDICTGENES_GENE_PREDICTOR_ARGNAME,
						help=gval.get_field_help(vm.PredictGenesValidator, ct.PREDICTGENES_GENE_PREDICTOR_ARGNAME))

	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME],
					    nargs="+", required=False, default=vm.PredictGenesValidator.model_fields[ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME].default,
						dest=ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.PredictGenesValidator, ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_TRANSLATION_TABLE_ARGNAME],
						required=False, default=vm.PredictGenesValidator.model_fields[ct.PREDICTGENES_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.PREDICTGENES_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.PredictGenesValidator, ct.PREDICTGENES_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_CPU_CORES_ARGNAME],
						required=False, default=vm.PredictGenesValidator.model_fields[ct.PREDICTGENES_CPU_CORES_ARGNAME].default,
						dest=ct.PREDICTGENES_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.PredictGenesValidator, ct.PREDICTGENES_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.PredictGenesValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	# Log command
	log_cli_command(gct.PREDICTGENES_MODULE_NAME)
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the PredictGenes module
	logger.info(f"Running the {gct.PREDICTGENES_MODULE_NAME} module...")
	try:
		predict_genes.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def run_create_schema():
	"""Create a schema based on genome assemblies or coding sequences."""

	from CHEWBBACA.CreateSchema import constants as ct
	from CHEWBBACA.CreateSchema import validation_models as vm
	from CHEWBBACA.CreateSchema import create_schema

	parser = argparse.ArgumentParser(prog=gct.CREATESCHEMA_MODULE_NAME,
									 description=gct.CREATESCHEMA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --input-path <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog="It is strongly advised to provide a training file to create a schema. "
											f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.CREATESCHEMA_MODULE_NAME)}")

	parser.add_argument("-i", "--input-path", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_INPUT_PATH_ARGNAME],
						required=True, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_INPUT_PATH_ARGNAME].default,
						dest=ct.CREATESCHEMA_INPUT_PATH_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_INPUT_PATH_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--sn", "--schema-name", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_SCHEMA_NAME_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_SCHEMA_NAME_ARGNAME].default,
						dest=ct.CREATESCHEMA_SCHEMA_NAME_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_SCHEMA_NAME_ARGNAME))

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME].default,
						dest=ct.CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME))

	parser.add_argument("--ml", "--minimum-length", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_MINIMUM_LENGTH_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_MINIMUM_LENGTH_ARGNAME].default,
						dest=ct.CREATESCHEMA_MINIMUM_LENGTH_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_MINIMUM_LENGTH_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_TRANSLATION_TABLE_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.CREATESCHEMA_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--st", "--size-threshold", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_SIZE_THRESHOLD_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_SIZE_THRESHOLD_ARGNAME].default,
						dest=ct.CREATESCHEMA_SIZE_THRESHOLD_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_SIZE_THRESHOLD_ARGNAME))

	parser.add_argument("--gp", "--gene-predictor", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_GENE_PREDICTOR_ARGNAME],
					 	required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_GENE_PREDICTOR_ARGNAME].default,
						dest=ct.CREATESCHEMA_GENE_PREDICTOR_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_GENE_PREDICTOR_ARGNAME))

	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME],
					    nargs="+", required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME].default,
						dest=ct.CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME))

	parser.add_argument("--co", "--clustering-options", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME],
					 	nargs="+", required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME].default,
						dest=ct.CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME))

	parser.add_argument("--bp", "--blast-path", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_BLAST_PATH_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_BLAST_PATH_ARGNAME].default,
						dest=ct.CREATESCHEMA_BLAST_PATH_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_BLAST_PATH_ARGNAME))

	parser.add_argument("--ci", "--cds-input", action="store_true",
					 	required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_CDS_INPUT_ARGNAME].default,
						dest=ct.CREATESCHEMA_CDS_INPUT_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_CDS_INPUT_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_CPU_CORES_ARGNAME],
						required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_CPU_CORES_ARGNAME].default,
						dest=ct.CREATESCHEMA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_CPU_CORES_ARGNAME))

	parser.add_argument("--no-cds-renaming", action="store_true",
					 	required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_NO_CDS_RENAMING_ARGNAME].default,
						dest=ct.CREATESCHEMA_NO_CDS_RENAMING_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_NO_CDS_RENAMING_ARGNAME))

	parser.add_argument("--no-cleanup", action='store_true',
					 	required=False, default=vm.CreateSchemaValidator.model_fields[ct.CREATESCHEMA_NO_CLEANUP_ARGNAME].default,
						dest=ct.CREATESCHEMA_NO_CLEANUP_ARGNAME,
						help=gval.get_field_help(vm.CreateSchemaValidator, ct.CREATESCHEMA_NO_CLEANUP_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.CreateSchemaValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.CREATESCHEMA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the CreateSchema module
	logger.info(f"Running the {gct.CREATESCHEMA_MODULE_NAME} module...")
	try:
		nloci = create_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
	
	logger.success(f"Process completed successfully!")


def run_allele_call():
	"""Determine allelic profiles."""

	from CHEWBBACA.AlleleCall import constants as ct
	from CHEWBBACA.AlleleCall import validation_models as vm
	from CHEWBBACA.AlleleCall import allele_call

	parser = argparse.ArgumentParser(prog=gct.ALLELECALL_MODULE_NAME,
									 description=gct.ALLELECALL_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --input-path <PATH> --schema-directory <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog="It is strongly advised to perform allele calling with the schema's "
									 		"parameters to ensure the consistency of the results. "
											f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.ALLELECALL_MODULE_NAME)}")

	parser.add_argument("-i", "--input-path", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_INPUT_PATH_ARGNAME],
						required=True, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_INPUT_PATH_ARGNAME].default,
						dest=ct.ALLELECALL_INPUT_PATH_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_INPUT_PATH_ARGNAME))

	parser.add_argument("-s", "--schema-directory", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.ALLELECALL_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.ALLELECALL_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--ll", "--loci-list", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_LOCI_LIST_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_LOCI_LIST_ARGNAME].default,
						dest=ct.ALLELECALL_LOCI_LIST_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_LOCI_LIST_ARGNAME))

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_BLAST_SCORE_RATIO_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_BLAST_SCORE_RATIO_ARGNAME].default,
						dest=ct.ALLELECALL_BLAST_SCORE_RATIO_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_BLAST_SCORE_RATIO_ARGNAME))

	parser.add_argument("--ml", "--minimum-length", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_MINIMUM_LENGTH_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_MINIMUM_LENGTH_ARGNAME].default,
						dest=ct.ALLELECALL_MINIMUM_LENGTH_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_MINIMUM_LENGTH_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_TRANSLATION_TABLE_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.ALLELECALL_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--st", "--size-threshold", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_SIZE_THRESHOLD_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_SIZE_THRESHOLD_ARGNAME].default,
						dest=ct.ALLELECALL_SIZE_THRESHOLD_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_SIZE_THRESHOLD_ARGNAME))

	parser.add_argument("-gp", "--gene-predictor", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_GENE_PREDICTOR_ARGNAME],
					 	required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_GENE_PREDICTOR_ARGNAME].default,
						dest=ct.ALLELECALL_GENE_PREDICTOR_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_GENE_PREDICTOR_ARGNAME))	

	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME],
					 	nargs="+", required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME].default,
						dest=ct.ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME))

	parser.add_argument("--co", "--clustering-options", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_CLUSTERING_OPTIONS_ARGNAME],
					 	nargs="+", required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_CLUSTERING_OPTIONS_ARGNAME].default,
						dest=ct.ALLELECALL_CLUSTERING_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_CLUSTERING_OPTIONS_ARGNAME))

	parser.add_argument("--bp", "--blast-path", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_BLAST_PATH_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_BLAST_PATH_ARGNAME].default,
						dest=ct.ALLELECALL_BLAST_PATH_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_BLAST_PATH_ARGNAME))

	parser.add_argument("--ci", "--cds-input", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_CDS_INPUT_ARGNAME].default,
						dest=ct.ALLELECALL_CDS_INPUT_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_CDS_INPUT_ARGNAME))

	parser.add_argument("--ni", "--no-inferred", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_NO_INFERRED_ARGNAME].default,
						dest=ct.ALLELECALL_NO_INFERRED_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_NO_INFERRED_ARGNAME))

	parser.add_argument("--oo", "--output-unclassified", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_OUTPUT_UNCLASSIFIED_ARGNAME].default,
						dest=ct.ALLELECALL_OUTPUT_UNCLASSIFIED_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_OUTPUT_UNCLASSIFIED_ARGNAME))

	parser.add_argument("--om", "--output-missing", action="store_true",
					 	required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_OUTPUT_MISSING_ARGNAME].default,
						dest=ct.ALLELECALL_OUTPUT_MISSING_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_OUTPUT_MISSING_ARGNAME))

	parser.add_argument("--on", "--output-novel", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_OUTPUT_NOVEL_ARGNAME].default,
						dest=ct.ALLELECALL_OUTPUT_NOVEL_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_OUTPUT_NOVEL_ARGNAME))

	parser.add_argument("--omp", "--output-masked-profiles", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_OUTPUT_MASKED_ARGNAME].default,
						dest=ct.ALLELECALL_OUTPUT_MASKED_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_OUTPUT_MASKED_ARGNAME))

	parser.add_argument("--fc", "--force-continue", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_FORCE_CONTINUE_ARGNAME].default,
						dest=ct.ALLELECALL_FORCE_CONTINUE_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_FORCE_CONTINUE_ARGNAME))

	parser.add_argument("--em", "--execution-mode", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_EXECUTION_MODE_ARGNAME],
					 	required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_EXECUTION_MODE_ARGNAME].default,
						dest=ct.ALLELECALL_EXECUTION_MODE_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_EXECUTION_MODE_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_CPU_CORES_ARGNAME],
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_CPU_CORES_ARGNAME].default,
						dest=ct.ALLELECALL_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_CPU_CORES_ARGNAME))

	parser.add_argument("--no-cds-renaming", action="store_true",
					 	required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_NO_CDS_RENAMING_ARGNAME].default,
						dest=ct.ALLELECALL_NO_CDS_RENAMING_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_NO_CDS_RENAMING_ARGNAME))

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, default=vm.AlleleCallValidator.model_fields[ct.ALLELECALL_NO_CLEANUP_ARGNAME].default,
						dest=ct.ALLELECALL_NO_CLEANUP_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallValidator, ct.ALLELECALL_NO_CLEANUP_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.AlleleCallValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.ALLELECALL_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the AlleleCall module
	logger.info(f"Running the {gct.ALLELECALL_MODULE_NAME} module...")
	try:
		allele_call.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def run_evaluate_schema():
	"""Build an interactive report for schema evaluation."""

	from CHEWBBACA.SchemaEvaluator import constants as ct
	from CHEWBBACA.SchemaEvaluator import validation_models as vm
	from CHEWBBACA.SchemaEvaluator import evaluate_schema

	parser = argparse.ArgumentParser(prog=gct.SCHEMAEVALUATOR_MODULE_NAME,
									 description=gct.SCHEMAEVALUATOR_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --schema-directory <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog="The module can evaluate schemas created with chewBBACA or from external platforms. "
									 		f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.SCHEMAEVALUATOR_MODULE_NAME)}")

	parser.add_argument("-s", "--schema-directory", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME],
					 	required=True, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME],
					 	required=True, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--ll", "--loci-list", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_LOCI_LIST_ARGNAME],
						required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_LOCI_LIST_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_LOCI_LIST_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_LOCI_LIST_ARGNAME))

	parser.add_argument("--af", "--annotations-file", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME],
					 	required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME],
						required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--st", "--size-threshold", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME],
						required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME))

	parser.add_argument("--ml", "--minimum-length", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME],
						required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_CPU_CORES_ARGNAME],
						required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_CPU_CORES_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_CPU_CORES_ARGNAME))

	parser.add_argument("--lr", "--loci-reports", action="store_true",
					 	required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_LOCI_REPORTS_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_LOCI_REPORTS_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_LOCI_REPORTS_ARGNAME))

	parser.add_argument("--light", action="store_true",
					 	required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_LIGTH_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_LIGTH_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_LIGTH_ARGNAME))

	parser.add_argument("--add-sequences", action="store_true",
					 	required=False, default=vm.SchemaEvaluatorValidator.model_fields[ct.SCHEMAEVALUATOR_ADD_SEQUENCES_ARGNAME].default,
						dest=ct.SCHEMAEVALUATOR_ADD_SEQUENCES_ARGNAME,
						help=gval.get_field_help(vm.SchemaEvaluatorValidator, ct.SCHEMAEVALUATOR_ADD_SEQUENCES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.SchemaEvaluatorValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)
	
	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.SCHEMAEVALUATOR_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the SchemaEvaluator module
	logger.info(f"Running the {gct.SCHEMAEVALUATOR_MODULE_NAME} module...")
	try:
		evaluate_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def run_evaluate_calls():
	"""Build an interactive report for allele calling evaluation."""

	from CHEWBBACA.AlleleCallEvaluator import constants as ct
	from CHEWBBACA.AlleleCallEvaluator import validation_models as vm
	from CHEWBBACA.AlleleCallEvaluator import evaluate_calls

	parser = argparse.ArgumentParser(prog=gct.ALLELECALLEVALUATOR_MODULE_NAME,
									 description=gct.ALLELECALLEVALUATOR_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --results-files <PATH> --schema-directory <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.ALLELECALLEVALUATOR_MODULE_NAME)}")

	parser.add_argument("-r", "--results-files", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME],
					 	required=True, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME))

	parser.add_argument("-s", "--schema-directory", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME],
					 	required=True, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME],
					 	required=True, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--af", "--annotations-file", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME],
					 	required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_CPU_CORES_ARGNAME],
						required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_CPU_CORES_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_CPU_CORES_ARGNAME))

	parser.add_argument("--npa", "--no-presence-absence", action="store_true",
					 	required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_NO_PRESENCE_ABSENCE_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_NO_PRESENCE_ABSENCE_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_NO_PRESENCE_ABSENCE_ARGNAME))

	parser.add_argument("--ndm", "--no-distance-matrix", action='store_true',
					 	required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_NO_DISTANCE_MATRIX_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_NO_DISTANCE_MATRIX_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_NO_DISTANCE_MATRIX_ARGNAME))

	parser.add_argument("--nnj", "--no-neighbor-joining", action="store_true",
					 	required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_NO_NEIGHBOR_JOINING_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_NO_NEIGHBOR_JOINING_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_NO_NEIGHBOR_JOINING_ARGNAME))

	parser.add_argument("--fcm", "--force-core-msa", action="store_true",
					 	required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_FORCE_CORE_MSA_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_FORCE_CORE_MSA_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_FORCE_CORE_MSA_ARGNAME))

	parser.add_argument("--light", action="store_true",
					 	required=False, default=vm.AlleleCallEvaluatorValidator.model_fields[ct.ALLELECALLEVALUATOR_LIGTH_ARGNAME].default,
						dest=ct.ALLELECALLEVALUATOR_LIGTH_ARGNAME,
						help=gval.get_field_help(vm.AlleleCallEvaluatorValidator, ct.ALLELECALLEVALUATOR_LIGTH_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.AlleleCallEvaluatorValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.ALLELECALLEVALUATOR_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the AlleleCallEvaluator module
	logger.info(f"Running the {gct.ALLELECALLEVALUATOR_MODULE_NAME} module...")
	try:
		evaluate_calls.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_compute_core():
	"""Determine the set of core loci based on allele calling results."""

	from CHEWBBACA.ComputeCore import constants as ct
	from CHEWBBACA.ComputeCore import validation_models as vm
	from CHEWBBACA.ComputeCore import compute_core

	parser = argparse.ArgumentParser(prog=gct.COMPUTECORE_MODULE_NAME,
									 description=gct.COMPUTECORE_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --profiles <FILE> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.COMPUTECORE_MODULE_NAME)}")

	parser.add_argument("-p", "--profiles", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_PROFILES_ARGNAME],
						required=True, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_PROFILES_ARGNAME].default,
						dest=ct.COMPUTECORE_PROFILES_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_PROFILES_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--lpt", "--loci-presence-threshold", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME], nargs="+",
						required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME].default,
						dest=ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME))

	parser.add_argument("--ss", "--sample-step", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_SAMPLE_STEP_ARGNAME],
					 	required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_SAMPLE_STEP_ARGNAME].default,
						dest=ct.COMPUTECORE_SAMPLE_STEP_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_SAMPLE_STEP_ARGNAME))

	parser.add_argument("--ca", "--compute-accessory", action="store_true",
					 	required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_COMPUTE_ACCESSORY_ARGNAME].default,
						dest=ct.COMPUTECORE_COMPUTE_ACCESSORY_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_COMPUTE_ACCESSORY_ARGNAME))

	parser.add_argument("--ra", "--rarefaction-analysis", action="store_true",
					 	required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_RAREFACTION_ANALYSIS_ARGNAME].default,
						dest=ct.COMPUTECORE_RAREFACTION_ANALYSIS_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_RAREFACTION_ANALYSIS_ARGNAME))

	parser.add_argument("--pn", "--permutation-number", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_PERMUTATION_NUMBER_ARGNAME],
					 	required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_PERMUTATION_NUMBER_ARGNAME].default,
						dest=ct.COMPUTECORE_PERMUTATION_NUMBER_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_PERMUTATION_NUMBER_ARGNAME))

	parser.add_argument("--ps", "--permutation-samples", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME],
					 	required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME].default,
						dest=ct.COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_CPU_CORES_ARGNAME],
						required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_CPU_CORES_ARGNAME].default,
						dest=ct.COMPUTECORE_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_CPU_CORES_ARGNAME))

	parser.add_argument("--el", "--exclude-loci", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_EXCLUDE_LOCI_ARGNAME],
						required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_EXCLUDE_LOCI_ARGNAME].default,
						dest=ct.COMPUTECORE_EXCLUDE_LOCI_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_EXCLUDE_LOCI_ARGNAME))

	parser.add_argument("--eg", "--exclude-genomes", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_EXCLUDE_GENOMES_ARGNAME],
						required=False, default=vm.ComputeCoreValidator.model_fields[ct.COMPUTECORE_EXCLUDE_GENOMES_ARGNAME].default,
						dest=ct.COMPUTECORE_EXCLUDE_GENOMES_ARGNAME,
						help=gval.get_field_help(vm.ComputeCoreValidator, ct.COMPUTECORE_EXCLUDE_GENOMES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.ComputeCoreValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.COMPUTECORE_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the ComputeCore module
	logger.info(f"Running the {gct.COMPUTECORE_MODULE_NAME} module...")
	try:
		compute_core.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_subset_results():
	"""Subset results files."""

	from CHEWBBACA.SubsetResults import constants as ct
	from CHEWBBACA.SubsetResults import validation_models as vm
	from CHEWBBACA.SubsetResults import subset_results

	parser = argparse.ArgumentParser(prog=gct.SUBSETRESULTS_MODULE_NAME,
									 description=gct.SUBSETRESULTS_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --results-files <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.SUBSETRESULTS_MODULE_NAME)}")

	parser.add_argument("-r", "--results-files", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_RESULTS_FILES_ARGNAME],
						required=True, default=vm.SubsetResultsValidator.model_fields[ct.SUBSETRESULTS_RESULTS_FILES_ARGNAME].default,
						dest=ct.SUBSETRESULTS_RESULTS_FILES_ARGNAME,
						help=gval.get_field_help(vm.SubsetResultsValidator, ct.SUBSETRESULTS_RESULTS_FILES_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.SubsetResultsValidator.model_fields[ct.SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.SubsetResultsValidator, ct.SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--ll", "--loci-list", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_LOCI_LIST_ARGNAME],
						required=False, default=vm.SubsetResultsValidator.model_fields[ct.SUBSETRESULTS_LOCI_LIST_ARGNAME].default,
						dest=ct.SUBSETRESULTS_LOCI_LIST_ARGNAME,
						help=gval.get_field_help(vm.SubsetResultsValidator, ct.SUBSETRESULTS_LOCI_LIST_ARGNAME))

	parser.add_argument("--sl", "--sample-list", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_SAMPLE_LIST_ARGNAME],
						required=False, default=vm.SubsetResultsValidator.model_fields[ct.SUBSETRESULTS_SAMPLE_LIST_ARGNAME].default,
						dest=ct.SUBSETRESULTS_SAMPLE_LIST_ARGNAME,
						help=gval.get_field_help(vm.SubsetResultsValidator, ct.SUBSETRESULTS_SAMPLE_LIST_ARGNAME))

	parser.add_argument("--invert-loci", action="store_true",
						required=False, default=vm.SubsetResultsValidator.model_fields[ct.SUBSETRESULTS_INVERT_LOCI_ARGNAME].default,
						dest=ct.SUBSETRESULTS_INVERT_LOCI_ARGNAME,
						help=gval.get_field_help(vm.SubsetResultsValidator, ct.SUBSETRESULTS_INVERT_LOCI_ARGNAME))

	parser.add_argument("--invert-samples", action="store_true",
						required=False, default=vm.SubsetResultsValidator.model_fields[ct.SUBSETRESULTS_INVERT_SAMPLES_ARGNAME].default,
						dest=ct.SUBSETRESULTS_INVERT_SAMPLES_ARGNAME,
						help=gval.get_field_help(vm.SubsetResultsValidator, ct.SUBSETRESULTS_INVERT_SAMPLES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.SubsetResultsValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.SUBSETRESULTS_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the SubsetResults module
	logger.info(f"Running the {gct.SUBSETRESULTS_MODULE_NAME} module...")
	try:
		subset_results.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_merge_results():
	"""Merge results files."""

	from CHEWBBACA.MergeResults import constants as ct
	from CHEWBBACA.MergeResults import validation_models as vm
	from CHEWBBACA.MergeResults import merge_results

	parser = argparse.ArgumentParser(prog=gct.MERGERESULTS_MODULE_NAME,
									 description=gct.MERGERESULTS_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --results-directories <PATH> --output-directory <PATH> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.MERGERESULTS_MODULE_NAME)}")

	parser.add_argument("-r", "--results-directories", nargs="+", type=ct.MERGERESULTS_ARGUMENT_TYPES[ct.MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME],
						required=True, default=vm.MergeResultsValidator.model_fields[ct.MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME].default,
						dest=ct.MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME,
						help=gval.get_field_help(vm.MergeResultsValidator, ct.MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.MERGERESULTS_ARGUMENT_TYPES[ct.MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.MergeResultsValidator.model_fields[ct.MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.MergeResultsValidator, ct.MERGERESULTS_COMMON_ARGNAME))

	parser.add_argument("--common", action="store_true",
						required=False, default=vm.MergeResultsValidator.model_fields[ct.MERGERESULTS_COMMON_ARGNAME].default,
						dest=ct.MERGERESULTS_COMMON_ARGNAME,
						help=gval.get_field_help(vm.MergeResultsValidator, ct.MERGERESULTS_COMMON_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.MergeResultsValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.MERGERESULTS_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the MergeResults module
	logger.info(f"Running the {gct.MERGERESULTS_MODULE_NAME} module...")
	try:
		merge_results.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_hash_profiles():
	"""Hash allelic profiles."""

	from CHEWBBACA.HashProfiles import constants as ct
	from CHEWBBACA.HashProfiles import validation_models as vm
	from CHEWBBACA.HashProfiles import hash_profiles

	parser = argparse.ArgumentParser(prog=gct.HASHPROFILES_MODULE_NAME,
									 description=gct.HASHPROFILES_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --allelic-profiles <FILE> --schema-directory <DIR> --output-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.HASHPROFILES_MODULE_NAME)}")

	parser.add_argument("-a", "--allelic-profiles", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_ALLELIC_PROFILES_ARGNAME],
						required=True, default=vm.HashProfilesValidator.model_fields[ct.HASHPROFILES_ALLELIC_PROFILES_ARGNAME].default,
						dest=ct.HASHPROFILES_ALLELIC_PROFILES_ARGNAME,
						help=gval.get_field_help(vm.HashProfilesValidator, ct.HASHPROFILES_ALLELIC_PROFILES_ARGNAME))

	parser.add_argument("-s", "--schema-directory", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.HashProfilesValidator.model_fields[ct.HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.HashProfilesValidator, ct.HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.HashProfilesValidator.model_fields[ct.HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.HashProfilesValidator, ct.HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--ha", "--hashing-algo", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_HASHING_ALGO_ARGNAME],
					 	required=False, default=vm.HashProfilesValidator.model_fields[ct.HASHPROFILES_HASHING_ALGO_ARGNAME].default,
						dest=ct.HASHPROFILES_HASHING_ALGO_ARGNAME,
						help=gval.get_field_help(vm.HashProfilesValidator, ct.HASHPROFILES_HASHING_ALGO_ARGNAME))

	parser.add_argument("--nc", "--nrow-chunk", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_NROW_CHUNK_ARGNAME],
					 	required=False, default=vm.HashProfilesValidator.model_fields[ct.HASHPROFILES_NROW_CHUNK_ARGNAME].default,
						dest=ct.HASHPROFILES_NROW_CHUNK_ARGNAME,
						help=gval.get_field_help(vm.HashProfilesValidator, ct.HASHPROFILES_NROW_CHUNK_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_CPU_CORES_ARGNAME],
						required=False, default=vm.HashProfilesValidator.model_fields[ct.HASHPROFILES_CPU_CORES_ARGNAME].default,
						dest=ct.HASHPROFILES_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.HashProfilesValidator, ct.HASHPROFILES_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.HashProfilesValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.HASHPROFILES_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the HashProfiles module
	logger.info(f"Running the {gct.HASHPROFILES_MODULE_NAME} module...")
	try:
		hash_profiles.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_get_alleles():
	"""Get allele sequences."""

	from CHEWBBACA.GetAlleles import constants as ct
	from CHEWBBACA.GetAlleles import validation_models as vm
	from CHEWBBACA.GetAlleles import get_alleles

	parser = argparse.ArgumentParser(prog=gct.GETALLELES_MODULE_NAME,
									 description=gct.GETALLELES_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --allelic-profiles <FILE> --schema-directory <DIR> --output-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.GETALLELES_MODULE_NAME)}")

	parser.add_argument("-a", "--allelic-profiles", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_ALLELIC_PROFILES_ARGNAME],
						required=True, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_ALLELIC_PROFILES_ARGNAME].default,
						dest=ct.GETALLELES_ALLELIC_PROFILES_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_ALLELIC_PROFILES_ARGNAME))

	parser.add_argument("-s", "--schema-directory", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.GETALLELES_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.GETALLELES_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_OUTPUT_DIRECTORY_ARGNAME))
	
	parser.add_argument("-ll", "--loci-list", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_LOCI_LIST_ARGNAME],
						required=False, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_LOCI_LIST_ARGNAME].default,
						dest=ct.GETALLELES_LOCI_LIST_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_LOCI_LIST_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_CPU_CORES_ARGNAME],
						required=False, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_CPU_CORES_ARGNAME].default,
						dest=ct.GETALLELES_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_CPU_CORES_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_TRANSLATION_TABLE_ARGNAME],
					 	required=False, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.GETALLELES_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--ta", "--translate-alleles", action="store_true",
						required=False, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_TRANSLATE_ALLELES_ARGNAME].default,
						dest=ct.GETALLELES_TRANSLATE_ALLELES_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_TRANSLATE_ALLELES_ARGNAME))

	parser.add_argument("--distinct", action="store_true",
						required=False, default=vm.GetAllelesValidator.model_fields[ct.GETALLELES_DISTINCT_ARGNAME].default,
						dest=ct.GETALLELES_DISTINCT_ARGNAME,
						help=gval.get_field_help(vm.GetAllelesValidator, ct.GETALLELES_DISTINCT_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.GetAllelesValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.GETALLELES_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the GetAlleles module
	logger.info(f"Running the {gct.GETALLELES_MODULE_NAME} module...")
	try:
		get_alleles.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_adapt_schema():
	"""Adapt schemas for use with chewBBACA."""

	from CHEWBBACA.AdaptSchema import constants as ct
	from CHEWBBACA.AdaptSchema import validation_models as vm
	from CHEWBBACA.AdaptSchema import adapt_schema

	parser = argparse.ArgumentParser(prog=gct.ADAPTSCHEMA_MODULE_NAME,
									 description=gct.ADAPTSCHEMA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --schema-directory <DIR> --output-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.ADAPTSCHEMA_MODULE_NAME)}")

	parser.add_argument("-s", "--schema-directory", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--ll", "--loci-list", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_LOCI_LIST_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_LOCI_LIST_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_LOCI_LIST_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_LOCI_LIST_ARGNAME))

	parser.add_argument("-gp", "--gene-predictor", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME))
	
	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME],
						nargs="+", required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME))

	parser.add_argument("--co", "--clustering-options", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME],
						nargs="+", required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME))

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME))

	parser.add_argument("--ml", "--minimum-length", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME],
					 	required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--st", "--size-threshold", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_CPU_CORES_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_CPU_CORES_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_CPU_CORES_ARGNAME))

	parser.add_argument("--bp", "--blast-path", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_BLAST_PATH_ARGNAME],
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_BLAST_PATH_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_BLAST_PATH_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_BLAST_PATH_ARGNAME))

	parser.add_argument("--size-filter", action="store_true",
						required=False, default=vm.AdaptSchemaValidator.model_fields[ct.ADAPTSCHEMA_SIZE_FILTER_ARGNAME].default,
						dest=ct.ADAPTSCHEMA_SIZE_FILTER_ARGNAME,
						help=gval.get_field_help(vm.AdaptSchemaValidator, ct.ADAPTSCHEMA_SIZE_FILTER_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.AdaptSchemaValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

########### Is this writing the log file into the schema folder???

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.ADAPTSCHEMA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the AdaptSchema module
	logger.info(f"Running the {gct.ADAPTSCHEMA_MODULE_NAME} module...")
	try:
		adapt_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")



def run_annotate_schema():
	"""Retrieve annotations for loci in a schema."""

	from CHEWBBACA.AnnotateSchema import constants as ct
	from CHEWBBACA.AnnotateSchema import validation_models as vm
	from CHEWBBACA.AnnotateSchema import annotate_schema

	parser = argparse.ArgumentParser(prog=gct.ANNOTATESCHEMA_MODULE_NAME,
									 description=gct.ANNOTATESCHEMA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --schema-directory <DIR> --output-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.ANNOTATESCHEMA_MODULE_NAME)}")

	parser.add_argument("-s", "--schema-directory", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--ll", "--loci-list", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_LOCI_LIST_ARGNAME],
						required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_LOCI_LIST_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_LOCI_LIST_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_LOCI_LIST_ARGNAME))

	parser.add_argument("--pt", "--protein-table", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME],
						required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME))

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME],
					 	required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME],
					 	required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_CPU_CORES_ARGNAME],
						required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_CPU_CORES_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_CPU_CORES_ARGNAME))

	parser.add_argument("--ptax", "--proteome-taxa", nargs="+", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME],
						required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME))

	parser.add_argument("--pm", "--proteome-matches", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME],
					 	required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME))

	parser.add_argument("--bp", "--blast-path", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_BLAST_PATH_ARGNAME],
						required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_BLAST_PATH_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_BLAST_PATH_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_BLAST_PATH_ARGNAME))

	parser.add_argument("--rapi", "--rest-api", action="store_true",
	                    required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_UNIPROT_REST_API_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_UNIPROT_REST_API_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_UNIPROT_REST_API_ARGNAME))

	parser.add_argument("--rm", "--rest-api-matches", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME],
					 	required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME))

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, default=vm.AnnotateSchemaValidator.model_fields[ct.ANNOTATESCHEMA_NO_CLEANUP_ARGNAME].default,
						dest=ct.ANNOTATESCHEMA_NO_CLEANUP_ARGNAME,
						help=gval.get_field_help(vm.AnnotateSchemaValidator, ct.ANNOTATESCHEMA_NO_CLEANUP_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.AnnotateSchemaValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.ANNOTATESCHEMA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the AnnotateSchema module
	logger.info(f"Running the {gct.ANNOTATESCHEMA_MODULE_NAME} module...")
	try:
		annotate_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
		
	logger.success(f"Process completed successfully!")


def run_compute_distances():
	"""Compute pairwise distances based on allele calling results."""

	from CHEWBBACA.ComputeDistances import constants as ct
	from CHEWBBACA.ComputeDistances import validation_models as vm
	from CHEWBBACA.ComputeDistances import compute_distances

	parser = argparse.ArgumentParser(prog=gct.COMPUTEDISTANCES_MODULE_NAME,
									 description=gct.COMPUTECORE_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --allelic-profiles <FILE> --output-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.COMPUTEDISTANCES_MODULE_NAME)}")

	parser.add_argument("-a", "--allelic-profiles", type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME],
						required=True, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--dcm", "--distance-computation-method", type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME],
					 	required=False, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME))

	parser.add_argument('--of', '--output-format', type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME],
					 	required=False, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME))

	parser.add_argument("--cs", "--compute-similarity", action="store_true",
					 	required=False, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_COMPUTE_SIMILARITY_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_COMPUTE_SIMILARITY_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_COMPUTE_SIMILARITY_ARGNAME))

	parser.add_argument("--no-mask", action="store_true",
					 	required=False, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_NO_MASK_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_NO_MASK_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_NO_MASK_ARGNAME))

	parser.add_argument('--cc', '--cpu-cores', type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_CPU_CORES_ARGNAME],
						required=False, default=vm.ComputeDistancesValidator.model_fields[ct.COMPUTEDISTANCES_CPU_CORES_ARGNAME].default,
						dest=ct.COMPUTEDISTANCES_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.ComputeDistancesValidator, ct.COMPUTEDISTANCES_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.ComputeDistancesValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.COMPUTEDISTANCES_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the ComputeDistances module
	logger.info(f"Running the {gct.COMPUTEDISTANCES_MODULE_NAME} module...")
	try:
		compute_distances.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def run_compute_msa():
	"""Compute Multiple Sequence Alignments."""

	from CHEWBBACA.ComputeMSA import constants as ct
	from CHEWBBACA.ComputeMSA import validation_models as vm
	from CHEWBBACA.ComputeMSA import compute_msa

	parser = argparse.ArgumentParser(prog=gct.COMPUTEMSA_MODULE_NAME,
									 description=gct.COMPUTEMSA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --input-path <PATH> --output-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.COMPUTEMSA_MODULE_NAME)}")

	parser.add_argument("-i", "--input-path", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_INPUT_PATH_ARGNAME],
						required=True, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_INPUT_PATH_ARGNAME].default,
						dest=ct.COMPUTEMSA_INPUT_PATH_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_INPUT_PATH_ARGNAME))

	parser.add_argument("-o", "--output-directory", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("-s", "--schema-directory", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME],
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("--tt", "--translation-table", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_TRANSLATION_TABLE_ARGNAME],
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_TRANSLATION_TABLE_ARGNAME].default,
						dest=ct.COMPUTEMSA_TRANSLATION_TABLE_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_TRANSLATION_TABLE_ARGNAME))

	parser.add_argument("--dm", "--dna-msa", action="store_true",
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_DNA_MSA_ARGNAME].default,
						dest=ct.COMPUTEMSA_DNA_MSA_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_DNA_MSA_ARGNAME))

	parser.add_argument("--ov", "--output-variable", action="store_true",
					 	required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_OUTPUT_VARIABLE_ARGNAME].default,
						dest=ct.COMPUTEMSA_OUTPUT_VARIABLE_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_OUTPUT_VARIABLE_ARGNAME))

	parser.add_argument("--gaps", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_GAPS_ARGNAME],
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_GAPS_ARGNAME].default,
						dest=ct.COMPUTEMSA_GAPS_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_GAPS_ARGNAME))

	parser.add_argument("--ambiguous", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_AMBIGUOUS_ARGNAME],
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_AMBIGUOUS_ARGNAME].default,
						dest=ct.COMPUTEMSA_AMBIGUOUS_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_AMBIGUOUS_ARGNAME))

	parser.add_argument("--only-loci-msas", action="store_true",
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_ONLY_LOCI_MSAS_ARGNAME].default,
						dest=ct.COMPUTEMSA_ONLY_LOCI_MSAS_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_ONLY_LOCI_MSAS_ARGNAME))

	parser.add_argument("--cmo", "--custom-mafft-options", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME],
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME].default,
						dest=ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME))

	parser.add_argument("--pi", "--protein-input", action="store_true",
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_PROTEIN_INPUT_ARGNAME].default,
						dest=ct.COMPUTEMSA_PROTEIN_INPUT_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_PROTEIN_INPUT_ARGNAME))

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_NO_CLEANUP_ARGNAME].default,
						dest=ct.COMPUTEMSA_NO_CLEANUP_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_NO_CLEANUP_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_CPU_CORES_ARGNAME],
						required=False, default=vm.ComputeMSAValidator.model_fields[ct.COMPUTEMSA_CPU_CORES_ARGNAME].default,
						dest=ct.COMPUTEMSA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.ComputeMSAValidator, ct.COMPUTEMSA_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.ComputeMSAValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.COMPUTEMSA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())
	
	# Run the ComputeMSA module
	logger.info(f"Running the {gct.COMPUTEMSA_MODULE_NAME} module...")
	try:
		compute_msa.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def run_download_schema():
	"""Download schemas from Chewie-NS."""

	from CHEWBBACA.CHEWBBACA_NS import constants as ct
	from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
	from CHEWBBACA.CHEWBBACA_NS import download_schema

	parser = argparse.ArgumentParser(prog=gct.DOWNLOADSCHEMA_MODULE_NAME,
									 description=gct.DOWNLOADSCHEMA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --output-directory <DIR> --species-id <INT/STR> --schema-id <INT/STR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.DOWNLOADSCHEMA_MODULE_NAME)}")

	parser.add_argument("-o", "--output-directory", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME))

	parser.add_argument("--spi", "--species-id", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_SPECIES_ID_ARGNAME],
						required=True, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_SPECIES_ID_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_SPECIES_ID_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_SPECIES_ID_ARGNAME))

	parser.add_argument("--sci", "--schema-id", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME],
						required=True, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME))

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME))

	parser.add_argument("--sd", "--schema-date", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME],
						required=False, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME))

	parser.add_argument("--lv", "--latest-version", action="store_true",
						required=False, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_LATEST_VERSION_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_LATEST_VERSION_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_LATEST_VERSION_ARGNAME))

	parser.add_argument("--bp", "--blast-path", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_BLAST_PATH_ARGNAME],
						required=False, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_BLAST_PATH_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_BLAST_PATH_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_BLAST_PATH_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_CPU_CORES_ARGNAME],
						required=False, default=vm.DownloadSchemaValidator.model_fields[ct.DOWNLOADSCHEMA_CPU_CORES_ARGNAME].default,
						dest=ct.DOWNLOADSCHEMA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.DownloadSchemaValidator, ct.DOWNLOADSCHEMA_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.DownloadSchemaValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.DOWNLOADSCHEMA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the DownloadSchema module
	logger.info(f"Running the {gct.DOWNLOADSCHEMA_MODULE_NAME} module...")
	try:
		download_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
	
	logger.success(f"Process completed successfully!")


def run_upload_schema():
	"""Upload schemas to Chewie-NS."""

	from CHEWBBACA.CHEWBBACA_NS import constants as ct
	from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
	from CHEWBBACA.CHEWBBACA_NS import upload_schema

	parser = argparse.ArgumentParser(prog=gct.UPLOADSCHEMA_MODULE_NAME,
									 description=gct.UPLOADSCHEMA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --schema-directory <DIR> --species-id <INT/STR> --schema-name <STR> --loci-prefix <STR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.UPLOADSCHEMA_MODULE_NAME)}")

	parser.add_argument("-s", "--schema-directory", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("--spi", "--species-id", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_SPECIES_ID_ARGNAME],
						required=True, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_SPECIES_ID_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_SPECIES_ID_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_SPECIES_ID_ARGNAME))

	parser.add_argument("--sn", "--schema-name", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_SCHEMA_NAME_ARGNAME],
						required=True, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_SCHEMA_NAME_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_SCHEMA_NAME_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_SCHEMA_NAME_ARGNAME))

	parser.add_argument("--lp", "--loci-prefix", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_LOCI_PREFIX_ARGNAME],
						required=True, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_LOCI_PREFIX_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_LOCI_PREFIX_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_LOCI_PREFIX_ARGNAME))

	parser.add_argument("--df", "--description-file", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME],
						required=False, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME))

	parser.add_argument("--af", "--annotations-file", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME],
						required=False, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME))

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME))

	parser.add_argument("--continue-upload", action="store_true",
						required=False, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_CONTINUE_UPLOAD_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_CONTINUE_UPLOAD_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_CONTINUE_UPLOAD_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_CPU_CORES_ARGNAME],
						required=False, default=vm.UploadSchemaValidator.model_fields[ct.UPLOADSCHEMA_CPU_CORES_ARGNAME].default,
						dest=ct.UPLOADSCHEMA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.UploadSchemaValidator, ct.UPLOADSCHEMA_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.UploadSchemaValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.UPLOADSCHEMA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the UploadSchema module
	logger.info(f"Running the {gct.UPLOADSCHEMA_MODULE_NAME} module...")
	try:
		upload_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)
	
	logger.success(f"Process completed successfully!")


def run_synchronize_schema():
	"""Synchronize local schemas with the remote versions on Chewie-NS."""

	from CHEWBBACA.CHEWBBACA_NS import constants as ct
	from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
	from CHEWBBACA.CHEWBBACA_NS import synchronize_schema

	parser = argparse.ArgumentParser(prog=gct.SYNCSCHEMA_MODULE_NAME,
									 description=gct.SYNCSCHEMA_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --schema-directory <DIR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.SYNCSCHEMA_MODULE_NAME)}")

	parser.add_argument("-s", "--schema-directory", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, default=vm.SynchronizeSchemaValidator.model_fields[ct.SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME].default,
						dest=ct.SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help=gval.get_field_help(vm.SynchronizeSchemaValidator, ct.SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME))

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, default=vm.SynchronizeSchemaValidator.model_fields[ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME].default,
						dest=ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help=gval.get_field_help(vm.SynchronizeSchemaValidator, ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME))

	parser.add_argument("--sa", "--submit-alleles", action="store_true",
						required=False, default=vm.SynchronizeSchemaValidator.model_fields[ct.SYNCSCHEMA_SUBMIT_ALLELES_ARGNAME].default,
						dest=ct.SYNCSCHEMA_SUBMIT_ALLELES_ARGNAME,
						help=gval.get_field_help(vm.SynchronizeSchemaValidator, ct.SYNCSCHEMA_SUBMIT_ALLELES_ARGNAME))

	parser.add_argument("--bp", "--blast-path", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_BLAST_PATH_ARGNAME],
						required=False, default=vm.SynchronizeSchemaValidator.model_fields[ct.SYNCSCHEMA_BLAST_PATH_ARGNAME].default,
						dest=ct.SYNCSCHEMA_BLAST_PATH_ARGNAME,
						help=gval.get_field_help(vm.SynchronizeSchemaValidator, ct.SYNCSCHEMA_BLAST_PATH_ARGNAME))

	parser.add_argument("--cc", "--cpu-cores", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_CPU_CORES_ARGNAME],
						required=False, default=vm.SynchronizeSchemaValidator.model_fields[ct.SYNCSCHEMA_CPU_CORES_ARGNAME].default,
						dest=ct.SYNCSCHEMA_CPU_CORES_ARGNAME,
						help=gval.get_field_help(vm.SynchronizeSchemaValidator, ct.SYNCSCHEMA_CPU_CORES_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.SynchronizeSchemaValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.SYNCSCHEMA_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the SynSchema module
	logger.info(f"Running the {gct.SYNCSCHEMA_MODULE_NAME} module...")
	try:
		synchronize_schema.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def run_stats_requests():
	"""Retrieve basic species and schema information from Chewie-NS."""

	from CHEWBBACA.CHEWBBACA_NS import constants as ct
	from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
	from CHEWBBACA.CHEWBBACA_NS import stats_requests

	parser = argparse.ArgumentParser(prog=gct.NSSTATS_MODULE_NAME,
									 description=gct.NSSTATS_MODULE_SHORT_DESCRIPTION,
									 usage="chewBBACA.py %(prog)s --stats-mode <STR> [OPTIONS]",
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(gct.NSSTATS_MODULE_NAME)}")

	parser.add_argument("--sm", "--stats-mode", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_STATS_MODE_ARGNAME],
						required=True, default=vm.NSStatsValidator.model_fields[ct.NSSTATS_STATS_MODE_ARGNAME].default,
						dest=ct.NSSTATS_STATS_MODE_ARGNAME,
						help=gval.get_field_help(vm.NSStatsValidator, ct.NSSTATS_STATS_MODE_ARGNAME))

	parser.add_argument("--spi", "--species-id", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_SPECIES_ID_ARGNAME],
						required=False, default=vm.NSStatsValidator.model_fields[ct.NSSTATS_SPECIES_ID_ARGNAME].default,
						dest=ct.NSSTATS_SPECIES_ID_ARGNAME,
						help=gval.get_field_help(vm.NSStatsValidator, ct.NSSTATS_SPECIES_ID_ARGNAME))

	parser.add_argument("--sci", "--schema-id", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_SCHEMA_ID_ARGNAME],
						required=False, default=vm.NSStatsValidator.model_fields[ct.NSSTATS_SCHEMA_ID_ARGNAME].default,
						dest=ct.NSSTATS_SCHEMA_ID_ARGNAME,
						help=gval.get_field_help(vm.NSStatsValidator, ct.NSSTATS_SCHEMA_ID_ARGNAME))

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, default=vm.NSStatsValidator.model_fields[ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME].default,
						dest=ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help=gval.get_field_help(vm.NSStatsValidator, ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME))

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	try:
		validated_args = vm.NSStatsValidator.model_validate(clean_args)
	except ValidationError as err:
		logger.error(f"Error during argument validation: {err}")
		sys.exit(gval.ExitCode.VALIDATION_ERROR)

	# Initialize file logging
	logfile_path = setup_file_logging(validated_args.output_directory)
	logger.info(f"File logging initialized at: {logfile_path}")
	# Log command
	# Use shlex to retain quotes and other special characters
	raw_cmd = " ".join(shlex.quote(arg) for arg in [sys.argv[0]]+[gct.NSSTATS_MODULE_NAME]+sys.argv[1:])
	logger.info(f"Command: {raw_cmd}")
	# Log list of validated arguments
	log_validated_args(validated_args.model_dump())

	# Run the NSStats module
	logger.info(f"Running the {gct.NSSTATS_MODULE_NAME} module...")
	try:
		stats_requests.main(validated_args)
	except Exception as err:
		logger.error(f"Process error: {err}")
		sys.exit(gval.ExitCode.GENERAL_ERROR)

	logger.success(f"Process completed successfully!")


def main():

	# Need the modules names in the global constants to be able to
	functions_info = {gct.PREDICTGENES_MODULE_NAME: [gct.PREDICTGENES_MODULE_SHORT_DESCRIPTION, run_predict_genes],
					  gct.CREATESCHEMA_MODULE_NAME: [gct.CREATESCHEMA_MODULE_SHORT_DESCRIPTION, run_create_schema],
					  gct.ALLELECALL_MODULE_NAME: [gct.ALLELECALL_MODULE_SHORT_DESCRIPTION, run_allele_call],
					  gct.SCHEMAEVALUATOR_MODULE_NAME: [gct.SCHEMAEVALUATOR_MODULE_SHORT_DESCRIPTION, run_evaluate_schema],
					  gct.ALLELECALLEVALUATOR_MODULE_NAME: [gct.ALLELECALLEVALUATOR_MODULE_SHORT_DESCRIPTION, run_evaluate_calls],
					  gct.COMPUTECORE_MODULE_NAME: [gct.COMPUTECORE_MODULE_SHORT_DESCRIPTION, run_compute_core],
					  gct.SUBSETRESULTS_MODULE_NAME: [gct.SUBSETRESULTS_MODULE_SHORT_DESCRIPTION, run_subset_results],
					  gct.ADAPTSCHEMA_MODULE_NAME: [gct.ADAPTSCHEMA_MODULE_SHORT_DESCRIPTION, run_adapt_schema],
					  gct.MERGERESULTS_MODULE_NAME: [gct.MERGERESULTS_MODULE_SHORT_DESCRIPTION, run_merge_results],
					  gct.HASHPROFILES_MODULE_NAME: [gct.HASHPROFILES_MODULE_SHORT_DESCRIPTION, run_hash_profiles],
					  gct.GETALLELES_MODULE_NAME: [gct.GETALLELES_MODULE_SHORT_DESCRIPTION, run_get_alleles],
					  gct.ANNOTATESCHEMA_MODULE_NAME: [gct.ANNOTATESCHEMA_MODULE_SHORT_DESCRIPTION, run_annotate_schema],
					  gct.COMPUTEDISTANCES_MODULE_NAME: [gct.COMPUTEDISTANCES_MODULE_SHORT_DESCRIPTION, run_compute_distances],
					  gct.COMPUTEMSA_MODULE_NAME: [gct.COMPUTEMSA_MODULE_SHORT_DESCRIPTION, run_compute_msa],
					  gct.DOWNLOADSCHEMA_MODULE_NAME: [gct.DOWNLOADSCHEMA_MODULE_SHORT_DESCRIPTION, run_download_schema],
					  gct.UPLOADSCHEMA_MODULE_NAME: [gct.UPLOADSCHEMA_MODULE_SHORT_DESCRIPTION, run_upload_schema],
					  gct.SYNCSCHEMA_MODULE_NAME: [gct.SYNCSCHEMA_MODULE_SHORT_DESCRIPTION, run_synchronize_schema],
					  gct.NSSTATS_MODULE_NAME: [gct.NSSTATS_MODULE_SHORT_DESCRIPTION, run_stats_requests]}

	if len(sys.argv) > 1 and sys.argv[1] in gct.VERSION_TRIGGERS:
		print(f"chewBBACA version: {gct.__version__}")
		# Exit after printing version
		sys.exit(0)

	print_initial_header()

	# Display help message if selected process is not valid
	if len(sys.argv) == 1 or sys.argv[1] not in functions_info or sys.argv[1] in gct.HELP_TRIGGERS:
		exit_code = 0
		# Detect if user passed module name that does not exist
		if len(sys.argv) > 1 and sys.argv[1] not in gct.HELP_TRIGGERS:
			print(f"No module named {sys.argv[1]}.\n")
			exit_code = 2
		print("USAGE: chewBBACA.py [module] -h, --help\n")
		print("Select one of the following modules:")
		for f in functions_info:
			print('{0}: {1}'.format(f, functions_info[f][0]))
		sys.exit(exit_code)

	# Check python version
	python_uptodate, python_version = gval.validate_python_version(gct.MIN_PYTHON)
	if not python_uptodate:
		sys.exit(gct.PYTHON_VERSION.format(python_version, gct.MIN_PYTHON))

	# Trigger module help message if no arguments are provided
	if len(sys.argv) == 2 and sys.argv[1] in functions_info:
		sys.argv.append(gct.HELP_TRIGGERS[0])

	process = sys.argv[1]
	# Pop the module name from the list of arguments because it is not defined as an argument in the parsers
	sys.argv.pop(1)
	# Call the selected module
	functions_info[process][1]()


if __name__ == "__main__":

	main()
