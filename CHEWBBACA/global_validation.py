#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions/classes related to the validation of the
arguments passed to chewBBACA's modules.

Code documentation
------------------
"""


import sys
import shutil
import platform
import multiprocessing
from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator)
from pydantic_core import PydanticUndefined

try:
	import global_constants as gct
	from utils import (blast_wrapper as bw,
					   file_operations as fo,
					   fasta_operations as fao,
					   iterables_manipulation as im,
					   pyrodigal_gene_prediction as pgp)
except ModuleNotFoundError:
	from CHEWBBACA import global_constants as gct
	from CHEWBBACA.utils import (blast_wrapper as bw,
	                             file_operations as fo,
								 fasta_operations as fao,
								 iterables_manipulation as im,
								 pyrodigal_gene_prediction as pgp)


def validate_python_version(minimum_version):
	"""Validate Python version used to run chewBBACA.

	Parameters
	----------
	minimum_version : tuple
		A tuple with the Python version as (MAJOR, MINOR, PATCH) according 
		to the rules of Semanting Versioning (https://semver.org/).

	Returns
	-------
	valid : bool
		True if the Python version is compatible with chewBBACA, False otherwise.
	version : str
		Python version in the format "MAJOR.MINOR.PATCH".

	Raises
	------
	SystemExit
		- If the Python version does not meet minimum requirements
		or it was not possible to determine/detect a version.
	"""
	version = platform.python_version()

	valid = tuple(map(int, version.split('.'))) >= minimum_version[0]

	return valid, version


def filter_inputs_extension(input_files, extensions=gct.FASTA_EXTENSIONS):
	"""Determine if files end with any of the supported file extensions.

	Parameters
	----------
	input_files : list
		List with paths to input files.
	extensions : list
		List of supported file extensions.

	Returns
	-------
	input_files : list
		List with paths to input files that have a supported file extension.
	
	Raises
	------
	SystemExit
		- If any of the input files has a non-supported file extension.
	"""
	# Need to verify if files end with any of the accepted file extensions
	valid_extension, invalid_extension = fo.filter_by_extension(input_files, extensions)
	if len(invalid_extension) > 0:
		# Join list of files with non-supported file extensions
		invalid_joined = "\n".join(invalid_extension)
		sys.exit(gct.INVALID_EXTENSION_EXCEPTION.format(invalid_joined, extensions))

	return input_files


def inputs_exist(input_files):
	"""Determine if all input files exist.

	Parameters
	----------
	input_files : list
		List with paths to input files.
	
	Returns
	-------
	input_files : list
		List with paths to input files that exist.

	Raises
	------
	SystemExit
		- If any of the input files does not exist.
	"""
	# Check that all files exist
	missing_inputs = [file for file in input_files if fo.exists(file) is False]
	if len(missing_inputs) > 0:
		missing_joined = "\n".join(missing_inputs)
		sys.exit(gct.MISSING_INPUTS_EXCEPTION.format(missing_joined))

	return input_files


def validate_inputs_fastas(input_files):
	"""Determine if all input files are in a valid FASTA format.

	Parameters
	----------
	input_files : list
		List with paths to input files.
	
	Returns
	-------
	input_files : list
		List with paths to input files that are in a valid FASTA format.

	Raises
	------
	SystemExit
		- If any of the input files is not in a valid FASTA format.
	"""
	# Input FASTA files must have a valid format
	fasta, non_fasta = fao.filter_non_fasta(input_files)
	if len(non_fasta) > 0:
		# Join files that are not in a valid FASTA format
		non_fasta_joined = "\n".join(non_fasta)
		sys.exit(gct.NON_FASTA_EXCEPTION.format(non_fasta_joined))

	return input_files


def check_unique_prefixes(input_files):
	"""Check if all input files have a unique identifier.

	Parameters
	----------
	input_list : list
		List with paths to input files.

	Returns
	-------
	input_files : list
		List with paths to input files that have a unique identifier.

	Raises
	------
	SystemExit
		- If there are multiple files sharing the same unique prefix.
	"""
	prefixes = fo.get_file_prefixes(input_files)
	# Detect if some inputs share the same unique prefix
	if len(set(prefixes)) < len(input_files):
		repeated_basenames = [f"{k}: {', '.join(v)}" for k, v in prefixes.items() if len(v) > 1]
		repeated_basenames = [",".join(l) for l in repeated_basenames]
		repeated_basenames = "\n".join(repeated_basenames)
		sys.exit(gct.INPUTS_SHARE_PREFIX.format(repeated_basenames))

	return input_files


def check_blanks(input_files):
	"""Check if input files do not include blank spaces in the filename.

	Parameters
	----------
	input_files : list
		List with paths to input files.

	Returns
	-------
	input_files : list
		List with paths to input files that do not include blank spaces in the filename.

	Raises
	------
	SystemExit
		- If there are blank spaces in any of the filenames.
	"""
	basenames = [fo.file_basename(file) for file in input_files]
	include_blanks = [name for name in basenames if ' ' in name]
	if len(include_blanks) > 0:
		include_blanks_joined = "\n".join(include_blanks)
		sys.exit(gct.INPUTS_INCLUDE_BLANKS.format(include_blanks_joined))

	return input_files


def check_blast(blast_path, blastp_alias=gct.BLASTP_ALIAS):
	"""Determine if BLAST executables are available in the system.

	Parameters
	----------
	blast_path : str
		Path to the folder containing BLAST executables.
	blastp_alias : str
		Alias used to call BLASTp.

	Returns
	-------
	blast_path : str
		Path to the folder containing BLAST executables.

	Raises
	------
	SystemExit
		- If BLAST executables are not available in the system.
	"""
	# User did not provide a path, try to determine path to BLAST executables using BLASTp alias
	if not blast_path:
		blast_path = fo.get_parent_directory(shutil.which(blastp_alias))

	if blast_path is None or fo.exists(blast_path) is False:
		sys.exit(gct.BLAST_MISSING)

	return blast_path


def check_blast_version(blast_path, blastp_alias=gct.BLASTP_ALIAS, major=gct.BLAST_MAJOR, minor=gct.BLAST_MINOR):
	"""Determine if the BLAST version meets minimum requirements.

	Parameters
	----------
	blast_path : str
		Path to the folder containing BLAST executables.
	blastp_alias : str
		Alias used to call BLASTp.
	major : int
		Minimum required MAJOR version of BLAST.
	minor : int
		Minimum required MINOR version of BLAST.
	
	Returns
	-------
	blast_path : str
		Path to the folder containing BLAST executables.

	Raises
	------
	SystemExit
		- If the BLAST version does not meet minimum requirements.
	"""
	# Create path to BLASTp executable
	blastp_path = fo.join_paths(blast_path, [blastp_alias])
	# Get BLAST version
	blast_version = bw.get_blast_version(blastp_path)
	# Determine if BLAST version meets minimum requirements
	if blast_version['MAJOR'] < major or (blast_version['MAJOR'] >= major and blast_version['MINOR'] < minor):
		sys.exit(gct.BLAST_INVALID_VERSION)

	return blast_path


def input_is_fasta(input_path):
	"""Determine if the provided input is a single FASTA file.

	Parameters
	----------
	input_path : str
		Path to the input file.

	Returns
	-------
	input_path : str
		Path to the input file if it is a single FASTA file.

	Raises
	------
	SystemExit
		- If the provided input is a FASTA file.
	"""
	# Check if input is a single FASTA file
	if fo.is_file(input_path)[0]:
		if fao.validate_fasta(input_path):
			sys.exit(gct.FASTA_INPUT_EXCEPTION.format(input_path))

	return input_path


def add_ptf_genetic_code(ptf_path, genetic_code):
	"""Determine if the provided Pyrodigal training file has a different genetic code than the one provided by the user.

	Parameters
	----------
	ptf_path : str
		Path to the Pyrodigal training file.
	genetic_code : int
		Genetic code provided by the user.

	Returns
	-------
	genetic_code : int
		Genetic code that will be used for Pyrodigal gene prediction.
	"""
	# Get translation table used to create training file
	ptf_genetic_code = pgp.read_training_file(ptf_path).translation_table
	if ptf_genetic_code != genetic_code:
		print(gct.FORCE_TRAINING_FILE.format(genetic_code, ptf_genetic_code))
		genetic_code = ptf_genetic_code

	return genetic_code


def create_training_file(training_reference, output_directory, translation_table):
	"""
	"""
	return pgp.create_training_file(training_reference, output_directory, translation_table)


def check_meta(pyrodigal_mode, pyrodigal_training_file):
	"""Determine if the provided Pyrodigal mode is compatible with the provided training file.

	Parameters
	----------
	pyrodigal_mode : str
		Pyrodigal mode provided by the user.
	pyrodigal_training_file : str
		Path to the Pyrodigal training file provided by the user.

	Returns
	-------
	pyrodigal_mode : str
		Pyrodigal mode that will be used for Pyrodigal gene prediction.

	Raises
	------
	SystemExit
		- If the provided Pyrodigal mode is "meta" and a training file was provided.
	"""
	if pyrodigal_mode == "meta" and pyrodigal_training_file is not None:
		sys.exit(gct.PYRODIGAL_META_NOPTF)

	return pyrodigal_mode


def validate_multiple_choices(choices, valid_options, exception_str):
	"""Validate if all provided choices are valid options.

	Parameters
	----------
	choices : list
		List of choices provided by the user.
	valid_options : list
		List of valid options.
	exception_str : str
		Exception string to be printed if any of the provided choices is not a valid option.

	Returns
	-------
	choices : list
		List of choices that are valid options.

	Raises
	------
	SystemExit
		- If any of the provided choices is not a valid option.
	"""
	if not all([option in valid_options for option in choices]):
		sys.exit(exception_str.format(valid_options))

	return choices


def verify_cpu_usage(cpu_value):
	"""Verify if the cores/threads value does not exceed available resources.

	Parameters
	----------
	cpu_value : int
		Value provided for the number of CPU cores/threads.

	Returns
	-------
	cpu_to_use : int
		Value of CPU cores/threads that will be used after
		determining if the provided value was safe.
	"""
	# Determine number of cores/threads available in the system
	system_ncores = multiprocessing.cpu_count()
	if cpu_value > system_ncores:
		sys.exit(gct.EXCEEDED_CPU.format(cpu_value, system_ncores))
	elif cpu_value == system_ncores:
		print(gct.CPU_VALUE_WARNING.format(cpu_value, system_ncores))

	return cpu_value


def parse_parameter_string(input_string):
	"""Parse a string of parameter=argument pairs.

	Parameters
	----------
	input_string : list
		List of strings with parameter=argument pairs.

	Returns
	-------
	config_args : dict
		Dictionary with parameter names as keys and argument values as values.
	"""
	config_args = {}
	for v in input_string:
		parameter, argument = v.split("=")
		parameter = parameter.replace("-", "_")
		# Split argument value in by "," to get multiple values
		if "," in argument:
			argument = argument.split(",")
		# Add parameter name and argument value to config dictionary
		config_args[parameter] = argument

	return config_args


def validate_choice(choice, valid_options, exception_str):
	"""Validate if the provided choice is a valid option.

	Parameters
	----------
	choice : str
		Choice provided by the user.
	valid_options : list
		List of valid options.
	exception_str : str
		Exception string to be printed if the provided choice is not a valid option.

	Returns
	-------
	choice : str
		The validated choice.

	Raises
	------
	SystemExit
		- If the provided choice is not a valid option.
	"""
	if choice not in valid_options:
		sys.exit(exception_str)

	return choice


def parentdir_exists(input_path):
	"""Determine if the parent directory of the provided path exists.
	
	Parameters
	----------
	input_path : str
		Path for which to check the parent directory.

	Returns
	-------
	input_path : str
		The input path if the parent directory exists.

	Raises
	------
	ValueError
		- If the parent directory for the provided path does not exist.
	"""
	parent_dir = fo.get_parent_directory(input_path)
	if not fo.is_directory(parent_dir)[0]:
		raise ValueError(f"Parent directory for '{input_path}' does not exist.")

	return input_path


def input_paths_exist(input_path, function, exception_str):
	"""Determine if the provided input path(s) exist.

	Parameters
	----------
	input_path : str or list
		Path(s) for which to check existence.
	function : callable
		Function to use for checking existence.
	exception_str : str
		Exception string to be printed if any path does not exist.

	Returns
	-------
	input_path : str or list
		The input path(s) if all exist.

	Raises
	------
	SystemExit
		- If any of the provided paths do not exist.
	"""
	# Convert into list if it is a single path/str
	path_list = input_path if type(input_path) is list else [input_path]
	for path in path_list:
		if not function(path):
			sys.exit(exception_str.format(path))

	return input_path


def create_output_directory(output_path):
	"""Create the output directory if it does not exist.
	
	Parameters
	----------
	output_path : str
		Path to the output directory.

	Returns
	-------
	output_path : str
		Path to the output directory.

	Raises
	------
	SystemExit
		- If the output directory already exists.
	"""
	created, output_path = fo.create_directory(output_path)
	if not created:
		sys.exit(gct.OUTPUT_DIRECTORY_EXISTS)

	return output_path


def check_ptf_tref_conflict(training_file, training_reference):
	"""Check if both a Pyrodigal training file and a training reference were provided.
	
	Parameters
	----------
	training_file : str
		Path to the Pyrodigal training file.
	training_reference : str
		Path to the training reference.

	Returns
	-------
	None

	Raises
	------
	SystemExit
		- If both a Pyrodigal training file and a training reference were provided.
	"""
	if training_file and training_reference:
		sys.exit(gct.CANNOT_PROVIDE_PTF_AND_TREF)


def schema_includes_fasta(schema_directory):
	"""Check if the schema directory includes FASTA files.

	Parameters
	----------
	schema_directory : str
		Path to the schema directory.

	Returns
	-------
	schema_directory : str
		Path to the schema directory.

	Raises
	------
	SystemExit
		- If no FASTA files are found in the schema directory.
	"""
	schema_files = fo.listdir(schema_directory)
	# Check if the folder includes FASTA files
	fasta_extension, _ = fo.filter_by_extension(schema_files, gct.FASTA_EXTENSIONS)
	if len(fasta_extension) == 0:
		sys.exit(gct.MISSING_SCHEMA_FASTAS)

	return schema_directory


def schema_check_genes_list(schema_directory):
	"""Check if the schema directory includes a genes list file and the loci IDs match the FASTA files.

	Parameters
	----------
	schema_directory : str
		Path to the schema directory.

	Returns
	-------
	schema_directory : str
		Path to the schema directory.

	Raises
	------
	SystemExit
		- If the genes list file is missing or the loci IDs do not match the FASTA files.
	"""
	schema_files = fo.listdir(schema_directory)
	# Check if schema includes the .genes_list file
	if gct.GENE_LIST_BASENAME not in schema_files:
		sys.exit(gct.MISSING_GENE_LIST)

	# Verify that loci FASTA files match list of loci in the .genes_list file
	schema_files_fasta = fo.filter_by_extension(schema_files, [".fasta"])[0]
	schema_files_basenames = [fo.file_basename(file, False) for file in schema_files_fasta]
	loci_list = fo.pickle_loader(fo.join_paths(schema_directory, [gct.GENE_LIST_BASENAME]))

	# Determine if the schema or the list of loci have different sets of loci IDs
	schema_exclusive = set(schema_files_basenames) - set(loci_list)
	loci_list_exclusive = set(loci_list) - set(schema_files_basenames)
	if len(schema_exclusive) > 0 or len(loci_list_exclusive) > 0:
		sys.exit(gct.SCHEMA_LOCI_GENE_LIST_NOMATCH.format(len(schema_exclusive), im.join_list(schema_exclusive, ", "), len(loci_list_exclusive), im.join_list(loci_list_exclusive, ", ")))

	return schema_directory


def schema_check_config(schema_directory):
	"""Check if the schema directory includes the configuration file.

	Parameters
	----------
	schema_directory : str
		Path to the schema directory.

	Returns
	-------
	schema_directory : str
		Path to the schema directory.

	Raises
	------
	SystemExit
		- If the configuration file is missing.
	"""
	schema_files = fo.listdir(schema_directory)
	# Check if schema includes the .schema_config file
	if gct.SCHEMA_CONFIG_BASENAME not in schema_files:
		sys.exit(gct.SCHEMA_MISSING_FILES.format(gct.SCHEMA_CONFIG_BASENAME))

	return schema_directory


def schema_check_short(schema_directory):
	"""Check if the schema directory includes the "short" directory.

	Parameters
	----------
	schema_directory : str
		Path to the schema directory.

	Returns
	-------
	schema_directory : str
		Path to the schema directory.

	Raises
	------
	SystemExit
		- If the "short" directory is missing.
	"""
	schema_files = fo.listdir(schema_directory)
	# Check if the "short" directory exists
	if "short" not in schema_files:
		sys.exit(gct.SCHEMA_MISSING_FILES.format("short"))

	return schema_directory


def check_argument_conflict(arg_name, arg_value, schema_config, force_continue):
	"""Check if the provided argument value conflicts with the values stored in the schema configuration.
	
	Parameters
	----------
	arg_name : str
		Name of the argument to check for conflicts.
	arg_value : any
		Value of the argument to check for conflicts.
	schema_config : dict
		Dictionary containing the schema configuration with the argument values used with the schema.
	force_continue : bool
		Flag indicating whether to force continuation despite conflicts.

	Returns
	-------
	arg_value : any
		The validated argument value.
	schema_config : dict
		Updated schema configuration with the argument value added if it was not already present.
	"""
	if arg_value not in schema_config[arg_name]:
		print(gct.CONFLICTING_ARGUMENTS.format(arg_name, schema_config[arg_name]))
		if not force_continue:
			answer = fo.input_timeout(gct.ARGS_DIFFER_PROMPT, gct.PROMPT_TIMEOUT)
		else:
			answer = 'yes'

		if answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			schema_config[arg_name].append(arg_value)

	return arg_value, schema_config


def contains_results(results_directory, results_files_basenames):
	"""Check if the results directory includes the expected results file.

	Parameters
	----------
	results_directory : str
		Path to the results directory.

	Returns
	-------
	results_directory : str
		Path to the results directory.

	Raises
	------
	SystemExit
		- If the expected results file is missing.
	"""
	# List files in input direcoty
	results_files = fo.listdir(results_directory)
	# Check if folder includes the results_alleles.tsv file with allelic profiles
	if any([basename not in results_files for basename in results_files_basenames]):
		sys.exit(gct.MISSING_RESULTS_FILES.format(", ".join(results_files_basenames)))

	return results_directory


# Define reusable annotated types/Fields for argument validation in multiple modules
OutputDirectory = Annotated[str,
                            Field(description=""),
							AfterValidator(parentdir_exists),
							AfterValidator(create_output_directory)]

InputPath = Annotated[str,
					   AfterValidator(partial(input_paths_exist,
					                          function=fo.exists,
											  exception_str=gct.MISSING_INPUT_ARG)),
					   AfterValidator(input_is_fasta),
					   AfterValidator(fo.list_input_files),
					   AfterValidator(filter_inputs_extension),
					   AfterValidator(inputs_exist),
					   AfterValidator(validate_inputs_fastas),
					   AfterValidator(check_unique_prefixes),
					   AfterValidator(check_blanks)]

ChewieSchema = Annotated[str,
						 AfterValidator(partial(input_paths_exist,
						                        function=fo.exists,
												exception_str=gct.SCHEMA_PATH_MISSING)),
						 AfterValidator(schema_includes_fasta),
						 AfterValidator(schema_check_genes_list),
						 AfterValidator(schema_check_config),
						 AfterValidator(schema_check_short)]

ExternalSchema = Annotated[str,
							AfterValidator(partial(input_paths_exist,
							                       function=fo.exists,
												   exception_str=gct.SCHEMA_PATH_MISSING)),
							AfterValidator(schema_includes_fasta)]

LociList = Annotated[str,
					 AfterValidator(partial(input_paths_exist,
					                        function=fo.exists,
											exception_str=gct.MISSING_LOCI_LIST)),
					 AfterValidator(input_is_fasta),
					 AfterValidator(fo.list_input_files)]

BLASTScoreRatio = Annotated[float,
							Field(default=gct.BSR_DEFAULT,
							      ge=gct.BSR_MIN,
								  le=gct.BSR_MAX)]

MinimumLength = Annotated[int,
						  Field(default=gct.MSL_DEFAULT,
						        ge=gct.MSL_MIN,
								le=gct.MSL_MAX)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT),
							 AfterValidator(partial(validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

SizeThreshold = Annotated[float | None,
						  Field(default=gct.ST_DEFAULT,
						        ge=gct.ST_MIN,
								le=gct.ST_MAX)]

MinimumLengthValues = Annotated[list,
                                Field(default=[])]

SizeThresholdValues = Annotated[list,
                                Field(default=[])]

TranslationTableValues = Annotated[list,
                                   Field(default=[])]

BlastScoreRatioValues = Annotated[list,
                                  Field(default=[])]

BLASTPath = Annotated[str | None, # Need to specify that the argument value may not be provided
                      Field(default=None, # Need to set default to None to really set arg value as None
							validate_default=True, # Needs to be validated even if it is None
							description=gct.BLAST_PATH_DESCRIPTION),
					  AfterValidator(check_blast),
					  AfterValidator(check_blast_version)]

CDSInput = Annotated[bool,
                     Field(default=False,
                           description=gct.CDS_INPUT_DESCRIPTION)]

NoCleanup = Annotated[bool,
                        Field(default=False,
                              description=gct.NO_CLEANUP_DESCRIPTION)]

NoCDSRemaining = Annotated[bool,
                        Field(default=False,
                              description=gct.NO_CDS_RENAMING_DESCRIPTION)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=multiprocessing.cpu_count()),
					 AfterValidator(verify_cpu_usage)]

WordSize = Annotated[int,
					 Field(default=gct.WORD_SIZE_DEFAULT,
					       ge=gct.WORD_SIZE_MIN,
						   le=gct.WORD_SIZE_MAX)]

WindowSize = Annotated[int,
					   Field(default=gct.WORD_SIZE_DEFAULT,
					         ge=gct.WINDOW_SIZE_MIN,
							 le=gct.WINDOW_SIZE_MAX)]

ClusteringSimilarity = Annotated[float,
								 Field(default=gct.CLUSTERING_SIMILARITY_DEFAULT,
								       ge=gct.CLUSTERING_SIMILARITY_MIN,
									   le=gct.CLUSTERING_SIMILARITY_MAX)]

RepresentativeFilter = Annotated[float,
								 Field(default=gct.REPRESENTATIVE_FILTER_DEFAULT,
								       ge=gct.REPRESENTATIVE_FILTER_MIN,
									   le=gct.REPRESENTATIVE_FILTER_MAX)]

IntraFilter = Annotated[float,
						Field(default=(gct.INTRA_CLUSTER_DEFAULT),
						      ge=gct.INTRA_CLUSTER_MIN,
							  le=gct.INTRA_CLUSTER_MAX)]

Annotations = Annotated[str,
                        AfterValidator(partial(input_paths_exist,
                                               function=fo.exists,
											   exception_str=gct.MISSING_ANNOTATIONS_FILE))]

BooleanFlag = Annotated[bool,
                        Field(default=False)]

NSCONFIG = Annotated[bool,
                     Field(default=False)]

# Nedd to add a AfterValidator to get and validate results files
ResultsFiles = Annotated[str,
                         AfterValidator(partial(input_paths_exist,
                                                function=fo.exists,
												exception_str=gct.MISSING_RESULTS)),
						 AfterValidator(partial(contains_results,
						                        results_files_basenames=[gct.ALLELIC_PROFILES_BASENAME]))]

# Should implement a validator for the profile filr format
AllelicProfiles = Annotated[str,
                            AfterValidator(partial(input_paths_exist,
                                                   function=fo.exists,
												   exception_str=gct.PROFILES_FILE_MISSING))]


# Use a separate model to validate clustering arguments passed to modules that use clustering
class ClusteringArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	word_size: WordSize
	window_size: WindowSize
	clustering_similarity: ClusteringSimilarity
	representative_filter: RepresentativeFilter
	intra_filter: IntraFilter


# Field definition for the validated clustering arguments
ValidatedClusteringOptions = Annotated[ClusteringArgs | None, Field(default=None)]


# Create str with the defautl values for the gene prediction parameter=argument pairs
# Added here to avoid circular imports
def create_default_clustering_str():
	default_clustering_config = ClusteringArgs.model_validate({}).model_dump()
	default_clustering_options = []
	for parameter, argument in default_clustering_config.items():
		if type(argument) == list:
			default_clustering_options.append(f'"{parameter}={','.join(argument)}"')
		else:
			default_clustering_options.append(f'"{parameter}={argument}"')

	default_clustering_str = ' '.join(default_clustering_options)

	return default_clustering_str


ClusteringOptions = Annotated[str,
                              Field(default=create_default_clustering_str(),
                                    description=gct.CLUSTERING_OPTIONS_DESCRIPTION),
                              AfterValidator(parse_parameter_string)]


def get_field_help(model, field_name):
	"""
	"""
	field = model.model_fields[field_name]
	description = field.description or ""

	# Append default if a valid defautl value is set
	if field.default is not PydanticUndefined and field.default is not None:
		return f"{description} (default: {field.default})."

	return description


def get_system_cpu_cores():
	""""""
	return multiprocessing.cpu_count()
