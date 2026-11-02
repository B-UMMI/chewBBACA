#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions/classes related to the validation
of the arguments passed to chewBABCA's modules.

Code documentation
------------------
"""


import os
import re
import sys
import zlib
import shutil
import hashlib
import argparse
import platform
import subprocess
import multiprocessing
from typing import Annotated
from functools import partial

import pandas as pd
from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  TypeAdapter,
					  AfterValidator,
					  BeforeValidator,
					  model_validator)

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   process_datetime as pd,
					   fasta_operations as fao,
					   chewiens_requests as cr,
					   iterables_manipulation as im,
					   pyrodigal_gene_prediction as pgp)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
	                             file_operations as fo,
	                             process_datetime as pd,
								 fasta_operations as fao,
								 chewiens_requests as cr,
								 iterables_manipulation as im,
								 pyrodigal_gene_prediction as pgp)


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

	# prog is the name of the program 'ex: chewBBACA.py'
	def __init__(self, prog, indent_increment=1, max_help_position=56, width=100):
		super().__init__(prog, indent_increment, max_help_position, width)

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


# def validate_ns_url(arg):
# 	"""Verify if the Chewie-NS URL passed to chewBBACA is valid.

# 	Parameters
# 	----------
# 	arg : str
# 		Identifier of the Chewie-NS instance or the URL
# 		to a instance of Chewie-NS.

# 	Returns
# 	-------
# 	ns_url : str
# 		URL to connect to the instance of Chewie-NS.

# 	Raises
# 	------
# 	SystemExit
# 		- If it is not possible to connect to the
# 		chewie-NS instance.
# 	"""
# 	if arg in ct.HOST_NS:
# 		ns_url = ct.HOST_NS[arg]
# 	else:
# 		ns_url = arg

# 	# sync schema has None by default to get ns_url in schema URI
# 	if ns_url is not None:
# 		# check if server is up
# 		conn = cr.check_connection(ns_url)
# 		if conn is False:
# 			sys.exit(ct.NS_CANNOT_CONNECT.format(ns_url))

# 	return ns_url


def validate_python_version(minimum_version):
	"""Validate Python version used to run chewBBACA.

	Parameters
	----------
	minimum_version : tuple
		A tuple with the Puthon version as (MAJOR, MINOR, PATCH).
		According to the rules of Semanting Versioning
		(https://semver.org/).

	Returns
	-------
	python_version : str
		Python version in format "MAJOR.MINOR.PATCH".

	Raises
	------
	SystemExit
		- If the Python version does not meet minimum requirements
		or it was not possible to determine/detect a version.
	"""
	version = platform.python_version()

	valid = tuple(map(int, version.split('.'))) >= minimum_version[0]

	return valid, version


def is_exe(fpath):
	"""Determine if path points to a file and if the file is an executable.

	Parameters
	----------
	fpath : str
		Path to a file.

	Returns
	-------
	True if the file exists and is executable, False otherwise.
	"""
	return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def get_blast_version(blastp_path):
	"""Determines BLAST version.

	Parameters
	----------
	blastp_path : str
		Path to the BLASTp executable.

	Returns
	-------
	version : dict or NoneType
		Dictionary with the BLAST MAJOR and MINOR versions or
		NoneType if it was not possible to determine the BLAST
		version.
	"""
	# Try to get BLAST version using BLASTp
	try:
		proc = subprocess.Popen([blastp_path, '-version'],
								stdout=subprocess.PIPE,
								stderr=subprocess.PIPE)
		stdout, stderr = proc.communicate()
		# Process the version string
		version_string = stdout.decode('utf8')
		version_pattern = r'^blastp:\s(?P<MAJOR>\d+).(?P<MINOR>\d+).(?P<REV>\d+).*'
		blast_version_pat = re.compile(version_pattern)
		match = blast_version_pat.search(version_string)
		# Got something
		if match is not None:
			version = {k: int(v) for k, v in match.groupdict().items()}
		# No matches
		else:
			version = None
	except:
		version = None

	return version


def get_file_prefixes(path_list):
	"""Determine the file prefix for each file in a list of file paths.

	Parameters
	----------
	path_list : list
		List with file paths.

	Returns
	-------
	prefixes :  dict
		Dictionary with file prefixes as keys and file basenames
		as values.
	"""
	basenames = {fo.file_basename(file, False): file for file in path_list}
	basename_counts = {}
	for k, v in basenames.items():
		if k not in basename_counts:
			basename_counts[k] = [v]
		else:
			basename_counts[k] += [v]

	return basename_counts


def list_input_files(input_path):
	""""""
	if os.path.isfile(input_path):
		# Read list of input files
		input_files = [line[0] for line in fo.read_tabular(input_path)]
	# Input path is for a directory
	elif os.path.isdir(input_path):
		# List absolute paths
		# Only keep paths to files
		input_files = [file for file in fo.listdir(input_path, full_path=True) if os.path.isdir(file) is False]

	return input_files


def filter_inputs_extension(input_files, extensions=ct.FASTA_EXTENSIONS):
	""""""
	# Need to verify if files end with any of the accepted file extensions
	valid_extension, invalid_extension = fo.filter_by_extension(input_files, extensions)
	if len(invalid_extension) > 0:
		# Join list of files with non-supported file extensions
		invalid_joined = "\n".join(invalid_extension)
		sys.exit(ct.INVALID_EXTENSION_EXCEPTION.format(invalid_joined, extensions))

	return input_files


def inputs_exist(input_files):
	""""""
	# Check that all files exist
	missing_inputs = [file for file in input_files if fo.exists(file) is False]
	if len(missing_inputs) > 0:
		missing_joined = "\n".join(missing_inputs)
		sys.exit(ct.MISSMISSING_INPUTS_EXCEPTIONING_INPUTS.format(missing_joined))

	return input_files


def validate_inputs_fastas(input_files):
	""""""
	# Input FASTA files must have a valid format
	fasta, non_fasta = fao.filter_non_fasta(input_files)
	if len(non_fasta) > 0:
		# Join files that are not in a valid FASTA format
		non_fasta_joined = "\n".join(non_fasta)
		sys.exit(ct.NON_FASTA_EXCEPTION.format(non_fasta_joined))

	return input_files


def check_unique_prefixes(input_files):
	"""Check if all input files have an unique identifier.

	Parameters
	----------
	input_list : str
		Path to file that contains the list of paths to input files.

	Returns
	-------
	False if there are no input files sharing the same identifier.

	Raises
	------
	SystemExit
		- If there are multiple files sharing the same prefix.
	"""
	prefixes = get_file_prefixes(input_files)
	# Detect if some inputs share the same unique prefix
	if len(set(prefixes)) < len(input_files):
		repeated_basenames = [f"{k}: {', '.join(v)}" for k, v in prefixes.items() if len(v) > 1]
		repeated_basenames = [",".join(l) for l in repeated_basenames]
		repeated_basenames = "\n".join(repeated_basenames)
		sys.exit(ct.INPUTS_SHARE_PREFIX.format(repeated_basenames))

	return input_files


def check_blanks(input_files):
	"""Check if input files do not include blank spaces in the filename.

	Parameters
	----------
	input_list : str
		Path to file that contains the list of paths to input files.

	Returns
	-------
	False if there are no blank spaces in the filenames.

	Raises
	------
	SystemExit
		- If there are blank spaces in any of the filenames.
	"""
	basenames = [fo.file_basename(file) for file in input_files]
	include_blanks = [name for name in basenames if ' ' in name]
	if len(include_blanks) > 0:
		include_blanks_joined = "\n".join(include_blanks)
		sys.exit(ct.INPUTS_INCLUDE_BLANKS.format(include_blanks_joined))

	return input_files


def check_blast(blast_path, blastp_alias=ct.BLASTP_ALIAS):
	"""
	"""
	# User did not provide a path, try to determine path to BLAST executables using BLASTp alias
	if not blast_path:
		blast_path = fo.get_parent_directory(shutil.which(blastp_alias))

	if blast_path is None or fo.exists(blast_path) is False:
		sys.exit(ct.BLAST_MISSING)

	return blast_path


def check_blast_version(blast_path, blastp_alias=ct.BLASTP_ALIAS, major=ct.BLAST_MAJOR, minor=ct.BLAST_MINOR):
	"""
	"""
	# Create path to BLASTp executable
	blastp_path = fo.join_paths(blast_path, [blastp_alias])
	# Get BLAST version
	blast_version = get_blast_version(blastp_path)
	# Determine if BLAST version meets minimum requirements
	if blast_version['MAJOR'] < major or (blast_version['MAJOR'] >= major and blast_version['MINOR'] < minor):
		sys.exit(ct.BLAST_INVALID_VERSION)

	return blast_path


def check_augustus(augustus_path, augustus_alias=ct.AUGUSTUS_ALIAS):
	"""
	"""
	# User did not provide a path, try to determine path to the AUGUSTUS executable using its alias
	if not augustus_path:
		augustus_path = fo.get_parent_directory(shutil.which(augustus_alias))

	if augustus_path is None or fo.exists(augustus_path) is False:
		sys.exit(ct.AUGUSTUS_MISSING)

	return augustus_path


def get_augustus_species_list(augustus_exe):
	"""Get the list of species' models supported by AUGUSTUS.

	Parameters
	----------
	augustus_alias : str
		Alias used to call AUGUSTUS.

	Returns
	-------
	version : str
	"""
	# Try to get the list os species' models supported by AUGUSTUS'
	proc = subprocess.Popen([augustus_exe, '--species=help'],
							stdout=subprocess.PIPE,
							stderr=subprocess.PIPE,
							text=True)
	stdout, stderr = proc.communicate()
	# Process the list of species printed to stdout
	species_list = {}
	# List of species is printed to stderr
	i = 0
	for line in stderr.split("\n"):
		if "|" in line:
			i += 1
			# Skip the table header and separator
			if i > 2:
				species_id, species_name = line.split("|")
				species_id = species_id.strip()
				species_name = species_name.strip()
				species_list[species_name] = species_id

	return species_list


def validate_augustus_species(species_id, augustus_path, augustus_alias=ct.AUGUSTUS_ALIAS):
	""""""
	augustus_exe = fo.join_paths(augustus_path, [augustus_alias])
	species_list = get_augustus_species_list(augustus_exe)

	if species_id not in species_list.values():
		# Create string with list of options
		species_list_str = [": ".join([k, v]) for k, v in species_list.items()]
		species_list_str = "\n".join(species_list_str)
		sys.exit(ct.AUGUSTUS_INVALID_SPECIES.format(species_list_str))

	return species_id


def input_is_fasta(input_path):
	""""""
	# Check if input is a single FASTA file
	if fo.is_file(input_path)[0]:
		if fao.validate_fasta(input_path):
			sys.exit(ct.FASTA_INPUT_EXCEPTION.format(input_path))

	return input_path


def add_ptf_genetic_code(ptf_path, genetic_code):
	"""
	"""
	# Get translation table used to create training file
	ptf_genetic_code = pgp.read_training_file(ptf_path).translation_table
	if ptf_genetic_code != genetic_code:
		print("Provided a Pyrodigal training file. Ignoring the translation table "
				f"value previously set ({genetic_code}) and using the genetic "
				f"code defined in the training file ({ptf_genetic_code}).")
		genetic_code = ptf_genetic_code

	return genetic_code


def check_meta(pyrodigal_mode, pyrodigal_training_file):
	"""
	"""
	if pyrodigal_mode == 'meta' and pyrodigal_training_file is not None:
		sys.exit(ct.PYRODIGAL_META_NOPTF)

	return pyrodigal_mode


def validate_multiple_choices(choices, valid_options, exception_str):
	""""""
	if not all([option in valid_options for option in choices]):
		sys.exit(exception_str.format(valid_options))

	return choices


def verify_cpu_usage(arg):
	"""Verify if the cores/threads value does not exceed available resources.

	Parameters
	----------
	cpu_to_use : int
		Value provided for the number of CPU cores/threads.

	Returns
	-------
	cpu_to_use : int
		Value of CPU cores/threads that will be used after
		determining if the provided value was safe.
	"""
	system_ncores = multiprocessing.cpu_count()
	if arg == system_ncores:
		print(ct.CPU_VALUE_WARNING.format(system_ncores))

	return arg


def parse_parameter_string(input_string):
	"""
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
	""""""
	if choice not in valid_options:
		sys.exit(exception_str)

	return choice


# Check if parent directory for the provided output_dir exists
# field_validator is used to check before the class is instantiated
def parentdir_exists(input_path):
	parent_dir = os.path.dirname(input_path)
	if not os.path.isdir(parent_dir):
		raise ValueError(f"Parent directory for '{input_path}' does not exist.")
	return input_path


def input_path_exists(input_path, function, exception_str):
	# Convert into list if it is a single path/str
	path_list = input_path if type(input_path) is list else [input_path]
	for path in path_list:
		if not function(path):
			sys.exit(exception_str.format(path))
	return input_path


def create_output_directory(output_path):
	""""""
	created, output_path = fo.create_directory(output_path)
	if not created:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)

	return output_path


def check_ptf_tref_conflict(training_file, training_reference):
	""""""
	if training_file and training_reference:
		sys.exit(ct.CANNOT_PROVIDE_PTF_AND_TREF)


def schema_includes_fasta(schema_directory):
	""""""
	schema_files = fo.listdir(schema_directory)
	# Check if the folder includes FASTA files
	fasta_extension, _ = fo.filter_by_extension(schema_files, ct.FASTA_EXTENSIONS)
	if len(fasta_extension) == 0:
		sys.exit(ct.MISSING_SCHEMA_FASTAS)

	return schema_directory


def schema_check_genes_list(schema_directory):
	""""""
	schema_files = fo.listdir(schema_directory)
	# Check if schema includes the .genes_list file
	if ct.GENE_LIST_BASENAME not in schema_files:
		sys.exit(ct.MISSING_GENE_LIST)

	# Verify that loci FASTA files match list of loci in the .genes_list file
	schema_files_fasta = fo.filter_by_extension(schema_files, [".fasta"])[0]
	schema_files_basenames = [fo.file_basename(file, False) for file in schema_files_fasta]
	loci_list = fo.pickle_loader(fo.join_paths(schema_directory, [ct.GENE_LIST_BASENAME]))

	# Determine if the schema or the list of loci have different sets of loci IDs
	schema_exclusive = set(schema_files_basenames) - set(loci_list)
	loci_list_exclusive = set(loci_list) - set(schema_files_basenames)
	if len(schema_exclusive) > 0 or len(loci_list_exclusive) > 0:
		sys.exit(ct.SCHEMA_LOCI_GENE_LIST_NOMATCH.format(len(schema_exclusive), im.join_list(schema_exclusive, ", "), len(loci_list_exclusive), im.join_list(loci_list_exclusive, ", ")))

	return schema_directory


def schema_check_config(schema_directory):
	""""""
	schema_files = fo.listdir(schema_directory)
	# Check if schema includes the .schema_config file
	if ct.SCHEMA_CONFIG_BASENAME not in schema_files:
		sys.exit(ct.SCHEMA_MISSING_FILES.format(ct.SCHEMA_CONFIG_BASENAME))

	return schema_directory


def schema_check_short(schema_directory):
	""""""
	schema_files = fo.listdir(schema_directory)
	# Check if the "short" directory exists
	if "short" not in schema_files:
		sys.exit(ct.SCHEMA_MISSING_FILES.format("short"))

	return schema_directory


def check_argument_conflict(arg_name, arg_value, schema_config, force_continue):
	""""""
	if arg_value not in schema_config[arg_name]:
		print(ct.CONFLICTING_ARGUMENTS.format(arg_name, schema_config[arg_name]))
		if not force_continue:
			answer = fo.input_timeout(ct.ARGS_DIFFER_PROMPT, ct.PROMPT_TIMEOUT)
		else:
			answer = 'yes'

		if answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			schema_config[arg_name].append(arg_value)

	return arg_value, schema_config


def contains_results(results_directory):
	""""""
	# List files in input direcoty
	results_files = fo.listdir(results_directory)
	# Check if folder includes the results_alleles.tsv file with allelic profiles
	if ct.RESULTS_ALLELES_BASENAME not in results_files:
		sys.exit(ct.MISSING_RESULTS_ALLELES)

	return results_directory


def ensure_list(input_val):
	"""
	"""
	if isinstance(input_val, str):
		return [input_val]
	return input_val


# Define reusable fields
OutputDirectory = Annotated[str,
							AfterValidator(parentdir_exists),
							AfterValidator(create_output_directory)]

InputFiles = Annotated[str,
					   AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_INPUT_ARG)),
					   AfterValidator(input_is_fasta),
					   AfterValidator(list_input_files),
					   AfterValidator(filter_inputs_extension),
					   AfterValidator(inputs_exist),
					   AfterValidator(validate_inputs_fastas),
					   AfterValidator(check_unique_prefixes),
					   AfterValidator(check_blanks)]

ChewieSchema = Annotated[str,
						 AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.SCHEMA_PATH_MISSING)),
						 AfterValidator(schema_includes_fasta),
						 AfterValidator(schema_check_genes_list),
						 AfterValidator(schema_check_config),
						 AfterValidator(schema_check_short)]

ExternalSchema = Annotated[str,
							AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.SCHEMA_PATH_MISSING)),
							AfterValidator(schema_includes_fasta)]

LociList = Annotated[str,
					 AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_LOCI_LIST)),
					 AfterValidator(input_is_fasta),
					 AfterValidator(list_input_files)]

BLASTScoreRatio = Annotated[float,
							Field(default=ct.BSR_DEFAULT,
							      ge=ct.BSR_MIN,
								  le=ct.BSR_MAX)]

MinimumLength = Annotated[int,
						  Field(default=ct.MSL_DEFAULT,
						        ge=ct.MSL_MIN,
								le=ct.MSL_MAX)]

TranslationTable = Annotated[int,
							 Field(default=ct.GENETIC_CODE_DEFAULT),
							 AfterValidator(partial(validate_choice, valid_options=ct.GENETIC_CODES, exception_str=ct.INVALID_GENETIC_CODE))]

SizeThreshold = Annotated[float | None,
						  Field(default=ct.ST_DEFAULT,
						        ge=ct.ST_MIN,
								le=ct.ST_MAX)]

GenePredictor = Annotated[str,
						  Field(default=ct.GENE_PREDICTOR_DEFAULT),
						  AfterValidator(partial(validate_choice, valid_options=ct.GENE_PREDICTORS, exception_str=ct.INVALID_GENE_PREDICTOR))]

# Validate if input is str, not if it is None
GenePredictionOptions = Annotated[list,
                                  Field(default={}),
                                  AfterValidator(parse_parameter_string)]

ClusteringOptions = Annotated[str,
                              Field(default={}),
                              AfterValidator(parse_parameter_string)]

BLASTPath = Annotated[str | None, # Need to specify that the argument value may not be provided
                      Field(default=None, # Need to set default to None to really set arg value as None
							validate_default=True), # Needs to be validated even if it is None
					  AfterValidator(check_blast),
					  AfterValidator(check_blast_version)]

CPUCores = Annotated[int,
					 Field(default=ct.CPU_CORES_DEFAULT,
					       ge=ct.CPU_CORES_MIN,
						   le=multiprocessing.cpu_count()),
					 AfterValidator(verify_cpu_usage)]

PyrodigalTrainingFile = Annotated[str,
                                  Field(default=None),
                                  AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.INVALID_PTF_PATH))]

PyrodigalMode = Annotated[str,
						  Field(default=ct.PYRODIGAL_DEFAULT_MODE),
						  AfterValidator(partial(validate_choice, valid_options=ct.PYRODIGAL_MODES, exception_str=ct.INVALID_PYRODIGAL_MODE))]

PyrodigalOutputFormats = Annotated[list,
								   Field(default=ct.PYRODIGAL_DEFAULT_OUTFMTS),
								   BeforeValidator(ensure_list),
								   AfterValidator(partial(validate_multiple_choices, valid_options=ct.PYRODIGAL_OUTFMTS, exception_str=ct.PYRODIGAL_INVALID_OUTFMT))]

PyrodigalMinimumConfidence = Annotated[float,
                                       Field(default=None,
                                             ge=ct.PYRODIGAL_MIN_CONFIDENCE,
											 le=ct.PYRODIGAL_MAX_CONFIDENCE)]

PyrodigalTrainingReference = Annotated[str,
                                       Field(default=None),
                                       AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.INVALID_PTF_PATH))]

PyrodigalJustTraining = Annotated[bool,
                                  Field(default=False)]

AugustusPath = Annotated[str | None,
                         Field(default=None,
                               validate_default=True),
						 AfterValidator(check_augustus)]

AugustusSpecies = Annotated[str,
                            Field(default=None)]

AugustusOutputFormat = Annotated[list,
								 Field(default=ct.AUGUSTUS_OUTFMT_DEFAULT),
								 BeforeValidator(ensure_list),
								 AfterValidator(partial(validate_multiple_choices, valid_options=ct.AUGUSTUS_OUTFMTS, exception_str=ct.AUGUSTUS_INVALID_OUTFMT))]

WordSize = Annotated[int,
					 Field(default=ct.WORD_SIZE_DEFAULT,
					       ge=ct.WORD_SIZE_MIN,
						   le=ct.WORD_SIZE_MAX)]

WindowSize = Annotated[int,
					   Field(default=ct.WORD_SIZE_DEFAULT,
					         ge=ct.WINDOW_SIZE_MIN,
							 le=ct.WINDOW_SIZE_MAX)]

ClusteringSimilarity = Annotated[float,
								 Field(default=ct.CLUSTERING_SIMILARITY_DEFAULT,
								       ge=ct.CLUSTERING_SIMILARITY_MIN,
									   le=ct.CLUSTERING_SIMILARITY_MAX)]

RepresentativeFilter = Annotated[float,
								 Field(default=ct.REPRESENTATIVE_FILTER_DEFAULT,
								       ge=ct.REPRESENTATIVE_FILTER_MIN,
									   le=ct.REPRESENTATIVE_FILTER_MAX)]

IntraFilter = Annotated[float,
						Field(default=(ct.INTRA_CLUSTER_DEFAULT),
						      ge=ct.INTRA_CLUSTER_MIN,
							  le=ct.INTRA_CLUSTER_MAX)]

Annotations = Annotated[str,
                        AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_ANNOTATIONS_FILE))]


BooleanFlag = Annotated[bool,
                        Field(default=False)]


# Use separate models to validate Pyrodigal and AUGUSTUS arguments, since they have different fields and validation requirements
class PyrodigalArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	pyrodigal_training_file: PyrodigalTrainingFile
	pyrodigal_mode: PyrodigalMode
	pyrodigal_output_formats: PyrodigalOutputFormats
	pyrodigal_minimum_confidence: PyrodigalMinimumConfidence
	pyrodigal_training_reference: PyrodigalTrainingReference
	pyrodigal_just_training: PyrodigalJustTraining


class AugustusArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	augustus_path: AugustusPath
	augustus_species: AugustusSpecies
	augustus_output_format: AugustusOutputFormat

	@model_validator(mode="after")
	def validate_species(self):
		self.augustus_species = validate_augustus_species(self.augustus_species, self.augustus_path)

		return self


# Field definition for the validated gene prediction arguments
ValidatedGenePredictionOptions = Annotated[PyrodigalArgs | AugustusArgs | None, Field(default=None)]


# Use a separate model to validate clustering arguments passed to modules that use clustering
class ClusteringArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	word_size: WordSize
	window_size: WindowSize
	clustering_similarity: ClusteringSimilarity
	representative_filter: RepresentativeFilter
	intra_filter: IntraFilter


ValidatedClusteringOptions = Annotated[ClusteringArgs | None, Field(default=None)]


class PredictGenesValidator(BaseModel):
	# Use `from_attributes=True` to inspect the attributes of the 
	# argparse.Namespace object directly without having to convert to a dictionary
	# Do not use `extra="forbid"` so that it accepts the arg with the module name without having to remove it explicitly
	model_config = ConfigDict(from_attributes=True)

	# Field order is preserved when serializing with model_dump()
	output_directory: OutputDirectory
	input_files: InputFiles
	gene_predictor: GenePredictor
	gene_prediction_options: GenePredictionOptions
	validated_gene_prediction_options: ValidatedGenePredictionOptions
	translation_table: TranslationTable
	cpu_cores: CPUCores

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options = PyrodigalArgs(**self.gene_prediction_options)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_options = AugustusArgs(**self.gene_prediction_options)

		return self

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_file is not None:
			self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_options.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options.pyrodigal_mode = check_meta(self.validated_gene_prediction_options.pyrodigal_mode, self.validated_gene_prediction_options.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def ptf_and_reference(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_reference:
			# Check if user provided path to Pyrodigal training file and to a training reference, which is not allowed
			check_ptf_tref_conflict(self.validated_gene_prediction_options.pyrodigal_training_file, self.validated_gene_prediction_options.pyrodigal_training_reference)
			# Create training file based on training reference
			print(f'Creating Pyrodigal training file based on {self.validated_gene_prediction_options.pyrodigal_training_reference}...')
			self.validated_gene_prediction_options.pyrodigal_training_file = pgp.create_training_file(self.validated_gene_prediction_options.pyrodigal_training_reference, self.output_directory, self.translation_table)
			print(f'Training file saved to {self.validated_gene_prediction_options.pyrodigal_training_file}')
			# Exit if user only requested to create training file
			if self.validated_gene_prediction_options.pyrodigal_just_training:
				sys.exit(ct.JUST_TRAINING)

		return self


SchemaName = Annotated[str,
                       Field(default=ct.SCHEMA_NAME_DEFAULT)]


class CreateSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_files: InputFiles
	schema_name: SchemaName
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	gene_predictor : GenePredictor
	gene_prediction_options: GenePredictionOptions
	validated_gene_prediction_options: ValidatedGenePredictionOptions
	clustering_options: ClusteringOptions
	validated_clustering_options: ValidatedClusteringOptions
	blast_path : BLASTPath
	cds_input: BooleanFlag
	no_cds_renaming: BooleanFlag
	cpu_cores: CPUCores
	no_cleanup: BooleanFlag

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options = PyrodigalArgs(**self.gene_prediction_options)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_options = AugustusArgs(**self.gene_prediction_options)

		return self

	# Further validation for clustering arguments
	@model_validator(mode="after")
	def validate_clustering_arguments(self):
		self.validated_clustering_options = ClusteringArgs(**self.clustering_options)

		return self

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_file is not None:
			self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_options.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options.pyrodigal_mode = check_meta(self.validated_gene_prediction_options.pyrodigal_mode, self.validated_gene_prediction_options.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def ptf_and_reference(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_reference:
			# Check if user provided path to Pyrodigal training file and to a training reference, which is not allowed
			check_ptf_tref_conflict(self.validated_gene_prediction_options.pyrodigal_training_file, self.validated_gene_prediction_options.pyrodigal_training_reference)
			# Create training file based on training reference
			print(f'Creating Pyrodigal training file based on {self.validated_gene_prediction_options.pyrodigal_training_reference}...')
			self.validated_gene_prediction_options.pyrodigal_training_file = pgp.create_training_file(self.validated_gene_prediction_options.pyrodigal_training_file, self.validated_gene_prediction_options.pyrodigal_training_reference, self.output_directory, self.translation_table)
			print(f'Training file saved to {self.validated_gene_prediction_options.pyrodigal_training_file}')
			# Exit if user only requested to create training file
			if self.validated_gene_prediction_options.pyrodigal_just_training:
				sys.exit(ct.JUST_TRAINING)

		return self


AllelecallMode = Annotated[int,
                           Field(default=ct.ALLELECALL_DEFAULT_MODE),
						   AfterValidator(partial(validate_choice, valid_options=ct.ALLELECALL_MODES, exception_str=ct.INVALID_ALLELECALL_MODE))]

NSCONFIG = Annotated[bool,
                     Field(default=False)]


class AlleleCallValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_files: InputFiles
	schema_directory: ChewieSchema
	# Do not run validation if loci_list is None
	loci_list: LociList | None = None
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	gene_predictor : GenePredictor
	gene_prediction_options: GenePredictionOptions
	validated_gene_prediction_options: ValidatedGenePredictionOptions
	clustering_options: ClusteringOptions
	validated_clustering_options: ValidatedClusteringOptions
	blast_path : BLASTPath
	cds_input: BooleanFlag
	no_inferred: BooleanFlag
	output_unclassified: BooleanFlag
	output_missing: BooleanFlag
	output_novel: BooleanFlag
	output_masked_profiles: BooleanFlag
	no_cds_renaming: BooleanFlag
	force_continue: BooleanFlag
	mode: AllelecallMode
	cpu_cores: CPUCores
	no_cleanup: BooleanFlag
	ns_config: NSCONFIG

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options = PyrodigalArgs(**self.gene_prediction_options)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_options = AugustusArgs(**self.gene_prediction_options)

		return self

	# Further validation for clustering arguments
	@model_validator(mode="after")
	def validate_clustering_arguments(self):
		self.validated_clustering_options = ClusteringArgs(**self.clustering_options)

		return self

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_ids = [fo.file_basename(file, False) for file in loci_files]
			extra_loci = set(self.loci_list) - set(loci_ids)
			if len(extra_loci) > 0:
				sys.exit(ct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(im.join_list(extra_loci, ", ")))
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# List all loci FASTA files in the schema and add it to the loci_list arg
	@model_validator(mode="after")
	def list_schema_loci(self):
		schema_files = fo.listdir(self.schema_directory, full_path=True)
		loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
		loci_ids = [fo.file_basename(file, False) for file in loci_files]
		self.loci_list = [loci_ids, loci_files]

		return self

	@model_validator(mode="after")
	def solve_conflicting_arguments(self):
		config_file = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
		config = fo.pickle_loader(config_file)
		# BLAST Score Ratio
		self.blast_score_ratio, config = check_argument_conflict(ct.BLAST_SCORE_RATIO_ARGNAME, self.blast_score_ratio, config, self.force_continue)
		# Minimum length
		self.minimum_length, config = check_argument_conflict(ct.MINIMUM_LENGTH_ARGNAME, self.minimum_length, config, self.force_continue)
		# Size threshold
		self.size_threshold, config = check_argument_conflict(ct.SIZE_THRESHOLD_ARGNAME, self.size_threshold, config, self.force_continue)
		# Pyrodigal training file
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_file is not None:
			self.validated_gene_prediction_options.pyrodigal_training_file, config = check_argument_conflict(ct.PYRODIGAL_TRAININGFILE_ARGNAME, self.validated_gene_prediction_options.pyrodigal_training_file, config, self.force_continue)
			# Compare current value for the translation table and the one in the training file
			self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_options.pyrodigal_training_file, self.translation_table)
			# Add translation table to config if it is different than the ones previously used
			if self.translation_table not in config[ct.TRANSLATION_TABLE_ARGNAME]:
				config[ct.TRANSLATION_TABLE_ARGNAME].append(self.translation_table)

		# Update schema config file
		fo.pickle_dumper(config, config_file)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options.pyrodigal_mode = check_meta(self.validated_gene_prediction_options.pyrodigal_mode, self.validated_gene_prediction_options.pyrodigal_training_file)

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [ct.NS_CONFIG_BASENAME])
		if fo.is_file(ns_config)[0]:
			self.ns_config = True

		return self


class SchemaEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: ChewieSchema
	loci_list: LociList | None = None
	annotations_file: Annotations | None = None
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	minimum_length: MinimumLength
	cpu_cores: CPUCores
	loci_reports: BooleanFlag
	light: BooleanFlag
	add_sequences: BooleanFlag
	chewie_schema: BooleanFlag
	minimum_length_values: Annotated[list, Field(default=[])]
	size_threshold_values: Annotated[list, Field(default=[])]
	translation_table_values: Annotated[list, Field(default=[])]
	blast_score_ratio_values: Annotated[list, Field(default=[])]
	training_file: BooleanFlag
	ns_config: NSCONFIG

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_ids = [fo.file_basename(file, False) for file in loci_files]
			extra_loci = set(self.loci_list) - set(loci_ids)
			if len(extra_loci) > 0:
				sys.exit(ct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(im.join_list(extra_loci, ", ")))
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# List all loci FASTA files in the schema and add it to the loci_list arg
	@model_validator(mode="after")
	def list_schema_loci(self):
		if not self.loci_list:
			schema_files = fo.listdir(self.schema_directory, full_path=True)
			loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_ids = [fo.file_basename(file, False) for file in loci_files]
			self.loci_list = [loci_ids, loci_files]

		return self

	@model_validator(mode="after")
	def set_schema_origin(self):
		schema_files = fo.listdir(self.schema_directory, full_path=False)
		# Schema was created with chewBBACA
		if ct.SCHEMA_CONFIG_BASENAME in schema_files:
			config_file  = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Set the argument that indicates that the schema was created with chewBBACA to True
			self.chewie_schema = config_values.get(ct.CHEWIE_VERSION_ARGNAME)[0]
		else:
			self.chewie_schema = False

		return self

	# Use the values stored in the schema config if it is a chewBBACA schema
	@model_validator(mode="after")
	def force_schema_args(self):
		if self.chewie_schema:
			print("Schema was created with chewBBACA. Using the values stored in the schema "
			      "config file for the minimum length, size threshold, and translation table.")
			config_file  = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Get the first minimum_length, size_threshold, and translation_table values used with the schema
			self.minimum_length = config_values.get(ct.MINIMUM_LENGTH_ARGNAME)[0]
			self.size_threshold = config_values.get(ct.SIZE_THRESHOLD_ARGNAME)[0]
			self.translation_table = config_values.get(ct.TRANSLATION_TABLE_ARGNAME)[0]

		return self

	@model_validator(mode="after")
	def set_schema_args(self):
		if self.chewie_schema:
			config_file  = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Get the first minimum_length, size_threshold, and translation_table values used with the schema
			self.minimum_length_values = config_values.get(ct.MINIMUM_LENGTH_ARGNAME)
			self.size_threshold_values = config_values.get(ct.SIZE_THRESHOLD_ARGNAME)
			self.translation_table_values = config_values.get(ct.TRANSLATION_TABLE_ARGNAME)
			self.blast_score_ratio_values = config_values.get(ct.BLAST_SCORE_RATIO_ARGNAME)
			self.training_file = True if config_values.get(ct.PYRODIGAL_TRAININGFILE_ARGNAME)[0] else False

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [ct.NS_CONFIG_BASENAME])
		if fo.is_file(ns_config)[0]:
			self.ns_config = True

		return self


# Nedd to add a AfterValidator to get and validate results files
ResultsFiles = Annotated[str,
                         AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_RESULTS)),
						 AfterValidator(contains_results)]


class AlleleCallEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	translation_table: TranslationTable
	output_directory: OutputDirectory
	schema_directory: ChewieSchema
	results_files: ResultsFiles
	annotations_file: Annotations | None = None
	cpu_cores: CPUCores
	light: BooleanFlag
	no_pa: BooleanFlag
	no_dm: BooleanFlag
	no_tree: BooleanFlag
	cg_alignment: BooleanFlag

	@model_validator(mode="after")
	def get_translation_table(self):
		config_file = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
		if fo.exists(config_file):
			config_values = fo.pickle_loader(config_file)
			self.translation_table = config_values[ct.TRANSLATION_TABLE_ARGNAME][0]

		return self


Profiles = Annotated[str,
                     AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_EXCLUDE_LOCI))]

LociPresenceThreshold = Annotated[list,
                            Field(default=ct.CGMLST_THRESHOLDS,
                                  ge=ct.CGMLST_THRESHOLD_MIN,
								  le=ct.CGMLST_THRESHOLD_MAX)]

SampleStep = Annotated[int,
                 Field(default=1)]

Permutation = Annotated[int,
                        Field(default=ct.PERMUTATION_NUMBER_DEFAULT)]

ExcludeLoci = Annotated[str,
                        AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_EXCLUDE_LOCI)),
						AfterValidator(input_is_fasta),
						AfterValidator(list_input_files)]

ExcludeGenomes = Annotated[str,
                           AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_EXCLUDE_GENOMES)),
						   AfterValidator(input_is_fasta),
						   AfterValidator(list_input_files)]


class ExtractCgMLSTValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	# Need to validate the profiles (header starting with FILE followed by loci IDs, header line followed by profiles, etc?)
	profiles: Profiles
	loci_presence_threshold: LociPresenceThreshold
	sample_step: SampleStep
	compute_accessory: BooleanFlag
	rarefaction_analysis: BooleanFlag
	permutation_number: Permutation
	permutation_samples: int | None = None
	cpu_cores: CPUCores
	exclude_loci: ExcludeLoci | None = None
	exclude_genomes: ExcludeGenomes | None = None

	@model_validator(mode="after")
	def check_exclude_loci(self):
		if self.exclude_loci:
			# Read list of loci in the profiles file
			results_loci = fo.read_lines(self.profiles, strip=True, num_lines=1)
			results_loci = results_loci[0].split("\t")[1:]
			extra_loci = set(self.exclude_loci) - set(results_loci)
			if len(extra_loci) > 0:
				print(ct.INVALID_EXCLUDE_LOCI.format("\n".join(extra_loci)))

		return self

	@model_validator(mode="after")
	def check_exclude_genomes(self):
		if self.exclude_genomes:
			# Read list of samples in the profiles file
			results_samples = fo.extract_column(self.profiles)
			extra_samples = set(self.exclude_genomes) - set(results_samples)
			if len(extra_samples) > 0:
				print(ct.INVALID_EXCLUDE_SAMPLES.format("\n".join(extra_samples)))

		return self


SampleList = Annotated[str,
                       AfterValidator(partial(input_path_exists, function=fo.exists, exception_str=ct.MISSING_EXCLUDE_GENOMES)),
					   AfterValidator(input_is_fasta),
					   AfterValidator(list_input_files)]


class SubsetResultsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_files: ResultsFiles
	loci_list: LociList | None = None
	sample_list: SampleList | None = None
	inverse_loci: BooleanFlag
	inverse_samples: BooleanFlag

	@model_validator(mode="after")
	def check_lists(self):
		# Do not proceed if user provided neither a list of loci nor of samples
		if not self.loci_list and not self.sample_list:
			sys.exit(ct.SUBSETRESULTS_MISSING_LISTS_EXCEPTION)

		return self

	@model_validator(mode="after")
	def check_files(self):
		# Do not proceed if the input folder does not include any of the supported output files
		input_files = fo.listdir(self.results_files)
		if len(set(ct.ALLELECALL_OUTFILES)-set(input_files)) == len(ct.ALLELECALL_OUTFILES):
			sys.exit("Input directory does not include any of the supported output files that can be subsetted."
			         "Please ensure that the input directory includes at least one of the following files:\n"
					 f"{'\n'.join(ct.ALLELECALL_OUTFILES)}")

	# Check if the loci IDs in the input list are all in the results
	@model_validator(mode="after")
	def check_loci(self):
		file_types = [ct.RESULTS_ALLELES_BASENAME, ct.RESULTS_COORDINATES_BASENAME, ct.PRESENCE_ABSENCE_BASENAME, ct.LOCI_STATS_BASENAME]
		if self.loci_list:
			for fname in file_types:
				# Determine if file exists in the input directory
				file_path = fo.join_paths(self.results_files, [fname])
				if fo.is_file(file_path)[0]:
					if fname in file_types[:3]:
						results_loci = fo.read_lines(file_path, strip=True, num_lines=1)
						results_loci = results_loci[0].split("\t")[1:]
					elif fname in file_types[3:]:
						results_loci = pd.read_csv(file_path, usecols=["Locus"], sep="\t", dtype=str)
						results_loci = results_loci["Locus"].tolist()

					non_matching = set(self.loci_list) - set(results_loci)
					if len(non_matching) > 0:
						sys.exit(f"The following loci IDs in the input loci list do not match "
							     f"loci IDs in the results files:\n{'\n'.join(non_matching)}\n"
								 "Please ensure that all loci IDs in the input loci list match "
								 "loci IDs in the resutls files.")
			# It is possible that the input directory does not include any of the files containing the loci IDs
			if not fo.is_file(file_path)[0]:
				print(f"Could not get the list of loci IDs from the results files in the input directory "
				      "to validate the loci IDs in the input loci list.")

		return self

	# Check if the sample IDs in the input list are all in the results
	@model_validator(mode="after")
	def check_samples(self):
		file_types = [ct.RESULTS_ALLELES_BASENAME, ct.RESULTS_COORDINATES_BASENAME, ct.PRESENCE_ABSENCE_BASENAME, ct.LOCI_STATS_BASENAME]
		if self.sample_list:
			for fname in file_types:
				file_path = fo.join_paths(self.results_files, [fname])
				if fo.is_file(file_path)[0]:
					results_samples = pd.read_csv(file_path, usecols=["FILE"], sep="\t", dtype=str)
					results_samples = results_samples["FILE"].tolist()
					non_matching = set(self.sample_list) - set(results_samples)
					if len(non_matching) > 0:
						sys.exit(f"The following sample IDs in the input sample list do not match "
							  	 f"sample IDs in the results files:\n{'\n'.join(non_matching)}\n"
							  	 "Please ensure that all sample IDs in the input sample list match "
							  	 "sample IDs in the resutls files.")
			if not fo.is_file(file_path)[0]:
				print(f"Could not get the list of sample IDs from the results files in the input directory "
				      "to validate the sample IDs in the input loci list.")

		return self


ResultsDirectories = Annotated[list,
								AfterValidator(partial(input_path_exists,
								                       function=fo.is_directory,
													   exception_str=ct.MERGERESULTS_INPUTFILE_EXCEPTION))]


class MergeResults(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_directories: ResultsDirectories
	common: BooleanFlag

	# Determine if user passed more than one input path
	@model_validator(mode="after")
	def check_multiple_input_folders(self):
		if len(self.results_directories) <= 1:
			sys.exit("Provided a single input path. Need multiple input paths to merge results files from different folders.")

		return self

	# Check if results have loci in common
	@model_validator(mode="after")
	def validate_loci(self):
		file_types = [ct.RESULTS_ALLELES_BASENAME, ct.RESULTS_COORDINATES_BASENAME, ct.PRESENCE_ABSENCE_BASENAME, ct.LOCI_STATS_BASENAME]
		results_files = {fname: [] for fname in file_types}
		for folder in self.results_directories:
			# List files in folder
			current_files = fo.listdir(folder)
			for fname, files in results_files.items():
				if fname in current_files:
					results_files[fname].append(fo.join_paths(folder, [fname]))

		for fname, files in results_files.items():
			if len(files) != len(self.results_directories):
				print(f"Cannot determine the set of shared loci based on the {fname} "
				      "file type because it is not present in all input directories.")
			# All input directories contain the file type
			# Can use the file type to check if results have common loci
			else:
				headers = []
				for f in files:
					if fname in file_types[:3]:
						header = fo.read_lines(f, strip=True, num_lines=1)
						results_loci = header[0].split("\t")[1:]
					elif fname in file_types[3:]:
						results_loci = pd.read_csv(f, usecols=["Locus"], sep="\t", dtype=str)
						results_loci = results_loci["Locus"].tolist()
					headers.append(results_loci)
	
				loci = set(headers[0])
				for h in headers[1:]:
					# Determine common loci ordered based on headers in first file
					loci = loci.intersection(set(h))

				# Sort list of common loci
				loci = sorted(list(loci))

				# Exit if there are no common loci between all results
				if len(loci) == 0:
					sys.exit(ct.MERGERESULTS_NOCOMMONLOCI_EXCEPTION)
				# Exit if results are not for the same set of loci
				elif len(loci) != len(headers[0]) and not self.common:
					sys.exit(ct.MERGERESULTS_LOCIDIFFER_EXCEPTION)

				self.common = [self.common, loci, len(headers[0])]

				# No need to check based on other file type
				break

			# Exit if it was not possible to determine if results have loci IDs in common
			if type(self.common) == bool:
				sys.exit("Could not determine if the results files included in the input "
				         "directories share loci IDs. Will not attempt to merge the results."
						 "Please make sure that all input directories include at least one "
						 f"of the following file types: {'\n'.join(file_types)}")

		return self


# Should implement a validator for the profile filr format
AllelicProfiles = Annotated[str,
                            AfterValidator(partial(input_path_exists,
                                                   function=fo.exists,
												   exception_str=ct.PROFILES_FILE_MISSING))]


def validate_hashing_algo(hashing_algo):
	""""""
	# Get hash function
	hashing_function = getattr(hashlib, hashing_algo, None)
	if hashing_function is None:
		hashing_function = getattr(zlib, hashing_algo, None)

	if hashing_function is None:
		sys.exit(ct.HASHPROFILES_INVALID_HASHING.format(hashing_algo))

	return [hashing_algo, hashing_function]


HashingAlgo = Annotated[str,
                        Field(default=ct.HASHPROFILES_DEFAULT_HASH_ALGO,
                              validate_default=True),
						AfterValidator(validate_hashing_algo)]

NrowChunk = Annotated[int,
                      Field(default=100)]


class HashProfilesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: AllelicProfiles
	schema_directory: ChewieSchema
	hashing_algo: HashingAlgo
	nrow_chunk: NrowChunk
	cpu_cores: CPUCores

	# Determine if all loci IDs in the profiles are in the schema
	@model_validator(mode="after")
	def profiles_loci_in_schema(self):
		# Read loci IDs in profiles file
		profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
		profiles_loci = profiles_header[0].split("\t")[1:]
		# Get list of loci in the schema
		schema_files = fo.listdir(self.schema_directory)
		schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
		schema_loci = [fo.file_basename(file, False) for file in schema_loci]

		extra_loci = set(profiles_loci) - set(schema_loci)
		if len(extra_loci) > 0:
			sys.exit("The following loci IDs in the input profiles are not in the schema:\n"
					 f"{'\n'.join(extra_loci)}\nPlease ensure that all loci IDs in the input profiles "
					 "are in the schema.")

		return self


class GetAllelesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: AllelicProfiles
	schema_directory: ChewieSchema
	loci_list: LociList | None = None
	cpu_cores: CPUCores
	translate_alleles: BooleanFlag
	translation_table: TranslationTable
	distinct: BooleanFlag

	# Use translation table defined in training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		# Determine if schema includes a training file
		schema_files = fo.listdir(self.schema_directory, full_path=True)
		training_file = fo.filter_by_extension(schema_files, [".trn"])[0]
		if len(training_file) > 0:
			self.translation_table = add_ptf_genetic_code(training_file[0], self.translation_table)

		return self

	# Determine if all loci IDs in the profiles are in the schema
	@model_validator(mode="after")
	def profiles_loci_in_schema(self):
		# Read loci IDs in profiles file
		profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
		profiles_loci = profiles_header[0].split("\t")[1:]
		# Get list of loci in the schema
		schema_files = fo.listdir(self.schema_directory)
		schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
		schema_loci = [fo.file_basename(file, False) for file in schema_loci]

		extra_loci = set(profiles_loci) - set(schema_loci)
		if len(extra_loci) > 0:
			sys.exit("The following loci IDs in the input profiles are not in the schema:\n"
				     f"{'\n'.join(extra_loci)}\nPlease ensure that all loci IDs in the input profiles "
					 "are in the schema.")

		return self

	# Verify that the profiles include the loci IDs in the list of loci
	@model_validator(mode="after")
	def loci_list_in_profiles(self):
		if self.loci_list:
			# Read loci IDs in profiles file
			profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
			profiles_loci = profiles_header[0].split("\t")[1:]
			extra_loci = set(self.loci_list)-set(profiles_loci)
			if len(extra_loci) > 0:
				sys.exit("The following loci IDs in the input loci list are not in the allelic profiles:\n"
				         f"{'\n'.join(extra_loci)}\nPlease ensure that are all loci IDs in the input loci "
						 "list are in the allelic profiles.")

		return self

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
			schema_loci = [fo.file_basename(file, False) for file in schema_loci]
			extra_loci = set(self.loci_list) - set(schema_loci)
			if len(extra_loci) > 0:
				sys.exit(ct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(im.join_list(extra_loci, ", ")))
			# Get full paths
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# Simply get the list of loci in the profiles if no loci list is provided
	@model_validator(mode="after")
	def get_profiles_loci(self):
		if not self.loci_list:
			# Read loci IDs in profiles file
			profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
			profiles_loci = profiles_header[0].split("\t")[1:]
			# Get full paths
			loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in profiles_loci]
			self.loci_list = [profiles_loci, loci_fpaths]

		return self


class PrepExternalSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: ExternalSchema
	loci_list: LociList | None = None
	gene_predictor: GenePredictor
	gene_prediction_options: GenePredictionOptions
	validated_gene_prediction_options: ValidatedGenePredictionOptions
	clustering_options: ClusteringOptions
	validated_clustering_options: ValidatedClusteringOptions
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	adaptation_minimum_length: Annotated[int, Field(default=0)]
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	adaptation_size_threshold: Annotated[float, Field(default=None)]
	cpu_cores: CPUCores
	blast_path : BLASTPath
	size_filter: BooleanFlag

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options = PyrodigalArgs(**self.gene_prediction_options)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_options = AugustusArgs(**self.gene_prediction_options)

		return self

	# Further validation for clustering arguments
	@model_validator(mode="after")
	def validate_clustering_arguments(self):
		self.validated_clustering_options = ClusteringArgs(**self.clustering_options)

		return self

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
			schema_loci = [fo.file_basename(file, False) for file in schema_loci]
			absent_files = set(self.loci_list) - set(schema_loci)
			if len(absent_files) > 0:
				sys.exit(ct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(im.join_list(absent_files, ", ")))
			# Get full paths
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	@model_validator(mode="after")
	def get_schema_loci(self):
		if not self.loci_list:
			# Get full paths
			schema_files = fo.listdir(self.schema_directory, full_path=True)
			loci_fpaths = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_basenames = [fo.file_basename(file, False) for file in loci_fpaths]
			self.loci_list = [loci_basenames, loci_fpaths]

		return self

	@model_validator(mode="after")
	def apply_size_filter(self):
		if self.size_filter:
			self.adaptation_minimum_length = self.minimum_length
			print(f"Using a minimum length value of {self.adaptation_minimum_length} for schema "
				  f"adaptation and {self.minimum_length} to store in the schema "
				  "config file.")
			self.adaptation_size_threshold = self.size_threshold
			print(f"Using a size threshold value of {self.adaptation_size_threshold} for schema "
		    	  f"adaptation and {self.size_threshold} to store in the schema "
		    	  "config file.")

		return self


ProteinTable = Annotated[str,
                         Field(default=None),
                         AfterValidator(partial(input_path_exists,
                                                function=fo.exists,
												exception_str=ct.MISSING_PROTEIN_TABLE))]

Taxa = Annotated[list,
                 Field(default=None)]

ProteomeMatches = Annotated[int,
                            Field(default=1, ge=1)]

RestApiMatches = Annotated[int,
                           Field(default=1, ge=1)]


class AnnotateSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: ExternalSchema
	loci_list: LociList | None = None
	protein_table: ProteinTable
	blast_score_ratio: BLASTScoreRatio
	translation_table: TranslationTable
	cpu_cores: CPUCores
	proteome_taxa: Taxa
	proteome_matches: ProteomeMatches
	blast_path: BLASTPath
	rest_api: BooleanFlag
	rest_api_matches: RestApiMatches
	no_cleanup: BooleanFlag

	# User must specify that loci should be annotated based on either proteomes or requests to UniProt's Proteins API or both
	@model_validator(mode="after")
	def check_annotation_method(self):
		if not self.proteome_taxa and not self.rest_api:
			sys.exit(ct.ANNOTATE_SCHEMA_NO_ANNOTATION_METHOD)

		return self

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_ids = [fo.file_basename(file, False) for file in loci_files]
			extra_loci = set(self.loci_list) - set(loci_ids)
			if len(extra_loci) > 0:
				sys.exit(ct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(im.join_list(extra_loci, ", ")))
			# Get full paths
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# Annotate all loci if user did not provide a list of loci
	@model_validator(mode="after")
	def get_schema_loci(self):
		if not self.loci_list:
			# Get full paths
			schema_files = fo.listdir(self.schema_directory, full_path=True)
			loci_fpaths = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_basenames = [fo.file_basename(file, False) for file in loci_fpaths]
			self.loci_list = [loci_basenames, loci_fpaths]

		return self

	# Use translation table stored in schema config file if it is a chewie schema
	@model_validator(mode="after")
	def get_translation_table(self):
		config_file = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
		if fo.exists(config_file):
			config_values = fo.pickle_loader(config_file)
			self.translation_table = config_values[ct.TRANSLATION_TABLE_ARGNAME][0]

		return self


DistanceComputationMethod = Annotated[str,
                   					  Field(default=ct.DEFAULT_DISTANCE_METHOD),
                   					  AfterValidator(partial(validate_choice, valid_options=ct.DISTANCE_METHODS, exception_str=ct.INVALID_DISTANCE_METHOD))]

OutputFormat = Annotated[str,
                   Field(default=ct.DEFAULT_OUTPUT_FORMAT),
				   AfterValidator(partial(validate_choice, valid_options=ct.OUTPUT_FORMATS, exception_str=ct.INVALID_OUTPUT_FORMAT))]


class ComputeDistancesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: AllelicProfiles
	distance_computation_method: DistanceComputationMethod
	output_format: OutputFormat
	compute_similarity: BooleanFlag
	no_mask: BooleanFlag
	cpu_cores: CPUCores


MSAInput = Annotated[str,
					 AfterValidator(partial(input_path_exists,
								    function=fo.exists,
									exception_str=ct.MISSING_INPUT_ARG))]

Gaps = Annotated[str,
                 Field(default=ct.DEFAULT_GAPS),
				 AfterValidator(partial(validate_choice, valid_options=ct.GAPS_CHOICES, exception_str=ct.COMPUTEMSA_INVALID_GAPS_OPTION))]

Ambiguous = Annotated[str,
                      Field(default=ct.DEFAULT_AMBIGUOUS),
					  AfterValidator(partial(validate_choice, valid_options=ct.AMBIGUOUS_CHOICES, exception_str=ct.COMPUTEMSA_INVALID_AMBIGUOUS_OPTION))]

CustomMafftOptions = Annotated[str,
                               Field(default=None)]

class ComputeMSAValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_path: MSAInput
	schema_directory: ExternalSchema
	dna_msa: BooleanFlag
	output_variable: BooleanFlag
	translation_table: TranslationTable
	cpu_cores: CPUCores
	only_loci_msas: BooleanFlag
	gaps: Gaps
	ambiguous: Ambiguous
	custom_mafft_options: CustomMafftOptions
	protein_input: BooleanFlag
	no_cleanup: BooleanFlag

	# Use the translation table defined in the training file if main input is a matrix with allelic profiles
	@model_validator(mode="after")
	def get_translation_table(self):
		if self.schema_directory:
			config_file = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
			# Allow for schemas that were not created with chewBBACA
			if fo.exists(config_file):
				config_values = fo.pickle_loader(config_file)
				self.translation_table = config_values[ct.TRANSLATION_TABLE_ARGNAME][0]

		return self


SpeciesId = Annotated[str,
                      Field(default=None)]

SchemaId = Annotated[str,
                     Field(default=None)]

NomenclatureServerInstance = Annotated[str | None,
                               		   Field(default=ct.HOST_NS[ct.DEFAULT_NOMENCLATURE_SERVER])]


SchemaDate = Annotated[str,
                       Field(default=None),
					   AfterValidator(pd.validate_date)] # Determine if user date is in a valid format


class DownloadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	species_id: SpeciesId
	schema_id: SchemaId
	nomenclature_server_instance: NomenclatureServerInstance
	schema_date: SchemaDate
	latest_version: BooleanFlag
	blast_path : BLASTPath
	cpu_cores: CPUCores

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self

	# Check if user provided a valid integer identifier or species name
	@model_validator(mode="after")
	def validate_species_id(self):
		species_info = cr.species_ids(self.species_id, self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		if isinstance(species_info, list):
			self.species_id = species_info
		else:
			sys.exit(f"There is no species with the identifier {self.species_id} in the specified "
			         f"Chewie-NS instance ({self.nomenclature_server_instance}).")

		return self

	@model_validator(mode="after")
	def validate_schema_id(self):
		# The function return a list with the schema ID, the schema URI, the schema name, and the schema parameteres
		schema_info = cr.get_species_schemas(self.schema_id, self.species_id, self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		self.schema_id = schema_info

		return self


DescriptionFile = Annotated[str,
                            Field(default=None),
                            AfterValidator(partial(input_path_exists,
                                                   function=fo.exists,
												   exception_str=ct.MISSING_INPUTS_EXCEPTION))]


class UploadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: ChewieSchema
	species_id: SpeciesId
	schema_name: str
	loci_prefix: str
	description_file: DescriptionFile
	annotations_file: Annotations | None = None
	nomenclature_server_instance: NomenclatureServerInstance
	continue_up: BooleanFlag
	cpu_cores: CPUCores

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self

	# Check if user provided a valid integer identifier or species name
	@model_validator(mode="after")
	def validate_species_id(self):
		species_info = cr.species_ids(self.species_id, self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		if isinstance(species_info, list):
			self.species_id = species_info
		else:
			sys.exit(f"There is no species with the identifier {self.species_id} in the specified "
			         f"Chewie-NS instance ({self.nomenclature_server_instance}).")

		return self

	# Read config file and validate values
	@model_validator(mode="after")
	def validate_schema_config(self):
		print("Validating the argument values stored in the schema config...")
		config_file = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
		config_values = fo.pickle_loader(config_file)
		if any([len(v) > 1 for k, v in config_values.items()]):
			sys.exit("Schemas can only be uploaded to Chewie-NS if they were used with a "
			         "single argument per parameter defined in the '.schema_config file'.")
		# Validate arguments if there is a single value per parameter
		# Pyrodigal training file
		ptf_hash = config_values["pyrodigal_training_file"][0]
		schema_files = fo.listdir(self.schema_directory, full_path=True)
		ptf_paths = fo.filter_by_extension(schema_files, [".trn"])[0]
		if len(ptf_paths) > 1:
			sys.exit(f"Schema includes multiple training files. Please ensure that the schema includes a single training file.")

		if fo.hash_file(ptf_paths[0], "blake2b") != ptf_hash:
			sys.exit(f"Training file in schema's directory does not match the one in the config file.")

		# Validate remaining argument values using Pydantic's TypeAdapter
		bsr_adapter = TypeAdapter(BLASTScoreRatio).validate_python(config_values["blast_score_ratio"][0])
		ml_adapter = TypeAdapter(MinimumLength).validate_python(config_values["minimum_length"][0])
		st_adapter = TypeAdapter(SizeThreshold).validate_python(config_values["size_threshold"][0])
		tt_adapter = TypeAdapter(TranslationTable).validate_python(config_values["translation_table"][0])
		wos_adapter = TypeAdapter(WordSize).validate_python(config_values["word_size"][0])
		win_adapter = TypeAdapter(WindowSize).validate_python(config_values["window_size"][0])
		cs_adapter = TypeAdapter(ClusteringSimilarity).validate_python(config_values["clustering_sim"][0])
		rf_adapter = TypeAdapter(RepresentativeFilter).validate_python(config_values["representative_filter"][0])
		if_adapter = TypeAdapter(IntraFilter).validate_python(config_values["intra_filter"][0])

		validated_config = {
			"bsr": config_values["blast_score_ratio"][0],
			"minimum_locus_length": config_values["minimum_length"][0],
			"size_threshold": config_values["size_threshold"][0],
			"translation_table": config_values["translation_table"][0],
			"prodigal_training_file": [ptf_paths[0], ptf_hash],
			"word_size": config_values["word_size"][0],
			"cluster_sim": config_values["clustering_sim"][0],
			"representative_filter": config_values["representative_filter"][0],
			"intraCluster_filter": config_values["intra_filter"][0],
			"chewBBACA_version": f"chewBBACA {config_values['chewie_version'][0]}"
		}

		self.schema_directory = [self.schema_directory, validated_config]

		return self


class SynchronizeSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: ChewieSchema
	nomenclature_server_instance: NomenclatureServerInstance
	submit_alleles: BooleanFlag
	blast_path : BLASTPath
	cpu_cores: CPUCores

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [ct.NS_CONFIG_BASENAME])
		if not fo.is_file(ns_config)[0]:
			sys.exit("Schema does not contain the '.ns_config' file.")

		return self


StatsMode = Annotated[str,
                      Field(default=None),
					  AfterValidator(partial(validate_choice,
					                         valid_options=ct.NSSTATS_MODE_CHOICES,
											 exception_str=ct.NS_STATS_INVALID_MODE))]


class NSStatsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	stats_mode: StatsMode
	species_id: SpeciesId
	schema_id: SchemaId
	nomenclature_server_instance: NomenclatureServerInstance

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, ct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self
