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
import shutil
import argparse
import platform
import subprocess
import multiprocessing
from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   chewiens_requests as cr,
					   fasta_operations as fao,
					   pyrodigal_gene_prediction as pgp)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (file_operations as fo,
								 chewiens_requests as cr,
								 fasta_operations as fao,
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


def validate_ptf_hash(ptf_hash, schema_ptfs, force_continue):
	""" Determines if the hash for the Prodigal training
		file matches any of the hashes from training files
		that have been used with the schema.

	Paramters
	---------
	ptf_hash : str
		BLAKE2b hash computed based on the contents of
		the training file.
	schema_ptfs : list
		List with the hashes of all training files that
		have been used with the schema.
	force_continue : bool
		True if the hash should be added to the list with
		all hashes from training files used with the
		schema without prompting the user. False otherwise.

	Returns
	-------
	unmatch : bool
		True if the hash is not in the list with all hashes
		from all training files used with the schema. False
		otherwise.

	Raises
	------
	SystemExit
		- If the user does not agree to add the hash from a
		new training file to the list with all hashes for
		training files that have been used with the schema.
	"""
	unmatch = False
	if ptf_hash not in schema_ptfs:
		ptf_num = len(schema_ptfs)
		if force_continue is False:
			if ptf_num == 1:
				ptf_answer = fo.input_timeout(ct.DIFFERENT_PTF_PROMPT, ct.PROMPT_TIMEOUT)
			if ptf_num > 1:
				ptf_answer = fo.input_timeout(ct.MULTIPLE_PTF_PROMPT.format(ptf_num), ct.PROMPT_TIMEOUT)
		else:
			ptf_answer = 'yes'

		if ptf_answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			unmatch = True

	return unmatch


def validate_ptf(ptf_path, schema_directory, schema_ptfs, force_continue):
	""" Validates the path to the Prodigal training file and
		its hash value.

	Parameters
	----------
	ptf_path : str or NoneType
		Path to the Prodigal training file or NoneType
		if no value was provided.
	schema_directory : str
		Path to the schema's directory.
	schema_ptfs : list
		List with the hashes of all training files that
		have been used with the schema.
	force_continue : bool
		True if the path and hash of the training file
		should be validated without prompting the user.
		False otherwise.

	Returns
	-------
	ptf_path : str or bool
		Path to the training file if the user provided a
		valid path or if no value was provided and the
		schema has a training file. False if the user
		passed 'False' or if no value was passed and the
		schema has no training file.
	ptf_hash : str
		BLAKE2b hash computed based on the contents of
		the training file.
	unmatch : bool
		True if the training file does not match any of
		the training files previously used with the schema.
	"""
	ptf_path = validate_ptf_path(ptf_path, schema_directory)

	# Determine PTF checksum
	if ptf_path is not None:
		ptf_hash = fo.hash_file(ptf_path, 'blake2b')
	else:
		ptf_hash = None

	unmatch = validate_ptf_hash(ptf_hash, schema_ptfs, force_continue)

	return [ptf_path, ptf_hash, unmatch]


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


def input_is_fasta(input_path):
	""""""
	# Check if input is a single FASTA file
	if fo.is_file(input_path)[0]:
		if fao.validate_fasta(input_path):
			sys.exit(ct.FASTA_INPUT_EXCEPTION)


def list_input_files(input_path):
	""""""
	if os.path.isfile(input_path):
		# Read list of input files
		input_files = [line[0] for line in fo.read_tabular(input_path)]
	# Input path is for a directory
	elif os.path.isdir(input_path):
		# List absolute paths
		# Only keep paths to files
		input_files = [file for file in fo.listdir_fullpath(input_path) if os.path.isdir(file) is False]

	return input_files


def filter_inputs_extension(input_files, extensions=ct.FASTA_EXTENSIONS):
	""""""
	# Need to verify if files end with any of the accepted file extensions
	valid_extension, invalid_extension = fo.filter_by_extension(input_files, extensions)
	if len(invalid_extension) > 0:
		sys.exit(ct.INVALID_EXTENSION_EXCEPTION)

	return input_files


def inputs_exist(input_files):
	""""""
	# Check that all files exist
	missing_inputs = [file for file in input_files if fo.exists(file) is False]
	if len(missing_inputs) > 0:
		sys.exit(ct.MISSING_INPUTS)

	return input_files


def validate_inputs_fastas(input_files):
	""""""
	# Input files must be vaid FASTA files
	fasta, non_fasta = fao.filter_non_fasta(input_files)
	if len(non_fasta) > 0:
		sys.exit(ct.NON_FASTA_EXCEPTION)

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
		repeated_basenames = [f'{k}: {", ".join(v)}' for k, v in prefixes.items() if len(v) > 1]
		repeated_basenames = [','.join(l) for l in repeated_basenames]
		sys.exit(ct.INPUTS_SHARE_PREFIX.format('\n'.join(repeated_basenames)))

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
		sys.exit(ct.INPUTS_INCLUDE_BLANKS.format('\n'.join(include_blanks)))

	return input_files


def translation_table_type(genetic_code, valid_genetic_codes):
	"""Validate the translation table value (TT) passed to chewBBACA.

	Parameters
	----------
	arg : int
		The TT value passed to chewBBACA. Must be of type int
		and match the identifier of one of the genetic codes.
	genetic_codes : dict
		Dictionary with genetic codes identifiers as keys and
		descriptions as values.

	Returns
	-------
	valid : int
		The TT value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the TT value cannot be converted to int type
		or if it does not match any of the acceptable genetic
		codes.
	"""
	# Set to default value if user did not provide a value
	if genetic_code not in valid_genetic_codes:
		sys.exit(ct.INVALID_GENETIC_CODE)

	return genetic_code


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
	for line in stderr.split("\n"):
		if "|" in line:
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
		sys.exit(ct.AUGUSTUS_INVALID_SPECIES)

	return species_id


def validate_augustus_outfmt(outfmt, valid_outfmts=ct.AUGUSTUS_OUTFMTS):
	""""""
	if outfmt not in valid_outfmts:
		sys.exit(ct.AUGUSTUS_INVALID_OUTFMT)

	return outfmt


def input_is_fasta(input_path):
	""""""
	# Check if input is a single FASTA file
	if fo.is_file(input_path)[0]:
		if fao.validate_fasta(input_path):
			sys.exit(ct.FASTA_INPUT_EXCEPTION)

	return input_path


def add_ptf_genetic_code(ptf_path, genetic_code):
	"""
	"""
	if ptf_path:
		print(f"Provided training file. Ignoring the translation table "
			  "value previously set ({genetic_code}) and using the genetic "
			  "code defined in the training file ({genetic_code}).")
		# Get translation table used to create training file
		genetic_code = pgp.read_training_file(ptf_path).translation_table

	return genetic_code


def validate_pyrodigal_mode(mode, valid_modes=ct.PYRODIGAL_MODES):
	""""""
	if mode not in valid_modes:
		sys.exit(ct.INVALID_PYRODIGAL_MODE)

	return mode


def check_meta(pyrodigal_mode, pyrodigal_training_file):
	"""
	"""
	if pyrodigal_mode == 'meta' and pyrodigal_training_file is not None:
		sys.exit(ct.PYRODIGAL_META_NOPTF)

	return pyrodigal_mode


def validate_pyrodigal_outfmt(outfmt, valid_outfmts=ct.PYRODIGAL_OUTFMTS):
	""""""
	if not all([of in valid_outfmts for of in outfmt]):
		sys.exit(ct.PYRODIGAL_INVALID_OUTFMT)

	return outfmt


def validate_pyrodigal_minimum_confidence(confidence, min_value=ct.PYRODIGAL_MIN_CONFIDENCE, max_value=ct.PYRODIGAL_MAX_CONFIDENCE):
	""""""
	if confidence:
		if confidence < min_value or confidence > max_value:
			sys.exit(ct.INVALID_PYRODIGAL_CONFIDENCE)

	return confidence


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
	if arg == multiprocessing.cpu_count():
		print(ct.CPU_VALUE_WARNING)

	return arg


def validate_word_size(word_size, min_value=ct.WORD_SIZE_MIN, max_value=ct.WORD_SIZE_MAX):
	"""Validate the word size value (WS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The WS value passed to chewBBACA.
	min_value : float
		Minimum acceptable WS value.
	max_value : float
		Maximum acceptable WS value.
	default_value : float
		The default WS value to use if none is provided.

	Returns
	-------
	valid : float
		The WS value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the WS value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	if word_size < min_value or word_size > max_value:
		sys.exit(ct.INVALID_WORD_SIZE)

	return word_size


def validate_window_size(window_size, min_value=ct.WINDOW_SIZE_MIN, max_value=ct.WINDOW_SIZE_MAX):
	"""Validate the window size value (WS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The WS value passed to chewBBACA.
	min_value : float
		Minimum acceptable WS value.
	max_value : float
		Maximum acceptable WS value.
	default_value : float
		The default WS value to use if none is provided.

	Returns
	-------
	valid : float
		The WS value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the WS value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	if window_size < min_value or window_size > max_value:
		sys.exit(ct.INVALID_WINDOW_SIZE)

	return window_size


def validate_clustering_similarity(clustering_similarity, min_value=ct.CLUSTERING_SIMILARITY_MIN, max_value=ct.CLUSTERING_SIMILARITY_MAX):
	"""Validate the clustering similarity value (WS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The WS value passed to chewBBACA.
	min_value : float
		Minimum acceptable WS value.
	max_value : float
		Maximum acceptable WS value.
	default_value : float
		The default WS value to use if none is provided.

	Returns
	-------
	valid : float
		The WS value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the WS value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	if clustering_similarity < min_value or clustering_similarity > max_value:
		sys.exit(ct.INVALID_CLUSTERING_SIMILARITY)

	return clustering_similarity


def validate_representative_filter(representative_filter, min_value=ct.REPRESENTATIVE_FILTER_MIN, max_value=ct.REPRESENTATIVE_FILTER_MAX):
	"""Validate the representative filter value (WS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The WS value passed to chewBBACA.
	min_value : float
		Minimum acceptable WS value.
	max_value : float
		Maximum acceptable WS value.
	default_value : float
		The default WS value to use if none is provided.

	Returns
	-------
	valid : float
		The WS value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the WS value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	if representative_filter < min_value or representative_filter > max_value:
		sys.exit(ct.INVALID_REPRESENTATIVE_FILTER)

	return representative_filter


def validate_intra_filter(intra_cluster, min_value=ct.INTRA_CLUSTER_MIN, max_value=ct.INTRA_CLUSTER_MAX):
	"""Validate the intra filter value (WS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The WS value passed to chewBBACA.
	min_value : float
		Minimum acceptable WS value.
	max_value : float
		Maximum acceptable WS value.
	default_value : float
		The default WS value to use if none is provided.

	Returns
	-------
	valid : float
		The WS value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the WS value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	if intra_cluster < min_value or intra_cluster > max_value:
		sys.exit(ct.INVALID_INTRA_CLUSTER_FILTER)

	return intra_cluster


def parse_parameter_string(input_string):
	"""
	"""
	config_args = {}
	if input_string is not None:
		for v in input_string:
			parameter, argument = v.split("=")
			parameter = parameter.replace("-", "_")
			# Split argument value in by "," to get multiple values
			if "," in argument:
				argument = argument.split(",")
			# Add parameter name and argument value to config dictionary
			config_args[parameter] = argument

	return config_args


def training_file_exists(training_file_path):
	""""""
	if training_file_path:
		if not fo.exists(training_file_path):
			sys.exit(ct.INVALID_INPUT_PATH)

	return training_file_path


def validate_gene_predictor(gene_predictor, valid_gene_predictors=ct.GENE_PREDICTORS):
	""""""
	if gene_predictor not in valid_gene_predictors:
		sys.exit(ct.INVALID_GENE_PREDICTOR)

	return gene_predictor


# Check if parent directory for the provided output_dir exists
# field_validator is used to check before the class is instantiated
def parentdir_exists(input_path):
	parent_dir = os.path.dirname(input_path)
	if not os.path.isdir(parent_dir):
		raise ValueError(f"Parent directory for '{input_path}' does not exist.")
	return input_path


def input_path_exists(input_path):
	if not fo.exists(input_path):
		sys.exit(f"Input path '{input_path}' does not exist.")
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


def validate_allelecall_mode(mode, valid_modes):
	""""""
	if mode not in valid_modes:
		sys.exit(ct.INVALID_ALLELECALL_MODE)

	return mode


def schema_includes_fasta(schema_directory):
	""""""
	schema_files = os.listdir(schema_directory)
	# Check if the folder includes FASTA files
	fasta_extension, _ = fo.filter_by_extension(schema_files, ct.FASTA_EXTENSIONS)
	if len(fasta_extension) == 0:
		sys.exit(ct.MISSING_SCHEMA_FASTAS)

	return schema_directory


def check_schema(schema_directory):
	""""""
	schema_files = os.listdir(schema_directory)
	# Check if the "short" directory exists
	if "short" not in schema_files:
		sys.exit(ct.SCHEMA_INVALID_PATH)
	# Check if schema includes the .schema_config file
	config_file = fo.join_paths(schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
	if not fo.is_file(config_file)[0]:
		sys.exit(ct.ADAPT_LEGACY_SCHEMA)

	return schema_directory


def check_bsr_conflict(user_bsr, schema_config, force_continue):
	""""""
	if user_bsr not in schema_config["bsr"]:
		print("The value provided for the BLAST Score Ratio does not match "
			  f"any of the values used with the schema ({schema_config["bsr"]}).")
		if not force_continue:
			proceed = fo.input_timeout(ct.ARGS_DIFFER_PROMPT, ct.PROMPT_TIMEOUT)
		else:
			params_answer = 'yes'

		if params_answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			schema_config["minimum_locus_bsrlength"].append(user_bsr)

	return user_bsr, schema_config


def check_ml_conflict(user_ml, schema_config, force_continue):
	""""""
	if user_ml not in schema_config["minimum_locus_length"]:
		print("The value provided for the minimum locus length does not match "
			  f"any of the values used with the schema ({schema_config["minimum_locus_length"]}).")
		if not force_continue:
			proceed = fo.input_timeout(ct.ARGS_DIFFER_PROMPT, ct.PROMPT_TIMEOUT)
		else:
			params_answer = 'yes'

		if params_answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			schema_config["minimum_locus_length"].append(user_ml)

	return user_ml, schema_config


def check_st_conflict(user_st, schema_config, force_continue):
	""""""
	if user_st not in schema_config["size_threshold"]:
		print("The value provided for the size threshold does not match "
			  f"any of the values used with the schema ({schema_config["size_threshold"]}).")
		if not force_continue:
			proceed = fo.input_timeout(ct.ARGS_DIFFER_PROMPT, ct.PROMPT_TIMEOUT)
		else:
			params_answer = 'yes'

		if params_answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			schema_config["size_threshold"].append(user_st)

	return user_st, schema_config


def check_ptf_conflict(user_ptf, translation_table, schema_config, force_continue):
	""""""
	user_ptf_hash = fo.hash_file(user_ptf, 'blake2b')
	if user_ptf_hash not in schema_config["prodigal_training_file"]:
		print("The Pyrodigal training file provided does not match "
				f"any of the training files used with the schema ({schema_config["prodigal_training_file"]}).")
		if not force_continue:
			proceed = fo.input_timeout(ct.ARGS_DIFFER_PROMPT, ct.PROMPT_TIMEOUT)
		else:
			params_answer = 'yes'

		if params_answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			schema_config["prodigal_training_file"].append(user_ptf_hash)
			# Get genetic code used with the training file
			genetic_code = pgp.read_training_file(user_ptf).translation_table
			translation_table = genetic_code
			# Add genetic code to schema config if it was never used
			if translation_table not in schema_config["translation_table"]:
				schema_config["translation_table"].append(translation_table)

		return user_ptf, translation_table, schema_config


def contains_results(results_directory):
	""""""
	# List files in input direcoty
	results_files = os.listdir(results_directory)
	# Check if folder includes the results_alleles.tsv file with allelic profiles
	if ct.RESULTS_ALLELES_BASENAME not in results_files:
		sys.exit(ct.MISSING_RESULTS_ALLELES)


def validate_cgmlst_thresholds(threshold_values):
	""""""
	for t in threshold_values:
		if t < ct.CGMLST_THRESHOLD_MIN or t > ct.CGMLST_THRESHOLD_MAX:
			sys.exit("Invalid loci presence threshold.")

	return threshold_values


def validate_distance_method(method, valid_methods=ct.DISTANCE_METHODS):
	""""""
	if method not in valid_methods:
		sys.exit(ct.INVALID_DISTANCE_METHOD)

	return method


def validate_output_format(output_format, valid_output_formats=ct.OUTPUT_FORMATS):
	""""""
	if output_format not in valid_output_formats:
		sys.exit(ct.INVALID_OUTPUT_FORMAT)

	return output_format

# Define reusable fields
OutputDirectory = Annotated[str,
							AfterValidator(parentdir_exists),
							AfterValidator(create_output_directory)]

InputFiles = Annotated[str,
					   AfterValidator(input_path_exists),
					   AfterValidator(input_is_fasta),
					   AfterValidator(list_input_files),
					   AfterValidator(filter_inputs_extension),
					   AfterValidator(inputs_exist),
					   AfterValidator(validate_inputs_fastas),
					   AfterValidator(check_unique_prefixes),
					   AfterValidator(check_blanks)]

SchemaDirectory = Annotated[str,
							AfterValidator(input_path_exists),
							AfterValidator(schema_includes_fasta)]

LociList = Annotated[str, 
					 AfterValidator(input_path_exists),
					 AfterValidator(input_is_fasta),
					 AfterValidator(list_input_files),
					 AfterValidator(filter_inputs_extension),
					 AfterValidator(inputs_exist),
					 AfterValidator(validate_inputs_fastas)]

BLASTScoreRatio = Annotated[float,
							Field(default=ct.BSR_DEFAULT, ge=ct.BSR_MIN, le=ct.BSR_MAX)]

MinimumLength = Annotated[int,
						  Field(default=ct.MSL_DEFAULT, ge=ct.MSL_MIN, le=ct.MSL_MAX)]

TranslationTable = Annotated[int,
							 Field(default=ct.GENETIC_CODE_DEFAULT),
							 AfterValidator(partial(translation_table_type, valid_genetic_codes=ct.GENETIC_CODES))]

SizeThreshold = Annotated[float | None,
						  Field(default=ct.ST_DEFAULT, ge=ct.ST_MIN, le=ct.ST_MAX)]

GenePredictor = Annotated[str,
						  Field(default=ct.GENE_PREDICTOR_DEFAULT),
						  AfterValidator(partial(validate_gene_predictor, valid_gene_predictors=ct.GENE_PREDICTORS))]

GenePredictionArguments = Annotated[str | None,
									AfterValidator(partial(parse_parameter_string, parameter_types=ct.ARGUMENT_TYPES))]

ValidatedGenePredictionArguments = Annotated[PyrodigalArgs | AugustusArgs | None,
											 Field(default=None)]

ClusteringArguments = Annotated[str | None, 
								AfterValidator(partial(parse_parameter_string, parameter_types=ct.ARGUMENT_TYPES))]

ValidatedClusteringArguments = Annotated[ClusteringArgs | None,
										 Field(default=None)]

BLASTPath = Annotated[str | None,
					  AfterValidator(check_blast),
					  AfterValidator(check_blast_version)]

CPUCores = Annotated[int,
					 Field(default=ct.CPU_CORES_DEFAULT, ge=ct.CPU_CORES_MIN, le=multiprocessing.cpu_count()),
					 AfterValidator(verify_cpu_usage)]

PyrodigalTrainingFile = Annotated[str | None,
								  AfterValidator(training_file_exists)]

PyrodigalMode = Annotated[str,
						  Field(default=ct.PYRODIGAL_DEFAULT_MODE),
						  AfterValidator(partial(validate_pyrodigal_mode, valid_modes=ct.PYRODIGAL_MODES))]

PyrodigalOutputFormats = Annotated[str,
								   Field(default=ct.PYRODIGAL_DEFAULT_OUTFMT),
								   AfterValidator(partial(validate_pyrodigal_outfmt, valid_outfmts=ct.PYRODIGAL_OUTFMTS))]

PyrodigalMinimumConfidence = Annotated[float | None,
									   AfterValidator(validate_pyrodigal_minimum_confidence)]

PyrodigalTrainingReference = Annotated[str | None,
									   AfterValidator(training_file_exists)]

AugustusPath = Annotated[str | None,
						 AfterValidator(check_augustus)]

AugustusOutputFormats = Annotated[str,
								  Field(default=ct.AUGUSTUS_OUTFMT_DEFAULT),
								  AfterValidator(partial(validate_augustus_outfmt, valid_gene_predictors=ct.GENE_PREDICTORS))]

WordSize = Annotated[int,
					 Field(default=ct.WORD_SIZE_DEFAULT),
					 AfterValidator(validate_word_size)]

WindowSize = Annotated[int,
					   Field(default=ct.WORD_SIZE_DEFAULT),
					   AfterValidator(validate_window_size)]

ClusteringSimilarity = Annotated[float,
								 Field(default=ct.CLUSTERING_SIMILARITY_DEFAULT),
								 AfterValidator(validate_clustering_similarity)]

RepresentativeFilter = Annotated[float,
								 Field(default=ct.REPRESENTATIVE_FILTER_DEFAULT),
								 AfterValidator(validate_representative_filter)]

IntraFilter = Annotated[float,
						Field(default=(ct.INTRA_CLUSTER_DEFAULT)),
						AfterValidator(validate_intra_filter)]

Annotations = Annotated[str,
						AfterValidator(input_path_exists)]


class PyrodigalArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	pyrodigal_training_file: PyrodigalTrainingFile
	pyrodigal_mode: PyrodigalMode
	pyrodigal_output_formats: PyrodigalOutputFormats
	pyrodigal_minimum_confidence: PyrodigalMinimumConfidence
	pyrodigal_training_reference: PyrodigalTrainingReference


class AugustusArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	augustus_path: AugustusPath
	augustus_species: str | None
	augustus_output_formats: AugustusOutputFormats

	@model_validator(mode="after")
	def validate_species(self):
		self.augustus_species = validate_augustus_species(self.augustus_species, self.augustus_path)

		return self


class ClusteringArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	word_size: WordSize
	window_size: WindowSize
	clustering_similarity: ClusteringSimilarity
	representative_filter: RepresentativeFilter
	intra_filter: IntraFilter


class PredictGenesValidator(BaseModel):
	# Use `from_attributes=True` to inspect the attributes of the 
	# argparse.Namespace object directly without having to convert to a dictionary
	model_config = ConfigDict(from_attributes=True)

	# Field order is preserved when serializing with model_dump()
	output_directory: OutputDirectory
	input_files: InputFiles
	gene_predictor: GenePredictor
	gene_prediction_arguments: GenePredictionArguments
	validated_gene_prediction_arguments: ValidatedGenePredictionArguments
	translation_table: TranslationTable
	cpu_cores: CPUCores

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_arguments = PyrodigalArgs(**self.gene_prediction_arguments)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_arguments = AugustusArgs(**self.gene_prediction_arguments)

		return self

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		self.validated_gene_prediction_arguments.pyrodigal_mode = check_meta(self.validated_gene_prediction_arguments.pyrodigal_mode, self.validated_gene_prediction_arguments.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def ptf_and_reference(self):
		# Check if user provided path to Pyrodigal training file and to a training reference, which is not allowed
		check_ptf_tref_conflict(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.validated_gene_prediction_arguments.training_reference)
		# Create training file based on training reference
		print(f'Creating Pyrodigal training file based on {self.validated_gene_prediction_arguments.training_reference}...')
		self.validated_gene_prediction_arguments.pyrodigal_training_file = pgp.create_training_file(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.validated_gene_prediction_arguments.training_reference, self.output_directory, self.translation_table)
		print(f'Training file saved to {self.validated_gene_prediction_arguments.pyrodigal_training_file}')


# I can use reusable fields like this one
# custom_option = Annotated[str, AfterValidator(validation_function)]
class CreateSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_files: InputFiles
	schema_name: Annotated[str, Field(default=ct.SCHEMA_NAME_DEFAULT)]
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	gene_predictor : GenePredictor
	gene_prediction_arguments: GenePredictionArguments
	validated_gene_prediction_arguments: ValidatedGenePredictionArguments
	clustering_parameters: ClusteringArguments
	validated_clustering_arguments: ValidatedClusteringArguments
	blast_path : BLASTPath
	cds_input: bool
	no_cds_renaming: bool
	cpu_cores: CPUCores
	no_cleanup: bool

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_arguments = PyrodigalArgs(**self.gene_prediction_arguments)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_arguments = AugustusArgs(**self.gene_prediction_arguments)

		return self

	# Further validation for clustering arguments
	@model_validator(mode="after")
	def validate_clustering_arguments(self):
		self.validated_clustering_arguments = ClusteringArgs(**self.clustering_parameters)

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		self.validated_gene_prediction_arguments.pyrodigal_mode = check_meta(self.validated_gene_prediction_arguments.pyrodigal_mode, self.validated_gene_prediction_arguments.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def ptf_and_reference(self):
		# Check if user provided path to Pyrodigal training file and to a training reference, which is not allowed
		check_ptf_tref_conflict(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.validated_gene_prediction_arguments.training_reference)
		# Create training file based on training reference
		print(f'Creating Pyrodigal training file based on {self.validated_gene_prediction_arguments.training_reference}...')
		self.validated_gene_prediction_arguments.pyrodigal_training_file = pgp.create_training_file(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.validated_gene_prediction_arguments.training_reference, self.output_directory, self.translation_table)
		print(f'Training file saved to {self.validated_gene_prediction_arguments.pyrodigal_training_file}')


class AlleleCallValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_files: InputFiles
	schema_directory: SchemaDirectory
	# Do not run validation if loci_list is None
	loci_list: LociList | None = None
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	gene_predictor : GenePredictor
	gene_prediction_arguments: GenePredictionArguments
	validated_gene_prediction_arguments: ValidatedGenePredictionArguments
	clustering_parameters: ClusteringArguments
	validated_clustering_arguments: ValidatedClusteringArguments
	blast_path : BLASTPath
	cds_input: bool
	no_inferred: bool
	output_unclassified: bool
	output_missing: bool
	output_novel: bool
	output_masked: bool
	no_cds_renaming: bool
	force_continue: bool
	mode: Annotated[int, Field(default=ct.ALLELECALL_DEFAULT_MODE), AfterValidator(partial(validate_allelecall_mode, ct.ALLELECALL_MODES))]
	cpu_cores: CPUCores
	no_cleanup: bool
	ns_config: Annotated[bool, Field(default=False)]

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_arguments = PyrodigalArgs(**self.gene_prediction_arguments)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_arguments = AugustusArgs(**self.gene_prediction_arguments)

		return self

	# Further validation for clustering arguments
	@model_validator(mode="after")
	def validate_clustering_arguments(self):
		self.validated_clustering_arguments = ClusteringArgs(**self.clustering_parameters)

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		self.validated_gene_prediction_arguments.pyrodigal_mode = check_meta(self.validated_gene_prediction_arguments.pyrodigal_mode, self.validated_gene_prediction_arguments.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def validate_schema(self):
		# Check that the schema includes the short directory and the .schema_config file necessary for allele calling
		self.schema_directory = check_schema(self.schema_directory)

	@model_validator(mode="after")
	def list_schema_loci(self):
		# List all loci FASTA files in the schema if no loci list was provided
		if not self.loci_list:
			self.loci_list = fo.listdir_fullpath(self.schema_directory, substring_filter=".fasta")

		return self

	@model_validator
	def solve_conflicting_arguments(self):
		config_file = fo.join_paths(self.schema_directory, [ct.SCHEMA_CONFIG_BASENAME])
		config = fo.pickle_loader(config_file)
		# BLAST Score Ratio
		self.blast_score_ratio, config = check_bsr_conflict(self.blast_score_ratio, config, self.force_continue)
		# Minimum length
		self.minimum_length, config = check_ml_conflict(self.minimum_length, config, self.force_continue)
		# Size threshold
		self.size_threshold, config = check_st_conflict(self.size_threshold, config, self.force_continue)
		# Pyrodigal training file
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_arguments.pyrodigal_training_file, self.translation_table, config = check_ptf_conflict(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.translation_table, config, self.force_continue)

		# Update schema config file
		fo.pickle_dumper(config, config_file)

		return self

	@model_validator
	def check_ns_config(self):
		ns_config = fo.join_paths(self.output_directory, [ct.NS_CONFIG_BASENAME])
		if fo.is_file(ns_config)[0]:
			self.ns_config = True


class SchemaEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: SchemaDirectory
	loci_list: LociList | None = None
	annotations: Annotations | None = None
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	minimum_length: MinimumLength
	cpu_cores: CPUCores
	loci_reports: Annotated[bool, Field(default=False)]
	light: Annotated[bool, Field(default=False)]
	add_sequences: Annotated[bool, Field(default=False)]


class AlleleCallEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_files: Annotated[str, AfterValidator(input_path_exists), AfterValidator(contains_results)]
	schema_directory: SchemaDirectory
	annotations: Annotations | None = None
	cpu_cores: CPUCores
	light: Annotated[bool, Field(default=False)]
	no_pa: Annotated[bool, Field(default=False)]
	no_dm: Annotated[bool, Field(default=False)]
	no_tree: Annotated[bool, Field(default=False)]
	cg_alignment: Annotated[bool, Field(default=False)]


class ExtractCgMLSTValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_files: Annotated[str, AfterValidator(input_path_exists), AfterValidator(contains_results)]
	threshold: Annotated[list, Field(default=ct.CGMLST_THRESHOLDS), AfterValidator(validate_cgmlst_thresholds)]
	step: Annotated[int, Field(default=1)]
	compute_accessory: Annotated[bool, Field(default=False)]
	rarefaction_analysis: Annotated[bool, Field(default=False)]
	permutation_number: Annotated[int, Field(default=ct.PERMUTATION_NUMBER_DEFAULT)]
	permutation_samples: int | None = None
	exclude_loci: Annotated[str, AfterValidator(input_path_exists)] | None = None
	exclude_genomes: Annotated[str, AfterValidator(input_path_exists)] | None = None
	cpu_cores: CPUCores


class SubsetResultsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_files: Annotated[str, AfterValidator(input_path_exists), AfterValidator(contains_results)]
	loci_list: LociList | None = None
	sample_list: Annotated[str, AfterValidator(input_path_exists)] | None = None
	inverse_loci: Annotated[bool, Field(default=False)]
	inverse_samples: Annotated[bool, Field(default=False)]


class MergeResults(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_files: Annotated[str, AfterValidator(input_path_exists), AfterValidator(contains_results)]
	common: Annotated[bool, Field(default=False)]


class HashProfilesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: Annotated[str, AfterValidator(input_path_exists)]
	schema_directory: SchemaDirectory
	hash_type: Annotated[str, Field(default="crc32")]
	nrows: Annotated[int, Field(default=100)]
	cpu_cores: CPUCores


class GetAllelesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: Annotated[str, AfterValidator(input_path_exists)]
	schema_directory: SchemaDirectory
	loci_list: LociList | None = None
	cpu_cores: CPUCores
	distinct: Annotated[bool, Field(default=False)]
	translate: Annotated[bool, Field(default=False)]
	translation_table: TranslationTable

	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		self.translation_table = add_ptf_genetic_code(self.validated_gene_prediction_arguments.pyrodigal_training_file, self.translation_table)

		return self


class PrepExternalSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: SchemaDirectory
	loci_list: LociList | None = None
	gene_predictor: GenePredictor
	gene_prediction_arguments: GenePredictionArguments
	validated_gene_prediction_arguments: ValidatedGenePredictionArguments
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	adaptation_minimum_length: Annotated[int, Field(default=0)]
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	adaptation_size_threshold: Annotated[float, Field(default=None)]
	cpu_cores: CPUCores
	blast_path : BLASTPath
	size_filter: Annotated[bool, Field(default=False)]

	@model_validator(mode="after")
	def apply_size_filter(self):
		if self.size_filter:
			self.adaptation_minimum_length = self.minimum_length
			self.adaptation_size_threshold = self.size_threshold


class UniprotFinderValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)
	
	output_directory: OutputDirectory
	schema_directory: SchemaDirectory
	loci_list: LociList | None = None
	protein_table: Annotated[str, AfterValidator(input_path_exists)]
	blast_score_ratio: BLASTScoreRatio
	cpu_cores: CPUCores
	taxa: Annotated[str, Field(default=None)]
	proteome_matches: Annotated[int, Field(default=1)]
	no_sparql: Annotated[bool, Field(default=False)]
	no_cleanup: Annotated[bool, Field(default=False)]
	no_cleanup: bool


class ComputeDistancesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: Annotated[str, AfterValidator(input_path_exists)]
	method: Annotated[str, Field(default=ct.DEFAULT_DISTANCE_METHOD), AfterValidator(validate_distance_method)]
	outfmt: Annotated[str, Field(default=ct.DEFAULT_OUTPUT_FORMAT), AfterValidator(validate_output_format)]
	no_mask: Annotated[bool, Field(default=False)]
	similarity: Annotated[bool, Field(default=False)]
	cpu_cores: CPUCores


class ComputeMSAValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_path: Annotated[str, AfterValidator(input_path_exists)]
	schema_directory: SchemaDirectory
	dna_msa: Annotated[bool, Field(default=False)]
	output_variable: Annotated[bool, Field(default=0)]
	translation_table: TranslationTable
	cpu_cores: CPUCores
	only_loci_msas: Annotated[bool, Field(default=False)]
	gaps: Annotated[str, Field(default=ct.DEFAULT_GAPS), AfterValidator(validate_choice)]
	ambiguous: Annotated[str, Field(default=ct.DEFAULT_AMBIGUOUS), AfterValidator(validate_choice)]
	custom_mafft_parameters: Annotated[str, AfterValidator()]
	protein_input: Annotated[bool, Field(default=False)]
	no_cleanup: Annotated[bool, Field(default=False)]


class DownloadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	species_id: Annotated[str, AfterValidator()]
	schema_id: Annotated[str, AfterValidator()]
	download_folder: Annotated[str, AfterValidator(input_path_exists)]
	cpu_cores: CPUCores
	nomenclature_server: Annotated[str, Field(default=ct.DEFAULT_NOMENCLATURE_SERVER), AfterValidator(validate_choice)]
	blast_path : BLASTPath
	date: Annotated[str, AfterValidator()] | None = None
	latest: Annotated[bool, Field(default=False)]


class UploadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: SchemaDirectory
	species_id: Annotated[str, AfterValidator()]
	schema_name: Annotated[str, AfterValidator()]
	loci_prefix: Annotated[str, AfterValidator()]
	description_file: Annotated[str, AfterValidator(input_path_exists)]
	annotations: Annotations | None = None
	cpu_cores: CPUCores
	nomenclature_server: Annotated[str, Field(default=ct.DEFAULT_NOMENCLATURE_SERVER), AfterValidator(validate_choice)]
	continue_up: Annotated[bool, Field(default=False)]


class SynchronizeSchemaValidator():
	model_config = ConfigDict(from_attributes=True)

	schema_directory: SchemaDirectory
	cpu_cores: CPUCores
	nomenclature_server: Annotated[str, Field(default=ct.DEFAULT_NOMENCLATURE_SERVER), AfterValidator(validate_choice)]
	blast_path : BLASTPath
	submit: Annotated[bool, Field(default=False)]


class NSStatsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	mode: Annotated[str, AfterValidator(validate_choice)]
	species_id:
	schema_id:
	nomenclature_server: Annotated[str, Field(default=ct.DEFAULT_NOMENCLATURE_SERVER), AfterValidator(validate_choice)]
