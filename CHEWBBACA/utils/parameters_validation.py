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
import hashlib
import argparse
import platform
import subprocess
import multiprocessing

try:
	from utils import (constants as ct,
					   blast_wrapper as bw,
					   gene_prediction as gp,
					   file_operations as fo,
					   chewiens_requests as cr,
					   fasta_operations as fao,
					   iterables_manipulation as im)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 blast_wrapper as bw,
								 gene_prediction as gp,
								 file_operations as fo,
								 chewiens_requests as cr,
								 fasta_operations as fao,
								 iterables_manipulation as im)


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


def arg_list(arg, arg_name):
	"""Determine if more than one value has been used for a single parameter.

	Parameter
	---------
	arg : list
		List with the values used for a sinlge parameter to perform
		allele calling.
	arg_name : str
		The name of the parameter to include in the exception
		message if more than one parameter value has been used
		to perform allele calling.

	Returns
	-------
	The single parameter value used to perform allele calling.

	Raises
	------
	SystemExit
		- If more than one parameter value has been used to
		perform allele calling.
	"""
	if isinstance(arg, list) is True:
		if len(arg) > 1:
			sys.exit('\nMultiple {0} values.'.format(arg_name))
		else:
			arg = arg[0]

	return arg


def bsr_type(arg, min_value=ct.BSR_MIN, max_value=ct.BSR_MAX):
	"""Validate the BLAST Score Ratio (BSR) value passed to chewBBACA.

	Parameters
	----------
	arg : float
		The BLAST Score Ratio (BSR) value passed to chewBBACA.
	min_value : float
		Minimum acceptable BSR value.
	max_value : float
		Maximum acceptable BSR value.

	Returns
	-------
	valid : float
		The BSR value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the BSR value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	arg = arg_list(arg, 'BLAST Score Ratio')

	try:
		schema_bsr = float(arg)
		if schema_bsr >= min_value and schema_bsr <= max_value:
			valid = schema_bsr
		elif schema_bsr < min_value or schema_bsr > max_value:
			sys.exit(ct.INVALID_BSR)
	except Exception:
		sys.exit(ct.INVALID_BSR_TYPE.format(arg))

	return valid


def minimum_sequence_length_type(arg, min_value=ct.MSL_MIN, max_value=ct.MSL_MAX):
	"""Validate the minimum sequence length value (MSL) passed to chewBBACA.

	Parameters
	----------
	arg : int
		The MSL value passed to chewBBACA.
	min_value : int
		Minimum acceptable MSL value.
	max_value : int
		Maximum acceptable MSL value.

	Returns
	-------
	valid : int
		The MSL value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the MSL value cannot be converted to int type
		or if it is not contained in the acceptable interval.
	"""
	arg = arg_list(arg, 'minimum sequence length')

	try:
		schema_ml = int(arg)
		if schema_ml >= min_value and schema_ml <= max_value:
			valid = schema_ml
		elif schema_ml < min_value or schema_ml > max_value:
			sys.exit(ct.INVALID_MINLEN)
	except Exception:
		sys.exit(ct.INVALID_MINLEN_TYPE)

	return valid


def size_threshold_type(arg, min_value=ct.ST_MIN, max_value=ct.ST_MAX):
	"""Validate the size threshold value (ST) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The ST value passed to chewBBACA. Must be of type float
		or NoneType if no size threshold filter should be applied.
	min_value : float
		Minimum acceptable ST value.
	max_value : float
		Maximum acceptable ST value.

	Returns
	-------
	valid : float
		The ST value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the ST value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	arg = arg_list(arg, 'size threshold')

	try:
		schema_st = float(arg)
		if schema_st >= min_value and schema_st <= max_value:
			valid = schema_st
		elif schema_st < min_value or schema_st > max_value:
			sys.exit(ct.INVALID_ST)
	except Exception:
		if arg in [None, 'None']:
			valid = None
		else:
			sys.exit(ct.INVALID_ST_TYPE)

	return valid


def translation_table_type(arg, genetic_codes=ct.GENETIC_CODES):
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
	arg = arg_list(arg, 'translation table')

	try:
		schema_gen_code = int(arg)
		if schema_gen_code in genetic_codes:
			valid = schema_gen_code
		else:
			valid = False
	except Exception:
		valid = False

	if valid is False:
		# Format available genetic codes into list
		lines = ['\t{0}: {1}'.format(k, v) for k, v in genetic_codes.items()]
		gc_table = '\n{0}\n'.format('\n'.join(lines))

		sys.exit(ct.INVALID_GENETIC_CODE.format(gc_table))

	return valid


def validate_ws(arg, min_value=ct.WORD_SIZE_MIN, max_value=ct.WORD_SIZE_MAX):
	"""Validate the word size value (WS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The WS value passed to chewBBACA.
	min_value : float
		Minimum acceptable WS value.
	max_value : float
		Maximum acceptable WS value.

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
	arg = arg_list(arg, 'word size')

	try:
		if arg is None:
			valid = 'None'
		else:
			word_size = int(arg)
			if word_size >= min_value and word_size <= max_value:
				valid = word_size
			else:
				sys.exit(ct.INVALID_WS.format(min_value, max_value))
	except Exception:
		sys.exit(ct.INVALID_WS_TYPE)

	return valid


def validate_cs(arg, min_value=ct.CLUSTERING_SIMILARITY_MIN,
				max_value=ct.CLUSTERING_SIMILARITY_MAX):
	"""Validate the clustering similarity value (CS) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The CS value passed to chewBBACA. Must be of type float.
	min_value : float
		Minimum acceptable CS value.
	max_value : float
		Maximum acceptable CS value.

	Returns
	-------
	valid : float
		The CS value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the CS value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	arg = arg_list(arg, 'clustering similarity')

	try:
		if arg is None:
			valid = 'None'
		else:
			cluster_sim = float(arg)
			if cluster_sim >= min_value and cluster_sim <= max_value:
				valid = cluster_sim
			else:
				sys.exit(ct.INVALID_CS)
	except Exception:
		sys.exit(ct.INVALID_CS_TYPE)

	return valid


def validate_rf(arg, min_value=ct.REPRESENTATIVE_FILTER_MIN,
				max_value=ct.REPRESENTATIVE_FILTER_MAX):
	"""Validate the representative filter value (RF) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The RF value passed to chewBBACA. Must be of type float.
	min_value : float
		Minimum acceptable RF value.
	max_value : float
		Maximum acceptable RF value.

	Returns
	-------
	valid : float
		The RF value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the RF value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	arg = arg_list(arg, 'representative filter')

	try:
		if arg is None:
			valid = 'None'
		else:
			representative_filter = float(arg)
			if representative_filter >= min_value and representative_filter <= max_value:
				valid = representative_filter
			else:
				sys.exit(ct.INVALID_RF)
	except Exception:
		sys.exit(ct.INVALID_RF_TYPE)

	return valid


def validate_if(arg, min_value=ct.INTRA_CLUSTER_MIN,
				max_value=ct.INTRA_CLUSTER_MAX):
	"""Validate the intra-cluster filter value (IF) passed to chewBBACA.

	Parameters
	----------
	arg : float
		The IF value passed to chewBBACA. Must be of type float.
	min_value : float
		Minimum acceptable IF value.
	max_value : float
		Maximum acceptable IF value.

	Returns
	-------
	valid : float
		The IF value passed to chewBBACA, if it is valid.

	Raises
	------
	SystemExit
		- If the IF value cannot be converted to float type
		or if it is not contained in the acceptable interval.
	"""
	arg = arg_list(arg, 'intra-cluster filter')

	try:
		if arg is None:
			valid = 'None'
		else:
			intraCluster_filter = float(arg)
			if intraCluster_filter >= min_value and intraCluster_filter <= max_value:
				valid = intraCluster_filter
			else:
				sys.exit(ct.INVALID_ICF)
	except Exception:
		sys.exit(ct.INVALID_ICF_TYPE)

	return valid


def validate_ns_url(arg):
	"""Verify if the Chewie-NS URL passed to chewBBACA is valid.

	Parameters
	----------
	arg : str
		Identifier of the Chewie-NS instance or the URL
		to a instance of Chewie-NS.

	Returns
	-------
	ns_url : str
		URL to connect to the instance of Chewie-NS.

	Raises
	------
	SystemExit
		- If it is not possible to connect to the
		chewie-NS instance.
	"""
	if arg in ct.HOST_NS:
		ns_url = ct.HOST_NS[arg]
	else:
		ns_url = arg

	# sync schema has None by default to get ns_url in schema URI
	if ns_url is not None:
		# check if server is up
		conn = cr.check_connection(ns_url)
		if conn is False:
			sys.exit(ct.NS_CANNOT_CONNECT.format(ns_url))

	return ns_url


def validate_python_version(minimum_version=ct.MIN_PYTHON):
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
	python_version = platform.python_version()

	try:
		assert tuple(map(int, python_version.split('.'))) >= minimum_version[0]
	except AssertionError:
		sys.exit(ct.PYTHON_VERSION.formta(python_version, minimum_version[1]))

	return python_version


def verify_cpu_usage(cpu_to_use):
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
	total_cpu = multiprocessing.cpu_count()

	cpu_to_use = int(cpu_to_use)

	# Do not allow a value greater than the number of cores
	if cpu_to_use >= total_cpu:
		# Define a value that is safe according to the number of
		# available cores/threads
		if total_cpu > 2:
			cpu_to_use = total_cpu - 2
		elif total_cpu == 2:
			cpu_to_use = 1
		print(ct.CPU_RESET_WARNING.format(cpu_to_use))
	elif cpu_to_use == (total_cpu - 1):
		print(ct.CPU_VALUE_WARNING.format(cpu_to_use, total_cpu))

	return cpu_to_use


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


def get_blast_path(blast_path):
	"""Determines if BLAST is in PATH.

	Parameters
	----------
	blast_path : str or NoneType
		Path to the directory with the BLAST executables or
		NoneType if user did not provide a value.

	Returns
	-------
	blast_path : str or NoneType
		Validated path to the directory that contains the BLAST
		executables or NoneType if it was not possible to validate
		the provided path or if the BLAST executables are not in PATH.
	"""
	# Search for BLAST in PATH
	if not blast_path:
		blastp_path = shutil.which(ct.BLASTP_ALIAS)
		if blastp_path is not None:
			blast_path = os.path.dirname(blastp_path)
		else:
			blast_path = None
	# Validate user-provided path
	else:
		blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
		executable = is_exe(blastp_path)
		if executable is False:
			blast_path = None

	return blast_path


def get_blast_version(blast_path):
	"""Determines BLAST version.

	Parameters
	----------
	blast_path : str
		Path to the directory that contains the BLAST executables.

	Returns
	-------
	version : dict or NoneType
		Dictionary with the BLAST MAJOR and MINOR versions or
		NoneType if it was not possible to determine the BLAST
		version.
	"""
	blastp_path = fo.join_paths(blast_path, [ct.BLASTP_ALIAS])
	# Check BLAST version
	try:
		proc = subprocess.Popen([blastp_path, '-version'],
								stdout=subprocess.PIPE,
								stderr=subprocess.PIPE)
		stdout, stderr = proc.communicate()
		version_string = stdout.decode('utf8')
		version_pattern = r'^blastp:\s(?P<MAJOR>\d+).(?P<MINOR>\d+).(?P<REV>\d+).*'
		blast_version_pat = re.compile(version_pattern)

		match = blast_version_pat.search(version_string)
		if match is not None:
			version = {k: int(v) for k, v in match.groupdict().items()}
		else:
			version = None
	except:
		version = None

	return version


def check_blast(blast_path, major=ct.BLAST_MAJOR, minor=ct.BLAST_MINOR):
	"""Determine if BLAST is installed and validates its version.

	Parameters
	----------
	blast_path : str or Nonetype
		Path to the directory with BLAST executables or
		NoneType if user did not provide a value.
	major : int
		BLAST minimun MAJOR version.
	minor : int
		BLAST minimum MINOR version.

	Returns
	-------
	blast_path : str
		Path to the directory with BLAST executables.

	Raises
	------
	SystemExit
		- If the user did not provide a value and BLAST is
		not in PATH.
		- If the user provided a value but that path does not
		contain BLAST executables.
		- If it is not possible to determine the BLAST
		version or if it does not match minimum requirements.
	"""
	# Validate BLAST path
	blast_path = get_blast_path(blast_path)
	# Exit if it is not possible to get BLAST path
	if blast_path is None:
		sys.exit(ct.BLAST_NO_PATH)
	# Get BLAST version
	blast_version = get_blast_version(blast_path)
	# Exit if it is not possible to get BLAST version
	if blast_version is None:
		sys.exit(ct.BLAST_NO_VERSION.format(major, minor))
	if blast_version['MAJOR'] < major or (blast_version['MAJOR'] >= major and blast_version['MINOR'] < minor):
		sys.exit(ct.BLAST_UPDATE.format(blast_version['MAJOR'],
										blast_version['MINOR'],
				 						major, minor))

	return blast_path


def hash_ptf(ptf_path):
	"""Determine hash value for a Prodigal training file.

	Parameters
	----------
	ptf_path : str or None
		Path to the Prodigal training file or None if no
		training file will be used.

	Returns
	-------
	ptf_hash : str
		Blake2b hash computed from file content.
	"""
	if ptf_path is not None:
		ptf_hash = fo.hash_file(ptf_path, hashlib.blake2b())
	else:
		ptf_hash = None

	return ptf_hash


def validate_ptf_path(ptf_path, schema_directory):
	""" Determines if the path to the Prodigal training file
		is valid. Gets the training file in the schema's
		directory if the input path is of type NoneType.

	Parameters
	----------
	ptf_path : str or NoneType
		Path to the Prodigal training file or NoneType
		if no value was provided.
	schema_directory : str
		Path to the schema's directory.

	Returns
	-------
	ptf_path : str or bool
		Path to the Prodigal training file or False if
		no training file should be used.

	Raises
	------
	SystemExit
		- If there is more than one training file in
		the schema's directory.
		- If a path was provided and it is not valid.
	"""
	if ptf_path is None:
		# Deal with multiple training files
		schema_ptfs = [file
					   for file in os.listdir(schema_directory)
					   if file.endswith('.trn')]
		if len(schema_ptfs) > 1:
			sys.exit(ct.MULTIPLE_PTFS)
		elif len(schema_ptfs) == 1:
			if schema_ptfs[0] is not None:
				ptf_path = os.path.join(schema_directory, schema_ptfs[0])
			else:
				print(ct.MISSING_PTF)
				ptf_path = None
	else:
		if os.path.isfile(ptf_path) is False:
			message = (ct.INVALID_PTF_PATH)
			sys.exit(message)

	return ptf_path


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
		ptf_hash = hash_ptf(ptf_path)
	else:
		ptf_hash = None

	unmatch = validate_ptf_hash(ptf_hash, schema_ptfs, force_continue)

	return [ptf_path, ptf_hash, unmatch]


def solve_conflicting_arguments(schema_params, ptf_path, blast_score_ratio,
								translation_table, minimum_length,
								size_threshold, force_continue, config_file,
								schema_directory):
	""" Compares schema parameters values stored in the config
		file with values provided by the user to solve conflicting
		cases. Adds/appends new values to the config file if the
		user wants to use values that do not match schema's
		default values.

	Parameters
	----------
	schema_params : dict
		Dictionary with the schema's config values.
	ptf_path : str or NoneType
		Path to the Prodigal training file or NoneType
		if no value was passed through the command line.
	blast_score_ratio : float or NoneType
		BLAST Score Ratio value. NoneType if no value was
		passed.
	translation_table : int
		Translation table value. NoneType if no value was
		passed.
	minimum_length : int
		Minimum sequence length value. NoneType if no value was
		passed.
	size_threshold : float
		Allele size variation threshold. NoneType if no value was
		passed.
	force_continue : bool
		True to validate parameters values without prompting users.
		False otherwise.
	config_file : str
		Path to the schema's configuration file.
	schema_directory : str
		Path to the schema's directory.

	Returns
	-------
	run_params : dict
		Dictionary with the arguments validated values that
		will be used for allele calling.
	"""
	# Parameter values for current run
	run_params = {'bsr': blast_score_ratio,
				  'minimum_locus_length': minimum_length,
				  'size_threshold': size_threshold}

	# Determine user provided values that differ from default
	unmatch_params = {k: v
					  for k, v in run_params.items()
					  if v not in schema_params[k] and v is not None}

	# Update run values equal to None
	for k, v in run_params.items():
		if v is None:
			run_params[k] = schema_params[k][0]

	if len(unmatch_params) > 0:
		print(ct.ARGS_DIFFER)
		params_diffs = [[p, ':'.join(map(str, schema_params[p])),
						 str(unmatch_params[p])]
						for p in unmatch_params]
		params_diffs_text = ['{:^20} {:^20} {:^10}'.format('Argument', 'Schema', 'Provided')]
		params_diffs_text += ['{:^20} {:^20} {:^10}'.format(p[0], p[1], p[2]) for p in params_diffs]
		print('\n'.join(params_diffs_text))
		if force_continue is False:
			params_answer = fo.input_timeout(ct.ARGS_DIFFER_PROMPT, ct.PROMPT_TIMEOUT)
		else:
			params_answer = 'yes'

		if params_answer.lower() not in ['y', 'yes']:
			sys.exit('Exited.')
		else:
			# Append new argument values to config values
			for p in unmatch_params:
				schema_params[p].append(unmatch_params[p])

	# Default is to get the training file in schema directory
	schema_ptfs = schema_params['prodigal_training_file']
	ptf_path, ptf_hash, unmatch = validate_ptf(ptf_path, schema_directory,
											   schema_ptfs, force_continue)

	run_params['ptf_path'] = ptf_path
	if unmatch is True:
		schema_params['prodigal_training_file'].append(ptf_hash)
		unmatch_params['prodigal_training_file'] = ptf_hash

	# Update translation table
	if ptf_path:
		# Get translation table used to create training file
		ptf_table = gp.read_training_file(ptf_path).translation_table
		run_params['translation_table'] = ptf_table
		if ptf_table not in schema_params['translation_table']:
			schema_params['translation_table'].append(ptf_table)
			unmatch_params['translation_table'] = ptf_table
	else:
		if not translation_table:
			run_params['translation_table'] = schema_params['translation_table'][0]
		else:
			run_params['translation_table'] = translation_table
			if translation_table not in schema_params['translation_table']:
				schema_params['translation_table'].append(translation_table)
				unmatch_params['translation_table'] = translation_table

	# Update schema config file
	if len(unmatch_params) > 0:
		fo.pickle_dumper(schema_params, config_file)

	return run_params


def write_gene_list(schema_dir):
	"""Save list of loci in a schema to the '.genes_list' file.

	Parameters
	----------
	schema_dir : str
		Path to the schema directory.

	Returns
	-------
	A list with two elements. A boolean value that
	is True if the file with the list of genes was
	created, False otherwise. The second element
	is the path to the created file.
	"""
	# Loci FASTA files must end with '.fasta' extension
	schema_files = os.listdir(schema_dir)
	loci_files, _ = fo.filter_by_extension(schema_files, ['.fasta'])
	output_file = fo.join_paths(schema_dir, [ct.GENE_LIST_BASENAME])
	fo.pickle_dumper(loci_files, output_file)

	return [os.path.isfile(output_file), output_file]


def write_schema_config(args, chewie_version, output_directory):
	""" Writes chewBBACA's parameter values used to create
		a schema to a file.

	Parameters
	----------
	args : dict
		Dictionary with the parameter values to store in the
		schema config file.
	chewie_version : str
		Version of the chewBBACA suite used to create
		the schema.
	output_directory : str
		Path to the output directory where the file with
		schema parameters values will be created.

	Returns
	-------
	A list with two elements. A boolean value that
	is True if the file with the parameters values was
	created and False otherwise. The second element
	is the path to the created file.
	"""

	# Deal with multiple names for the same parameter
	size_threshold = None if args['size_threshold'] in [None, 'None'] else float(args['size_threshold'])
	bsr = float(args.get('blast_score_ratio')) if 'blast_score_ratio' in args else float(args['bsr'])
	minimum_locus_length = int(args.get('minimum_length')) if 'minimum_length' in args else int(args['minimum_locus_length'])
	cluster_sim = args.get('clustering_sim') if 'clustering_sim' in args else args['cluster_sim']
	intraCluster_filter = args.get('intra_filter') if 'intra_filter' in args else args['intraCluster_filter']

	params = {}
	params['bsr'] = [bsr]
	params['prodigal_training_file'] = [args['ptf_path']]
	params['translation_table'] = [int(args['translation_table'])]
	params['minimum_locus_length'] = [minimum_locus_length]
	params['chewBBACA_version'] = [chewie_version]
	params['size_threshold'] = [size_threshold]
	params['word_size'] = [args['word_size']]
	params['window_size'] = [args['window_size']]
	params['cluster_sim'] = [cluster_sim]
	params['representative_filter'] = [args['representative_filter']]
	params['intraCluster_filter'] = [intraCluster_filter]

	config_file = os.path.join(output_directory, ct.SCHEMA_CONFIG_BASENAME)
	fo.pickle_dumper(params, config_file)

	return [os.path.isfile(config_file), config_file]


def read_configs(schema_path, filename):
	""" Reads file with schema config values.

	Parameters
	----------
	schema_path : str
		Path to the schema's directory.
	filename : str
		Name of the file that contains the config values.

	Returns
	-------
	configs : dict
		Dictionary with config names as keys and config
		values as values.
	"""

	config_file = os.path.join(schema_path, filename)
	if os.path.isfile(config_file):
		# Load configs dictionary
		configs = fo.pickle_loader(config_file)
	else:
		sys.exit(ct.MISSING_CONFIG)

	return configs


def check_input_type(input_path, output_file):
	"""Validate input and create list of files to use.

	Parameters
	----------
	input_path : str
		Path to a file or directory.
	output_file : str
		Path to the output file created to store the paths
		to valid FASTA files.

	Returns
	-------
	output_file : str
		Path to the output file created to store the paths
		to valid FASTA files.

	Raises
	------
	SystemExit
		- If the input path is not a valid path for a file or
		  directory.
	"""
	# Input path is for a file
	if os.path.isfile(input_path):
		output_file, total_inputs = validate_input_file(input_path, output_file)
	# Input path is for a directory
	elif os.path.isdir(input_path):
		output_file, total_inputs = validate_input_dir(input_path, output_file)
	else:
		sys.exit(ct.INVALID_INPUT_PATH)

	return output_file, total_inputs


def validate_input_file(input_path, output_file):
	"""Validate a file with a list of paths to input files.

	Parameters
	----------
	input_path : str
		Path to a file with a list of paths.
	output_file : str
		Path to the output file created to store the paths
		to valid FASTA files.

	Returns
	-------
	output_file : str
		Path to the output file created to store the paths
		to valid FASTA files.

	Raises
	------
	SystemExit
		- If the input path is for a FASTA file.
		- If any of the provided paths does not exist.
		- If any of the file basenames does not end with one of
		  the accepted file extensions.
		- If the format of any of the files is not FASTA.
	"""
	# Check if it is a single FASTA file
	if fao.validate_fasta(input_path) is True:
		# Exit if input is a single FASTA file
		sys.exit(ct.FASTA_INPUT_EXCEPTION)

	# Read list of input files
	files = [line[0] for line in fo.read_tabular(input_path)]

	invalid_files = []
	# Need to verify if files end with any of the accepted file
	# extensions, not only '.fasta'
	valid_extension, invalid_extension = fo.filter_by_extension(files, ct.FASTA_EXTENSIONS)
	if len(invalid_extension) > 0:
		invalid_files.append([invalid_extension, ct.INVALID_EXTENSION_EXCEPTION])

	# Check that all files exist
	missing = [file for file in files if os.path.exists(file) is False]
	if len(missing) > 0:
		invalid_files.append([missing, ct.MISSING_INPUTS_EXCEPTION])

	# Only keep files whose content is typical of a FASTA file
	fasta_files, non_fasta = fao.filter_non_fasta(files)
	if len(non_fasta) > 0:
		invalid_files.append([non_fasta, ct.NON_FASTA_EXCEPTION])

	# Exit if list of input files contained invalid files
	if len(invalid_files) > 0:
		exception_messages = [e[1].format(im.join_list(e[0], '\n')) for e in invalid_files]
		sys.exit(im.join_list(exception_messages, '\n'))
	# Save file paths to output file
	else:
		fo.write_lines(files, output_file)

	return output_file, len(files)


def validate_input_dir(input_path, output_file):
	"""List and validate input files in a directory.

	Parameters
	----------
	input_path : str
		Path to the directory that contains the input files.
	output_file : str
		Path to the output file created to store the paths
		to valid FASTA files.

	Returns
	-------
	output_file : str
		Path to the output file created to store the paths
		to valid FASTA files.

	Raises
	------
	SystemExit
		- If there are no valid FASTA files in the input directory.
	"""
	# List absolute paths
	# Only keep paths to files
	files = [file for file in fo.listdir_fullpath(input_path)
			 if os.path.isdir(file) is False]

	# Filter based on file extension
	valid_extension, invalid_extension = fo.filter_by_extension(files, ct.FASTA_EXTENSIONS)

	# Only keep files whose content is typical of a FASTA file
	fasta_files, non_fasta = fao.filter_non_fasta(valid_extension)

	# If there are FASTA files
	if len(fasta_files) > 0:
		# Save file paths to output file
		fo.write_lines(fasta_files, output_file)
	else:
		sys.exit(ct.MISSING_FASTAS_EXCEPTION)

	return output_file, len(files)


def validate_loci_list(input_path, output_file, parent_dir=None):
	"""Validate a list of paths to loci FASTA files or loci IDs.

	Parameters
	----------
	input_path : str
		Path to a file with a list of paths.
	output_file : str
		Path to the output file created to store the paths
		to valid loci FASTA files.
	parent_dir : str
		Parent directory to add to construct paths
		to input files when users provide a file
		with file names.

	Returns
	-------
	Path to the output file created to store the paths
	to valid loci FASTA files.

	Raises
	------
	SystemExit
		- If the input path is for a FASTA file.
		- If any of the provided paths does not exist.
		- If the format of any of the files is not FASTA.
	"""
	# Check if it is a single FASTA file
	if fao.validate_fasta(input_path) is True:
		# Exit if input is a single FASTA file
		sys.exit(ct.FASTA_LOCI_LIST_EXCEPTION)

	# Read list of input files
	files = [line[0] for line in fo.read_tabular(input_path)]

	# List must have full paths
	if parent_dir is not None:
		# Add parent directory path if necessary
		files = [os.path.join(parent_dir, fo.file_basename(file))
				 if parent_dir not in file
				 else file
				 for file in files]

	# Add '.fasta' extension if it is missing from IDs
	# Loci files use the '.fasta' extension, anything else might
	# mean there is an issue with the schema or it is an external schema
	files = [file+'.fasta'
			 if any([file.endswith(ext) for ext in ct.FASTA_EXTENSIONS]) is False
			 else file
			 for file in files]

	# Check that all files exist
	invalid_files = []
	missing = [file for file in files if os.path.exists(file) is False]
	if len(missing) > 0:
		invalid_files.append([missing, ct.MISSING_LOCI_EXCEPTION])

	# Only keep files whose content is typical of a FASTA file
	fasta_files, non_fasta = fao.filter_non_fasta(files)
	if len(non_fasta) > 0:
		invalid_files.append([non_fasta, ct.NON_FASTA_LOCI_EXCEPTION])

	# Exit if list of input files contained invalid files
	if len(invalid_files) > 0:
		exception_messages = [e[1].format(im.join_list(e[0], '\n')) for e in invalid_files]
		sys.exit(im.join_list(exception_messages, '\n'))
	# Save file paths to output file
	else:
		fo.write_lines(files, output_file)

	return output_file


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
	basenames = [fo.file_basename(file) for file in path_list]
	prefixes = {}
	for name in basenames:
		prefix = fo.split_joiner(name, [0], '.')
		prefixes.setdefault(prefix, []).append(name)

	return prefixes


def check_unique_prefixes(input_list):
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
	input_paths = fo.read_lines(input_list)
	prefixes = get_file_prefixes(input_paths)

	# Detect if some inputs share the same unique prefix
	if len(set(prefixes)) < len(input_paths):
		repeated_basenames = [v for k, v in prefixes.items() if len(v) > 1]
		repeated_basenames = [','.join(l) for l in repeated_basenames]
		sys.exit(ct.INPUTS_SHARE_PREFIX.format('\n'.join(repeated_basenames)))

	return False


def check_blanks(input_list):
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
	input_paths = fo.read_lines(input_list)
	basenames = [fo.file_basename(file) for file in input_paths]
	include_blanks = [name for name in basenames if ' ' in name]

	if len(include_blanks) > 0:
		sys.exit(ct.INPUTS_INCLUDE_BLANKS.format('\n'.join(include_blanks)))

	return False


def check_prefix_length(input_list):
	"""Check if input file prefixes are not longer than 30 characters.

	Parameters
	----------
	input_list : str
		Path to file that contains the list of paths to input files.

	Returns
	-------
	False if there are no input files with a unique identifier with
	more than 30 characters.

	Raises
	------
	SystemExit
		- If any input file has a unique identifier with more than 30 characters.
	"""
	input_paths = fo.read_lines(input_list)
	prefixes = get_file_prefixes(input_paths)
	long_prefixes = {k: v for k, v in prefixes.items() if len(k) > ct.PREFIX_MAXLEN}
	# Exit if any file prefix is longer than 30 characters
	if len(long_prefixes) > 0:
		long_prefixes_msg = [f'{v[0]} ({k}, {len(k)} chars)' for k, v in long_prefixes.items()]
		long_prefixes_msg = '\n'.join(long_prefixes_msg)
		sys.exit(ct.INPUTS_LONG_PREFIX.format(long_prefixes_msg))

	return False


def check_prefix_pdb(input_list, output_directory, makeblastdb_path, blastdbcmd_path):
	"""Check if the BLAST database includes the expected sequence IDs.

	Parameters
	----------
	input_list : str
		Path to file that contains the list of paths to input files.
	output_directory : str
		Path to the directory where dummy data will be created.
	makeblastdb_path : str
		Path to the makeblastdb executable.
	blastdbcmd_path : str
		Path to the blastdbcmd executable.

	Returns
	-------
	None
	"""
	input_paths = fo.read_lines(input_list)
	prefixes = get_file_prefixes(input_paths)
	# Create directory to store dummy data
	dummy_dir = os.path.join(output_directory, ct.DUMMY_DIR)
	fo.create_directory(dummy_dir)
	# Create dummy FASTA records
	dummy_seqids = [f'{i}-protein1' for i in (prefixes)]
	dummy_records = [fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [i, ct.DUMMY_PROT]) for i in dummy_seqids]
	dummy_fasta = os.path.join(dummy_dir, ct.DUMMY_FASTA)
	fo.write_lines(dummy_records, dummy_fasta)
	# Create BLAST db
	dummy_blastdb = os.path.join(dummy_dir, ct.DUMMY_BLASTDB)
	bw.make_blast_db(makeblastdb_path, dummy_fasta, dummy_blastdb, 'prot')
	# Get sequence from BLAST db
	dummy_blastdbcmd_fasta = os.path.join(dummy_dir, ct.DUMMY_BLASTDBCMD_FASTA)
	blastdbcmd_std = bw.run_blastdbcmd(blastdbcmd_path, dummy_blastdb, dummy_blastdbcmd_fasta)
	# Check if the sequence IDs in the BLAST db match the expected IDs
	records = fao.sequence_generator(dummy_blastdbcmd_fasta)
	recids = [record.id for record in records]
	modified_seqids = set(dummy_seqids) - set(recids)
	# Delete dummy data
	fo.delete_directory(dummy_dir)
	# Exit if BLAST db includes sequence IDs that do not match the expected IDs
	if len(modified_seqids) > 0:
		pdb_prefix = [prefixes[p.split('-protein')[0]][0] for p in modified_seqids]
		sys.exit(ct.INPUTS_PDB_PREFIX.format('\n'.join(pdb_prefix)))
