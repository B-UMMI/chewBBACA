#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions/classes related to parameter
and argument validation.

Code documentation
------------------
"""


import os
import re
import sys
import csv
import zlib
import shutil
import hashlib
import logging
import argparse
import platform
import subprocess

try:
    from utils import (constants as ct,
                       file_operations as fo,
                       chewiens_requests as cr,
                       fasta_operations as fao)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 chewiens_requests as cr,
                                 fasta_operations as fao)


logger = logging.getLogger('PV')


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

    # prog is the name of the program 'ex: chewBBACA.py'
    def __init__(self, prog, indent_increment=2, max_help_position=56, width=100):
        super().__init__(prog, indent_increment, max_help_position, width)

    # override split lines method
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

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return '%s %s' % (', '.join(parts), args_string)

            return ', '.join(parts)

    def _get_default_metavar_for_optional(self, action):
        return action.dest.upper()


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
            logger.error('BSR value is not contained in the '
                         f'{[min_value, max_value]} interval.')
            sys.exit(1)
    except Exception:
        logger.error(f'Invalid BSR value: {arg}. BSR value must be '
                     f'contained in the {[min_value, max_value]} interval.')
        sys.exit(1)

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
            logger.error(f'Invalid minimum sequence length value: {arg}. '
                         f'Must be equal or greater than {min_value}.')
            sys.exit(1)
    except Exception:
        logger.error(f'Invalid minimum sequence length value: {arg}. '
                     f'Must be an integer equal or greater than {min_value}.')
        sys.exit(1)

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
            logger.error(f'Invalid size threshold value: {arg}. Must be None '
                         f'or be contained in the {[min_value, max_value]} '
                         'interval.')
            sys.exit(1)
    except Exception:
        if arg in [None, 'None']:
            valid = None
        else:
            logger.error(f'Invalid size threshold value: {arg}. Must be None '
                         f'or be contained in the {[min_value, max_value]} '
                         'interval.')
            sys.exit(1)

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
        # format available genetic codes into list
        lines = ['{0}: {1}'.format(k, v) for k, v in genetic_codes.items()]
        gc_table = '{0}\n'.format('\n'.join(lines))
        logger.error(f'Invalid genetic code value: {arg}. Value must '
                     'correspond to one of the accepted genetic codes. '
                     f'Accepted genetic codes:\n{gc_table}')
        sys.exit(1)

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
                sys.exit('\nWord size for the clustering step '
                         'must be equal or greater than {0} and '
                         'equal or smaller than {1}.'.format(min_value, max_value))
    except Exception:
        sys.exit('\nSchema created with invalid clustering word '
                 'size value.')

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
                sys.exit('\nClustering similarity threshold value '
                         'must be contained in the [0.0, 1.0] '
                         'interval.')
    except Exception:
        sys.exit('\nSchema created with invalid clustering '
                 'threshold value.')

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
                sys.exit('\nRepresentative filter threshold value '
                         'must be contained in the [0.0, 1.0] '
                         'interval.')
    except Exception:
        sys.exit('\nSchema created with invalid representative filter value.')

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
                sys.exit('\nIntra-cluster filter value '
                         'must be contained in the [0.0, 1.0] '
                         'interval.')
    except Exception:
        sys.exit('\nSchema created with invalid intra-cluster filter '
                 'value.')

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
            sys.exit('Failed to establish a connection to the Chewie-NS '
                     'instance at {0}.'.format(ns_url))

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
        print('Python version found: {} '.format(python_version))
        print('Please use version Python >= {0}'.format(minimum_version[1]))
        sys.exit(0)

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
    # get number of logical CPU cores
    total_cpu = os.cpu_count()
    cpu_to_use = int(cpu_to_use)

    # do not allow a value greater than the number of cores
    if cpu_to_use >= total_cpu:
        logger.warning(f'CPU core count value provided ({cpu_to_use}) is '
                       'equal to or exceeds the number of available CPU '
                       f'cores ({total_cpu}).')
        # define a value that is safe according to the number of
        # available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1
        logger.warning(f'Resetting CPU core count to: {cpu_to_use}')
    elif cpu_to_use == (total_cpu - 1):
        logger.warning(f'CPU core count value provided ({cpu_to_use}) is '
                       'close to the number of available CPU cores '
                       f'({total_cpu}). This may affect system '
                       'responsiveness.')

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


def which(program):
    """Determine if a program is in PATH.

    Parameters
    ----------
    program : str
        The name used to call the program.

    Returns
    -------
    program_path : str
        Program path added to PATH.
    """
    proc = subprocess.Popen(['which', program],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    stdout, stderr = proc.communicate()

    program_path = stdout.decode('utf8').strip()

    return program_path


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
    # search for BLAST in PATH
    if blast_path == '':
        blastp_path = which(ct.BLASTP_ALIAS)
        if not blastp_path:
            logger.error('Could not find BLAST executables in PATH.')
            sys.exit(1)
        else:
            blast_path = os.path.dirname(blastp_path)
    # validate user-provided path
    else:
        blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
        executable = is_exe(blastp_path)
        if executable is False:
            logger.error(f'Could not find BLAST executables in {blast_path}')
            sys.exit(1)

    # check BLAST version
    try:
        proc = subprocess.Popen([blastp_path, '-version'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except Exception as e:
        logger.error('Could not determine BLAST version. Please verify '
                     'that BLAST is installed and added to PATH.')
        sys.exit(1)

    stdout, stderr = proc.communicate()

    version_string = stdout.decode('utf8')
    version_pattern = r'^blastp:\s(?P<MAJOR>\d+).(?P<MINOR>\d+).(?P<REV>\d+).*'
    blast_version_pat = re.compile(version_pattern)

    match = blast_version_pat.search(version_string)
    if match is None:
        logger.error('Could not determine BLAST version. Please verify '
                     'that BLAST is installed and added to PATH.')
        sys.exit(1)

    version = {k: v for k, v in match.groupdict().items()}
    if int(version['MAJOR']) < major or (int(version['MAJOR']) >= major and int(version['MINOR']) < minor):
        logger.error(f'Please update BLAST to version >= {major}.{minor}.0')
        sys.exit(1)

    version_string = ".".join(list(version.values()))
    logger.info(f'BLAST version: {version_string}')

    return blast_path


def check_prodigal(prodigal_path):
    """Determine if Prodigal is installed and added to the PATH.

    Parameters
    ----------
    prodigal_path : str
        Name or path used to call Prodigal.

    Returns
    -------
    True if Prodigal is installed.

    Raises
    ------
    SystemExit
        - If it is not possible to run the command to determine
        the Prodigal version.
    """
    # check Prodigal version
    try:
        proc = subprocess.Popen([prodigal_path, '-v'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except Exception as e:
        logger.error('Could not determine Prodigal version '
                     f'({prodigal_path} -v). Please verify that Prodigal '
                     'is installed and added to PATH.')
        # automatically captures and logs traceback
        logger.exception(e)
        sys.exit(1)

    # prodigal version is captured in stderr
    stdout, stderr = proc.communicate()

    return stderr.decode().strip()


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


def prompt_arguments(ptf_path, blast_score_ratio, translation_table,
                     minimum_length, size_threshold, version):
    """ Detects if a valid argument value was passed and asks
        for a value if the value of any argument is of type
        NoneType. Raises SystemExit if a provided value is not
        valid.
    """
    prompt = ('It seems that your schema was created with an older '
              'version or missing some files.\nIt is highly recommended '
              'that you run the PrepExternalSchema process to guarantee full '
              'compatibility with the new chewBBACA version.\nIf you '
              'wish to continue, the AlleleCall process will convert '
              'the schema to v{0}, but will not determine if schema '
              'structure respects configuration values.\nDo you wish '
              'to proceed?\n'.format(version))
    proceed = fo.input_timeout(prompt, ct.prompt_timeout)
    if proceed.lower() not in ['y', 'yes']:
        logger.info(f'Input: {proceed}. Exited.')
        sys.exit(1)

    # determine parameters default values to include in config file
    if ptf_path is None:
        prompt = ('Full path to the Prodigal training file:\n')
        ptf_path = fo.input_timeout(prompt, ct.prompt_timeout)
        if ptf_path in ['', 'None']:
            ptf_path = None
        else:
            if os.path.isfile(ptf_path) is False:
                logger.error('Provided path does not exist.')
                sys.exit(1)
    else:
        if os.path.isfile(ptf_path) is False:
            logger.error('Provided path for Prodigal training file does '
                         'not exist.')
            sys.exit(1)

    if blast_score_ratio is None:
        prompt = (f'BLAST score ratio value ([{ct.BSR_MIN}, '
                  f'{ct.BSR_MAX}]):\n')
        blast_score_ratio = fo.input_timeout(prompt, ct.prompt_timeout)

    blast_score_ratio = bsr_type(blast_score_ratio)

    if translation_table is None:
        table_ids = ','.join(map(str, list(ct.GENETIC_CODES.keys())))
        prompt = (f'Translation table value ({{{table_ids}}}):\n')
        translation_table = fo.input_timeout(prompt, ct.prompt_timeout)

    translation_table = translation_table_type(translation_table)

    if minimum_length is None:
        prompt = (f'Minimum length value (>={ct.MSL_MIN}):\n')
        minimum_length = fo.input_timeout(prompt, ct.prompt_timeout)

    minimum_length = minimum_sequence_length_type(minimum_length)

    if size_threshold is None:
        prompt = (f'Size threshold value ([{ct.ST_MIN}, {ct.ST_MAX}]):\n')
        size_threshold = fo.input_timeout(prompt, ct.prompt_timeout)

    size_threshold = size_threshold_type(size_threshold)

    return [ptf_path, blast_score_ratio, translation_table,
            minimum_length, size_threshold]


def auto_arguments(ptf_path, blast_score_ratio, translation_table,
                   minimum_length, size_threshold):
    """ Detects if a valid argument value was passed and
        selects default config values if the provided value
        is of type NoneType.
    """
    if ptf_path is not None:
        if os.path.isfile(ptf_path) is False:
            logger.error('Provided path to Prodigal training file '
                         'does not exist.')
            sys.exit(1)
    else:
        ptf_path = None

    blast_score_ratio = (blast_score_ratio
                         if blast_score_ratio is not None
                         else ct.DEFAULT_BSR)
    blast_score_ratio = bsr_type(blast_score_ratio)

    translation_table = (translation_table
                         if translation_table is not None
                         else ct.GENETIC_CODES_DEFAULT)
    translation_table = translation_table_type(translation_table)

    minimum_length = (minimum_length
                      if minimum_length is not None
                      else ct.MSL_MIN)
    minimum_length = minimum_sequence_length_type(minimum_length)

    size_threshold = (size_threshold
                      if size_threshold is not None
                      else ct.SIZE_THRESHOLD_DEFAULT)
    size_threshold = size_threshold_type(size_threshold)

    return [ptf_path, blast_score_ratio, translation_table,
            minimum_length, size_threshold]


def upgrade_legacy_schema(ptf_path, schema_directory, blast_score_ratio,
                          translation_table, minimum_length, version,
                          size_threshold, force_continue):
    """Upgrade a legacy schema to current version.

    Determines configuration values and adds them to the
    configuration file. Creates the file with the list of loci in
    the schema.

    Parameters
    ----------
    ptf_path : str or None
        Path to the Prodigal training file or NoneType
        if no value was provided.
    schema_directory : str
        Path to the schema's directory.
    blast_score_ratio : float
        BLAST Score Ratio value.
    translation_table : int
        Translation table integer identifier.
    minimum_length : int
        Minimum sequence length value.
    version : str
        chewBBACA's version.
    size_threshold : float
        Allele size variation threshold.
    force_continue : bool
        True if argument values should be automatically
        validated and set to default values if their type
        is NoneType. False to validate provided values and
        ask for a value if the user did not provide one.

    Returns
    -------
    Valid/selected values for ptf_path, blast_score_ratio,
    translation_table, minimum_length and size_threshold.
    """
    if force_continue is False:
        values = prompt_arguments(ptf_path, blast_score_ratio,
                                  translation_table, minimum_length,
                                  size_threshold, version)
    # forced method, create with default args defined in constants
    else:
        values = auto_arguments(ptf_path, blast_score_ratio,
                                translation_table, minimum_length,
                                size_threshold)

    ptf_path, blast_score_ratio,\
        translation_table, minimum_length, size_threshold = values

    # copy training file to schema directory
    if ptf_path is not None:
        fo.copy_file(ptf_path, schema_directory)
        logger.info(f'Copied {ptf_path} to {schema_directory}')

    # determine PTF hash
    ptf_hash = hash_ptf(ptf_path)

    # write schema config file
    schema_config = write_schema_config(blast_score_ratio, ptf_hash,
                                        translation_table, minimum_length,
                                        version, size_threshold,
                                        ct.WORD_SIZE_DEFAULT,
                                        ct.WINDOW_SIZE_DEFAULT,
                                        ct.CLUSTERING_SIMILARITY_DEFAULT,
                                        ct.REPRESENTATIVE_FILTER_DEFAULT,
                                        ct.INTRA_CLUSTER_DEFAULT,
                                        schema_directory)

    # create hidden file with genes/loci list
    genes_list_file = write_gene_list(schema_directory)

    return [ptf_path, blast_score_ratio, translation_table,
            minimum_length, size_threshold]


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
        # deal with multiple training files
        schema_ptfs = [file
                       for file in os.listdir(schema_directory)
                       if file.endswith('.trn')]
        if len(schema_ptfs) > 1:
            logger.error('Found more than one Prodigal training '
                         'file in the schema directory.\nPlease maintain '
                         'only the training file used in the schema '
                         'creation process.')
            sys.exit(1)
        elif len(schema_ptfs) == 1:
            if schema_ptfs[0] is not None:
                ptf_path = os.path.join(schema_directory, schema_ptfs[0])
            else:
                logger.warning('There is no Prodigal training file in schema '
                               'directory.')
                ptf_path = None
    else:
        if os.path.isfile(ptf_path) is False:
            logger.error('Cannot find specified Prodigal training file.'
                         '\nPlease provide a valid training file.\n\nYou '
                         'can create a training file for a species of '
                         'interest with the following command:\n  prodigal '
                         '-i <reference_genome> -t <training_file.trn> -p '
                         'single\n\nIt is strongly advised to provide a '
                         'high-quality and closed genome for the training '
                         'process.')
            sys.exit(1)

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
                logger.warning('Prodigal training file is not the one '
                               'used to create the schema.')
                prompt = ('Using this training file might lead to '
                          'results not consistent with previous runs '
                          'and invalidate the schema for usage with '
                          'Chewie-NS.\nContinue?\n')
                ptf_answer = fo.input_timeout(prompt, ct.prompt_timeout)
            if ptf_num > 1:
                logger.warning('Prodigal training file is not any of the '
                               f'{ptf_num} used in previous runs.')
                prompt = ('Continue?\n')
                ptf_answer = fo.input_timeout(prompt, ct.prompt_timeout)
        else:
            ptf_answer = 'yes'

        if ptf_answer.lower() not in ['y', 'yes']:
            sys.exit(1)
        else:
            unmatch = True

    return unmatch


def validate_ptf(ptf_path, schema_directory, schema_ptfs, force_continue):
    """Validate the path to the Prodigal training file an its hash value.

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

    # determine PTF checksum
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

    # run parameters values
    run_params = {'bsr': blast_score_ratio,
                  'translation_table': translation_table,
                  'minimum_locus_length': minimum_length,
                  'size_threshold': size_threshold}

    # determine user provided arguments values that differ from default
    unmatch_params = {k: v
                      for k, v in run_params.items()
                      if v not in schema_params[k] and v is not None}
    # determine arguments values not provided by user
    default_params = {k: schema_params[k][0]
                      for k, v in run_params.items()
                      if v is None}

    # update arguments for current run
    for k in run_params:
        if k in default_params:
            run_params[k] = default_params[k]

    if len(unmatch_params) > 0:
        params_diffs = [[p, ':'.join(map(str, schema_params[p])),
                         str(unmatch_params[p])]
                        for p in unmatch_params]
        params_diffs_text = ['{:^20} {:^20} {:^10}'.format('Argument', 'Schema', 'Provided')]
        params_diffs_text += ['{:^20} {:^20} {:^10}'.format(p[0], p[1], p[2]) for p in params_diffs]
        logger.warning('Provided parameter values differ from values used for '
                       'schema creation:\n{"\n".join(params_diffs_text)}')
        if force_continue is False:
            prompt = ('\nContinuing might lead to results not '
                      'consistent with previous runs.\nProviding '
                      'parameter values that differ from the ones '
                      'used for schema creation will also invalidate '
                      'the schema for uploading and synchronization '
                      'with Chewie-NS.\nContinue? (yes/no)\n')
            params_answer = fo.input_timeout(prompt, ct.prompt_timeout)
        else:
            params_answer = 'yes'

        if params_answer.lower() not in ['y', 'yes']:
            sys.exit(1)
        else:
            # append new arguments values to configs values
            for p in unmatch_params:
                schema_params[p].append(unmatch_params[p])

    # default is to get the training file in schema directory
    schema_ptfs = schema_params['prodigal_training_file']
    ptf_path, ptf_hash, unmatch = validate_ptf(ptf_path, schema_directory,
                                               schema_ptfs, force_continue)

    run_params['ptf_path'] = ptf_path
    if unmatch is True:
        schema_params['prodigal_training_file'].append(ptf_hash)
        unmatch_params['prodigal_training_file'] = ptf_hash

    # save updated schema config file
    if len(unmatch_params) > 0:
        fo.pickle_dumper(schema_params, config_file)

    return run_params


def write_gene_list(schema_dir):
    """Create file with list of loci in a schema.

    Parameters
    ----------
    schema_dir : str
        Path to the directory with schema files.

    Returns
    -------
    A list with two elements. A boolean value that
    is True if the file with the list of genes was
    created and False otherwise. The second element
    is the path to the created file.
    """
    loci_files = [file
                  for file in os.listdir(schema_dir)
                  if '.fasta' in file]
    loci_list_file = fo.join_paths(schema_dir, ['.genes_list'])
    fo.pickle_dumper(loci_files, loci_list_file)
    logger.debug(f'Wrote list of loci to {loci_list_file}')

    return loci_list_file


def write_schema_config(blast_score_ratio, ptf_hash, translation_table,
                        minimum_sequence_length, chewie_version, size_threshold,
                        word_size, window_size, clustering_sim, representative_filter,
                        intra_filter, output_directory):
    """Write argument values used to create a schema to a file.

    Parameters
    ----------
    blast_score_ratio : float
        BLAST Score Ratio value used to create the
        schema.
    ptf_hash : str
        BLAKE2 hash of the Prodigal training file
        content.
    translation_table : int
        Genetic code used to predict and translate
        coding sequences.
    minimum_sequence_length : int
        Minimum sequence length, sequences with a
        length value lower than this value are not
        included in the schema.
    chewie_version : str
        Version of the chewBBACA suite used to create
        the schema.
    size_threshold : float
        Sequence size variation percentage threshold,
        new alleles cannot have a length value that
        deviates +/- than this value in relation to the
        locus's representative sequence.
    word_size : int
        Word/k value used to cluster protein sequences
        during schema creation and allele calling.
    clustering_sim : float
        Proportion of k-mers/minimizers that two proteins
        need to have in common to be clustered together.
    representative_filter : float
        Proportion of k-mers/minimizers that a clustered
        protein has to have in common with the representative
        protein of the cluster to be considered the same gene.
    intra_filter : float
        Proportion of k-mers/minimizers that clustered
        proteins have to have in common to be considered
        of the same gene.
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
    size_threshold = None if size_threshold in [None, 'None'] else float(size_threshold)

    params = {}
    params['bsr'] = [float(blast_score_ratio)]
    params['prodigal_training_file'] = [ptf_hash]
    params['translation_table'] = [int(translation_table)]
    params['minimum_locus_length'] = [int(minimum_sequence_length)]
    params['chewBBACA_version'] = [chewie_version]
    params['size_threshold'] = [size_threshold]
    params['word_size'] = [word_size]
    params['window_size'] = [window_size]
    params['cluster_sim'] = [clustering_sim]
    params['representative_filter'] = [representative_filter]
    params['intraCluster_filter'] = [intra_filter]

    config_file = os.path.join(output_directory, '.schema_config')
    fo.pickle_dumper(params, config_file)
    logger.debug(f'Wrote schema config to {config_file}')

    return config_file


def read_configs(schema_path, filename):
    """Read file with schema config values.

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
        sys.exit('Could not find a valid config file.')

    return configs


def check_input_type(input_path, output_file, parent_dir=None):
    """ Checks if the input path is for a file or for a
        directory. If the path is for a directory, the
        function creates a file with the list of paths
        to FASTA files in the directory.

    Parameters
    ----------
    input_path : str
        Path to file or directory.
    output_file : str
        Path to the output file with the list of FASTA
        files.
    parent_dir : str
        Parent directory to add to construct paths
        to input files when users provide a file
        with file names.

    Returns
    -------
    output_file : str
        Path to a file with a list of paths for FASTA
        files.

    Raises
    ------
    SystemExit
        - If there were no FASTA files in the directory.
        - If input path is not a valid path for a file or
          for a directory.
    """

    # check if input argument is a file or a directory
    if os.path.isfile(input_path):
        # check if it's not a FASTA file
        if fao.validate_fasta(input_path) is True:
            logger.error('Input file is a FASTA file. Please provide '
                         'the path to the parent directory that contains '
                         'the FASTA files or a file with the list of full '
                         'paths to the FASTA files (one per line).')
            sys.exit(1)

        # read list of input files
        with open(input_path, 'r') as infile:
            lines = list(csv.reader(infile))
            lines = [f[0] for f in lines]

        # list of input genes must have full paths
        if parent_dir is not None:
            # add parent directory path if necessary
            lines = [os.path.join(parent_dir, f)
                     if parent_dir not in f
                     else f
                     for f in lines]
            # add FASTA extension if it's missing
            lines = [file+'.fasta'
                     if file.endswith('.fasta') is False
                     else file
                     for file in lines]

        # check that all files exist
        missing = [file for file in lines
                   if os.path.exists(file) is False]
        if len(missing) > 0:
            missing_str = '\n'.join(missing)
            logger.error('Could not find some of the files provided in '
                         'the input list. Please verify that you\'ve '
                         'provided the full paths to valid input '
                         f'files. Files missing: \n{missing_str}')
            sys.exit(1)
        # save file paths to output file
        else:
            with open(output_file, 'w') as outfile:
                outfile.write('\n'.join(lines))

    elif os.path.isdir(input_path):
        # we need to get only files with FASTA extension
        files = os.listdir(input_path)
        files = fo.filter_files(files, ct.FASTA_SUFFIXES)
        # get absolute paths
        files = [os.path.join(input_path, file) for file in files]
        # filter any directories that might end with FASTA extension
        files = [file for file in files if os.path.isdir(file) is False]

        # only keep files whose content is typical of a FASTA file
        fasta_files = fao.filter_non_fasta(files)

        # if there are FASTA files
        if len(fasta_files) > 0:
            # store full paths to FASTA files
            with open(output_file, 'w') as f:
                for file in fasta_files:
                    f.write(file + '\n')
        else:
            logger.error('Could not get input files. Please '
                         'provide a directory with FASTA files '
                         'or a file with the list of full paths '
                         'to the FASTA files and ensure that '
                         'filenames end with one of the '
                         f'following suffixes: {ct.FASTA_SUFFIXES}.')
            sys.exit(1)
    else:
        logger.error('Input argument is not a valid directory or '
                     'file with a list of paths. Please provide a '
                     'valid input, either a folder with FASTA files '
                     'or a file with the list of full paths to FASTA '
                     'files (one per line).')
        sys.exit(1)

    return output_file


def check_hash_type(hash_function):
    """Check if hash type is available in hashlib or zlib libraries.

    Parameters
    ----------
    hash_function : str
        String used to identify the hash function in hashlib or zlib.

    Returns
    -------
    hash_function : str
        String used to identify the hash function in hashlib or zlib.
    """
    # get hash function
    in_hashlib = getattr(hashlib, hash_function, None)
    if in_hashlib is None:
        in_zlib = getattr(zlib, hash_function, None)
        if in_zlib is None:
            logger.error(f'{hash_function} hash function is not available in '
                         'hashlib or zlib modules.')
            sys.exit(1)

    return hash_function


def check_schema_directory(schema_path):
    """Check if path to schema directory includes the necessary files.

    Parameters
    ----------
    schema_path : str
        Path to the schema directory.

    Returns
    -------
    schema_path : str
        Path to the schema directory.
    """
    try:
        schema_files = os.listdir(schema_path)
        # exit if there is no 'short' directory or if there are no FASTA files
        if 'short' not in schema_files or len(fo.filter_files(schema_files, ['.fasta'])) == 0:
            logger.error('Path provided as schema directory does not include '
                         'all the necessary files. Please verify that you '
                         'have passed the correct path to the schema.')
            sys.exit(1)
    except Exception as e:
        logger.exception(e)
        sys.exit(1)

    return schema_path
