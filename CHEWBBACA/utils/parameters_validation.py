#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import sys
import argparse

try:
    from utils import constants as cnts
    from utils import auxiliary_functions as aux
except:
    from CHEWBBACA.utils import constants as cnts
    from CHEWBBACA.utils import auxiliary_functions as aux


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

    # prog is the name of the program 'ex: chewBBACA.py'
    def __init__(self, prog, indent_increment=2, max_help_position=56, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    # override split lines method
    def _split_lines(self, text, width):
        lines = super()._split_lines(text, width) + ['']
        return lines


def arg_list(arg, arg_name):
    """
    """

    if isinstance(arg, list) is True:
        if len(arg) > 1:
            sys.exit('\nMultiple {0} values.'.format(arg_name))
        else:
            arg = arg[0]

    return arg


# custom functions to validate arguments type and value
def bsr_type(arg, min_value=cnts.BSR_MIN, max_value=cnts.BSR_MAX):
    """
    """

    arg = arg_list(arg, 'BLAST Score Ratio')

    try:
        schema_bsr = float(arg)
        if schema_bsr >= min_value and schema_bsr <= max_value:
            valid = schema_bsr
        elif schema_bsr < min_value or schema_bsr > max_value:
            sys.exit('\nBSR value is not contained in the '
                     '[0.0, 1.0] interval.')
    except Exception:
        sys.exit('\nInvalid BSR value of {0}. BSR value must be contained'
                 ' in the [0.0, 1.0] interval.'.format(arg[0]))

    return valid


def minimum_sequence_length_type(arg, min_value=cnts.MSL_MIN, max_value=cnts.MSL_MAX):
    """
    """

    arg = arg_list(arg, 'minimum sequence length')

    try:
        schema_ml = int(arg)
        if schema_ml >= min_value and schema_ml <= max_value:
            valid = schema_ml
        elif schema_ml < min_value or schema_ml > max_value:
            sys.exit('\nInvalid minimum sequence length value. '
                     'Must be equal or greater than 0.')
    except Exception:
        sys.exit('\nInvalid minimum sequence length value used to '
                 'create schema. Value must be a positive integer.')

    return valid


def size_threshold_type(arg, min_value=cnts.ST_MIN, max_value=cnts.ST_MAX):
    """
    """

    arg = arg_list(arg, 'size threshold')

    try:
        schema_st = float(arg)
        if schema_st >= min_value and schema_st <= max_value:
            valid = schema_st
        elif schema_st < min_value or schema_st > max_value:
            sys.exit('\nInvalid size threshold value. '
                     'Must be contained in the [0.0, 1.0] interval.')
    except Exception:
        if arg in [None, 'None']:
            valid = None
        else:
            sys.exit('\nInvalid size threshold value used to '
                     'create schema. Value must be None or a '
                     'positive float in the [0.0, 1.0] interval.')

    return valid


def translation_table_type(arg, genetic_codes=cnts.GENETIC_CODES):
    """
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
        lines = ['\t{0}: {1}'.format(k, v) for k, v in genetic_codes.items()]
        gc_table = '\n{0}\n'.format('\n'.join(lines))

        sys.exit('\nInvalid genetic code value.\nValue must correspond to '
                 'one of the accepted genetic codes\n\nAccepted genetic '
                 'codes:\n{0}'.format(gc_table))

    return valid


def validate_cv(arg, chewie_versions=cnts.CHEWIE_VERSIONS):
    """
    """

    arg = arg_list(arg, 'chewBBACA version')

    if arg in chewie_versions:
        valid = arg
    else:
        sys.exit('\nSchema created with chewBBACA version that '
                 'is not suitable to work with the NS.')

    return valid


def validate_ws(arg, min_value=cnts.WS_MIN, max_value=cnts.WS_MAX):
    """
    """

    arg = arg_list(arg, 'word size')

    try:
        if arg == None:
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


def validate_cs(arg, min_value=cnts.CS_MIN, max_value=cnts.CS_MAX):
    """
    """

    arg = arg_list(arg, 'clustering similarity')

    try:
        if arg == None:
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


def validate_rf(arg, min_value=cnts.RF_MIN, max_value=cnts.RF_MAX):
    """
    """

    arg = arg_list(arg, 'representative filter')

    try:
        if arg == None:
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


def validate_if(arg, min_value=cnts.IF_MIN, max_value=cnts.IF_MAX):
    """
    """

    arg = arg_list(arg, 'intra-cluster filter')

    try:
        if arg == None:
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


def validate_ptf(arg, input_path):
    """
    """

    arg = arg_list(arg, 'Prodigal training file')

    if arg == '':
        sys.exit('Cannot upload a schema that was created '
                 'without a Prodigal training file.')

    schema_ptfs = [os.path.join(input_path, file)
                   for file in os.listdir(input_path) if '.trn' in file]
    if len(schema_ptfs) == 1:
        schema_ptf = schema_ptfs[0]
        ptf_hash = aux.hash_file(schema_ptf, 'rb')
        if ptf_hash == arg:
            valid = [schema_ptf, arg]
        else:
            sys.exit('Training file in schema directory is not the original.')
    elif len(schema_ptfs) > 1:
        sys.exit('More than one training file in schema directory.')
    elif len(schema_ptfs) == 0:
        sys.exit('Could not find a valid training file in schema directory.')

    return valid


def validate_ns_url(arg):
    """
    """

    if arg in cnts.HOST_NS:
        ns_url = cnts.HOST_NS[arg]
    else:
        ns_url = arg

    # sync schema has None by default to get ns_url in schema URI
    if ns_url is not None:
        # check if server is up
        conn = aux.check_connection(ns_url)
        if conn is False:
            sys.exit('Failed to establish a connection to the Chewie-NS '
                     'at {0}.'.format(ns_url))

    return ns_url
