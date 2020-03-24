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


import argparse

from utils import constants as cnts


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

    # prog is the name of the program 'ex: chewBBACA.py'
    def __init__(self, prog, indent_increment=2, max_help_position=56, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    # override split lines method
    def _split_lines(self, text, width):
        lines = super()._split_lines(text, width) + ['']
        return lines


# custom functions to validate arguments type and value
def bsr_type(arg, min_value=cnts.BSR_MIN, max_value=cnts.BSR_MAX):
    """
    """

    try:
        farg = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid float value: {0}'.format(arg))

    if farg < min_value or farg > max_value:
        raise argparse.ArgumentTypeError('value must be > {0} '
                                         'and <= {1}.'.format(min_value, max_value))

    return farg


def minimum_sequence_length_type(arg, min_value=cnts.MSL_MIN, max_value=cnts.MSL_MAX):
    """
    """

    try:
        iarg = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid int value: {0}'.format(arg))

    if iarg < min_value or iarg > max_value:
        raise argparse.ArgumentTypeError('value must be > {0} '
                                         'and <= {1}.'.format(min_value, max_value))

    return iarg


def size_threshold_type(arg, min_value=cnts.ST_MIN, max_value=cnts.ST_MAX):
    """
    """

    if arg is not None:
        try:
            farg = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError('invalid float value: {0}'.format(arg))

        if farg < min_value:
            raise argparse.ArgumentTypeError('value must be > {0} '
                                             'and <= {1}.'.format(min_value, max_value))
    else:
        farg = None

    return farg


def translation_table_type(arg, genetic_codes=cnts.GENETIC_CODES):
    """
    """

    try:
        iarg = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid int value: {0}'.format(arg))

    if iarg not in genetic_codes:
        # format available genetic codes into list
        lines = []
        for k, v in genetic_codes.items():
            new_line = '\t{0}: {1}'.format(k, v)
            lines.append(new_line)

        gc_table = '\n{0}\n'.format('\n'.join(lines))

        raise argparse.ArgumentTypeError('value must correspond to '
                                         'one of the accepted genetic '
                                         'codes\n\nAccepted genetic codes:\n{0}'.format(gc_table))

    return iarg
