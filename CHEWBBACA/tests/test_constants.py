#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


Code documentation
------------------
"""


# Valid FASTA file extensions
FASTA_EXTENSIONS = [".fasta", ".fna", ".ffn", ".fa", ".fas"]

INVALID_BSR = ('\nBSR value is not contained in the [0.0, 1.0] interval.')
INVALID_BSR_TYPE = ('\nInvalid BSR value of {0}. BSR value must be contained in the [0.0, 1.0] interval.')

INVALID_MINLEN = ('\nInvalid minimum sequence length value. Must be equal or greater than 0.')
INVALID_MINLEN_TYPE = ('\nInvalid minimum sequence length value. Value must be a positive integer.')

INVALID_ST = ('\nInvalid size threshold value. Must be contained in the [0.0, 1.0] interval.')
INVALID_ST_TYPE = ('\nInvalid size threshold value used to create schema. Value must be None or a positive float in the [0.0, 1.0] interval.')

# Input path is neither a file nor a directory
INVALID_INPUT_PATH = ('Input argument is not a valid directory or '
					  'file with a list of paths to FASTA files. Please '
					  'provide a valid input, either a folder with FASTA '
					  'files or a file with the list of full paths to FASTA '
					  'files (one per line and ending with one of the '
					  f'following file extensions: {FASTA_EXTENSIONS}).')

# Input directory does not contain FASTA files
MISSING_FASTAS_EXCEPTION = ('Could not get input files. Please provide '
							'a directory with FASTA files or a file with '
							'the list of full paths to the FASTA files '
							'and ensure that filenames end with one of '
							f'the following extensions: {FASTA_EXTENSIONS}.')

INVALID_PTF_PATH = ("Cannot find the specified Prodigal training file in {0}."
					"\nPlease provide a path to a valid training file.")
