#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


MERGERESULTS_MODULE_NAME = "MergeResults"

MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME = "results_directories"
MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
MERGERESULTS_COMMON_ARGNAME = "common"

MERGERESULTS_ARGUMENT_TYPES = {
	"results_directories": str,
	"output_directory": str,
	"common": bool
}

# Exceptions

MERGERESULTS_INPUTFILE_EXCEPTION = ("Input paths should all be folders. The following are not folders:\n{0}")

MERGERESULTS_NOCOMMONLOCI_EXCEPTION = ('Results do not have loci in common.')

MERGERESULTS_LOCIDIFFER_EXCEPTION = ('Files have different sets of loci. Please provide '
				 					 'files with the results for the same set of loci or '
				 					 'provide the "--common" parameter to create merged files '
				 					 'for the set of loci shared by all results folders.')
