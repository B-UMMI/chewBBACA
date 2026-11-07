#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the MergeResults module.
"""


# Parameter names
MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME = "results_directories"
MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
MERGERESULTS_COMMON_ARGNAME = "common"

# Argument types
MERGERESULTS_ARGUMENT_TYPES = {
	"results_directories": str,
	"output_directory": str,
	"common": bool
}

# Parameter descriptions
MERGERESULTS_RESULTS_DIRECTORIES_DESCRIPTION = ("Paths to the directories containing the results files created by chewBBACA. The "
							 "results must have been determined with the same schema and share all the loci "
							 "or a subset of the loci if using the --common parameter.")
MERGERESULTS_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory.")
MERGERESULTS_COMMON_DESCRIPTION = ("Merge the results based on the subset of loci shared between all inputs.")

# Exceptions
MERGERESULTS_INPUTFILE_EXCEPTION = ("Input paths should all be folders. The following are not folders:\n{0}")

MERGERESULTS_NOCOMMONLOCI_EXCEPTION = ('Results do not have loci in common.')

MERGERESULTS_LOCIDIFFER_EXCEPTION = ('Files have different sets of loci. Please provide '
				 					 'files with the results for the same set of loci or '
				 					 'provide the "--common" parameter to create merged files '
				 					 'for the set of loci shared by all results folders.')
