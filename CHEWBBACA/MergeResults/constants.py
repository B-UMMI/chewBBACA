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

# List of valid files for the SubsetResutls and MergeResults modules
# Includes all standard files created by the AlleleCall module plus the ones created with specific parameters
SUPPORTED_RESULTS_FILES = [ALLELIC_PROFILES_BASENAME, COORDINATES_PROFILES_BASENAME, SAMPLE_STATISTICS_BASENAME,
                           LOCI_STATS_BASENAME, INVALID_CDS_BASENAME, GENE_COORDINATES_BASENAME, MISSING_FASTA_BASENAME,
						   MISSING_TSV_BASENAME, UNCLASSIFIED_BASENAME, PARALOGOUS_COUNTS_BASENAME, PARALOGOUS_TABLE_BASENAME,
						   INFERRED_ALLELES_BASENAME, PRESENCE_ABSENCE_PROFILES_BASENAME]

# Exceptions

MERGERESULTS_INPUTFILE_EXCEPTION = ("Input paths should all be folders. The following are not folders:\n{0}")

MERGERESULTS_NOCOMMONLOCI_EXCEPTION = ('Results do not have loci in common.')

MERGERESULTS_LOCIDIFFER_EXCEPTION = ('Files have different sets of loci. Please provide '
				 					 'files with the results for the same set of loci or '
				 					 'provide the "--common" parameter to create merged files '
				 					 'for the set of loci shared by all results folders.')
