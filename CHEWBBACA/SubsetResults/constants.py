#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the SubsetResults module.
"""


# Parameter names
SUBSETRESULTS_RESULTS_FILES_ARGNAME = "results_files"
SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
SUBSETRESULTS_LOCI_LIST_ARGNAME = "loci_list"
SUBSETRESULTS_SAMPLE_LIST_ARGNAME = "sample_list"
SUBSETRESULTS_INVERT_LOCI_ARGNAME = "invert_loci"
SUBSETRESULTS_INVERT_SAMPLES_ARGNAME = "invert_samples"

# Argument types
SUBSETRESULTS_ARGUMENT_TYPES = {
	"results_files": str,
	"output_directory": str,
	"loci_list": str,
	"sample_list": str,
	"invert_loci": bool,
	"invert_samples": bool
}

# Parameter descriptions
SUBSETRESULTS_RESULTS_FILES_DESCRIPTION = ("Path to the directory containing the files to be subsetted.")
SUBSETRESULTS_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory.")
SUBSETRESULTS_LOCI_LIST_DESCRIPTION = ("Path to a TXT/TSV file containing a list of loci to select, one locus identifier "
							 "per line. If the file contains multiple columns, the loci identifiers must be in "
							 "the first column.")
SUBSETRESULTS_SAMPLE_LIST_DESCRIPTION = ("Path to a TXT/TSV file containing a list of samples to select, one sample identifier "
							 "per line. If the file contains multiple columns, the sample identifiers must be in "
							 "the first column.")
SUBSETRESULTS_INVERT_LOCI_DESCRIPTION = ("If provided, the process will select the loci that are not in the input loci list.")
SUBSETRESULTS_INVERT_SAMPLES_DESCRIPTION = ("If provided, the process will select the samples that are not in the input samples list.")

# Exceptions
SUBSETRESULTS_MISSING_LISTS_EXCEPTION = ('Did not provide a list of loci or samples. Please '
										 'provide at least one of those lists to subset results.')

SUBSETRESULTS_MISSING_FILES_EXCEPTION = ('Did not find any of the valid file types to subset. Please make '
		   		 						 'sure that the input directory contains at least the file type '
				 						 'created by chewBBACA to store the allelic profiles (the basename '
				 						 'must be `results_alleles.tsv`).')

SUBSETRESULTS_MISSING_PROFILES_EXCEPTION = ('Did not find the file type used by chewBBACA to store the allelic '
		   		 						    'profiles (`results_alleles.tsv`). Please make sure that the input '
				 						    'directory contains at least that file type.')

SUBSETRESULTS_ABSENT_LOCI_EXCEPTION = ('Please make sure that all loci IDs in the list provided to '
		 		  	 				   '`--loci-list` match loci IDs in the input results.')

SUBSETRESULTS_ABSENT_SAMPLES_EXCEPTION = ('Please make sure that all sample IDs in the list provided to '
				  	 					  '`--sample-list` match sample IDs in the input results.')

MISSING_EXCLUDE_LOCI = ("Path to file with list of loci is not valid.")

MISSING_EXCLUDE_GENOMES = ("Path to file with list of genomes is not valid.")
