#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


SUBSETRESULTS_MODULE_NAME = "SubsetResults"

SUBSETRESULTS_RESULTS_FILES_ARGNAME = "results_files"
SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
SUBSETRESULTS_LOCI_LIST_ARGNAME = "loci_list"
SUBSETRESULTS_SAMPLE_LIST_ARGNAME = "sample_list"
SUBSETRESULTS_INVERT_LOCI_ARGNAME = "invert_loci"
SUBSETRESULTS_INVERT_SAMPLES_ARGNAME = "invert_samples"

SUBSETRESULTS_ARGUMENT_TYPES = {
	"results_files": str,
	"output_directory": str,
	"loci_list": str,
	"sample_list": str,
	"invert_loci": bool,
	"invert_samples": bool
}

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
