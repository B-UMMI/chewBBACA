#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the PredictGenes module.
"""


# Parameter names
COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME = "allelic_profiles"
COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME = "distance_computation_method"
COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME = "output_format"
COMPUTEDISTANCES_COMPUTE_SIMILARITY_ARGNAME = "compute_similarity"
COMPUTEDISTANCES_NO_MASK_ARGNAME = "no_mask"
COMPUTEDISTANCES_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types
COMPUTEDISTANCES_ARGUMENT_TYPES = {
	"allelic_profiles": str,
	"output_directory": str,
	"distance_computation_method": str,
	"output_format": str,
	"compute_similarity": bool,
	"no_mask": bool,
	"cpu_cores": int
}

# Parameter descriptions
COMPUTEDISTANCES_ALLELIC_PROFILES_DESCRIPTION = ("Path to a TSV file containing allelic profiles determined by the AlleleCall "
							 "module.")
COMPUTEDISTANCES_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory where the process will store intermediate "
							 "and final results.")
COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_DESCRIPTION = ("Distance method used to compute the distance matrix. The module supports "
							 "the hamming, jaccard, loci (number of loci not shared), and core (number of "
							 "different alleles for core loci) methods.")
COMPUTEDISTANCES_OUTPUT_FORMAT_DESCRIPTION = ("Output format for the distance matrix (upper_triangular, lower_triangular, "
							 "symmetric, table).")
COMPUTEDISTANCES_COMPUTE_SIMILARITY_DESCRIPTION = ("Compute similarity values instead of distance values.")
COMPUTEDISTANCES_NO_MASK_DESCRIPTION = ("Do not mask missing data when computing the distance matrix. This option "
							 "is useful when the input profiles are already masked.")
COMPUTEDISTANCES_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be "
							 "used to run the process (chewie resets to a "
							 "lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads).")

# Default distance method used to compute pairwise distances
DEFAULT_DISTANCE_METHOD = "hamming"

# List of supported distance methods
DISTANCE_METHODS = ["hamming", "jaccard", "loci", "core"]

# Output format options
DEFAULT_OUTPUT_FORMAT = "lower_triangular"
OUTPUT_FORMATS = ["upper_triangular", "lower_triangular", "symmetric", "table"]

# Exceptions
INVALID_DISTANCE_METHOD = ("Invalid distance method specified.")

INVALID_OUTPUT_FORMAT = ("Invalid output format specified.")
