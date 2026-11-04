#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


COMPUTEDISTANCES_MODULE_NAME = "ComputeDistances"

COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME = "allelic_profiles"
COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME = "distance_computation_method"
COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME = "output_format"
COMPUTEDISTANCES_COMPUTE_SIMILARITY_ARGNAME = "compute_similarity"
COMPUTEDISTANCES_NO_MASK_ARGNAME = "no_mask"
COMPUTEDISTANCES_CPU_CORES_ARGNAME = "cpu_cores"

COMPUTEDISTANCES_ARGUMENT_TYPES = {
	"allelic_profiles": str,
	"output_directory": str,
	"distance_computation_method": str,
	"output_format": str,
	"compute_similarity": bool,
	"no_mask": bool,
	"cpu_cores": int
}

DEFAULT_DISTANCE_METHOD = "hamming"

DISTANCE_METHODS = ["hamming", "jaccard", "loci", "core"]

DEFAULT_OUTPUT_FORMAT = "lower_triangular"

OUTPUT_FORMATS = ["upper_triangular", "lower_triangular", "symmetric", "table"]

# Exceptions

INVALID_DISTANCE_METHOD = ("Invalid distance method specified.")

INVALID_OUTPUT_FORMAT = ("Invalid output format specified.")
