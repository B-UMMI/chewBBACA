#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


COMPUTECORE_MODULE_NAME = "ComputeCore"

COMPUTECORE_PROFILES_ARGNAME = "profiles"
COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME = "loci_presence_threshold"
COMPUTECORE_SAMPLE_STEP_ARGNAME = "sample_step"
COMPUTECORE_COMPUTE_ACCESSORY_ARGNAME = "compute_accessory"
COMPUTECORE_RAREFACTION_ANALYSIS_ARGNAME = "rarefaction_analysis"
COMPUTECORE_PERMUTATION_NUMBER_ARGNAME = "permutation_number"
COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME = "permutation_samples"
COMPUTECORE_CPU_CORES_ARGNAME = "cpu_cores"
COMPUTECORE_EXCLUDE_LOCI_ARGNAME = "exclude_loci"
COMPUTECORE_EXCLUDE_GENOMES_ARGNAME = "exclude_genomes"

COMPUTECORE_ARGUMENT_TYPES = {
	"profiles": str,
	"output_directory": str,
	"loci_presence_threshold": float,
	"sample_step": int,
	"compute_accessory": bool,
	"rarefaction_analysis": bool,
	"permutation_number": int,
	"permutation_samples": int,
	"cpu_cores": int,
	"exclude_loci": str,
	"exclude_genomes": str
}

# Default loci presence thresholds used to compute the cgMLST
CGMLST_THRESHOLDS = [0.95, 0.99, 1]
CGMLST_THRESHOLD_MIN = 0
CGMLST_THRESHOLD_MAX = 1

PRESENCE_ABSENCE_PROFILES_BASENAME = "presence_absence_profiles.tsv"

# Header for TSV file with sample loci presence stats
SAMPLE_LOCUS_PRESENCE_HEADER = ["Sample", "Locus presence count", "Locus presence proportion"]
SAMPLE_LOCUS_PRESENCE_BASENAME = "sample_presence_stats.tsv"
# Header for TSV file with locus sample presence stats
LOCUS_SAMPLE_PRESENCE_HEADER = ["Locus", "Sample presence count", "Sample presence proportion"]
LOCUS_SAMPLE_PRESENCE_BASENAME = "locus_presence_stats.tsv"

PERMUTATION_NUMBER_DEFAULT = 100
PERMUTATION_SAMPLES_DEFAULT = None

# Exceptions

INVALID_EXCLUDE_LOCI = ("The following loci IDs in the list provided to `--exclude-loci` do not match "
                        "any of the loci IDs in the input results:\n{0}\nPlease ensure that all loci "
						"IDs in the input list match loci IDs in the profiles.")

INVALID_EXCLUDE_SAMPLES = ("The following sample IDs in the list provided to `--exclude-samples` do not match "
                           "any of the sample IDs in the input results:\n{0}\nPlease ensure that all loci "
						   "IDs in the input list match loci IDs in the profiles.")

MISSING_EXCLUDE_LOCI = ("Path to file with list of loci is not valid.")
