#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the ComputeCore module.
"""


# Parameter names
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

# Argument types
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

# Parameter descriptions
COMPUTECORE_PROFILES_DESCRIPTION = ("Path to the TSV file that contains the allelic profiles determined by "
							 "the AlleleCall module.")
COMPUTECORE_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the directory where the process will store the output files.")
COMPUTECORE_LOCI_PRESENCE_THRESHOLD_DESCRIPTION = ("Loci/genes that constitute the core genome must be in a proportion of "
							 "genomes that is at least equal to this value. Provide multiple values "
							 "to compute the core genome for multiple threshold values.")
COMPUTECORE_SAMPLE_STEP_DESCRIPTION = ("The allele calling results are processed iteratively to evaluate the "
							 "impact of adding subsets of the results in computing the core genome. "
							 "The step value controls the number of allelic profiles added in each "
							 "iteration until all profiles are included.")
COMPUTECORE_COMPUTE_ACCESSORY_DESCRIPTION = ("Determine the set of accessory loci. The accessory genome corresponds "
							 "to all the loci not included in the core genome. The accessory genome "
							 "is determined for each core genome threshold.")
COMPUTECORE_RAREFACTION_ANALYSIS_DESCRIPTION = ("Perform rarefaction analysis to evaluate the stability of the core genome "
							 "are randomly selected and the number of core loci (for the core genome) or "
							 "total loci (for the pangenome) is computed for each subset from 1 to the "
							 "total number of samples. The average number of core loci (for the core genome) "
							 "or total loci (for the pangenome) is computed for each sample size subset "
							 "using the values of all permutations to plot the rarefaction curve. The "
							 "rarefaction analysis is performed for each threshold. A power law model "
							 "is fitted to the rarefaction curve to estimate the stability of the core genome "
							 "and the openness of the pangenome.")
COMPUTECORE_PERMUTATION_NUMBER_DESCRIPTION = ("Number of permutations for the rarefaction analysis. The rarefaction analysis "
							 "is repeated a number of times equal to the value provided to this parameter.")
COMPUTECORE_PERMUTATION_SAMPLES_DESCRIPTION = ("Number of samples randomly selected for each permutation. All samples will be "
							 "used if this value is not provided.")
COMPUTECORE_CPU_CORES_DESCRIPTION = ("Maximum number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads).")
COMPUTECORE_EXCLUDE_LOCI_DESCRIPTION = ("Path to a file with a list of loci identifiers to exclude from the analysis "
							 "(one locus identifier per line).")
COMPUTECORE_EXCLUDE_GENOMES_DESCRIPTION = ("Path to a file with a list of genome identifiers to exclude from the analysis "
							 "(one genome identifier per line).")

# Default loci presence thresholds used to compute the cgMLST
CGMLST_THRESHOLDS = [0.95, 0.99, 1]
CGMLST_THRESHOLD_MIN = 0
CGMLST_THRESHOLD_MAX = 1

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

MISSING_EXCLUDE_GENOMES = ("Path to file with list of genomes is not valid.")
