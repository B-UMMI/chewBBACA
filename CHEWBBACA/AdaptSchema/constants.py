#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the AdaptSchema module.
"""


# Parameter names
ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
ADAPTSCHEMA_LOCI_LIST_ARGNAME = "loci_list"
ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME = "gene_predictor"
ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME = "gene_prediction_options"
ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME = "clustering_options"
ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME = "blast_score_ratio"
ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME = "minimum_length"
ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME = "translation_table"
ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME = "size_threshold"
ADAPTSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"
ADAPTSCHEMA_BLAST_PATH_ARGNAME = "blast_path"
ADAPTSCHEMA_SIZE_FILTER_ARGNAME = "size_filter"

# Argument types
ADAPTSCHEMA_ARGUMENT_TYPES = {
	"schema_directory": str,
	"output_directory": str,
	"loci_list": str,
	"gene_predictor": str,
	"gene_prediction_options": str,
	"clustering_options": str,
	"blast_score_ratio": float,
	"minimum_length": int,
	"translation_table": int,
	"size_threshold": float,
	"cpu_cores": int,
	"blast_path": str,
	"size_filter": bool
}

# Parameter descriptions
ADAPTSCHEMA_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the directory of the schema to adapt. The schema must contain one "
							 "FASTA file per gene/locus.")
ADAPTSCHEMA_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory where the adapted schema will be created.")
ADAPTSCHEMA_LOCI_LIST_DESCRIPTION = ("Path to a file with the list of loci in the schema that the process should "
							 "adapt (one per line, full paths or loci IDs).")
ADAPTSCHEMA_GENE_PREDICTOR_DESCRIPTION = ("Specify which gene prediction software to use. Default is Pyrodigal "
								 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
								 "for prokaryotic and eukaryotic genomes.")
ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_DESCRIPTION = ("List of arguments passed to configure the gene prediction. When providing "
							 "genome assemblies in FASTA format, the list of arguments for each parameter "
							 "used to configure the gene prediction can be passed as the long format of "
							 "the parameter name followed by the argument value (e.g., pyrodigal-training"
							 "-file=/path/to/file).")
ADAPTSCHEMA_CLUSTERING_OPTIONS_DESCRIPTION = ("List of arguments passed to configure the clustering of predicted "
								"coding sequences (CDSs). The list of arguments for each parameter "
								"used to configure the clustering can be passed as the long format "
								"of the parameter name followed by the argument value (e.g., word-size=).")
ADAPTSCHEMA_BLAST_SCORE_RATIO_DESCRIPTION = ("BLAST Score Ratio (BSR) value. The process selects representative alleles "
							 "for each locus based on this value. Representative alleles are selected "
							 "until all alleles in a locus align against one of the representatives "
							 "with a BSR >= than the specified value.")
ADAPTSCHEMA_MINIMUM_LENGTH_DESCRIPTION = ("Minimum sequence length value stored in the schema config file. The "
							 "schema adaptation process will only discard sequences smaller than "
							 "this value if the --size-filter parameter is provided.")
ADAPTSCHEMA_TRANSLATION_TABLE_DESCRIPTION = ("Genetic code used for allele translation. This value is ignored if "
							 "a valid training file is passed to `--ptf`, `--training-file`.")
ADAPTSCHEMA_SIZE_THRESHOLD_DESCRIPTION = ("Allele size variation threshold value stored in the schema config "
							 "file. The schema adaptation process will only discard alleles with "
							 "a size that deviates from the locus length mode +- the size theshold "
							 "value if the --size-filter parameter is provided.")
ADAPTSCHEMA_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")
ADAPTSCHEMA_SIZE_FILTER_DESCRIPTION = ("Apply the minimum length and size threshold values to filter out "
							 "alleles during schema adaptation.")

# File header for file with summary statistics created by the AdaptSchema module
ADAPTSCHEMA_STATS_HEADER = ["Locus", "Total_alleles", "Valid_alleles", "No. Representatives"]
