#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the CreateSchema module.
"""


# Parameter names
CREATESCHEMA_INPUT_PATH_ARGNAME = "input_path"
CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
CREATESCHEMA_SCHEMA_NAME_ARGNAME = "schema_name"
CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME = "blast_score_ratio"
CREATESCHEMA_MINIMUM_LENGTH_ARGNAME = "minimum_length"
CREATESCHEMA_TRANSLATION_TABLE_ARGNAME = "translation_table"
CREATESCHEMA_SIZE_THRESHOLD_ARGNAME = "size_threshold"
CREATESCHEMA_GENE_PREDICTOR_ARGNAME = "gene_predictor"
CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME = "gene_prediction_options"
CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME = "clustering_options"
CREATESCHEMA_BLAST_PATH_ARGNAME = "blast_path"
CREATESCHEMA_CDS_INPUT_ARGNAME = "cds_input"
CREATESCHEMA_CPU_CORES_ARGNAME = "cpu_cores"
CREATESCHEMA_NO_CDS_RENAMING_ARGNAME = "no_cds_renaming"
CREATESCHEMA_NO_CLEANUP_ARGNAME = "no_cleanup"

# Argument types
CREATESCHEMA_ARGUMENT_TYPES = {
	"input_path": str,
	"output_directory": str,
	"schema_name": str,
	"blast_score_ratio": float,
	"minimum_length": int,
	"translation_table": int,
	"size_threshold": float,
	"gene_predictor": str,
	"gene_prediction_options": str,
	"clustering_options": str,
	"blast_path": str,
	"cds_input": bool,
	"cpu_cores": int,
	"no_cds_renaming": bool,
	"no_cleanup": bool,
}

# Parameter descriptions
CREATESCHEMA_INPUT_PATH_DESCRIPTION = ("Path to the directory that contains the input FASTA files or to a file"
							 " with a list of full paths to FASTA files, one per line.")
CREATESCHEMA_OUTPUT_DIRECTORY_DESCRIPTION = ("Output directory where the process will store intermediate files and "
							 "create the schema's directory.")
CREATESCHEMA_SCHEMA_NAME_DESCRIPTION = ("Name given to the schema folder.")
CREATESCHEMA_BLAST_SCORE_RATIO_DESCRIPTION = ("BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp "
							 "alignment and aligned sequences with a BSR >= than the defined value "
							 "are considered to be alleles of the same gene.")
CREATESCHEMA_MINIMUM_LENGTH_DESCRIPTION = ("Minimum sequence length value. Predicted coding sequences (CDSs) "
							 "shorter than this value are excluded.")
CREATESCHEMA_TRANSLATION_TABLE_DESCRIPTION = ("The translation table, matching a specific genetic code, used for gene "
							 "prediction and to translate CDSs. This value is ignored if a valid Pyrodigal "
							 "training file is passed to `-gpo`, `--gene-prediction-options`. In that case, "
							 "the translation table specified in the Pyrodigal training file is used.")
CREATESCHEMA_SIZE_THRESHOLD_DESCRIPTION = ("CDS size variation threshold. This value is used to define how much allele "
							 "size can vary for each locus (e.g. setting this value to 0.2 will define a "
							 "threshold that limits allele size variation to +/-20%% the allele size mode "
							 "for each locus). This value is not used in the CreateSchema process, but is "
							 "added to the schema's config file to be used during allele calling to check "
							 "if the sequence size of the potential new alleles is within the defined size "
							 "threshold.")
CREATESCHEMA_CPU_CORES_DESCRIPTION = ("Number of CPU cores that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores).")


# Proteins to cluster are divided into a maximum of 40 smaller groups in CreateSchema
# Since we do not have pre-defined clusters during schema creation, dividing based on
# the number of CPU cores would lead to variable results. Therefore, we use a fixed
# number of groups to divide the proteins to ensure results are consistent across
# different runs
CREATESCHEMA_CLUSTERING_NGROUPS = 40

# Default name for the folder containing the schema seeds created by the CreateSchema module
SCHEMA_NAME_DEFAULT = "schema_seed"
