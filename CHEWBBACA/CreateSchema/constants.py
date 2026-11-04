#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


CREATESCHEMA_MODULE_NAME = "CreateSchema"

CREATESCHEMA_INPUT_FILES_ARGNAME = "input_files"
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

# Define expected types for the arguments
CREATESCHEMA_ARGUMENT_TYPES = {
	"input_files": str,
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


# Proteins to cluster are divided into a maximum of 40 smaller groups in CreateSchema
# Since we do not have pre-defined clusters during schema creation, dividing based on
# the number of CPU cores would lead to variable results. Therefore, we use a fixed
# number of groups to divide the proteins to ensure results are consistent across
# different runs
CREATESCHEMA_CLUSTERING_NGROUPS = 40

# CreateSchema exclusive section headers
FINAL_BLASTp = "Final BLASTp"

SCHEMA_NAME_DEFAULT = "schema_seed"
