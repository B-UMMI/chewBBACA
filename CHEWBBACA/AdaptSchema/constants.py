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

# File header for file with summary statistics created by the AdaptSchema module
ADAPTSCHEMA_STATS_HEADER = ["Locus", "Total_alleles", "Valid_alleles", "No. Representatives"]
