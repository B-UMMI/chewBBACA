#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""

PREDICTGENES_MODULE_NAME = "PredictGenes"

PREDICTGENES_INPUT_FILES_ARGNAME = "input_files"
PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
PREDICTGENES_GENE_PREDICTOR_ARGNAME = "gene_predictor"
PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME = "gene_prediction_options"
PREDICTGENES_TRANSLATION_TABLE_ARGNAME = "translation_table"
PREDICTGENES_CPU_CORES_ARGNAME = "cpu_cores"

# Define expected types for the arguments
PREDICTGENES_ARGUMENT_TYPES = {
	"input_files": str,
	"output_directory": str,
	"gene_predictor": str,
	"gene_prediction_options": str,
	"translation_table": int,
	"cpu_cores": int,
}

GENE_PREDICTORS = ["pyrodigal", "augustus"]

GENE_PREDICTOR_DEFAULT = "pyrodigal"

VALID_PARAMETERS = {"augustus": ["augustus_species", "augustus_output_formats", "augustus_path"],
					"pyrodigal": ["pyrodigal_training_file", "pyrodigal_mode",
							  	  "pyrodigal_output_formats", "pyrodigal_minimum_confidence",
							  	  "pyrodigal_training_reference", "pyrodigal_just_training"]}

PYRODIGAL_MODES = ["single", "meta"]

PYRODIGAL_OUTFMTS = ["genes", "translations", "gff", "genbank", "scores"]

PYRODIGAL_DEFAULT_OUTFMTS = ["genes"]

PYRODIGAL_MIN_CONFIDENCE = 0.0

PYRODIGAL_MAX_CONFIDENCE = 100.0

PYRODIGAL_DEFAULT_MODE = "single"

AUGUSTUS_ALIAS = "augustus"

AUGUSTUS_OUTFMTS = ["genes", "gff"]

AUGUSTUS_OUTFMT_DEFAULT = ["genes"]

AUGUSTUS_GFF_FILTERS = ["# start gene", "# coding sequence"]

GENE_COORDINATES_BASENAME = "gene_coordinates.tsv" # TSV file containing the coordinates of predicted CDSs
GENE_COORDINATES_EXCLUDED_BASENAME = "gene_coordinates_excluded.tsv" # TSV file containing the list of predicted CDSs that were excluded

# COLUMNS for the DataFrame with assembly statistics created by the PredictGenes module
ASSEMBLY_STATS_COLUMNS = ["FILE", "Number of contigs", "Average contig size",
						  "N50", "Assembly size", "GC content", "Missing values (Ns)",
						  "Total genes", "Excluded genes"]

# Exceptions

PREDICTGENES_INVALID_PYRODIGAL_MODE = ("Invalid execution mode specified for Pyrodigal. Please choose one of the "
                          			   f"following supported execution modes: {PYRODIGAL_MODES}.")

PREDICTGENES_INVALID_GENE_PREDICTOR = ("Specified gene predictor is not valid. Please provide one of the following "
                          			   f"supported gene predictors: {', '.join(GENE_PREDICTORS)}.")

AUGUSTUS_MISSING = ("Could not find AUGUSTUS executables. Please ensure that AUGUSTUS is installed "
                    "and added to the PATH environment variable.")

AUGUSTUS_INVALID_SPECIES = ("Specified species ID is not in the list of supported species for "
                            "AUGUSTUS. Please provide one of the following supported species IDs:\n{0}.")

AUGUSTUS_INVALID_OUTFMT = ("Invalid output format specified for AUGUSTUS results. Please select a single "
                            "or multiple of the following supported output formats: {0}")

PYRODIGAL_INVALID_OUTFMT = ("Invalid output format specified for Pyrodigal results. Please select a single "
                            "or multiple of the following supported output formats: {0}")

INVALID_PTF_PATH = ("Cannot find the specified Prodigal training file in {0}."
					"\nPlease provide a path to a valid training file.")
