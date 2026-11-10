#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the PredictGenes module.
"""


# Parameter names
PREDICTGENES_INPUT_PATH_ARGNAME = "input_path"
PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
PREDICTGENES_GENE_PREDICTOR_ARGNAME = "gene_predictor"
PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME = "gene_prediction_options"
PREDICTGENES_TRANSLATION_TABLE_ARGNAME = "translation_table"
PREDICTGENES_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types
PREDICTGENES_ARGUMENT_TYPES = {
	"input_path": str,
	"output_directory": str,
	"gene_predictor": str,
	"gene_prediction_options": str,
	"translation_table": int,
	"cpu_cores": int,
}

# Parameter descriptions
PREDICTGENES_INPUT_PATH_DESCRIPTION = ("Path to the directory that contains the input files or to a file "
							 "with a list of full paths to the input files, one per line. Input "
							 "files must be in FASTA format.")
PREDICTGENES_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory where the process will store the files with the predicted genes.")
PREDICTGENES_GENE_PREDICTOR_DESCRIPTION = ("Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")
PREDICTGENES_GENE_PREDICTION_OPTIONS_DESCRIPTION = ("List of arguments passed to configure the gene prediction. When providing "
							 "genome assemblies in FASTA format, the list of arguments for each parameter "
							 "used to configure the gene prediction can be passed as the long format of "
							 "the parameter name followed by the argument value.")
PREDICTGENES_TRANSLATION_TABLE_DESCRIPTION = ("The translation table, matching a specific genetic code, used for gene "
							 "prediction. This value is ignored if a Pyrodigal training file is passed to "
							 "'--gpo, --gene-prediction-options'. In that case, the translation table "
							 "specified in the Pyrodigal training file is used.")
PREDICTGENES_CPU_CORES_DESCRIPTION = ("Number of CPU cores that will be used to run the process.")

# Gene predictor options
GENE_PREDICTORS = ["pyrodigal", "augustus"]
GENE_PREDICTOR_DEFAULT = "pyrodigal"
# List of valid parameter names for each gene predictor
VALID_PARAMETERS = {"augustus": ["augustus_species", "augustus_output_formats", "augustus_path"],
					"pyrodigal": ["pyrodigal_training_file", "pyrodigal_mode",
							  	  "pyrodigal_output_formats", "pyrodigal_minimum_confidence",
							  	  "pyrodigal_training_reference", "pyrodigal_just_training"]}

# Pyrodigal option
PYRODIGAL_MODES = ["single", "meta"]
PYRODIGAL_OUTFMTS = ["genes", "translations", "gff", "genbank", "scores"]
PYRODIGAL_DEFAULT_OUTFMTS = ["genes"]
PYRODIGAL_MIN_CONFIDENCE = 0.0
PYRODIGAL_MAX_CONFIDENCE = 100.0
PYRODIGAL_DEFAULT_MODE = "single"

# AUGUSTUS options
AUGUSTUS_ALIAS = "augustus"
AUGUSTUS_OUTFMTS = ["genes", "gff"]
AUGUSTUS_OUTFMT_DEFAULT = ["genes"]
AUGUSTUS_GFF_FILTERS = ["# start gene", "# coding sequence"]

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

# Could not predict CDSs for input FASTA files
# e.g. files only contain sequence headers, contain invalid
# sequences/chars or pyrodigal cannot predict any genes
PREDICTGENES_NO_PREDICTIONS = ("Could not predict CDSs from any of the input files. Please verify the format of the input files and the argument values passed to chewBBACA.")
