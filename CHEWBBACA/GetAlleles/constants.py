#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the GetAlleles module.
"""


# Parameter names
GETALLELES_ALLELIC_PROFILES_ARGNAME = "allelic_profiles"
GETALLELES_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
GETALLELES_LOCI_LIST_ARGNAME = "loci_list"
GETALLELES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
GETALLELES_CPU_CORES_ARGNAME = "cpu_cores"
GETALLELES_TRANSLATION_TABLE_ARGNAME = "translation_table"
GETALLELES_TRANSLATE_ALLELES_ARGNAME = "translate_alleles"
GETALLELES_DISTINCT_ARGNAME = "distinct"

# Argument types
GETALLELES_ARGUMENT_TYPES = {
	"allelic_profiles": str,
	"schema_directory": str,
	"loci_list": str,
	"output_directory": str,
	"cpu_cores": int,
	"translation_table": int,
	"translate_alleles": bool,
	"distinct": bool
}

# Parameter descriptions
GETALLELES_ALLELIC_PROFILES_DESCRIPTION = ("Path to the TSV file containing the allelic profiles.")
GETALLELES_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the schema directory.")
GETALLELES_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory.")
GETALLELES_LOCI_LIST_DESCRIPTION = ("Path to a file with the list of genes/loci to create FASTA files for. The "
							 "file must include the identifiers of the loci, one per line, without the "
							 ".fasta extension.")
GETALLELES_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")
GETALLELES_TRANSLATION_TABLE_DESCRIPTION = ("Genetic code used to translate coding DNA sequences (CDSs). If no value"
							 " is specified, the process tries to get the value stored in the schema "
							 "config file. If the schema does not include a config file, the process "
							 "uses the default translation table (11).")
GETALLELES_TRANSLATE_ALLELES_DESCRIPTION = ("Create FASTA files with the translated alleles.")
GETALLELES_DISTINCT_DESCRIPTION = ("Only include the distinct alleles in the output FASTA files.")

# Header for TSV file created by the GetAlleles module
GETALLELES_LOCI_STATS_HEADER = ["Locus", "Total alleles in schema", "Samples with locus", "Distinct alleles in dataset"]
