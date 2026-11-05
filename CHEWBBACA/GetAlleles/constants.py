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

# Header for TSV file created by the GetAlleles module
GETALLELES_LOCI_STATS_HEADER = ["Locus", "Total alleles in schema", "Samples with locus", "Distinct alleles in dataset"]
