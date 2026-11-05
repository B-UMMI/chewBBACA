#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the ComputeMSA module.
"""


# Parameter names
COMPUTEMSA_INPUT_FILES_ARGNAME = "input_path"
COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
COMPUTEMSA_TRANSLATION_TABLE_ARGNAME = "translation_table"
COMPUTEMSA_DNA_MSA_ARGNAME = "dna_msa"
COMPUTEMSA_OUTPUT_VARIABLE_ARGNAME = "output_variable"
COMPUTEMSA_GAPS_ARGNAME = "gaps"
COMPUTEMSA_AMBIGUOUS_ARGNAME = "ambiguous"
COMPUTEMSA_ONLY_LOCI_MSAS_ARGNAME = "only_loci_msas"
COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME = "custom_mafft_options"
COMPUTEMSA_PROTEIN_INPUT_ARGNAME = "protein_input"
COMPUTEMSA_NO_CLEANUP_ARGNAME = "no_cleanup"
COMPUTEMSA_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types
COMPUTEMSA_ARGUMENT_TYPES = {
	"input_path": str,
	"output_directory": str,
	"schema_directory": str,
	"translation_table": int,
	"dna_msa": bool,
	"output_variable": bool,
	"gaps": str,
	"ambiguous": str,
	"only_loci_msas": bool,
	"custom_mafft_options":str,
	"protein_input": bool,
	"no_cleanup": bool,
	"cpu_cores": int
}

# Filenames for the files created by the ComputeMSA module
COMPUTEMSA_PROTEIN_CONCAT = "protein_concat.fasta"
COMPUTEMSA_DNA_CONCAT = "dna_concat.fasta"
COMPUTEMSA_PROTEIN_MSA = "protein_msa.fasta"
COMPUTEMSA_PROTEIN_MSA_VARIABLE = "protein_msa_variable.fasta"
COMPUTEMSA_DNA_MSA = "dna_msa.fasta"
COMPUTEMSA_DNA_MSA_VARIABLE = "dna_msa_variable.fasta"

# Gap character used in the MSAs
GAP_CHAR = "-"
# Default GAP behaviour (exclude all positions/columns with Gaps)
DEFAULT_GAPS = "exclude"
# GAP behaviour options
GAPS_CHOICES = ["ignore", "exclude"]

# Default Ambiguous behaviour (exclude all positions/columns with ambiguous bases)
DEFAULT_AMBIGUOUS = "exclude"
# Ambiguous behaviour options
AMBIGUOUS_CHOICES = ["ignore", "exclude"]

# Exceptions
COMPUTEMSA_INVALID_GAPS_OPTION = (f"Invalid option provided to '--gaps'. Please provide one of the following option: {DEFAULT_GAPS}")

COMPUTEMSA_INVALID_AMBIGUOUS_OPTION = (f"Invalid option provided to '--ambiguous'. Please provide one of the following options: {DEFAULT_AMBIGUOUS}")

COMPUTEMSA_NO_SCHEMA = ('Schema directory must be provided when input is a TSV file with allelic profiles.')

NO_MSAS_CREATED = ('Could not compute the MSA for any of the files. Exiting...')
