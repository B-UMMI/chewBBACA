#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the ComputeMSA module.
"""


# Parameter names
COMPUTEMSA_INPUT_PATH_ARGNAME = "input_path"
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

# Parameter descriptions
COMPUTEMSA_INPUT_PATH_DESCRIPTION = ("Path to a TSV file containing allelic profiles or to a folder containing "
							 "FASTA files. If a TSV file containing allelic profiles is provided, it is "
							 "necessary to provide the path to the schema to the `--schema-directory` "
							 "parameter. The module will create a FASTA file with the alleles identified "
							 "in the samples for each schema locus and compute a MSA. The loci MSAs are "
							 "joined to create the complete MSA based on the allele calling results. If "
							 "a path to a folder is provided, the module computes a MSA for each FASTA "
							 "file in the folder, but will not attempt to join the MSAs as it does not "
							 "have the sample information (in this case, it is not necessary to pass the "
							 "schema path).")
COMPUTEMSA_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory where the process will store intermediate and "
							 "final results.")
COMPUTEMSA_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the schema\'s directory. This parameter is only required if the "
							 "input is a TSV file with allelic profiles.")
COMPUTEMSA_TRANSLATION_TABLE_DESCRIPTION = ("Genetic code used for sequence translation.")
COMPUTEMSA_DNA_MSA_DESCRIPTION = ("Converts the protein MSA back to DNA to create an additional output file "
							 "with the DNA MSA.")
COMPUTEMSA_OUTPUT_VARIABLE_DESCRIPTION = ("Output a reduced MSA including only the variable positions. If the "
							 "`--dna-msa` parameter is provided, the process will output a reduced MSA "
							 "for both the protein and DNA MSAs.")
COMPUTEMSA_GAPS_DESCRIPTION = ("How to treat gaps when determining the reduced MSA for the variable "
							 "positions. The default value, `exclude`, removes variable positions if "
							 "any of the aligned sequences contain a gap. The `ignore` option allows "
							 "to consider variable positions that include gaps in some sequences as "
							 "long as other sequences include variable non-gap characters. The character "
							 "used to represent gaps is `-`.")
COMPUTEMSA_AMBIGUOUS_DESCRIPTION = ("How to treat ambiguous amino acids or nucleotides when determining the "
							 "reduced MSA for the variable positions. The default value, `exclude`, "
							 "removes variable positions if any of the aligned sequences contain an "
							 "ambiguous amino acid or nucleotide. The `ignore` option allows to consider "
							 "variable positions that include ambiguous amino acids or nucleotides in "
							 "some sequences as long as other sequences include variable non-ambiguous "
							 "characters. The characters interpreted as ambiguous amino acids are "
							 "[B, Z, X, J]. The characters interpreted as ambiguous nucleotides are"
                             " [R, Y, S, W, K, M, B, D, H, V, N].")
COMPUTEMSA_ONLY_LOCI_MSAS_DESCRIPTION = ("Do not compute the full MSA when the input file is a TSV file containing "
							 "allelic profiles (this is already the default when the input is a path to "
							 "a folder with FASTA files).")
COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_DESCRIPTION = ("Custom parameters to pass to MAFFT when computing the loci MSAs. The "
							 "value must be a single string with all parameters enclosed in quotes "
							 "(e.g. '--retree 1 --maxiterate 0').")
COMPUTEMSA_PROTEIN_INPUT_DESCRIPTION = ("Input files contain protein sequences. This option skips sequence translation "
							 "and is only valid for cases when users provide a path to a directory containing "
							 "FASTA files.")
COMPUTEMSA_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

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
