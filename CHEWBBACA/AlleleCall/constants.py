#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the AlleleCall module.
"""


# Parameter names
ALLELECALL_INPUT_FILES_ARGNAME = "input_path"
ALLELECALL_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
ALLELECALL_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
ALLELECALL_LOCI_LIST_ARGNAME = "loci_list"
ALLELECALL_BLAST_SCORE_RATIO_ARGNAME = "blast_score_ratio"
ALLELECALL_MINIMUM_LENGTH_ARGNAME = "minimum_length"
ALLELECALL_TRANSLATION_TABLE_ARGNAME = "translation_table"
ALLELECALL_SIZE_THRESHOLD_ARGNAME = "size_threshold"
ALLELECALL_GENE_PREDICTOR_ARGNAME = "gene_predictor"
ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME = "gene_prediction_options"
ALLELECALL_CLUSTERING_OPTIONS_ARGNAME = "clustering_options"
ALLELECALL_BLAST_PATH_ARGNAME = "blast_path"
ALLELECALL_CDS_INPUT_ARGNAME = "cds_input"
ALLELECALL_NO_INFERRED_ARGNAME = "no_inferred"
ALLELECALL_OUTPUT_UNCLASSIFIED_ARGNAME = "output_unclassified"
ALLELECALL_OUTPUT_MISSING_ARGNAME = "output_missing"
ALLELECALL_OUTPUT_NOVEL_ARGNAME = "output_novel"
ALLELECALL_OUTPUT_MASKED_ARGNAME = "output_masked_profiles"
ALLELECALL_FORCE_CONTINUE_ARGNAME = "force_continue"
ALLELECALL_EXECUTION_MODE_ARGNAME = "execution_mode"
ALLELECALL_CPU_CORES_ARGNAME = "cpu_cores"
ALLELECALL_NO_CDS_RENAMING_ARGNAME = "no_cds_renaming"
ALLELECALL_NO_CLEANUP_ARGNAME = "no_cleanup"

# Argument types
ALLELECALL_ARGUMENT_TYPES = {
	"input_path": str,
	"schema_directory": str,
	"output_directory": str,
	"loci_list": str,
	"blast_score_ratio": float,
	"minimum_length": int,
	"translation_table": int,
	"size_threshold": float,
	"gene_predictor": str,
	"gene_prediction_options": str,
	"clustering_options": str,
	"blast_path": str,
	"cds_input": bool,
	"no_inferred": bool,
	"output_unclassified": bool,
	"output_missing": bool,
	"output_novel": bool,
	"output_masked": bool,
	"force_continue": bool,
	"execution_mode": int,
	"cpu_cores": int,
	"no_cds_renaming": bool,
	"no_cleanup": bool
}

# Headers used in output files
# Headers for TSV files with paralogous loci count and per genome
PARALOGOUS_COUNTS_HEADER = ["Locus", "Count"]
PARALOGOUS_TABLE_HEADER = ["Genome", "Loci", "CDS"]
# Header for TSV file with loci stats
LOCI_STATS_HEADER = ["Locus", "EXC", "INF", "PLOT3", "PLOT5", "LOTSC", "NIPH",
					 "NIPHEM", "ALM", "ASM", "PAMA", "LNF", "Total CDSs Classified"]
INFERRED_ALLELES_TABLE_HEADER = ["Allele ID", "CDS ID", "Representative"]
# Header for TSV file with information about CDSs classified as ambiguous
MISSING_TSV_HEADER = ["Index", "Genome", "Locus", "Locus_classification", "CDS", "CDS_classification"]

# Maximum number of values stored while creating the 'allelic_profiles.tsv' and 'coordinates_profiles.tsv' files
MAX_STORED_VALUES = 300000

# Maximum number of allele hashes per pre-computed file
HASH_TABLE_MAXIMUM_ALLELES = 200000

# Dictionary template to map variables returned by the AlleleCall module
ALLELECALL_DICT = {"classification_files": None,
				   "basename_map": None,
				   "cds_coordinates": None,
				   "cds_counts": None,
				   "dna_fasta": None,
				   "protein_fasta": None,
				   "dna_hashtable": None,
				   "protein_hashtable": None,
				   "invalid_alleles": None,
				   "unclassified_ids": None,
				   "self_scores": None,
				   "representatives": None}

# AlleleCall section headers
PRECOMPUTED_DATA = "Pre-computed data"
CDS_PREDICTION = "CDS prediction"
CDS_DEDUPLICATION = "CDS deduplication"
CDS_EXACT = "CDS exact matching"
CDS_TRANSLATION = "CDS translation"
PROTEIN_DEDUPLICATION = "Protein deduplication"
PROTEIN_EXACT = "Protein exact matching"
PROTEIN_CLUSTERING = "Protein clustering"
REPRESENTATIVE_DETERMINATION = "Representative determination"
WRAPPING_UP = "Wrapping up"

# Execution modes
ALLELECALL_MODES = {1: "Searches for exact matches at the DNA level between the schema alleles and the input CDSs.",
                    2: "Searches for exact matches at the DNA and protein levels between the schema alleles and the input CDSs.",
					3: "Searches for exact matches at the DNA and protein levels and clusters CDSs against all representative alleles followed by intracluster BLASTp.",
					4: "The default mode. Runs the whole process (exact matching, clustering followed by intracluster BLASTp, and a final all-vs-all BLASTp between all representative alleles and the remaining unclassified CDSs)."}

ALLELECALL_DEFAULT_MODE = 4

# Exceptions
ALLELECALL_INVALID_MODE = ("Specified invalid mode for allele calling. Please provide one of the "
                           f"following supported modes: {'\n'.join([str(k)+': '+v for k, v in ALLELECALL_MODES.items()])}")

