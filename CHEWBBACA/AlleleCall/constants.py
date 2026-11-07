#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the AlleleCall module.
"""


# Parameter names
ALLELECALL_INPUT_PATH_ARGNAME = "input_path"
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

# Parameter descriptions
ALLELECALL_INPUT_PATH_DESCRIPTION = ("Path to the directory that contains the input FASTA files or to "
							 "a file with a list of full paths to FASTA files, one per line.")
ALLELECALL_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the schema directory. The schema directory contains the "
							 "loci FASTA files and a folder named `short` that contains the FASTA "
							 "files with the loci representative alleles.")
ALLELECALL_OUTPUT_DIRECTORY_DESCRIPTION = ("Output directory where the process will store intermediate files "
							 "and allele calling results.")
ALLELECALL_LOCI_LIST_DESCRIPTION = ("Path to a file with the list of schema loci to perform allele "
							 "calling for. The file must include the full paths to the loci FASTA "
							 "files or the loci IDs, one per line. The process will perform "
							 "allele calling only for the subset of loci provided in the file.")
ALLELECALL_BLAST_SCORE_RATIO_DESCRIPTION = ("BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp "
							 "alignment and aligned sequences with a BSR >= than the defined value "
							 "are considered to be alleles of the same locus.")
ALLELECALL_MINIMUM_LENGTH_DESCRIPTION = ("Minimum sequence length value. Predicted coding sequences (CDSs) "
							 "shorter than this value are excluded.")
ALLELECALL_TRANSLATION_TABLE_DESCRIPTION = ("Genetic code used to predict genes and to translate coding DNA "
							 "sequences (CDSs). This value will be ignored if a training file "
							 "is used.")
ALLELECALL_SIZE_THRESHOLD_DESCRIPTION = ("Coding sequence (CDS) size variation threshold. At the default "
							 "value of 0.2, CDSs with a size that deviates +-20 percent from "
							 "the locus length mode are classified as ASM/ALM.")
ALLELECALL_NO_INFERRED_DESCRIPTION = ("If provided, the process will not add the sequences of "
							 "inferred alleles (INF) to the schema. Allelic profiles "
							 "will still include the allele identifiers attributed to "
							 "the inferred alleles. Use this parameter if the schema "
							 "is being accessed by multiple processes/users simultaneously.")
ALLELECALL_OUTPUT_UNCLASSIFIED_DESCRIPTION = ("Create a Fasta file with the coding sequences (CDSs) that "
							 "were not classified.")
ALLELECALL_OUTPUT_MISSING_DESCRIPTION = ("Create a Fasta file with coding sequences (CDSs) classified "
							 "as NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5 and LOTSC.")
ALLELECALL_OUTPUT_NOVEL_DESCRIPTION = ("Create a Fasta file with the novel alleles inferred during "
							 "allele calling. The sequence headers include the locus and "
							 "allele identifiers attributed by chewBBACA based on the "
							 "allele calling results.")
ALLELECALL_OUTPUT_MASKED_DESCRIPTION = ("Create a TSV file with the masked allelic profiles. The "
							 "masking process removes the `INF-` prefix from inferred "
							 "alleles and substitutes all special classes (NIPH, NIPHEM, "
							 "ASM, ALM, PLOT3, PLOT5, LOTSC, PAMA) with `0`.")
ALLELECALL_FORCE_CONTINUE_DESCRIPTION = ("If provided, chewie will not warn users and ask for "
							 "permission to continue if any of the provided argument "
							 "values does not match the values in the config file.")
ALLELECALL_EXECUTION_MODE_DESCRIPTION = ("Execution mode (1: only exact matches at DNA level; 2: "
							 "exact matches at DNA and Protein level; 3: exact matches "
							 "and minimizer-based clustering to find similar alleles "
							 "based on BSR+0.1; 4: run the full process to find exact "
							 "matches and similar matches based on BSR value, including "
							 "the determination of new representative alleles to add to "
							 "the schema).")
ALLELECALL_CPU_CORES_DESCRIPTION = ("Number of CPU cores that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores).")
ALLELECALL_NS_CONFIG_DESCRIPTION = ("")

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

