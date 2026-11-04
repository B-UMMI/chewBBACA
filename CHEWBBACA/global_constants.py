#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains default argument values and messages used by chewBBACA's modules.

Code documentation
------------------
"""


import sys
import shutil

##########################################
## Minimum, Maximum, and Default values ##
##########################################

__version__ = "4.0.0-beta"

# Authors, GitHub repository, documentation, tutorial and contacts
AUTHORS = 'Rafael Mamede, Pedro Cerqueira, Mickael Silva, João Carriço, Mário Ramirez'
REPOSITORY = 'https://github.com/B-UMMI/chewBBACA'
DOCUMENTATION = 'https://chewbbaca.readthedocs.io/en/latest/index.html'
CONTACTS = 'imm-bioinfo@medicina.ulisboa.pt'

# BLAST Score Ratio default values
# Value must be contained in the [0.0, 1.0] interval
BSR_MIN = 0.0
BSR_MAX = 1.0
BSR_DEFAULT = 0.6

# Minimum sequence length defaults
MSL_MIN = 0
# Large value to ensure that all sequences above minimum value are accepted
MSL_MAX = sys.maxsize
# Default minimum sequence length value of 201 nucleotides (67 aminoacids)
MSL_DEFAULT = 201

# Allele size variation threshold defaults
ST_MIN = 0.0
ST_MAX = 1.0
# New alleles are inferred if their length value does not deviate more than 
# this value from the locus sequence length mode
# CDSs are classified as ASM or ALM if their length value is below or above
# the computed threshold values, respectively
ST_DEFAULT = 0.2

# Word size/k-mer size value used for minimizer clustering
WORD_SIZE_MIN = 5
WORD_SIZE_MAX = 5
WORD_SIZE_DEFAULT = 5

# Number of adjacent k-mers to consider when selecting minimizers
WINDOW_SIZE_MIN = 5
WINDOW_SIZE_MAX = 5
WINDOW_SIZE_DEFAULT = 5

# Minimum decimal proportion of shared distinct minimizers for a sequence to be added to a cluster
CLUSTERING_SIMILARITY_MIN = 0.20
CLUSTERING_SIMILARITY_MAX = 0.20
CLUSTERING_SIMILARITY_DEFAULT = 0.20

# Decimal proportion of shared distinct minimizers with cluster representative
REPRESENTATIVE_FILTER_MIN = 0.9
REPRESENTATIVE_FILTER_MAX = 0.9
# In the CreateSchema process, clustered sequences are excluded if they share
# this proportion of distinct minimizers with the cluster representative
# This value is not used in the AlleleCall module
REPRESENTATIVE_FILTER_DEFAULT = 0.9

# Decimal proportion of shared distinct minimizers with other clustered sequences
INTRA_CLUSTER_MIN = 0.9
INTRA_CLUSTER_MAX = 0.9
# In the CreateSchema process, clustered sequences are excluded if they share
# this proportion of distinct minimizers with another clustered sequence of
# equal or greater length
# This value is not used in the AlleleCall module
INTRA_CLUSTER_DEFAULT = 0.9

# Genetic codes/translation tables
GENETIC_CODES = {1: "The Standard Code",
				 2: "The Vertebrate Mitochondrial Code",
				 3: "The Yeast Mitochondrial Code",
				 4: "The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
				 5: "The Invertebrate Mitochondrial Code",
				 6: "The Ciliate, Dasycladacean and Hexamita Nuclear Code",
				 9: "The Echinoderm and Flatworm Mitochondrial Code",
				 10: "The Euplotid Nuclear Code",
				 11: "The Bacterial, Archaeal and Plant Plastid Code",
				 12: "The Alternative Yeast Nuclear Code",
				 13: "The Ascidian Mitochondrial Code",
				 14: "The Alternative Flatworm Mitochondrial Code",
				 15: "Blepharisma Nuclear Code",
				 16: "Chlorophycean Mitochondrial Code",
				 21: "Trematode Mitochondrial Code",
				 22: "Scenedesmus obliquus Mitochondrial Code",
				 23: "Thraustochytrium Mitochondrial Code",
				 24: "Rhabdopleuridae Mitochondrial Code",
				 25: "Candidate Division SR1 and Gracilibacteria Code"}

GENETIC_CODE_DEFAULT = 11

CPU_CORES_DEFAULT = 1
CPU_CORES_MIN = 1

# Basenames for the output files
ALLELIC_PROFILES_BASENAME = "allelic_profiles.tsv" # TSV file containing the allelic profiles
COORDINATES_PROFILES_BASENAME = "coordinates_profiles.tsv" # TSV file containing the allele genomic coordinates
PARALOGOUS_COUNTS_BASENAME = "paralogous_counts.tsv" # TSV file containing the list and count per paralogous locus
PARALOGOUS_TABLE_BASENAME = "paralogous_loci.tsv" # TSV file containing each ocurrence of a paralogous locus
ALLELIC_PROFILES_MASKED_BASENAME = "allelic_profiles_masked.tsv" # Masked allelic profiles where all special classes are converted to '0'
SAMPLE_STATISTICS_BASENAME = "sample_stats.tsv" # TSV file containing sample summary statistics
LOCI_STATS_BASENAME = "loci_stats.tsv" # TSV file containing loci summary statistics
UNCLASSIFIED_BASENAME = "unclassified.fasta" # FASTA file containing the unclassified CDSs
MISSING_FASTA_BASENAME = "missing_classes.fasta" # FASTA file containing the CDSs/loci assigned special classifications
MISSING_TSV_BASENAME = "missing_classes.tsv" # TSV file containing information about the CDSs/loci assigned special classifications
INFERRED_ALLELES_BASENAME = "inferred_alleles.fasta" # FASTA file containing the alleles inferred through allele calling
INFERRED_ALLELES_TABLE_BASENAME = "inferred_alleles.tsv" # TSV file containing data about the inferred alleles

# Valid FASTA file extensions
FASTA_EXTENSIONS = [".fasta", ".fna", ".ffn", ".fa", ".fas"]

# Character replacements for genome and loci identifiers
CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""), (")", ""), ("'", ""), ("\"", ""), (":", "")]

# Minimum Python version
MIN_PYTHON = [(3, 10, 0), "3.10.0"]

## GET and POST requests

# Chewie-NS related constants
HEADERS_GET = {"Authorization": None,
			   "accept": "application/octet-stream"}

HEADERS_GET_JSON = {"Authorization": None,
					"accept": "application/json"}

HEADERS_POST = {"Authorization": None,
				"user_id": None}

HEADERS_POST_JSON = {"Authorization": None,
					 "Content-type": "application/json",
					 "accept": "application/json",
					 "user_id": None}

# Timeout, in seconds, to wait for user input
PROMPT_TIMEOUT = 30

# Ambiguous characters for DNA and protein sequences
PROTEIN_AMBIGUOUS_CHARS = ['B', 'Z', 'X', 'J']
DNA_AMBIGUOUS_CHARS = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']

INVALID_CDS_BASENAME = "invalid_cds.txt" # TXT file with the list of invalid CDSs
SCHEMA_CONFIG_BASENAME = ".schema_config" # Pickle file containing schema config values
GENE_LIST_BASENAME = ".genes_list" # Pickle file containing the list of genes in a schema
NS_CONFIG_BASENAME = ".ns_config" # Pickle file added to schemas downloaded to Chewie-NS

# Allele classes
ALLELE_CLASSES = ["EXC", "INF", "PLOT3", "PLOT5", "LOTSC", "NIPH", "NIPHEM", "ALM", "ASM", "PAMA", "LNF"]
# PLNF class for modes {1,2,3}
PROBABLE_LNF = "PLNF"

# String template for a standard single line FASTA record
FASTA_RECORD_TEMPLATE = ">{0}\n{1}"

NUCLEOTIDES = "AGCT"

# Links to Documentation pages
DOCS_URL_TEMPLATE = "https://chewbbaca.readthedocs.io/en/latest/user/modules/{0}.html"

## BLAST

# Minimum MAJOR and MINOR BLAST versions
BLAST_MAJOR = 2
BLAST_MINOR = 9

# Paths to BLASTp and makeblastdb executables in Linux and Windows
BLASTP_ALIAS = "blastp"
MAKEBLASTDB_ALIAS = shutil.which("makeblastdb")
BLASTDB_ALIASTOOL_ALIAS = shutil.which("blastdb_aliastool")
BLASTDBCMD_ALIAS = shutil.which("blastdbcmd")

# Define minimum sequence threshold to change BLASTn and BLASTp task type
BLAST_TASK_THRESHOLD = {"blastn": 50, "blastp": 30}

# BLAST outfmt
BLAST_DEFAULT_OUTFMT = "6 qseqid qstart qend qlen sseqid slen score"

# Define sequence identifier prefix to avoid issues where makeblastdb modifies the IDs
# of the sequences when creating the database (e.g. when IDs are interpreted as PDB chain IDs)
BLASTDB_LCL_PREFIX = "lcl|SEQ"
BLASTDB_SEQ_PREFIX = "SEQ"

## MAFFT

# Path to MAFFT executable
MAFFT_ALIAS = shutil.which("mafft")
# MAFFT default parameters
MAFFT_DEFAULT_PARAMETERS = ["--thread", "1", "--retree", "1", "--maxiterate", "0", "--treeout"]

############################
## Warning and Exceptions ##
############################

MISSING_GENE_LIST = ("Could not find the '.genes_list' file in the schema's directory.")

SCHEMA_LOCI_GENE_LIST_NOMATCH = ("The identifiers of the loci FASTA files in the schema do not match "
                                 "the loci identifiers stored in the '.genes_list' file. The schema "
								 "includes {0} exclusive loci IDs ({1}) and the '.genes_list' file "
								 "includes {2} exclusive loci IDs ({3}).")

LOCI_LIST_SCHEMA_LOCI_NOMATCH = ("The file containing the list of loci includes loci IDs that are not in the "
                                 "schema. Please ensure that all loci IDs in the file match IDs of loci in the "
								 "schema. The following loci are not in the schema: {0}.")

# Input file is a FASTA file but chewBBACA expects a file with a list of
# file paths
FASTA_INPUT_EXCEPTION = ("Input path {0} points to a FASTA file. Please provide a path to "
                         "a folder that contains FASTA files or a file with the list of "
						 "full paths to the FASTA files (one per line).")

# Some input files have an invalid file extension
INVALID_EXTENSION_EXCEPTION = ("The following input files do not have any of the supported "
                               " file extensions:\n{0}\nPlease provide FASTA files with one "
							   "of the following file extensions: {1}.")

# Some of the file paths provided do not exist
MISSING_INPUTS_EXCEPTION = ("Could not find the following input files:\n{0}\n"
                            "Please verify that all input paths exist.")

# Files that do not have the expected format of a FASTA file
NON_FASTA_EXCEPTION = ("The following input files are not valid FASTA files:\n{0}\n."
                       "Please make sure that all provided input files are valid FASTA files.")

MISSING_SCHEMA_FASTAS = ("Input path does not include FASTA files. Please provide "
						 "a valid path for a folder containing only FASTA files ending "
						 f"in one of the following file extensions: {FASTA_EXTENSIONS}")

# Path to schema does not exist
SCHEMA_PATH_MISSING = ("Path to input schema does not exist. Please provide a valid path.")

# User provided legacy schema. Tell user to adapt with the PrepExternalSchema module
SCHEMA_MISSING_FILES = ("Schema does not include a `{0}` file. Please ensure that "
					    "the path to the schema is correct and that it includes all mandatory "
					    "files. If the schema was not created by chewBBACA, you can adapt it "
					    "to make it compatible with chewBBACA using the PrepExternalSchema module.")

# Output directory exists
OUTPUT_DIRECTORY_EXISTS = ("Output directory already exists. Please provide a non-existent path "
                           "to create a new directory where results will be stored.")

# Could not predict CDSs for input FASTA files
# e.g. files only contain sequence headers, contain invalid
# sequences/chars or pyrodigal cannot predict any genes
CANNOT_PREDICT = ('Could not predict CDSs from any of the input files.'
				  '\nPlease verify the format of the input files.')

INVALID_GENETIC_CODE = ("Invalid genetic code value. Value must match one of the following accepted genetic codes"
                        f"codes:\n{'\n'.join([str(k)+': '+v for k, v in GENETIC_CODES.items()])}")

PYTHON_VERSION = ('Python version found: {0}\nPlease use Python >= {1}')

CPU_VALUE_WARNING = ("Provided a CPU core count value that is equal to the number of CPU cores "
                     "in the system {0}. This may affect system responsiveness.")

BLAST_MISSING = ("Could not find the BLAST executables. Please make sure that BLAST "
				 "is installed and added to the PATH environment variable.")

BLAST_INVALID_VERSION = ("Could not determine BLAST version or the version is not valid. Please make "
				 		 f"sure that BLAST>={BLAST_MAJOR}.{BLAST_MINOR} is installed.")

JUST_TRAINING = ('User specified that the process should only create a '
				 'training file (--just-training). Exited.')

ARGS_DIFFER = ('Provided argument values differ from the values '
			   'used for schema creation:\n')

ARGS_DIFFER_PROMPT = ('\nContinuing might lead to results not '
					  'consistent with previous runs.\nProviding '
					  'parameter values that differ from the values '
					  'used for schema creation will also invalidate '
					  'the schema for uploading and synchronization '
					  'with Chewie-NS.\nContinue? (yes/no)\n')

MISSING_CONFIG = ('Could not find a valid config file.')

INPUTS_SHARE_PREFIX = ("The following input files share the same unique filename prefix:\n{0}\n"
                       "Please ensure that every input file has a unique filename prefix. This "
					   "is necessary to unambiguously identify each input during intermediate "
					   "steps and in the output files.")

INPUTS_INCLUDE_BLANKS = ("The following input files include blank spaces in the filename:\n{0}\n"
                         "Please ensure that filenames do not include blank spaces or special "
						 "characters (e.g. !@#?$^*()+). This is necessary to avoid issues related "
						 "to how these characters are recognized by chewBBACA and its dependencies.")

MISSING_INPUT_ARG = ("Input path '{0}' does not exist. Please provide a valid path.")

MISSING_LOCI_LIST = ("Path to the file containing the loci list is invalid.")

ALM_MSG = ('allele greater than {0}% locus length mode ({1}>{2})')
ASM_MSG = ('allele smaller than {0}% locus length mode ({1}<{2})')

CONFLICTING_ARGUMENTS = ("The value provided for the {0} does not match any of the values used with the schema ({1}).")

PYRODIGAL_META_NOPTF = ("Cannot use a training file when running Pyrodigal in `meta` mode. "
						"Please do not provide a training file if setting Pyrodigal's running mode to `meta`.")

CANNOT_PROVIDE_PTF_AND_TREF = ("Cannot provide a Pyrodigal training file and a training reference "
							   "to create a training file. Please provide either a training "
							   "file or a training reference to create a new training file, "
							   "not both.")

MISSING_RESULTS_FILES = (f"Input directory with results data must include the following file types: {0}.")

MISSING_ANNOTATIONS_FILE = ("Path provided for the annotations file {0} is not valid.")

MISSING_RESULTS = ("Path provided for the folder with the allele calling results is not valid.")

PROFILES_FILE_MISSING = ("The path provided to the TSV file with allelic profiles is not valid.")

FORCE_TRAINING_FILE = ("Provided a Pyrodigal training file. Ignoring the translation table value previously set "
                       "({0}) and using the genetic code defined in the training file ({1}).")
