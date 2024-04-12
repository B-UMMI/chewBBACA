#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains default values for chewBBACA's
parameters.

Code documentation
------------------
"""


import sys
import shutil
import platform


# BLAST Score Ratio default values
# Value must be contained in the [0.0,1.0] interval
BSR_MIN = 0.0
BSR_MAX = 1.0
DEFAULT_BSR = 0.6

# Minimum sequence length defaults
MSL_MIN = 0
# Large value to ensure that all sequences above minimum value are accepted
MSL_MAX = sys.maxsize
# Default minimum sequence length value of 201 nucleotides (67 aminoacids)
MINIMUM_LENGTH_DEFAULT = 201

# Size variation threshold defaults
ST_MIN = 0.0
ST_MAX = 1.0
# New alleles are inferred if their length value does
# not deviate more than this value from the locus sequence
# length mode
# ASM if below threshold and ALM if above
SIZE_THRESHOLD_DEFAULT = 0.2

# Word size/k value used for minimizer clustering
# this value should not be modified
WORD_SIZE_MIN = 5
WORD_SIZE_MAX = 5
WORD_SIZE_DEFAULT = 5

# Number of adjacent kmers to consider when selecting minimizers
# this value should not be modified
WINDOW_SIZE_MIN = 5
WINDOW_SIZE_MAX = 5
WINDOW_SIZE_DEFAULT = 5

# Minimum decimal proportion of shared distinct minimizers for
# a sequence to be added to a cluster
# this value should not be modified
CLUSTERING_SIMILARITY_MIN = 0.20
CLUSTERING_SIMILARITY_MAX = 0.20
CLUSTERING_SIMILARITY_DEFAULT = 0.20

# Decimal proportion of shared distinct minimizers with cluster
# representative
REPRESENTATIVE_FILTER_MIN = 0.9
REPRESENTATIVE_FILTER_MAX = 0.9
# In the CreateSchema process, clustered sequences are excluded
# if they share this proportion of distinct minimizers with the
# cluster representative
REPRESENTATIVE_FILTER_DEFAULT = 0.9

# Decimal proportion of shared distinct minimizers with other
# clustered sequences
INTRA_CLUSTER_MIN = 0.9
INTRA_CLUSTER_MAX = 0.9
# In the CreateSchema process, clustered sequences are excluded
# if they share this proportion of distinct minimizers with another
# clustered sequence of equal or greater length
INTRA_CLUSTER_DEFAULT = 0.9

# Genetic codes/translation tables
GENETIC_CODES = {1: 'Standard',
                 4: 'The mold, protozoan, and coelenterate mitochondrial '
                    'code and the mycoplasma/spiroplasma code',
                 11: 'The Bacterial, Archaeal and Plant Plastid code',
                 25: 'Candidate division SR1 and gracilibacteria code'}

GENETIC_CODES_DEFAULT = 11

# Proteins to cluster are divided into a maximum
# of 40 smaller groups in CreateSchema
# Dividing based on the number of CPU cores can lead to
# variable results because we do not have pre-defined clusters
# in CreateSchema.
CREATESCHEMA_CLUSTERING_NGROUPS = 40

# Valid FASTA file extensions
FASTA_EXTENSIONS = ['.fasta', '.fna', '.ffn', '.fa', '.fas']

# Chewie-NS related constants
HEADERS_GET_ = {'Authorization': None,
                'accept': 'application/octet-stream'}

HEADERS_GET_JSON = {'Authorization': None,
                    'accept': 'application/json'}

HEADERS_POST = {'Authorization': None,
                'user_id': None}

HEADERS_POST_JSON = {'Authorization': None,
                     'Content-type': 'application/json',
                     'accept': 'application/json',
                     'user_id': None}

# List of Chewie-NS instance identifiers and URLs
HOST_NS = {'main': 'https://chewbbaca.online/api/NS/api/',
           'tutorial': 'https://tutorial.chewbbaca.online/api/NS/api/',
           'local': 'http://127.0.0.1:5000/NS/api/'}

# Authors, GitHub repository, documentation, tutorial and contacts
authors = 'Rafael Mamede, Pedro Cerqueira, Mickael Silva, João Carriço, Mário Ramirez'
repository = 'https://github.com/B-UMMI/chewBBACA'
documentation = 'https://chewbbaca.readthedocs.io/en/latest/index.html'
contacts = 'imm-bioinfo@medicina.ulisboa.pt'

# Timeout, in seconds, to wait for user input
prompt_timeout = 30

# Minimum MAJOR and MINOR BLAST versions
BLAST_MAJOR = 2
BLAST_MINOR = 9

# Paths to BLASTp and makeblastdb executables in Linux and Windows
BLASTP_ALIAS = 'blastp.exe' if platform.system() == 'Windows' else shutil.which('blastp')
MAKEBLASTDB_ALIAS = 'makeblastdb.exe' if platform.system() == 'Windows' else shutil.which('makeblastdb')
BLASTDB_ALIASTOOL_ALIAS = 'blastdb_aliastool.exe' if platform.system() == 'Windows' else shutil.which('blastdb_aliastool')

# BLAST warnings to be ignored
# This warning is raised in BLAST>=2.10 when passing a TXT file with sequence identifiers to -seqidlist
# To avoid this warning the TXT file must be converted to binary with the blastdb_aliastool
# Performance can be severely affected if the TXT is not converted to binary
IGNORE_RAISED = ['Warning: [blastp] To obtain better run time '
                 'performance, please run blastdb_aliastool '
                 '-seqid_file_in <INPUT_FILE_NAME> -seqid_file_out '
                 '<OUT_FILE_NAME> and use <OUT_FILE_NAME> as the '
                 'argument to -seqidlist']

# Path to MAFFT executable
MAFFT_ALIAS = shutil.which('mafft')

# Replacements for genome and loci identifiers
CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""),
                     (")", ""), ("'", ""), ("\"", ""),
                     (":", "")]

# Minimum Python version
MIN_PYTHON = [(3, 6, 0), '3.6.0']

# UniProt SPARQL endpoint
UNIPROT_SPARQL = 'https://sparql.uniprot.org/sparql'
UNIPROT_SPARQL_THREADS = 4
# Maximum number of retries if querying the SPARQL endpoint fails
MAX_RETRIES = 2
# Maximum number of sequences used to query the SPARQL endpoint
MAX_QUERIES = 20

# FTP to get UniProt's reference proteomes
UNIPROT_PROTEOMES_FTP = ('https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/')

# List of UniProt's uninformative terms
UNIPROT_UNINFORMATIVE = ['uncharacterized', 'hypothetical', 'duf']

# AlleleCall logfile basename and content template
LOGFILE_BASENAME = 'logging_info.txt'
LOGFILE_TEMPLATE = ('Started script at: {0}\n'
                    'Finished script at: {1}\n'
                    'Number of inputs: {2}\n'
                    'Number of loci: {3}\n'
                    'Used this number of CPU cores: {4}\n'
                    'Used a BSR of: {5}\n')

# Basename for files created by the AlleleCall module
RESULTS_COORDINATES_BASENAME = 'results_contigsInfo.tsv'
PARALOGOUS_COUNTS_BASENAME = 'paralogous_counts.tsv'
PARALOGOUS_LOCI_BASENAME = 'paralogous_loci.tsv'
RESULTS_ALLELES_BASENAME = 'results_alleles.tsv'
RESULTS_STATISTICS_BASENAME = 'results_statistics.tsv'
LOCI_STATS_BASENAME = 'loci_summary_stats.tsv'
UNCLASSIFIED_BASENAME = 'unclassified_sequences.fasta'
MISSING_FASTA_BASENAME = 'missing_classes.fasta'
MISSING_TSV_BASENAME = 'missing_classes.tsv'
NOVEL_BASENAME = 'novel_alleles.fasta'
CDS_COORDINATES_BASENAME = 'cds_coordinates.tsv'
INVALID_CDS_BASENAME = 'invalid_cds.txt'
SCHEMA_CONFIG_BASENAME = '.schema_config'
GENE_LIST_BASENAME = '.genes_list'
# Header for TSV file with loci stats
LOCI_STATS_HEADER = ('Locus\tEXC\tINF\tPLOT3\tPLOT5\tLOTSC\tNIPH\t'
                     'NIPHEM\tALM\tASM\tPAMA\tLNF\tTotal_CDS')
# Header for TSV file with information about extracted CDSs
CDS_TABLE_HEADER = 'Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n'
# Headers for TSV files with paralogous loci count and per genome
PARALOGOUS_COUNTS_HEADER = 'Locus\tCount'
PARALOGOUS_LIST_HEADER = 'Genome\tLoci\tCDS'
# Header for TSV file with information about CDSs classified as ambiguous
MISSING_HEADER = 'Index\tGenome\tLocus\tLocus_classification\tCDS\tCDS_classification'

# Allele calling classifications
ALLELECALL_CLASSIFICATIONS = ['EXC', 'INF', 'PLOT3', 'PLOT5',
                              'LOTSC', 'NIPH', 'NIPHEM', 'ALM',
                              'ASM', 'PAMA', 'LNF']

# PLNF classificaton for modes {1,2,3}
PROBABLE_LNF = 'PLNF'

# Maximum number of values stored while creating the 'results_contigsInfo.tsv' file
RESULTS_MAXVALS = 300000

# String template for a standard single line FASTA record
FASTA_RECORD_TEMPLATE = '>{0}\n{1}'
FASTA_CDS_TEMPLATE = '>{0}-protein{1}\n{2}'

DNA_BASES = 'AGCT'

# Define default BLASTp task
BLAST_TASK_THRESHOLD = {'blastn': 50, 'blastp': 30}

# BLAST outfmt
BLAST_DEFAULT_OUTFMT = '6 qseqid qstart qend qlen sseqid slen score'

# Input file prefix maximum length
PREFIX_MAXLEN = 30

# Dictionary template to map variables returned by AlleleCall
ALLELECALL_DICT = {'classification_files': None,
                   'basename_map': None,
                   'cds_coordinates': None,
                   'cds_counts': None,
                   'dna_fasta': None,
                   'protein_fasta': None,
                   'dna_hashtable': None,
                   'protein_hashtable': None,
                   'invalid_alleles': None,
                   'unclassified_ids': None,
                   'self_scores': None,
                   'representatives': None}

GENOME_LIST = 'listGenomes2Call.txt'
LOCI_LIST = 'listGenes2Call.txt'

# Maximum number of allele hashes per pre-computed file
HASH_TABLE_MAXIMUM_ALLELES = 200000

# AlleleCall section headers
CONFIG_VALUES = 'Configuration values'
PRECOMPUTED_DATA = 'Pre-computed data'
CDS_PREDICTION = 'CDS prediction'
CDS_DEDUPLICATION = 'CDS deduplication'
CDS_EXACT = 'CDS exact matching'
CDS_TRANSLATION = 'CDS translation'
PROTEIN_DEDUPLICATION = 'Protein deduplication'
PROTEIN_EXACT = 'Protein exact matching'
PROTEIN_CLUSTERING = 'Protein clustering'
REPRESENTATIVE_DETERMINATION = 'Representative determination'
WRAPPING_UP = 'Wrapping up'

# CreateSchema exclusive section headers
EXCLUDE_SMALL = 'Short CDS removal'
FINAL_BLASTp = 'Final BLASTp'

# File header for file with summary statistics created by PrepExternalSchema
PREPEXTERNAL_SUMMARY_STATS_HEADER = ('Gene\tTotal_alleles\tValid_alleles\t'
                                     'Number_representatives')

# Default loci presence thresholds used to compute the cgMLST
CGMLST_THRESHOLDS = [0.95, 0.99, 1]

# HTML template to create Schema Report
# need to include '.' at start to work properly when referencing local files
SCHEMA_REPORT_HTML = ("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>Schema Evaluator - React Edition</title>
</head>
<body style="background-color: #f6f6f6">
    <noscript> You need to enable JavaScript to run this app. </noscript>
    <div id="root"></div>
    <script> preComputedData = {0} </script>
    <script src="./report_bundle.js"></script>
</body>
</html>
""")

# HTML template to create Loci Reports
# need to include '.' at start to work properly when referencing local files
LOCUS_REPORT_HTML = ("""
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="UTF-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        <title>Schema Evaluator - Individual Analysis</title>
    </head>
    <body style="background-color: #f6f6f6">
        <noscript> You need to enable JavaScript to run this app. </noscript>
        <div id="root"></div>
        <script src="https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.js"></script>
        <link type=text/css rel=stylesheet href=https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.css />
        <script> preComputedDataInd = {0} </script>
        <script src="./report_bundle.js"></script>
    </body>
</html>
""")

# HTML template to create main AlleleCall report
ALLELECALL_REPORT_HTML = ("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>AlleleCall Report - React Edition</title>
</head>
<body style="background-color: #f6f6f6">
    <noscript> You need to enable JavaScript to run this app. </noscript>
    <div id="root"></div>
    <script> preComputedData = {0} </script>
    <script src="./report_bundle.js"></script>
</body>
</html>
""")

# Basename for files created by the SchemaEvaluator module
SCHEMA_REPORT_BASENAME = 'schema_report.html'

# Basename for files created by the AlleleCallEvaluator module
DISTANCE_MATRIX_BASENAME = 'distance_matrix.tsv'
CORE_MSA_BASENAME = 'cgMLST_MSA.fasta'
MASKED_PROFILES_BASENAME = 'masked_profiles.tsv'
CGMLST_PROFILES_BASENAME = 'cgMLST_profiles.tsv'
ALLELECALL_REPORT_BASENAME = 'allelecall_report.html'

# Basename for files created by ExtractCgMLST module
PRESENCE_ABSENCE_BASENAME = 'presence_absence.tsv'
MISSING_LOCI_BASENAME = 'missing_loci_stats.tsv'

# Relative path to the JS bundles used by SchemaEvaluator
# Main page
SCHEMA_EVALUATOR_SCHEMA_BUNDLE = 'report_template_components/src/bundles/SchemaEvaluator/schema_report/report_bundle.js'
# Loci pages
SCHEMA_EVALUATOR_LOCI_BUNDLE = 'report_template_components/src/bundles/SchemaEvaluator/loci_reports/report_bundle.js'

# Relative path to the JS bundle used by AlleleCallEvaluator
ALLELECALL_EVALUATOR_BUNDLE = 'report_template_components/src/bundles/AlleleCallEvaluator/report_bundle.js'

# Do not use list of strings as constants if the strings include formatting
# placeholders. Multiple references to the list of strings will have the same
# id and altering the strings with format will not change the list id. In
# multiprocessing it can reference the same list/id in different processess
# and use the latest changes to a string in the list/id when those changes
# might not refer to the current process (returning an incorrect value).

# Table header for the Schema report summary data
SCHEMA_SUMMARY_TABLE_HEADERS = ('Loci\tAlleles\tValid Alleles\tInvalid '
                                'alleles\tIncomplete ORF\tAmbiguous Bases\t'
                                'Missing Start/Stop Codon\tIn-frame Stop '
                                'Codon\tAlleles < {0}bp\tAlleles below '
                                'threshold\tAlleles above threshold')

# Column headers for the Loci Analysis Table in the Schema report
LOCI_ANALYSIS_COLUMNS = ('Locus\tTotal Alleles\tValid Alleles\tInvalid '
                         'Alleles\tProportion of Validated Alleles\tDistinct '
                         'Protein Alleles\tIncomplete '
                         'ORF\tAmbiguous Bases\tMissing '
                         'Start/Stop Codon\tIn-frame Stop Codon\tAlleles '
                         '< {0}bp\tAlleles below threshold\tAlleles above '
                         'threshold\tMissing Allele IDs')

# Column headers for the Summary Table in the Loci reports
LOCUS_COLUMNS = ('Locus\tTotal Alleles\tValid Alleles\tInvalid '
                 'Alleles\tProportion of Validated Alleles\tDistinct '
                 'Protein Alleles\t'
                 'Incomplete ORF\tAmbiguous Bases\tMissing Start/Stop '
                 'Codon\tIn-frame Stop Codon\tAlleles < {0}bp\tSize Range '
                 '(bp)\tLength Median (bp)\tLength Mode (bp)\tAlleles below '
                 'threshold ({1}bp)\tAlleles above threshold ({2}bp)\t'
                 'Missing Allele IDs')

# Column headers for the Invalid Alleles table in the loci reports
INVALID_ALLELES_COLUMNS = ['Allele ID', 'Exception Category',
                           'Exception Description']

TRANSLATION_EXCEPTIONS = ['Extra in frame stop codon',
                          'is not a start codon',
                          'is not a stop codon',
                          'sequence length is not a multiple of 3',
                          'ambiguous or invalid characters']

DISTINCT_ALLELES_COLUMNS = ['Protein Allele ID',
                            'Count',
                            'List of Distinct Alleles']

# Column headers for the Sample Stats table in the allele calling report
SAMPLE_STATS_COLUMNS = ['Sample', 'Total Contigs', 'Total CDSs',
                        'Proportion of Classified CDSs', 'Identified Loci',
                        'Proportion of Identified Loci',
                        'Valid Classifications', 'Invalid Classifications']

# Column headers for the Loci Stats table in the allele calling report
LOCI_STATS_COLUMNS = ['Locus', 'Total CDSs', 'Valid Classes',
                      'Invalid Classes', 'Proportion Samples']

# Column headers for the Summary Stats table in the allele calling report
SUMMARY_STATS_COLUMNS = ['Total Samples', 'Total Loci', 'Total CDSs',
                         'Total CDSs Classified', 'EXC', 'INF',
                         'PLOT3', 'PLOT5', 'LOTSC', 'NIPH',
                         'NIPHEM', 'ALM', 'ASM', 'PAMA', 'LNF']

# Exception messages

# Input file is a FASTA file but chewBBACA expects a file with a list of
# file paths
FASTA_INPUT_EXCEPTION = ('Input file is a FASTA file. Please provide '
                         'the path to the parent directory that contains '
                         'the FASTA files or a file with the list of full '
                         'paths to the FASTA files (one per line).')

# Some input files have an invalid file extension
INVALID_EXTENSION_EXCEPTION = ('The following input files do not have a '
                               'valid file extension:\n{0}\nPlease ensure '
                               'that the filenames end with one of the '
                               f'following extensions: {FASTA_EXTENSIONS}.')

# Some of the file paths provided do not exist
MISSING_INPUTS_EXCEPTION = ('Could not find some of the files provided in '
                            'the input list. Please verify that you\'ve '
                            'provided valid paths to the following input '
                            'files.\n{0}')

# Files that do not have the expected format of a FASTA file
NON_FASTA_EXCEPTION = ('The following input files are not in FASTA format:\n{0}')

# Input directory does not contain FASTA files
MISSING_FASTAS_EXCEPTION = ('Could not get input files. Please provide '
                            'a directory with FASTA files or a file with '
                            'the list of full paths to the FASTA files '
                            'and ensure that filenames end with one of '
                            f'the following extensions: {FASTA_EXTENSIONS}.')

# Input path is neither a file nor a directory
INVALID_INPUT_PATH = ('Input argument is not a valid directory or '
                      'file with a list of paths to FASTA files. Please '
                      'provide a valid input, either a folder with FASTA '
                      'files or a file with the list of full paths to FASTA '
                      'files (one per line and ending with one of the '
                      f'following file extensions: {FASTA_EXTENSIONS}).')

# Path to schema does not exist
SCHEMA_PATH_MISSING = ('Path to input schema does not exist. Please provide '
                       'a valid path.')

# Path to schema does not include expected files
SCHEMA_INVALID_PATH = ('Provided path does not include all the necessary '
                       'schema files. Please verify that you have passed '
                       'the correct path to the schema.')

# User provided legacy schema. Tell user to adapt with the PrepExternalSchema module
ADAPT_LEGACY_SCHEMA = ('Schema does not include a config file. Probably because '
                       'it was created with chewBBACA<=2.1.0. Please adapt schema '
                       'with the PrepExternalSchema module to add a config file.')

# Output directory exists
OUTPUT_DIRECTORY_EXISTS = ('Output directory already exists. Please '
                           'provide a path to a directory that will be '
                           'created to store the results.')

# Input file is a FASTA file but chewBBACA expects a file with a list of
# loci file paths
FASTA_LOCI_LIST_EXCEPTION = ('Path provided to --gl is for a FASTA file. Please provide '
                             'a file with the list of locus identifiers or full '
                             'paths to the loci FASTA files (one per line).')

# Invalid paths to loci FASTA files
MISSING_LOCI_EXCEPTION = ('Could not find some of the loci FASTA files provided in '
                          'the input list. Please verify that you\'ve '
                          'provided valid paths to the following input '
                          'files.\n{0}')

# Invalid format for loci files
NON_FASTA_LOCI_EXCEPTION = ('The following loci files are not in FASTA format:\n{0}')

# User does not have permissions to upload schemas to Chewie-NS
LOADSCHEMA_NO_PERMISSIONS = ('Current user has no Administrator or Contributor '
                             'permissions.\nNot allowed to upload schemas.')

# PTF is missing from schema's directory
LOADSCHEMA_MISSING_PTF = ('Please ensure that the schema\'s directory includes the '
                          'Prodigal training file used to create the schema.')

# Path for PTF does not exist
INVALID_PTF_PATH = 'Invalid path for Prodigal training file.'

# Could not predict CDSs for input FASTA files
# e.g. files only contain sequence headers, contain invalid
# sequences/chars or pyrodigal cannot predict any genes
CANNOT_PREDICT = ('Could not predict CDSs from any of the input files.'
				  '\nPlease provide input files in the accepted FASTA format.')

INVALID_BSR = ('\nBSR value is not contained in the [0.0, 1.0] interval.')
INVALID_BSR_TYPE = ('\nInvalid BSR value of {0}. BSR value must be contained in the [0.0, 1.0] interval.')

INVALID_MINLEN = ('\nInvalid minimum sequence length value. Must be equal or greater than 0.')
INVALID_MINLEN_TYPE = ('\nInvalid minimum sequence length value. Value must be a positive integer.')

INVALID_ST = ('\nInvalid size threshold value. Must be contained in the [0.0, 1.0] interval.')
INVALID_ST_TYPE = ('\nInvalid size threshold value used to create schema. Value must be None or a positive float in the [0.0, 1.0] interval.')

INVALID_GENETIC_CODE = ('\nInvalid genetic code value.\nValue must correspond to '
                 		'one of the accepted genetic codes\n\nAccepted genetic '
                 		'codes:\n{0}')

INVALID_WS = ('\nWord size for the clustering step '
                     'must be equal or greater than {0} and '
                     'equal or smaller than {1}.')
INVALID_WS_TYPE = ('\nSchema created with invalid clustering word size value.')

INVALID_CS = ('\nClustering similarity threshold value '
              'must be contained in the [0.0, 1.0] '
              'interval.')
INVALID_CS_TYPE = ('\nSchema created with invalid clustering threshold value.')

INVALID_RF = ('\nRepresentative filter threshold value '
              'must be contained in the [0.0, 1.0] '
              'interval.')
INVALID_RF_TYPE = ('\nSchema created with invalid representative filter value.')

INVALID_ICF = ('\nIntra-cluster filter value '
               'must be contained in the [0.0, 1.0] '
               'interval.')
INVALID_ICF_TYPE = ('\nSchema created with invalid intra-cluster filter value.')

NS_CANNOT_CONNECT = ('Failed to establish a connection to the Chewie-NS instance at {0}.')

PYTHON_VERSION = ('Python version found: {0}\nPlease use Python >= {1}')

CPU_RESET_WARNING = ('Warning! You have provided a CPU core count value '
             		 'that is equal to or exceeds the number of CPU '
             		 'cores in your system! Resetting to: {0}')
CPU_VALUE_WARNING = ('Warning! You have provided a CPU core count value '
              		 'that is close to the maximum core count of your '
              		 'machine ({0}/{1}). This may affect your system '
              		 'responsiveness.')

BLAST_NO_PATH = ('Could not find BLAST executables.')
BLAST_NO_VERSION = ('Could not determine BLAST version. Please make '
                 	'sure that BLAST>={0}.{1} is installed.')
BLAST_UPDATE = ('Found BLAST {0}.{1}. Please update BLAST to version >={2}.{3}')

MULTIPLE_PTFS = ('Found more than one Prodigal training '
                 'file in the schema directory.\nPlease maintain '
                 'only the training file used in the schema '
                 'creation process.')
MISSING_PTF = ('Could not find a Prodigal training file in the schema directory.')

INVALID_PTF_PATH = ('Cannot find specified Prodigal training file.'
                    '\nPlease provide a valid training file.\nYou '
                    'can create a training file for a species of '
                    'interest with the following command:\n\n  prodigal '
                    '-i <reference_genome> -t <training_file.trn> -p '
                    'single\n\nIt is strongly advised to provide a '
                    'high-quality and closed genome for the training '
                    'process.')

DIFFERENT_PTF_PROMPT = ('Prodigal training file is not the one '
              			'used to create the schema. Using this training '
			  			'file might lead to results not consistent with '
			  			'previous runs and invalidate the schema for '
			  			'usage with Chewie-NS.\nContinue process?\n')
MULTIPLE_PTF_PROMPT = ('Prodigal training file is not any of the {0} '
                       'used in previous runs.\nContinue?\n')

ARGS_DIFFER = ('Provided argument values differ from the values '
			   'used for schema creation:\n')

ARGS_DIFFER_PROMPT = ('\nContinuing might lead to results not '
                      'consistent with previous runs.\nProviding '
                      'parameter values that differ from the values '
                      'used for schema creation will also invalidate '
                      'the schema for uploading and synchronization '
                      'with Chewie-NS.\nContinue? (yes/no)\n')

MISSING_CONFIG = ('Could not find a valid config file.')

INPUTS_SHARE_PREFIX = ('The following input files share the same filename prefix '
					   '(substring before the first "." in the filename):\n{0}\n'
					   'Please ensure that every input file has a unique '
					   'filename prefix.')

INPUTS_INCLUDE_BLANKS = ('The following input files include blank spaces '
						 'in the filename:\n{0}\nPlease ensure that filenames '
						 'do not include blank spaces or special characters '
						 '(e.g. !@#?$^*()+)')

INPUTS_LONG_PREFIX = ('The following input files have a prefix longer than '
					  '30 characters:\n{0}\nPlease make sure that input '
					  'files have a shorter and unique prefix.')

MISSING_INPUT_ARG = ('Path to input files does not exist. Please provide a valid path.')
