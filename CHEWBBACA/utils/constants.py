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

# UniProt SPARQL endpoint
UNIPROT_SPARQL = 'http://sparql.uniprot.org/sparql'
MAX_QUERIES = 10

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

# BLASTp and makeblastdb aliases for Linux and Windows
BLASTP_ALIAS = 'blastp.exe' if platform.system() == 'Windows' else 'blastp'
MAKEBLASTDB_ALIAS = 'makeblastdb.exe' if platform.system() == 'Windows' else 'makeblastdb'

# Prodigal alias (must be in PATH)
PRODIGAL_PATH = 'prodigal'

# BLAST warnings to be ignored
IGNORE_RAISED = ['Warning: [blastp] To obtain better run time '
                 'performance, please run blastdb_aliastool '
                 '-seqid_file_in <INPUT_FILE_NAME> -seqid_file_out '
                 '<OUT_FILE_NAME> and use <OUT_FILE_NAME> as the '
                 'argument to -seqidlist']

# Replacements for genome and loci identifiers
CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""),
                     (")", ""), ("'", ""), ("\"", ""),
                     (":", "")]

# Minimum Python version
MIN_PYTHON = [(3, 6, 0), '3.6.0']

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

# Basename for files created by AlleleCall
RESULTS_ALLELES_BASENAME = 'results_alleles.tsv'
RESULTS_STATISTICS_BASENAME = 'results_statistics.tsv'
UNCLASSIFIED_BASENAME = 'unclassified_sequences.fasta'
NOVEL_BASENAME = 'novel_alleles.fasta'
PARALOGS_BASENAME = 'RepeatedLoci.tsv'
LOCI_STATS_BASENAME = 'loci_summary_stats.tsv'
CDS_COORDINATES_BASENAME = 'cds_coordinates.tsv'
PARALOGOUS_COUNTS_BASENAME = 'paralogous_counts.tsv'
PARALOGOUS_LOCI_BASENAME = 'paralogous_loci.tsv'
# Header for TSV file with loci stats
LOCI_STATS_HEADER = ('Locus\tEXC\tINF\tPLOT3\tPLOT5\tLOTSC\tNIPH\t'
                     'NIPHEM\tALM\tASM\tPAMA\tLNF\tTotal_CDS')
# Header for TSV file with information about extracted CDSs
CDS_TABLE_HEADER = 'Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n'
# Headers for TSV files with paralogous loci count and per genome
PARALOGOUS_COUNTS_HEADER = 'Locus\tCount'
PARALOGOUS_LIST_HEADER = 'Genome\tLoci\tCDS'

# Allele calling classifications
ALLELECALL_CLASSIFICATIONS = ['EXC', 'INF', 'PLOT3', 'PLOT5',
                              'LOTSC', 'NIPH', 'NIPHEM', 'ALM',
                              'ASM', 'PAMA', 'LNF']

# PLNF classificaton for modes {1,2,3}
PROBABLE_LNF = 'PLNF'

# Regex pattern to match locus identifier
LOCUS_ID_PATTERN = r'.*-protein[0-9]+'

# String template for a standard single line FASTA record
FASTA_RECORD_TEMPLATE = '>{0}\n{1}'
FASTA_CDS_TEMPLATE = '>{0}-protein{1}\n{2}'

DNA_BASES = 'AGCT'

# Define default BLASTp task
BLAST_TASK_THRESHOLD = {'blastn': 50, 'blastp': 30}

# BLAST outfmt
BLAST_DEFAULT_OUTFMT = '6 qseqid qstart qend qlen sseqid slen score'

# Dictionary template to map variables returned by AlleleCall
ALLELECALL_DICT = {'classification_files': None,
                   'basename_map': None,
                   'cds_coordinates': None,
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

LOCI_LIST_FILE = '.genes_list'

# Maximum number of allele hashes per pre-computed file
HASH_TABLE_MAXIMUM_ALLELES = 200000

# File header for file with summary statistics created by PrepExternalSchema
PREPEXTERNAL_SUMMARY_STATS_HEADER = ('Gene\tTotal_alleles\tValid_alleles\t'
                                     'Number_representatives')

UNIPROT_SPARQL_THREADS = 4

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

TRANSLATION_EXCEPTIONS = ['Extra in frame stop codon found',
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

# Invalid path to schema
SCHEMA_INVALID_PATH = ('Provided path does not include all the necessary '
                       'schema files. Please verify that you have passed '
                       'the correct path to the schema.')

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
