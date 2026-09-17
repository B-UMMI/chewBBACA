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


# BLAST Score Ratio default values
# Value must be contained in the [0.0,1.0] interval
BSR_MIN = 0.0
BSR_MAX = 1.0
BSR_DEFAULT = 0.6

# Minimum sequence length defaults
MSL_MIN = 0
# Large value to ensure that all sequences above minimum value are accepted
MSL_MAX = sys.maxsize
# Default minimum sequence length value of 201 nucleotides (67 aminoacids)
MSL_DEFAULT = 201

# Size variation threshold defaults
ST_MIN = 0.0
ST_MAX = 1.0
# New alleles are inferred if their length value does
# not deviate more than this value from the locus sequence
# length mode
# ASM if below threshold and ALM if above
ST_DEFAULT = 0.2

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
GENETIC_CODES = {1: 'The Standard Code',
				 2: 'The Vertebrate Mitochondrial Code',
				 3: 'The Yeast Mitochondrial Code',
				 4: 'The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code',
				 5: 'The Invertebrate Mitochondrial Code',
				 6: 'The Ciliate, Dasycladacean and Hexamita Nuclear Code',
				 9: 'The Echinoderm and Flatworm Mitochondrial Code',
				 10: 'The Euplotid Nuclear Code',
				 11: 'The Bacterial, Archaeal and Plant Plastid Code',
				 12: 'The Alternative Yeast Nuclear Code',
				 13: 'The Ascidian Mitochondrial Code',
				 14: 'The Alternative Flatworm Mitochondrial Code',
				 15: 'Blepharisma Nuclear Code',
				 16: 'Chlorophycean Mitochondrial Code',
				 21: 'Trematode Mitochondrial Code',
				 22: 'Scenedesmus obliquus Mitochondrial Code',
				 23: 'Thraustochytrium Mitochondrial Code',
				 24: 'Rhabdopleuridae Mitochondrial Code',
				 25: 'Candidate Division SR1 and Gracilibacteria Code'}

GENETIC_CODE_DEFAULT = 11

# Proteins to cluster are divided into a maximum
# of 40 smaller groups in CreateSchema
# Dividing based on the number of CPU cores can lead to
# variable results because we do not have pre-defined clusters
# in CreateSchema.
CREATESCHEMA_CLUSTERING_NGROUPS = 40

# Valid FASTA file extensions
FASTA_EXTENSIONS = ['.fasta', '.fna', '.ffn', '.fa', '.fas']

# Chewie-NS related constants
HEADERS_GET = {'Authorization': None,
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
AUTHORS = 'Rafael Mamede, Pedro Cerqueira, Mickael Silva, João Carriço, Mário Ramirez'
REPOSITORY = 'https://github.com/B-UMMI/chewBBACA'
DOCUMENTATION = 'https://chewbbaca.readthedocs.io/en/latest/index.html'
CONTACTS = 'imm-bioinfo@medicina.ulisboa.pt'

# Timeout, in seconds, to wait for user input
PROMPT_TIMEOUT = 30

# Minimum MAJOR and MINOR BLAST versions
BLAST_MAJOR = 2
BLAST_MINOR = 9

# Paths to BLASTp and makeblastdb executables in Linux and Windows
BLASTP_ALIAS = 'blastp'
MAKEBLASTDB_ALIAS = shutil.which('makeblastdb')
BLASTDB_ALIASTOOL_ALIAS = shutil.which('blastdb_aliastool')
BLASTDBCMD_ALIAS = shutil.which('blastdbcmd')

# Protein to create dummy FASTA records used to check if sequence IDs are interpreted as PDB IDs
DUMMY_PROT = 'MKFFYRPTGLAISINDAYQKVNFSTDGSSLRVDNPTPYFITYDQIKINGKSVKNVDMVAPYSQQTYPFKGARANETVQWTVVNDYGGDQKGESILH'
DUMMY_FASTA = 'dummy.fasta'
DUMMY_BLASTDB = 'dummy_db'
DUMMY_DIR = 'dummy_dir'
DUMMY_BLASTDBCMD_FASTA = 'dummy_blastdbcmd.fasta'

# BLAST warnings to be ignored
# This warning is raised in BLAST>=2.10 when passing a TXT file with sequence identifiers to -seqidlist
# To avoid this warning the TXT file must be converted to binary with the blastdb_aliastool
# Performance can be severely affected if the TXT is not converted to binary
# Since the TXT file is converted with blastdb_aliastool, it is no longer necessary to ignore this warning
IGNORE_RAISED = ['Warning: [blastp] To obtain better run time '
				 'performance, please run blastdb_aliastool '
				 '-seqid_file_in <INPUT_FILE_NAME> -seqid_file_out '
				 '<OUT_FILE_NAME> and use <OUT_FILE_NAME> as the '
				 'argument to -seqidlist']

# Path to MAFFT executable
MAFFT_ALIAS = shutil.which('mafft')
# MAFFT default parameters
MAFFT_DEFAULT_PARAMETERS = ['--thread', '1', '--retree', '1', '--maxiterate', '0', '--treeout']

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
RESULTS_ALLELES_MASKED_BASENAME = 'results_alleles_masked.tsv'
RESULTS_STATISTICS_BASENAME = 'results_statistics.tsv'
LOCI_STATS_BASENAME = 'loci_summary_stats.tsv'
UNCLASSIFIED_BASENAME = 'unclassified_sequences.fasta'
MISSING_FASTA_BASENAME = 'missing_classes.fasta'
MISSING_TSV_BASENAME = 'missing_classes.tsv'
NOVEL_BASENAME = 'novel_alleles.fasta'
GENE_COORDINATES_BASENAME = 'gene_coordinates.tsv'
GENE_COORDINATES_EXCLUDED_BASENAME = 'gene_coordinates_excluded.tsv'
INVALID_CDS_BASENAME = 'invalid_cds.txt'
SCHEMA_CONFIG_BASENAME = '.schema_config'
NS_CONFIG_BASENAME = ".ns_config"
GENE_LIST_BASENAME = '.genes_list'

# Header for TSV file with loci stats
LOCI_STATS_HEADER = ('Locus\tEXC\tINF\tPLOT3\tPLOT5\tLOTSC\tNIPH\t'
					 'NIPHEM\tALM\tASM\tPAMA\tLNF\tTotal CDSs Classified')
# COLUMNS for the DataFrame with assembly statistics created by the PredictGenes module
ASSEMBLY_STATS_COLUMNS = ['FILE', 'Number of contigs', 'Average contig size',
						  'N50', 'Assembly size', 'GC content', 'Missing values (Ns)',
						  'Total genes', 'Excluded genes']
# Filename and Headers for the file with the list of novel alleles
NOVEL_ALLELES_FILENAME = "new_alleles.tsv"
NOVEL_ALLELES_LIST_HEADER = 'Allele ID\tCDS ID\tRepresentative'
# Header for TSV file with information about extracted CDSs
GENE_TABLE_HEADER = 'CDS_ID\tGenome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\tConfidence\tSHA256\n'
# Headers for TSV files with paralogous loci count and per genome
PARALOGOUS_COUNTS_HEADER = 'Locus\tCount'
PARALOGOUS_LIST_HEADER = 'Genome\tLoci\tCDS'
# Header for TSV file with information about CDSs classified as ambiguous
MISSING_HEADER = 'Index\tGenome\tLocus\tLocus_classification\tCDS\tCDS_classification'
# Header for TSV file created by the GetAlleles module
GETALLELES_LOCI_STATS_HEADER = 'Locus\tTotal alleles in schema\tSamples with locus\tDistinct alleles in dataset'
# Filenames for the files created by the ComputeMSA module
COMPUTEMSA_PROTEIN_CONCAT = 'protein_concat.fasta'
COMPUTEMSA_DNA_CONCAT = 'dna_concat.fasta'
COMPUTEMSA_PROTEIN_MSA = 'protein_msa.fasta'
COMPUTEMSA_PROTEIN_MSA_VARIABLE = 'protein_msa_variable.fasta'
COMPUTEMSA_DNA_MSA = 'dna_msa.fasta'
COMPUTEMSA_DNA_MSA_VARIABLE = 'dna_msa_variable.fasta'
# Ambiguous characters for DNA and protein sequences
PROTEIN_AMBIGUOUS_CHARS = ['B', 'Z', 'X', 'J']
DNA_AMBIGUOUS_CHARS = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
# Gap character used in MSAs
GAP_CHAR = '-'

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

GENOME_LIST = 'input_files.txt'
LOCI_LIST = 'loci_list.txt'

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
CGMLST_THRESHOLD_MIN = 0
CGMLST_THRESHOLD_MAX = 1

GENOMES_MISSING_COLUMNS = ['Sample', 'Loci presence count', 'Loci presence proportion']
GENOMES_MISSING_BASENAME = 'sample_presence_stats.tsv'
LOCI_MISSING_COLUMNS = ['Locus', 'Sample presence count', 'Sample presence proportion']
LOCI_MISSING_BASENAME = 'loci_presence_stats.tsv'

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
NON_FASTA_EXCEPTION = ('Some of the input files are not in FASTA format.')

# Input directory does not contain FASTA files
MISSING_FASTAS_EXCEPTION = ('Could not get input files. Please provide '
							'a directory with FASTA files or a file with '
							'the list of full paths to the FASTA files '
							'and ensure that filenames end with one of '
							f'the following extensions: {FASTA_EXTENSIONS}.')

MISSING_SCHEMA_FASTAS = ("Input path does not include FASTA files. Please provide "
						 "a valid path for a folder containing only FASTA files ending "
						 f"in one of the following file extensions: {FASTA_EXTENSIONS}")

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
ADAPT_LEGACY_SCHEMA = ('Schema does not include a config file or includes files that'
					   ' are no longer used or supported. Please use the PrepExternalSchema '
					   'module to adapt the schema and make it compatible with chewBBACA\'s '
					   'latest version.')

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

# Could not predict CDSs for input FASTA files
# e.g. files only contain sequence headers, contain invalid
# sequences/chars or pyrodigal cannot predict any genes
CANNOT_PREDICT = ('Could not predict CDSs from any of the input files.'
				  '\nPlease verify the format of the input files.')

INVALID_BSR = ('\nBSR value is not contained in the [0.0, 1.0] interval.')
INVALID_BSR_TYPE = ('\nInvalid BSR value of {0}. BSR value must be contained in the [0.0, 1.0] interval.')

INVALID_MINLEN = ('\nInvalid minimum sequence length value. Must be equal or greater than 0.')
INVALID_MINLEN_TYPE = ('\nInvalid minimum sequence length value. Value must be a positive integer.')

INVALID_ST = ('\nInvalid size threshold value. Must be contained in the [0.0, 1.0] interval.')
INVALID_ST_TYPE = ('\nInvalid size threshold value used to create schema. Value must be None or a positive float in the [0.0, 1.0] interval.')

INVALID_GENETIC_CODE = ('\nInvalid genetic code value.\nValue must correspond to '
				 		'one of the accepted genetic codes\n\nAccepted genetic '
				 		f'codes:\n{"\n".join([str(k)+": "+v for k, v in GENETIC_CODES.items()])}')

INVALID_WORD_SIZE = ("Invalid clustering word size value.")

INVALID_WINDOW_SIZE = ("Invalid clustering window size value.")

INVALID_CLUSTERING_SIMILARITY = ("Invalid clustering similarity threshold value.")

INVALID_REPRESENTATIVE_FILTER = ("Invalid clustering representative filter threshold value.")

INVALID_INTRA_CLUSTER_FILTER = ("Invalid clustering intra-cluster filter threshold value.")

NS_CANNOT_CONNECT = ('Failed to establish a connection to the Chewie-NS instance at {0}.')

PYTHON_VERSION = ('Python version found: {0}\nPlease use Python >= {1}')

CPU_VALUE_WARNING = ('You have provided a CPU core count value that is equal to the '
					 'number of CPU cores in your system! This may affect your system '
					 'responsiveness.')

BLAST_MISSING = ("Could not find the BLAST executables. Please make sure that BLAST "
				 "is installed and added to the PATH environment variable.")
BLAST_INVALID_VERSION = ("Could not determine BLAST version or the version is not valid. Please make "
				 		 f"sure that BLAST>={BLAST_MAJOR}.{BLAST_MINOR} is installed.")

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

INPUTS_SHARE_PREFIX = ('The following input files share the same filename prefix '
					   '(substring before the first "." in the filename):\n{0}\n'
					   'Please ensure that every input file has a unique '
					   'filename prefix. This is necessary to unambiguously '
					   'identify each input during intermediate steps and in '
					   'output files.')

INPUTS_INCLUDE_BLANKS = ('The following input files include blank spaces '
						 'in the filename:\n{0}\nPlease ensure that filenames '
						 'do not include blank spaces or special characters '
						 '(e.g. !@#?$^*()+). This is necessarry to avoid issues '
						 'related to how these characters are recognized by '
						 'chewBBACA and its dependencies.')

INPUTS_LONG_PREFIX = ('The following input files have a prefix longer than '
					  '30 characters:\n{0}\nPlease make sure that input '
					  'files have a shorter and unique prefix (substring before '
					  'the first "." in the filename). The prefixes are used as '
					  'unique identifiers and long prefixes might lead to issues '
					  '(e.g. BLAST does not accept sequence IDs longer than 50 '
					  'characters when creating a database).')

INPUTS_PDB_PREFIX = ('The following input files have prefixes that are interpreted by BLAST '
					 'as chain PDB IDs:\n{0}\nBLAST modifies the '
					 'IDs of the CDSs that include these prefixes when creating a database, '
					 'which leads to issues when chewBBACA cannot find the original '
					 'IDs in the results. Please ensure that the file prefixes (substring '
					 'before the first "." in the filename) cannot be interpreted as chain PDB IDs.')

MISSING_INPUT_ARG = ('Path to input files does not exist. Please provide a valid path.')

MISSING_ALLELES = ('\nCould not create the FASTA files for {0} loci.'
				   'Some alleles are not in the schema\'s FASTA files. Alleles are not '
				   'added to the schema if the allele calling process did not '
				   'complete successfully or if the --no-inferred option is used.')

MISSING_LOCI_LIST = ("Path for the loci list is invalid.")

ALM_MSG = ('allele greater than {0}% locus length mode ({1}>{2})')
ASM_MSG = ('allele smaller than {0}% locus length mode ({1}<{2})')

NO_MSAS_CREATED = ('Could not compute the MSA for any of the files. Exiting...')

COMPUTEMSA_NO_SCHEMA = ('Schema directory must be provided when input is a TSV file with allelic profiles.')

# Define sequence identifier prefix to avoid issues where makeblastdb modifies the IDs
# of the sequences when creating the database (e.g. when IDs are interpreted as PDB chain IDs)
BLASTDB_LCL_PREFIX = 'lcl|SEQ'
BLASTDB_SEQ_PREFIX = 'SEQ'

# List of valid files for the MergeResults module
# Includes all standard files created by the AlleleCall module plus the ones created with specific parameters
ALLELECALL_OUTFILES = ['results_alleles.tsv', 'results_contigsInfo.tsv',
					   'results_statistics.tsv', 'loci_summary_stats.tsv',
					   'invalid_cds.txt', 'cds_coordinates.tsv',
					   'missing_classes.fasta', 'missing_classes.tsv',
					   'unclassified_sequences.fasta', 'paralogous_counts.tsv',
					   'paralogous_loci.tsv', 'novel_alleles.fasta',
					   'presence_absence.tsv']

SUBSETRESULTS_MISSING_LISTS_EXCEPTION = ('Did not provide a list of loci or samples. Please '
										 'provide at least one of those lists to subset results.')

SUBSETRESULTS_MISSING_FILES_EXCEPTION = ('Did not find any of the valid file types to subset. Please make '
		   		 						 'sure that the input directory contains at least the file type '
				 						 'created by chewBBACA to store the allelic profiles (the basename '
				 						 'must be `results_alleles.tsv`).')

SUBSETRESULTS_MISSING_PROFILES_EXCEPTION = ('Did not find the file type used by chewBBACA to store the allelic '
		   		 						    'profiles (`results_alleles.tsv`). Please make sure that the input '
				 						    'directory contains at least that file type.')

SUBSETRESULTS_ABSENT_LOCI_EXCEPTION = ('Please make sure that all loci IDs in the list provided to '
		 		  	 				   '`--loci-list` match loci IDs in the input results.')

SUBSETRESULTS_ABSENT_SAMPLES_EXCEPTION = ('Please make sure that all sample IDs in the list provided to '
				  	 					  '`--sample-list` match sample IDs in the input results.')

MERGERESULTS_INPUTFILE_EXCEPTION = ('Input paths should all be folders. The following are not folders:\n{0}')

MERGERESULTS_NOCOMMONLOCI_EXCEPTION = ('Results do not have loci in common.')

MERGERESULTS_LOCIDIFFER_EXCEPTION = ('Files have different sets of loci. Please provide '
				 					 'files with the results for the same set of loci or '
				 					 'provide the "--common" parameter to create merged files '
				 					 'for the set of loci shared by all results folders.')

HASHPROFILES_INVALID_HASHING = ('{0} hash function is not available in the hashlib '
				 				'(https://docs.python.org/3/library/hashlib.html) and '
				 				'zlib (https://docs.python.org/3/library/zlib.html) modules.')

GENE_PREDICTORS = ["pyrodigal", "augustus"]

INVALID_GENE_PREDICTOR = ("Specified gene predictor is not valid.")

AUGUSTUS_ALIAS = 'augustus'

AUGUSTUS_MISSING = ('Could not find AUGUSTUS executables. Please ensure that AUGUSTUS is installed and the executables were added to PATH.')

AUGUSTUS_INVALID_SPECIES = ("Specified species ID is not in the list of supported species for AUGUSTUS. Please provide a valid species ID.")

AUGUSTUS_OUTFMTS = ['genes', 'gff']

AUGUSTUS_OUTFMT_DEFAULT = "genes"

AUGUSTUS_INVALID_OUTFMT = ("Invalid output format specified for AUGUSTUS results.")

PYRODIGAL_MODES = ["single", "meta"]

PYRODIGAL_DEFAULT_MODE = "single"

PYRODIGAL_OUTFMTS = ['genes', 'translations', 'gff', 'genbank', 'scores']

PYRODIGAL_DEFAULT_OUTFMT = ['genes']

PYRODIGAL_MIN_CONFIDENCE = 0.0

PYRODIGAL_MAX_CONFIDENCE = 100.0

GENE_PREDICTOR_DEFAULT = "pyrodigal"

GENE_PREDICTION_DEFAULT_ARGUMENTS = {"pyrodigal_training_file": None,
									 "pyrodigal_mode": PYRODIGAL_DEFAULT_MODE,
									 "pyrodigal_output_formats": PYRODIGAL_DEFAULT_OUTFMT,
									 "pyrodigal_minimum_confidence": None,
									 "pyrodigal_training_reference": None,
									 "pyrodigal_just_training": False,
									 "augustus_species": None,
							  		 "augustus_output_formats": ["genes"],
							  		 "augustus_path": None}

CLUSTERING_DEFAULT_ARGUMENTS = {"word_size": WORD_SIZE_DEFAULT,
								"window_size": WINDOW_SIZE_DEFAULT,
								"clustering_sim": CLUSTERING_SIMILARITY_DEFAULT,
								"representative_filter": REPRESENTATIVE_FILTER_DEFAULT,
								"intra_filter": INTRA_CLUSTER_DEFAULT}

PYRODIGAL_INVALID_OUTFMT = ("Invalid output format specified for Pyrodigal results.")

INVALID_PYRODIGAL_MODE = ("Invalid mode specified for Pyrodigal. Please provide a valid mode.")

INVALID_PYRODIGAL_CONFIDENCE = ("Invalid confidence value specified for Pyrodigal. Please provide a value between 0.0 and 100.0.")

AUGUSTUS_GFF_FILTERS = ["# start gene", "# coding sequence"]

VALID_PARAMETERS = {"augustus": ["augustus_species", "augustus_output_formats", "augustus_path"],
					"pyrodigal": ["pyrodigal_training_file", "pyrodigal_mode",
							  	  "pyrodigal_output_formats", "pyrodigal_minimum_confidence",
							  	  "pyrodigal_training_reference", "pyrodigal_just_training"]}

MISSING_INPUTS = ("Some of the paths to input files are not valid.")

CPU_CORES_DEFAULT = 1
CPU_CORES_MIN = 1

# List of parameter names used by chewBBACA
INPUT_FILES_ARGNAME = "input_files"
OUTPUT_DIRECTORY_ARGNAME = "output_directory"
GENE_PREDICTOR_ARGNAME = "gene_predictor"
GENE_PREDICTION_OPTIONS_ARGNAME = "gene_prediction_options"
TRANSLATION_TABLE_ARGNAME = "translation_table"
CPU_CORES_ARGNAME = "cpu_cores"
SCHEMA_NAME_ARGNAME = "schema_name"
BLAST_SCORE_RATIO_ARGNAME = "blast_score_ratio"
MINIMUM_LENGTH_ARGNAME = "minimum_length"
SIZE_THRESHOLD_ARGNAME = "size_threshold"
CLUSTERING_OPTIONS_ARGNAME = "clustering_options"
BLAST_PATH_ARGNAME = "blast_path"
CDS_INPUT_ARGNAME = "cds_input"
NO_CDS_RENAMING_ARGNAME = "no_cds_renaming"
NO_CLEANUP_ARGNAME = "no_cleanup"
SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
LOCI_LIST_ARGNAME = "loci_list"
NO_INFERRED_ARGNAME = "no_inferred"
OUTPUT_UNCLASSIFIED_ARGNAME = "output_unclassified"
OUTPUT_MISSING_ARGNAME = "output_missing"
OUTPUT_NOVEL_ARGNAME = "output_novel"
OUTPUT_MASKED_ARGNAME = "output_masked"
FORCE_CONTINUE_ARGNAME = "force_continue"
EXECUTION_MODE_ARGNAME = "execution_mode"
ANNOTATIONS_FILE_ARGNAME = "annotations_file"
LOCI_REPORTS_ARGNAME = "loci_reports"
LIGTH_ARGNAME = "light"
ADD_SEQUENCES_ARGNAME = "add_sequences"
RESULTS_FILES_ARGNAME = "results_files"
NO_PRESENCE_ABSENCE_ARGNAME = "no_presence_absence"
NO_DISTANCE_MATRIX_ARGNAME = "no_distance_matrix"
NO_NEIGHBOR_JOINING_ARGNAME = "no_neighbor_joining"
FORCE_CORE_MSA_ARGNAME = "force_core_msa"
LOCI_PRESENCE_THRESHOLD_ARGNAME = "loci_presence_threshold"
SAMPLE_STEP_ARGNAME = "sample_step"
COMPUTE_ACCESSORY_ARGNAME = "compute_accessory"
RAREFACTION_ANALYSIS_ARGNAME = "rarefaction_analysis"
PERMUTATION_NUMBER_ARGNAME = "permutation_number"
PERMUTATION_SAMPLES_ARGNAME = "permutation_samples"
EXCLUDE_LOCI_ARGNAME = "exclude_loci"
EXCLUDE_GENOMES_ARGNAME = "exclude_genomes"
SAMPLE_LIST_ARGNAME = "sample_list"
INVERT_LOCI_ARGNAME = "invert_loci"
INVERT_SAMPLES_ARGNAME = "invert_samples"
COMMON_ARGNAME = "common"
ALLELIC_PROFILES_ARGNAME = "allelic_profiles"
HASHING_ALGO_ARGNAME = "hashing_algo"
NROW_CHUNK_ARGNAME = "nrow_chunk"
TRANSLATE_ALLELES_ARGNAME = "translate_alleles"
DISTINCT_ARGNAME = "distinct"
SIZE_FILTER_ARGNAME = "size_filter"
PROTEIN_TABLE_ARGNAME = "protein_table"
PROTEOME_TAXA_ARGNAME = "proteome_taxa"
PROTEOME_MATCHES_ARGNAME = "proteome_matches"
NO_SPARQL_ARGNAME = "no_sparql"
DISTANCE_COMPUTATION_METHOD_ARGNAME = "distance_computation_method"
OUTPUT_FORMAT_ARGNAME = "output_format"
COMPUTE_SIMILARITY_ARGNAME = "compute_similarity"
NO_MASK_ARGNAME = "no_mask"
DNA_MSA_ARGNAME = "dna_msa"
OUTPUT_VARIABLE_ARGNAME = 'output_variable'
GAPS_ARGNAME = "gaps"
AMBIGUOUS_ARGNAME = "ambiguous"
ONLY_LOCI_MSAS_ARGNAME = "only_loci_msas"
CUSTOM_MAFFT_OPTIONS_ARGNAME = "custom_mafft_options"
PROTEIN_INPUT_ARGNAME = "protein_input"
SPECIES_ID_ARGNAME = "species_id"
SCHEMA_ID_ARGNAME = "schema_id"
DOWNLOAD_FOLDER_ARGNAME = "download_folder"
NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
SCHEMA_DATE_ARGNAME = "schema_date"
LATEST_VERSION_ARGNAME = "latest_version"
LOCI_PREFIX_ARGNAME = "loci_prefix"
DESCRIPTION_FILE_ARGNAME = "description_file"
CONTINUE_UPLOAD_ARGNAME = "continue_upload"
SUBMIT_ALLELES_ARGNAME = "submit_alleles"
STATS_MODE_ARGNAME = "stats_mode"

# List of argument default values

PYRODIGAL_MODE_ARGNAME = "pyrodigal_mode"
PYRODIGAL_TRAININGFILE_ARGNAME = "pyrodigal_training_file"
AUGUSTUS_PATH_ARGNAME = "augustus_path"
CLUSTERING_WORD_ARGNAME = "word_size"
CLUSTERING_WINDOW_ARGNAME = "window_size"
CLUSTERING_SIMILARITY_ARGNAME = "clustering_sim"
CLUSTERING_REPRESENTATIVEFILTER_ARGNAME = "representative_filter"
CLUSTERING_INTRAFILTER_ARGNAME = "intra_filter"

NSSTATS_MODE_CHOICES = ['species', 'schemas']

DEFAULT_NOMENCLATURE_SERVER = "main"

DEFAULT_GAPS = "exclude"
GAPS_CHOICES = ['ignore', 'exclude']

DEFAULT_AMBIGUOUS = "exclude"
AMBIGUOUS_CHOICES = ['ignore', 'exclude']

DEFAULT_DISTANCE_METHOD = "hamming"
DISTANCE_METHODS = ['hamming', 'jaccard', 'loci', 'core']
INVALID_DISTANCE_METHOD = ("Invalid distance method specified.")

DEFAULT_OUTPUT_FORMAT = "upper_triangular"
OUTPUT_FORMATS = ['upper_triangular', 'lower_triangular', 'symmetric', 'table']
INVALID_OUTPUT_FORMAT = ("Invalid output format specified.")

PERMUTATION_NUMBER_DEFAULT = 100
PERMUTATION_SAMPLES_DEFAULT = None

INVALID_GENEPREDICTOR_PARAMETER = ("{0} is not a valid parameter=argument pair used to configure gene prediction with {1}.")
INVALID_PARAMETER_STR = ("{0} is not a valid parameter=argument pair.")

PYRODIGAL_META_NOPTF = ("Cannot use a training file when running Pyrodigal in meta mode. "
						"Please do not provide a training file if setting Pyrodigal's running mode to meta.")

# Define expected types for the arguments
ARGUMENT_TYPES = {
	"input_files": str,
	"output_directory": str,
	"gene_predictor": str,
	"gene_prediction_options": str,
	"translation_table": int,
	"cpu_cores": int,
	"schema_name": str,
	"blast_score_ratio": float,

	"minimum_length": int,
	"size_threshold": float,
	"clustering_parameters": str,
	"blast_path": str,
	"cds_input": bool,
	"no_cds_renaming": bool,
	"no_cleanup": bool,
	"augustus_species": str,
	"augustus_output_formats": str,
	"augustus_path": str,
	"pyrodigal_training_file": str,
	"pyrodigal_mode": str,
	"pyrodigal_output_formats": str,
	"pyrodigal_minimum_confidence": float,
	"pyrodigal_training_reference": str,
	"pyrodigal_just_training": bool,
	"word_size": int,
	"window_size": int,
	"clustering_sim": float,
	"representative_filter": float,
	"intra_filter": float,
	"annotations": str,
	"threshold": float,
	"step":int,
	"exclude_loci": str,
	"exclude_genomes": str,
	"allelic_profiles": str,
	"hash_type": str,
	"size_filter": bool,
	"protein_table": str,
	"taxa":str,
	"proteome_matches": int,
	"no_sparql": bool,
	"output_format": str,
	'output_variable': bool,
}

CANNOT_PROVIDE_PTF_AND_TREF = ("Cannot provide a training file and a training reference "
							   "to create a training file. Please provide either a training "
							   "file or a training reference to create a new training file, "
							   "not both.")

SCHEMA_NAME_DEFAULT = "schema_seed"

# Links to Documentation pages
PredictGenesDocs = "https://chewbbaca.readthedocs.io/en/latest/user/modules/PredictGenes.html"
CreateSchemaDocs = "https://chewbbaca.readthedocs.io/en/latest/user/modules/CreateSchema.html"
AlleleCallDocs = "https://chewbbaca.readthedocs.io/en/latest/user/modules/AlleleCall.html"

ALLELECALL_MODES = [1, 2, 3, 4]
ALLELECALL_DEFAULT_MODE = 1

INVALID_ALLELECALL_MODE = ("Specified invalid mode for allele calling.")

MISSING_RESULTS_ALLELES = (f"Input directory with results data does not include the {RESULTS_ALLELES_BASENAME} file.")
