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


import platform


# BLAST Score Ratio default values
# value must be contained in the [0.0,1.0] interval
BSR_MIN = 0.0
BSR_MAX = 1.0
DEFAULT_BSR = 0.6

# minimum sequence length defaults
MSL_MIN = 0
# large value to ensure that all sequences above minimum value are accepted
MSL_MAX = 99999
# default value of 201 nucleotides (67 aminoacids)
MINIMUM_LENGTH_DEFAULT = 201

# size variation threshold defaults
ST_MIN = 0.0
ST_MAX = 1.0
# new alleles are inferred if their length value does
# not deviate more than this value from the locus sequence
# length mode
# ASM if below threshold and ALM if above
SIZE_THRESHOLD_DEFAULT = 0.2

# word size/k size used for minimizer clustering
# this value should not be modified
WORD_SIZE_MIN = 5
WORD_SIZE_MAX = 5
WORD_SIZE_DEFAULT = 5

# number of adjacent kmers to consider when selecting minimizers
# this value should not be modified
WINDOW_SIZE_MIN = 5
WINDOW_SIZE_MAX = 5
WINDOW_SIZE_DEFAULT = 5

# minimum decimal proportion of shared distinct minimizers for
# a sequence to be added to a cluster
CLUSTERING_SIMILARITY_MIN = 0.20
CLUSTERING_SIMILARITY_MAX = 0.20
CLUSTERING_SIMILARITY_DEFAULT = 0.20

# decimal proportion of shared distinct minimizers with cluster
# representative
REPRESENTATIVE_FILTER_MIN = 0.9
REPRESENTATIVE_FILTER_MAX = 0.9
# in the CreateSchema process, clustered sequences are excluded
# if they share this proportion of distinct minimizers with the
# cluster representative
REPRESENTATIVE_FILTER_DEFAULT = 0.9

# decimal proportion of shared distinct minimizers with other
# clustered sequences
INTRA_CLUSTER_MIN = 0.9
INTRA_CLUSTER_MAX = 0.9
# in the CreateSchema process, clustered sequences are excluded
# if they share this proportion of distinct minimizers with another
# clustered sequence of equal or greater length
INTRA_CLUSTER_DEFAULT = 0.9

# genetic codes/translation tables
GENETIC_CODES = {1: 'Standard',
                 4: 'The mold, protozoan, and coelenterate mitochondrial '
                    'code and the mycoplasma/spiroplasma code',
                 11: 'The Bacterial, Archaeal and Plant Plastid code',
                 25: 'Candidate division SR1 and gracilibacteria code'}

GENETIC_CODES_DEFAULT = 11

# valid FASTA file extensions
FASTA_SUFFIXES = ['.fasta', '.fna', '.ffn', '.fa']

# NS related constants
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

# list of Chewie-NS instance identifiers and URLs
HOST_NS = {'main': 'https://chewbbaca.online/api/NS/api/',
           'tutorial': 'https://tutorial.chewbbaca.online/api/NS/api/',
           'local': 'http://127.0.0.1:5000/NS/api/'}

# UniProt SPARQL endpoint
UNIPROT_SPARQL = 'http://sparql.uniprot.org/sparql'
MAX_QUERIES = 10

# authors, GitHub repository, documentation,
# tutorial repository and contacts
authors = 'Mickael Silva, Pedro Cerqueira, Rafael Mamede'
repository = 'https://github.com/B-UMMI/chewBBACA'
wiki = 'https://github.com/B-UMMI/chewBBACA/wiki'
tutorial = 'https://github.com/B-UMMI/chewBBACA_tutorial'
contacts = 'imm-bioinfo@medicina.ulisboa.pt'

# timeout, in seconds, for user input
prompt_timeout = 30

# minimum MAJOR and MINOR BLAST versions
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

# replacements for genome and loci identifiers
CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""),
                     (")", ""), ("'", ""), ("\"", ""),
                     (":", "")]

# minimum Python version
MIN_PYTHON = [(3, 6, 0), '3.6.0']

# FTP to get UniProt's reference proteomes
UNIPROT_PROTEOMES_FTP = ('ftp://ftp.uniprot.org/pub/databases/'
                         'uniprot/current_release/knowledgebase/'
                         'reference_proteomes/')

# list of UniProt's uninformative terms
UNIPROT_UNINFORMATIVE = ['uncharacterized', 'hypothetical', 'duf']

# AlleleCall logfile basename and content template
LOGFILE_BASENAME = 'logging_info.txt'
LOGFILE_TEMPLATE = ('Started script at: {0}\n'
                    'Finished script at: {1}\n'
                    'Number of inputs: {2}\n'
                    'Number of loci: {3}\n'
                    'Used this number of CPU cores: {4}\n'
                    'Used a BSR of: {5}\n')

# basename for files created by AlleleCall
RESULTS_ALLELES_BASENAME = 'results_alleles.tsv'
UNCLASSIFIED_BASENAME = 'unclassified_sequences.fasta'
PARALOGS_BASENAME = 'RepeatedLoci.tsv'
LOCI_STATS_BASENAME = 'loci_summary_stats.tsv'
CDS_COORDINATES_BASENAME = 'cds_coordinates.tsv'
PARALOGOUS_COUNTS_BASENAME = 'paralogous_counts.tsv'
PARALOGOUS_LOCI_BASENAME = 'paralogous_loci.tsv'
# header for TSV file with loci stats
LOCI_STATS_HEADER = ('Locus\tEXC\tINF\tPLOT3\tPLOT5\tLOTSC\tNIPH\t'
                     'NIPHEM\tALM\tASM\tLNF\tTotal_CDS')
# header for TSV file with information about extracted CDSs
CDS_TABLE_HEADER = 'Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n'
# header for TSV file with paralogous counts
PARALOGOUS_COUNTS_HEADER = 'Locus\tParalog_Count'
# header for TSV file with paralogous loci per genome
PARALOGOUS_LIST_HEADER = 'Genome\tLoci\tCDS'

# allele calling classifications
ALLELECALL_CLASSIFICATIONS = ['EXC', 'INF', 'PLOT3', 'PLOT5',
                              'LOTSC', 'NIPH', 'NIPHEM', 'ALM',
                              'ASM', 'LNF']

# regex pattern to match locus identifier
LOCUS_ID_PATTERN = r'.*-protein[0-9]+'

# string template for a standard single line FASTA record
FASTA_RECORD_TEMPLATE = '>{0}\n{1}'
FASTA_CDS_TEMPLATE = '>{0}-protein{1}\n{2}'

DNA_BASES = 'AGCT'

# define default BLASTp task
BLAST_TASK_THRESHOLD = {'blastn': 50, 'blastp': 30}

# BLAST outfmt
BLAST_DEFAULT_OUTFMT = '6 qseqid qstart qend qlen sseqid slen score'

# dictionary template to map variables returned by AlleleCall
ALLELECALL_DICT = {'classification_files': None,
                   'basename_map': None,
                   'cds_coordinates': None,
                   'dna_fasta': None,
                   'protein_fasta': None,
                   'dna_hashtable': None,
                   'protein_hashtable': None,
                   'invalid_inputs': None,
                   'invalid_alleles': None,
                   'unclassified_ids': None,
                   'self_scores': None,
                   'representatives': None}

GENOME_LIST = 'listGenomes2Call.txt'
LOCI_LIST = 'listGenes2Call.txt'
