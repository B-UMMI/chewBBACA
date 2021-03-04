#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This modules contains constants used to define default parameters
or variables values.

Code documentation
------------------
"""


import platform


# BLAST Score Ratio defaults
BSR_MIN = 0.0
BSR_MAX = 1.0
DEFAULT_BSR = 0.6

# minimum sequence length defaults
MSL_MIN = 0
# large value to ensure that all sequences above minimum value are accepted
MSL_MAX = 99999
# 201 nucleotides (67 aminoacids)
MINIMUM_LENGTH_DEFAULT = 201

# size variation threshold defaults
ST_MIN = 0.0
ST_MAX = 1.0
# new alleles are inferred if their length value does
# not deviate more than this value from the locus sequence
# length mode
SIZE_THRESHOLD_DEFAULT = 0.2

# word size used in clustering
WORD_SIZE_MIN = 5
WORD_SIZE_MAX = 5
WORD_SIZE_DEFAULT = 5

# number of kmers in a window
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
# clustered sequences are excluded if they share this proportion
# of distinct minimizers with the cluster representative
REPRESENTATIVE_FILTER_DEFAULT = 0.9

# decimal proportion of shared distinct minimizers with other
# clustered sequences
INTRA_CLUSTER_MIN = 0.9
INTRA_CLUSTER_MAX = 0.9
# clustered sequences are excluded if they share this proportion
# of distinct minimizers with another clustered sequence of equal
# or greater length
INTRA_CLUSTER_DEFAULT = 0.9

# genetic codes/translation tables
GENETIC_CODES = {1: 'Standard',
                 4: 'The mold, protozoan, and coelenterate mitochondrial code and the mycoplasma/spiroplasma code',
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

# list with chewie-NS instances identifiers and URLs
HOST_NS = {'main': 'https://chewbbaca.online/api/NS/api/',
           'tutorial': 'https://tutorial.chewbbaca.online/api/NS/api/',
           'local': 'http://127.0.0.1:5000/NS/api/'}

# UniProt SPARQL endpoint
UNIPROT_SPARQL = 'http://sparql.uniprot.org/sparql'
MAX_QUERIES = 10

# authors names, GitHub repository, documentation,
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
IGNORE_RAISED = ['Warning: [blastp] To obtain better run time performance, please run '
                 'blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> -seqid_file_out '
                 '<OUT_FILE_NAME> and use <OUT_FILE_NAME> as the argument to -seqidlist']

# replacements for genome and loci identifiers
CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""),
                     (")", ""), ("'", ""), ("\"", ""), (":", "")]

# minimum Python version
MIN_PYTHON = [(3, 4, 0), '3.4.0']
