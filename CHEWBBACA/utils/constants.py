#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import inspect
import platform

# remove this!
CHEWIE_VERSIONS = ['2.5.0', '2.5.1', '2.5.2', '2.5.3',
				   '2.5.4', '2.5.5', '2.5.6', '2.6.0']

# BSR
BSR_MIN = 0.0
BSR_MAX = 1.0
DEFAULT_BSR = 0.6

# minimum sequence length
MSL_MIN = 0
MSL_MAX = 99999
MINIMUM_LENGTH_DEFAULT = 201

# size variation threshold
ST_MIN = 0.0
ST_MAX = 1.0
SIZE_THRESHOLD_DEFAULT = 0.2

# clustering config
WORD_SIZE_MIN = 5
WORD_SIZE_MAX = 5
WORD_SIZE_DEFAULT = 5

WINDOW_SIZE_MIN = 5
WINDOW_SIZE_MAX = 5
WINDOW_SIZE_DEFAULT = 5

CLUSTERING_SIMILARITY_MIN = 0.20
CLUSTERING_SIMILARITY_MAX = 0.20
CLUSTERING_SIMILARITY_DEFAULT = 0.20

REPRESENTATIVE_FILTER_MIN = 0.8
REPRESENTATIVE_FILTER_MAX = 0.8
REPRESENTATIVE_FILTER_DEFAULT = 0.8

INTRA_CLUSTER_MIN = 0.8
INTRA_CLUSTER_MAX = 0.8
INTRA_CLUSTER_DEFAULT = 0.8

# genetic codes/translation tables
GENETIC_CODES = {1: 'Standard',
				 4: 'The mold, protozoan, and coelenterate mitochondrial code and the mycoplasma/spiroplasma code',
				 11: 'The Bacterial, Archaeal and Plant Plastid code',
				 25: 'Candidate division SR1 and gracilibacteria code'}
GENETIC_CODES_DEFAULT = 11

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

HOST_NS = {'main': 'https://chewbbaca.online/api/NS/api/',
           'tutorial': 'https://tutorial.chewbbaca.online/api/NS/api/',
           'local': 'http://127.0.0.1:5000/NS/api/'}

# UniProt SPARQL endpoint
UNIPROT_SPARQL = 'http://sparql.uniprot.org/sparql'
MAX_QUERIES = 10

authors = 'Mickael Silva, Pedro Cerqueira, Rafael Mamede'
repository = 'https://github.com/B-UMMI/chewBBACA'
wiki = 'https://github.com/B-UMMI/chewBBACA/wiki'
tutorial = 'https://github.com/B-UMMI/chewBBACA_tutorial'
contacts = 'imm-bioinfo@medicina.ulisboa.pt'

# timeout when the process asks users for input
prompt_timeout = 30

BLAST_MAJOR = 2
BLAST_MINOR = 9
BLASTP_ALIAS = 'blastp.exe' if platform.system() == 'Windows' else 'blastp'
MAKEBLASTDB_ALIAS = 'makeblastdb.exe' if platform.system() == 'Windows' else 'makeblastdb'

PRODIGAL_PATH = 'prodigal'

IGNORE_RAISED = ['Warning: [blastp] To obtain better run time performance, please run '
				 'blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> -seqid_file_out '
				 '<OUT_FILE_NAME> and use <OUT_FILE_NAME> as the argument to -seqidlist']

SPACED_MATRIX = os.path.join(os.path.dirname(inspect.getfile(inspect.currentframe())), 'spaced_matrix')
