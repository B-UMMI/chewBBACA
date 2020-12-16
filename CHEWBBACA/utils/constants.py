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


import platform


CHEWIE_VERSIONS = ['2.5.0', '2.5.1', '2.5.2', '2.5.3',
				   '2.5.4', '2.5.5', '2.5.6']

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

WS_MIN = 0.0
WS_MAX = 1.0

CS_MIN = 0.0
CS_MAX = 1.0

RF_MIN = 0.0
RF_MAX = 1.0

IF_MIN = 0.0
IF_MAX = 1.0

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
