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

CHEWIE_VERSIONS = ['2.5.0']

# BSR
BSR_MIN = 0.0
BSR_MAX = 1.0

# minimum sequence length
MSL_MIN = 0
MSL_MAX = 99999

# size variation threshold
ST_MIN = 0.0
ST_MAX = 1.0

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

FASTA_SUFFIXES = ['.fasta', '.fna', '.ffn']

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

HOST_NS = 'https://chewbbaca.online/api/NS/api/'
#HOST_NS = 'http://127.0.0.1:5000/NS/api/'

# UniProt SPARQL endpoint
UNIPROT_SPARQL = 'http://sparql.uniprot.org/sparql'
MAX_QUERIES = 10
