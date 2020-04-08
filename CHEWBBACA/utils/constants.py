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

CHEWIE_VERSIONS = ['2.1.0']

# BSR
BSR_MIN = 0.0
BSR_MAX = 1.0

# minimum sequence length
MSL_MIN = 0
MSL_MAX = 99999

# size variation threshold
ST_MIN = 0.0
ST_MAX = 1.0

# genetic codes/translation tables
GENETIC_CODES = {1: 'Standard',
				 4: 'The mold, protozoan, and coelenterate mitochondrial code and the mycoplasma/spiroplasma code',
				 11: 'The Bacterial, Archaeal and Plant Plastid code',
				 25: 'Candidate division SR1 and gracilibacteria code'}

FASTA_SUFFIXES = ['.fasta', '.fna', '.ffn']

# NS related constants

#HOST_NS = 'https://194.210.120.209/api/NS/api/'
HOST_NS = 'http://127.0.0.1:5000/NS/api/'

SFTP_PTF_TEMP = '/prodigal_training_files/upload_temp'
SFTP_PTF_AVAI = '/prodigal_training_files/available'
