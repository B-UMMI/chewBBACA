#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


HASHPROFILES_MODULE_NAME = "HashProfiles"

HASHPROFILES_ALLELIC_PROFILES_ARGNAME = "allelic_profiles"
HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
HASHPROFILES_HASHING_ALGO_ARGNAME = "hashing_algo"
HASHPROFILES_NROW_CHUNK_ARGNAME = "nrow_chunk"
HASHPROFILES_CPU_CORES_ARGNAME = "cpu_cores"

HASHPROFILES_ARGUMENT_TYPES = {
	"allelic_profiles": str,
	"schema_directory": str,
	"output_directory": str,
	"hashing_algo": str,
	"nrow_chunk": int,
	"cpu_cores": int
}

HASHPROFILES_DEFAULT_HASH_ALGO = "crc32"

# Exceptions

HASHPROFILES_INVALID_HASHING = ('{0} hash function is not available in the hashlib '
				 				'(https://docs.python.org/3/library/hashlib.html) and '
				 				'zlib (https://docs.python.org/3/library/zlib.html) modules.')
