#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the HashProfiles module.
"""


# Parameter names
HASHPROFILES_ALLELIC_PROFILES_ARGNAME = "allelic_profiles"
HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
HASHPROFILES_HASHING_ALGO_ARGNAME = "hashing_algo"
HASHPROFILES_NROW_CHUNK_ARGNAME = "nrow_chunk"
HASHPROFILES_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types
HASHPROFILES_ARGUMENT_TYPES = {
	"allelic_profiles": str,
	"schema_directory": str,
	"output_directory": str,
	"hashing_algo": str,
	"nrow_chunk": int,
	"cpu_cores": int
}

# Parameter descriptions
HASHPROFILES_ALLELIC_PROFILES_DESCRIPTION = ("Path to the TSV file that contains the allelic profiles determined by "
							 "the AlleleCall module.")
HASHPROFILES_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the schema's directory to get the allele sequences and compute the hashes.")
HASHPROFILES_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory.")
HASHPROFILES_HASHING_ALGO_DESCRIPTION = ("Hashing algorithm used to hash the allele sequences. The hashing algorithms "
							 "implemented in the hashlib and zlib Python libraries are supported.")
HASHPROFILES_NROW_CHUNK_DESCRIPTION = ("Divide the input file into chunks of this many rows to process larger files"
						 	 " more efficiently.")
HASHPROFILES_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

# Default hashing algorithm
HASHPROFILES_DEFAULT_HASH_ALGO = "crc32"

# Exceptions
HASHPROFILES_INVALID_HASHING = ('{0} hash function is not available in the hashlib '
				 				'(https://docs.python.org/3/library/hashlib.html) and '
				 				'zlib (https://docs.python.org/3/library/zlib.html) modules.')
