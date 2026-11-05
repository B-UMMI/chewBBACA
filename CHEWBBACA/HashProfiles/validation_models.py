#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the HashProfiles module.
"""


import sys
import zlib
import hashlib
from typing import Annotated

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

import CHEWBBACA.global_validation as gval
import CHEWBBACA.HashProfiles.constants as ct
from CHEWBBACA.utils import file_operations as fo


# Function used to validate the hashing algorithm specified by the user
def validate_hashing_algo(hashing_algo):
	"""Validate the hashing algorithm.

	Parameters
	----------
	hashing_algo : str
		The name of the hashing algorithm to validate.

	Returns
	-------
	tuple
		A tuple containing the validated hashing algorithm and its corresponding function.

	Raises
	------
	SystemExit
		- If the hashing algorithm is not supported.
	"""
	# Get hash function
	hashing_function = getattr(hashlib, hashing_algo, None)
	if hashing_function is None:
		hashing_function = getattr(zlib, hashing_algo, None)

	if hashing_function is None:
		sys.exit(ct.HASHPROFILES_INVALID_HASHING.format(hashing_algo))

	return [hashing_algo, hashing_function]


# Annotated Fields
HashingAlgo = Annotated[str,
                        Field(default=ct.HASHPROFILES_DEFAULT_HASH_ALGO,
                              validate_default=True),
						AfterValidator(validate_hashing_algo)]

NrowChunk = Annotated[int,
                      Field(default=100)]


# Main validation class
class HashProfilesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	allelic_profiles: gval.AllelicProfiles
	schema_directory: gval.ChewieSchema
	hashing_algo: HashingAlgo
	nrow_chunk: NrowChunk
	cpu_cores: gval.CPUCores

	# Determine if all loci IDs in the profiles are in the schema
	@model_validator(mode="after")
	def profiles_loci_in_schema(self):
		# Read loci IDs in profiles file
		profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
		profiles_loci = profiles_header[0].split("\t")[1:]
		# Get list of loci in the schema
		schema_files = fo.listdir(self.schema_directory)
		schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
		schema_loci = [fo.file_basename(file, False) for file in schema_loci]

		extra_loci = set(profiles_loci) - set(schema_loci)
		if len(extra_loci) > 0:
			sys.exit("The following loci IDs in the input profiles are not in the schema:\n"
					 f"{'\n'.join(extra_loci)}\nPlease ensure that all loci IDs in the input profiles "
					 "are in the schema.")

		return self
