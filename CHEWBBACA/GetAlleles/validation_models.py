#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the GetAlleles module.
"""


import sys
from pydantic import (BaseModel, ConfigDict, model_validator)


import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Main validation class
class GetAllelesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	allelic_profiles: gval.AllelicProfiles
	schema_directory: gval.ChewieSchema
	loci_list: gval.LociList | None = None
	cpu_cores: gval.CPUCores
	translate_alleles: gval.BooleanFlag
	translation_table: gval.TranslationTable
	distinct: gval.BooleanFlag

	# Use translation table defined in training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		# Determine if schema includes a training file
		schema_files = fo.listdir(self.schema_directory, full_path=True)
		training_file = fo.filter_by_extension(schema_files, [".trn"])[0]
		if len(training_file) > 0:
			self.translation_table = gval.add_ptf_genetic_code(training_file[0], self.translation_table)

		return self

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

	# Verify that the profiles include the loci IDs in the list of loci
	@model_validator(mode="after")
	def loci_list_in_profiles(self):
		if self.loci_list:
			# Read loci IDs in profiles file
			profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
			profiles_loci = profiles_header[0].split("\t")[1:]
			extra_loci = set(self.loci_list)-set(profiles_loci)
			if len(extra_loci) > 0:
				sys.exit("The following loci IDs in the input loci list are not in the allelic profiles:\n"
				         f"{'\n'.join(extra_loci)}\nPlease ensure that are all loci IDs in the input loci "
						 "list are in the allelic profiles.")

		return self

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
			schema_loci = [fo.file_basename(file, False) for file in schema_loci]
			extra_loci = set(self.loci_list) - set(schema_loci)
			if len(extra_loci) > 0:
				sys.exit(gct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(", ".join(extra_loci)))
			# Get full paths
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# Simply get the list of loci in the profiles if no loci list is provided
	@model_validator(mode="after")
	def get_profiles_loci(self):
		if not self.loci_list:
			# Read loci IDs in profiles file
			profiles_header = fo.read_lines(self.allelic_profiles, strip=True, num_lines=1)
			profiles_loci = profiles_header[0].split("\t")[1:]
			# Get full paths
			loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in profiles_loci]
			self.loci_list = [profiles_loci, loci_fpaths]

		return self
