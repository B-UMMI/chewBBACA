#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the GetAlleles module.
"""


import sys
from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
                      Field,
                      ConfigDict,
                      AfterValidator,
					  model_validator)


import CHEWBBACA.GetAlleles.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
AllelicProfiles = Annotated[str,
                            Field(description=ct.GETALLELES_ALLELIC_PROFILES_DESCRIPTION),
                            AfterValidator(partial(gval.input_paths_exist,
                                                   function=fo.exists,
												   exception_str=gct.PROFILES_FILE_MISSING))]

OutputDirectory = Annotated[str,
                            Field(description=ct.GETALLELES_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

ChewieSchema = Annotated[str,
                         Field(description=ct.GETALLELES_SCHEMA_DIRECTORY_DESCRIPTION),
						 AfterValidator(partial(gval.input_paths_exist,
						                        function=fo.exists,
												exception_str=gct.SCHEMA_PATH_MISSING)),
						 AfterValidator(gval.schema_includes_fasta),
						 AfterValidator(gval.schema_check_genes_list),
						 AfterValidator(gval.schema_check_config),
						 AfterValidator(gval.schema_check_short)]

LociList = Annotated[str,
                     Field(default=None, description=ct.GETALLELES_LOCI_LIST_DESCRIPTION),
					 AfterValidator(partial(gval.input_paths_exist,
					                        function=fo.exists,
											exception_str=gct.MISSING_LOCI_LIST)),
					 AfterValidator(gval.input_is_fasta),
					 AfterValidator(fo.list_input_files)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.GETALLELES_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

TranslateAlleles = Annotated[bool,
                        Field(default=False,
                              description=ct.GETALLELES_TRANSLATE_ALLELES_DESCRIPTION)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT,
							       description=ct.GETALLELES_TRANSLATION_TABLE_DESCRIPTION),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

Distinct = Annotated[bool,
                        Field(default=False,
                              description=ct.GETALLELES_DISTINCT_DESCRIPTION)]


# Main validation class
class GetAllelesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: AllelicProfiles
	schema_directory: ChewieSchema
	loci_list: LociList
	cpu_cores: CPUCores
	translate_alleles: TranslateAlleles
	translation_table: TranslationTable
	distinct: Distinct

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
