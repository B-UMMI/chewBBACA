#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the ComputeMSA module.
"""


from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  field_validator,
					  model_validator)

import CHEWBBACA.ComputeMSA.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
MSAInput = Annotated[str,
                     Field(description=ct.COMPUTEMSA_INPUT_PATH_DESCRIPTION),
					 AfterValidator(partial(gval.input_paths_exist,
								    function=fo.exists,
									exception_str=gct.MISSING_INPUT_ARG))]

OutputDirectory = Annotated[str,
                            Field(description=ct.COMPUTEMSA_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

ExternalSchema = Annotated[str | None,
                           Field(default=None, description=ct.COMPUTEMSA_SCHEMA_DIRECTORY_DESCRIPTION)]

DNAMSA = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTEMSA_DNA_MSA_DESCRIPTION)]

OutputVariable = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTEMSA_OUTPUT_VARIABLE_DESCRIPTION)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT,
							       description=ct.COMPUTEMSA_TRANSLATION_TABLE_DESCRIPTION),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.COMPUTEMSA_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

OnlyLociMSAs = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTEMSA_ONLY_LOCI_MSAS_DESCRIPTION)]

Gaps = Annotated[str,
                 Field(default=ct.DEFAULT_GAPS,
                       description=ct.COMPUTEMSA_GAPS_DESCRIPTION),
				 AfterValidator(partial(gval.validate_choice, valid_options=ct.GAPS_CHOICES, exception_str=ct.COMPUTEMSA_INVALID_GAPS_OPTION))]

Ambiguous = Annotated[str,
                      Field(default=ct.DEFAULT_AMBIGUOUS,
                            description=ct.COMPUTEMSA_AMBIGUOUS_DESCRIPTION),
					  AfterValidator(partial(gval.validate_choice, valid_options=ct.AMBIGUOUS_CHOICES, exception_str=ct.COMPUTEMSA_INVALID_AMBIGUOUS_OPTION))]

CustomMafftOptions = Annotated[str,
                               Field(default=None,
                                     description=ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_DESCRIPTION)]

ProteinInput = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTEMSA_PROTEIN_INPUT_DESCRIPTION)]


# Main validation class
class ComputeMSAValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_path: MSAInput
	schema_directory: ExternalSchema
	dna_msa: DNAMSA
	output_variable: OutputVariable
	translation_table: TranslationTable
	cpu_cores: CPUCores
	only_loci_msas: OnlyLociMSAs
	gaps: Gaps
	ambiguous: Ambiguous
	custom_mafft_options: CustomMafftOptions
	protein_input: ProteinInput
	no_cleanup: gval.NoCleanup

	@field_validator("schema_directory")
	@classmethod
	def schema_exists(cls, value):
		if value is not None:
			return gval.input_paths_exist(value, fo.exists, gct.SCHEMA_PATH_MISSING)

		return value

	@field_validator("schema_directory")
	@classmethod
	def schema_includes_fasta(cls, value):
		if value is not None:
			return gval.schema_includes_fasta(value)

		return value

	# Use the translation table defined in the training file if main input is a matrix with allelic profiles
	@model_validator(mode="after")
	def get_translation_table(self):
		if self.schema_directory:
			config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
			# Allow for schemas that were not created with chewBBACA
			if fo.exists(config_file):
				config_values = fo.pickle_loader(config_file)
				self.translation_table = config_values[gct.SCHEMA_CONFIG_TRANSLATION_TABLE_ARGNAME][0]

		return self
