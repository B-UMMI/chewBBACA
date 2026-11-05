#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

try:
	import ComputeMSA.constants as ct
	import global_constants as gct
	import global_validation as gval
	from utils import file_operations as fo
except ModuleNotFoundError:
	import CHEWBBACA.ComputeMSA.constants as ct
	import CHEWBBACA.global_constants as gct
	import CHEWBBACA.global_validation as gval
	from CHEWBBACA.utils import file_operations as fo


MSAInput = Annotated[str,
					 AfterValidator(partial(gval.input_paths_exist,
								    function=fo.exists,
									exception_str=gct.MISSING_INPUT_ARG))]

Gaps = Annotated[str,
                 Field(default=ct.DEFAULT_GAPS),
				 AfterValidator(partial(gval.validate_choice, valid_options=ct.GAPS_CHOICES, exception_str=ct.COMPUTEMSA_INVALID_GAPS_OPTION))]

Ambiguous = Annotated[str,
                      Field(default=ct.DEFAULT_AMBIGUOUS),
					  AfterValidator(partial(gval.validate_choice, valid_options=ct.AMBIGUOUS_CHOICES, exception_str=ct.COMPUTEMSA_INVALID_AMBIGUOUS_OPTION))]

CustomMafftOptions = Annotated[str,
                               Field(default=None)]


class ComputeMSAValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	input_path: MSAInput
	schema_directory: gval.ExternalSchema
	dna_msa: gval.BooleanFlag
	output_variable: gval.BooleanFlag
	translation_table: gval.TranslationTable
	cpu_cores: gval.CPUCores
	only_loci_msas: gval.BooleanFlag
	gaps: Gaps
	ambiguous: Ambiguous
	custom_mafft_options: CustomMafftOptions
	protein_input: gval.BooleanFlag
	no_cleanup: gval.BooleanFlag

	# Use the translation table defined in the training file if main input is a matrix with allelic profiles
	@model_validator(mode="after")
	def get_translation_table(self):
		if self.schema_directory:
			config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
			# Allow for schemas that were not created with chewBBACA
			if fo.exists(config_file):
				config_values = fo.pickle_loader(config_file)
				self.translation_table = config_values[gct.TRANSLATION_TABLE_ARGNAME][0]

		return self
