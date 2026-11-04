#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


from pydantic import (BaseModel, ConfigDict, model_validator)

try:
	import global_constants as gct
	import global_validation as gval
	from utils import file_operations as fo
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct
	import CHEWBBACA.global_validation as gval
	from CHEWBBACA.utils import file_operations as fo


class AlleleCallEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	translation_table: gval.TranslationTable
	output_directory: gval.OutputDirectory
	schema_directory: gval.ChewieSchema
	results_files: gval.ResultsFiles
	annotations_file: gval.Annotations | None = None
	cpu_cores: gval.CPUCores
	light: gval.BooleanFlag
	no_pa: gval.BooleanFlag
	no_dm: gval.BooleanFlag
	no_tree: gval.BooleanFlag
	cg_alignment: gval.BooleanFlag

	@model_validator(mode="after")
	def get_translation_table(self):
		config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
		if fo.exists(config_file):
			config_values = fo.pickle_loader(config_file)
			self.translation_table = config_values[gct.TRANSLATION_TABLE_ARGNAME][0]

		return self
