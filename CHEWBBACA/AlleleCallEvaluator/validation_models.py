#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the AlleleCallEvaluator module.
"""


from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
                      Field,
                      ConfigDict,
                      AfterValidator,
					  model_validator)

import CHEWBBACA.AlleleCallEvaluator.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
# Nedd to add a AfterValidator to get and validate results files
ResultsFiles = Annotated[str,
                         Field(description=ct.ALLELECALLEVALUATOR_RESULTS_FILES_DESCRIPTION),
                         AfterValidator(partial(gval.input_paths_exist,
                                                function=fo.exists,
												exception_str=gct.MISSING_RESULTS)),
						 AfterValidator(partial(gval.contains_results,
						                        results_files_basenames=[gct.ALLELIC_PROFILES_BASENAME]))]

ChewieSchema = Annotated[str,
                         Field(description=ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_DESCRIPTION),
						 AfterValidator(partial(gval.input_paths_exist,
						                        function=fo.exists,
												exception_str=gct.SCHEMA_PATH_MISSING)),
						 AfterValidator(gval.schema_includes_fasta),
						 AfterValidator(gval.schema_check_genes_list),
						 AfterValidator(gval.schema_check_config),
						 AfterValidator(gval.schema_check_short)]

OutputDirectory = Annotated[str,
                            Field(description=ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

Annotations = Annotated[str,
                        Field(default=None, description=ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_DESCRIPTION),
                        AfterValidator(partial(gval.input_paths_exist,
                                               function=fo.exists,
											   exception_str=gct.MISSING_ANNOTATIONS_FILE))]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.ALLELECALLEVALUATOR_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

NoPA = Annotated[bool,
                 Field(default=False,
                       description=ct.ALLELECALLEVALUATOR_NO_PRESENCE_ABSENCE_DESCRIPTION)]

NoDM = Annotated[bool,
                 Field(default=False,
                       description=ct.ALLELECALLEVALUATOR_NO_DISTANCE_MATRIX_DESCRIPTION)]

NoTree = Annotated[bool,
                   Field(default=False,
                         description=ct.ALLELECALLEVALUATOR_NO_NEIGHBOR_JOINING_DESCRIPTION)]

CGAlignment = Annotated[bool,
                         Field(default=False,
                               description=ct.ALLELECALLEVALUATOR_FORCE_CORE_MSA_DESCRIPTION)]

Light = Annotated[bool,
                   Field(default=False,
                         description=ct.ALLELECALLEVALUATOR_LIGTH_DESCRIPTION)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]


# Main validation class
class AlleleCallEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	results_files: ResultsFiles
	schema_directory: ChewieSchema
	output_directory: OutputDirectory
	annotations_file: Annotations
	cpu_cores: CPUCores
	no_presence_absence: NoPA
	no_distance_matrix: NoDM
	no_neighbor_joining: NoTree
	force_core_msa: CGAlignment
	light: Light
	translation_table: TranslationTable

	@model_validator(mode="after")
	def get_translation_table(self):
		config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
		if fo.exists(config_file):
			config_values = fo.pickle_loader(config_file)
			self.translation_table = config_values[ct.TRANSLATION_TABLE_ARGNAME][0]

		return self
