#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the SchemaEvaluator module.
"""


import sys
from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
                      Field,
                      ConfigDict,
                      AfterValidator,
					  model_validator)

import CHEWBBACA.global_constants as gct
import CHEWBBACA.SchemaEvaluator.constants as ct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
OutputDirectory = Annotated[str,
                            Field(description=ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

ExternalSchema = Annotated[str,
                           Field(description=ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_DESCRIPTION),
							AfterValidator(partial(gval.input_paths_exist,
							                       function=fo.exists,
												   exception_str=gct.SCHEMA_PATH_MISSING)),
							AfterValidator(gval.schema_includes_fasta)]

LociList = Annotated[str,
                     Field(default=None, description=ct.SCHEMAEVALUATOR_LOCI_LIST_DESCRIPTION),
					 AfterValidator(partial(gval.input_paths_exist,
					                        function=fo.exists,
											exception_str=gct.MISSING_LOCI_LIST)),
					 AfterValidator(gval.input_is_fasta),
					 AfterValidator(fo.list_input_files)]

Annotations = Annotated[str,
                        Field(default=None, description=ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_DESCRIPTION),
                        AfterValidator(partial(gval.input_paths_exist,
                                               function=fo.exists,
											   exception_str=gct.MISSING_ANNOTATIONS_FILE))]

MinimumLength = Annotated[int,
						  Field(default=gct.MSL_DEFAULT,
						        ge=gct.MSL_MIN,
								le=gct.MSL_MAX,
								description=ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_DESCRIPTION)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT,
							       description=ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_DESCRIPTION),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

SizeThreshold = Annotated[float | None,
						  Field(default=gct.ST_DEFAULT,
						        ge=gct.ST_MIN,
								le=gct.ST_MAX,
								description=ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_DESCRIPTION)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.SCHEMAEVALUATOR_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

MinimumLengthValues = Annotated[list,
                                Field(default=[])]

SizeThresholdValues = Annotated[list,
                                Field(default=[])]

TranslationTableValues = Annotated[list,
                                   Field(default=[])]

BlastScoreRatioValues = Annotated[list,
                                  Field(default=[])]

LociReports = Annotated[bool,
                     Field(default=False,
                           description=ct.SCHEMAEVALUATOR_LOCI_REPORTS_DESCRIPTION)]

Light = Annotated[bool,
                     Field(default=False,
                           description=ct.SCHEMAEVALUATOR_LIGTH_DESCRIPTION)]

AddSequences = Annotated[bool,
                     Field(default=False,
                           description=ct.SCHEMAEVALUATOR_ADD_SEQUENCES_DESCRIPTION)]

ChewieSchema = Annotated[bool,
                     Field(default=False)]

TrainingFile = Annotated[bool,
                     Field(default=False)]

NSConfig = Annotated[bool,
                     Field(default=False)]


# Main validation class
class SchemaEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: ExternalSchema
	loci_list: LociList
	annotations_file: Annotations
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	minimum_length: MinimumLength
	cpu_cores: CPUCores
	loci_reports: LociReports
	light: Light
	add_sequences: AddSequences
	chewie_schema: ChewieSchema
	minimum_length_values: MinimumLengthValues
	size_threshold_values: SizeThresholdValues
	translation_table_values: TranslationTableValues
	blast_score_ratio_values: BlastScoreRatioValues
	training_file: TrainingFile
	ns_config: NSConfig

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_ids = [fo.file_basename(file, False) for file in loci_files]
			extra_loci = set(self.loci_list) - set(loci_ids)
			if len(extra_loci) > 0:
				sys.exit(gct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(", ".join(extra_loci)))
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# List all loci FASTA files in the schema and add it to the loci_list arg
	@model_validator(mode="after")
	def list_schema_loci(self):
		if not self.loci_list:
			schema_files = fo.listdir(self.schema_directory, full_path=True)
			loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_ids = [fo.file_basename(file, False) for file in loci_files]
			self.loci_list = [loci_ids, loci_files]

		return self

	@model_validator(mode="after")
	def set_schema_origin(self):
		schema_files = fo.listdir(self.schema_directory, full_path=False)
		# Schema was created with chewBBACA
		if gct.SCHEMA_CONFIG_BASENAME in schema_files:
			config_file  = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Set the argument that indicates that the schema was created with chewBBACA to True
			self.chewie_schema = config_values.get(gct.SCHEMA_CONFIG_CHEWIE_VERSION_ARGNAME)[0]
		else:
			self.chewie_schema = False

		return self

	# Use the values stored in the schema config if it is a chewBBACA schema
	@model_validator(mode="after")
	def force_schema_args(self):
		if self.chewie_schema:
			print("Schema was created with chewBBACA. Using the values stored in the schema "
			      "config file for the minimum length, size threshold, and translation table.")
			config_file  = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Get the first minimum_length, size_threshold, and translation_table values used with the schema
			self.minimum_length = config_values.get(gct.SCHEMA_CONFIG_MINIMUM_LENGTH_ARGNAME)[0]
			self.size_threshold = config_values.get(gct.SCHEMA_CONFIG_SIZE_THRESHOLD_ARGNAME)[0]
			self.translation_table = config_values.get(gct.SCHEMA_CONFIG_TRANSLATION_TABLE_ARGNAME)[0]

		return self

	@model_validator(mode="after")
	def set_schema_args(self):
		if self.chewie_schema:
			config_file  = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Get the first minimum_length, size_threshold, and translation_table values used with the schema
			self.minimum_length_values = config_values.get(gct.SCHEMA_CONFIG_TRANSLATION_TABLE_ARGNAME)
			self.size_threshold_values = config_values.get(gct.SCHEMA_CONFIG_SIZE_THRESHOLD_ARGNAME)
			self.translation_table_values = config_values.get(gct.SCHEMA_CONFIG_TRANSLATION_TABLE_ARGNAME)
			self.blast_score_ratio_values = config_values.get(gct.SCHEMA_CONFIG_BLAST_SCORE_RATIO_ARGNAME)
			self.training_file = True if config_values.get(gct.SCHEMA_CONFIG_PYRODIGAL_TRAININGFILE_ARGNAME)[0] else False

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [gct.NS_CONFIG_BASENAME])
		if fo.is_file(ns_config)[0]:
			self.ns_config = True

		return self
