#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the CreateSchema module.
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
import CHEWBBACA.global_validation as gval
import CHEWBBACA.utils.file_operations as fo
import CHEWBBACA.CreateSchema.constants as ct
import CHEWBBACA.PredictGenes.validation_models as pval


# Annotated Fields
OutputDirectory = Annotated[str,
                            Field(description=ct.CREATESCHEMA_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

# Need to include list type to avoid PydanticSerializationUnexpectedValue
# Input is a str but AfterValidator converts it into a list of valid file paths
InputPath = Annotated[str | list,
                      Field(description=ct.CREATESCHEMA_INPUT_PATH_DESCRIPTION),
					   AfterValidator(partial(gval.input_paths_exist,
					                          function=fo.exists,
											  exception_str=gct.MISSING_INPUT_ARG)),
					   AfterValidator(gval.input_is_fasta),
					   AfterValidator(fo.list_input_files),
					   AfterValidator(gval.filter_inputs_extension),
					   AfterValidator(gval.inputs_exist),
					   AfterValidator(gval.validate_inputs_fastas),
					   AfterValidator(gval.check_unique_prefixes),
					   AfterValidator(gval.check_blanks)]

SchemaName = Annotated[str,
                       Field(default=ct.SCHEMA_NAME_DEFAULT,
                             description=ct.CREATESCHEMA_SCHEMA_NAME_DESCRIPTION)]

BLASTScoreRatio = Annotated[float,
							Field(default=gct.BSR_DEFAULT,
							      ge=gct.BSR_MIN,
								  le=gct.BSR_MAX,
								  description=ct.CREATESCHEMA_BLAST_SCORE_RATIO_DESCRIPTION)]

MinimumLength = Annotated[int,
						  Field(default=gct.MSL_DEFAULT,
						        ge=gct.MSL_MIN,
								le=gct.MSL_MAX,
								description=ct.CREATESCHEMA_MINIMUM_LENGTH_DESCRIPTION)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT,
							       description=ct.CREATESCHEMA_TRANSLATION_TABLE_DESCRIPTION),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

SizeThreshold = Annotated[float | None,
						  Field(default=gct.ST_DEFAULT,
						        ge=gct.ST_MIN,
								le=gct.ST_MAX,
								description=ct.CREATESCHEMA_SIZE_THRESHOLD_DESCRIPTION)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.CREATESCHEMA_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]


# Main validation class
class CreateSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	input_path: InputPath
	schema_name: SchemaName
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	gene_predictor : pval.GenePredictor
	gene_prediction_options: pval.GenePredictionOptions
	validated_gene_prediction_options: pval.ValidatedGenePredictionOptions
	clustering_options: gval.ClusteringOptions
	validated_clustering_options: gval.ValidatedClusteringOptions
	blast_path: gval.BLASTPath
	cds_input: gval.CDSInput
	no_cds_renaming: gval.NoCDSRenaming
	cpu_cores: CPUCores
	no_cleanup: gval.NoCleanup

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options = pval.PyrodigalArgs(**self.gene_prediction_options)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_options = pval.AugustusArgs(**self.gene_prediction_options)

		return self

	# Further validation for clustering arguments
	@model_validator(mode="after")
	def validate_clustering_arguments(self):
		self.validated_clustering_options = gval.ClusteringArgs(**self.clustering_options)

		return self

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_file is not None:
			self.translation_table = gval.add_ptf_genetic_code(self.validated_gene_prediction_options.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options.pyrodigal_mode = gval.check_meta(self.validated_gene_prediction_options.pyrodigal_mode, self.validated_gene_prediction_options.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def ptf_and_reference(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_reference:
			# Check if user provided path to Pyrodigal training file and to a training reference, which is not allowed
			gval.check_ptf_tref_conflict(self.validated_gene_prediction_options.pyrodigal_training_file, self.validated_gene_prediction_options.pyrodigal_training_reference)
			# Create training file based on training reference
			print(f'Creating Pyrodigal training file based on {self.validated_gene_prediction_options.pyrodigal_training_reference}...')
			self.validated_gene_prediction_options.pyrodigal_training_file = gval.create_training_file(self.validated_gene_prediction_options.pyrodigal_training_reference, self.output_directory, self.translation_table)
			print(f'Training file saved to {self.validated_gene_prediction_options.pyrodigal_training_file}')
			# Exit if user only requested to create training file
			if self.validated_gene_prediction_options.pyrodigal_just_training:
				sys.exit(gct.JUST_TRAINING)

		return self
