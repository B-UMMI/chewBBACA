#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the CreateSchema module.
"""


import sys
from typing import Annotated

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  model_validator)

import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
import CHEWBBACA.CreateSchema.constants as ct
import CHEWBBACA.PredictGenes.validation_models as pval


# Annotated Fields
SchemaName = Annotated[str,
                       Field(default=ct.SCHEMA_NAME_DEFAULT)]


# Main validation class
class CreateSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	input_files: gval.InputFiles
	schema_name: SchemaName
	blast_score_ratio: gval.BLASTScoreRatio
	minimum_length: gval.MinimumLength
	translation_table: gval.TranslationTable
	size_threshold: gval.SizeThreshold
	gene_predictor : pval.GenePredictor
	gene_prediction_options: pval.GenePredictionOptions
	validated_gene_prediction_options: pval.ValidatedGenePredictionOptions
	clustering_options: gval.ClusteringOptions
	validated_clustering_options: gval.ValidatedClusteringOptions
	blast_path: gval.BLASTPath
	cds_input: gval.BooleanFlag
	no_cds_renaming: gval.BooleanFlag
	cpu_cores: gval.CPUCores
	no_cleanup: gval.BooleanFlag

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
