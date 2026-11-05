#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the AdaptSchema module.
"""

import sys
from typing import Annotated

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  model_validator)

try:
	import global_constants as gct
	import global_validation as gval
	import PredictGenes.validation_models as pval
	from utils import file_operations as fo
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct
	import CHEWBBACA.global_validation as gval
	import CHEWBBACA.PredictGenes.validation_models as pval
	from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
AdaptationMinimumLength = Annotated[int,
                                    Field(default=0)]

AdaptationSizeThreshold = Annotated[float,
                                    Field(default=None)]


# Main validation class
class AdaptSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	schema_directory: gval.ExternalSchema
	loci_list: gval.LociList | None = None
	gene_predictor: pval.GenePredictor
	gene_prediction_options: pval.GenePredictionOptions
	validated_gene_prediction_options: pval.ValidatedGenePredictionOptions
	clustering_options: gval.ClusteringOptions
	validated_clustering_options: gval.ValidatedClusteringOptions
	blast_score_ratio: gval.BLASTScoreRatio
	minimum_length: gval.MinimumLength
	adaptation_minimum_length: AdaptationMinimumLength
	translation_table: gval.TranslationTable
	size_threshold: gval.SizeThreshold
	adaptation_size_threshold: AdaptationSizeThreshold
	cpu_cores: gval.CPUCores
	blast_path: gval.BLASTPath
	size_filter: gval.BooleanFlag

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

	# Verify that the schema contains FASTA files for all loci in the list of loci
	@model_validator(mode="after")
	def loci_list_in_schema(self):
		if self.loci_list:
			# Compare list of loci IDs
			schema_files = fo.listdir(self.schema_directory)
			schema_loci = fo.filter_by_extension(schema_files, [".fasta"])[0]
			schema_loci = [fo.file_basename(file, False) for file in schema_loci]
			absent_files = set(self.loci_list) - set(schema_loci)
			if len(absent_files) > 0:
				sys.exit(gct.LOCI_LIST_SCHEMA_LOCI_NOMATCH.format(", ".join(absent_files)))
			# Get full paths
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	@model_validator(mode="after")
	def get_schema_loci(self):
		if not self.loci_list:
			# Get full paths
			schema_files = fo.listdir(self.schema_directory, full_path=True)
			loci_fpaths = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_basenames = [fo.file_basename(file, False) for file in loci_fpaths]
			self.loci_list = [loci_basenames, loci_fpaths]

		return self

	@model_validator(mode="after")
	def apply_size_filter(self):
		if self.size_filter:
			self.adaptation_minimum_length = self.minimum_length
			print(f"Using a minimum length value of {self.adaptation_minimum_length} for schema "
				  f"adaptation and {self.minimum_length} to store in the schema "
				  "config file.")
			self.adaptation_size_threshold = self.size_threshold
			print(f"Using a size threshold value of {self.adaptation_size_threshold} for schema "
		    	  f"adaptation and {self.size_threshold} to store in the schema "
		    	  "config file.")

		return self
