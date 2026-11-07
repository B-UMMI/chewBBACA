#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the AdaptSchema module.
"""

import sys
from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

import CHEWBBACA.AdaptSchema.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
import CHEWBBACA.PredictGenes.validation_models as pval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
OutputDirectory = Annotated[str,
                            Field(description=ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

ExternalSchema = Annotated[str,
                           Field(description=ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_DESCRIPTION),
							AfterValidator(partial(gval.input_paths_exist,
							                       function=fo.exists,
												   exception_str=gct.SCHEMA_PATH_MISSING)),
							AfterValidator(gval.schema_includes_fasta)]

LociList = Annotated[str,
                     Field(description=ct.ADAPTSCHEMA_LOCI_LIST_DESCRIPTION),
					 AfterValidator(partial(gval.input_paths_exist,
					                        function=fo.exists,
											exception_str=gct.MISSING_LOCI_LIST)),
					 AfterValidator(gval.input_is_fasta),
					 AfterValidator(fo.list_input_files)]

BLASTScoreRatio = Annotated[float,
							Field(default=gct.BSR_DEFAULT,
							      ge=gct.BSR_MIN,
								  le=gct.BSR_MAX,
								  description=ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_DESCRIPTION)]

MinimumLength = Annotated[int,
						  Field(default=gct.MSL_DEFAULT,
						        ge=gct.MSL_MIN,
								le=gct.MSL_MAX,
								description=ct.ADAPTSCHEMA_MINIMUM_LENGTH_DESCRIPTION)]

AdaptationMinimumLength = Annotated[int,
                                    Field(default=0)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT,
							       description=ct.ADAPTSCHEMA_TRANSLATION_TABLE_DESCRIPTION),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

SizeThreshold = Annotated[float | None,
						  Field(default=gct.ST_DEFAULT,
						        ge=gct.ST_MIN,
								le=gct.ST_MAX,
								description=ct.ADAPTSCHEMA_SIZE_THRESHOLD_DESCRIPTION)]

AdaptationSizeThreshold = Annotated[float,
                                    Field(default=None)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.ADAPTSCHEMA_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

SizeFilter = Annotated[bool,
                        Field(default=False,
                              description=ct.ADAPTSCHEMA_SIZE_FILTER_DESCRIPTION)]


# Main validation class
class AdaptSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	schema_directory: ExternalSchema
	loci_list: LociList
	gene_predictor: pval.GenePredictor
	gene_prediction_options: pval.GenePredictionOptions
	validated_gene_prediction_options: pval.ValidatedGenePredictionOptions
	clustering_options: gval.ClusteringOptions
	validated_clustering_options: gval.ValidatedClusteringOptions
	blast_score_ratio: BLASTScoreRatio
	minimum_length: MinimumLength
	adaptation_minimum_length: AdaptationMinimumLength
	translation_table: TranslationTable
	size_threshold: SizeThreshold
	adaptation_size_threshold: AdaptationSizeThreshold
	cpu_cores: CPUCores
	blast_path: gval.BLASTPath
	size_filter: SizeFilter

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
