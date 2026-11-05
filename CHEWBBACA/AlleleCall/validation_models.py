#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the AlleleCall module.
"""


import sys

from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

try:
	import AlleleCall.constants as ct
	import global_constants as gct
	import global_validation as gval
	import PredictGenes.validation_models as pval
	from utils import file_operations as fo
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct
	import CHEWBBACA.global_validation as gval
	from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
AllelecallMode = Annotated[int,
                           Field(default=ct.ALLELECALL_DEFAULT_MODE),
						   AfterValidator(partial(gval.validate_choice, valid_options=ct.ALLELECALL_MODES, exception_str=ct.ALLELECALL_INVALID_MODE))]


# Main validation class
class AlleleCallValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	input_files: gval.InputFiles
	schema_directory: gval.ChewieSchema
	# Do not run validation if loci_list is None
	loci_list: gval.LociList | None = None
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
	no_inferred: gval.BooleanFlag
	output_unclassified: gval.BooleanFlag
	output_missing: gval.BooleanFlag
	output_novel: gval.BooleanFlag
	output_masked_profiles: gval.BooleanFlag
	no_cds_renaming: gval.BooleanFlag
	force_continue: gval.BooleanFlag
	mode: AllelecallMode
	cpu_cores: gval.CPUCores
	no_cleanup: gval.BooleanFlag
	ns_config: gval.NSCONFIG

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
		schema_files = fo.listdir(self.schema_directory, full_path=True)
		loci_files = fo.filter_by_extension(schema_files, [".fasta"])[0]
		loci_ids = [fo.file_basename(file, False) for file in loci_files]
		if self.loci_list:
			self.loci_list = [self.loci_list, [loci_ids, loci_files]]
		else:
			self.loci_list = [[loci_ids, loci_files], [loci_ids, loci_files]]

		return self

	@model_validator(mode="after")
	def solve_conflicting_arguments(self):
		config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
		config = fo.pickle_loader(config_file)
		# BLAST Score Ratio
		self.blast_score_ratio, config = gval.check_argument_conflict(gct.SCHEMA_CONFIG_BLAST_SCORE_RATIO_ARGNAME, self.blast_score_ratio, config, self.force_continue)
		# Minimum length
		self.minimum_length, config = gval.check_argument_conflict(gct.SCHEMA_CONFIG_MINIMUM_LENGTH_ARGNAME, self.minimum_length, config, self.force_continue)
		# Size threshold
		self.size_threshold, config = gval.check_argument_conflict(gct.SCHEMA_CONFIG_SIZE_THRESHOLD_ARGNAME, self.size_threshold, config, self.force_continue)
		# Pyrodigal training file
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_file is not None:
			self.validated_gene_prediction_options.pyrodigal_training_file, config = gval.check_argument_conflict(gct.SCHEMA_CONFIG_PYRODIGAL_TRAININGFILE_ARGNAME, self.validated_gene_prediction_options.pyrodigal_training_file, config, self.force_continue)
			# Compare current value for the translation table and the one in the training file
			self.translation_table = gval.add_ptf_genetic_code(self.validated_gene_prediction_options.pyrodigal_training_file, self.translation_table)
			# Add translation table to config if it is different than the ones previously used
			if self.translation_table not in config[gct.SCHEMA_CONFIG_TRANSLATION_TABLE_ARGNAME]:
				config[gct.SCHEMA_CONFIG_TRANSLATION_TABLE_ARGNAME].append(self.translation_table)

		# Update schema config file
		fo.pickle_dumper(config, config_file)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options.pyrodigal_mode = gval.check_meta(self.validated_gene_prediction_options.pyrodigal_mode, self.validated_gene_prediction_options.pyrodigal_training_file)

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [gct.NS_CONFIG_BASENAME])
		if fo.is_file(ns_config)[0]:
			self.ns_config = True

		return self
