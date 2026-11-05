#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the SchemaEvaluator module.
"""


import sys

from pydantic import (BaseModel, ConfigDict, model_validator)

import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Main validation class
class SchemaEvaluatorValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	schema_directory: gval.ChewieSchema
	loci_list: gval.LociList | None = None
	annotations_file: gval.Annotations | None = None
	translation_table: gval.TranslationTable
	size_threshold: gval.SizeThreshold
	minimum_length: gval.MinimumLength
	cpu_cores: gval.CPUCores
	loci_reports: gval.BooleanFlag
	light: gval.BooleanFlag
	add_sequences: gval.BooleanFlag
	chewie_schema: gval.BooleanFlag
	minimum_length_values: gval.MinimumLengthValues
	size_threshold_values: gval.SizeThresholdValues
	translation_table_values: gval.TranslationTableValues
	blast_score_ratio_values: gval.BlastScoreRatioValues
	training_file: gval.BooleanFlag
	ns_config: gval.NSCONFIG

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
			self.chewie_schema = config_values.get(gct.CHEWIE_VERSION_ARGNAME)[0]
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
			self.minimum_length = config_values.get(gct.MINIMUM_LENGTH_ARGNAME)[0]
			self.size_threshold = config_values.get(gct.SIZE_THRESHOLD_ARGNAME)[0]
			self.translation_table = config_values.get(gct.TRANSLATION_TABLE_ARGNAME)[0]

		return self

	@model_validator(mode="after")
	def set_schema_args(self):
		if self.chewie_schema:
			config_file  = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
			config_values = fo.pickle_loader(config_file)
			# Get the first minimum_length, size_threshold, and translation_table values used with the schema
			self.minimum_length_values = config_values.get(gct.MINIMUM_LENGTH_ARGNAME)
			self.size_threshold_values = config_values.get(gct.SIZE_THRESHOLD_ARGNAME)
			self.translation_table_values = config_values.get(gct.TRANSLATION_TABLE_ARGNAME)
			self.blast_score_ratio_values = config_values.get(gct.BLAST_SCORE_RATIO_ARGNAME)
			self.training_file = True if config_values.get(gct.PYRODIGAL_TRAININGFILE_ARGNAME)[0] else False

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [gct.NS_CONFIG_BASENAME])
		if fo.is_file(ns_config)[0]:
			self.ns_config = True

		return self
