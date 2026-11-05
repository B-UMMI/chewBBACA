#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the AnnotateSchema module.
"""


import sys
from typing import Annotated
from functools import partial

import pandas as pd
from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

import CHEWBBACA.AnnotateSchema.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
ProteinTable = Annotated[str,
                         Field(default=None),
                         AfterValidator(partial(gval.input_paths_exist,
                                                function=fo.exists,
												exception_str=ct.MISSING_PROTEIN_TABLE))]

Taxa = Annotated[list,
                 Field(default=None)]

ProteomeMatches = Annotated[int,
                            Field(default=1, ge=1)]

RestApiMatches = Annotated[int,
                           Field(default=1, ge=1)]


# Main validation class
class AnnotateSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	schema_directory: gval.ExternalSchema
	loci_list: gval.LociList | None = None
	protein_table: ProteinTable
	blast_score_ratio: gval.BLASTScoreRatio
	translation_table: gval.TranslationTable
	cpu_cores: gval.CPUCores
	proteome_taxa: Taxa
	proteome_matches: ProteomeMatches
	blast_path: gval.BLASTPath
	rest_api: gval.BooleanFlag
	rest_api_matches: RestApiMatches
	no_cleanup: gval.BooleanFlag

	# User must specify that loci should be annotated based on either proteomes or requests to UniProt's Proteins API or both
	@model_validator(mode="after")
	def check_annotation_method(self):
		if not self.proteome_taxa and not self.rest_api:
			sys.exit(gct.ANNOTATE_SCHEMA_NO_ANNOTATION_METHOD)

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
			# Get full paths
			else:
				loci_fpaths = [fo.join_paths(self.schema_directory, [locus_id+".fasta"]) for locus_id in self.loci_list]
				self.loci_list = [self.loci_list, loci_fpaths]

		return self

	# Annotate all loci if user did not provide a list of loci
	@model_validator(mode="after")
	def get_schema_loci(self):
		if not self.loci_list:
			# Get full paths
			schema_files = fo.listdir(self.schema_directory, full_path=True)
			loci_fpaths = fo.filter_by_extension(schema_files, [".fasta"])[0]
			loci_basenames = [fo.file_basename(file, False) for file in loci_fpaths]
			self.loci_list = [loci_basenames, loci_fpaths]

		return self

	# Use translation table stored in schema config file if it is a chewie schema
	@model_validator(mode="after")
	def get_translation_table(self):
		config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
		if fo.exists(config_file):
			config_values = fo.pickle_loader(config_file)
			self.translation_table = config_values[gct.TRANSLATION_TABLE_ARGNAME][0]

		return self
