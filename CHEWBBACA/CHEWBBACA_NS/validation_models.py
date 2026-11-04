#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


import sys

from typing import Annotated
from functools import partial

import pandas as pd
from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  TypeAdapter,
					  AfterValidator,
					  model_validator)

try:
	import global_constants as gct
	import global_validation as gval
	from utils import file_operations as fo
	from utils import chewiens_requests as cr
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct
	import CHEWBBACA.global_validation as gval
	from CHEWBBACA.utils import file_operations as fo
	from CHEWBBACA.utils import chewiens_requests as cr


SpeciesId = Annotated[str,
                      Field(default=None)]

SchemaId = Annotated[str,
                     Field(default=None)]

NomenclatureServerInstance = Annotated[str | None,
                               		   Field(default=gct.HOST_NS[gct.DEFAULT_NOMENCLATURE_SERVER])]

SchemaDate = Annotated[str,
                       Field(default=None),
					   AfterValidator(pd.validate_date)] # Determine if user date is in a valid format

DescriptionFile = Annotated[str,
                            Field(default=None),
                            AfterValidator(partial(gval.input_paths_exist,
                                                   function=fo.exists,
												   exception_str=gct.MISSING_INPUTS_EXCEPTION))]

StatsMode = Annotated[str,
                      Field(default=None),
					  AfterValidator(partial(gval.validate_choice,
					                         valid_options=gct.NSSTATS_MODES,
											 exception_str=gct.NS_STATS_INVALID_MODE))]


class DownloadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	species_id: gval.SpeciesId
	schema_id: gval.SchemaId
	nomenclature_server_instance: gval.NomenclatureServerInstance
	schema_date: gval.SchemaDate
	latest_version: gval.BooleanFlag
	blast_path : gval.BLASTPath
	cpu_cores: gval.CPUCores

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self

	# Check if user provided a valid integer identifier or species name
	@model_validator(mode="after")
	def validate_species_id(self):
		species_info = cr.species_ids(self.species_id, self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if isinstance(species_info, list):
			self.species_id = species_info
		else:
			sys.exit(f"There is no species with the identifier {self.species_id} in the specified "
			         f"Chewie-NS instance ({self.nomenclature_server_instance}).")

		return self

	@model_validator(mode="after")
	def validate_schema_id(self):
		# The function return a list with the schema ID, the schema URI, the schema name, and the schema parameteres
		schema_info = cr.get_species_schemas(self.schema_id, self.species_id, self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		self.schema_id = schema_info

		return self


class UploadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: gval.ChewieSchema
	species_id: gval.SpeciesId
	schema_name: str
	loci_prefix: str
	description_file: gval.DescriptionFile
	annotations_file: gval.Annotations | None = None
	nomenclature_server_instance: gval.NomenclatureServerInstance
	continue_up: gval.BooleanFlag
	cpu_cores: gval.CPUCores

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self

	# Check if user provided a valid integer identifier or species name
	@model_validator(mode="after")
	def validate_species_id(self):
		species_info = cr.species_ids(self.species_id, self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if isinstance(species_info, list):
			self.species_id = species_info
		else:
			sys.exit(f"There is no species with the identifier {self.species_id} in the specified "
			         f"Chewie-NS instance ({self.nomenclature_server_instance}).")

		return self

	# Read config file and validate values
	@model_validator(mode="after")
	def validate_schema_config(self):
		print("Validating the argument values stored in the schema config...")
		config_file = fo.join_paths(self.schema_directory, [gct.SCHEMA_CONFIG_BASENAME])
		config_values = fo.pickle_loader(config_file)
		if any([len(v) > 1 for k, v in config_values.items()]):
			sys.exit("Schemas can only be uploaded to Chewie-NS if they were used with a "
			         "single argument per parameter defined in the '.schema_config file'.")
		# Validate arguments if there is a single value per parameter
		# Pyrodigal training file
		ptf_hash = config_values["pyrodigal_training_file"][0]
		schema_files = fo.listdir(self.schema_directory, full_path=True)
		ptf_paths = fo.filter_by_extension(schema_files, [".trn"])[0]
		if len(ptf_paths) > 1:
			sys.exit(f"Schema includes multiple training files. Please ensure that the schema includes a single training file.")

		if fo.hash_file(ptf_paths[0], "blake2b") != ptf_hash:
			sys.exit(f"Training file in schema's directory does not match the one in the config file.")

		# Validate remaining argument values using Pydantic's TypeAdapter
		bsr_adapter = TypeAdapter(gval.BLASTScoreRatio).validate_python(config_values["blast_score_ratio"][0])
		ml_adapter = TypeAdapter(gval.MinimumLength).validate_python(config_values["minimum_length"][0])
		st_adapter = TypeAdapter(gval.SizeThreshold).validate_python(config_values["size_threshold"][0])
		tt_adapter = TypeAdapter(gval.TranslationTable).validate_python(config_values["translation_table"][0])
		wos_adapter = TypeAdapter(gval.WordSize).validate_python(config_values["word_size"][0])
		win_adapter = TypeAdapter(gval.WindowSize).validate_python(config_values["window_size"][0])
		cs_adapter = TypeAdapter(gval.ClusteringSimilarity).validate_python(config_values["clustering_sim"][0])
		rf_adapter = TypeAdapter(gval.RepresentativeFilter).validate_python(config_values["representative_filter"][0])
		if_adapter = TypeAdapter(gval.IntraFilter).validate_python(config_values["intra_filter"][0])

		validated_config = {
			"bsr": config_values["blast_score_ratio"][0],
			"minimum_locus_length": config_values["minimum_length"][0],
			"size_threshold": config_values["size_threshold"][0],
			"translation_table": config_values["translation_table"][0],
			"prodigal_training_file": [ptf_paths[0], ptf_hash],
			"word_size": config_values["word_size"][0],
			"cluster_sim": config_values["clustering_sim"][0],
			"representative_filter": config_values["representative_filter"][0],
			"intraCluster_filter": config_values["intra_filter"][0],
			"chewBBACA_version": f"chewBBACA {config_values['chewie_version'][0]}"
		}

		self.schema_directory = [self.schema_directory, validated_config]

		return self


class SynchronizeSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: gval.ChewieSchema
	nomenclature_server_instance: gval.NomenclatureServerInstance
	submit_alleles: gval.BooleanFlag
	blast_path: gval.BLASTPath
	cpu_cores: gval.CPUCores

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self

	# Check if schema was downloaded from Chewie-NS
	@model_validator(mode="after")
	def check_ns_config(self):
		ns_config = fo.join_paths(self.schema_directory, [gct.NS_CONFIG_BASENAME])
		if not fo.is_file(ns_config)[0]:
			sys.exit("Schema does not contain the '.ns_config' file.")

		return self


class NSStatsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	stats_mode: gval.StatsMode
	species_id: gval.SpeciesId
	schema_id: gval.SchemaId
	nomenclature_server_instance: gval.NomenclatureServerInstance

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self
