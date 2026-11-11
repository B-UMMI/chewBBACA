#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the modules that allow for the integration with Chewie-NS.
"""


import sys

from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  TypeAdapter,
					  AfterValidator,
					  model_validator)

import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
import CHEWBBACA.CHEWBBACA_NS.constants as ct
from CHEWBBACA.utils import file_operations as fo
from CHEWBBACA.utils import chewiens_requests as cr
from CHEWBBACA.utils.process_datetime import validate_date


# DownloadSchema Annotated Fields
DownloadSchemaOutputDirectory = Annotated[str,
                            Field(description=ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

DownloadSchemaSpeciesId = Annotated[str | list,
                      Field(default=None,
                            description=ct.DOWNLOADSCHEMA_SPECIES_ID_DESCRIPTION)]

DownloadSchemaSchemaId = Annotated[str | list,
                     Field(default=None,
                           description=ct.DOWNLOADSCHEMA_SCHEMA_ID_DESCRIPTION)]

DownloadSchemaNomenclatureServerInstance = Annotated[str | None,
                               		   Field(default=ct.HOST_NS[ct.DEFAULT_NOMENCLATURE_SERVER],
                                                description=ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION)]

DownloadSchemaSchemaDate = Annotated[str,
                       Field(default=None,
                             description=ct.DOWNLOADSCHEMA_SCHEMA_DATE_DESCRIPTION),
					   AfterValidator(validate_date)] # Determine if user date is in a valid format

DownloadSchemaLatestVersion = Annotated[bool,
                        Field(default=False,
                              description=ct.DOWNLOADSCHEMA_LATEST_VERSION_DESCRIPTION)]

DownloadSchemaCPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.DOWNLOADSCHEMA_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]


# Main validation class for the DownloadSchema module
class DownloadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: DownloadSchemaOutputDirectory
	species_id: DownloadSchemaSpeciesId
	schema_id: DownloadSchemaSchemaId
	nomenclature_server_instance: DownloadSchemaNomenclatureServerInstance
	schema_date: DownloadSchemaSchemaDate
	latest_version: DownloadSchemaLatestVersion
	blast_path : gval.BLASTPath
	cpu_cores: DownloadSchemaCPUCores

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


# UploadSchema Annotated Fields
UploadSchemaSchemaDirectory = Annotated[str,
                                        Field(description=ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_DESCRIPTION),
						 AfterValidator(partial(gval.input_paths_exist,
						                        function=fo.exists,
												exception_str=gct.SCHEMA_PATH_MISSING)),
						 AfterValidator(gval.schema_includes_fasta),
						 AfterValidator(gval.schema_check_genes_list),
						 AfterValidator(gval.schema_check_config),
						 AfterValidator(gval.schema_check_short)]

UploadSchemaSpeciesId = Annotated[str,
                      Field(default=None,
                            description=ct.UPLOADSCHEMA_SPECIES_ID_DESCRIPTION)]

UploadSchemaSchemaName = Annotated[str,
                       Field(description=ct.UPLOADSCHEMA_SCHEMA_NAME_DESCRIPTION)]

UploadSchemaLociPrefix = Annotated[str,
                       Field(description=ct.UPLOADSCHEMA_LOCI_PREFIX_DESCRIPTION)]

UploadSchemaDescriptionFile = Annotated[str,
                            Field(default=None,
                                  description=ct.UPLOADSCHEMA_DESCRIPTION_FILE_DESCRIPTION),
                            AfterValidator(partial(gval.input_paths_exist,
                                                   function=fo.exists,
												   exception_str=gct.MISSING_INPUTS_EXCEPTION))]

UploadSchemaAnnotations = Annotated[str,
                                    Field(description=ct.UPLOADSCHEMA_ANNOTATIONS_FILE_DESCRIPTION),
                        			AfterValidator(partial(gval.input_paths_exist,
                                               function=fo.exists,
											   exception_str=gct.MISSING_ANNOTATIONS_FILE))]

UploadSchemaNomenclatureServerInstance = Annotated[str | None,
                               		   Field(default=ct.HOST_NS[ct.DEFAULT_NOMENCLATURE_SERVER],
                                                description=ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION)]

UploadSchemaContinueUp = Annotated[bool,
                       Field(description=ct.UPLOADSCHEMA_CONTINUE_UPLOAD_DESCRIPTION)]

UploadSchemaCPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.UPLOADSCHEMA_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]


# Main validation class for the UploadSchema module
class UploadSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: UploadSchemaSchemaDirectory
	species_id: UploadSchemaSpeciesId
	schema_name: UploadSchemaSchemaName
	loci_prefix: UploadSchemaLociPrefix
	description_file: UploadSchemaDescriptionFile
	annotations_file: UploadSchemaAnnotations
	nomenclature_server_instance: UploadSchemaNomenclatureServerInstance
	continue_upload: UploadSchemaContinueUp
	cpu_cores: UploadSchemaCPUCores

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


# SyncSchema Annotated Fields
SyncSchemaSchemaDirectory = Annotated[str,
                                        Field(description=ct.SYNCSCHEMA_SCHEMA_DIRECTORY_DESCRIPTION),
						 AfterValidator(partial(gval.input_paths_exist,
						                        function=fo.exists,
												exception_str=gct.SCHEMA_PATH_MISSING)),
						 AfterValidator(gval.schema_includes_fasta),
						 AfterValidator(gval.schema_check_genes_list),
						 AfterValidator(gval.schema_check_config),
						 AfterValidator(gval.schema_check_short)]

SyncSchemaNomenclatureServerInstance = Annotated[str | None,
                               		   Field(default=ct.HOST_NS[ct.DEFAULT_NOMENCLATURE_SERVER],
                                                description=ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION)]

SyncSchemaSubmitAlleles = Annotated[bool,
                          Field(description=ct.SYNCSCHEMA_SUBMIT_ALLELES_DESCRIPTION)]

SyncSchemaCPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.SYNCSCHEMA_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]


# Main validation class for the SyncSchema module
class SynchronizeSchemaValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	schema_directory: SyncSchemaSchemaDirectory
	nomenclature_server_instance: SyncSchemaNomenclatureServerInstance
	submit_alleles: SyncSchemaSubmitAlleles
	blast_path: gval.BLASTPath
	cpu_cores: SyncSchemaCPUCores

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


# NSStats Annotated Fields
NSStatsStatsMode = Annotated[str,
                      Field(default=None,
                            description=ct.NSSTATS_STATS_MODE_DESCRIPTION),
					  AfterValidator(partial(gval.validate_choice,
					                         valid_options=ct.NSSTATS_MODES,
											 exception_str=ct.NSSTATS_INVALID_MODE))]

NSStatsSpeciesId = Annotated[str,
                      Field(default=None,
                            description=ct.NSSTATS_SPECIES_ID_DESCRIPTION)]

NSStatsSchemaId = Annotated[str,
                     Field(default=None,
                           description=ct.NSSTATS_SCHEMA_ID_DESCRIPTION)]

NSStatsNomenclatureServerInstance = Annotated[str | None,
                               		   Field(default=ct.HOST_NS[ct.DEFAULT_NOMENCLATURE_SERVER],
                                                description=ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION)]


# Main validation class for the NSStats module
class NSStatsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	stats_mode: NSStatsStatsMode
	species_id: NSStatsSpeciesId
	schema_id: NSStatsSchemaId
	nomenclature_server_instance: NSStatsNomenclatureServerInstance

	# Check that it is possible to connect Chewie-NS
	@model_validator(mode="after")
	def check_chewiens_connection(self):
		connection = cr.check_connection(self.nomenclature_server_instance, gct.HEADERS_GET_JSON)
		if not connection:
			sys.exit(f"Could not connect to the Chewie-NS instance at {self.nomenclature_server_instance} (returned status code {connection[1]}).")

		return self


# Define custom exceptions
class DownloadSchemaInvalidDateFormat(Exception):
	"""Raised when the date format is invalid"""
	pass

class DownloadSchemaInvalidDate(Exception):
	"""Raised when the date is not valid due to matching a date prior or later than the schema dates."""
	pass

class DownloadSchemaFailedLociFastaDownload(Exception):
	"""Raised when the download of the loci FASTA files fails."""
	pass

class SchemaLocked(Exception):
	"""Raised when a remote schema is locked."""
	pass

class UploadSchemaCannotContinueUpload(Exception):
	""""""
	pass

class UploadSchemaCannotDetermineSchemaStatus(Exception):
	""""""
	pass

class UploadSchemaCannotCreateSchema(Exception):
	""""""
	pass

class UploadSchemaFailedUpload(Exception):
	""""""
	pass

class SyncSchemaMultipleArgumentValues(Exception):
	""""""
	pass

class SyncSchemaSchemaConfigNoMatch(Exception):
	""""""
	pass

class NSStatsCouldNotRetrieveSpeciesData(Exception):
	""""""
	pass

class NSStatsInvalidSpeciesID(Exception):
	""""""
	pass

class NSStatsCouldNotRetrieveSchemas(Exception):
	""""""
	pass

class NSStatsInvalidSchemaID(Exception):
	""""""
	pass

class NoPermissions(Exception):
	""""""
	pass
