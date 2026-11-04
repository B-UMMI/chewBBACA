#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


DOWNLOADSCHEMA_MODULE_NAME = "DownloadShema"

DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
DOWNLOADSCHEMA_SPECIES_ID_ARGNAME = "species_id"
DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME = "schema_id"
DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME = "schema_date"
DOWNLOADSCHEMA_LATEST_VERSION_ARGNAME = "latest_version"
DOWNLOADSCHEMA_BLAST_PATH_ARGNAME = "blast_path"
DOWNLOADSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"

DOWNLOADSCHEMA_ARGUMENT_TYPES = {
	"output_directory": str,
	"species_id": str,
	"schema_id": str,
	"nomenclature_server_instance": str,
	"schema_date": str,
	"latest_version": bool,
	"blast_path": str,
	"cpu_cores": int,
}

UPLOADSCHEMA_MODULE_NAME = "UploadShema"

UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
UPLOADSCHEMA_SPECIES_ID_ARGNAME = "species_id"
UPLOADSCHEMA_SCHEMA_NAME_ARGNAME = "schema_name"
UPLOADSCHEMA_LOCI_PREFIX_ARGNAME = "loci_prefix"
UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME = "description_file"
UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME = "annotations_file"
UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
UPLOADSCHEMA_CONTINUE_UPLOAD_ARGNAME = "continue_upload"
UPLOADSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"

UPLOADSCHEMA_ARGUMENT_TYPES = {
	"schema_directory": str,
	"species_id": str,
	"schema_name": str,
	"loci_prefix": str,
	"description_file": str,
	"annotations_file": str,
	"nomenclature_server_instance": str,
	"continue_upload": bool,
	"cpu_cores": int,
}

SYNCSCHEMA_MODULE_NAME = "SyncShema"

SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
SYNCSCHEMA_SUBMIT_ALLELES_ARGNAME = "submit_alleles"
SYNCSCHEMA_BLAST_PATH_ARGNAME = "blast_path"
SYNCSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"

SYNCSCHEMA_ARGUMENT_TYPES = {
	"schema_directory": str,
	"nomenclature_server_instance": str,
	"submit_alleles": bool,
	"blast_path": str,
	"cpu_cores": int,
}

NSSTATS_MODULE_NAME = "NSStats"

NSSTATS_STATS_MODE_ARGNAME = "stats_mode"
NSSTATS_SPECIES_ID_ARGNAME = "species_id"
NSSTATS_SCHEMA_ID_ARGNAME = "schema_id"
NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"

NSSTATS_ARGUMENT_TYPES = {
	"stats_mode": str,
	"species_id": str,
	"schema_id": str,
	"nomenclature_server_instance": str,
}

# List of Chewie-NS instance identifiers and URLs
HOST_NS = {"main": "https://chewbbaca.online/api/NS/api/",
		   "tutorial": "https://tutorial.chewbbaca.online/api/NS/api/",
		   "local": "http://127.0.0.1:5000/NS/api/"}

DEFAULT_NOMENCLATURE_SERVER = "main"

NSSTATS_MODES = ["species", "schemas"]

LOCI_TO_ADAPT_BASENAME = 'loci_to_adapt.txt'

# Exceptions

NSSTATS_INVALID_MODE = (f"Invalid value passed to '--sm, --stats-mode'. Please select one of the following options: {NSSTATS_MODES}")

# User does not have permissions to upload schemas to Chewie-NS
LOADSCHEMA_NO_PERMISSIONS = ('Current user has no Administrator or Contributor '
							 'permissions.\nNot allowed to upload schemas.')

NS_CANNOT_CONNECT = ('Failed to establish a connection to the Chewie-NS instance at {0}.')
