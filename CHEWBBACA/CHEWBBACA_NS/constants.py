#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the modules taht allow for the integration with Chewie-NS.
"""


# Parameter names for the DownloadSchema module
DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
DOWNLOADSCHEMA_SPECIES_ID_ARGNAME = "species_id"
DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME = "schema_id"
DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME = "schema_date"
DOWNLOADSCHEMA_LATEST_VERSION_ARGNAME = "latest_version"
DOWNLOADSCHEMA_BLAST_PATH_ARGNAME = "blast_path"
DOWNLOADSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types for the DownloadSchema module
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

# Parameter descriptions for the DownloadSchema module
DOWNLOADSCHEMA_OUTPUT_DIRECTORY_DESCRIPTION = ("Output folder to which the schema will be saved.")
DOWNLOADSCHEMA_SPECIES_ID_DESCRIPTION = ("The integer identifier or name of the species that the schema is "
							 "associated to in Chewie-NS.")
DOWNLOADSCHEMA_SCHEMA_ID_DESCRIPTION = ("The URI, integer identifier or name of the schema to download from "
							 "Chewie-NS.")
DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION = ("The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https://"
							 "tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")
DOWNLOADSCHEMA_SCHEMA_DATE_DESCRIPTION = ("Download schema with state from specified date. Must be in the format "
							 "`Y-m-dTH:M:S`.")
DOWNLOADSCHEMA_LATEST_VERSION_DESCRIPTION = ("If the compressed version that is available is not the latest, downloads "
							 "all loci FASTA files and constructs schema locally using the AdaptSchema module.")
DOWNLOADSCHEMA_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads). This value is only used if it is "
							 "necessary to construct the schema locally.")

# Parameter names for the UploadSchema module
UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
UPLOADSCHEMA_SPECIES_ID_ARGNAME = "species_id"
UPLOADSCHEMA_SCHEMA_NAME_ARGNAME = "schema_name"
UPLOADSCHEMA_LOCI_PREFIX_ARGNAME = "loci_prefix"
UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME = "description_file"
UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME = "annotations_file"
UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
UPLOADSCHEMA_CONTINUE_UPLOAD_ARGNAME = "continue_upload"
UPLOADSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types for the UploadSchema module
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

# Parameter descriptions for the UploadSchema module
UPLOADSCHEMA_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the directory of the schema to upload.")
UPLOADSCHEMA_SPECIES_ID_DESCRIPTION = ("The integer identifier or name of the species that the schema will be "
							 "associated to in Chewie-NS.")
UPLOADSCHEMA_SCHEMA_NAME_DESCRIPTION = ("A brief and meaningful name that should help understand the type and "
							 "content of the schema.")
UPLOADSCHEMA_LOCI_PREFIX_DESCRIPTION = ("Prefix included in the name of each locus of the schema.")
UPLOADSCHEMA_DESCRIPTION_FILE_DESCRIPTION = ("Path to a text file with a description about the schema. Markdown "
							 "syntax is supported in order to offer greater customizability of "
							 "the rendered description in the Frontend. Will default to the schema's "
							 "name if the user does not provide a valid path for a file.")
UPLOADSCHEMA_ANNOTATIONS_FILE_DESCRIPTION = ("Path to a TSV file with loci annotations. The first column has loci "
							 "identifiers (w/o .fasta extension), the second has UniProt protein "
							 "names, the third has UniProt gene names, the fourth has UniProt URIs, "
							 "the fifth has user annotations, and the sixth has custom annotations.")
UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION = ("The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https://"
							 "tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")
UPLOADSCHEMA_CONTINUE_UPLOAD_DESCRIPTION = ("Check if the schema upload was interrupted and attempt to continue upload.")
UPLOADSCHEMA_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads). This value is used to accelerate the quality "
							 "control step that checks all alleles in the schema.")

# Parameter names for the SyncSchema module
SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"
SYNCSCHEMA_SUBMIT_ALLELES_ARGNAME = "submit_alleles"
SYNCSCHEMA_BLAST_PATH_ARGNAME = "blast_path"
SYNCSCHEMA_CPU_CORES_ARGNAME = "cpu_cores"

# Argument types for the SyncSchema module
SYNCSCHEMA_ARGUMENT_TYPES = {
	"schema_directory": str,
	"nomenclature_server_instance": str,
	"submit_alleles": bool,
	"blast_path": str,
	"cpu_cores": int,
}

# Parameter descriptions for the SyncSchema module
SYNCSCHEMA_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the directory with the schema to be synced.")
SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION = ("The base URL for the Chewie-NS instance. The default option will get the base "
							 "URL from the schema's URI. It is also possible to specify other options that "
							 "are available in chewBBACA\'s configs, such as: `main` will establish a connection"
							 " to `https://chewbbaca.online/`, `tutorial` to `https://tutorial.chewbbaca.online/` "
							 "and `local` to `http://127.0.0.1:5000/NS/api/` (localhost). Users may also provide "
							 "the IP address to other Chewie-NS instances.")
SYNCSCHEMA_SUBMIT_ALLELES_DESCRIPTION = ("If the process should identify new alleles in the local schema and send them "
							 "to the Chewie-NS instance. (only authorized users can submit new alleles).")
SYNCSCHEMA_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads). This value is only used if the process retrieves "
							 "novel alleles from the remote schema and needs to redetermine the set of "
							 "representative alleles for the local schema.")

# Parameter names for the NSStats module
NSSTATS_STATS_MODE_ARGNAME = "stats_mode"
NSSTATS_SPECIES_ID_ARGNAME = "species_id"
NSSTATS_SCHEMA_ID_ARGNAME = "schema_id"
NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME = "nomenclature_server_instance"

# Argument types for the NSStats module
NSSTATS_ARGUMENT_TYPES = {
	"stats_mode": str,
	"species_id": str,
	"schema_id": str,
	"nomenclature_server_instance": str,
}

# Parameter descriptions for the NSStats module
NSSTATS_STATS_MODE_DESCRIPTION = ("The process can retrieve the list of species (`species` option) in Chewie-NS "
							 "or the list of schemas for a species (`schemas` option).")
NSSTATS_SPECIES_ID_DESCRIPTION = ("The integer identifier of a species in Chewie-NS.")
NSSTATS_SCHEMA_ID_DESCRIPTION = ("The integer identifier of a schema in Chewie-NS.")
NSSTATS_NOMENCLATURE_SERVER_INSTANCE_DESCRIPTION = ("The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https:"
							 "//tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")

# List of Chewie-NS instance identifiers and URLs
HOST_NS = {"main": "https://chewbbaca.online/api/NS/api/",
		   "tutorial": "https://tutorial.chewbbaca.online/api/NS/api/",
		   "local": "http://127.0.0.1:5000/NS/api/"}

# Default instance
DEFAULT_NOMENCLATURE_SERVER = "main"

# Data request modes for the NSStats mode
NSSTATS_MODES = ["species", "schemas"]

LOCI_TO_ADAPT_BASENAME = 'loci_to_adapt.txt'

# Exceptions
NSSTATS_INVALID_MODE = (f"Invalid value passed to '--sm, --stats-mode'. Please select one of the following options: {NSSTATS_MODES}")

# User does not have permissions to upload schemas to Chewie-NS
LOADSCHEMA_NO_PERMISSIONS = ('Current user has no Administrator or Contributor '
							 'permissions.\nNot allowed to upload schemas.')

NS_CANNOT_CONNECT = ('Failed to establish a connection to the Chewie-NS instance at {0}.')
