#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the AnnotateSchema module.
"""


# Parameter names
ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
ANNOTATESCHEMA_LOCI_LIST_ARGNAME = "loci_list"
ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME = "protein_table"
ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME = "blast_score_ratio"
ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME = "translation_table"
ANNOTATESCHEMA_CPU_CORES_ARGNAME = "cpu_cores"
ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME = "proteome_taxa"
ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME = "proteome_matches"
ANNOTATESCHEMA_BLAST_PATH_ARGNAME = "blast_path"
ANNOTATESCHEMA_UNIPROT_REST_API_ARGNAME = "rest_api"
ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME = "rest_api_matches"
ANNOTATESCHEMA_NO_CLEANUP_ARGNAME = "no_cleanup"

# Argument types
ANNOTATESCHEMA_ARGUMENT_TYPES = {
	"schema_directory": str,
	"output_directory": str,
	"loci_list": str,
	"protein_table": str,
	"blast_score_ratio": float,
	"translation_table": int,
	"cpu_cores": int,
	"proteome_taxa":str,
	"proteome_matches": int,
	"blast_path": str,
	"rest_api": bool,
	"rest_api_matches": int,
	"no_cleanup": bool
}

# UniProt Proteins API URL
UNIPROT_PROTEINS_API = "https://www.ebi.ac.uk/proteins/api/proteins"
# UniProt entry template URL
UNIPROT_ENTRY_URL = "https://www.uniprot.org/uniprotkb/{0}/entry"
# FTP to get UniProt's reference proteomes
UNIPROT_PROTEOMES_FTP = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/")

# Querying UniProt
UNIPROT_REST_API_THREADS = 4
UNIPROT_RESTAPI_MAX_RETRIES = 3
# Time allowed to establish TCP connection
# Set to 3.05 based on the standard practices to set slightly greater than TCP packet retransmission cycles
UNIPROT_RESTAPI_SERVER_TIMEOUT = 3.05
UNIPROT_RESTAPI_RESPONSE_TIMEOUT = 5
UNIPROT_RESTAPI_RETRY_SLEEP = 1
# Define list of status codes after which is allowed to retry
UNIPROT_RESTAPI_STATUS_FORCELIST = [500, 502, 503, 504]

# Column headers/names for the output file
ANNOTATE_SCHEMA_FIRST_COLUMN = "Locus"
ANNOTATE_SCHEMA_PROTEOME_COLUMNS = ["Proteome_ID", "Proteome_Product", "Proteome_Gene_Name", "Proteome_Species", "Proteome_BSR"]
ANNOTATE_SCHEMA_RESTAPI_COLUMNS = ["RestAPI_Accession", "RestAPI_Source", "RestAPI_TaxID", "RestAPI_Taxon", "RestAPI_Protein_Name", "RestAPI_Gene_Name", "RestAPI_URL"]

# Exceptions
MISSING_PROTEIN_TABLE = ("Path for protein table is not valid.")

ANNOTATE_SCHEMA_NO_ANNOTATION_METHOD = ("No annotation method was specified. Please provide either a single "
                                        "or multiple taxon names to '--ptax, --proteome-taxa' or '--rapi, "
										"--rest-api' to query UniProt's Proteins API or both options.")
