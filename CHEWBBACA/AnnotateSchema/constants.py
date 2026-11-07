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

# Parameter descriptions
ANNOTATESCHEMA_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the schema's directory.")
ANNOTATESCHEMA_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory where the process will store intermediate "
							 "files and save the final TSV file with the loci annotations.")
ANNOTATESCHEMA_LOCI_LIST_DESCRIPTION = ("Path to a file with the list of loci in the schema that the process "
							 "should find annotations for (one per line, full paths or loci IDs).")
ANNOTATESCHEMA_PROTEIN_TABLE_DESCRIPTION = ("Path to the TSV file with coding sequence (CDS) coordinate data, "
							 "`cds_coordinates.tsv`, created by the CreateSchema process.")
ANNOTATESCHEMA_BLAST_SCORE_RATIO_DESCRIPTION = ("BLAST Score Ratio value. The BSR is only used when taxa names are "
							 "provided to the --taxa parameter and local sequences are aligned "
							 "against reference proteomes downloaded from UniProt. Annotations "
							 "are selected based on a BSR >= than the specified value.")
ANNOTATESCHEMA_TRANSLATION_TABLE_DESCRIPTION = ("Genetic code used to translate coding DNA sequences (CDSs). If no value"
							 " is specified, the process tries to get the value stored in the schema "
							 "config file. If the schema does not include a config file, the process "
							 "uses the default translation table (11).")
ANNOTATESCHEMA_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")
ANNOTATESCHEMA_PROTEOME_TAXA_DESCRIPTION = ("List of scientific names for a set of taxa. The process will download "
							 "reference proteomes from UniProt associated to taxa names that contain "
							 "any of the provided terms. The schema representative alleles are aligned "
							 "against the reference proteomes to assign annotations based on high-BSR "
							 "matches.")
ANNOTATESCHEMA_PROTEOME_MATCHES_DESCRIPTION = ("Maximum number of proteome matches to report.")
ANNOTATESCHEMA_UNIPROT_REST_API_DESCRIPTION = ("Search for annotations using UniProt's Proteins API.")
ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_DESCRIPTION = ("Maximum number of matches retrieved from UniProt's Proteins API to report.")
ANNOTATESCHEMA_NO_CLEANUP_DESCRIPTION = ("If provided, intermediate files generated during process execution "
							 "are not removed at the end.")

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
