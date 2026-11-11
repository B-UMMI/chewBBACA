#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module allows authorized users to upload chewBBACA's schemas
to a Chewie-NS instance. The process requests the credentials of
the user trying to upload a schema and determines if the user is
allowed to upload the schema (only Admin and Contributor level
users can upload schemas). The schema config file is read to
validate the argument values used to create the schema. Only
schemas that have been used with a single valid value per parameter
can be uploaded. Invalid or multiple values for a single parameter
can lead to inconsistent results; thus, it is strongly advised to
always perform allele calling with the same set of parameters
and refrain from altering the initial set of parameters values
defined in the schema creation or adaptation processes.

Code documentation
------------------
"""


import os
import json
import time
import requests
import itertools
import multiprocessing

from loguru import logger
from urllib3.exceptions import InsecureRequestWarning

import CHEWBBACA.global_constants as gct
from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
import CHEWBBACA.CHEWBBACA_NS.constants as ct
from CHEWBBACA.utils import (file_operations as fo,
							 chewiens_requests as cr,
							 fasta_operations as fao,
							 sequence_manipulation as sm)


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


def import_annotations(tsv_file):
	"""Read a TSV file to get loci annotations.

	This function imports loci annotations from
	a user-provided TSV file. Users can provide
	annotations for the full set of loci in the
	schema or only for a subset.

	Parameters
	----------
	tsv_file : str
		Path to a TSV file with loci annotations.
		The TSV file has the following three columns:

		- Loci identifiers (schema's filenames stripped of
		  the '.fasta' suffix).
		- User annotation, a term that the user attributes
		  to each locus.
		- Custom annotation, an additional term that the
		  user might want to provide.

		It is not mandatory to provide terms for both
		annotation fields.

	Returns
	-------
	annotations : dict
		Dictionary with loci identifiers as keys and lists
		with user and custom annotations as values.
	"""
	# Read TSV file containing annotations
	lines = fo.read_tabular(tsv_file)
	annotations = {line[0]: line[1:] for line in lines}

	# Substitute blank fields by N/A
	for l, a in annotations.items():
		annotations[l] = [t if t != "" else "N/A" for t in a]

	return annotations


def search_schema(schemas, schema_name):
	"""Determine if any of the schemas in the input list has the provided name.

	Parameters
	----------
	schemas : list of dict
		A list that contains dictionaries.
		Each dictionary has information about a schema.
	schema_name : str
		Name of the schema to look for in the input list
		of schemas.

	Returns
	-------
	list
		A list with the schema URI and schema identifier
		if the schema is found or a list with the boolean
		value of False if it is not found.
	"""
	match = [s for s in schemas if s["name"]["value"] == schema_name]

	if len(match) > 0:
		schema = match[0]
		schema_uri = schema["schemas"]["value"]
		schema_id = schema_uri.split("/")[-1]

		return [schema_uri, schema_id]
	else:
		return [False]


def schema_status(base_url, headers_get, schema_name, species_id, continue_up):
	"""Determine if it is possible to create a schema or resume its upload.

	Parameters
	----------
	base_url : str
		Base URL of the Chewie Nomenclature Server.
	headers_get : dict
		HTTP headers for GET requests.
	schema_name : str
		Name of the schema that the user wants to create
		or resume the upload of.
	species_id : str
		The identifier of the schema's species in the
		Chewie-NS.
	continue_up : bool
		False if the schema is to be uploaded as a new schema
		or True if the schema already exists.

	Returns
	-------
	list
		A list with the upload type ('novel' for new schemas,
		'incomplete' for upload resuming) and the schema
		identifier (NoneType if the schema still does not
		exist).

	Raises
	------
	SystemExit
		- If `continue_up` is True and the schema does not exist.
		- If `continue_up` is False and a schema with name value
		  equal to `schema_name` already exists.
		- If `continue_up` is True and the schema has been fully
		  uploaded.
		- If `continue_up` is True and the current user is not the
		  one that administrates the schema.
	"""
	schema_id = None
	schemas_url, schemas_get = cr.simple_get_request(base_url, headers_get, ["species", species_id, "schemas"])
	schemas_status = schemas_get.status_code
	if schemas_status in [200, 201]:
		species_schemas = schemas_get.json()
		# Determine if there is a schema for current species with provided name
		schema_info = search_schema(species_schemas, schema_name)
		# There is no schema with provided name
		if schema_info[0] is False:
			if continue_up is False:
				upload_type = "novel"
			elif continue_up is True:
				logger.error("Cannot continue uploading a schema that does not exist.")
				raise vm.UploadSchemaCannotContinueUpload
		# There is a schema with provided name
		else:
			if continue_up is False:
				logger.error("Cannot upload schema. A schema with provided name already exists.")
				raise vm.UploadSchemaCannotContinueUpload
			elif continue_up is True:
				schema_url, schema_id = schema_info
				# Determine if schema upload was complete
				schema_get = cr.simple_get_request(base_url, headers_get, ["species", species_id, "schemas", schema_id])[1]
				current_schema = schema_get.json()[0]
				schema_date = current_schema["dateEntered"]["value"]
				if schema_date != "singularity":
					logger.error("Schema finished uploading. Cannot proceed.")
					raise vm.UploadSchemaCannotContinueUpload
				# determine if user was the one that started the upload
				ask_admin = cr.simple_get_request(base_url, headers_get, ["species", species_id, "schemas", schema_id, "administrated"])[1]
				schema_admin = ask_admin.json()
				if schema_admin is False:
					logger.error("Current user is not the user that started schema upload. Exited.")
					raise vm.UploadSchemaCannotContinueUpload
				upload_type = 'incomplete'
	else:
		if continue_up is True:
			logger.error("Cannot continue uploading a schema that does not exist.")
			raise vm.UploadSchemaCannotContinueUpload
		else:
			upload_type = "novel"

	return [upload_type, schema_id]


def quality_control(locus_input):
	"""Translate DNA sequences and report invalid alleles.

	Parameters
	----------
	locus_input : list
		A list with the following elements:

		- Path to the FASTA file with the locus DNA
		  sequences (str).
		- Locus identifier (str).
		- Genetic code that should be used to translate
		  the DNA sequences (int).
		- Minimum sequence length value (int). Alleles
		  with a length value that is lower than this value
		  will be considered invalid.
		- Sequence size threshold value (NoneType). Set
		  to None to avoid removal of alleles that deviate
		  from the sequence length mode based on the current
		  set of sequences, but that did not deviate when they
		  were added.

	Returns
	-------
	list
		A list with the following elements:

		- Path to the FASTA file with the locus DNA
		  sequences (str).
		- Path to a pickled file that contains a dictionary
		  with protein identifiers as keys and protein sequences
		  as values.
		- A list with the identifiers of invalid alleles
		  (list of str).
	"""
	fasta_path, locus_id, output_directories, table_id, min_len, size_threshold = locus_input

	# Import DNA sequences
	dna_seqs = fao.import_sequences(fasta_path)
	# Translate sequences
	_, protein_file, _, invalid = fao.translate_fasta(fasta_path, output_directories[1], table_id)
	prot_seqs = fao.import_sequences(protein_file)

	if size_threshold is not None and len(prot_seqs) > 0:
		# Remove alleles based on length mode and size threshold
		modes, alm, asm, allele_sizes = sm.mode_filter(dna_seqs, size_threshold)
		excluded = set(asm + alm)
		# Remove excluded alleles from dictionaries
		dna_seqs = {seqid: seq
					for seqid, seq in dna_seqs.items()
					if seqid not in excluded}
		prot_seqs = {seqid: seq
					for seqid, seq in prot_seqs.items()
					if seqid not in excluded}

		modes_concat = ":".join(map(str, modes))
		st_percentage = int(size_threshold*100)
		invalid += [[s, gct.ALM_MSG.format(st_percentage, allele_sizes[s], modes_concat)] for s in alm]
		invalid += [[s, gct.ASM_MSG.format(st_percentage, allele_sizes[s], modes_concat)] for s in asm]

	# Rewrite FASTA files with valid sequences only
	if len(invalid) > 0:
		# Create FASTA file with valid DNA sequences
		valid_fasta = os.path.join(output_directories[0], "{0}.fasta".format(locus_id))
		dna_records = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, [(seqid, seq) for seqid, seq in dna_seqs.items()])
		fo.write_lines(dna_records, valid_fasta)
		# Create FASTA file with valid protein sequences
		valid_prot_file = protein_file
		prot_records = fao.fasta_lines(gct.FASTA_RECORD_TEMPLATE, [(seqid, seq) for seqid, seq in prot_seqs.items()])
		fo.write_lines(prot_records, valid_prot_file)
	else:
		valid_fasta = fasta_path
		valid_prot_file = protein_file

	return [valid_fasta, valid_prot_file, invalid]


def create_lengths_files(loci_files, output_directory):
	"""Determine the length of the sequences in a set of FASTA files.

	Parameters
	----------
	loci_files : list of str
		List with paths to FASTA files.
	out_dir : str
		Directory where the output files with the length
		values of the sequences from each input FASTA file
		will be created.

	Returns
	-------
	length_files : list of str
		List with paths to the pickled files that contain
		the sequence length values determined for each input
		FASTA file. Each pickled file contains a dictionary
		with hashes (sequence hashes computed with sha256)
		as keys and sequence lengths as values.
	"""
	length_files = []
	for file in loci_files:
		locus_basename = fo.file_basename(file)
		locus_lengths = {locus_basename: fao.sequence_lengths(file, True)}
		lengths_file = os.path.join(output_directory, "{0}_lengths".format(locus_basename.split(".fasta")[0]))

		fo.pickle_dumper(locus_lengths, lengths_file)
		length_files.append(lengths_file)

	return length_files


def schema_completedness(base_url, species_id, schema_id, headers_get,
						 hashed_files):
	"""Determine a schema's upload status.

	Returns information that indicates which schema data has
	been or hasn't been fully uploaded and inserted into the
	Chewie-NS database.

	Parameters
	----------
	base_url : str
		Base URL of the Chewie Nomenclature server.
	species_id : int
		The identifier of the schema's species in the NS.
	schema_id : int
		The identifier of the schema in the NS.
	headers_get : dict
		HTTP headers for GET requests.
	hashed_files : dict
		Dictionary with file hashes (obtained through BLAKE2
		hashing of file content) as keys and file paths as
		values.

	Returns
	-------
	list
		A list with the following elements:

		- A dictionary with paths to loci FASTA files as keys
		  and boolean values that indicate if loci data has been
		  uploaded and inserted.
		- A dictionary with the same structure as the previous
		  one but that only contains entries for loci that were
		  not fully inserted.
		- A list with the paths to FASTA files of the loci whose
		  set of alleles was not fully uploaded.

		The first and second elements, both dictionaries, have
		lists with four elements as values. The first, second
		and third elements on those lists are Boolean values
		that indicate if the locus has been created, linked to
		the species and linked to the schema, respectively. The
		last element is the BLAKE2 hash of the locus file.
	"""
	# Get info about loci and alleles upload
	schema_loci = cr.simple_get_request(base_url, headers_get, ["species", species_id, "schemas", schema_id, "loci", "data"])[1]
	schema_loci = schema_loci.json()
	if "Not found" in schema_loci:
		logger.error("{0}".format(schema_loci["Not found"]))
		raise vm.UploadSchemaCannotDetermineSchemaStatus
	else:
		schema_loci = schema_loci["hashes"]

	loci_info = {hashed_files[k]: v[1]+[k] for k, v in schema_loci.items()}
	# Determine loci that were not fully inserted
	absent_loci = {hashed_files[k]: v[1]+[k] for k, v in schema_loci.items() if all(v[1]) is not True}
	# Determine loci without alleles
	absent_alleles = [hashed_files[k] for k, v in schema_loci.items() if v[0] is False]

	return [loci_info, absent_loci, absent_alleles]


def create_schema(base_url, headers_post, species_id, params):
	"""Create a schema in Chewie-NS.

	Parameters
	----------
	base_url : str
		Base URL of the Nomenclature server.
	headers_post : dict
		HTTP headers for POST requests.
	species_id : int
		The identifier of the schema's species in the
		Chewie-NS.
	params : dict
		A dictionary with the following key/values pairs:

		- bsr (str): BLAST Score Ratio value used to create the
		  schema and perform allele calling.
		- prodigal_training_file (str): BLAKE2 hash of the Prodigal
		  training file content.
		- translation_table (str): genetic code used to predict
		  and translate coding sequences.
		- minimum_locus_length (str): minimum sequence length,
		  sequences with a length value lower than this value
		  are not included in the schema.
		- chewBBACA_version (str): version of the chewBBACA suite
		  used to create the schema and perform allele calling.
		- size_threshold (str): sequence size variation percentage
		  threshold, new alleles cannot have a length value that
		  deviates +/- than this value in relation to the locus's
		  representative sequence.
		- word_size (str): word/k value used to cluster protein
		  sequences during schema creation and allele calling.
		- cluster_sim (str): proportion of k-mers that two proteins
		  need to have in common to be clustered together.
		- representative_filter (str): proportion of k-mers that a
		  clustered protein has to have in common with the representative
		  protein of the cluster to be considered the same gene.
		- intraCluster_filter (str): proportion of k-mers that clustered
		  proteins have to have in common to be considered of the same
		  gene.
		- name (str): name of the new schema. Must be unique among
		  the schemas' of the species.
		- SchemaDescription (str): BLAKE2 of the contents of the file
		  with the schema description.
		- schema_hashes (list of str): list that contains the BLAKE2
		  hashes of the schema files.

	Returns
	-------
	list
		A list with the following elements: the schema URI
		(str) and the schema identifier (str) in the Chewie-NS.
	"""
	schema_post = cr.simple_post_request(base_url, headers_post, ["species", species_id, "schemas"], data=json.dumps(params))[1]
	schema_status = schema_post.status_code

	# check status code
	if schema_status not in [200, 201]:
		logger.error(f"Could not create new schema in the NS. Returned {schema_status} status code ({schema_post.text}).")
		raise vm.UploadSchemaCannotCreateSchema

	schema_url = schema_post.json()["url"]
	schema_id = schema_url.split("/")[-1]

	return [schema_url, schema_id]


def create_loci_file(schema_files, output_directory, annotations, species_id, schema_id, loci_prefix, absent_loci=None):
	"""Create a file with the essential data to insert loci in Chewie-NS.

	Parameters
	----------
	schema_files : list of str
		List with paths to schema's FASTA files.
	annotations : dict
		Dictionary with loci identifiers as keys and lists
		with annotation terms as values.
	species_id : int
		The identifier of the schema's species in the NS.
	schema_id : int
		The identifier of the schema in the NS.
	loci_prefix : str
		Prefix to include in loci identifiers.
	absent_loci : dict, optional
		Loci that have not been created or are incomplete.
		Default value of None collects and saves data
		from all schema's loci.

	Returns
	-------
	loci_file : str
		Path to the file with the necessary data to insert loci
		and associate with species and schema.
	"""
	loci_file = os.path.join(output_directory, "{0}_{1}_loci".format(species_id, schema_id))
	loci_data = [loci_prefix, []]
	for file in schema_files:
		if absent_loci is None or file in absent_loci:

			locus = os.path.basename(file).split(".fasta")[0]
			locus_hash = fo.hash_file(file, "blake2b")

			locus_annotations = annotations[locus]

			loci_data[1].append([locus, locus_hash]+locus_annotations)

	if len(loci_data[1]) > 0:
		fo.pickle_dumper(loci_data, loci_file)

	return loci_file


def upload_loci_data(loci_file, base_url, species_id, schema_id, headers_post, headers_get, hashed_files):
	"""Send data to insert and associate loci to a species and schema.

	Parameters
	----------
	loci_file : str
		Path to the file with the data to upload to the
		Chewie-NS.
	base_url : str
		Base URL of the Nomenclature server.
	species_id : int
		The identifier of the schema's species in the NS.
	schema_id : int
		The identifier of the schema in the NS.
	headers_post : dict
		HTTP headers for POST requests.
	hashed_files : dict
		Dictionary with file hashes as keys and file paths
		as values.

	Returns
	-------
	response_data : dict
		Dictionary with paths to loci FASTA files as keys and
		lists with four elements as values. The first, second
		and third elements on those lists are Boolean values
		that indicate if the locus has been created, linked to
		the species and linked to the schema, respectively. The
		last element is the BLAKE2 hash of the locus file.

	Raises
	------
	SystemExit

		- If the data sent to the NS does not match the data that
		  was determined to be in the schema.
		- If the response returned from the NS indicates that it
		  was not possible to create and associate any/some loci.
	"""
	# Compress file with loci data to reduce upload size
	loci_zip_file = "{0}.zip".format(loci_file)
	fo.file_zipper(loci_file, loci_zip_file)

	zip_url = cr.make_url(base_url, "species", species_id, "schemas", schema_id, "loci", "data")

	# Upload file as multipart-encoded file
	filename = "{0}_{1}_loci.zip".format(species_id, schema_id)
	response = cr.upload_file(loci_zip_file, filename, zip_url, headers_post, False)

	response = response.json()
	start_nr_loci = int(response["nr_loci"])
	start_sp_loci = int(response["sp_loci"])
	start_sc_loci = int(response["sc_loci"])

	time_limit = 2100
	current_time = 0
	status = "Inserting"
	while status != "Complete" and (current_time < time_limit):
		insertion_status = cr.simple_get_request(base_url, headers_get, ["species", species_id, "schemas", schema_id, "loci", "data"])[1]
		insertion_status = insertion_status.json()

		if insertion_status["status"] == "complete":
			status = "Complete"
			response_data = insertion_status["hashes"]

		nr_loci = int(insertion_status["nr_loci"])
		sp_loci = int(insertion_status["sp_loci"])
		sc_loci = int(insertion_status["sc_loci"])

		inserted_loci = nr_loci - start_nr_loci
		linked_sp = sp_loci - start_sp_loci
		linked_sc = sc_loci - start_sc_loci

		logger.info(f"Inserted {inserted_loci} loci; Linked {linked_sp} to species; Linked {linked_sc} to schema.")
		time.sleep(2)
		current_time += 2

	response_data = {hashed_files[k]: v+[k] for k, v in response_data.items()}

	return response_data


def create_alleles_files(schema_files, output_directory, loci_responses, invalid_alleles,
						 species_name, base_url, species_id,
						 schema_id, user_id):
	"""Create files with data to insert alleles in Chewie-NS.

	Parameters
	----------
	schema_files : list
		List with paths to schema's FASTA files.
	loci_responses : dict
		Dictionary with files paths as keys and
	invalid_alleles : list
		List with the identifiers of the alleles
		that were determined to be invalid.
	species_name : str
		Name of the species.
	base_url : str
		Base URL of the Nomenclature server.
	species_id : int
		The identifier of the schema's species in the NS.
	schema_id : int
		The identifier of the schema in the NS.
	user_id : int
		Current user identifier.

	Returns
	-------
	list of list
		List with the following elements:

		- List with paths to files that contain the data that
		  is necessary to upload to insert all alleles.
		- List with the loci identifiers in the NS.
		- List with file hashes for loci files.
		- List with loci basenames.
	"""
	loci_ids = []
	loci_names = []
	loci_hashes = []
	alleles_files = []
	user_uri = f"{base_url}users/{user_id}"
	for file in schema_files:
		locus_uri = loci_responses[file][1][0]

		alleles_sequences = [str(rec.seq)
							 for rec in fao.sequence_generator(file)
							 if rec.id not in invalid_alleles]

		post_inputs = [locus_uri, species_name,
					   user_uri, tuple(alleles_sequences)]

		locus_id = locus_uri.split("/")[-1]
		loci_ids.append(locus_id)

		locus_hash = loci_responses[file][-1]
		loci_hashes.append(locus_hash)

		locus_name = os.path.basename(file).split(".fasta")[0]
		loci_names.append(locus_name)

		alleles_file = os.path.join(output_directory, f"{species_id}_{schema_id}_{locus_id}")
		alleles_files.append(alleles_file)

		fo.pickle_dumper(post_inputs, alleles_file)

	return [alleles_files, loci_ids, loci_hashes, loci_names]


def upload_alleles_data(alleles_data, length_files, base_url, headers_post, headers_post_bytes, species_id, schema_id):
	"""Upload files with data to insert alleles.

	Parameters
	----------
	alleles_data : list
		List with tuples, one per locus, that contain the path
		to the ZIP archive with the data to insert alleles,
		the identifier of the locus, the locus file hash and
		the basename of the locus file.
	length_files : list
		List with paths to the pickled files that contain a
		dictionary with sequences hashes as keys and sequence
		length as values.
	base_url : str
		Base URL of the Nomenclature server.
	headers_post : dict
		HTTP headers for POST requests that accept JSON
		formatted data.
	headers_post_bytes : dict
		HTTP headers for POST requests that support file
		upload.
	species_id : int
		The identifier of the schema's species in the NS.
	schema_id : int
		The identifier of the schema in the NS.

	Returns
	-------
	failed : list of str
		List with the identifiers of the loci whose alleles
		data could not be fully uploaded.
	"""
	failed = []
	uploaded = 0
	for i, a in enumerate(alleles_data):
		current_locus = a[1]
		locus_hash = a[2]

		# get length of alleles from current locus
		current_len = length_files[i]
		data = fo.pickle_loader(current_len)
		data = {current_locus: data[next(iter(data))]}
		data = {"content": data, "locus_hash": locus_hash}

		# send data to the NS
		send_url = cr.make_url(base_url, "species", species_id, "schemas", schema_id, "loci", current_locus, "lengths")

		lengths_res = cr.simple_post_request(send_url, headers_post, data=json.dumps(data))
		length_status = lengths_res[1].status_code

		# get path to ZIP archive with data to insert alleles
		current_zip = a[0]

		# send data to insert alleles in the NS
		zip_url = cr.make_url(base_url, "species", species_id, "schemas", schema_id, "loci", current_locus, "data")

		zip_res = cr.upload_file(current_zip, locus_hash, zip_url, headers_post_bytes, False)
		zip_status = zip_res.status_code

		# determine if upload was successful
		if length_status not in [200, 201] or zip_status not in [200, 201]:
			failed.append(current_locus)
		elif length_status in [200, 201] and zip_status in [200, 201]:
			uploaded += 1
			logger.info("Sent data for alleles of {uploaded} loci.")

	return failed


def main(config):
	"""Upload a schema to a Chewie-NS module.

	Parameters
	----------
	schema_directory : str
		Path to the directory of the schema to upload.
	species_id : int
		The integer identifier or name of the species that
		the schema will be associated to in Chewie-NS.
	schema_name : str
		A brief and meaningful name that should help
		understand the type and content of the schema.
	loci_prefix : str
		A short prefix included in the name of each locus.
	description_file : str
		Path to a text file with a description about the schema. Markdown
		syntax is supported in order to allow greater customizability of
		the rendered description in the Frontend. Will default to the
		schema's name if the user does not provide a valid path for a
		file.
	annotations : str
		Path to a TSV file with loci annotations. The first column has
		loci identifiers (w/o .fasta extension), the second has UniProt protein names,
		the third has UniProt gene names, the fourth has UniProt URIs, the fifth has
		user annotations, and the sixth has custom annotations.
	cpu_cores : int
		Number of CPU cores that will be used in the pre-processing steps.
	nomenclature_server : str
		The base URL for the Chewie-NS instance.
	continue_up : bool
		If the process should check if the schema upload was interrupted
		and try to resume it.
	"""
	# Verify user
	if config.nomenclature_server_instance != ct.HOST_NS["tutorial"]:
		token = cr.capture_login_credentials(config.nomenclature_server_instance)
	else:
		token = ""

	# GET request headers
	headers_get = gct.HEADERS_GET_JSON
	headers_get["Authorization"] = token

	# Determine current user ID and Role
	if "tutorial" not in config.nomenclature_server_instance:
		user_id, user_role, user_auth = cr.user_info(config.nomenclature_server_instance, headers_get)
	else:
		user_id, user_role, user_auth = ["", "", True]

	logger.info(f"User id: {user_id}")
	logger.info(f"User role: {user_role}")
	logger.info(f"Authorized: {user_auth}")

	# Only Admin or Contributor type users can upload schemas
	if not user_auth:
		logger.error(ct.LOADSCHEMA_NO_PERMISSIONS)
		raise vm.NoPermissions

	# POST requests headers
	headers_post = gct.HEADERS_POST_JSON
	headers_post["Authorization"] = token
	headers_post["user_id"] = user_id
	# POST headers to send binary data
	headers_post_bytes = gct.HEADERS_POST
	headers_post_bytes["Authorization"] = token
	headers_post_bytes["user_id"] = user_id

	logger.info(f"Local schema: {config.schema_directory[0]}")

	# Get schema files from genes list file
	genes_list = fo.join_paths(config.schema_directory[0], [gct.GENE_LIST_BASENAME])
	genes = fo.pickle_loader(genes_list)
	fasta_paths = [os.path.join(config.schema_directory[0], file) for file in genes]
	fasta_paths.sort()

	# Total number of loci
	total_loci = len(fasta_paths)
	# Total number of alelles
	total_alleles = sum(list(map(fao.count_sequences, fasta_paths)))

	logger.info(f"Number of loci: {total_loci}")
	logger.info(f"Number of alleles: {total_alleles}\n")

	# Add name of the schema to the schema config
	params = config.schema_directory[1]
	params["name"] = config.schema_name
	# Print list of schema config
	params_lines = [f"{k}: {str(v)}" for k, v in params.items()]
	params_text = "\n".join(params_lines)
	logger.info(params_text)

	# Determine schema status
	upload_type = schema_status(config.nomenclature_server_instance, headers_get,
	                            config.schema_name, config.species_id[0], config.continue_up)

	if upload_type[0] == "novel":
		logger.info("New schema name: {0}".format(config.schema_name))
	else:
		schema_url = "{0}species/{1}/schemas/{2}".format(config.nomenclature_server_instance,
														 config.species_id[0],
														 upload_type[1])
		schema_id = schema_url.split("/")[-1]
		logger.info(f"Schema exists and upload was not completed ({config.schema_name}, id={schema_id}).")

	# Create folder to store temporary files
	upload_temp_dir = os.path.join(config.schema_directory[0], "upload_temp")
	fo.create_directory(upload_temp_dir)

	# Get schema description
	if config.continue_up is False:
		if description_file is not None and os.path.isfile(description_file) is True:
			# Determine file hash
			description_hash = fo.hash_file(description_file, "blake2b")
			logger.info("Schema description: {0}".format(description_file))
		else:
			logger.info("Could not get a description from a file. Will use the name of the schema as description.")
			description_file = fo.join_paths(upload_temp_dir, ["schema_description.txt"])
			with open(description_file, "w") as sd:
				sd.write(config.schema_name)
			description_hash = fo.hash_file(description_file, "blake2b")

		params["SchemaDescription"] = description_hash

	# Hash schema files to get unique identifiers based on content
	hashed_files = {fo.hash_file(file, "blake2b"): file for file in fasta_paths}
	logger.info("Determining data to upload...")
	absent_loci = fasta_paths
	if upload_type[0] == "incomplete":
		loci_info, absent_loci, fasta_paths = schema_completedness(config.nomenclature_server_instance,
		                                                           config.species_id[0],
																   upload_type[1],
																   headers_get,
																   hashed_files)

	logger.info(f"Loci to create and associate with species and schema: {len(absent_loci)}")
	logger.info(f"Loci without the full set of alleles: {len(fasta_paths)}")

	# Create temporary directory to store FASTA files containing only the DNA dna protein sequences for the valid alleles
	valid_fasta_dir = os.path.join(upload_temp_dir, "valid_fasta")
	fo.create_directory(valid_fasta_dir)
	translated_dir = os.path.join(upload_temp_dir, "translated_fasta")
	fo.create_directory(translated_dir)

	# Create inputs for QC step
	inputs = [(file,
			   file.split("/")[-1].split(".fasta")[0],
			   [valid_fasta_dir, translated_dir],
			   int(params["translation_table"]),
			   0,
			   None) for file in fasta_paths]

	# Validate schema data and create files with translated sequences
	logger.info("Translating sequences based on schema config...")
	qc_results = []
	genes_pools = multiprocessing.Pool(processes=config.cpu_cores)
	rawr = genes_pools.map_async(quality_control, inputs,
								 callback=qc_results.extend)
	rawr.wait()

	# Get invalid alleles
	invalid_alleles = [r[2] for r in qc_results]
	invalid_alleles = list(itertools.chain.from_iterable(invalid_alleles))
	invalid_identifiers = set([r[0] for r in invalid_alleles])

	logger.info(f"Found a total of {len(invalid_identifiers)} invalid alleles.")

	# List translated sequences files
	dna_files = [r[0] for r in qc_results]
	prot_files = [r[1] for r in qc_results]

	# Determine loci missing annotations
	miss_annotation = [os.path.basename(pf.split("_protein")[0]) for pf in prot_files
					   if pf.split("_protein")[0] + ".fasta" in absent_loci]

	logger.info(f"Missing annotations for {len(miss_annotation)} loci")

	if len(miss_annotation) > 0:
		# Get locus annotations
		loci_annotations = {locus_id: [] for locus_id in miss_annotation}
		# Import user-provided annotations
		if os.path.isfile(str(config.annotations_file)) is True:
			user_annotations = import_annotations(config.annotations_file)
		else:
			if config.annotations_file is None:
				logger.info("User did not provide a path to a file with annotations.")
			else:
				# Value passed by user is not a valid file path
				logger.info("Could not find the provided file with annotations.")
			user_annotations = []

		provided = 0
		for l in loci_annotations:
			if l in user_annotations:
				loci_annotations[l].extend(user_annotations[l])
				provided += 1
			else:
				loci_annotations[l].extend(["N/A", "N/A", "N/A", "N/A", "N/A"])
		logger.info(f"User provided valid annotations for {provided} loci.")

	# Convert parameters to string type because Chewie-NS expects strings
	params = {k: str(v) for k, v in params.items()}
	params["schema_hashes"] = list(hashed_files.keys())

	# Start sending data
	# Build the new schema URL and POST to NS
	if config.continue_up is False:
		schema_url, schema_id = create_schema(config.nomenclature_server_instance,
											  headers_post,
											  config.species_id[0],
											  params)
		# Send file with description
		description_uri = cr.make_url(config.nomenclature_server_instance, "species",
									  config.species_id[0], "schemas", schema_id,
									  "description")

		desc_res = cr.upload_file(description_file, description_hash,
								  description_uri, headers_post_bytes,
								  False)

		logger.info(f"Created schema with name {config.schema_name} (id={schema_id}).")
	else:
		logger.info(f"Continuing upload of schema with name {config.schema_name} (id={schema_id}).")

	# Start creating new loci and adding/linking alleles
	if upload_type[0] == "incomplete":
		if len(absent_loci) == 0:
			response_data = {k: [v[0], True, True, True, v[3]]
							 for k, v in loci_info.items()}
			logger.info("All loci were inserted in a previous process.")
		else:
			loci_data_dir = os.path.join(upload_temp_dir, "loci_data")
			logger.info("Collecting loci data...")
			loci_file = create_loci_file(dna_files, loci_data_dir, loci_annotations,
										 config.species_id[0],
										 schema_id, config.loci_prefix,
										 absent_loci)
			logger.info("Sending data to chewie-NS...")
			absent_data = upload_loci_data(loci_file, config.nomenclature_server_instance,
										   config.species_id[0], schema_id,
										   headers_post_bytes, headers_get,
										   hashed_files)
			response_data = {k: [v[0], True, True, True, v[3]]
							 for k, v in loci_info.items() if k not in absent_data}
			response_data = {**absent_data, **response_data}
			logger.info(f"The NS completed the insertion of {len(absent_loci)} loci.")
	else:
		# Create directory to store loci data files
		loci_data_dir = os.path.join(upload_temp_dir, 'loci_data')
		fo.create_directory(loci_data_dir)
		logger.info("Collecting loci data...")
		loci_file = create_loci_file(dna_files, loci_data_dir, loci_annotations,
									 config.species_id[0],
									 schema_id, config.loci_prefix)
		logger.info("Sending data to chewie-NS...")
		response_data = upload_loci_data(loci_file, config.nomenclature_server_instance,
										 config.species_id[0], schema_id,
										 headers_post_bytes, headers_get,
										 hashed_files)

		logger.info("The NS completed the insertion of {len(response_data)} loci.")

	# Create files with info for posting alleles
	allele_data_dir = os.path.join(upload_temp_dir, "allele_data")
	fo.create_directory(allele_data_dir)
	logger.info("Collecting alleles data...")
	(alleles_files, loci_ids, loci_hashes, loci_names) = create_alleles_files(dna_files, allele_data_dir, response_data,
										invalid_identifiers, config.species_id[1],
										config.nomenclature_server_instance, config.species_id[0],
										schema_id, user_id)
	# Determine length of all alleles per locus
	length_files = create_lengths_files(dna_files, allele_data_dir)

	# Compress all files to reduce upload size
	logger.info("Compressing files with alleles data...")
	zipped_files = ["{0}.zip".format(file) for file in alleles_files]
	list(map(fo.file_zipper, alleles_files, zipped_files))
	alleles_data = list(zip(zipped_files, loci_ids, loci_hashes, loci_names))

	logger.info("Sending alleles data to chewie-NS...")
	# send POST with file contents and process each file in the NS
	failed = upload_alleles_data(alleles_data, length_files,
								 config.nomenclature_server_instance, headers_post,
								 headers_post_bytes, config.species_id[0],
								 schema_id)

	if len(failed) > 0:
		logger.error("Could not upload data for alleles of following loci:\n"
				 	"{0}\n Please retry with the '--continue_up' option and "
				 	"contact the NS Admin if the problem persists.".format(','.join(failed)))
		raise vm.UploadSchemaFailedUpload
	else:
		# Upload training file to Chewie-NS
		logger.info("Uploading Prodigal training file...")
		ptf_url = cr.make_url(config.nomenclature_server_instance, "species", config.species_id[0],
							  "schemas", schema_id, "ptf")
		ptf_res = cr.upload_file(params["prodigal_training_file"][0], params["prodigal_training_file"][1],
								 ptf_url, headers_post_bytes,
								 False)
		logger.info(list(ptf_res.json().values())[0])
		logger.info("Chewie-NS has received the data and will insert the alleles into the database.")
		logger.info("Schema will be available for download as soon as the process has completed.")
		logger.info("Schema information will also be available on the chewie-NS website.")

	# Delete all intermediate files
	logger.info("Removing intermediate files...")
	fo.delete_directory(upload_temp_dir)
