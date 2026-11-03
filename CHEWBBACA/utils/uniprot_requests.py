#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to perform requests to UniProt.

Code documentation
------------------
"""


import time
import hashlib
import requests
from urllib3.util import Retry
from urllib.request import urlopen
from requests.adapters import HTTPAdapter

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   fasta_operations as fao,
					   sequence_manipulation as sm)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 file_operations as fo,
								 fasta_operations as fao,
								 sequence_manipulation as sm)


def website_availability(url):
	"""Determine if a website is reachable.

	Parameters
	----------
	url : str
		URL to a website.

	Returns
	-------
	reponse : http.client.HTTPResponse or urllib.error.HTTPError
		Response with status code 200 if website is reachable
		or a different status code if the website was not reachable.
	"""
	try:
		response = urlopen(url)
	except Exception as e:
		response = e

	return response


def get_proteomes(proteome_ids, output_dir):
	"""Download reference proteomes from UniProt's FTP.

	Parameters
	----------
	proteomes : list
		List with a sublist per proteome to download.
		Each sublist has the information about a proteome
		that was contained in the README file with the list
		of UniProt's reference proteomes.
	output_dir : str
		Path to the output directory where downloaded
		proteomes will be saved to.

	Returns
	-------
	proteomes_files : list
		Local paths to the downloaded proteomes.
	"""
	print('Downloading reference proteomes...')
	# construct FTP URLs for each proteome
	downloaded = 0
	proteomes_files = []
	for pid in proteome_ids:
		domain = '{0}{1}'.format(pid[3][0].upper(), pid[3][1:])
		proteome_id = '{0}_{1}'.format(pid[0], pid[1])
		proteome_file = '{0}.fasta.gz'.format(proteome_id)
		local_proteome_file = fo.join_paths(output_dir, [proteome_file])
		proteome_url = fo.join_paths(ct.UNIPROT_PROTEOMES_FTP, [domain, pid[0], proteome_file])
		res = fo.download_file(proteome_url, local_proteome_file)
		proteomes_files.append(local_proteome_file)
		downloaded += 1
		print('\r', 'Downloaded {0}/{1}'.format(downloaded, len(proteome_ids)), end='')
		time.sleep(0.1)

	return proteomes_files


def create_requests_session(retries, backoff_factor, status_forcelist):
	""""""
	session = requests.Session()
	# Configure retry strategy
	retry_config = Retry(total=retries, # Total number of retries
	                     backoff_factor=backoff_factor, # Delay between retries {backoff_factor} * (2 ** (number of previous retries))
						 status_forcelist=status_forcelist, # Retry on these HTTP status codes
						 raise_on_status=False, # Return final response instead of raising exception if request fails
						 allowed_methods=["GET"] # Restrict to specific HTTP methods
						)

	adapter = HTTPAdapter(max_retries=retry_config)

	# Mount adapter to both HTTP and HTTPS protocols
	session.mount("http://", adapter)
	session.mount("https://", adapter)

	return session


def send_request(session, url, params, headers):
	""""""
	try:
		response = session.get(url,
		                       params=params,
		                       headers=headers,
		                       # Set timeout to connect to server and to return response
							   timeout=(ct.UNIPROT_RESTAPI_SERVER_TIMEOUT, ct.UNIPROT_RESTAPI_RESPONSE_TIMEOUT)
							)
		response.raise_for_status()
		return response
	except requests.exceptions.RetryError:
		return "Maximum retries exceeded."
	except requests.exceptions.Timeout:
		return "The request timed out."
	except requests.exceptions.HTTPError as err:
		return f"HTTP error occurred {err}"


def get_annotation(input_file, translation_table, maximum_matches):
	"""Retrieve and select best annotation for a locus.

	Parameters
	----------
	input_file : str
		Path to a FASTA file with the DNA sequences for
		the alleles of the locus.
	translation_table : int
		Translation table used to translate DNA sequences.

	Returns
	-------
	gene : str
		Path to a FASTA file with the DNA sequences for
		the alleles of the locus.
	selected_name : str
		Product name selected from the terms retrieved from UniProt.
	selected_url : str
		URL for the page with information about the selected terms.
	"""
	queried = []
	exceptions = []
	accession = None
	# Import locus alleles
	sequences = fao.import_sequences(input_file)
	annotations = []
	# Create session object to perform requests
	session = create_requests_session(ct.UNIPROT_RESTAPI_MAX_RETRIES, ct.UNIPROT_RESTAPI_RETRY_SLEEP, ct.UNIPROT_RESTAPI_STATUS_FORCELIST)
	for seqid, sequence in sequences.items():
		# Translate allele and determine MD5 hash
		protein_sequence = str(sm.translate_sequence(sequence, table_id=translation_table))
		md5_hash = hashlib.md5(protein_sequence.encode("utf-8")).hexdigest()
		# Do not query UniProt if another allele encoding the same protein has already been seen
		if md5_hash not in queried:
			params = {"md5": md5_hash}
			response = send_request(session, ct.UNIPROT_PROTEINS_API, params=params, headers=ct.HEADERS_GET_JSON)
			# Extract results if request was successful
			if response.status_code == 200:
				data = response.json()
				if len(data) > 0:
					# Extract annotation terms
					for i in range(min(maximum_matches, len(data))):
						current_data = data[i]
						# Get UniProt accession number
						accession = current_data.get("accession", "")
						# Get the source of the record (TrEMBL or Swiss-Prot)
						source = current_data.get("info", {}).get("type", "")
						# Get the taxonomic identifier of the organism the protein was identified in
						taxid = str(current_data.get("organism", {}).get("taxonomy", ""))
						# Get teh scientific name of the organism
						taxon = current_data.get("organism", {}).get("names", [{}])[0].get("value", "")
						# Get the protein name/product
						protein_name = current_data.get("protein", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
						# Get the name of the gene encoding the protein
						gene_name = current_data.get("gene", [{}])[0].get("name", {}).get("value", "")
						# Construct the link to the record's page on UniProt
						gene_url = ct.UNIPROT_ENTRY_URL.format(accession)

						annotations.append([accession, source, taxid, taxon, protein_name, gene_name, gene_url])

				queried.append(md5_hash)
			else:
				exceptions.append(response)

		# Stop iterating when one sequence matched a record on UniProt
		if accession:
			return [input_file, annotations]

	return [input_file, []]


def extract_proteome_terms(header_items):
	"""Extract information from the header of a proteome record.

	Extracts the sequence identifier, product name, gene name
	and species name fields from the sequence header of a
	reference proteome from Uniprot.

	Parameters
	----------
	header_items : dict
		Dictionary with the keys and values from a Biopython
		Bio.SeqRecord.SeqRecord object. Created by passing the
		Biopython Bio.SeqRecord.SeqRecord object to the vars
		function.

	Returns
	-------
	seqid : str
		Sequence identifier of the record. Composed
		of the db, UniqueIdentifier and EntryName fields.
	record_product : str
		ProteinName field value in the sequence header.
	record_gene_name : str
		GeneName field (GN) in the sequence header.
	record_species : str
		OrganismName field (OS) in the sequence header.
	"""
	# some tags might be missing
	seqid = header_items.get('id', 'not_found')
	record_description = header_items.get('description', 'not_found')
	if record_description != '':
		if 'OS=' in record_description:
			# get organism name
			record_species = (record_description.split('OS=')[1]).split(' OX=')[0]
			record_product = (record_description.split(seqid+' ')[1]).split(' OS=')[0]
		else:
			record_species = 'not_found'
			record_product = 'not_found'

		if 'GN=' in record_description:
			record_gene_name = (record_description.split('GN=')[1]).split(' PE=')[0]
		else:
			record_gene_name = 'not_found'

	return [seqid, record_product, record_gene_name, record_species]
