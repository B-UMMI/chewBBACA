#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module allows authorized users to upload chewBBACA's schemas
to the Chewie-NS.

The process for schema upload has four stages:

    - User Permissions: Determines if the current user has permission
      to upload schemas. Only Admin or Contributor level users can
      upload schemas to the Chewie-NS.

    - Parameters Validation: Validation of the set of parameters associated
      with the schema. Only schemas that have been used with a single valid
      value per parameter can be uploaded. Invalid or multiple values
      for a single parameter can lead to inconsistent results; thus,
      it is strongly advised to always perform allele calling with
      the same set of parameters and refrain from altering the initial
      set of parameters values defined in the schema creation or
      adaptation processes.

    - Schema Pre-processing: Applies quality control measures to identify
      and exclude invalid alleles. Searches for annotations on UniProt
      and imports annotations provided by users.

    - Schema Upload: Collects essential data and sends it to the Chewie-NS
      for schema creation and data insertion. The process finishes when all
      the necessary data has been uploaded. The Chewie-NS automatically
      detects that all data has been received and finishes data insertion.

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_files`` : Path to the directory of the schema to upload.

    - e.g.: ``/home/user/schemas/ypestis_schema``

- ``-sp``, ``species_id`` : The integer identifier or name of the species that
  the schema will be associated to in the Chewie-NS.

    - e.g.: ``1`` or ``'Yersinia pestis'``

- ``-sn``, ``schema_name`` : A brief and meaningful name that should help
  understand the type and content of the schema.

    - e.g.: ``ypestis_cgMLST`` or ``ypestis cgMLST``

- ``-lp``, ``loci_prefix`` : Prefix included in the name of each locus of the
  schema.

    - e.g.: ``ypestis``

- ``--df``, ``description_file`` : Path to a text file with a description
  about the schema. Markdown syntax is supported in order to allow greater
  customizability of the rendered description in the Frontend. Will default
  to the schema's name if the user does not provide a valid path for a
  file (default=None).

    - e.g.: ``/home/user/schemas/ypestis_description``

- ``--a``, ``annotations`` : Path to a TSV file with loci annotations. The
  first column has loci identifiers (w/o .fasta extension), the second has
  user annotations and the third has custom annotations (default=None).

    - e.g.: ``/home/user/schemas/ypestis_annotations``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores that will be used in the
  Schema Pre-processing step (default=1).

    - e.g.: ``4``

- ``--thr``, ``threads`` : Number of threads to use to search for annotations
  on UniProt (default=20).

    - e.g.: ``20``

- ``--ns_url``, ``nomenclature_server_url`` : The base URL for the Nomenclature
  Server. The default value, "main", will establish a connection to
  "https://chewbbaca.online/", "tutorial" to "https://tutorial.chewbbaca.online/"
  and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). Users may also
  provide the IP address to other Chewie-NS instances.

    - e.g.: ``http://127.0.0.1:5000/NS/api/`` (localhost)

- ``--continue_up`` : If the process should check if the schema upload was
  interrupted and try to resume it. ``True`` if provided, ``False`` otherwise.

Code documentation
------------------
"""


import os
import sys
import csv
import time
import json
import pickle
import argparse
import requests
import itertools
import datetime as dt
import multiprocessing
import concurrent.futures

from Bio import SeqIO
from SPARQLWrapper import SPARQLWrapper, JSON
from urllib3.exceptions import InsecureRequestWarning

try:
    from utils import constants as cnst
    from utils import auxiliary_functions as aux
    from utils import parameters_validation as pv
except:
    from CHEWBBACA.utils import constants as cnst
    from CHEWBBACA.utils import auxiliary_functions as aux
    from CHEWBBACA.utils import parameters_validation as pv


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


uniprot_sparql = SPARQLWrapper(cnst.UNIPROT_SPARQL)


def import_annotations(tsv_file):
    """ Reads a TSV file to get loci annotations.

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

    # read file and fill missing values with N/A
    with open(tsv_file, 'r') as tf:
        lines = csv.reader(tf, delimiter='\t')
        lines = [l + ['N/A']*(3-len(l)) for l in lines]

    annotations = {}
    for l in lines:
        locus = l[0]
        # do not keep blank values
        a = ['N/A' if i.strip() == '' else i for i in l[1:]]

        # only keep cases that have at least one non-NA field
        if a[0] != 'N/A' or a[1] != 'N/A':
            annotations[locus] = a

    return annotations


def search_schema(schemas, schema_name):
    """ Determines if any of the schemas in the input list
        has the provided name.

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

    match = [s for s in schemas if s['name']['value'] == schema_name]

    if len(match) > 0:
        schema = match[0]
        schema_uri = schema['schemas']['value']
        schema_id = schema_uri.split('/')[-1]

        return [schema_uri, schema_id]
    else:
        return [False]


def schema_status(base_url, headers_get, schema_name, species_id, continue_up):
    """ Determines if it is possible to create a new schema or
        resume the upload of an incomplete schema.

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
    schemas_get = aux.simple_get_request(base_url, headers_get,
                                         ['species', species_id, 'schemas'])
    schemas_status = schemas_get.status_code
    if schemas_status in [200, 201]:
        species_schemas = schemas_get.json()
        # determine if there is a schema for current
        # species with provided name
        schema_info = search_schema(species_schemas, schema_name)

        # there is no schema with provided name
        if schema_info[0] is False:
            if continue_up is False:
                upload_type = 'novel'
            elif continue_up is True:
                sys.exit('Cannot continue uploading a schema that '
                         'does not exist.')
        # there is a schema with provided name
        else:
            if continue_up is False:
                sys.exit('Cannot upload schema. A schema with '
                         'provided name already exists.')
            elif continue_up is True:
                schema_url, schema_id = schema_info
                # determine if schema upload was complete
                schema_get = aux.simple_get_request(base_url, headers_get,
                                                    ['species', species_id,
                                                     'schemas', schema_id])
                current_schema = schema_get.json()[0]
                schema_date = current_schema['dateEntered']['value']
                if schema_date != 'singularity':
                    sys.exit('Schema finished uploading. Cannot proceed.')
                # determine if user was the one that started the upload
                ask_admin = aux.simple_get_request(base_url, headers_get,
                                                   ['species', species_id,
                                                    'schemas', schema_id,
                                                    'administrated'])
                schema_admin = ask_admin.json()
                if schema_admin is False:
                    sys.exit('Current user is not the user that '
                             'started schema upload. Exited.')
                upload_type = 'incomplete'
    else:
        if continue_up is True:
            sys.exit('Cannot continue uploading a schema that '
                     'does not exist.')
        else:
            upload_type = 'novel'

    return [upload_type, schema_id]


def create_uniprot_queries(file, max_queries=cnst.MAX_QUERIES):
    """ Creates SPARQL queries to search for protein annotations
        on UniProt.

        Parameters
        ----------
        file : str
            Path to a pickled file that contains a dictionary
            with protein identifiers as keys and protein
            sequences as values.
        max_queries : int
            Maximum number of queries that will be created.
            Only the first `max_queries` unique proteins will
            be used to create queries.

        Returns
        -------
        queries_file : str
            Path to the pickled file that contains a list
            with the SPARQL queries that were created.
    """

    # import dictionary with protein sequences as values
    protein_seqs = aux.pickle_loader(file)

    # create queries based on unique sequences only
    unique_prots = set(list(protein_seqs.values()))
    selected = unique_prots if len(unique_prots) <= max_queries \
                               else list(unique_prots)[0:max_queries]

    queries = [aux.uniprot_query(prot) for prot in selected]
    # save SPARQL queries with pickle
    queries_file = '{0}_up'.format(file)
    aux.pickle_dumper(queries_file, queries)

    return queries_file


def get_annotation(queries_file, max_queries=cnst.MAX_QUERIES):
    """ Queries the UniProt SPARQL endpoint to retrieve
        protein annotations.

        This function uses the SPARQL queries imported from
        a pickled file to search for protein annotations
        through UniProt's SPARQL endpoint. It searches until
        it retrieves terms that are considered informative or
        until it reaches the maximum number of searches that is
        allowed.

        Parameters
        ----------
        queries_file : str
            Path to a pickled file that contains a list with
            SPARQL queries that can be used to search for protein
            annotations.

        Returns
        -------
        annotation_info : list of str
            A list with the following elements:

            - Locus identifier.
            - UniProt annotation name.
            - UniProt annotation label.
            - UniProt annotation URL.
    """

    locus = queries_file.split('/')[-1].split('_prots_up')[0]

    # load queries for locus
    queries = aux.pickle_loader(queries_file)

    uniprot_sparql.setReturnFormat(JSON)
    uniprot_sparql.setTimeout(10)

    prev_url = 'N/A'
    prev_name = 'N/A'
    prev_label = 'N/A'
    uninformative = ['Uncharacterized protein',
                     'hypothetical protein',
                     'DUF']

    tries = 0
    # define maximum number of tries
    max_tries = max_queries
    found = False
    while found is False:

        uniprot_sparql.setQuery(queries[tries])

        try:
            result = uniprot_sparql.query().convert()
            name, url, label = aux.select_name(result)

            if name != '':
                if prev_name == 'N/A':
                    prev_name = name
                    prev_label = label
                    prev_url = url
                # keep latest term if previous term was uninformative
                elif any([n in prev_name for n in uninformative]) is True \
                        and any([n in name for n in uninformative]) is False:
                    prev_name = name
                    prev_label = label
                    prev_url = url

            if prev_name not in uninformative:
                found = True

        except Exception:
            pass

        tries += 1
        if tries == max_tries or tries == len(queries):
            found = True

    annotation_info = [locus, prev_name, prev_label, prev_url]

    return annotation_info


def quality_control(locus_input):
    """ Translates DNA sequences based on a set of quality criteria
        and reports invalid alleles.

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

    res = aux.get_seqs_dicts(*locus_input, max_proteins=cnst.MAX_QUERIES)

    prots_file = '{0}_prots'.format(locus_input[0].split('.fasta')[0])
    aux.pickle_dumper(prots_file, res[1])

    if len(res[2]) > 0:
        print('  Found {0} invalid alleles for '
              'locus {1}.'.format(len(res[2]), locus_input[1]))

    return [locus_input[0], prots_file, res[2]]


def create_lengths_files(loci_files, out_dir):
    """ Determines the length of the sequences in a set of
        FASTA files.

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
        locus = os.path.basename(file).split('.fasta')[0]
        locus_lengths = aux.sequences_lengths(file)
        lengths_file = os.path.join(out_dir,
                                    '{0}_lengths'.format(locus))

        aux.pickle_dumper(lengths_file, locus_lengths)
        length_files.append(lengths_file)

    return length_files


def schema_completedness(base_url, species_id, schema_id, headers_get,
                         hashed_files):
    """ Returns information that indicates which schema data has
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

    # get info about loci and alleles upload
    schema_loci = aux.simple_get_request(base_url, headers_get,
                                         ['species', species_id,
                                          'schemas', schema_id,
                                          'loci', 'data'])
    schema_loci = schema_loci.json()
    if 'Not found' in schema_loci:
        sys.exit('{0}'.format(schema_loci['Not found']))
    else:
        schema_loci = schema_loci['hashes']

    loci_info = {hashed_files[k]: v[1]+[k] for k, v in schema_loci.items()}
    # determine loci that were not fully inserted
    absent_loci = {hashed_files[k]: v[1]+[k]
                   for k, v in schema_loci.items() if all(v[1]) is not True}
    # determine loci without alleles
    absent_alleles = [hashed_files[k]
                      for k, v in schema_loci.items() if v[0] is False]

    return [loci_info, absent_loci, absent_alleles]


def create_schema(base_url, headers_post, species_id, params):
    """ Creates a schema in the Chewie-NS.

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

    schema_post = aux.simple_post_request(base_url, headers_post,
                                          ['species', species_id, 'schemas'],
                                          params)
    schema_status = schema_post.status_code

    # check status code
    if schema_status not in [200, 201]:
        print(schema_status, schema_post.text)
        sys.exit('Could not create new schema in the NS.')

    schema_url = schema_post.json()['url']
    schema_id = schema_url.split('/')[-1]

    return [schema_url, schema_id]


def create_loci_file(schema_files, annotations, schema_dir,
                     species_id, schema_id, loci_prefix,
                     absent_loci=None):
    """ Creates a file with the essential data to insert loci in
        the Chewie-NS.

        Parameters
        ----------
        schema_files : list of str
            List with paths to schema's FASTA files.
        annotations : dict
            Dictionary with loci identifiers as keys and lists
            with annotation terms as values.
        schema_dir : str
            Path to the directory with schema files.
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

    loci_file = os.path.join(schema_dir,
                             '{0}_{1}_loci'.format(species_id, schema_id))
    loci_data = [loci_prefix, []]
    for file in schema_files:
        if absent_loci is None or file in absent_loci:

            locus = os.path.basename(file).split('.fasta')[0]
            locus_hash = aux.hash_file(file, 'rb')

            locus_annotations = annotations[locus]

            loci_data[1].append([locus, locus_hash]+locus_annotations)

    if len(loci_data[1]) > 0:
        with open(loci_file, 'wb') as lf:
            pickle.dump(loci_data, lf)

    return loci_file


def upload_loci_data(loci_file, base_url, species_id,
                     schema_id, headers_post, headers_get, hashed_files):
    """ Sends the necessary data to create loci and associate
        those loci with species and schema.

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

    # compress file with loci data to reduce upload size
    loci_zip_file = '{0}.zip'.format(loci_file)
    aux.file_zipper(loci_file, loci_zip_file)

    zip_url = aux.make_url(base_url, 'species', species_id,
                           'schemas', schema_id, 'loci',
                           'data')

    # upload file as multipart-encoded file
    filename = '{0}_{1}_loci.zip'.format(species_id, schema_id)
    response = aux.upload_file(loci_zip_file, filename,
                               zip_url, headers_post,
                               False)

    response = response.json()
    start_nr_loci = int(response['nr_loci'])
    start_sp_loci = int(response['sp_loci'])
    start_sc_loci = int(response['sc_loci'])

    time_limit = 2100
    current_time = 0
    status = 'Inserting'
    while status != 'Complete' and (current_time < time_limit):
        insertion_status = aux.simple_get_request(
            base_url, headers_get, ['species', species_id,
                                    'schemas', schema_id,
                                    'loci', 'data'])
        insertion_status = insertion_status.json()

        if insertion_status['status'] == 'complete':
            status = 'Complete'
            response_data = insertion_status['hashes']

        nr_loci = int(insertion_status['nr_loci'])
        sp_loci = int(insertion_status['sp_loci'])
        sc_loci = int(insertion_status['sc_loci'])

        inserted_loci = nr_loci - start_nr_loci
        linked_sp = sp_loci - start_sp_loci
        linked_sc = sc_loci - start_sc_loci

        print('\r', '   Inserted {0} loci; '
              'Linked {1} to species; '
              'Linked {2} to schema.'.format(inserted_loci, linked_sp, linked_sc), end='')
        time.sleep(2)
        current_time += 2

    response_data = {hashed_files[k]: v+[k] for k, v in response_data.items()}

    return response_data


def create_alleles_files(schema_files, loci_responses, invalid_alleles,
                         species_name, base_url, species_id,
                         schema_id, user_id):
    """ Creates files with the essential data to insert alleles in the NS.

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
    schema_dir = os.path.dirname(schema_files[0])
    user_uri = '{0}users/{1}'.format(base_url, user_id)
    for file in schema_files:
        locus_uri = loci_responses[file][1][0]

        alleles_sequences = [str(rec.seq)
                             for rec in SeqIO.parse(file, 'fasta')
                             if rec.id not in invalid_alleles]

        post_inputs = [locus_uri, species_name,
                       user_uri, tuple(alleles_sequences)]

        locus_id = locus_uri.split('/')[-1]
        loci_ids.append(locus_id)

        locus_hash = loci_responses[file][-1]
        loci_hashes.append(locus_hash)

        locus_name = os.path.basename(file).split('.fasta')[0]
        loci_names.append(locus_name)

        alleles_file = os.path.join(schema_dir,
                                    '{0}_{1}_{2}'.format(species_id,
                                                         schema_id,
                                                         locus_id))
        alleles_files.append(alleles_file)

        aux.pickle_dumper(alleles_file, post_inputs)

    return [alleles_files, loci_ids, loci_hashes, loci_names]


def upload_alleles_data(alleles_data, length_files, base_url,
                        headers_post, headers_post_bytes,
                        species_id, schema_id):
    """ Uploads files with the data to insert alleles and the
        length values for the sequences of each locus.

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
        data = aux.pickle_loader(current_len)
        data = {current_locus: data[next(iter(data))]}
        data = json.dumps({'content': data, 'locus_hash': locus_hash})

        # send data to the NS
        send_url = aux.make_url(base_url, 'species', species_id,
                                'schemas', schema_id, 'loci',
                                current_locus, 'lengths')

        lengths_res = aux.upload_data(data, send_url, headers_post, False)
        length_status = lengths_res.status_code

        # get path to ZIP archive with data to insert alleles
        current_zip = a[0]

        # send data to insert alleles in the NS
        zip_url = aux.make_url(base_url, 'species', species_id,
                               'schemas', schema_id, 'loci',
                               current_locus, 'data')

        zip_res = aux.upload_file(current_zip, locus_hash,
                                  zip_url, headers_post_bytes,
                                  False)
        zip_status = zip_res.status_code

        # determine if upload was successful
        if length_status not in [200, 201] or zip_status not in [200, 201]:
            failed.append(current_locus)
        elif length_status in [200, 201] and zip_status in [200, 201]:
            uploaded += 1
            print('\r', '    Sent data for alleles of '
                  '{0} loci.'.format(uploaded), end='')

    return failed


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory of the schema to upload.')

    parser.add_argument('-sp', type=str, required=True,
                        dest='species_id',
                        help='The integer identifier or name of the species '
                             'that the schema will be associated to in '
                             'the NS.')

    parser.add_argument('-sn', type=str, required=True,
                        dest='schema_name',
                        help='A brief and meaningful name that '
                             'should help understand the type and content '
                             'of the schema.')

    parser.add_argument('-lp', type=str, required=True,
                        dest='loci_prefix',
                        help='Prefix included in the name of each locus of '
                             'the schema.')

    parser.add_argument('--df', type=str, required=False,
                        dest='description_file', default=None,
                        help='Path to a text file with a description '
                             'about the schema. Markdown syntax is supported '
                             'in order to offer greater customizability of '
                             'the rendered description in the Frontend. '
                             'Will default to the schema\'s name if the user '
                             'does not provide a valid path for a file.')

    parser.add_argument('--a', type=str, required=False,
                        dest='annotations', default=None,
                        help='Path to a TSV file with loci annotations. '
                             'The first column has loci identifiers '
                             '(w/o .fasta extension), the second has user '
                             'annotations and the third has custom '
                             'annotations.')

    parser.add_argument('--cpu', type=int, required=False,
                        dest='cpu_cores', default=1,
                        help='Number of CPU cores that will '
                             'be used in the Schema Pre-processing step.')

    parser.add_argument('--thr', type=int, required=False,
                        default=20, dest='threads',
                        help='Number of threads to use to search for '
                             'annotations on UniProt.')

    parser.add_argument('--ns', type=pv.validate_ns_url, required=False,
                        default='main',
                        dest='nomenclature_server',
                        help='The base URL for the Nomenclature Server. '
                             'The default value, "main", will establish a '
                             'connection to "https://chewbbaca.online/", '
                             '"tutorial" to "https://tutorial.chewbbaca.online/" '
                             'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
                             'Users may also provide the IP address to other '
                             'Chewie-NS instances.')

    parser.add_argument('--continue_up', required=False, action='store_true',
                        dest='continue_up',
                        help='If the process should check if the schema '
                             'upload was interrupted and try to resume it.')

    args = parser.parse_args()

    schema_directory = args.schema_directory
    species_id = args.species_id
    schema_name = args.schema_name
    loci_prefix = args.loci_prefix
    description_file = args.description_file
    annotations = args.annotations
    cpu_cores = args.cpu_cores
    threads = args.threads
    nomenclature_server = args.nomenclature_server
    continue_up = args.continue_up

    return [schema_directory, species_id, schema_name,
            loci_prefix, description_file, annotations,
            cpu_cores, threads, nomenclature_server,
            continue_up]


def main(input_files, species_id, schema_name, loci_prefix, description_file,
         annotations, cpu_cores, threads, base_url, continue_up):

    if 'tutorial' not in base_url:
        token = aux.capture_login_credentials(base_url)
    else:
        token = ''

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
    print('Started at: {0}\n'.format(start_date_str))

    # verify user
    print('-- User Permissions --')
    # GET request headers
    headers_get = cnst.HEADERS_GET_JSON
    headers_get['Authorization'] = token

    # determine current user ID and Role
    if 'tutorial' not in base_url:
        user_id, user_role, user_auth = aux.user_info(base_url, headers_get)
    else:
        user_id, user_role, user_auth = ['', '', True]

    print('User id: {0}'.format(user_id))
    print('User role: {0}'.format(user_role))
    print('Authorized: {0}\n'.format(user_auth))

    # only Admin or Contributor type users can upload schemas
    if not user_auth:
        sys.exit('Current user has no Administrator '
                 'or Contributor permissions.\n'
                 'Not allowed to upload schemas.')

    # POST requests headers
    headers_post = cnst.HEADERS_POST_JSON
    headers_post['Authorization'] = token
    headers_post['user_id'] = user_id
    # POST headers to send binary data
    headers_post_bytes = cnst.HEADERS_POST
    headers_post_bytes['Authorization'] = token
    headers_post_bytes['user_id'] = user_id

    print('-- Parameters Validation --')
    print('Local schema: {0}'.format(input_files))

    # Get schema files from genes list file
    genes_list = os.path.join(input_files, '.genes_list')
    with open(genes_list, 'rb') as gl:
        genes = pickle.load(gl)

    fasta_paths = [os.path.join(input_files, file) for file in genes]
    fasta_paths.sort()

    # total number of loci
    total_loci = len(fasta_paths)
    # total number of alelles
    total_alleles = sum(list(map(aux.count_sequences, fasta_paths)))

    # Get the name of the species from the provided id
    # or vice-versa
    species_info = aux.species_ids(species_id, base_url, headers_get)
    if isinstance(species_info, list):
        species_id, species_name = species_info
        print("Schema's species: {0} (id={1})".format(species_name, species_id))
    else:
        sys.exit('There is no species with the provided identifier in the NS.')

    print('Number of loci: {0}'.format(total_loci))
    print('Number of alleles: {0}\n'.format(total_alleles))

    # verify schema configs
    print('Verifying schema configs...')
    # load schema config file
    configs = aux.read_configs(input_files, '.schema_config')

    # validate arguments values
    ptf_val = pv.validate_ptf(configs.get('prodigal_training_file', ''), input_files)
    bsr_val = pv.bsr_type(configs.get('bsr', ''))
    msl_val = pv.minimum_sequence_length_type(configs.get('minimum_locus_length', ''))
    tt_val = pv.translation_table_type(configs.get('translation_table', ''))
    st_val = pv.size_threshold_type(configs.get('size_threshold', ''))
    cv_val = pv.validate_cv(configs.get('chewBBACA_version', ''))
    ws_val = pv.validate_ws(configs.get('word_size', None))
    cs_val = pv.validate_cs(configs.get('cluster_sim', None))
    rf_val = pv.validate_rf(configs.get('representative_filter', None))
    if_val = pv.validate_if(configs.get('intraCluster_filter', None))

    # dictionary with schema parameters values to send to the NS
    params = {'bsr': bsr_val, 'prodigal_training_file': ptf_val[1],
              'translation_table': tt_val, 'minimum_locus_length': msl_val,
              'chewBBACA_version': cv_val, 'size_threshold': st_val,
              'word_size': ws_val, 'cluster_sim': cs_val,
              'representative_filter': rf_val, 'intraCluster_filter': if_val}

    params_values = list(params.values())

    if '' in params_values:
        sys.exit('Found invalid parameters values and exited.')
    else:
        params_lines = ['  {0}: {1}'.format(k, str(v))
                        for k, v in params.items()
                        if k not in ['prodigal_training_file', 'name']]
        params_text = '\n'.join(params_lines)
        print(params_text)
        print('All configurations successfully validated.\n')
        params['name'] = schema_name
        ptf_file = ptf_val[0]
        ptf_hash = ptf_val[1]

    # determine schema status
    upload_type = schema_status(base_url, headers_get,
                                schema_name, species_id,
                                continue_up)

    if upload_type[0] == 'novel':
        print('New schema name: "{0}" '.format(schema_name))
    else:
        schema_url = '{0}species/{1}/schemas/{2}'.format(base_url, species_id,
                                                         upload_type[1])
        schema_id = schema_url.split('/')[-1]
        print('Schema exists and is incomplete '
              '("{0}", id={1})'.format(schema_name, schema_id))

    # get schema description
    if continue_up is False:
        if description_file is not None and os.path.isfile(description_file) is True:
            # determine file hash
            description_hash = aux.hash_file(description_file, 'rb')
            print('Schema description: {0}'.format(description_file))
        else:
            print('Could not get a description from a file. '
                  'Reset to schema name.')
            description_file = 'schema_description.txt'
            with open(description_file, 'w') as sd:
                sd.write(schema_name)
            description_hash = aux.hash_file(description_file, 'rb')

        params['SchemaDescription'] = description_hash

    print('\n-- Schema Pre-processing --')

    # hash schema files to get unique identifiers based on content
    hashed_files = {aux.hash_file(file, 'rb'): file for file in fasta_paths}
    print('Determining data to upload...')
    absent_loci = fasta_paths
    if upload_type[0] == 'incomplete':

        loci_info, absent_loci, fasta_paths = schema_completedness(base_url, species_id,
                                                                   upload_type[1], headers_get,
                                                                   hashed_files)

    print('  Loci to create and associate with species '
          'and schema: {0}'.format(len(absent_loci)))
    print('  Loci without the full set of alleles: '
          '{0}\n'.format(len(fasta_paths)))

    # create inputs for QC step
    inputs = [(file,
               file.split('/')[-1].split('.fasta')[0],
               int(params['translation_table']),
               0,
               None) for file in fasta_paths]

    # validate schema data and create files with translated sequences
    print('Translating sequences based on schema configs...')
    qc_results = []
    genes_pools = multiprocessing.Pool(processes=cpu_cores)
    rawr = genes_pools.map_async(quality_control, inputs,
                                 callback=qc_results.extend)
    rawr.wait()

    # get invalid alleles
    invalid_alleles = [r[2] for r in qc_results]
    invalid_alleles = list(itertools.chain.from_iterable(invalid_alleles))
    invalid_identifiers = set([r[0] for r in invalid_alleles])

    print('  Found a total of {0} invalid '
          'alleles.\n'.format(len(invalid_identifiers)))

    # list translated sequences files
    dna_files = [r[0] for r in qc_results]
    prot_files = [r[1] for r in qc_results]

    # determine loci missing annotations
    miss_annotation = [pf for pf in prot_files
                       if pf.split('_prots')[0] + '.fasta' in absent_loci]

    print('Loci missing UniProt annotation: {0}'.format(len(miss_annotation)))

    queries_files = []
    if len(miss_annotation) > 0:
        # create SPARQL queries to query UniProt SPARQL endpoint
        print('Creating SPARQL queries to search UniProt for annotations...')
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            for res in executor.map(create_uniprot_queries, miss_annotation):
                queries_files.append(res)

        # multithreaded annotation searching
        print('Searching for annotations on UniProt...')
        loci_annotations = {}
        total_found = 0
        total_loci = len(queries_files)
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            # Start the load operations and mark each future with its URL
            for res in executor.map(get_annotation, queries_files):
                loci_annotations[res[0]] = res[1:]
                total_found += 1
                print('\r', 'Searched annotations for '
                      '{0}/{1} loci'.format(total_found, total_loci), end='')

        # get user and custom annotations
        if os.path.isfile(str(annotations)) is True:
            user_annotations = import_annotations(annotations)
            valid = 0
            for a in loci_annotations:
                if a in user_annotations:
                    loci_annotations[a].extend(user_annotations[a])
                    valid += 1
                else:
                    loci_annotations[a].extend(['N/A', 'N/A'])
            print('\nUser provided valid annotations for {0} '
                  'loci.'.format(valid))
        else:
            loci_annotations = {locus: loci_annotations[locus]+['N/A', 'N/A'] for locus in loci_annotations}
            if annotations is None:
                print('\nUser did not provide file with annotations.')
            else:
                print('\nInvalid annotations file value.')

    # convert parameters to string type because the Chewie-NS
    # expects strings
    params = {k: str(v) for k, v in params.items()}
    params['schema_hashes'] = list(hashed_files.keys())

    # start sending data
    print('\n-- Schema Upload --')
    # Build the new schema URL and POST to NS
    if continue_up is False:
        schema_url, schema_id = create_schema(base_url,
                                              headers_post,
                                              species_id,
                                              params)
        # send file with description
        description_uri = aux.make_url(base_url, 'species', species_id,
                                       'schemas', schema_id, 'description')

        desc_res = aux.upload_file(description_file, description_hash,
                                   description_uri, headers_post_bytes,
                                   False)

        print('Created schema with name {0} '
              '(id={1}).\n'.format(schema_name, schema_id))
    else:
        print('Continuing upload of schema with '
              'name {0} (id={1})\n'.format(schema_name, schema_id))

    # start creating new loci and adding/linking alleles
    print('Loci data:')

    if upload_type[0] == 'incomplete':
        if len(absent_loci) == 0:
            response_data = {k: [v[0], True, True, True, v[3]]
                             for k, v in loci_info.items()}
            print('All loci were inserted in a previous process.\n')
        else:
            print('  Collecting loci data...')
            loci_file = create_loci_file(dna_files, loci_annotations,
                                         input_files, species_id,
                                         schema_id, loci_prefix,
                                         absent_loci)
            print('  Sending data to the NS...')
            absent_data = upload_loci_data(loci_file, base_url,
                                           species_id, schema_id,
                                           headers_post_bytes, headers_get,
                                           hashed_files)
            response_data = {k: [v[0], True, True, True, v[3]]
                             for k, v in loci_info.items() if k not in absent_data}
            response_data = {**absent_data, **response_data}
            print('\n  The NS completed the insertion of {0} '
                  'loci.\n'.format(len(absent_loci)))
    else:
        print('  Collecting loci data...')
        loci_file = create_loci_file(dna_files, loci_annotations,
                                     input_files, species_id,
                                     schema_id, loci_prefix)
        print('  Sending data to the NS...')
        response_data = upload_loci_data(loci_file, base_url,
                                         species_id, schema_id,
                                         headers_post_bytes, headers_get,
                                         hashed_files)

        print('\n  The NS completed the insertion of {0} '
              'loci.\n'.format(len(response_data)))

    # create files with info for posting alleles
    print('Alleles data:')
    print('  Collecting alleles data...')
    (alleles_files, loci_ids, loci_hashes,
     loci_names) = create_alleles_files(dna_files, response_data,
                                        invalid_identifiers, species_name,
                                        base_url, species_id,
                                        schema_id, user_id)
    # determine length of all alleles per locus
    length_files = create_lengths_files(dna_files, input_files)

    # zip all files to reduce upload size
    print('  Compressing files with alleles data...')
    zipped_files = ['{0}.zip'.format(file) for file in alleles_files]
    list(map(aux.file_zipper, alleles_files, zipped_files))
    alleles_data = list(zip(zipped_files, loci_ids, loci_hashes, loci_names))

    print('  Sending alleles data to the NS...')
    # send POST with file contents and process each file in the NS
    failed = upload_alleles_data(alleles_data, length_files, base_url,
                                 headers_post, headers_post_bytes,
                                 species_id, schema_id)

    if len(failed) > 0:
        sys.exit('Could not upload data for alleles of following loci:\n'
                 '{0}\n Please retry with the "--continue_up" option and '
                 'contact the NS Admin if the problem '
                 'persists.'.format(','.join(failed)))
    else:
        # send training file to NS
        print('\n\nUploading Prodigal training file...')
        ptf_url = aux.make_url(base_url, 'species', species_id,
                               'schemas', schema_id, 'ptf')
        ptf_res = aux.upload_file(ptf_file, ptf_hash,
                                  ptf_url, headers_post_bytes,
                                  False)
        print(list(ptf_res.json().values())[0])
        print('\nThe NS has received the data and will insert '
              'the alleles into the database.')
        print('Schema will be available for download as soon as '
              'the process has completed.')
        print('Schema information will also be available on the '
              'NS website.\n')

    # delete all intermediate files
    print('Removing intermediate files...')
    aux.remove_files(prot_files)
    aux.remove_files(queries_files)
    aux.remove_files(length_files)
    aux.remove_files(alleles_files)
    aux.remove_files(zipped_files)

    if len(absent_loci) > 0:
        os.remove(loci_file)
        os.remove('{0}.zip'.format(loci_file))

    end_date = dt.datetime.now()
    end_date_str = dt.datetime.strftime(end_date, '%Y-%m-%dT%H:%M:%S')

    delta = end_date - start_date
    minutes, seconds = divmod(delta.total_seconds(), 60)

    print('\nFinished at: {0}'.format(end_date_str))
    print('Elapsed time: {0:.0f}m{1:.0f}s'.format(minutes, seconds))


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4],
         args[5], args[6], args[7], args[8], args[9])
