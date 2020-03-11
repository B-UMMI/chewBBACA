#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import json
import time
import pickle
import hashlib
import argparse
import requests
import concurrent.futures
from getpass import getpass
from collections import Counter
from SPARQLWrapper import SPARQLWrapper, JSON

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

from utils import constants as cnst
from utils import auxiliary_functions as aux


virtuoso_server = SPARQLWrapper('http://sparql.uniprot.org/sparql')


def binary_file_hash(binary_file):
    """
    """

    with open(binary_file, 'rb') as bf:
        file_hash = hashlib.blake2b()
        file_text = bf.read()
        file_hash.update(file_text)
        file_hash = file_hash.hexdigest()

    return file_hash


def retrieve_schema_info(schemas_list, schema_desc):
    """
    """

    schema_exists = False
    for s in schemas_list:
        current_desc = s['name']['value']
        if current_desc == schema_desc:
            schema_exists = True
            schema_url = s['schemas']['value']
            schema_id = schema_url.split('/')[-1]

    if schema_exists:
        return [schema_url, schema_id]
    else:
        return 404


def determine_upload(local_schema_loci, ns_schema_loci,
                     ns_schema_locid_map, local_path,
                     base_url, headers_get):
    """
    """

    missing = []
    incomplete = []
    for locus in local_schema_loci:
        local_locus = locus

        if local_locus in ns_schema_loci:
            local_file = os.path.join(local_path, locus)
            local_sequences = [str(rec.seq)
                               for rec in SeqIO.parse(local_file, 'fasta')]

            ns_uri = ns_schema_locid_map[locus][0]
            ns_locus_id = ns_uri.split('/')[-1]
            ns_info = aux.simple_get_request(base_url, headers_get,
                                             ['loci', ns_locus_id, 'fasta'])
            ns_sequences = [seq['nucSeq']['value']
                            for seq in ns_info.json()['Fasta']]

            local_set = set(local_sequences)
            ns_set = set(ns_sequences)

            if len(ns_set) > len(local_set):
                message = ('A locus in the NS has more sequences '
                           'than the local locus.\nLocal schema '
                           'is not the original.')
                return (400, message)

            ns_diff = ns_set - local_set
            if len(ns_diff) > 0:
                message = ('A locus in the NS has sequences '
                           'that are not in the local locus.'
                           '\nLocal schema is not the original.')
                return (401, message)

            local_diff = local_set - ns_set
            if len(local_diff) > 0:
                incomplete.append(locus)

        else:
            missing.append(locus)

    upload = missing + incomplete
    # incomplete_text = ', '.join(incomplete)
    # print('Incomplete: {0}'.format(incomplete_text))
    print('Incomplete: {0}'.format(len(upload)))

    return upload


def create_uniprot_queries(fasta_paths):
    """ Create queries to search for protein annotations
        on UniProt and save the queries in a binary file.

        Args:
            fasta_paths (list): list of paths to fasta files.
        Returns:
            queries_files (list): list of paths to binary files
            with queries to search for protein annotations on
            UniProt.
    """

    queries_files = []
    for file in fasta_paths:

        dna_seqs = [(record.id, str(record.seq))
                    for record in SeqIO.parse(file, 'fasta')]
        protein_seqs = [str(translate_sequence(rec[1], 11))
                        for rec in dna_seqs]

        # create queries
        # do not create duplicated queries with same protein
        unique_prots = []
        for prot in protein_seqs:
            if prot not in unique_prots:
                unique_prots.append(prot)

        queries = [aux.uniprot_query(prot) for prot in unique_prots]
        # save in binary file
        binary_file = '{0}_up'.format(file.split('.fasta')[0])
        with open(binary_file, 'wb') as bup:
            pickle.dump(queries, bup)

        queries_files.append((file, binary_file))

    return queries_files


def check_configs(file_path, input_path, schema_desc):
    """ Checks and validates each parameter used to create
        a schema.

        Args:
            file_path (str): path to the binary file
            with the parameters used to create the schema.
        Returns:
            A list with following elements:
                - a boolean indicating if all
                parameters were validated (True)
                or not (False);
                - a list with the reason/reasons
                for not validating the schema
                parameters or a dictionary with
                all parameters if they were all
                valid.
    """

    params = {}
    messages = []

    if os.path.isfile(file_path):
        print('Found config file. Loading configs...')
        # Load configs dictionary
        with open(file_path, 'rb') as cf:
            configs = pickle.load(cf)

        params['name'] = schema_desc
        # Determine if the configs are valid
        # Prodigal training file
        schema_ptf_name = configs.get('prodigal_training_file', 'na')
        schema_ptf_path = os.path.join(input_path, schema_ptf_name)
        if os.path.isfile(schema_ptf_path):
            print('Found valid training file in schema directory.')
            params['prodigal_training_file'] = schema_ptf_path
        else:
            message = 'Could not find valid training file in schema directory.'
            messages.append(message)

        # BSR value
        schema_bsr = configs.get('bsr', 'na')
        try:
            schema_bsr = float(schema_bsr)
            if schema_bsr > 0.0 and schema_bsr < 1.0:
                print('Schema created with BSR value of '
                      '{0}.'.format(schema_bsr))
                params['bsr'] = str(schema_bsr)
            else:
                raise ValueError('Value is not contained in the '
                                 '[0.0, 1.0] range.')
        except ValueError:
            message = ('Invalid BSR value of {0}. BSR value must be contained'
                       ' in the [0.0, 1.0] interval.'.format(schema_bsr))
            messages.append(message)

        # Minimum sequence length
        schema_ml = int(configs.get('minimum_locus_length', 'na'))
        try:
            schema_ml = int(schema_ml)
            if schema_ml >= 0:
                print('Schema created with a minimum sequence length '
                      'parameter of {0}.'.format(schema_ml))
                params['minimum_locus_length'] = str(schema_ml)
            else:
                raise ValueError('Invalid minimum sequence length value. '
                                 'Must be equal or greater than 0.')
        except ValueError:
            message = ('Invalid minimum sequence length value used to '
                       'create schema. Value must be a positive integer.')
            messages.append(message)

        # translation table
        schema_gen_code = int(configs.get('translation_table', 'na'))
        if schema_gen_code in cnst.GENETIC_CODES:
            genetic_code_desc = cnst.GENETIC_CODES[schema_gen_code]
            print('Schema genes were predicted with genetic code '
                  '{0} ({1}).'.format(schema_gen_code, genetic_code_desc))
            params['translation_table'] = str(schema_gen_code)
        else:
            message = ('Genetic code used to create schema is not valid.')
            messages.append(message)

        # chewie version
        schema_chewie_version = configs.get('chewBBACA_version', 'na')
        if schema_chewie_version in cnst.CHEWIE_VERSIONS:
            chewie_version = cnst.CHEWIE_VERSIONS[cnst.CHEWIE_VERSIONS.index(schema_chewie_version)]
            print('Schema created with chewBBACA v{0}.'.format(chewie_version))
            params['chewBBACA_version'] = schema_chewie_version
        else:
            message = ('Schema created with chewBBACA version that '
                       'is not suitable to work with the NS.')
            messages.append(message)

        # clustering word_size
        word_size = configs.get('word_size', 'na')
        try:
            if word_size == 'na':
                print('Schema created with a chewBBACA version '
                      'that did not use clustering.')
            else:
                word_size = int(word_size)
                if word_size >= 4:
                    print('Schema created with a clustering word size '
                          'value of {0}.'.format(word_size))
                else:
                    raise ValueError('Word size for the clustering step '
                                     'must be equal or greater than 4.')
            params['word_size'] = str(word_size)
        except ValueError:
            message = ('Schema created with invalid clustering word '
                       'size value.')
            messages.append(message)

        # clustering similarity threshold
        cluster_sim = configs.get('cluster_sim', 'na')
        try:
            if cluster_sim == 'na':
                print('Schema created with a chewBBACA version '
                      'that did not use clustering.')
            else:
                cluster_sim = float(cluster_sim)
                if cluster_sim >= 0.0:
                    print('Schema created with a clustering threshold '
                          'value of {0}.'.format(cluster_sim))
                else:
                    raise ValueError('Clustering similarity threshold value '
                                     'must be contained in the [0.0, 1.0] '
                                     'interval.')
            params['cluster_sim'] = str(cluster_sim)
        except ValueError:
            message = ('Schema created with invalid clustering '
                       'threshold value.')
            messages.append(message)

        # representative_filter
        representative_filter = configs.get('representative_filter', 'na')
        try:
            if representative_filter == 'na':
                print('Schema created with a chewBBACA version '
                      'that did not use clustering.')
            else:
                representative_filter = float(representative_filter)
                if representative_filter >= 0.0 and representative_filter <= 1.0:
                    print('Schema created with a representative filter '
                          'value of {0}.'.format(representative_filter))
                else:
                    raise ValueError('Representative filter threshold value '
                                     'must be contained in the [0.0, 1.0] '
                                     'interval.')
            params['representative_filter'] = str(representative_filter)
        except ValueError:
            message = ('Schema created with invalid representative filter value.')
            messages.append(message)

        # intraCluster_filter
        intraCluster_filter = configs.get('intraCluster_filter', 'na')
        try:
            if intraCluster_filter == 'na':
                print('Schema created with a chewBBACA version '
                      'that did not use clustering.')
            else:
                intraCluster_filter = float(intraCluster_filter)
                if intraCluster_filter >= 0.0 and intraCluster_filter <= 1.0:
                    print('Schema created with a intraCluster filter '
                          'value of {0}.'.format(intraCluster_filter))
                else:
                    raise ValueError('intraCluster filter threshold value '
                                     'must be contained in the [0.0, 1.0] '
                                     'interval.')
            params['intraCluster_filter'] = str(intraCluster_filter)
        except ValueError:
            message = ('Schema created with invalid intraCluster filter '
                       'value.')
            messages.append(message)

        if 'message' not in locals():
            message = 'All configurations successfully validated.'
            print(message)
            return [True, params]
        else:
            for m in messages:
                print(m)
            return [False, messages]
    else:
        message = ('Could not find a valid config file. Cannot upload'
                   ' schema without checking for valid parameters values.')
        print(message)

        return [False, messages]


def check_schema_status(status_code, species_name, upload_type):
    """ Checks the schema post status and determines
        if the schema was successfully created in the NS.

        Args:
            status_code (int): schema post status code.
            species_name (str): name of the species.
        Returns:
            message (str): message indicating if the
            schema post was successful or not and why.
    """

    if status_code in [200, 201]:
        if upload_type == 'de novo':
            message = ('A new schema for {0} was created '
                       'succesfully.'.format(species_name))
        elif upload_type == 'continue':
            message = ('Schema exists. Will try to continue upload.')
    else:
        if upload_type == 'de novo':
            if status_code == 403:
                message = ('{0}: No permission to load '
                           'schema.'.format(status_code))
            elif status_code == 404:
                message = ('{0}: Cannot upload a schema for a species '
                           'that is not in NS.'.format(status_code))
            elif status_code == 409:
                message = ('{0}: Cannot upload a schema with the same '
                           'description as a schema that is in the '
                           'NS.'.format(status_code))
            else:
                message = '{0}: Could not insert schema.'.format(status_code)
        elif upload_type == 'continue':
            message = ('{0}: Cannot continue uploading data for a schema that '
                       'does not exist.'.format(status_code))

    return message


def post_locus(base_url, headers_post, locus_prefix, keep_file_name, gene):
    """ Adds a new locus to the NS.

        Args:
            base_url (str): the base URI for the NS, used to concatenate
            with a list of elements and obtain endpoints URIs.
            headers_post (dict): headers for the POST method used to
            insert data into the NS.
            locus_prefix (str): prefix for the locus identifier.
            keep_file_name (bool): boolean value indicating if the original
            schema file identifier should be stored in the NS.
            gene (str): identifier of the original schema file.
        Returns:
            loci_url (str): API endpoint for the new locus.
    """

    # Build the url for loci/list
    url_loci = aux.make_url(base_url, 'loci', 'list')

    # Define POST request parameters
    params = {}
    params['prefix'] = locus_prefix

    if keep_file_name:
        params['locus_ori_name'] = gene

    # Add locus to species
    res = requests.post(url_loci, data=json.dumps(params),
                        headers=headers_post, timeout=30)

    res_status = res.status_code
    if res_status == 409:
        message = '{0}: Locus already exists on NS.'.format(res_status)
    elif res_status == 404:
        message = '{0}: Species not found.'.format(res_status)
    elif res_status == 403:
        message = ('{0}: Unauthorized. No permission to add '
                   'new locus.'.format(res_status))
    elif res_status == 400:
        message = ('{0}: Please provide a valid locus '
                   'prefix.'.format(res_status))

    if 'message' in locals():
        return [False, message]
    else:
        loci_url = res.json()['uri']
        return [True, loci_url]


def get_annotation(sparql_queries):
    """ Queries the UniProt SPARQL endpoint to retrieve
        protein annotations.

        Args:
            sparql_queries (list): a list with a tuple.
            The first element of each tuple is the locus
            identifier and the second element is a list with
            the SPARQL queries and DNA sequences for each
            allele of the locus.
        Returns:
            A tuple with the following elements:
                - locus (str): the locus identifier;
                - prev_name (str): the annotation description.
                - label (str): a descriptive label;
                - url (str): the URL to the UniProt page about
                the protein;
                - dna_seq_to_ns (list): the list of DNA sequences
                that belong to the locus and should be added to the NS.
    """

    locus = sparql_queries[0]
    queries_file = sparql_queries[1]

    # load queries for locus
    with open(queries_file, 'rb') as bup:
        queries = pickle.load(bup)

    virtuoso_server.setReturnFormat(JSON)
    virtuoso_server.setTimeout(10)

    prev_url = ''
    prev_name = ''
    prev_label = ''
    found = False
    unpreferred_names = ['Uncharacterized protein',
                         'hypothetical protein',
                         'DUF',
                         '']

    a = 0
    # define maximum number of tries
    max_tries = 10
    while found is False:

        virtuoso_server.setQuery(queries[a])

        try:
            result = virtuoso_server.query().convert()

            name, url, label = aux.select_name(result)

            if prev_name == '' and name != '':
                prev_name = name
                prev_label = label
                prev_url = url
                if prev_name not in unpreferred_names:
                    found = True

            elif prev_name in unpreferred_names and name not in unpreferred_names:
                prev_name = name
                prev_label = label
                prev_url = url
                found = True

        # retry if the first query failed
        except Exception:
            pass

        a += 1
        if a == max_tries or a == len(queries):
            found = True

    if prev_name == '':
        prev_name = 'not found'
    if prev_label == '':
        prev_label = 'not found'
    if prev_url == '':
        # virtuoso needs a string that looks like an URL
        prev_url = 'http://not.found.org'

    return (locus, prev_name, prev_label, prev_url)


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.
        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    try:
        protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

        return protseq

    except TranslationError as e:
        return e


def post_species_loci(url, species_id, locus_id, headers_post):
    """ Adds a new loci to an existing species.

        Args:
            url (str): NS base URI.
            species_id (int): integer identifier for the
            species in the NS.
            locus_id (int): identifier of the locus that
            will be inserted.
            headers_post (dict): headers for the POST method used to
            insert data into the NS.
        Returns:
            True if the POST was successful.
    """

    # Define POST request parameters
    params = {}
    params['locus_id'] = locus_id

    # Build the url for the loci of the new schema
    url_species_loci = aux.make_url(url, 'species', species_id, 'loci')

    # insert new locus
    res = requests.post(url_species_loci, data=json.dumps(params),
                        headers=headers_post, timeout=30)

    if res.status_code > 201:
        message = ('{0}: Failed to link locus to '
                   'species.'.format(res.status_code))

    if 'message' in locals():
        return [False, message]
    else:
        species_loci_url = '{0}/{1}'.format(url_species_loci, locus_id)
        return [True, species_loci_url]


def post_schema_loci(loci_url, schema_url, headers_post):
    """ Adds a new loci to an existing schema.

        Args:
            loci_url (str): URI for the loci that will be added.
            schema_url (str): URI for the schema that we want to
            add the locus to.
            headers_post (dict): headers for the POST method used to
            insert data into the NS.
        Returns:
            True if the POST was successful.
    """

    # Get the new loci id from the new loci url
    new_loci_id = str(int(loci_url.split('/')[-1]))

    # Define POST request parameters
    params = {}
    params['loci_id'] = new_loci_id

    # Build the url for the loci of the new schema
    url_schema_loci = aux.make_url(schema_url, 'loci')

    res = requests.post(url_schema_loci, data=json.dumps(params),
                        headers=headers_post, timeout=30)

    if res.status_code > 201:
        message = ('{0}: Failed to link locus to '
                   'schema.'.format(res.status_code))

    if 'message' in locals():
        return [False, message]
    else:
        schema_loci_url = '{0}/{1}'.format(url_schema_loci, new_loci_id)
        return [True, schema_loci_url]


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, dest='input_files', required=True,
                        help='Path to folder containing the schema fasta '
                        'files (alternatively, a file with a list of paths '
                        'to fasta files).')

    parser.add_argument('-sp', type=str, dest='species_id', required=True,
                        help='The identifier for the schemas species in '
                        'the NS. Can be the species name of the integer '
                        'identifier for that species in the NS.')

    parser.add_argument('-sdesc', type=str, dest='schema_desc', required=True,
                        help='A brief description for the schema. It is '
                        'important to pass a meaningful value')

    parser.add_argument('-lprefix', type=str, dest='loci_prefix',
                        required=True,
                        help='Prefix added to the identifier of each schema '
                        'gene. For instance, ACIBA will produce '
                        'ACIBA00001.fasta.')

    parser.add_argument('--thr', type=int, required=False, dest='threads',
                        default=30, help='The number of threads to use. '
                        'The process will use multithreading to search for '
                        'annotations on UniProt and to send data to the NS '
                        '(default=30).')

    # change this so that it gets default from config file
    parser.add_argument('--url', type=str, required=False, dest='base_url',
                        default='http://127.0.0.1:5000/NS/api/',
                        help='The base URL for the NS server. Will be used '
                        'to create endpoints URLs.')

    parser.add_argument('--cont', type=str, dest='continue_up', required=False,
                        default='no', choices=['no', 'yes'],
                        help='Flag used to indicate if the process should '
                        'try to continue a schema upload that crashed '
                        '(default=no.')

    args = parser.parse_args()

    return [args.input_files, args.species_id, args.schema_desc,
            args.loci_prefix, args.threads, args.base_url,
            args.continue_up]


def main(input_files, species_id, schema_desc, loci_prefix, threads,
         base_url, continue_up):

    check_cds = True

    # login with master key
    login_key = False
    if login_key:
        pass
    # if the login key is not found ask for credentials
    else:
        print('\nCould not find private key.')
        print('\nPlease provide login credentials:')
        user = input('USERNAME: ')
        password = getpass('PASSWORD: ')
        print()
        # get token
        token = aux.login_user_to_NS(base_url, user, password)
        # if login was not successful, stop the program
        if token is False:
            message = '403: Invalid credentials.'
            print(message)
            return message

    start = time.time()

    # Define the headers of the requests
    headers_get = {'Authorization': token,
                   'accept': 'application/json'}

    # verify user role to check permission
    user_info = aux.simple_get_request(base_url, headers_get,
                                       ['user', 'current_user'])
    user_info = user_info.json()
    user_role = any(role in user_info['roles']
                    for role in ['Admin', 'Contributor'])

    if not user_role:
        print('\n403: Current user has no Administrator '
              'or Contributor permissions.\n'
              'Not allowed to upload schemas.')
        return 403

    user_id = str(user_info['id'])
    headers_post = {'Authorization': token,
                    'Content-type': 'application/json',
                    'accept': 'application/json',
                    'user_id': user_id}

    # check if there is config file and load it
    config_file = os.path.join(input_files, '.schema_config')
    configs_validation = check_configs(config_file, input_files, schema_desc)
    if configs_validation[0] is not True:
        print(configs_validation)
        return configs_validation
    else:
        params = configs_validation[1]

    # Check if user provided a list of genes or a folder
    fasta_paths = [os.path.join(input_files, file)
                   for file in os.listdir(input_files) if '.fasta' in file]
    fasta_paths.sort()

    # Get the name of the species from the provided id
    # or vice-versa
    species_info = aux.species_ids(species_id, base_url, headers_get)
    if isinstance(species_info, list):
        species_id, species_name = species_info
        print('\nNS species with identifier {0} is {1}.'.format(species_id,
                                                                species_name))
    else:
        print('\nThere is no species with the provided identifier in the NS.')
        return 1

    # translate loci sequences and contruct SPARQL queries to query UniProt
    # create queries for each file and save in binary with pickle
    print('\nCreating queries to search UniProt for annotations...')
    # we have to verify if all sequences can be translated!!!
    queries_files = create_uniprot_queries(fasta_paths)
    print('Done.')

    # find annotations for all loci with multithreading
    print('\nSearching for annotations on UniProt...')
    results = []
    total_found = 0
    total_loci = len(queries_files)
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(get_annotation, queries_files):
            results.append(res)
            total_found += 1
            print('\r', 'Found annotations for '
                  '{0}/{1} loci.'.format(total_found, total_loci), end='')

    # delete binary files with queries
    for file in queries_files:
        os.remove(file[1])

    # start sending data
    print('\n\nSending data to NS...')

    # schema parameters to send to NS have sftp path
    # for the Prodigal training file
    params['prodigal_training_file'] = os.path.join(cnst.SFTP_PTF_AVAI,
                                                    params['prodigal_training_file'])

    # Build the new schema URL and POST to NS
    if continue_up == 'no':
        print('\nCreating new schema...')
        # schema is created in locked state
        schema_post = aux.simple_post_request(base_url, headers_post,
                                              ['species', species_id, 'schemas'],
                                              params)
        schema_status = schema_post.status_code
    elif continue_up == 'yes':
        print('\nChecking if schema already exists...')
        schema_get = aux.simple_get_request(base_url, headers_get,
                                            ['species', species_id, 'schemas'])
        schema_get_status = schema_get.status_code
        species_schemas = schema_get.json()
        if schema_get_status in [200, 201]:
            # determine if there is a schema for current
            # species with same description
            schema_info = retrieve_schema_info(species_schemas, schema_desc)

            if isinstance(schema_info, int):
                print('\nCannot continue uploading to a schema that '
                      'does not exist.')
                return 404
            else:
                schema_url, schema_id = schema_info
        else:
            print('\nCould not retrieve schemas for current species.')
            return 404

        # compare list of genes, if they do not intersect, halt process
        # get list of loci for schema in NS
        ns_loci_get = aux.simple_get_request(base_url, headers_get,
                                             ['species', species_id,
                                              'schemas', schema_id,
                                              'loci'])
        # get loci files names from response
        print(ns_loci_get.json())
        ns_schema_loci = []
        ns_schema_locid_map = {}
        if 'message' not in ns_loci_get.json():
            for l in ns_loci_get.json()['Loci']:
                locus_file = l['original_name']['value']
                ns_schema_loci.append(locus_file)
                locus_uri = l['locus']['value']
                locus_name = l['name']['value']
                ns_schema_locid_map[locus_file] = (locus_uri, locus_name)

        # get list of loci for schema to upload
        local_schema_loci = [l[0].split('/')[-1] for l in results]

        local_loci_set = set(local_schema_loci)
        ns_loci_set = set(ns_schema_loci)

        # verify that the number of loci in NS is not greater than
        # in the local schema
        ns_loci_diff = ns_loci_set - local_loci_set
        if len(ns_loci_diff) > 0:
            print('NS schema has loci that are not in the local schema.')
            return 400

        if len(ns_loci_set) > len(local_loci_set):
            print('NS schema has more loci than the local schema.')
            return 400
        elif len(ns_loci_set) < len(local_loci_set):
            print('NS schema has less loci than local schema.')
            absent_loci = list(local_loci_set-ns_loci_set)
            absent_text = [absent_loci[i:i+4] for i in range(0, len(absent_loci), 4)]
            absent_text = [('{:30} '*len(g)).format(*g) for g in absent_text]
            absent_text = '\n'.join(absent_text)
            print('Absent loci: \n{0}'.format(absent_text))
        elif len(ns_loci_set) == len(local_loci_set):
            print('NS and local schemas have the same number of loci.')

        # if the set of loci is equal, check sequences in each locus
        # if a locus in the NS has more sequences than one in the local
        # set, halt process
        upload = determine_upload(local_schema_loci, ns_schema_loci,
                                  ns_schema_locid_map, input_files,
                                  base_url, headers_get)

        if isinstance(upload, tuple):
            return upload[0]
        elif len(upload) == 0:
            print('Local and NS schemas are identical. Nothing left to do!')
            return 'I am such a happy tato'

        results = [list(res) for res in results if
                   any(locus in res[0] for locus in upload)]

        for r in range(len(results)):
            lfile = (results[r][0]).split('/')[-1]
            if lfile in ns_schema_locid_map:
                results[r].append(ns_schema_locid_map[lfile][0])

        schema_status = 200

    # check status code
    # add other prints for cases that fail so that users see a print explaining
    upload_type = 'de novo' if continue_up == 'no' else 'continue'
    schema_insert = check_schema_status(schema_status,
                                        species_name,
                                        upload_type)
    if schema_status in [200, 201]:
        print('{0}'.format(schema_insert))
    else:
        print('{0}'.format(schema_insert))
        return schema_insert

    if continue_up == 'no':
        schema_url = schema_post.json()['url']
        schema_id = schema_url.split('/')[-1]

    # Get the new schema url from the response
    print('Schema description: {0}'.format(schema_desc))
    print('Schema URI: {0}\n'.format(schema_url))

    # start creating new loci and adding/linking alleles
    print('Creating loci and adding/linking alleles...\n')

    # create new loci, link them to species and to new schema and add alleles
    failed = 0
    repeated = 0
    invalid_cds = 0
    new_alleles = 0
    inserted_loci = 0
    linked_alleles = 0
    hash_collisions = 0
    total_loci = len(results)
    for locus in results:

        locus_file = locus[0]
        name = locus[1]
        label = locus[2]
        url = locus[3]
        allele_seq_list = [str(rec.seq)
                           for rec in SeqIO.parse(locus_file, 'fasta')]
        locus_basename = os.path.basename(locus_file)

        if len(locus) == 5:
            # re-upload alleles to existing locus
            new_loci_url = locus[4]
            print('Re-uploading locus: {0}'.format(new_loci_url))

        else:
            # Create a new locus
            new_loci_status, new_loci_url = post_locus(base_url, headers_post,
                                                       loci_prefix, True,
                                                       locus_basename)
            if new_loci_status is False:
                print('{0}'.format(new_loci_url))
                continue
            elif new_loci_status is True:
                print('Created new locus: {0}'.format(new_loci_url))

            # Get the new loci ID
            new_loci_id = new_loci_url.split('/')[-1]

            # Associate the new loci id to the species
            species_link_status, species_link_url = post_species_loci(base_url,
                                                                      species_id,
                                                                      new_loci_id,
                                                                      headers_post)
            if species_link_status is False:
                print('{0}'.format(species_link_url))
                continue
            elif species_link_status is True:
                print('Linked new locus to species: '
                      '{0}'.format(species_link_url))

            # Associate the new loci id to the new schema
            schema_loci_status, schema_link_url = post_schema_loci(new_loci_url,
                                                                   schema_url,
                                                                   headers_post)
            if schema_loci_status is False:
                print('{0}'.format(schema_link_url))
                continue
            elif schema_loci_status is True:
                print('Linked new locus to schema: '
                      '{0}'.format(schema_link_url))

        # create inputs
        post_inputs = aux.create_allele_data(allele_seq_list, new_loci_url,
                                             name, label, url, species_name,
                                             check_cds, headers_post, user_id,
                                             1)

        post_results = []
        total_inserted = 0
        total_alleles = len(allele_seq_list)
        # if the locus has a big number of sequences, use multithreading
        # do not use too much workers as it will make the process hang
        # between the end of each
        if total_alleles > 20:
            workers = threads
        else:
            workers = 1

        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            # Start the load operations and mark each future with its URL
            for res in executor.map(aux.post_allele, post_inputs):
                post_results.append(res)
                total_inserted += 1
                print('\r', 'Processed {0}/{1} '
                      'alleles.'.format(total_inserted, total_alleles), end='')

        requests_res = [r.json() for r in post_results]
        res_definitions = [list(r.keys())[0] for r in requests_res]
        res_collisions = [r for r in requests_res if 'HASH COLLISION' in r]
        if len(res_collisions) > 0:
            print(res_collisions)
        res_counts = Counter(res_definitions)

        inserted_loci += 1
        failed += res_counts.get('FAIL', 0)
        new_alleles += res_counts.get('ADD', 0)
        linked_alleles += res_counts.get('LINK', 0)
        invalid_cds += res_counts.get('INVALID CDS', 0)
        repeated += res_counts.get('REPEATED ALLELE', 0)
        hash_collisions += res_counts.get('HASH COLLISION', 0)

        successful = res_counts.get('ADD', 0)+res_counts.get('LINK', 0)
        print('\nSuccessfully sent {0} alleles.'.format(successful))
        print('Could not insert {0} alleles.'.format(total_alleles-successful))
        print('Inserted {0}/{1} loci.\n'.format(inserted_loci, total_loci))

    # print final statistics
    print('\tNew alleles: {0}\n'
          '\tLinked alleles: {1}\n'
          '\tFailed: {2}\n'
          '\tInvalid CDS: {3}\n'
          '\tRepeated allele: {4}\n'
          '\tHash collision: {5}\n'.format(new_alleles,
                                           linked_alleles,
                                           failed,
                                           invalid_cds,
                                           repeated,
                                           hash_collisions),
          end='\r')

    # send training file to sftp folder
    print('\nUploading Prodigal training file...')
    # add schema identifier to training file name
    # same training file can be uploaded more than once with different name
    # but this redundancy simplifies the process of maintaining training files
#    ptf_nsid = '{0}_{1}.trn'.format(ptf_basename.rstrip('.trn'), schema_id)
#    aux.upload_sftp(cnst.HOST_NS,
#                    user.split('@')[0],
#                    password,
#                    os.path.join(input_files, ptf_basename),
#                    os.path.join(cnst.SFTP_PTF_TEMP, ptf_nsid))

    # unlock schema
    schema_unlock = aux.simple_post_request(base_url, headers_post,
                                            ['species', species_id,
                                             'schemas', schema_id, 'lock'],
                                            {'action': 'unlock'})

    end = time.time()
    delta = end - start

    # determine elapsed time in minutes
    minutes = int(delta / 60)
    seconds = int(delta % 60)
    print('\nElapsed time: {0}m{1}s'.format(minutes, seconds))


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2],
         args[3], args[4], args[5], args[6])
