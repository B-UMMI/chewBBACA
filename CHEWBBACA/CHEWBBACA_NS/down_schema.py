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
import re
import time
import pickle
import shutil
import requests
import argparse
import urllib.request
import concurrent.futures
from getpass import getpass
from itertools import repeat

from utils import auxiliary_functions as aux
from PrepExternalSchema import PrepExternalSchema


def species_list(base_url, headers_get, endpoint_list):
    """
    """

    res = simple_get_request(base_url, headers_get, endpoint_list)
    res = res.json()
    species_lst = {}
    for sp in res:
        species = sp['name']['value']
        species_url = sp['species']['value']
        species_id = species_url.split('/')[-1]

        species_lst[species] = species_id

    return species_lst


def simple_get_request(base_url, headers, endpoint_list):
    """ Constructs an endpoint URI and uses a GET method to retrive
        information from the endpoint.

        Args:
            base_url (str): the base URI for the NS, used to concatenate
            with a list of elements and obtain endpoints URL.
            headers (dict): headers for the GET method used to
            get data from the API endpoints.
            endpoint_list (list): list with elements that will be
            concatenated to the base URL to obtain the URL for
            the API endpoint.
        Returns:
            res (requests.models.Response): response object from
            the GET method.
    """

    # unpack list of sequential endpoints and pass to create URI
    url = aux.make_url(base_url, *endpoint_list)
    res = requests.get(url, headers=headers, timeout=30)

    return res


def species_ids(species_id, base_url, headers_get):
    """
    """

    try:
        int(species_id)
        species_info = simple_get_request(base_url, headers_get,
                                          ['species', species_id])
        if species_info.status_code == 200:
            species_name = species_info.json()[0]['name']['value']
            return [species_id, species_name]
        else:
            return 404
    except ValueError:
        species_name = species_id
        ns_species = species_list(base_url, headers_get, ['species', 'list'])
        species_id = ns_species.get(species_name, 'not_found')
        if species_id != 'not_found':
            return [species_id, species_name]
        else:
            return 404


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


def build_fasta_files(response_dict, download_folder):
    """ Write fasta files from get request responses

        Args:
            response_dict (dict): Contains the names as keys,
            get request response as values.

            path2down (str): Path to the download directory

        Results:
            failed_downloads (list): Contains the names.
    """

    failed_downloads = []
    for locus, response in response_dict.items():

        # if the fasta already exists, stop beep boop
        locus_file = os.path.join(download_folder, locus)
        if os.path.exists(locus_file):
            continue

        if response.status_code > 200:
            failed_downloads.append([locus, response])
            continue
        else:
            locus_name = locus.rstrip('.fasta')
            ns_data = response.json()['Fasta']

            # allele identifier to DNA sequence
            locus_alleles = []
            for allele in ns_data:
                allele_id = int(allele['allele_id']['value'])
                allele_seq = f"{allele['nucSeq']['value']}"
                locus_alleles.append((allele_id, allele_seq))

            # write sequences to FASTA file
            records = []
            for allele in locus_alleles:
                record = '>{0}_{1}\n{2}'.format(locus_name,
                                                allele[0],
                                                allele[1])
                records.append(record)

            locus_file = os.path.join(download_folder, locus)
            with open(locus_file, 'w') as lf:
                concat_records = '\n'.join(records)
                lf.write(concat_records)

    return failed_downloads


def get_fasta_seqs(headers_get, url):
    """ Perform get requests to the NS

        Args:
            headers_get (dict): headers for the GET request
            url (str): url to perform the request

        Return:
            res (response): response of the request containing
            fasta sequences
    """

    res = requests.get(url, headers=headers_get, timeout=30)

    return (url.rstrip('/fasta'), res)


def get_schema(schema_uri, download_folder, headers_get, schema_desc):
    """ Downloads, builds a writes a schema from NS

        Args:
            schema_uri (str): uri of the schema to download
            path2down (str): Path to the download directory
            cpu2use (int): Number of cpu allowed
            maxBsrShort (float): Maximum BSR allowed for the
            representative selection
            headers_get (dict): headers for the GET request
        Returns:
    """

    schema_loci_uri = aux.make_url(schema_uri, 'loci')

    # get list of loci and build dictionary locus_id --> gene_name
    loci_res = requests.get(schema_loci_uri, headers=headers_get)
    loci_list = loci_res.json()
    loci_list = loci_list['Loci']
    server_time = loci_res.headers['Server-Date']

    # get server time and save on config before starting to download
    # useful for further sync function
    ns_config = os.path.join(download_folder, '.ns_config')
    if not os.path.exists(ns_config):
        with open(ns_config, 'wb') as nc:
            download_info = [server_time, schema_uri]
            pickle.dump(download_info, nc)

    # Total number of loci
    total_loci = len(loci_list)
    loci_names = {}
    print(f"Number of loci to download: {str(total_loci)}")
    for locus in loci_list:
        loci_names[str(locus['locus']['value'])] = f"{locus['name']['value']}.fasta"

    # build the list of urls to get
    fasta_urls = [aux.make_url(locus, 'fasta') for locus in loci_names]

    # multithread the requests
    # must write files or this will go into RAM and explode
    print('Downloading schema files...\n')
    responses = {}
    downloaded = 0
    total_files = len(fasta_urls)
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        for result in executor.map(get_fasta_seqs, repeat(headers_get), fasta_urls):
            responses[result[0]] = result[1]
            downloaded += 1
            print('\r', 'Downloaded: '
                  '{0}/{1}'.format(downloaded, total_files), end='')

    # Build a dict with fasta_file ----> response
    res_dict = {loci_names[locus]: responses[locus] for locus in responses}

    # Build the fasta files
    print('\nWriting fasta files...\n')
    failed_loci = build_fasta_files(res_dict, download_folder)
    successful_loci = len(res_dict) - len(failed_loci)
    print('Downloaded and wrote FASTA files for '
          '{0}/{1} loci'.format(successful_loci, len(res_dict)))
    print('Failed for {0} loci.'.format(len(failed_loci)))

    # output dir for the result of PrepExternalSchema
    local_schema_name = '{0}_schema'.format(schema_desc)
    local_schema_path = os.path.join(download_folder, local_schema_name)

    ns_files = list(res_dict.keys())

    return [local_schema_path, ns_files, ns_config]


def main(schema_id, species_id, download_folder, core_num, base_url):

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

    # Define the headers of the requests
    headers_get = {'Authorization': token,
                   'accept': 'application/json'}

    # Get the name of the species from the provided id
    # or vice-versa
    species_info = species_ids(species_id, base_url, headers_get)
    if isinstance(species_info, list):
        species_id, species_name = species_info
        print('\nNS species with identifier {0} is {1}.'.format(species_id,
                                                                species_name))
    else:
        print('\nThere is no species with the provided identifier in the NS.')
        return 1

    print('\nChecking if schema exists...')

    # check if schema identifier provided by user is the schema URI
    schema_uri = re.findall('{0}species/{1}/schemas/[0-9]*'.format(base_url,
                                                                   species_id),
                            schema_id)

    if len(schema_uri) > 0:
        schema_response = requests.get(schema_uri[0],
                                       headers=headers_get,
                                       timeout=5)
        if schema_response.status_code > 201:
            print('There is no schema with URI "{0}" for '
                  '{1}.'.format(species_name))
        else:
            schema_uri = schema_uri[0]
            schema_id = schema_uri.split('/')[0]
    # check if user provided schema identifier or schema description
    elif len(schema_uri) == 0:
        # get info about all the species schemas
        schema_get = simple_get_request(base_url, headers_get,
                                        ['species', species_id, 'schemas'])
        schema_get_status = schema_get.status_code
        if schema_get_status in [200, 201]:
            species_schemas = schema_get.json()

            schemas_info = {}
            for s in species_schemas:
                suri = s['schemas']['value']
                sid = suri.split('/')[-1]
                sdesc = s['name']['value']
                schemas_info[sid] = [suri, sdesc]

            # determine if schema identifier provided
            # by user is an integer identifier or a description
            found = False
            for k, v in schemas_info.items():
                if schema_id == k or schema_id == v[1]:
                    schema_id = k
                    schema_uri = v[0]
                    schema_desc = v[1]
                    schema_file_desc = '{0}{1}_{2}'.format(species_name[0].lower(),
                                                           species_name.split(' ')[1],
                                                           schema_desc)
                    found = True

            if found is False:
                print('\nCould not find a schema with such description.')
                return 404
        else:
            print('\nCould not retrieve schemas for current species.')

    # create download folder if it does not exist
    if not os.path.exists(download_folder):
        os.mkdir(download_folder)
    else:
        # verify that folder is empty and abort if it is not
        print('Download folder already exists...')
        download_folder_files = os.listdir(download_folder)
        if len(download_folder_files) > 0:
            print('Download folder is not empty.\n'
                  'Please ensure that folder is empty to guarantee proper\n'
                  'schema creation or provide a valid path for a new folder '
                  'that will be created.')
            return 1

    print('Schema info:\nID: {0}\nURI: {1}\nSpecies: {2}\nDownload directory: {3}'.format(schema_id,
                                                                                          schema_uri,
                                                                                          species_name,
                                                                                          download_folder))

    # get schema parameters
    schema_params = requests.get(schema_uri, headers=headers_get)
    schema_params = schema_params.json()[0]

    # create parameters dict
    schema_params_dict = {k: schema_params[k]['value']
                          for k in schema_params.keys()
                          if k != 'name'}

    print('Downloading schema...')

    start = time.time()

    schema_path, ns_files, ns_config = get_schema(schema_uri,
                                                  download_folder,
                                                  headers_get,
                                                  schema_file_desc)

    # determine representatives and create schema
    PrepExternalSchema.main(download_folder,
                            schema_path,
                            core_num,
                            float(schema_params_dict['bsr']),
                            int(schema_params_dict['minimum_locus_length']),
                            int(schema_params_dict['translation_table']))

    # copy ns_config file
    shutil.copy(ns_config, schema_path)
    if os.path.isfile(ns_config):
        os.remove(ns_config)

    # remove FASTA files with sequences from the NS
    for file in ns_files:
        os.remove(os.path.join(download_folder, file))

    # write hidden schema config file
    schema_config = os.path.join(schema_path, '.schema_config')
    with open(schema_config, 'wb') as scf:
        pickle.dump(schema_params_dict, scf)

    # Download prodigal training file
#    ptf_url = schema_params_dict['prodigal_training_file']
#    ptf_name = ptf_url.split('/')[-1]
#    local_ptf = os.path.join(schema_path, ptf_name)
#    urllib.request.urlretrieve(ptf_url, local_ptf)

    print('Schema is now available at: {0}'.format(schema_path))

    end = time.time()
    delta = end - start

    # determine elapsed time in minutes
    minutes = int(delta / 60)
    seconds = int(delta % 60)
    print('\nElapsed time: {0}m{1}s'.format(minutes, seconds))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-schema', type=str, dest='schema_id', required=True,
                        help='')

    parser.add_argument('-species', type=str, dest='species_id', required=True,
                        help='The identifier for the schemas species in '
                        'the NS. Can be the species name or the integer '
                        'identifier for that species in the NS.')

    parser.add_argument('-out', type=str, required=True,
                        dest='download_folder',
                        help='')

    parser.add_argument('--cores', type=int, required=False, dest='core_num',
                        default=1,
                        help='')

    parser.add_argument('--ns_url', type=str, required=False, dest='ns_url',
                        default='http://127.0.0.1:5000/NS/api/',
                        help='')

    args = parser.parse_args()

    return [args.schema_id, args.species_id, args.download_folder,
            args.core_num, args.ns_url]


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2],
         args[3], args[4])
