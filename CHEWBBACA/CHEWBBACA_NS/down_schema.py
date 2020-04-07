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
import sys
import time
import pickle
import shutil
import requests
import argparse
import datetime as dt
import concurrent.futures
from getpass import getpass
from itertools import repeat

from utils import constants as cnst
from utils import auxiliary_functions as aux
from PrepExternalSchema import PrepExternalSchema


def build_fasta(locus_id, locus_info, download_folder):
    """ Write fasta files from get request responses

        Args:
            response_dict (dict): Contains the names as keys,
            get request response as values.

            path2down (str): Path to the download directory

        Results:
            failed_downloads (list): Contains the names.
    """

    # if the fasta already exists, stop beep boop
    locus_file = os.path.join(download_folder, locus_id)
    if locus_info.status_code > 200:
        result = [locus_id, 1]
    else:
        locus_name = locus_id.rstrip('.fasta')
        ns_data = locus_info.json()['Fasta']

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

        locus_file = os.path.join(download_folder, locus_id)
        with open(locus_file, 'w') as lf:
            concat_records = '\n'.join(records)
            lf.write(concat_records)

        result = [locus_id, 0]

    return result


def get_fasta_seqs(headers_get, schema_date, url):
    """ Perform get requests to the NS

        Args:
            headers_get (dict): headers for the GET request
            url (str): url to perform the request

        Return:
            res (response): response of the request containing
            fasta sequences
    """

    payload = {'date':schema_date}
    res = requests.get(url, headers=headers_get, timeout=30, params=payload)

    return (url.rstrip('/fasta'), res)


def get_schema(schema_uri, download_folder, headers_get, schema_desc, schema_date):
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

    # Total number of loci
    total_loci = len(loci_list)
    loci_names = {}
    print(f"Number of loci to download: {str(total_loci)}")
    for locus in loci_list:
        loci_names[str(locus['locus']['value'])] = f"{locus['name']['value']}.fasta"

    # build the list of urls to get
    fasta_urls = [aux.make_url(locus, 'fasta') for locus in loci_names]

    # multithread the requests
    print('Downloading schema files...')
    total = 0
    failed = 0
    downloaded = 0
    ns_files = []
    total_files = len(fasta_urls)
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        for result in executor.map(get_fasta_seqs, repeat(headers_get), repeat(schema_date), fasta_urls):
            locus_id = loci_names[result[0]]
            locus_info = result[1]
            locus_file = build_fasta(locus_id, locus_info, download_folder)
            if locus_file[1] == 0:
                downloaded += 1
                ns_files.append(locus_id)
            else:
                failed += 1
            total += 1
            print('\r', 'Downloaded: '
                  '{0}/{1}'.format(downloaded, total_files), end='')

    print('\nDownloaded and wrote FASTA files for '
          '{0}/{1} loci'.format(downloaded, total))
    print('Failed download for {0} loci.\n'.format(failed))

    # output dir for the result of PrepExternalSchema
    local_schema_name = '{0}_schema'.format(schema_desc)
    local_schema_path = os.path.join(download_folder, local_schema_name)

    return [local_schema_path, ns_files]


def main(schema_id, species_id, download_folder, core_num, base_url, date, latest):

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
            sys.exit(message)

    # Define the headers of the requests
    headers_get = {'Authorization': token,
                   'accept': 'application/json'}

    # Get the name of the species from the provided id
    # or vice-versa
    species_info = aux.species_ids(species_id, base_url, headers_get)
    if isinstance(species_info, list):
        species_id, species_name = species_info
        print('NS species with identifier {0} is {1}.'.format(species_id,
                                                                species_name))
    else:
        sys.exit('There is no species with the provided identifier in the NS.')

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
            schema_id = schema_uri.split('/')[-1]
    # check if user provided schema identifier or schema description
    elif len(schema_uri) == 0:
        # get info about all the species schemas
        schema_get = aux.simple_get_request(base_url, headers_get,
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

    # get schema parameters
    schema_params = requests.get(schema_uri, headers=headers_get)
    schema_params = schema_params.json()[0]

    # create parameters dict
    schema_params_dict = {k: schema_params[k]['value']
                          for k in schema_params.keys()
                          if k != 'name'}
    schema_desc = schema_params['name']['value']

    # check if schema is locked
    lock_status = schema_params_dict['Schema_lock']
    if lock_status != 'Unlocked':
        sys.exit('403: schema is locked. This might be because it '
                 'is being uploaded, updated or compressed.'
                 ' Please try again later and contact the Administrator '
                 'if the schema stays locked for a long period of time.')

    # get zip information
    zip_uri = '{0}/zip'.format(schema_uri)
    zip_response = requests.get(zip_uri, headers=headers_get, params={'request_type': 'check'})
    zip_info = zip_response.json()
    if 'zip' in zip_info:
        zip_file = zip_info['zip'][0]
        zip_date = zip_file.split('_')[-1].split('.zip')[0]
    else:
        zip_date = None

    if date is not None and latest is not True:
        insertion_date = schema_params_dict['dateEntered']
        modification_date = schema_params_dict['last_modified']
        insertion_date_obj = dt.datetime.strptime(insertion_date, '%Y-%m-%dT%H:%M:%S.%f')
        modification_date_obj = dt.datetime.strptime(modification_date, '%Y-%m-%dT%H:%M:%S.%f')
        user_date = dt.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S')
        if user_date >= insertion_date_obj and user_date <= modification_date_obj:
            schema_date = date
        elif user_date < insertion_date_obj:
            sys.exit('Provided date is prior to the date of schema insertion. '
                     'Please provide a date later than schema insertion date ({0}).'.format(insertion_date))
        elif user_date > modification_date_obj:
            sys.exit('Provided date is greater than the last modification date. '
                     'Please provide a date equal or prior to the schema modification date ({0}).'.format(insertion_date))
    elif latest is True:
        schema_date = schema_params_dict['last_modified']
    else:
        schema_date = None

    if schema_date is None and zip_date is not None:
        schema_date = zip_date
    elif schema_date is None and zip_date is None:
        schema_date = schema_params_dict['last_modified']

    # create download folder if it does not exist
    if not os.path.exists(download_folder):
        os.mkdir(download_folder)
    else:
        # verify that folder is empty and abort if it is not
        print('Download folder already exists...')
        download_folder_files = os.listdir(download_folder)
        if len(download_folder_files) > 0:
            sys.exit('Download folder is not empty. '
                     'Please ensure that folder is empty to guarantee proper '
                     'schema creation or provide a valid path for a new folder '
                     'that will be created.')

    print('\nSchema info:\n  ID: {0}\n  URI: {1}\n  Species: {2}\n  Download directory: {3}'.format(schema_id,
                                                                                                    schema_uri,
                                                                                                    species_name,
                                                                                                    download_folder))

    start = time.time()

    # download compressed version if date matches compression version
    if schema_date == zip_date:
        print('\nDownloading schema archive...')
        zip_name = '{0}{1}_{2}.zip'.format(species_name[0].lower(), species_name.split(' ')[-1], schema_desc)
        schema_path = os.path.join(download_folder, zip_name.split('.zip')[0])
        os.mkdir(schema_path)
        zip_response = requests.get(zip_uri, headers=headers_get, params={'request_type': 'download'})
        zip_path = os.path.join(schema_path, zip_name)
        open(zip_path, 'wb').write(zip_response.content)
        shutil.unpack_archive(zip_path, extract_dir=schema_path)
        os.remove(zip_path)
    # user wants schema at specific time point or latest
    # version that has yet to be compressed
    else:
        print('\nDownloading loci files for local construction...')
        schema_path, ns_files = get_schema(schema_uri,
                                           download_folder,
                                           headers_get,
                                           schema_file_desc,
                                           schema_date)

        # download Prodigal training file
        ptf_hash = schema_params_dict['prodigal_training_file']

        #payload = {'species_id':species_id, 'schema_id': schema_id}

        ptf_uri = aux.make_url(base_url, *['species', species_id, 'schemas', schema_id, 'ptf'])

        ptf_response = requests.get(ptf_uri, headers=headers_get)#, params=payload)

        ptf_file = os.path.join(download_folder, '{0}.trn'.format('_'.join(species_name.split(' '))))

        open(ptf_file ,'wb').write(ptf_response.content)

        # determine representatives and create schema
        PrepExternalSchema.main(download_folder,
                                schema_path,
                                core_num,
                                float(schema_params_dict['bsr']),
                                int(schema_params_dict['minimum_locus_length']),
                                int(schema_params_dict['translation_table']),
                                ptf_file,
                                None)

        # copy Prodigal training file to schema directory
        shutil.copy(ptf_file, schema_path)
        os.remove(ptf_file)

        # remove FASTA files with sequences from the NS
        for file in ns_files:
            os.remove(os.path.join(download_folder, file))

        # write hidden schema config file
        del(schema_params_dict['Schema_lock'])
        schema_config = os.path.join(schema_path, '.schema_config')
        with open(schema_config, 'wb') as scf:
            pickle.dump(schema_params_dict, scf)

        # create ns_config file
        ns_config = os.path.join(schema_path, '.ns_config')
        if not os.path.exists(ns_config):
            with open(ns_config, 'wb') as nc:
                download_info = [schema_date, schema_uri]
                pickle.dump(download_info, nc)

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

    parser.add_argument('--d', type=str, required=False,
                        default=None,
                        dest='date',
                        help='Download schema with state from specified date. '
                             'Must be in the format "Y-m-dTH:M:S".')

    parser.add_argument('--latest', required=False, action='store_true',
                        dest='latest',
                        help='If the compressed version that is available is '
                             'not the latest, downloads all loci and constructs '
                             'schema locally.')

    args = parser.parse_args()

    return [args.schema_id, args.species_id, args.download_folder,
            args.core_num, args.ns_url, args.date, args.latest]


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2],
         args[3], args[4], args[5],
         args[6])
