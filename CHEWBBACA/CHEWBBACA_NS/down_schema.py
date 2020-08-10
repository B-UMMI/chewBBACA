#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the download of chewBBACA's schemas from the
Chewie-NS.

The process enables the download of ZIP archives that contain ready-to-use
versions of any schema in the Chewie-NS. It also allows users to download
any schema with the structure it had at a specific time point. It is also
possible to download the latest version of the schema through requests to
the Chewie-NS API, if the compressed version that is available does not
match the latest version of the schema. An alternative approach that can
be applied to get the latest version of the schema, if the compressed version
does not provide it, is to download the compressed version that is available
and run the SyncSchema process to retrieve the alleles that were added to the
schema after the compression date.

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-sp``, ``species_id`` : The integer identifier or name of the species
  that the schema will be associated to in the Chewie-NS.

    - e.g.: ``1`` or ``'Yersinia pestis'``

- ``-sc``, ``schema_id`` : The schema identifier in the Chewie-NS.

    - e.g.: ``1``

- ``-o``, ``download_folder`` : Path to the parent directory of the folder
  that will store the downloaded schema. The process will create a folder
  with the schema's name inside the directory specified through this argument.

    - e.g.: ``/home/user/chewie_schemas``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores that will be used to
  construct the schema if the process downloads FASTA files instead of
  the compressed version.

    - e.g.: ``4``

- ``--ns_url``, ``nomenclature_server_url`` : The base URL for the Nomenclature
  Server. The default value, "main", will establish a connection to
  "https://chewbbaca.online/", "tutorial" to "https://tutorial.chewbbaca.online/"
  and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). Users may also
  provide the IP address to other Chewie-NS instances.

    - e.g.: ``http://127.0.0.1:5000/NS/api/`` (local host)

- ``--d``, ``date`` : Download schema with structure it had at
  specified date. Must be in the format "Y-m-dTH:M:S" or "Y-m-dTH:M:S.f".

    - e.g.: ``2020-03-27T11:38:00`` or ``2020-03-27T11:38:01.100``

- ``--latest`` : If the compressed version that is available is not the
  latest, downloads all loci and constructs schema locally.

Code documentation
------------------
"""


import os
import sys
import pickle
import shutil
import requests
import argparse
import datetime as dt
import concurrent.futures
from itertools import repeat
from urllib3.exceptions import InsecureRequestWarning

try:
    from utils import constants as cnst
    from utils import auxiliary_functions as aux
    from utils import parameters_validation as pv
    from PrepExternalSchema import PrepExternalSchema
except:
    from CHEWBBACA.utils import constants as cnst
    from CHEWBBACA.utils import auxiliary_functions as aux
    from CHEWBBACA.utils import parameters_validation as pv
    from CHEWBBACA.PrepExternalSchema import PrepExternalSchema


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


def check_compressed(schema_uri, headers_get):
    """ Determines if there is a compressed version of
        a schema.

        Parameters
        ----------
        schema_uri : str
            The URI of the schema in the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        list
            A list with the following elements:

            - The URI for the compressed version of the schema (str).
            - The timestamp of the compressed version. Indicates the
              last modification date of the schema at time of compression.
    """

    zip_uri = '{0}/zip'.format(schema_uri)
    zip_response = requests.get(zip_uri, headers=headers_get,
                                params={'request_type': 'check'}, verify=False)
    zip_info = zip_response.json()
    if 'zip' in zip_info:
        zip_file = zip_info['zip'][0]
        zip_date = zip_file.split('_')[-1].split('.zip')[0]
    else:
        zip_date = None

    return [zip_uri, zip_date]


def download_date(user_date, zip_date, latest, insertion_date,
                  modification_date):
    """ Determines the date that will be used as timepoint for
        schema download.

        Parameters
        ----------
        user_date : str or None
            Date value that was received by :py:func:`main`.
            If the user provided a value, it should be a date
            in the format %Y-%m-%dT%H:%M:%S or %Y-%m-%dT%H:%M:%S.%f.
            If the user did not provide a value for the argument, it
            will be None by default.
        zip_date : str or None
            Date value of the compressed version of the schema that
            is available for download. Will be in the format
            %Y-%m-%dT%H:%M:%S or %Y-%m-%dT%H:%M:%S.%f if there is
            a compressed version. Will be None if there is no compressed
            version.
        latest : bool
            Boolean value that indicates if the user wants to download
            the latest version of the schema.
        insertion_date : str
            Date on which the initial schema upload/insertion was
            completed.
        modification_date : str
            Date on which the schema was last modified.

        Returns
        -------
        schema_date : str
            Date selected as timepoint. Schema will be downloaded with
            the data structure it had on this date.

        Raises
        ------
        SystemExit
            - If the value received for the `user_date` is not a
              valid date in the format %Y-%m-%dT%H:%M:%S or
              %Y-%m-%dT%H:%M:%S.%f.
            - If the value received for the `user_date` is a date
              that is prior to `insertion_date` or greater than the
              `modification_date`.
    """

    # user did not provide a value for 'date' or 'latest' arguments
    if user_date is None and latest is False:
        if zip_date is not None:
            schema_date = zip_date
        # there is no compressed version
        else:
            schema_date = modification_date
    # user wants latest schema version
    elif latest is True:
        schema_date = modification_date
    # user wants schema at a particular time point
    elif user_date is not None:
        # get schema insertion and last modification date
        insertion_date_obj = dt.datetime.strptime(insertion_date,
                                                  '%Y-%m-%dT%H:%M:%S.%f')
        modification_date_obj = dt.datetime.strptime(modification_date,
                                                     '%Y-%m-%dT%H:%M:%S.%f')
        # determine if date given by user is valid
        user_date_obj = aux.validate_date(user_date)

        if user_date_obj is False:
            sys.exit('Provided date is invalid. Please provide a date '
                     'in the format "%Y-%m-%dT%H:%M:%S" or '
                     '"%Y-%m-%dT%H:%M:%S.%f"')
        if user_date_obj >= insertion_date_obj and \
           user_date_obj <= modification_date_obj:
            schema_date = dt.datetime.strftime(user_date_obj,
                                               '%Y-%m-%dT%H:%M:%S.%f')
        elif user_date_obj < insertion_date_obj:
            sys.exit('Provided date is prior to the date of schema '
                     'insertion. Please provide a date later than '
                     'schema insertion date ({0}).'.format(insertion_date))
        elif user_date_obj > modification_date_obj:
            sys.exit('Provided date is greater than the last '
                     'modification date. Please provide a date '
                     'equal or prior to the schema modification '
                     'date ({0}).'.format(insertion_date))

    return schema_date


def build_fasta(locus_id, locus_info, download_folder):
    """ Writes DNA sequences from a response object into
        a FASTA file.

        Parameters
        ----------
        locus_id : str
            Name of the locus in the Chewie-NS.
        locus_info : dict
            Response object with the list of DNA
            sequences of the locus.
        download_folder : str
            Path to the directory where the FASTA
            file with the locus sequences will be
            created.

        Returns
        -------
        locus_file : str
            Path to the FASTA file with the locus
            sequences.
    """

    locus_name = locus_id.rstrip('.fasta')
    locus_file = os.path.join(download_folder, locus_id+'.fasta')
    ns_data = locus_info.json()['Fasta']

    # extract identifiers and DNA sequences from response
    locus_alleles = []
    for allele in ns_data:
        allele_id = int(allele['allele_id']['value'])
        allele_seq = allele['nucSeq']['value']
        locus_alleles.append((allele_id, allele_seq))

    # create FASTA records
    records = ['>{0}_{1}\n{2}'.format(locus_name, a[0], a[1])
               for a in locus_alleles]

    # write sequences to FASTA file
    with open(locus_file, 'w') as lf:
        concat_records = '\n'.join(records)
        lf.write(concat_records)

    return locus_file


def get_fasta_seqs(url, headers_get, schema_date):
    """ Retrieves the DNA sequences of a locus in the
        Chewie-NS.

        Parameters
        ----------
        url : str
            Endpoint URL to make the request.
        headers_get : dict
            HTTP headers for GET requests.
        schema_date : str
            The function will only retrieve alleles
            that were inserted up to this date.

        Returns
        -------
        tuple
            Tuple with the following elements:
            - URI of the locus.
            - Response object with the DNA sequences
              that were downloaded.
    """

    payload = {'date': schema_date}
    tries = 0
    max_tries = 3
    downloaded = False
    while downloaded is False:
        res = requests.get(url, headers=headers_get, timeout=180,
                           params=payload, verify=False)

        tries += 1
        if res.status_code in [200, 201] or tries == max_tries:
            downloaded = True

    return (url.rstrip('/fasta'), res)


def schema_loci(schema_uri, headers_get):
    """ Retrieves the list of loci for a schema.

        Parameters
        ----------
        schema_uri : str
            The URI of the schema in the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        loci : dict
            A dictionary with loci URIs as keys and
            loci names as values.
    """

    # get the list of loci
    loci_uri = aux.make_url(schema_uri, 'loci')
    loci_res = requests.get(loci_uri, headers=headers_get, verify=False)
    loci_res = loci_res.json()['Loci']

    # locus URI to locus name
    loci = {}
    for locus in loci_res:
        loci[str(locus['locus']['value'])] = locus['name']['value']

    return loci


def download_fastas(loci, download_folder, headers_get, schema_date):
    """ Downloads and writes FASTA files for the loci of a
        schema in the Chewie-NS

        Parameters
        ----------
        loci : dict
            A dictionary with loci URIs as keys and
            loci names as values.
        download_folder : str
            Path to the directory where the FASTA files
            will be created.
        headers_get : dict
            HTTP headers for GET requests.
        schema_date : str
            The function will only retrieve alleles
            that were inserted up to this date.

        Returns
        -------
        ns_files : list
            List with the paths to the schema's FASTA
            files that were created.
    """

    # Total number of loci
    total_loci = len(loci)
    print('Number of loci to download: {0}'.format(total_loci))

    # build the list of urls to get
    fasta_urls = [aux.make_url(locus, 'fasta') for locus in loci]

    # multithread the requests
    print('Downloading schema files...')
    total = 0
    failed = []
    downloaded = 0
    ns_files = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        for result in executor.map(get_fasta_seqs, fasta_urls,
                                   repeat(headers_get), repeat(schema_date)):
            locus_id = loci[result[0]]
            locus_info = result[1]
            if locus_info.status_code in [200, 201]:
                locus_file = build_fasta(locus_id, locus_info, download_folder)
                ns_files.append(locus_file)
                downloaded += 1
            else:
                failed.append(locus_id)
            total += 1
            print('\r', 'Downloaded: '
                  '{0}/{1}'.format(downloaded, total_loci), end='')

    print('\nDownloaded and wrote FASTA files for '
          '{0}/{1} loci'.format(downloaded, total))
    print('Failed download for {0} loci.\n'.format(len(failed)))
    if len(failed) > 0:
        sys.exit('Failed download for following loci: {0}\n'
                 'Please download files for failed loci '
                 'through the API or retry full schema '
                 'download'.format(','.join(failed)))

    return ns_files


def download_compressed(zip_uri, species_name, schema_name,
                        download_folder, headers_get):
    """ Downloads and extracts a ZIP archive with a ready-to-use
        version of a schema in the Chewie-NS.

        Parameters
        ----------
        zip_uri : str
            Endpoint URL to make the request to download
            the compressed schema.
        species_name : str
            Scientific name of the schema species.
        schema_name : str
            Name of the schema in the Chewie-NS.
        download_folder : str
            Path to the directory to which the ZIP archive
            will be saved.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        schema_path : str
            ZIP archive contents will be extracted to this
            directory.
    """

    zip_name = '{0}{1}_{2}.zip'.format(species_name[0].lower(),
                                       species_name.split(' ')[-1],
                                       schema_name)
    schema_path = os.path.join(download_folder,
                               zip_name.split('.zip')[0])
    os.mkdir(schema_path)

    # download ZIP archive
    zip_response = requests.get(zip_uri, headers=headers_get,
                                params={'request_type': 'download'},
                                verify=False)
    zip_path = os.path.join(schema_path, zip_name)
    open(zip_path, 'wb').write(zip_response.content)
    # uncompress
    print('Decompressing schema...')
    shutil.unpack_archive(zip_path, extract_dir=schema_path)
    # delete ZIP
    os.remove(zip_path)

    return schema_path


def download_ptf(ptf_hash, download_folder, schema_id,
                 species_id, species_name, headers_get, base_url):
    """ Downloads the Prodigal training file for a schema.

        Parameters
        ----------
        ptf_hash : str
            Unique identifier of the Prodigal training file
            (BLAKE2 hash).
        download_folder : str
            Path to the directory to which the Prodigal
            training file should be saved.
        schema_id : str
            The identifier of the schema in the Chewie-NS.
        species_id : str
            The identifier of the schema's species in the
            Chewie-NS.
        species_name : str
            Scientific name of the schema species.
        headers_get : dict
            HTTP headers for GET requests.
        base_url : str
            Base URL of the Chewie Nomenclature server.

        Returns
        -------
        ptf_file : str
            Path to the Prodigal training file.
    """

    ptf_uri = aux.make_url(base_url, *['species', species_id,
                                       'schemas', schema_id, 'ptf'])

    ptf_response = requests.get(ptf_uri, headers=headers_get, verify=False)

    ptf_file = os.path.join(download_folder,
                            '{0}.trn'.format(species_name.replace(' ', '_')))

    open(ptf_file, 'wb').write(ptf_response.content)

    return ptf_file


def main(species_id, schema_id, download_folder, core_num,
         base_url, date, latest):

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
    print('Started at: {0}\n'.format(start_date_str))

    # GET request headers
    headers_get = cnst.HEADERS_GET_JSON

    # Get the name of the species from the provided id
    # or vice-versa
    species_info = aux.species_ids(species_id, base_url, headers_get)
    if isinstance(species_info, list):
        species_id, species_name = species_info
    else:
        sys.exit('There is no species with the provided '
                 'identifier in the Chewie-NS.')

    # check if user provided schema identifier or schema description
    # get info about all the species schemas
    schema_id, schema_uri,\
        schema_name, schema_params = aux.get_species_schemas(schema_id,
                                                             species_id,
                                                             base_url,
                                                             headers_get)

    print('Schema id: {0}'.format(schema_id))
    print('Schema name: {0}'.format(schema_name))
    print("Schema's species: {0} "
          "(id={1})".format(species_name, species_id))

    # create parameters dict
    schema_params_dict = {k: schema_params[k]['value']
                          for k in schema_params.keys()
                          if k != 'name'}

    # check if schema is locked
    lock_status = schema_params_dict['Schema_lock']
    if lock_status != 'Unlocked':
        sys.exit('Schema is locked. This might be because it '
                 'is being uploaded, updated or compressed.'
                 ' Please try again later and contact the Administrator '
                 'if the schema stays locked for a long period of time.')

    # get zip information
    zip_uri, zip_date = check_compressed(schema_uri, headers_get)

    schema_date = download_date(date, zip_date, latest,
                                schema_params_dict['dateEntered'],
                                schema_params_dict['last_modified'])

    # check output folder
    if not os.path.exists(download_folder):
        os.mkdir(download_folder)
    else:
        # verify that folder is empty and abort if it is not
        download_folder_files = os.listdir(download_folder)
        if len(download_folder_files) > 0:
            sys.exit('Download folder is not empty. Please ensure '
                     'that folder is empty to guarantee proper '
                     'schema creation or provide a valid path for '
                     'a new folder that will be created.')

    if schema_date == zip_date:
        print('\nDownloading compressed version...')
        schema_path = download_compressed(zip_uri, species_name, schema_name,
                                          download_folder, headers_get)
    else:
        print('\nDownloading schema FASTA files...')
        # download FASTA files
        loci = schema_loci(schema_uri, headers_get)
        ns_files = download_fastas(loci, download_folder, headers_get,
                                   schema_date)

        # download Prodigal training file
        ptf_hash = schema_params_dict['prodigal_training_file']
        ptf_file = download_ptf(ptf_hash, download_folder, schema_id,
                                species_id, species_name, headers_get,
                                base_url)

        # use PrepExternalSchema main to determine representatives
        genus, epithet = species_name.split(' ')
        schema_name = '{0}{1}_{2}'.format(genus[0].lower(), epithet, schema_name)
        schema_path = os.path.join(download_folder, schema_name)

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
            os.remove(file)

        # write hidden schema config file
        del(schema_params_dict['Schema_lock'])
        schema_config = aux.write_schema_config(schema_params_dict['bsr'],
                                                ptf_hash,
                                                schema_params_dict['translation_table'],
                                                schema_params_dict['minimum_locus_length'],
                                                schema_params_dict['chewBBACA_version'],
                                                schema_params_dict['size_threshold'],
                                                schema_path)

        # create ns_config file
        ns_config = os.path.join(schema_path, '.ns_config')
        if not os.path.exists(ns_config):
            with open(ns_config, 'wb') as nc:
                download_info = [schema_date, schema_uri]
                pickle.dump(download_info, nc)

        genes_list_file = aux.write_gene_list(schema_path)

    print('Schema is now available at: {0}'.format(schema_path))

    end_date = dt.datetime.now()
    end_date_str = dt.datetime.strftime(end_date, '%Y-%m-%dT%H:%M:%S')

    delta = end_date - start_date
    minutes, seconds = divmod(delta.total_seconds(), 60)

    print('\nFinished at: {0}'.format(end_date_str))
    print('Elapsed time: {0:.0f}m{1:.0f}s'.format(minutes, seconds))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-sp', type=str, required=True,
                        dest='species_id',
                        help='The integer identifier or name of the '
                             'species that the schema is associated '
                             'to in the Chewie-NS.')

    parser.add_argument('-sc', type=str, required=True,
                        dest='schema_id',
                        help='The schema identifier in the Chewie-NS.')

    parser.add_argument('-o', type=str, required=True,
                        dest='download_folder',
                        help='Path to the parent directory of the '
                             'folder that will store the downloaded schema.')

    parser.add_argument('--cpu', type=int, required=False,
                        dest='cpu_cores', default=1,
                        help='Number of CPU cores that will '
                             'be used to construct the schema '
                             'if the process downloads FASTA '
                             'files instead of a compressed version.')

    parser.add_argument('--ns', type=pv.validate_ns_url, required=False,
                        dest='nomenclature_server',
                        default='main',
                        help='The base URL for the Nomenclature Server. '
                             'The default value, "main", will establish a '
                             'connection to "https://chewbbaca.online/", '
                             '"tutorial" to "https://tutorial.chewbbaca.online/" '
                             'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
                             'Users may also provide the IP address to other '
                             'Chewie-NS instances.')

    parser.add_argument('--d', type=str, required=False,
                        default=None,
                        dest='date',
                        help='Download schema with structure it had at '
                             'specified date. Must be in the format '
                             '"Y-m-dTH:M:S" or "Y-m-dTH:M:S.f".')

    parser.add_argument('--latest', required=False,
                        action='store_true', dest='latest',
                        help='If the compressed version that is available '
                             'is not the latest, downloads all loci and '
                             'constructs schema locally.')

    args = parser.parse_args()

    return [args.species_id, args.schema_id, args.download_folder,
            args.cpu_cores, args.nomenclature_server, args.date,
            args.latest]


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3],
         args[4], args[5], args[6])
