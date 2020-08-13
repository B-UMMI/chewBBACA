#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the retrieval of information/stats from the
Chewie-NS. Its main objective is to provide information about
the list of species and schemas in the Chewie-NS, so that users
can quickly identify a schema of interest and download it (this
process generates tables with species and schemas identifiers that
can be passed to the `-sc` and `-sp` arguments of DownloadSchema).

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-m``, ``stats_mode`` : The process can retrieve the list of species
  ("species" option) in the Chewie-NS, the list of schemas for a species
  ("schemas" option and valid value for `--sp`) or information about a
  single schema ("schemas" option and valid values for `--sp` and `--sc`).

    - e.g.: ``species`` or ``schemas``

- ``--ns_url``, ``nomenclature_server_url`` : The base URL for the Nomenclature Server.
  The default value, "main", will establish a connection to "https://chewbbaca.online/",
  "tutorial" to "https://tutorial.chewbbaca.online/"" and "local" to
  "http://127.0.0.1:5000/NS/api/" (localhost). Users may also provide the IP address to
  other Chewie-NS instances.

    - e.g.: ``http://127.0.0.1:5000/NS/api/`` (localhost)

- ``--sp``, ``species_id`` : The integer identifier of a species
  in the Chewie-NS. The process will retrieve the list of schemas
  for the species with specified identifier.

    - e.g.: ``2``

- ``--sc``, ``schema_id`` : The integer identifier of a schema in
  the Chewie-NS. The process will retrieve information about the
  schema with specified identifier.

    - e.g.: ``4``

Code documentation
------------------
"""


import sys
import requests
import argparse
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


def species_stats(base_url, headers_get):
    """ Retrieves the list of species in the Chewie-NS
        and the total number of schemas, loci and alleles
        per species.

        Parameters
        ----------
        base_url : str
            Base URL of the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        species_stats : list of str
            A list with formatted strings that can be joined
            and printed to display a table with the information
            retrieved from the Chewie-NS.
    """

    # get species and number of schemas
    species = species_schemas_count(base_url, headers_get)

    species_stats = []
    # table header
    species_stats.append('-'*78)
    species_stats.append('{:<30}  {:^10}  {:^10}  {:^10}  '
                         '{:^10}'.format('Species', 'id', '#schemas',
                                         '#loci', '#alleles'))
    species_stats.append('-'*78)

    for sp in species:
        # get totla number of loci and alleles per species
        stats = schema_stats(sp[0], base_url, headers_get)
        total_loci = 0
        total_alleles = 0
        for s in stats:
            total_loci += int(s['nr_loci'])
            total_alleles += int(s['nr_alleles'])

        species_stats.append('{:<30}  {:^10}  {:^10}  {:^10}  '
                             '{:^10}'.format(sp[1], sp[0], sp[2],
                                             total_loci, total_alleles))
    species_stats.append('-'*78)

    return species_stats


def schema_stats(species_id, base_url, headers_get):
    """ Retrieves schema properties, number of loci and
        number of alleles for all schemas of a species in
        the Chewie-NS.

        Parameters
        ----------
        species_id : str
            The integer identifier of the species in
            the Chewie-NS.
        base_url : str
            Base URL of the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        res : list of dict or None
            List with one dict per schema ot NoneType
            if it was not possible to retrieve information.
    """

    endpoint_list = ['stats', 'species', species_id, 'totals']
    # unpack list of sequential endpoints and pass to create URI
    res = aux.simple_get_request(base_url, headers_get, endpoint_list)
    status_code = res.status_code
    if status_code not in [200, 201]:
        res = None
    else:
        res = res.json()['message']

    return res


def species_schemas_count(base_url, headers_get):
    """ Returns the number of schemas per species in
        the Chewie-NS.

        Parameters
        ----------
        base_url : str
            Base URL of the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        info : list of list
            A list with a sublist per species.
            Each sublist contains the species
            identifier, the name of the species
            and the total number of schemas.
    """

    endpoint_list = ['stats', 'species']
    # unpack list of sequential endpoints and pass to create URI
    res = aux.simple_get_request(base_url, headers_get, endpoint_list)
    res = res.json()

    if 'message' in res:
        res = res['message']
    else:
        sys.exit('Could not retrieve species info.')

    if len(res) == 0:
        sys.exit('Could not retrieve species info.')
    else:
        info = []
        for s in res:
            sid = s['species']['value'].split('/')[-1]
            name = s['name']['value']
            num_schemas = s['schemas']['value']
            info.append([sid, name, num_schemas])

    # sort by species identifier
    info = sorted(info, key=lambda x: int(x[0]))

    return info


def species_schemas(species_id, base_url, headers_get):
    """ Retrieves the species and all the schemas for that
        species.

        Parameters
        ----------
        species_id : str
            The integer identifier of the species in
            the Chewie-NS.
        base_url : str
            Base URL of the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        res : list of dict
            The first dictionary contains the species
            URI and name and the following dictionaries
            contain the URI and name for all schemas
            associated with the species.
    """

    endpoint_list = ['species', species_id]
    # unpack list of sequential endpoints and pass to create URI
    res = aux.simple_get_request(base_url, headers_get, endpoint_list)
    res = res.json()

    return res


def single_species(species_id, base_url, headers_get):
    """ Retrieves the list of schemas for a species in
        the Chewie-NS and the total number of loci and
        alleles per schema.

        Parameters
        ----------
        species_id : str
            The integer identifier of the species in
            the Chewie-NS.
        base_url : str
            Base URL of the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        schemas_stats : list of str
            A list with formatted strings that can be joined
            and printed to display a table with the information
            retrieved from the Chewie-NS.

        Raises
        ------
        SystemExit
            - If the process cannot retrieve a species with
              specified identifier.
            - If the process could not retrieve schemas for the
              species.
    """

    # get all schemas for species
    schemas = species_schemas(species_id, base_url, headers_get)
    if 'NOT FOUND' in schemas:
        sys.exit('Could not find information about a species with '
                 'provided identifier (id={0}).'.format(species_id))
    species_name = schemas[0]['name']['value']
    schemas = {s['schemas']['value'].split('/')[-1]: s['schemaName']['value']
               for s in schemas[1:]}

    # get more info for all schema
    stats = schema_stats(species_id, base_url, headers_get)
    if stats is None:
        sys.exit('Could not retrieve schemas for {0} '
                 '(id={1}).'.format(species_name, species_id))
    schemas_stats = []
    schemas_stats.append('{0} (id={1})'.format(species_name, species_id))
    schemas_stats.append('-'*66)
    schemas_stats.append('{:<30}  {:^10}  {:^10}  {:^10}'
                         ''.format('Schema_name', 'id', '#loci', '#alleles'))
    schemas_stats.append('-'*66)
    for s in stats:
        schema_id = s['uri'].split('/')[-1]
        total_loci = s['nr_loci']
        total_alleles = s['nr_alleles']
        schema_name = schemas[schema_id]
        schemas_stats.append('{:<30}  {:^10}  {:^10}  {:^10}'
                             ''.format(schema_name, schema_id,
                                       total_loci, total_alleles))
    schemas_stats.append('-'*66)

    return schemas_stats


def single_schema(species_id, schema_id, base_url, headers_get):
    """ Retrieves information about a schema in the Chewie-NS.

        Parameters
        ----------
        species_id : str
            The integer identifier of the species in
            the Chewie-NS.
        schema_id : str
            The integer identifier of the schema in
            the Chewie-NS.
        base_url : str
            Base URL of the Chewie-NS.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        schema_info : list of str
            A list with formatted strings that can be joined
            and printed to display a table with the information
            retrieved from the Chewie-NS.

        Raises
        ------
        SystemExit
            - If the process cannot retrieve a species with
              specified identifier.
            - If the process cannot retrieve schemas for the
              species.
            - If the process cannot retrieve information about
              the schema with specified identifier.
    """

    # get all schemas for species
    species = species_schemas(species_id, base_url, headers_get)
    if 'NOT FOUND' in species:
        sys.exit('Could not find information about a species with '
                 'provided identifier (id={0}).'.format(species_id))
    species_name = species[0]['name']['value']

    schemas = schema_stats(species_id, base_url, headers_get)

    if schemas is None:
        sys.exit('Could not retrieve schemas for {0} '
                 '(id={1}).'.format(species_name, species_id))
    else:
        schema = [s for s in schemas if s['uri'].split('/')[-1] == schema_id]

    if len(schema) == 0:
        sys.exit('Could not find information about schema with '
                 'specified identifier (id={0}).'.format(schema_id))

    schema = schema[0]
    schema_info = []
    line_length = len(species_name)+len(schema['name'])+3
    schema_info.append('-'*line_length)
    schema_info.append('{0} - {1}'.format(species_name, schema['name']))
    schema_info.append('-'*line_length)

    schema_info.append('\nID: {0}'
                       ''.format(schema['uri'].split('/')[-1]))
    schema_info.append('Created by: {0}'
                       ''.format(schema['user']))
    schema_info.append('Total loci: {0}'
                       ''.format(schema['nr_loci']))
    schema_info.append('Total alleles: {0}'
                       ''.format(schema['nr_alleles']))
    schema_info.append('BLAST Score Ratio: {0}'
                       ''.format(schema['bsr']))
    schema_info.append('chewBBACA version: {0}'
                       ''.format(schema['chewBBACA_version']))
    schema_info.append('Genetic code: {0}'
                       ''.format(schema['translation_table']))
    schema_info.append('Minimum sequence length: {0}'
                       ''.format(schema['minimum_locus_length']))
    schema_info.append('Sequence length variation threshold: {0}'
                       ''.format(schema['size_threshold']))
    schema_info.append('Clustering word size: {0}'
                       ''.format(schema['word_size']))
    schema_info.append('Clustering similarity: {0}'
                       ''.format(schema['cluster_sim']))
    schema_info.append('Representative similarity filter: {0}'
                       ''.format(schema['representative_filter']))
    schema_info.append('Intracluster similarity filter: {0}'
                       ''.format(schema['intraCluster_filter']))
    schema_info.append('Creation date: {0}'
                       ''.format(schema['dateEntered']))
    schema_info.append('Last modified: {0}'
                       ''.format(schema['last_modified']))

    return schema_info


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', type=str, required=True,
                        dest='stats_mode', choices=['species', 'schemas'],
                        help='The process can retrieve the list of species '
                             '("species" option) in the Chewie-NS or the '
                             'list of schemas for a species '
                             '("schemas" option).')

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

    parser.add_argument('--sp', type=str, required=False,
                        dest='species_id', default=None,
                        help='The integer identifier of a '
                             'species in the Chewie-NS.')

    parser.add_argument('--sc', type=str, required=False,
                        dest='schema_id', default=None,
                        help='The integer identifier of a schema in '
                             'the Chewie-NS.')

    args = parser.parse_args()

    return [args.stats_mode, args.nomenclature_server,
            args.species_id, args.schema_id]


def main(mode, base_url, species_id, schema_id):

    headers_get = cnst.HEADERS_GET_JSON

    print('\nRetrieving data...')
    if mode == 'species':
        stats = species_stats(base_url, headers_get)
    elif mode == 'schemas':
        if species_id is not None:
            if schema_id is not None:
                stats = single_schema(species_id, schema_id,
                                      base_url, headers_get)
            else:
                stats = single_species(species_id, base_url,
                                       headers_get)
        else:
            sys.exit('\nPlease provide a valid species identifier '
                     'to get the list of available schemas.\n')

    # print stats
    stats_text = '\n'.join(stats)
    print('\n{0}\n'.format(stats_text))


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2], args[3])
