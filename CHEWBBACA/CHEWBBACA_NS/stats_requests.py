#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 15:37:16 2019

@author: rfm

DESCRIPTION

    Script used by users to get useful stats:
        -

"""


import sys
import argparse
import requests
import concurrent.futures
from getpass import getpass

from SPARQLWrapper import SPARQLWrapper

from utils import auxiliary_functions as aux

# way of importing scripts or modules that are not in the current directory
#sys.path.insert(1, '/home/rfm/Desktop/rfm/Lab_Software/Chewie_NS/NS_tests/ns_security_backup_20191009')
#
#
#virtuoso_server = SPARQLWrapper(app.config['LOCAL_SPARQL'])
#virtuoso_server = SPARQLWrapper('http://localhost:8890/sparql')
#url_send_local_virtuoso = 'http://localhost:8890/DAV/test_folder/data'


def simple_request(base_url, headers_get, endpoint_list):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    url = aux.make_url(base_url, *endpoint_list)
    res = requests.get(url, headers=headers_get, timeout=30)

    return res


def species_stats(base_url, headers_get, endpoint_list):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    species_namids = species_list(base_url, headers_get, endpoint_list)

    species_stats = []
    species_stats.append('{:<30}  {:^10}  {:^10}  {:^10}'.format('Species',
                                                                         'id',
                                                                         '#schemas',
                                                                         '#loci'))
                                                                         #'#alleles'))
    for sp in species_namids:
        species = sp[0]
        species_id = sp[1]

        # get number of schemas
        species_schemas = count_schemas(base_url,
                                        headers_get,
                                        ['species', species_id, 'schemas'])

        # get number of loci
        species_total_loci, species_loci_ids = count_loci(base_url,
                                                          headers_get,
                                                          ['species', species_id, 'loci'])

        # get total number of alleles
#        species_alleles = 0
#        for locus in species_loci_ids:
#            species_alleles += count_loci_alleles(base_url,
#                                                  headers_get,
#                                                  ['loci', locus, 'alleles'])
        
        # get total number of isolates



        species_stats.append('{:<30}  {:^10}  {:^10}  {:^10}'.format(species,
                                                                     species_id,
                                                                     str(species_schemas),
                                                                     str(species_total_loci)))
                                                                             #str(species_alleles)))

    return species_stats


def species_list(base_url, headers_get, endpoint_list):
    """
    """

    res = simple_request(base_url, headers_get, endpoint_list)
    res = res.json()
    species_lst = []
    for sp in res:
        species = sp['name']['value']
        species_url = sp['species']['value']
        species_id = species_url.split('/')[-1]

        species_lst.append((species, species_id))

    return species_lst


def count_loci(base_url, headers_get, endpoint_list):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    res = simple_request(base_url, headers_get, endpoint_list)
    res = res.json()

    loci_ids = []
    if 'message' not in res:
        loci = res['Loci']
        for locus in loci:
            locus_id = locus['locus']['value'].split('/')[-1]
            loci_ids.append(locus_id)

        total_loci = len(res['Loci'])

    else:
        total_loci = 0

    return [total_loci, loci_ids]


def count_schemas(base_url, headers_get, endpoint_list):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    res = simple_request(base_url, headers_get, endpoint_list)
    res = res.json()

    if 'message' not in res:
        total_schemas = len(res)
    else:
        total_schemas = 0

    return total_schemas


def schemas_info(base_url, headers_get, endpoint_list):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    res = simple_request(base_url, headers_get, endpoint_list)
    res = res.json()

    sc_info = {}
    if 'message' not in res:
        for s in res:
            schema_id = s['schemas']['value'].split('/')[-1]
            schema_desc = s['name']['value']
            sc_info[schema_id] = schema_desc

    return sc_info


def count_loci_alleles(loci_info):
    """
    """
    base_url = loci_info[0]
    headers_get = loci_info[1]
    endpoint_list = loci_info[2]

    # unpack list of sequential endpoints and pass to create URI
    res = simple_request(base_url, headers_get, endpoint_list)
    res = res.json()

    total_alleles = len(res)

    return total_alleles


def species_schemas_stats(base_url, headers_get, endpoint_list, species_id, workers):
    """
    """

    given = species_id

    # unpack list of sequential endpoints and pass to create URI
    species_namids = species_list(base_url, headers_get, endpoint_list)

    species_names = {t[0]:t[1] for t in species_namids}
    species_ids = {t[1]:t[0] for t in species_namids}

    sp = species_ids.get(species_id, 0)
    if sp == 0:
        sp = species_id
        species_id = species_names.get(species_id, 0)
        if species_id == 0:
            return 'There is no taxon name or identifier "{0}" in server.'.format(given)

    schema_stats = []
    schema_stats.append('{0} available schemas\n'.format(sp))
    schema_stats.append('{:<10}  {:^10}  {:^20}'.format('Schema_id',
                                                                '#loci',
                                                                '#alleles',
                                                                'description'))

    # get schemas for selected species
    schemas = schemas_info(base_url, headers_get, ['species', species_id, 'schemas'])

    for sid, sdesc in schemas.items():

        # get number of loci in schema
        schema_total_loci, schema_loci_ids = count_loci(base_url,
                                                        headers_get,
                                                        ['species', species_id,
                                                         'schemas', sid, 'loci'])

        # get number of alleles
        schema_alleles = 0
        schema_loci_ids = [[base_url, headers_get, ['loci', locus, 'alleles']] for locus in schema_loci_ids]
        #print(schema_loci_ids)
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            for res in executor.map(count_loci_alleles, schema_loci_ids):
                schema_alleles += res
            print(schema_alleles)

        schema_stats.append('{:^10}  {:^10}  {:^10}  {:^20}'.format(sid,
                                                                    schema_total_loci,
                                                                    schema_alleles,
                                                                    sdesc))

    return schema_stats


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', type=str, required=True, dest='stats_mode',
                        help='')

    parser.add_argument('--url', type=str, required=False, dest='base_url',
                        default='http://127.0.0.1:5000/NS/api/',
                        help='')

    parser.add_argument('--taxon', type=str, required=False, dest='taxon',
                        default='0', help='')

    args = parser.parse_args()

    return [args.stats_mode, args.base_url, args.taxon]

# add multithreading to GET number of alleles
def main(mode, base_url, taxid):
    """
    """

    # request username and password
    print('\nPlease provide login credentials:')
    user = input('USERNAME: ')
    password = getpass('PASSWORD: ')
    print()

    # get token
    token = aux.login_user_to_NS(base_url, user, password)

    headers_get = {'X-API-KEY': token,
                   'accept': 'application/json'}

    if mode == 'species_list':
        endpoint_list = ['species', 'list']
        stats = species_stats(base_url, headers_get, endpoint_list)

        # print stats
        stats_text = '\n'.join(stats)
        print(stats_text)

    elif mode == 'schema_info':

        species_id = taxid

        if species_id != '0':
            endpoint_list = ['species', 'list']
            stats = species_schemas_stats(base_url, headers_get,
                                          endpoint_list, species_id, 30)

            # print stats
            stats_text = '\n'.join(stats)
            print(stats_text)
        else:
            print('\nPlease provide a valid taxon name or identifier '
                  'to get the list of available schemas.\n')


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2])
