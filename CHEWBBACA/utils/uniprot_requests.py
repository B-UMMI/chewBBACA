#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to perform requests to
UniProts's SPARQL endpoint and process retrieved data.

Code documentation
------------------
"""


import time
from SPARQLWrapper import SPARQLWrapper, JSON

try:
    from utils import (constants as ct,
                       file_operations as fo)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo)


UNIPROT_SERVER = SPARQLWrapper(ct.UNIPROT_SPARQL)


def select_name(result):
    """ Extracts the annotation description from the result
        of a query to UniProt's SPARQL endpoint.

        Parameters
        ----------
        result : dict
            A dictionary with the results from querying
            the UniProt SPARQL endpoint.

        Returns
        -------
        A list with the following elements:
            name : str
                The annotation descrition.
            url : str
                The URL to the UniProt page about the protein record.
            label : str
                A label that has descriptive value.
    """

    url = ''
    name = ''
    label = ''

    i = 1
    found = False
    # get the entries with results
    aux = result['results']['bindings']
    total_res = len(aux)
    # only check results that are not empty
    if total_res > 0:
        # iterate over all results to find suitable
        while found is False:
            current_res = aux[i]
            res_keys = aux[i].keys()

            # annotation name can be associated
            # to different keys
            if 'fname' in res_keys:
                name = str(current_res['fname']['value'])
                found = True
            elif 'sname2' in res_keys:
                name = str(current_res['sname2']['value'])
                found = True
            elif 'label' in res_keys:
                name = str(current_res['label']['value'])
                found = True

            if 'label' in res_keys:
                label = str(current_res['label']['value'])
            else:
                label = name

            # get UniProt URL
            if 'uri' in res_keys:
                url = str(current_res['seq']['value'])
            elif 'seq' in res_keys:
                url = str(current_res['seq']['value'])

            if i == total_res:
                found = True

    return [name, url, label]


def uniprot_query(sequence):
    """ Constructs a SPARQL query to search for exact matches in
        UniProt's SPARQL endpoint.

        Parameters
        ----------
        sequence : str
            The Protein sequence that will be added to
            the query/searched for.

        Returns
        -------
        query : str
            The SPARQL query that will allow to search for
            exact matches in the UniProt database.
    """

    query = ('PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>  '
             'PREFIX up: <http://purl.uniprot.org/core/> '
             'select ?seq ?fname ?sname2 ?label  where {'
             '{?b a up:Simple_Sequence; rdf:value '
             '"'+sequence+'". ?seq up:sequence ?b. '
             'OPTIONAL{?seq up:submittedName ?sname. ?sname up:fullName ?sname2} '
             'OPTIONAL{?seq up:recommendedName ?rname.?rname up:fullName ?fname} }'
             'UNION{?seq a up:Sequence; rdf:value "'+sequence+'"; '
             'rdfs:label ?label. }}')

    return query


def get_data(sparql_query):
    """ Sends requests to query UniProts's SPARQL endpoint.

        Parameters
        ----------
        sparql_query : str
            SPARQL query.

        Returns
        -------
        result : dict
            Dictionary with data retrieved from UniProt.
    """

    tries = 0
    max_tries = 5
    success = False
    while success is False and tries < max_tries:
        try:
            UNIPROT_SERVER.setQuery(sparql_query)
            UNIPROT_SERVER.setReturnFormat(JSON)
            UNIPROT_SERVER.setTimeout(60)
            result = UNIPROT_SERVER.query().convert()
            success = True
        except Exception as e:
            tries += 1
            result = e
            time.sleep(1)

    return result


def get_proteomes(proteome_ids, output_dir):
    """ Downloads reference proteomes from UniProt's FTP.
    
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
