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
                       file_operations as fo,
                       fasta_operations as fao,
                       sequence_manipulation as sm)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_manipulation as sm)


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
    selected_name = ''
    selected_url = ''
    selected_label = ''

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
            elif 'sname2' in res_keys:
                name = str(current_res['sname2']['value'])
            elif 'label' in res_keys:
                name = str(current_res['label']['value'])

            if 'label' in res_keys:
                label = str(current_res['label']['value'])
            else:
                label = name

            # get UniProt URL
            if 'uri' in res_keys:
                url = str(current_res['seq']['value'])
            elif 'seq' in res_keys:
                url = str(current_res['seq']['value'])

            if any([term in name for term in ct.UNIPROT_UNINFORMATIVE]) is False:
                selected_name = name
                selected_url = url
                selected_label = label
                found=True
            else:
                if selected_name == '':
                    selected_name = name
                    selected_url = url
                    selected_label = label

            i += 1
            if i == total_res:
                found = True

    return [selected_name, selected_url, selected_label]


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


def get_annotation(gene, translation_table):
    """ Retrieves and selects annotation terms for a set of
        alleles from the same locus.

        Parameters
        ----------
        gene : str
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

    selected_name = ''
    selected_url = ''
    sequences = fao.import_sequences(gene)
    for seqid, sequence in sequences.items():
        protein_sequence = str(sm.translate_sequence(sequence,
                                                     table_id=translation_table))

        protein_sequence = protein_sequence.replace("*", "")

        query = uniprot_query(protein_sequence)
        result = get_data(query)

        name, url, label = select_name(result)

        lowercase_name = name.lower()
        if any([term in lowercase_name for term in ct.UNIPROT_UNINFORMATIVE]) is True:
            if selected_name == '':
                selected_name = name
                selected_url = url
            continue
        elif name == '':
            continue
        else:
            selected_name = name
            selected_url = url
            break

    return [gene, selected_name, selected_url]


def extract_proteome_terms(header_items):
    """ Extracts the sequence identifier, product name,
        gene name and species name fields from the sequence
        header of a reference proteome from Uniprot.

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
