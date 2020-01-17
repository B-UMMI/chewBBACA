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
import shutil
import pickle
import random
import argparse
import requests
from copy import deepcopy
import concurrent.futures
from getpass import getpass
from itertools import repeat
import extra_scripts.utils as ut
from collections import defaultdict
from SPARQLWrapper import SPARQLWrapper, JSON
from extra_scripts import PrepExternalSchema

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

virtuoso_server = SPARQLWrapper('http://sparql.uniprot.org/sparql')


def determine_upload(local_schema_loci, ns_schema_loci,
                     ns_schema_locid_map, schema_path,
                     temp_path, base_url, headers_get):
    """
    """

    local_uniq = {}
    comp = []
    for locus in local_schema_loci:
        #updated_locus = []
        local_locus = locus
        locus_id = local_locus.rstrip('.fasta')
        pickled_file = os.path.join(temp_path, '{0}_pickled'.format(locus_id))
        with open(pickled_file, 'rb') as pf:
            locus_sequences = pickle.load(pf)
        local_sequences = local_schema_loci[locus]

        if local_locus in ns_schema_loci:
            local_uniq[locus] = []

            # get uniprot info for locus
            ns_uri = ns_schema_locid_map[local_locus][0]
            ns_locus_id = ns_uri.split('/')[-1]

            ns_uniprot = simple_get_request(base_url, headers_get,
                                            ['loci', ns_locus_id, 'uniprot'])
            locus_name = ns_uniprot.json()['UniprotInfo'][0]['UniprotSName']['value']
            locus_label = ns_uniprot.json()['UniprotInfo'][0]['UniprotLabel']['value']
            locus_annotation = ns_uniprot.json()['UniprotInfo'][0]['UniprotURI']['value']

            # after sync, the allele identifiers of the sequences
            # exclusive to the local schema should be ok to use
            # as identifiers for new alleles in the NS
            for rec in local_sequences:
                seqid = rec[0]
                seqid = seqid.replace('*', '')
                seq = rec[1]

                local_uniq[local_locus].append([ns_uri,
                                               locus_name,
                                               locus_label,
                                               locus_annotation,
                                               seqid,
                                               seq])

                # delete entry with '*'
                allele_id = int(seqid.split('_')[-1])
                locus_sequences[allele_id] = (str(allele_id), locus_sequences[allele_id][1])
            
            # save updated pickle
            os.remove(pickled_file)
            with open(pickled_file, 'wb') as pf:
                pickle.dump(locus_sequences, pf)
            
            comp.append(locus)

    print('Found a total of {0} new sequences in'
          ' local schema.'.format(len(local_uniq)))

    return [local_uniq, comp]


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
    url = ut.make_url(base_url, *endpoint_list)

    res = requests.get(url, headers=headers, timeout=30)

    return res


def select_name(result):
    """
    """
    
    name = ''
    url = ''
    label = ''

    aux = result["results"]["bindings"]
    
    for elem in aux:
        if 'fname' in elem.keys():
            name = str(elem['fname']['value'])
        elif 'sname2' in elem.keys():
            name = str(elem['sname2']['value'])
        elif 'label' in elem.keys():
            name = str(elem['label']['value'])
            
        if 'label' in elem.keys():
            label = str(elem['label']['value'])

        url = str(elem['seq']['value'])
        
        break
        
    return [name, url, label]
    

def get_data(sparql_query):
    """ Retrieve annotations from uniprot
    """
    
    virtuoso_server.setReturnFormat(JSON)
    virtuoso_server.setTimeout(10)

    url = ''
    name = ''
    prev_name = ''
    label = ''
    found = False
    unpreferred_names = ['Uncharacterized protein', 'hypothetical protein', 'DUF']

    # implement more than 1 retry!!!
    alleles = len(sparql_query)
    a = 0
    while found is False:
        virtuoso_server.setQuery(sparql_query)

        try:
            result = virtuoso_server.query().convert()
            # Slowing the requests down so that Uniprot doesn't blacklist us :)
            time.sleep(random.randint(1, 3))

            name, url, label = select_name(result)

            if prev_name == '' and not any([n in name for n in unpreferred_names]):
                prev_name = name
                found = True

        except:
            #print("A request to uniprot timed out, trying new request")
            time.sleep(5)
            result = virtuoso_server.query().convert()
            name, url, label = select_name(result)
            if prev_name == '' and not any([n in name for n in unpreferred_names]):
                prev_name = name
                found = True

        a += 1
        if a == alleles:
            found = True
        

    return (prev_name, label, url)


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
        


def uniprot_query(sequence):
    """
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


def check_seq_req(headers_get, url):
    """
    """
    
    url_req = url[0]
    
    res = requests.get(url_req, headers = headers_get, timeout = 30)
    
    return (res, url[1])


def check_seq2(fasta, URL, headers_get, cpu):
    """ Checks if a sequence already exists in NS
    """
    
    cpu = 30
    
    # Get the sequences of each record from the fasta file
    sequences = [(str(rec.seq), rec.id) for rec in SeqIO.parse(fasta, "fasta")]
    
    responses = {}
    
    responses[fasta] = []
    
    urls = []
    
    for seq in sequences:
        
        url_seq_info = ut.make_url(URL, "sequences", "seq_info", sequence=seq[0])
        urls.append((url_seq_info, seq[1]))
    
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=cpu) as executor:
        for result in executor.map(check_seq_req, repeat(headers_get), urls):
            responses[fasta].append(result)
    
    return responses


def build_fasta_files(new_allele_seq_dict, path2schema, path_new_alleles):
    """
    """
    
    for name, allele_dict in new_allele_seq_dict.items():
        
        for allele_id, seq in allele_dict.items():
            
            auxDict = {}
            auxname = name.split(".")
            listIndexAux = []

		# create a dictionary with the info from fasta to build the fasta file; {id:sequence}
        listIndexAux.append(int(allele_id))
        auxDict[int(allele_id)] = f"{seq}"

		# sort by allele id, create the fasta string and write the file
        listIndexAux.sort()
        auxString = ''
        for index in listIndexAux:
            auxString += f">{auxname[0]}_{index}\n{auxDict[index]}\n"
        
        # if the fasta already exists, add new allele
        if os.path.exists(os.path.join(path2schema, name)):
            
            with open(os.path.join(path2schema, name), 'a') as f:
                f.write(auxString)
            
            # Copy the updated file to new_alleles dir 
            new_alleles_path = shutil.copy(os.path.join(path2schema, name), path_new_alleles)
        
        # write new file to schema dir and new_alleles dir
        else:
            
            # Schema dir
            with open(os.path.join(path2schema, name), 'w') as f1:
                f1.write(auxString)
            
            # New alleles dir
            with open(os.path.join(path_new_alleles, name), 'w') as f2:
                f2.write(auxString)
                
    
    print(f"{len(new_allele_seq_dict)} files have been written")

    return True


def retrieve_latest_alleles(loci_new_alleles, server_time, schema_uri, new_count, headers_get):
    """
    """

    # request the new alleles starting on the date given
    uri = ut.make_url(schema_uri, 'loci', date=server_time)

    # get the new alleles
    response = requests.get(uri, headers=headers_get)

    response_content = response.json()
    # get headers info with 
    response_headers = response.headers
    # get info about sequences that were added since last date
    new_ns_alleles = response_content['newAlleles']

    if len(new_ns_alleles) > 0:
        previous_server_time = server_time
        # get date of last added allele
        server_time = response_headers['Last-Allele']
    else:
        # get current server date if no sequence were added
        # since last server time
        server_time = response_headers['Server-Date']
        return (loci_new_alleles, True, server_time, new_count)

    # group retrieved alleles by locus
    already = 0
    for allele in new_ns_alleles:

        locus = '{0}{1}'.format(allele['locus_name']['value'], '.fasta')
        allele_id = allele['allele_id']['value']
        sequence_uri = allele['sequence']['value']

        if locus in loci_new_alleles:
            if allele_id in loci_new_alleles[locus]:
                if loci_new_alleles[locus][allele_id] == sequence_uri:
                    print('This was already here!')
                    already += 1

        loci_new_alleles[locus][allele_id] = sequence_uri

    # keep count of the number of retrieved alleles
    new_count += len(response_content['newAlleles']) - already

    # if we got the maximum number of alleles that can be returned
    # in one pass, return False to re-run the function until we get all
    # new alleles
    if 'All-Alleles-Returned' in response_headers:
        returned = response_headers['All-Alleles-Returned']
        all_retrieved = False if returned == 'False' else True
        if all_retrieved == True:
            server_time = response_headers['Server-Date']
        return (loci_new_alleles, all_retrieved, server_time, new_count)
    # return True to stop since there are no more new alleles
    elif 'Server-Date' in response_headers:
        return (loci_new_alleles, True, server_time, new_count)


def get_allele_seq(headers_get, uri):
    """
    """

    seq_hash = uri.split('/sequences/')[1]

    url = ut.make_url(uri.split(seq_hash)[0], 'seq_info', seq_id=seq_hash)

    r = requests.get(url, headers=headers_get, timeout=30)

    r_content = r.json()

    allele_seq = r_content['sequence']['value']

    return [uri, allele_seq]


def get_new_alleles_seqs(alleles_info, headers_get):
    """
    """
    
    new_alleles = {}
    for locus in alleles_info.keys():
        allele_dict = alleles_info[locus]
        # invert locus dictionary to get mapping from allele
        # sequence hash to allele identifier
        allele_hash_map = {v:k for k, v in allele_dict.items()}
        # get all endpoints for the new sequences
        sequence_uri_list = list(allele_dict.values())
        
        new_alleles[locus] = {}

        # get allele DNA sequences
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            for result in executor.map(get_allele_seq, repeat(headers_get), sequence_uri_list):
                new_alleles[locus][allele_hash_map[result[0]]] = result[1]

    return new_alleles


def update_loci_files(new_alleles, local_loci, schema_dir, temp_dir):
    """
    """

    # Check if new alleles are already on schema
    not_in_ns = {}
    pickled_loci = {}
    #print(new_alleles.keys())
    for locus in new_alleles:
        #print(locus)
        not_in_ns[locus] = []
        pickled_loci[locus] = ''
        if locus in local_loci:
            # get locus file identifier without '.fasta'
            locus_id = locus.rstrip('.fasta')
            # get latest locus alleles retrieved from the NS
            locus_ns_seqs = {seq: seqid
                             for seqid, seq in new_alleles[locus].items()}

            # full path for local locus file
            locus_file = os.path.join(schema_dir, locus)
            # full path for temporary file that will
            # substitute the current file
            temp_file = os.path.join(temp_dir, '{0}_pickled'.format(locus_id))

            # get local locus sequences
            with open(locus_file, 'r') as sf:
                # read sequences in local file
                records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                           for rec in SeqIO.parse(sf, 'fasta')}

            # deepcopy records to alter and sync with the NS locus
            altered_locus = deepcopy(records)
            id_map = {}
            for seqid, seq in altered_locus.items():
                if '*' in seq[0]:
                    id_map[int(seq[0][1:])] = seqid
                else:
                    id_map[int(seq[0])] = seqid

            # determine max integer identifier to later increment
            max_id = max(list(id_map.keys()))

            # initalize list
            switched = []
            # for each local record
            for seqid, seq in records.items():
                # get full sequence idetifier
                current_seqid = seqid
                # get integer identifier
                allele_id = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
                # get DNA sequence
                current_seq = seq[1]
                # if the DNA sequence is in the NS but was added to the
                # local schema through allele call
                if current_seq in locus_ns_seqs and '*' in current_seqid:
                    # get allele identifier in the NS
                    ns_id = int(locus_ns_seqs[current_seq])
                    # allele has the same integer identifier
                    # in local and NS schemas
                    if ns_id == allele_id:
                        # simple '*' replace suffices
                        new_id = current_seqid.replace('*', '')
                        # add record without '*'
                        altered_locus[new_id] = [str(allele_id), seq[1]]
                        # delete record with '*'
                        del altered_locus[current_seqid]
                        # delete from latest NS sequences
                        del locus_ns_seqs[current_seq]
                    # allele may have a different identifier between the local
                    # and NS schemas
                    elif ns_id != allele_id:
                        # the same allele identifier might have been attributed
                        # to different allele sequences
                        if ns_id in id_map:
                            # NS identifier will be the new local identifier
                            new_id = '{0}_{1}'.format(locus.rstrip('.fasta'),
                                                      ns_id)
                            # the local allele identifier will replace the
                            # identifier of the local sequence that currently has
                            # the same identifier as the NS allele
                            old_id = current_seqid
                            # where is the '*'
                            old_seq = altered_locus['{0}_*{1}'.format(locus_id, ns_id)][1]
                            altered_locus[new_id] = [str(ns_id), current_seq]
                            switched.append(current_seq)
                            if old_seq not in switched:
                                altered_locus[old_id] = [allele_id, old_seq]

                            del altered_locus['{0}_*{1}'.format(locus_id, ns_id)]
                            del locus_ns_seqs[current_seq]
                        # the identifier attributed to the allele in the NS
                        # might be greater than the number of alleles in the
                        # local schema
                        elif ns_id not in id_map:
                            new_id = '{0}_{1}'.format(locus_id, ns_id)
                            altered_locus[new_id] = [str(ns_id), current_seq]
                            # simply delete the old record with the same allele
                            del altered_locus[current_seqid]
                            del locus_ns_seqs[current_seq]

                # local allele may not be in the NS
                elif current_seq not in locus_ns_seqs and '*' in current_seqid:
                    # local sequence may have an identifier
                    # that was attributed to another allele in the NS
                    if str(allele_id) in list(locus_ns_seqs.values()):
                        # attribute identifier greater than total number
                        # of alleles
                        max_id += 1
                        new_id = '{0}_*{1}'.format(locus_id, max_id)
                        altered_locus[new_id] = ['*{0}'.format(max_id), current_seq]

                        # get sequence in the NS that has that allele
                        ns_seq = [k for k, v in locus_ns_seqs.items() if int(v) == allele_id]
                        updated_seqid = current_seqid.replace('*', '')
                        # attribute the NS sequence to the allele identifier
                        altered_locus[updated_seqid] = [str(allele_id), ns_seq[0]]

                        # delete identifier with '*'
                        del altered_locus[current_seqid]
                        del locus_ns_seqs[ns_seq[0]]

                        # store info about sequences that are not in the NS
                        not_in_ns[locus].append(('{0}_*{1}'.format(locus_id,
                                                                  max_id),
                                                                  current_seq))
                    # local sequence has an identifier that is not attributed
                    # to any sequence in the NS
                    else:
                        not_in_ns[locus].append((current_seqid, current_seq))
                # local sequence was previously synced
                elif current_seq in locus_ns_seqs and '*' not in current_seqid:
                    del locus_ns_seqs[current_seq]

        # get local records that were synced in terms of identifiers
        updated_records = {}
        for seqid, seq in altered_locus.items():
            int_seqid = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
            updated_records[int_seqid] = (seq[0], seq[1])

        # check if there are any completely new sequences that are only present in the NS
        if len(locus_ns_seqs) > 0:
            ns_rest = {int(seqid): (seqid, seq) for seq, seqid in locus_ns_seqs.items()}
            updated_records = {**updated_records, **ns_rest}

        # sort by integer identifier
        #updated_records = sorted(updated_records, key=lambda x: x[0])

        with open(temp_file, 'wb') as pl:
            pickle.dump(updated_records, pl)

        pickled_loci[locus] = temp_file

    return [not_in_ns, pickled_loci]


def update_local_schema(last_sync_date, schema_uri, schema_dir, temp_dir,
                        bsr, ns_url, headers_get, headers_post):
    """
    """

    # get list of sequences that are new considering the last date
    new_count = 0
    new_seqs = False
    server_time = last_sync_date
    loci_new_alleles = defaultdict(dict)

	# get all sequences until the number of new sequences is less than 100k, the maximum the server return is 100k
    # start getting new sequences that were added to the server and that are not in local schema
    print('Determining if there are new sequences in the NS...')
    print('Getting new alleles added to the NS since {0}'.format(str(server_time)))
    while not new_seqs:
        loci_new_alleles, new_seqs, server_time, new_count = retrieve_latest_alleles(loci_new_alleles,
                                                                                     server_time,
                                                                                     schema_uri,
                                                                                     new_count,
                                                                                     headers_get)
    print('Retrieved {0} new alleles added since {1}'.format(new_count, last_sync_date))

    # if new_alleles dict has no entries
    # there are no new alleles
    if not bool(loci_new_alleles):
        return 'There were no new alleles to retrieve...'

    # get DNA sequences for retrieved alleles
    print('Getting DNA sequences of latest alleles...')
    fasta_response = get_new_alleles_seqs(loci_new_alleles, headers_get)

    # get list of local schema files
    schema_fasta_files = [file for file in os.listdir(schema_dir) if '.fasta' in file]

    not_in_ns, pickled_loci = update_loci_files(fasta_response,
                                                schema_fasta_files,
                                                schema_dir,
                                                temp_dir)

    return [not_in_ns, pickled_loci, server_time]


def load_binary(parent_dir, file_name):
    """
    """

    binary_file = os.path.join(parent_dir,
                               file_name)

    if os.path.exists(binary_file):
        with open(binary_file, 'rb') as bf:
            content = pickle.load(bf)
        return content
    else:
        return False


def create_allele_data(allele_seq_list, species_name, check_cds,
                       headers_post, user_id):
    """
    """

    post_inputs = []
    for allele in allele_seq_list:
        allele_uri = '{0}/alleles/{1}'.format(allele[0], (allele[4]).split('_')[-1])
        post_inputs.append((allele[-1], allele[1], allele[2], allele[3],
                            allele[0], species_name,
                            check_cds, headers_post, allele_uri, user_id))

    return post_inputs


def post_allele(input_stuff):
    """ Adds a new allele to the NS.

        Args:
            A tuple with 8 elements:
                - sequence (str): the DNA sequence to send to NS.
                - name (str): protein annotation name.
                - label (str): protein annotation label.
                - uniprot_url (str): URL to the UniProt entry.
                - loci_url (str): URI of the locus in NS.
                - species_name (str): name of the species the allele
                belongs to.
                - cds_check (bool): if the sequence must be a complete CDS.
                - headers_post (dict): headers for the POST method used to
                insert data into the NS.
        Returns:
            response (requests.models.Response): response object from
            the POST method.
    """

    # getting inputs from multithreading
    sequence = input_stuff[0]
    name = input_stuff[1]
    label = input_stuff[2]
    uniprot_url = input_stuff[3]
    loci_url = input_stuff[4]
    species_name = input_stuff[5]
    cds_check = input_stuff[6]
    headers_post = input_stuff[7]
    allele_uri = input_stuff[8]
    user_id = input_stuff[9]

    # Build the url for loci/loci_id/alleles
    url = ut.make_url(loci_url, 'alleles')

    params = {}
    params['sequence'] = sequence
    params['species_name'] = species_name
    params['enforceCDS'] = cds_check
    params['uniprot_url'] = uniprot_url
    params['uniprot_label'] = label
    params['uniprot_sname'] = name
    params['input'] = 'auto'
    params['sequence_uri'] = allele_uri
    params['user_id'] = user_id

    response = requests.post(url, data=json.dumps(params),
                             headers=headers_post, timeout=30)

    return response


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-schema', type=str, dest='schema_dir', required=True,
                        help='Path to the directory with the schema to be'
                             'synced.')

    parser.add_argument('--cores', type=int, required=False, dest='core_num', 
                        default=1, help='Number of CPU cores to use to construct'
                        ' the complete schema based on the synced/updated loci'
                        ' files.')

    parser.add_argument('--ns_url', type=str, required=False, dest='ns_url', 
                        default='http://127.0.0.1:5000/NS/api/',
                        help='The base URI for the NS endpoints.')

    parser.add_argument('--sub', type=str, required=False, dest='submit', 
                        default='no',
                        help='If the process should detect new alleles'
                        'in the local schema and send them to the NS. (only'
                        ' authorized users can submit new alleles).')

    args = parser.parse_args()

    return [args.schema_dir, args.core_num, args.ns_url,
            args.submit]


#bsr = 0.6
#core_num = 6
#submit = True
#schema_dir = '/home/rfm/Desktop/rfm/Lab_Software/Chewie_NS/NS_tests/test_sync_schema/test_schema/ypestis_testsync2_schema'
#temp_dir = '/home/rfm/Desktop/rfm/Lab_Software/Chewie_NS/NS_tests/test_sync_schema/test_schema/temp'
#ns_url = 'http://127.0.0.1:5000/NS/api/'


def main(schema_dir, core_num, ns_url, submit):

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
        token = ut.login_user_to_NS(ns_url, user, password)
        # if login was not successful, stop the program
        if token is False:
            message = '403: Invalid credentials.'
            print(message)
            return message

    # Define the headers of the requests
    headers_get = {'Authorization': token,
                   'accept': 'application/json'}

    headers_post = {'Authorization': token,
                    'Content-type': 'application/json',
                    'accept': 'application/json'}

    # get schema configs
    ns_vars = load_binary(schema_dir, '.ns_config')
    if ns_vars is not False:
        last_sync_date, schema_uri = ns_vars
    else:
        print('Could not find schema NS variables.\n'
              'Cannot sync schema.')
        return 1

    schema_params = load_binary(schema_dir, '.schema_config')
    if schema_params is False:
        print('Could not find schema configurations file.\n'
              'Cannot sync schema')
        return 1

    # check if schema exists in the NS
    schema_response = requests.get(schema_uri,
                                   headers=headers_get,
                                   timeout=5)

    if schema_response.status_code > 201:
        print('There is no schema with URI "{0}" in the NS.')
        return 1

    # check if list of loci is the same in local and NS schemas

    # update local schema
    # download new sequences and change local identifiers
    bsr = float(schema_params['bsr'])

    #### change last_sync_date just to be able to test
    last_sync_date = '2019-09-07T13:31:04.677308'

    # Create an intermediate dir for the new alleles
    temp_dir = os.path.join(os.path.dirname(schema_dir), 'temp')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    not_in_ns, pickled_loci, server_time = update_local_schema(last_sync_date,
                                                               schema_uri,
                                                               schema_dir,
                                                               temp_dir,
                                                               bsr,
                                                               ns_url,
                                                               headers_get,
                                                               headers_post)

    if submit == 'yes':

        # verify user role to check permission
        user_info = simple_get_request(ns_url, headers_get,
                                       ['user', 'current_user'])
        user_info = user_info.json()
        user_role = any(role in user_info['roles']
                        for role in ['Admin', 'Contributor'])

        if not user_role:
            print('\n403: Current user has no Administrator '
                  'or Contributor permissions.\n'
                  'Not allowed to upload schemas.')
            #return 403

        user_id = str(user_info['id'])
        headers_post['user_id'] = user_id

        species_id = schema_uri.split('/')[-3]

        # Get the name of the species from the provided id
        # or vice-versa
        species_info = species_ids(species_id, ns_url, headers_get)
        if isinstance(species_info, list):
            species_id, species_name = species_info
            print('\nNS species with identifier {0} is {1}.'.format(species_id,
                                                                    species_name))
        else:
            print('\nThere is no species with the provided identifier in the NS.')
            #return 1

        # compare list of genes, if they do not intersect, halt process
        # get list of loci for schema in the NS
        schema_id = schema_uri.split('/')[-1]
        ns_loci_get = simple_get_request(ns_url, headers_get,
                                         ['species', species_id,
                                          'schemas', schema_id,
                                          'loci'])
        # get loci files names from response
        ns_schema_loci = []
        ns_schema_locid_map = {}
        for l in ns_loci_get.json()['Loci']:
            locus_file = l['original_name']['value']
            locus_name = l['name']['value'] + '.fasta'
            ns_schema_loci.append(locus_name)
            locus_uri = l['locus']['value']
            ns_schema_locid_map[locus_name] = (locus_uri, locus_file)

        completed = []
        incomplete = []
        for locus in not_in_ns:
            if locus in ns_schema_loci:
                if len(not_in_ns[locus]) > 0:
                    incomplete.append(locus)
                else:
                    completed.append(locus)
        
        if len(incomplete) > 0:

            uniq_local = {k:v for k, v in not_in_ns.items() if k in incomplete}
            loci_uris = {k:v for k, v in ns_schema_locid_map.items() if k in incomplete}
    
            # determine sequences that are not in the NS
            # we synced before so we just need to go to each file and
            # identify the alleles with '*' as the alleles to be added!
            # we have to submit them and alter locally too...
            upload, comp = determine_upload(uniq_local, incomplete,
                                            loci_uris, schema_dir,
                                            temp_dir, ns_url, headers_get)

            completed.extend(comp)

            # add sequences to the NS
            # use multiprocessing
            alleles_data = []
            for locus, info in upload.items():
                allele_seq_list = info
                post_data = create_allele_data(allele_seq_list,
                                               species_name,
                                               True,
                                               headers_post,
                                               user_id)
                alleles_data.append(post_data)

            for inlist in alleles_data:
                post_results = []
                total_inserted = 0
                total_alleles = len(inlist)

                workers = 10

                with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
                    # Start the load operations and mark each future with its URL
                    for res in executor.map(post_allele, inlist):
                        post_results.append(res)
                        total_inserted += 1
                        print('\r',
                              'Processed {0}/{1} alleles.'.format(total_inserted,
                                                                  total_alleles),
                              end='')
        else:
            print('There are no new alleles in the local schema that '
                  'need to be sent to the NS.')

    # change profiles in the master file so that local
    # profiles have updated allele identifiers

    # change pickled files to FASTA files
    for locus in completed:
        locus_id = locus.rstrip('.fasta')
        pickled_file = pickled_loci[locus]
        with open(pickled_file, 'rb') as pf:
            locus_sequences = pickle.load(pf)

        natsorted_locus = sorted(locus_sequences)

        fasta_path = os.path.join(temp_dir, locus)
        records = ['>{0}_{1}\n{2}'.format(locus_id,
                                          locus_sequences[seqid][0],
                                          locus_sequences[seqid][1])
                   for seqid in natsorted_locus]

        fasta_text = '\n'.join(records)

        with open(fasta_path, 'w') as fp:
            fp.write(fasta_text)

        os.remove(pickled_file)

    # Re-determine the representative sequences
    PrepExternalSchema.main(temp_dir, schema_dir,
                            core_num, bsr)

    # delete invalid alleles and genes files
    parent_dir = os.path.dirname(schema_dir)
    files = [os.path.join(parent_dir, file)
             for file in os.listdir(parent_dir)
             if 'invalid' in file]

    for f in files:
        os.remove(f)

    # delete temp directory
    shutil.rmtree(temp_dir)

    # update NS config file with latest server time
    ns_configs = os.path.join(schema_dir, '.ns_config')
    with open(ns_configs, 'wb') as nc:
        pickle.dump([server_time, schema_uri], nc)

    print('Done!')


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3])
