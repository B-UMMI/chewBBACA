#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:20:29 2019

@author: pcerqueira
"""

import os
import json
import time
import random
import argparse
import requests
import concurrent.futures
from itertools import repeat
import extra_scripts.utils as ut
from SPARQLWrapper import SPARQLWrapper, JSON

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

NS_BASE_URL = 'http://127.0.0.1:5000/NS/api/'

virtuoso_server = SPARQLWrapper('http://sparql.uniprot.org/sparql')


def simple_get_request(base_url, headers, endpoint_list):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    url = ut.make_url(base_url, *endpoint_list)
    res = requests.get(url, headers=headers, timeout=30)

    return res


def simple_post_request(base_url, headers, endpoint_list, data):
    """
    """

    # unpack list of sequential endpoints and pass to create URI
    url = ut.make_url(base_url, *endpoint_list)
    res = requests.post(url, data=json.dumps(data), headers=headers)

    return res


def check_seq(headers_get, url):
    """ Checks if a DNA sequence is in the NS.

        Args:
            headers_get (dict): headers for the GET method used to
            get data from the API endpoints.
            url (str): the endpoint url that will be queried to know
            if the DNA sequence is in the NS.
        Returns:
            response_tup (tup): a tuple with the query response and
            the sequence identifier.
    """

    url_req = url[0]
    seqid = url[1]
    response = requests.get(url_req, headers=headers_get, timeout=30)
    response_tup = (response, seqid)

    return response_tup
    

def check_gene_seqs(fasta, species, base_url, headers_get, threads):
    """ Checks if the sequences in a FASTA file are in the NS.

        Args:
            fasta (str): path to the FASTA file with the sequences.
            species (str): species identifier in the NS.
            base_url (str): the base URL for the NS, used to concatenate
            with a list of elements and obtain endpoints URL.
            headers_get (dict): headers for the GET method used to
            get data from the API endpoints.
            cpu (int): number of workers to pass to multi-threading.
        Returns:
            responses (dict): a dictionary with the FASTA file path as key
            and tuples with the query response and sequence identifiers
            as values.
    """

    # Get the sequences of each record from the fasta file
    sequences = [(str(rec.seq), rec.id) for rec in SeqIO.parse(fasta, "fasta")]

    responses = {}
    responses[fasta] = []

    # create url for every sequence query
    urls = []
    for seq in sequences:
        url_seq_info = ut.make_url(base_url, 'sequences', 'seq_info',
                                   sequence=seq[0], species_id=species)
        urls.append((url_seq_info, seq[1]))

    # check if sequences are in the server with multithreading
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        for result in executor.map(check_seq, repeat(headers_get), urls):
            responses[fasta].append(result)

    return responses


def post_loci(base_url, species, headers_post, locus_prefix, keepFileName, gene):
    """
    """
    

    #### THE LOCI MUST EXIST IN ORDER TO ADD SCHEMA LOCI
    
    # Build the url for loci/list
    url_sp_loci = ut.make_url(base_url, "loci", "list")
    
    # Define POST request parameters
    params = {}
    params['prefix'] = locus_prefix
    
    if keepFileName:
        params['locus_ori_name'] = gene
        
    
#    sucess_send = False
#    waitFactor = 4
#    while not sucess_send:
    
    # Add locus to species
    r_sp_loci = requests.post(url_sp_loci, data = json.dumps(params), headers = headers_post, timeout = 30)
    
    if r_sp_loci.status_code == 409:
        print("Locus already exists on species")
        #return
    elif r_sp_loci.status_code == 404:
        print("species not found")
        #return
    
#    elif r_sp_loci.status_code > 201:
#        print("Server returned code " + str(r_sp_loci.status_code))
#        print("Retrying in seconds " + str(waitFactor))
#        time.sleep(waitFactor)
#        waitFactor = waitFactor * 2
    
#    else:
#        sucess_send = True
            
    
    loci_url = r_sp_loci.json()['url']
    
    return loci_url


def select_name(result):
    """
    """
    name = ''
    url = ''
    label = ''

    aux = result["results"]["bindings"]
    
#    try:
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
    

def get_data(sparql_queries):
    """ Retrieve annotations from uniprot
    """

    locus = sparql_queries[0]
    queries = [qq[0] for qq in sparql_queries[1]]
    dna_seq_to_ns = [dna[1] for dna in sparql_queries[1]]
    
    virtuoso_server.setReturnFormat(JSON)
    virtuoso_server.setTimeout(10)

    url = ''
    name = ''
    prev_name = ''
    label = ''
    found = False
    unpreferred_names = ['Uncharacterized protein', 'hypothetical protein', 'DUF']

    # implement more than 1 retry!!!
    alleles = len(queries)
    a = 0
    while found is False:
        virtuoso_server.setQuery(queries[a])

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
        

    return (locus, prev_name, label, url, dna_seq_to_ns)


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


def create_protein_files(response_dict, parent_dir):
    """ Creates fasta file with the protein sequence of
        the allele to be uploaded to NS.
    """
    
    protein_files = []
    dna_files = []
    invalid_prots = []
    for file, v in list(response_dict.items()):
        dna_files.append(file)
        indexed_fasta = SeqIO.index(file, 'fasta')
        prots = []
        locus_id = os.path.basename(file)
        locus_id = locus_id.split('.fasta')[0]
        for allele in v:
            allele_seq = indexed_fasta.get(allele[1]).seq
            protein = translate_sequence(str(allele_seq), 11)
            if isinstance(protein, TranslationError):
                invalid_prots.append(allele_seq)
            else:
                prots.append((allele[1], protein))
        
        protein_lines = []
        for prot in prots:
            header = '>{0}\n'.format(prot[0])
            sequence = '{0}\n'.format(prot[1])
            protein_lines.append(header)
            protein_lines.append(sequence)
        
        protein_lines = ''.join(protein_lines)
        protein_file = '{0}/{1}{2}'.format(parent_dir, locus_id, '_protein.fasta')
        with open(protein_file, 'w') as pf:
            pf.write(protein_lines)
        protein_files.append(protein_file)
        
    return protein_files, dna_files


def create_uniprot_queries(protein_files, dna_files):
    """ Build SPARQL queries for Uniprot
    """
    # construct queries
    queries = {}
    for pfile, dnafile in zip(protein_files, dna_files):
        pfile_id = os.path.basename(pfile).split('.fasta')[0]
        queries[pfile_id] = []
        #index dna file
        indexed_dna_fasta = SeqIO.index(dnafile, 'fasta')
        for record in SeqIO.parse(pfile, 'fasta'):
            # get the dna seq
            dna_seq = str(indexed_dna_fasta.get(record.id).seq)
            queries[pfile_id].append((uniprot_query(str(record.seq)), dna_seq))
            
    return queries
    


def post_alleles3(sequence, name, label, uniprot_url, loci_url, species_name, noCDSCheck, headers_post):
    """ Creates new alleles in the endpoint. Multiple alleles can be created.
    """    
    
    # Build the url for loci/loci_id/alleles
    url = ut.make_url(loci_url, "alleles")

    params = {}
    params['sequence'] = sequence
    params['species_name'] = species_name
    
    if not uniprot_url == '':
        params['uniprot_url'] = uniprot_url
        
    if not label == '':
        params['uniprot_label'] = label
        
    if not name == '':
        params['uniprot_sname'] = name
    
    if not noCDSCheck:
        params['enforceCDS'] = "False"
        
    params["input"] = "auto"
         
    r = requests.post(url, data = json.dumps(params), headers = headers_post, timeout = 30)

    return r
    

def post_species_loci(url, species_id, locus_id, headers_post):
    """
    """
    
    # Define POST request parameters
    params = {}
    params['locus_id'] = locus_id
    
    # Build the url for the loci of the new schema        
    url_species_loci = ut.make_url(url, "species", species_id, "loci")
    
#    req_success = False
#    sleepfactor = 4
#    while not req_success:
        
    r = requests.post(url_species_loci, data = json.dumps(params), headers = headers_post, timeout = 30)
    if r.status_code > 201:
        print("post_species_loci, " + r.status_code)
#        print("failed, retrying in seconds " + str(sleepfactor))
#            time.sleep(sleepfactor)
#            sleepfactor = sleepfactor * 2
#        else:
#            req_success = True
            
    return True

    


def post_schema_loci(loci_url, schema_url, headers_post):
    """
    """
    
    # Get the new loci id from the new loci url
    new_loci_id = str(int(loci_url.split("/")[-1]))
    
    # Define POST request parameters
    params = {}
    params['loci_id'] = new_loci_id
    
    # Build the url for the loci of the new schema        
    url_schema_loci = ut.make_url(schema_url, "loci")
    
#    req_success = False
#    sleepfactor = 4
#    while not req_success:
        
    r = requests.post(url_schema_loci, data = json.dumps(params), headers = headers_post, timeout = 30)
    if r.status_code > 201:
        print("post_schema_loci, "+ r.status_code)
#        print("failed, retrying in seconds " + str(sleepfactor))
#        time.sleep(sleepfactor)
#        sleepfactor = sleepfactor * 2
#    else:
#        req_success = True
            
    return True


def link_allele(sequence_uri, loci_url, species_name, headers_post):
    """ Links an existing allele uri (in NS) 
        to a newly created locus.
    """
    
    # Build the url for loci/loci_id/alleles
    url = ut.make_url(loci_url, "alleles")

    params = {}
    params['sequence_uri'] = sequence_uri
    params['species_name'] = species_name
            
    params["input"] = "link"
         
    r = requests.post(url, data = json.dumps(params), headers = headers_post, timeout = 30)

    return r


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, dest='input_files', required=True,
                        help='Path to folder containg the schema fasta files '
                        '(alternatively, a file with a list of paths to '
                        'fasta files).')

    parser.add_argument('-sp', type=str, dest='species_id', required=True,
                        help='The identifier for the schemas species in '
                        'the NS.')

    parser.add_argument('--sdesc', type=str, dest='schema_desc', required=True,
                        help='A brief description for the schema. It is '
                        'important to pass a meaningful value')

    parser.add_argument('--lprefix', type=str, dest='loci_prefix', required=True,
                        help='Prefix added to the identifier of each schema gene.'
                        ' For instance, ACIBA will produce ACIBA00001.fasta.')

    parser.add_argument('--cpu', type=int, required=False, dest='cpu_threads', 
                        default=1, help='The number of CPU threads to use '
                        '(default=1).')

    parser.add_argument('--url', type=str, required=False, dest='base_url',
                        default='http://127.0.0.1:5000/NS/api/',
                        help='The base URL for the NS server. Will be used '
                        'to create endpoints URLs.')

    parser.add_argument('--cont', type=bool, dest='continue_up', required=False,
                        default=False, help='Flag used to indicate if the process '
                        'should try to continue a schema upload that crashed '
                        '(default=False.')

    args = parser.parse_args()

    return [args.input_files, args.species_id, args.schema_desc,
            args.loci_prefix, args.cpu_threads, args.base_url,
            args.continue_up]


def main():
    
    # Get the parent directory of the schema directory
    parent_dir = os.path.dirname(input_files)

    # Check if the user provided a valid number of cpu
    cpu_threads = ut.verify_cpu_usage(cpu_threads)

    # Check if user provided a list of genes or a folder
    input_files = ut.check_if_list_or_folder(input_files)

    # Create list of genes and sort it by name
    with open(input_files, 'r') as gf:
        fasta_paths = [path.strip() for path in gf]
    fasta_paths.sort()

    # remove file with list of paths
    if os.path.isfile('listGenes.txt'):
        os.remove('listGenes.txt')

    print('Processing the fastas')

    # login with master key
    login_key = False
    if login_key:
        pass
    # if the login key is not found ask for credentials
    else:
        print('\nPlease provide login credentials:')
        user = input('USERNAME: ')
        password = getpass('PASSWORD: ')
        print()
        # get token
        token = ut.login_user_to_NS(base_url, user, password)

    # Define the headers of the requests
    headers_get = {'X-API-KEY': token,
                   'accept': 'application/json'}

    headers_post = {'X-API-KEY': token,
                    'Content-type': 'application/json',
                    'accept': 'application/json'}

    # verify user role to check permission
    user_info = simple_get_request(base_url, headers_get,
                                   ['user', 'current_user'])
    user_info = user_info.json()
    user_role = any(role in user_info['roles']
                    for role in ['Admin', 'Contributor'])

    if not user_role:
        print('Current user has no Administrator or Contributor permissions.\n'
              'Not allowed to upload schemas.')

        return 403

     # Get the name of the species from the provided id
    species_info = simple_get_request(base_url, headers_get,
                                  ['species', species_id])
    species_name = species_info.json()[0]['name']['value']

    # Create the parameters for the GET request
    params = {}
    params['name'] = schema_name
    params['bsr'] = schema_bsr
    params['ptf'] = schema_ptf 
    params['translation_table'] = schema_translation_table 
    params['min_locus_len'] = schema_min_locus_len 
    params['chewBBACA_version'] = schema_chewBBACA_version
  
    # Build the new schema URL and POST to NS
    print('Posting to /schemas endpoint...')
    schema_post = simple_post_request(base_url, headers_post,
                                      ['species', species_id, 'schemas'],
                                      params)

     # check status code
    # add other prints for cases that fail so that users see a print explaining
    schema_post_status = schema_post.status_code
    if schema_post_status in [200, 201]:
        message = ('A new schema for {0} was created '
                   'succesfully.'.format(species_name))
        print(message)
    else:
        if schema_post_status == 403:
            message = ('{0}: No permission to load '
                       'schema.'.format(schema_post_status))

        elif schema_post_status == 404:
            message = ('{0}: Cannot load a schema for species '
                       'that is not in NS.'.format(schema_post_status))

        elif schema_post_status == 409:
            message = ('{0}: Cannot load a schema with the same '
                       'description as a schema that is in the '
                       'NS.'.format(schema_post_status))

        else:
            message = '{0}: Could not insert schema.'.format(schema_post_status)

        print(message)
        
        return message
    
    # Get the new schema url from the response
    schema_url = r_schema.json()['url']
    
    start = time.time()
    # 1 LOCI PER FILE, MULTIPLE ALLELES PER FILE
#    conflict_loci_ids = []
    print('Checking if sequences in each FASTA file are in the NS...\n')   
#    response_200 = []
#    response_404 = {}
    fp = 0
    testing = []
    for gene in listGenes:        
        
        # Check if sequences are already on NS
        check_seq_responses = check_seq2(gene, species, NS_BASE_URL, headers_get, cpu2use)
        
        testing.append(check_seq_responses)
        
        fp += 1
        print("\r", "Files checked: {0}/{1}".format(fp, len(listGenes)), end='')
    
    end_check = time.time()
    delta_check = end_check - start
    print("\nCheck files", delta_check)
        
    
    for fiel in testing:
#        print(fiel.keys())
        response_404_file = {}
        response_200_file = []
        for k, v in fiel.items():
            response_404_file[k] = [res for res in v if res[0].status_code == 404]
            response_200_file = [res[0].json()["sequence_uri"] for res in v if res[0].status_code == 200]

            # Create a new loci
            new_loci_url = post_loci(NS_BASE_URL, species, headers_post, schema_prefix, keepFileName, os.path.basename(k))      
            
            # Get the new loci ID
            new_loci_id = new_loci_url.split("/")[-1]
                        
            if not ut.isListEmpty(response_200_file):
                # Link existing alleles
                a = 0
                print("Linking existing alleles to new schema...")
                for uri in response_200_file:
                    rr = link_allele(uri, new_loci_url, species_name, headers_post)
                    a += 1
                    print("\r", "Alleles linked {0}/{1}\n".format(a, len(response_200_file)), end='')
                    
                    
                # Associate the new loci id to the species
                new_species_loci = post_species_loci(NS_BASE_URL, species, new_loci_id, headers_post)
                
                # Associate the new loci id to the new schema
                new_schema_loci = post_schema_loci(new_loci_url, schema_url, headers_post)   

            
            #if any([v for v in response_404.values() if v != []]):
            if not bool([v for v in response_404_file.values() if v == []]):
                
                #novelty_loci = {k: v for k, v in response_404.items() if len(v) > 0}
                #print(novelty_loci)

                # Insert new alleles
                total_alleles = 0
                for i in response_404_file.values():
                    total_alleles += len(i)
                    
                print("\n Total alleles to insert:", str(total_alleles))
                
                # Create protein fasta files for uniprot queries    
                print("\nCreating protein files...\n")
                protein_files, dna_files = create_protein_files(response_404_file, parent_dir)
                
                # Create queries Uniprot
                print("Building Uniprot queries...\n")
                queries = create_uniprot_queries(protein_files, dna_files)
        
                # Remove protein files
                print("Removing protein files...\n")
                for pf in protein_files:
                    os.remove(pf)
            
                # Get the Uniprot annotations
                print("Retrieving Uniprot annotations...\n")
                results = []
                with concurrent.futures.ThreadPoolExecutor(max_workers=30) as executor:
                    # Start the load operations and mark each future with its URL
                    for res in executor.map(get_data, list(queries.items())):
                        results.append(res)
                        
                # Insert alleles and annotations
                       
                print("Inserting alleles in NS...\n")
                for r in results:
                    
                    loci = r[0]
                    name = r[1]
                    label = r[2]
                    url = r[3]
                    allele_seq_list = r[4]
                    
                    p = 0
                    for al in allele_seq_list:
                        # Create a new allele and add the sequences of the file
                        new_alleles = post_alleles3(al, name, label, url, new_loci_url, species_name, noCDSCheck, headers_post)
                        
                        p += 1
                        print("\r", "Alleles sent {0}/{1}\n".format(p, len(allele_seq_list)), end='')

                        
                    # Associate the new loci id to the species
                    new_species_loci = post_species_loci(NS_BASE_URL, species, new_loci_id, headers_post)
                    
                    # Associate the new loci id to the new schema
                    new_schema_loci = post_schema_loci(new_loci_url, schema_url, headers_post)   
                    
    
    print("\nDone")    
    end = time.time()
    delta = end - start

    
if __name__ == "__main__":
    main()        
    
    
