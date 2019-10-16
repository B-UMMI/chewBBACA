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

from app import app

virtuoso_server = SPARQLWrapper('http://sparql.uniprot.org/sparql')


def check_seq_req(headers_get, url):
    """
    """
    
    url_req = url[0]
    
    res = requests.get(url_req, headers = headers_get, timeout = 30)
    
    return (res, url[1])
    

def check_seq2(fasta, species, URL, headers_get, cpu):
    """ Checks if a sequence already exists in NS
    """
    
    cpu = 30
    
    # Get the sequences of each record from the fasta file
    sequences = [(str(rec.seq), rec.id) for rec in SeqIO.parse(fasta, "fasta")]
    
    responses = {}
    
    responses[fasta] = []
    
    urls = []
    
    for seq in sequences:
        
        url_seq_info = ut.make_url(URL, "sequences", "seq_info", sequence=seq[0], species_id=species)
        urls.append((url_seq_info, seq[1]))
    
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=cpu) as executor:
        for result in executor.map(check_seq_req, repeat(headers_get), urls):
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
    
#sparql_queries = ('GCA-001223005-protein100_protein', queries['GCA-001223005-protein100_protein'])

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
    params['uniprot_url'] = uniprot_url
    params['uniprot_label'] = label
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




def main():
    
    parser = argparse.ArgumentParser(
            description="This program loads a schema to the nomenclature server, given the fasta files")
    parser.add_argument('-i', nargs='?', type=str, help='path to folder containg the schema fasta files ( alternative a list of fasta files)', required=True)
    parser.add_argument('-sp', nargs='?', type=str, help='species id', required=True)
    parser.add_argument('-t', nargs='?', type=str, help='token', required=True)
    parser.add_argument('--sname', nargs='?', type=str, help='schema name', required=True)
    parser.add_argument('--sprefix', nargs='?', type=str, help='loci prefix, for instance ACIBA will produce ACIBA00001.fasta', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument('--keep', help='store original fasta name too', required=False, default=False, action='store_true')
    parser.add_argument('--notCDS', help='dont enforce the sequences to be cds', required=False,default=False,action='store_false')
    parser.add_argument('--cont', help='use this flag to continue a schema upload that crashed in between', required=False, default=False, action='store_true')
    
    args = parser.parse_args()
    geneFiles = args.i
    species = args.sp
    token = args.t
    schema_name = args.sname
    schema_prefix = args.sprefix
    cpu2Use=args.cpu
    keepFileName=args.keep
    continue_previous_upload=args.cont
    noCDSCheck=args.notCDS
    
    #### DEBUG ONLY!!
    geneFiles = "/home/pcerqueira/Lab_Software/refactored_ns/ns_security/extra_scripts/adapted_senterica_20_2nd_samples"
    species = "13"
    schema_name = "test_debugging_script_senterica8"
    schema_prefix = "test_senterica8"
    cpu2Use = 6
    keepFileName = True
    continue_previous_upload = False
    noCDSCheck = False

    # Get the parent directory of the schema directory
    parent_dir = os.path.dirname(geneFiles)
    
    # Check if the user provided a valid number of cpu
    cpu2use = ut.verify_cpu_usage(cpu2Use)
    
    
    # Check if user provided a list of genes or a folder
    geneFiles = ut.check_if_list_or_folder(geneFiles)
	
	# Create list of genes and sort it by name
    listGenes = []
    with open(geneFiles, "r") as gf:
        for gene in gf:
            gene = gene.rstrip('\n')
            listGenes.append(gene)
    listGenes.sort()
    
    try:
        os.remove("listGenes.txt")
    except:
        pass
    
    print("Processing the fastas")
    
    
    # Obtain the base url from the app's config file
    URL = app.config['BASE_URL']
    
    ## FOR DEBUG ONLY! ###########################################################################
    # Log the user in
    token = ut.login_user_to_NS(URL, "test@refns.com", "mega_secret")    
    #############################################################################################

    # Define the headers of the requests
    headers_get = {'X-API-KEY': token,
                   'accept' : 'application/json'}
    
    headers_post = {'X-API-KEY': token,
                    'Content-type': 'application/json',
                    'accept' : 'application/json'}

    
    # Get the name of the species from the provided id
    species_url = ut.make_url(URL, "species", str(species))
    
    r_species_name = requests.get(species_url, headers = headers_get)
    
    species_name = r_species_name.json()[0]["name"]["value"]
        
    
    # Create the parameters for the GET request
    params = {}
    params['name'] = schema_name
    
    # ON A GET REQUEST DO NOT ADD CONTENT TYPE APPLICATION JSON
    # IT WILL TRY TO DECODE THE (EMPTY) BODY OF THE REQUEST AND
    # RETURN 400
    
    # Build the Schemas url
    url_schema = ut.make_url(URL, "species", species, "schemas")
    
    print("Posting to /schemas")
    
    # POST a new schema for the provided species ID  
    r_schema = requests.post(url_schema, data = json.dumps(params), headers = headers_post)

    if r_schema.status_code == 409:
        return "schema already exists"
    
    elif r_schema.status_code > 201:
        print(r_schema.status_code)
        return "something went wrong, species probably does not exist"
    
    # Get the new schema url from the response
    schema_url = r_schema.json()['url']
    
    
    start = time.time()
    # 1 LOCI PER FILE, MULTIPLE ALLELES PER FILE
    conflict_loci_ids = []
    print("Checking files...\n")   
    response_200 = []
    response_404 = {}
    fp = 0
    for gene in listGenes:        
        
        # Check if sequences are already on NS
        check_seq_responses = check_seq2(gene, species, URL, headers_get, cpu2use)

        # Sequence already exists, associate existing loci id with new schema 
        response_200.append([response for response in check_seq_responses[gene] if response[0].status_code == 200])
        # New sequence
        response_404[gene] = [response for response in check_seq_responses[gene] if response[0].status_code == 404]
        
        fp += 1
        print("\r", "Files checked: {0}/{1}".format(fp, len(listGenes)), end='')
    
    end_check = time.time()
    delta_check = end_check - start
    print("Check files", delta_check)
    
    ###### LINK EXISTING ALLELES ##################################################################################
    
    if not ut.isListEmpty(response_200):
                    
        # Get the responses' content
        response_contents = [res1[0][0].json() for res1 in response_200]
        response_contents_flat = ut.flatten_list(response_contents)
        
        # Obtain the loci_ids from the responses
        loci_ids = {res3["locus"]["value"].split("/")[-1] for res3 in response_contents_flat}
                    
        # Build the schema/loci url
        schema_loci_url = ut.make_url(schema_url, "loci")
              
        for loci_id in loci_ids:
            
            # Define the post data
            schema_loci_params = {}
            schema_loci_params["loci_id"] = str(loci_id)
                
            # Associate existing loci ID with new schema
            r_schema_loci = requests.post(schema_loci_url, data = json.dumps(schema_loci_params), headers = headers_post)
                
            if r_schema_loci.status_code > 201:
                print("Error occurred", r_schema_loci.status_code)
                
            if r_schema_loci.status_code == 409:
                print(gene)
                conflict_loci_ids.append(r_schema_loci.json())
                
    ###### INSERT NEW ALLELES #####################################################################################
    
    # Check if response_404's values are empty lists
    if not bool([v for v in response_404.values() if v == []]):
        
        total_alleles = 0
        for i in response_404.values():
            total_alleles += len(i)
    
        print("\n Total alleles to insert: ", str(total_alleles))

    
        # Create protein fasta files for uniprot queries    
        print("\nCreating protein files...\n")
        protein_files, dna_files = create_protein_files(response_404, parent_dir)
            
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
        
        end_uniprot = time.time()
        delta_uniprot = end_uniprot - start 
        print("Uniprot annotations", delta_uniprot)
        
        # Insert alleles and annotations
        p = 0        
        print("Inserting alleles in NS...\n")
        for r in results:
            
            loci = r[0]
            name = r[1]
            label = r[2]
            url = r[3]
            allele_seq_list = r[4]
            
            # Create a new loci
            loci_id = loci.replace("_protein", ".fasta")
            new_loci_url = post_loci(URL, species, headers_post, schema_prefix, keepFileName, loci_id)
            
            # Get the new loci ID
            new_loci_id = new_loci_url.split("/")[-1]
            
            for al in allele_seq_list:
                # Create a new allele and add the sequences of the file
                new_alleles = post_alleles3(al, name, label, url, new_loci_url, species_name, noCDSCheck, headers_post)
                
            # Associate the new loci id to the species
            new_species_loci = post_species_loci(URL, species, new_loci_id, headers_post)
            
            # Associate the new loci id to the new schema
            new_schema_loci = post_schema_loci(new_loci_url, schema_url, headers_post)   
            
            p += 1
            print("\r", "Loci inserted {0}/{1}".format(p, len(results)), end='')
        
    print("\nDone")    
    end = time.time()
    delta = end - start
    


if __name__ == "__main__":
    main()        
    
    
