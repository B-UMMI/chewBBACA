#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:20:29 2019

@author: pcerqueira
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:56:24 2019

@author: pcerqueira
"""

import os
import json
import argparse
import requests
import concurrent.futures
from itertools import repeat
import extra_scripts.utils as ut

from Bio import SeqIO
from app import app


def check_seq(fasta, species, URL, headers_get):
    """ Checks if a sequence already exists in NS
    """
    
    # Get the sequences of each record from the fasta file
    sequences = [str(rec.seq) for rec in SeqIO.parse(fasta, "fasta")]
    
    responses = {}
    
    #responses[fasta] = {}
    
    for seq in sequences:
        
#        allele_counter = 1
    
        # Build the url for seq_info endpoint, to check if the sequence
        # is already attributed to a locus
        url_seq_info = ut.make_url(URL, "sequences", "seq_info", sequence=seq, species_id=species)
        
        # Execute the GET request
        res = requests.get(url_seq_info, headers = headers_get, timeout = 30)
        
        # Save the response
        responses[fasta] = res
            
    return responses


def check_seq_req(headers_get, url):
    """
    """
    
    res = requests.get(url, headers = headers_get, timeout = 30)
    
    return res
    
    
    

def check_seq2(fasta, species, URL, headers_get, cpu):
    """ Checks if a sequence already exists in NS
    """
    
    # Get the sequences of each record from the fasta file
    sequences = [str(rec.seq) for rec in SeqIO.parse(fasta, "fasta")]
    
    responses = {}
    
    responses[fasta] = []
    
    urls = []
    
    for seq in sequences:
        
        url_seq_info = ut.make_url(URL, "sequences", "seq_info", sequence=seq, species_id=species)
        urls.append(url_seq_info)
    
    
    with concurrent.futures.ThreadPoolExecutor(max_workers = cpu) as executor:
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
        params['locus_ori_name'] = os.path.basename(gene)
        
    
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


def post_alleles(gene, token, loci_url, species_name, noCDSCheck, headers_post):
    """ Creates new alleles in the endpoint. Multiple alleles can be created.
    """
    
    # Get the sequences of each record from the fasta file
    sequences = [str(rec.seq) for rec in SeqIO.parse(gene, "fasta")]
    
    # Build the url for loci/loci_id/alleles
    url = ut.make_url(loci_url, "alleles")

    for allele_seq in sequences:

        params = {}
        params['sequence'] = allele_seq
        params['species_name'] = species_name
    
        if not noCDSCheck:
            params['enforceCDS'] = "False"
        
#        headers = {'X-API-KEY': token,
#                   'Content-type': 'application/json',
#                   'accept': 'application/json'}
        
#        req_success = False
#        sleepfactor = 4
#        while not req_success:
#            try:

        r = requests.post(url, data = json.dumps(params), headers = headers_post, timeout = 30)
        
        if r.status_code > 201:
            print(r)
            print("failed sending sequence, retrying in seconds ")
#                  + str(sleepfactor))
#            time.sleep(sleepfactor)
#            sleepfactor = sleepfactor * 2
#                else:
#                    req_success = True
#            except:
#                time.sleep(sleepfactor)
#                sleepfactor = sleepfactor * 2
#                pass
    
        req_code = int(r.status_code)
        # allele_url=((r.content).decode("utf-8")).replace('"', '').strip()
    
    return req_code



def post_alleles2(file, loci_url, species_name, noCDSCheck, headers_post):
    """ Creates new alleles in the endpoint. Multiple alleles can be created.
    """
    
    # Get the sequences of each record from the fasta file
    sequences = [str(rec.seq) for rec in SeqIO.parse(file, "fasta")]
    
    # Build the url for loci/loci_id/alleles
    url = ut.make_url(loci_url, "alleles")

    for allele_seq in sequences:

        params = {}
        params['sequence'] = allele_seq
        params['species_name'] = species_name
    
        if not noCDSCheck:
            params['enforceCDS'] = "False"

     
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
        print(r.status_code)
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
        print(r.status_code)
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
    geneFiles = "/home/pcerqueira/Lab_Software/chewie_NS/chewBBACA/CHEWBBACA_NS/schema_pyogenes_test"
    species = "13"
    schema_name = "test_pyogenes_sync"
    schema_prefix = "test_sync_pyogenes"
    cpu2Use = 6
    keepFileName = True
    continue_previous_upload = False
    noCDSCheck = False

    
    
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
    
    auth_params = {}
    auth_params["email"] = "test@refns.com" 
    auth_params["password"] = "mega_secret"
    
    auth_headers = {}
    auth_headers["Content-Type"] = "application/json"
    auth_headers["accepts"] = "application/json"
    
    auth_url = URL + "auth/login"
    
    auth_r = requests.post(auth_url, data=json.dumps(auth_params), headers=auth_headers)
    
    auth_result = auth_r.json() 
    
    token = auth_result["X-API-KEY"]
    
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


    
    # 1 LOCI PER FILE, MULTIPLE ALLELES PER FILE
    print("Looping through files...\n")
    for gene in listGenes:
        print(gene)
        
        # Check if sequences are already on NS
        print("Checking sequences...")
#        check_seq_responses = check_seq(gene, species, URL, headers_get)
        check_seq_responses = check_seq2(gene, species, URL, headers_get, cpu2use)
        
        # Responses are alleles
        for response in check_seq_responses[gene]:


            # Sequence already exists, associate existing loci id with new schema
            if response.status_code == 200:
                
                print("Sequence already exists, assigning existing loci_id to new schema...\n ")
                
                # Get the response content
                response_content = response.json()
                
                # Get ID of first result (may need to add ORDER BY to SPARQL query)
                loci_id = response_content[0]["locus"]["value"].split("/")[-1]
                
                # Build the schema/loci url
                schema_loci_url = ut.make_url(schema_url, "loci")
                
                # Define the post data
                schema_loci_params = {}
                schema_loci_params["loci_id"] = str(loci_id)
                
                # Associate existing loci ID with new schema
                r_schema_loci = requests.post(schema_loci_url, data = json.dumps(schema_loci_params), headers = headers_post)
                
                if r_schema_loci.status_code > 201:
                    print("Error occurred", r_schema_loci.status_code)
                    #return
            
            # New sequence
            if response.status_code == 404:
                
                print("New sequences detected, adding new locus...\n ")
                
                # Create a new loci
                new_loci_url = post_loci(URL, species, headers_post, schema_prefix, keepFileName, gene)
                
                # Get the new loci ID
                new_loci_id = new_loci_url.split("/")[-1]
                
                # Create a new allele and add the sequences of the file
                new_alleles = post_alleles2(gene, new_loci_url, species_name, noCDSCheck, headers_post)
                
                # Associate the new loci id to the species
                new_species_loci = post_species_loci(URL, species, new_loci_id, headers_post)
                
                # Associate the new loci id to the new schema
                new_schema_loci = post_schema_loci(new_loci_url, schema_url, headers_post)   
    
    print("Done")
        

if __name__ == "__main__":
    main()        
    
    