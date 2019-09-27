#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:54:54 2019

@author: pcerqueira
"""

import os
import time
import shutil
import requests
import concurrent.futures
from itertools import repeat
import extra_scripts.utils as ut
from extra_scripts import chewBBACA_PrepExternalSchema_optimization_tests as PrepExternalSchema


def build_fasta_files(response_dict, path2down):
    """ Write fasta files from get request responses
    
        Args:
            response_dict (dict): Contains the names as keys, get request response
            as values.
            
            path2down (str): Path to the download directory
            
        Results:
            failed_downloads (list): Contains the names 
    """

    failed_downloads = []

    for name, response in response_dict.items():
        
        #if the fasta already exist don't do anything
        if os.path.exists(os.path.join(path2down, name)):
            return True
        
        if response.status_code > 200:
            failed_downloads.append([name, response])
        
        res = response.json()
        res = res["Fasta"]
        auxDict = {}
        auxname = name.split(".")
        listIndexAux = []
        
        if os.path.isfile(os.path.join(path2down, name)):
            print (f"{name} already exists")
            return True

		#create a dictionary with the info from fasta to build the fasta file. {id:sequence}
        for allele in res:
            listIndexAux.append(int(allele['allele_id']['value']))
            auxDict[int(allele['allele_id']['value'])] = f"{allele['nucSeq']['value']}"

		# sort by allele id, create the fasta string and write the file
        listIndexAux.sort()
        auxString = ''
        for index in listIndexAux:
            auxString += f">{auxname[0]}_{index}\n{auxDict[index]}\n"
            
        with open(os.path.join(path2down, name), 'w') as f:
            f.write(auxString)
    
    print(f"{len(failed_downloads)} files were not downloaded\n")
    print(f"{len(response_dict) - len(failed_downloads)} files have been written")
    
    return failed_downloads


def get_fasta_seqs(headers_get, url):
    """ Perform get requests to the NS
    
        Args:
            headers_get (dict): headers for the GET request
            url (str): url to perform the request
            
        Return:
            res (response): response of the request containing
            fasta sequences
    """
    
    res = requests.get(url, headers = headers_get, timeout = 30)
    
    return res



def get_schema(schema_uri, path2down, cpu2use, maxBsrShort, headers_get):
    """ Downloads, builds a writes a schema from NS
    
        Args:
            schema_uri (str): uri of the schema to download
            path2down (str): Path to the download directory
            cpu2use (int): Number of cpu allowed
            maxBsrShort (float): Maximum BSR allowed for the representative selection
            headers_get (dict): headers for the GET request
    """

	#return error if folder to down schema doesnt exist        
    if not os.path.exists(path2down):
        os.mkdir(path2down)
       
    schemaLociList = ut.make_url(schema_uri, "loci")

	#get list of loci and build dictionary locus_id --> gene_name
    dictLoci = {}
    r = requests.get(schemaLociList, headers = headers_get)
    result = r.json()
    result = result["Loci"]
    serverTime = r.headers['Server-Date']


	#get server time and save on config before starting to down stuff, useful for further sync function
    if not os.path.exists(os.path.join(path2down, "config.py")):
        with open(os.path.join(path2down, ".config.txt"), 'w') as f:
            f.write(serverTime +'\n' + schema_uri)
       
    
    # Total number of loci
    total_len = len(result)
    
    print (f"Number of loci to Down: {str(total_len)}")
    for locus in result:
        dictLoci[str(locus['locus']['value'])] = f"{locus['name']['value']}.fasta"
    
    
    # build the list of urls to get
    ordered_keys = sorted(dictLoci.keys())
    
    fasta_urls = [ut.make_url(key, "fasta") for key in ordered_keys]    

    
    # Multithread the requests    
    print("Downloading schema...\n")
    responses = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers = 10) as executor:
        for result in executor.map(get_fasta_seqs, repeat(headers_get), fasta_urls):
            responses.append(result)
            
    # order dictLoci values by keys
    dictLoci_values_sorted = [dictLoci[key] for key in ordered_keys]
            
    # Build a dict with fasta_file ----> response
    res_dict = dict(zip(dictLoci_values_sorted, responses))
    
    # Build the fasta files 
    print("Writing fasta files...\n")
    fasta_files = build_fasta_files(res_dict, path2down)
    
    # Build list with path to fasta files for PrepExternalSchema
#    fasta_files = [os.path.join(path2down, fasta_file) for fasta_file in dictLoci.values()]
    
    # Select the representative sequences of the schema
    # calculate the BSR values and the short files

    # output dir for the result of PrepExternalSchema
    output = f"{path2down}_down"
    
    PrepExternalSchema.main(path2down, output, cpu2use, maxBsrShort)
    
    #copy .config.py to output
    shutil.copy(os.path.join(path2down, ".config.txt"), output)
    
    # remove path2down
    shutil.rmtree(path2down)
    
    # rename output to the name provided by the user (path2down)
    os.rename(output, path2down)
    
    print ("\n################\n" + str(len(ordered_keys)) + " loci downloaded")
    print ("Schema is now available at " + os.path.abspath(path2down))


def main(uri2schema, path2down, cpu2use, maxBsrShort):
    
#    uri2schema = "url_to_schema_in_NS"
#    
#    server_url =  uri2schema.split("/species")[0]
#
#    
#    path2down = "/path/to/schema/downloaded/from/NS"
#    cpu2use = 6
#    maxBsrShort = 0.6
    
    ## FOR DEBUG ONLY! ###########################################################################
    # Log the user in
    token = ut.login_user_to_NS(server_url, "test@refns.com", "mega_secret")
    #############################################################################################

    # Define the headers of the requests
    headers_get = {'X-API-KEY': token,
                   'accept' : 'application/json'}

    
    r = requests.get(uri2schema, headers = headers_get, timeout=5)
    if r.status_code > 200:
        print("server returned error with code " + str(
                r.status_code) + ", program will stop since there seems to be a problem contacting the server")
        return
    else:
        print("server is running, will start the program")
        
    #maxBsrShort=maxBsrShort

    
    start = time.time()
    
    schema_fasta = get_schema(uri2schema, path2down, cpu2use, maxBsrShort, headers_get)

    end = time.time()
    delta = end - start
    print(delta/60)



if __name__ == "__main__":
	main()
