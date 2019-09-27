#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:15:29 2019

@author: pcerqueira
"""


import os
import shutil
import requests
import concurrent.futures
from itertools import repeat
import extra_scripts.utils as ut
from collections import defaultdict
from extra_scripts import chewBBACA_PrepExternalSchema_optimization_tests as PrepExternalSchema

from Bio import SeqIO





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


def getNewsequences(newDict, lastSyncServerDate, schemaUri, numberSequences, headers_get):
    
    print("Start trying to get new sequences from server")
	
    #request the new alleles starting on the date given
    params = {}
    params['date'] = lastSyncServerDate
    uri = ut.make_url(schemaUri, "loci", date=lastSyncServerDate)	
    # Request the new alleles
    r = requests.get(uri, headers = headers_get)
    
    result = r.json()
    result = result["newAlleles"]

	#if no Last-Allele on headers is because there are no alleles
    try:
        servertime = r.headers['Last-Allele']
    except:
        servertime = r.headers['Server-Date']
        print("The schema is already up to date at: " + str(lastSyncServerDate))
        return (newDict, True, servertime, numberSequences)
    
    print("Getting NS new alleles since " + str(lastSyncServerDate))



	#for each new allele add info to the dictionary
    for newAllele in result:
        
        geneFasta = (newAllele['locus_name']['value'])+".fasta"
        alleleid = newAllele['allele_id']['value']
        sequenceUri = newAllele['sequence']['value']
        newDict[geneFasta][alleleid] = sequenceUri
        
    numberSequences += len(result)
	
    #if the result is the maximum number of new alleles the server can return, re-do the function, else stop doing the function
    if r.headers['All-Alleles-Returned'] == 'False':
        return (newDict, False, servertime, numberSequences)
    else:
        return (newDict, True, servertime, numberSequences)



def get_new_allele_seqs(headers_get, uri):
    """
    """
    
    seq_hash = uri.split("/sequences/")[1]
    
    url = ut.make_url(uri.split(seq_hash)[0], "seq_info", seq_id=seq_hash)
    
    r = requests.get(url, headers = headers_get, timeout=30)
    
    r_content = r.json()
    
    allele_seq = r_content["sequence"]["value"]
    
    return allele_seq
    


def syncSchema(lastSyncServerDate, schemaUri, path2schema, path_new_alleles, cpu2use, bsrThresh, headers_get):

	#get list of sequences that are new considering the last date
    allNewSequences = False
    newDict = defaultdict(dict)
    servertime = lastSyncServerDate
    numberSequences = 0

	#get all sequences until the number of new sequences is less than 100k, the maximum the server return is 100k
    while not allNewSequences:
        newDict, allNewSequences, servertime, numberSequences = getNewsequences(newDict, servertime, schemaUri, numberSequences, headers_get)
        
    fasta_response = {}
    
    for fasta_file_name in newDict.keys():
        
        allele_dict = newDict[fasta_file_name]
        sequence_uri_list = list(allele_dict.values())
        
        # Multithread requests    
        with concurrent.futures.ThreadPoolExecutor(max_workers = 10) as executor:
            for result in executor.map(get_new_allele_seqs, repeat(headers_get), sequence_uri_list):
                fasta_response[fasta_file_name] = dict(zip(list(allele_dict.keys()), [result]))
    
    
    # Create an intermediate dir for the new alleles
    if not os.path.exists(path_new_alleles):
        os.mkdir(path_new_alleles)
              
    # Check if new alleles are already on schema  
    schema_fasta_files = {file for file in os.listdir(path2schema) if file not in ["short", ".config.txt"]}

    #fasta_response2 = copy.deepcopy(fasta_response)
    for res in fasta_response:
        
        if res in schema_fasta_files:
            print("Locus exists on schema, checking if allele exists...")
            
            #os.path.join(path2schema, res)
            #remove the * from the local alleles that were on the new alleles synced
            with open(os.path.join(path2schema, res), "r") as test, open(os.path.join(path_new_alleles, os.path.splitext(res)[0] + "_corrected.fasta"), "w") as corrected:
                records = SeqIO.parse(test, "fasta")
                for record in records:
                    if str(record.seq) in list(fasta_response[res].values()) and "*" in record.id:
                            record.id = record.id.replace("*", "")
                            record.description = record.id.replace("*", "")  
                            del fasta_response[res][record.id.split("_")[-1]]                                                       
                    SeqIO.write(record, corrected, "fasta")
    
    
    # Write new alleles in existing schema files
    keys_to_remove = []             
    for res in fasta_response:
        
        if res in schema_fasta_files:
            
            with open(os.path.join(path_new_alleles, os.path.splitext(res)[0] + "_corrected.fasta"), "a") as f:
                for k, v in fasta_response[res].items():
                    fasta_str = f">{os.path.splitext(res)[0]}_{k}\n{v}"
                    f.write(fasta_str)
        keys_to_remove.append(res)

    
    # remove already existing loci files
    fasta_response_final = {key: fasta_response[key] for key in fasta_response if key not in keys_to_remove}
           

    if not fasta_response_final == {}:
        # Create the fasta files for the new alleles
        build_fasta_files(fasta_response_final, path2schema, path_new_alleles)                
    
    # Re-determine the representative sequences    
    PrepExternalSchema.main(path_new_alleles, path2schema, cpu2use, bsrThresh)
        
    # remove the intermediate folder of the new alleles
    shutil.rmtree(path_new_alleles)

    return servertime


def get_schema_vars(schemapath):

	schemapath_config = os.path.join(schemapath,".config.txt")
	try:
		with open(schemapath_config) as fp:
			lines = []
			for line in fp:
				lines.append(line.strip())

	except:
		print ("your schema has no config file")
		raise

	lastSyncServerDate = lines[0]
	schemaUri = lines[1]

	return lastSyncServerDate,schemaUri


def main(path2schema, cpu2use, bsrThresh):
    
    #path2schema = "/path/to/schema/downloaded/from/NS"
    #path_new_alleles = "/path/to/intermediate/dir"
    cpu2use = 6
    bsrThresh = 0.6
    
    lastSyncServerDate, schemaUri = get_schema_vars(path2schema)
    
    server_url = schemaUri.split("/species")[0]

    
    ## FOR DEBUG ONLY! ###########################################################################
    # Log the user in
    token = ut.login_user_to_NS(server_url, "test@refns.com", "mega_secret")
    #############################################################################################

    # Define the headers of the requests
    headers_get = {'X-API-KEY': token,
                   'accept' : 'application/json'}


    # Request to see if server is responding correctly
    r = requests.get(schemaUri, headers = headers_get, timeout=5)
    if r.status_code > 200:
        print("server returned error with code " + str(
                r.status_code) + ", program will stop since there seems to be a problem contacting the server")
        return
    else:
        print("server is running, will start the program")
        
    
    newserverTime = syncSchema(lastSyncServerDate, schemaUri, path2schema, path_new_alleles, cpu2use, bsrThresh, headers_get)
    
    if newserverTime == lastSyncServerDate:
        return
    
    schemapath_config = os.path.join(path2schema, ".config.txt")
    try:
        with open(schemapath_config, "w") as fp:
            fp.write(newserverTime + "\n" + schemaUri)
    except:
        print ("your schema has no config file")
        raise


if __name__ == "__main__":
	main()

