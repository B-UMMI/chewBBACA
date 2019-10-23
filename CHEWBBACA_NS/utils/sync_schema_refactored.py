#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:15:29 2019

@author: pcerqueira
"""


import os
import json
import time
import shutil
import random
import requests
import concurrent.futures
from itertools import repeat
import extra_scripts.utils as ut
from collections import defaultdict
from SPARQLWrapper import SPARQLWrapper, JSON
from extra_scripts import chewBBACA_PrepExternalSchema_optimization_tests as PrepExternalSchema

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

virtuoso_server = SPARQLWrapper('http://sparql.uniprot.org/sparql')


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
    


def syncSchema(lastSyncServerDate, schemaUri, path2schema, path_new_alleles, cpu2use, 
               bsrThresh, server_url, headers_get, headers_post, submit=False):

	#get list of sequences that are new considering the last date
    allNewSequences = False
    newDict = defaultdict(dict)
    servertime = lastSyncServerDate
    numberSequences = 0

	#get all sequences until the number of new sequences is less than 100k, the maximum the server return is 100k
    while not allNewSequences:
        newDict, allNewSequences, servertime, numberSequences = getNewsequences(newDict, servertime, schemaUri, numberSequences, headers_get)
    
    # If newDict is empty, schema is already updated
    if not bool(newDict):
        return 
    
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
    
    # submit option
    if submit:
    
        species_url = ut.make_url(server_url, "species", schemaUri.split("/")[-3])
        
        r_species_name = requests.get(species_url, headers = headers_get)
        
        species_name = r_species_name.json()[0]["name"]["value"]
    
        # loop over path2schema
        path2schema_files = [os.path.join(path2schema, fi) for fi in os.listdir(path2schema) if ".fasta" in fi]
       
        files_res = []
        # ask NS if they exist (check_seq2 from load_schema)
        for f in path2schema_files:
            # Check if sequences are already on NS
            check_seq_responses = check_seq2(f, server_url, headers_get, cpu2use)
                    
            files_res.append(check_seq_responses)
            
        for fiel in files_res:
            response_404_file = []
            response_200_file = []
            for k, v in fiel.items():
                # Get alleles that are not present in NS
                response_404_file = [res for res in v if res[0].status_code == 404]
                # Get alleles already present on NS
                response_200_file = [res for res in v if res[0].status_code == 200]
                
                # Proceed only if response_200_file list is not empty
                if not ut.isListEmpty(response_200_file):
                    for r in response_200_file:
                        
                        for locus in r[0].json()["result"]:
                            ns_locus_id = locus["locus"]["value"].split("/")[-1]
                        
                            # if it belongs to another locus in NS, add the allele to the current locus
                            if ns_locus_id not in os.path.basename(k):
                                
                                print("Detected allele belonging to another locus in NS, adding to the current locus...\n")
                                # get local locus id
                                local_locus_id_temp = k.split("-")[1].replace(".fasta", "")
                                if "_corrected" in local_locus_id_temp:
                                    local_locus_id = int(local_locus_id_temp.replace("_corrected", ""))
                                else:
                                    local_locus_id = int(local_locus_id_temp)
                                # get allele sequenece to upload
                                allele_seq = ut.get_sequence_from_url(r[0].url)
                                print("Translating allele...\n")
                                # translate sequence
                                prot_allele_seq = str(translate_sequence(allele_seq, 11))
                                # build uniprot query
                                unip_query = uniprot_query(prot_allele_seq)
                                print("Retrieving Uniprot annotation...\n")
                                # obtain uniprot results
                                name, label, url = get_data(unip_query)
                                # build loci url using local loci id
                                loci_url = ut.make_url(server_url, "loci", str(local_locus_id))
                                print("Adding to NS...\n")
                                # upload allele
                                new_alleles = post_alleles3(allele_seq, name, label, url, loci_url, species_name, headers_post)
                                
                if not ut.isListEmpty(response_404_file):
                    
                    print("Detected new alleles, adding them to NS...\n")
                    
                    # get local locus id
                    local_locus_id2_temp = os.path.basename(k).split("-")[1].replace(".fasta", "")
                    if "_corrected" in local_locus_id2_temp:
                        local_locus_id2 = int(local_locus_id2_temp.replace("_corrected", ""))
                    else:
                        local_locus_id2 = int(local_locus_id2_temp)
    
                    
                    for r2 in response_404_file:
                    
                        # get allele sequenece to upload
                        allele_seq2 = ut.get_sequence_from_url(r2[0].url)
                        print("Translating allele...\n")
                        # translate sequence
                        prot_allele_seq2 = str(translate_sequence(allele_seq2, 11))
                        # build uniprot query
                        unip_query2 = uniprot_query(prot_allele_seq2)
                        print("Retrieving Uniprot annotation...\n")
                        # obtain uniprot results
                        name2, label2, url2 = get_data(unip_query2)
                        # build loci url using local loci id
                        loci_url2 = ut.make_url(server_url, "loci", str(local_locus_id2))
                        print("Adding to NS...\n")
                        # upload allele
                        new_alleles2 = post_alleles3(allele_seq2, name2, label2, url2, loci_url2, species_name, headers_post)



    # if 404, nice, add the new allele to that locus
    # if 200, ask if they belong to the same locus on the NS and locally
    ### for example, locally its test_senterica15034336 and in NS its test_senterica15034336
    # if it belongs to the same locus do nothing
    # if it belongs to another locus in NS, add the allele to the current locus

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


def post_alleles3(sequence, name, label, uniprot_url, loci_url, species_name, headers_post, noCDSCheck=""):
    """ Creates new alleles in the endpoint. Multiple alleles can be created.
    """    
    
    # Build the url for loci/loci_id/alleles
    url = ut.make_url(loci_url, "alleles")
    
    # Build the request parameters
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


def main(path2schema, cpu2use, bsrThresh):
    
    path2schema = "/home/pcerqueira/Lab_Software/refactored_ns/ns_security/extra_scripts/test_down_senterica4"
    path_new_alleles = "/home/pcerqueira/Lab_Software/refactored_ns/ns_security/extra_scripts/senterica_new_alleles"
    cpu2use = 6
    bsrThresh = 0.6
    submit = True
    
#    blastPath = 'blastp'
    lastSyncServerDate, schemaUri = get_schema_vars(path2schema)
    
    server_url = schemaUri.split("/species")[0]

    
    ## FOR DEBUG ONLY! ###########################################################################
    # Log the user in
    token = ut.login_user_to_NS(server_url, "test@refns.com", "mega_secret")
    #############################################################################################

    # Define the headers of the requests
    headers_get = {'X-API-KEY': token,
                   'accept' : 'application/json'}
    
    headers_post = {'X-API-KEY': token,
                    'Content-type': 'application/json',
                    'accept' : 'application/json'}

    # Request to see if server is responding correctly
    r = requests.get(schemaUri, headers = headers_get, timeout=5)
    if r.status_code > 200:
        print("server returned error with code " + str(
                r.status_code) + ", program will stop since there seems to be a problem contacting the server")
        return
    else:
        print("server is running, will start the program")
        
    
    newserverTime = syncSchema(lastSyncServerDate, schemaUri, path2schema, path_new_alleles, 
                               cpu2use, bsrThresh, server_url, headers_get, headers_post, submit)
    
    
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

