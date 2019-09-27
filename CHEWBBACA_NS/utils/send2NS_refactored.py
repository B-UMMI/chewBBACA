#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:45:07 2019

@author: pcerqueira
"""

import os
import csv
import json
import time
import requests
from collections import defaultdict
import extra_scripts.utils as ut

def process_genome(genome, schemaURI, profile, listHeaders, headers_post):

    # schema uri is like http://137.205.69.51/app/NS/species/1/schemas/1 remove last 2 and get species uri
    server_ip = ut.make_url(schemaURI.split("/schemas")[0], "profiles")
#    server_ip.pop()
#    server_ip.pop()
#    server_ip = (str.join("/", server_ip)) + "/profiles"
#    server_ip = "http://127.0.0.1:5000/NS/api/species/13/profiles"

    event_data = {}
    event_data['profile'] = {genome: profile}
    event_data['headers'] = listHeaders
    
#    event_data['profile'] = {ola[0]: ola2}
    
#    r = requests.post(server_ip, headers = headers_post, data=json.dumps(event_data))
    

    # ~ server_return = requests.post(server_ip, headers=headers, json=event_data)
    
    
    server_return = requests.post(server_ip, headers=headers_post, data=json.dumps(event_data))
    
    if server_return.status_code == 409:
        print(genome + " already exists on NS, continuing...")
        return 

    
    
#    sucess_send = False
#    waitFactor = 4
#    attempts=0
#    while not sucess_send and attempts<5:
#
#        server_return = requests.post(server_ip, headers=headers_post, data=json.dumps(event_data))
#
##        if server_return.status_code>201:
##            time.sleep(waitFactor)
##            waitFactor=waitFactor*2
##            attempts+=1
#        
#        if server_return.status_code == 409:
#            
#
#        else:
#            sucess_send=True

    print(server_return.content.decode("utf-8"))

    try:
        time.sleep(1)
        genome_new_uri=((server_return.content.decode("utf-8").split("at "))[-1]).strip()
    except:
        genome_new_uri=False

    return [genome, genome_new_uri]



def collect_new_alleles_seq(profileFile, pathSchema, schemaUri, cpu2Use, 
                            percentMDallowed, metadataFile, headers_get, headers_post):

    #read profile to upload and put into a dictionary {genome:[full profile]}
    profileDict2={}
    with open(profileFile) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            profileDict2[row[0]] = row

    #build dictionary with {locus:[alleles per locus]}
    profileDict = defaultdict(list)
    listHeaders = profileDict2.pop("FILE")
    genomes = profileDict2.keys()
    i = 0
    while i < len(listHeaders):

        for genome in genomes:
            profileDict[listHeaders[i]].append(profileDict2[genome][i])
        i += 1


    dict_new_alleles = defaultdict(list)


    listGenomes2discard = []

    #re-read the profile and count the missing data per genome
    with open(profileFile) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        #~ headers = reader.next()
        headers = next(reader)

        for row in reader:

            missing_data_count = 0
            missing_data_count += row.count("LNF")
            missing_data_count += row.count("NIPH")
            missing_data_count += row.count("PLOT")
            missing_data_count += row.count("ASM")
            missing_data_count += row.count("ALM")

            i = 1

            #if number of missing data is larger than the missing data allowed the genome is discarded
            if int(len(listHeaders)*(percentMDallowed)) <= missing_data_count:
                i += 1
                print (row[0] + " discarded from further analysis because missing data count is over " + str(percentMDallowed))
                listGenomes2discard.append(row[0])
                continue

#            hasAlleles = False

            #if number of missing data is acceptable, get the alleles that are only local and store them on a dictionary {locus:[new alleles]}
            while i<len(row):

                gene = headers[i]
                if "*" in row[i]:
                    new_allele= int(row[i].replace('*',''))
                    dict_new_alleles[headers[i]].append(new_allele)

                i+=1

    lists2print = []
    modifiedProfileDict = {}
    i = 0

    #create the string of the new profile
    while i < len(profileDict["FILE"]):
        aux=[]
        for header in listHeaders:
            aux.append(profileDict[header][i])
        i += 1
        modifiedProfileDict[aux.pop(0)] = aux
        lists2print.append(aux)

    newProfileStr=('\t'.join(map(str,listHeaders))) + "\n"
    i = 0
    for elem in lists2print:
        newProfileStr += str(list(genomes)[i]) + "\t" + ('\t'.join(map(str, elem))) + "\n"
        i += 1
        
    for genome in sorted(modifiedProfileDict):

        if genome in listGenomes2discard:
            print (genome + " profile has a low locus count, probably a different species")
            continue
        
        else:
            profile = modifiedProfileDict[genome]
            ola = process_genome(genome, schemaUri, profile, listHeaders, headers_post)
                
    print ("DONE")





def get_schema_vars(schemapath):

    #get schema .config file and its variables, schema uri and last server data
    schemapath_config=os.path.join(schemapath,".config.txt")
    try:
        with open(schemapath_config) as fp:
            lines=[]
            for line in fp:
                lines.append(line.strip())

    except:
        print ("your schema has no config file")
        raise

    lastSyncServerDate=lines[0]
    schemaUri=lines[1]

    return lastSyncServerDate,schemaUri


def main(profileFile,pathSchema,token,metadataFile,percentMDallowed,cpu):
    
    cpu2Use=int(cpu)
    percentMDallowed = 0.5
    percentMDallowed=float(percentMDallowed)
    # profileFile = "/home/pcerqueira/Lab_Software/refactored_ns/ns_security/extra_scripts/toxic_schema_test/results_alleles_teste.tsv"
    # schemaUri = "http://127.0.0.1:5000/NS/api/species/13/schemas/6"

    #get schema uri and date from schema .config file
    lastSyncServerDate, schemaUri = get_schema_vars(pathSchema)
    
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
    
    collect_new_alleles_seq(profileFile, pathSchema, schemaUri, cpu2Use, 
                            percentMDallowed, metadataFile, headers_get, headers_post)



if __name__ == "__main__":
    main()

