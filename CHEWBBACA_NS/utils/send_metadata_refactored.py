#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:49:13 2019

@author: pcerqueira
"""

import csv
import json
import time
import requests
import multiprocessing
from ast import literal_eval

from extra_scripts import utils as ut
from app import app


def read_metadata(metadataFile, mustHaveHeaders, genomesURI=None):
    """
    """
    
    metadata={}
    with open(metadataFile) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')

        firstrow=next(reader)
        headers=set(firstrow)
        mustHaveHeaders=set(mustHaveHeaders)
        
        if not mustHaveHeaders == headers:
            print("header not correct, columns to remove/add (use NA if no attribute to be set for the isolate): " 
                     + str(mustHaveHeaders.symmetric_difference(headers)))
            raise

        for row in reader:
            if genomesURI is None:
                metadata[row[0]] = row
            else:
                try:
                    realURI=genomesURI[row[0]]
                    metadata[realURI] = row
                except:
                    pass

    return metadata, firstrow


def send_metadata(headers, metadata, map_headers, url, headers_post):

	i=0
	params = {}
	metadatacount = 0
	header=''
	
    #build the params to add to the request, if no params the request will not be send
	try:
		for elem in metadata:
			header=headers[i]
			if header=="FILE":
				i += 1
				continue
			param_name=map_headers[header]

			#check if the metadata cell is empty or with NA
			if not elem=='' and not elem.lower()=='na':
				params[param_name]=elem
				metadatacount += 1
			i+=1
	except Exception as e:
		pass


	if metadatacount<1:
		print("No metadata to upload")
		return
	else:
		print(url+" has "+ str(metadatacount)+ " metadata fields to add, sending to server...\n")

	sucess_send = False
	attempts = 0
	waitFactor = 4
	while not sucess_send and attempts < 5:
		r = requests.post(url, data=json.dumps(params), headers=headers_post, timeout=60)
		req_code = int(r.status_code)

		if req_code == 409:
			print(url)
			print("Genome already has all metadata")
			return

		elif req_code == 403:
			print(url)
			print("Genome is not your property")
			return

		elif req_code > 201:
			print((r.content).decode("utf-8"))
			print("Server returned code " + str(req_code))
			time.sleep(waitFactor)
			waitFactor = waitFactor * 2
			attempts += 1
		# raise
		else:
			sucess_send = True
			results=r.json()
			print(url)
			print("\nnumber of uploaded: " +str(results['Uploaded_total'][0]))
			notuploaded=results['Not_uploaded'][0]
			print("\nnot uploaded data:")
			print(str(notuploaded))


	return "Done"




def main(metadataFile, cpu2Use, genomesURI=None):
    
    server_url = app.config['BASE_URL']
    
    # Log user in
    token = ut.login_user_to_NS(server_url, "test@refns.com", "mega_secret")
    
    # Define the headers of the requests
    headers_get = {'X-API-KEY': token,
                   'accept' : 'application/json'}
    
    headers_post = {'X-API-KEY': token,
                    'Content-type': 'application/json',
                    'accept' : 'application/json'}


    #turn string to dict
    if not genomesURI is None:
        genomesURI = literal_eval(genomesURI)

    necessaryHeadersList=['FILE',
                          'ACCESSION',
                          'COUNTRY',
                          'ST',
                          'collection_date',
                          'host',
                          'host_disease',
                          'lat',
                          'long',
                          'isol_source',
                          'STRAIN'
                          ]

    #map headers from the tsv with the params name from the request
    map_headers={'FILE':'uri',
                          'ACCESSION':'accession',
                          'COUNTRY':'country',
                          'collection_date':'collection_date',
                          'host':'host',
                          'host_disease':'host_disease',
                          'lat':'lat',
                          'long':'long',
                          'isol_source':'isol_source',
                          'ST':'ST',
                          'STRAIN':'strainId'
    }

    metadataFile, headers = read_metadata(metadataFile, necessaryHeadersList, genomesURI)

    print("metadata file read, continue")
    listUploadedGenomes = sorted(list(metadataFile.keys()))
    
    for genome in listUploadedGenomes:
        genome = genome.replace('"', '')
        print(genome)
        
        send_metadata(headers, metadataFile[genome], map_headers, genome, headers_post)




#    pool = multiprocessing.Pool(cpu2Use)
#    for genome in listUploadedGenomes:
#        genome = genome.replace('"', '')
#        p = pool.apply_async(send_metadata, args=[token, headers, metadataFile[genome],map_headers,genome])
#
#    pool.close()
#    pool.join()

if __name__ == "__main__":
    main()

