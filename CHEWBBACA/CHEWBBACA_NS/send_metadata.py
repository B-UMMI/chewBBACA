#!/usr/bin/env python3

import requests
import multiprocessing
import time
import csv
from ast import literal_eval

try:
	from utils import CommonFastaFunctions
except ImportError:
	from CHEWBBACA_NS.utils import CommonFastaFunctions

try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO

def sanitize_input(mystring):
	print ("sanitizing")
	mystring=mystring.replace("'", "")
	mystring=mystring.encode('ascii', 'ignore')
	mystring=mystring.decode("utf-8")
	mystring=mystring.replace("\\", "")
	return mystring

def read_metadata(metadataFile,mustHaveHeaders,genomesURI=None):
	metadata={}
	with open(metadataFile) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')

		firstrow=next(reader)
		headers=set(firstrow)
		mustHaveHeaders=set(mustHaveHeaders)
		if not mustHaveHeaders==headers:
		#if "FILE" not in firstrow or "ACCESSION" not in firstrow or  "COUNTRY" not in firstrow or  "ST" not in firstrow or  "STRAIN" not in firstrow :
			print( "header not correct, columns to remove/add (use NA if no attribute to be set for the isolate): "+str(mustHaveHeaders.symmetric_difference(headers)))
			raise
		#metadata["headers"]=firstrow
		#print(genomesURI)

		for row in reader:
			#print(row[0])
			if genomesURI is None:
				metadata[row[0]] = row
			else:
				try:
					realURI=genomesURI[row[0]]
					metadata[realURI] = row
				except:
					#metadata[row[0]]=row
					pass
	#print(metadata)
	return 	metadata,firstrow

#send the request to add metadata, per genome
def send_metadata(token, headers,metadata,map_headers,url):

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
		#print(e)
		pass

	headers = {'Authentication-Token': token}

	if metadatacount<1:
		print("No metadata to upload")
		return
	else:
		print(url+" has "+ str(metadatacount)+ " metadata fields to add, sending to server...\n")

	sucess_send = False
	attempts = 0
	waitFactor = 4
	while not sucess_send and attempts < 5:
		r = requests.post(url, data=params, headers=headers, timeout=60)
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
			print(notuploaded)


	return "Done"



def main(metadataFile, cpu2Use,token,genomesURI=None):

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

	metadataFile,headers = read_metadata(metadataFile,necessaryHeadersList,genomesURI)

	print("metadata file read, continue")
	listUploadedGenomes=sorted(list(metadataFile.keys()))
	#listUploadedGenomes.remove("headers")

	#genome=(listUploadedGenomes[0]).replace('"', '')

	#print(listUploadedGenomes)
	pool = multiprocessing.Pool(cpu2Use)
	for genome in listUploadedGenomes:
		genome = genome.replace('"', '')
		p = pool.apply_async(send_metadata, args=[token, headers, metadataFile[genome],map_headers,genome])

	pool.close()
	pool.join()

if __name__ == "__main__":
	main()
