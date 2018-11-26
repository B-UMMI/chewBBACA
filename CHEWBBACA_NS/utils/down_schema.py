#!/usr/bin/env python3

import requests
import os
import multiprocessing
import time

try:
	from utils import init_schema_4_bbaca
except ImportError:
	from CHEWBBACA_NS.utils import init_schema_4_bbaca

class Result():
	def __init__(self):
		self.val = []

	def update_result(self, val):
		lala = val
		self.val.append(lala)

	def get_result(self):
		return self.val


def down_fasta(uri,name,path2down,auxBar):

	#if the fasta already exist don't do anything
	if os.path.exists(os.path.join(path2down, name)):
		return True

	# request get fasta data for locus
	uriFasta=uri+"/fasta"

	sucess_send = False
	waitFactor = 4
	attempts=0
	while not sucess_send and attempts<5:
		r = requests.get(uriFasta, timeout=10)
		if r.status_code > 201:
			time.sleep(waitFactor)
			waitFactor=waitFactor*2
			attempts+=1
			print("Couldn't download the fasta from " + uri)
			print("Retrying in second " + str(waitFactor))
		else:
			sucess_send=True


	# get request was sucessfull, build the fasta file
	try:
		result=r.json()
		result = result["Fasta"]
		auxDict={}
		auxname=name.split(".")
		#auxString={}
		listIndexAux=[]

		if os.path.isfile(os.path.join(path2down, name)):
			print (name+ "already exists")
			return True

		#create a dictionary with the info from fasta to build the fasta file; {id:sequence}
		for allele in result:
			listIndexAux.append(int(allele['allele_id']['value']))
			auxDict[int(allele['allele_id']['value'])]=str(allele['nucSeq']['value'])

		# sort by allele id, create the fasta string and write the file
		listIndexAux.sort()
		auxString=''
		for index in listIndexAux:
			auxString+=">"+auxname[0]+"_"+str(index)+"\n"+str(auxDict[index])+"\n"

		with open(os.path.join(path2down, name), 'w') as f:
			f.write(auxString)

		if uri in auxBar:
			auxlen=len(auxBar)
			index=auxBar.index(uri)
			print ( "["+"="*index+">"+" "*(auxlen-index)+"] Downloading fastas "+str(int((index/auxlen)*100))+"%")


		return True

	except Exception as e:
		print (e)
		return uri

def get_schema(schema_uri,path2down,cpu2use,maxBsrShort):

	#return error if folder to down schema doesnt exist
	if not os.path.exists(path2down):
		return "create the dir"

	schemaLociList=schema_uri+"/loci"

	#get list of loci and build dictionary locus_id --> gene_name
	dictLoci={}
	r = requests.get(schemaLociList)
	result=r.json()
	result = result["Loci"]
	serverTime = r.headers['Server-Date']


	#get server time and save on config before starting to down stuff, usefull for further sync function
	if not os.path.exists(os.path.join(path2down, "config.py")):
		with open(os.path.join(path2down, ".config.txt"), 'w') as f:
			f.write(serverTime+'\n'+schema_uri)

	totallen=len(result)

	print ("Number of loci to Down: "+str(totallen))
	for locus in result:
		dictLoci[str(locus['locus']['value'])]=str(locus['name']['value'])+".fasta"

	#test down bar
	auxBar=[]
	orderedkeys=sorted(dictLoci.keys())
	step=int((len(orderedkeys))/10)+1
	counter=0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter+=step

	new_result = Result()

	#get the sequences per locus and build the fasta with it
	pool = multiprocessing.Pool(cpu2use)
	for uri in orderedkeys:

		#~ p=pool.apply_async(down_fasta,args=[uri,name,path2down])
		p=pool.apply_async(down_fasta,args=[uri,dictLoci[uri],path2down,auxBar],callback=new_result.update_result)
	#~ result=p.get()
	#~ counter.append(result)
	#~ print ("Downloaded "+str(len(counter))+"/"+str(totallen))

	pool.close()
	pool.join()

	listNotGotFastas = new_result.get_result()
	listNotGotFastas=list(set(listNotGotFastas))
	if len(listNotGotFastas)>1:
		print("some fastas were not properly downloaded")
		raise

	#~ print (str(len(counter))+" locus downloaded successfully") 

	print ("processing the fastas for chewBBACA usability")

	#calculate the BSR values and the short files, takes a lot of time!!!111
	init_schema_4_bbaca.main(path2down, cpu2use,maxBsrShort)

	print ("\n################\n" +str(len(orderedkeys))+" loci downloaded")
	print ("Schema is now available at "+os.path.abspath(path2down))

def main(uri2schema,path2down,cpu2use,maxBsrShort):

	server_url = (uri2schema.split("/species"))[0]

	r = requests.get(server_url, timeout=5)
	if r.status_code > 200:
		print("server returned error with code " + str(
			r.status_code) + ", program will stop since there seems to be a problem contacting the server")
		return
	else:
		print("server is running, will start the program")

	maxBsrShort=maxBsrShort

	lala=get_schema(uri2schema,path2down,cpu2use,maxBsrShort)




if __name__ == "__main__":
	main()
