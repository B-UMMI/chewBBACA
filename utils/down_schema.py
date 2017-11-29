#!/usr/bin/env python3

import requests,json
import csv
import argparse
from collections import defaultdict
import os
import subprocess
import multiprocessing

def down_fasta(uri,name,path2down,auxBar):
	
	if os.path.exists(os.path.join(path2down, name)):
		return True
	
	uriFasta=uri+"/fasta"
	r = requests.get(uriFasta, timeout=10)
	result=r.json()
	auxDict={}
	auxname=name.split(".")
	auxString={}
	listIndexAux=[]
	
	if os.path.isfile(os.path.join(path2down, name)):
		print (name+ "already exists")
		return True
	
	for allele in result:
		listIndexAux.append(int(allele['allele_id']['value']))
		auxDict[int(allele['allele_id']['value'])]=str(allele['nucSeq']['value'])
	
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

def get_schema(schema_uri,path2down,cpu2use):
	
	
	if not os.path.exists(path2down):
		return "create the dir"
	
	schemaLociList=schema_uri+"/loci"
	
	#build dictionary locus_id --> gene_name
	dictLoci={}
	r = requests.get(schemaLociList)
	result=r.json()
	
	serverTime=result.pop()
	#get server time before starting to down stuff, usefull for further sync function
	if not os.path.exists(os.path.join(path2down, "config.py")):
		with open(os.path.join(path2down, ".config.txt"), 'w') as f:
			f.write(serverTime+'\n'+schema_uri)

	totallen=len(result)
	
	print ("Number of loci to Down: "+str(totallen))
	for locus in result:
		dictLoci[str(locus['locus']['value'])]=str(locus['name']['value'])
	
	#test down bar
	auxBar=[]
	orderedkeys=sorted(dictLoci.keys())
	step=int((len(orderedkeys))/10)+1
	counter=0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter+=step

	
	#get the sequences per locus and build the fasta with it
	pool = multiprocessing.Pool(cpu2use)
	for uri in orderedkeys:
				
		#~ p=pool.apply_async(down_fasta,args=[uri,name,path2down])
		p=pool.apply_async(down_fasta,args=[uri,dictLoci[uri],path2down,auxBar])
		#~ result=p.get()
		#~ counter.append(result)
		#~ print ("Downloaded "+str(len(counter))+"/"+str(totallen))
	
	pool.close()
	pool.join()
				


	#~ print (str(len(counter))+" locus downloaded successfully") 
	
	print ("processing the fastas for chewBBACA usability")
	
	ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'init_schema_4_bbaca.py')
	
	args = [ScriptPath, '-i', path2down, '--cpu', str(cpu2use)]
	proc = subprocess.Popen(args)
	proc.wait()
	
	print ("\n################\n" +str(len(orderedkeys))+" loci downloaded")
	print ("Schema is now available at "+path2down)

def main():

	parser = argparse.ArgumentParser(description="Download a schema from NS")
	parser.add_argument('-i', nargs='?', type=str, help='uri to NS schema', required=True)
	parser.add_argument('-p', nargs='?', type=str, help='path where to down', required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

	args = parser.parse_args()

	uri2schema = args.i
	path2down = args.p
	cpu2use = args.cpu
	#~ uri2schema="http://137.205.69.51/app/NS/species/1/schemas/1"
	
	lala=get_schema(uri2schema,path2down,cpu2use)
	
	
	  

if __name__ == "__main__":
    main()
