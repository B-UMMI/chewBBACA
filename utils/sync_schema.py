#!/usr/bin/env python3

import requests,json
import csv
import argparse
from collections import defaultdict
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import subprocess
import multiprocessing

def get_schema_vars(schemapath):
	
	schemapath_config=os.path.join(schemapath,".config.txt")
	try:
		with open(schemapath_config) as fp:
			lines=[]
			for line in fp:
				lines.append(line.strip())
	
	except:
		print ("your schema has no config file")
		raise
		return False 
	
	lastSyncServerDate=lines[0]
	schemaUri=lines[1]
	
	return lastSyncServerDate,schemaUri
	
def process_locus(gene,path2schema,newDict,auxBar):
	
	
	fastapath=os.path.join((os.path.abspath(path2schema)),gene)
	localnewAlleles={}
	localAlleles={}
	auxDict={}
	newAlleles=False
	shortGeneName=gene.replace(".fasta","")
	
	#read local fasta and put in memory necessary info
	for allele in SeqIO.parse(fastapath, "fasta",generic_dna):
		alleleI=((allele.name).split("_"))[-1]
		if alleleI.startswith("*"):
			localnewAlleles[str(allele.seq)]=alleleI				
		localAlleles[alleleI]=str(allele.seq)
		auxDict[alleleI]=int(alleleI.replace("*",""))
	
	#for each gene that has new alleles check if the new NS alleles are present on the local alleles
	for allele in newDict[gene]:
		alleleid=allele[0]
		sequenceuri=allele[1]
		r = requests.get(sequenceuri)
		result=r.json()
		dnaSequence=str(result[0]['sequence']['value'])
		
		# check if sequence is in one of the * local alleles
		try:
			localindex=localnewAlleles[dnaSequence]
			localAlleles[dnaSequence]=alleleid
			auxDict[alleleid]=int(alleleid.replace("*",""))
			auxDict.pop(localindex,None)
			localAlleles[alleleid]=dnaSequence
			newAlleles=True
		
		#sequence not in one of the local alleles	
		except Exception as e:
			
			#if id not attributed on local fasta is new local allele
			if not alleleid in localAlleles.keys():
				localAlleles[alleleid]=str(dnaSequence)
				auxDict[alleleid]=int(alleleid.replace("*",""))
				localAlleles[alleleid]=dnaSequence
				newAlleles=True
						

	newSortedIds=sorted(auxDict, key=auxDict.get)
	towrite=''
	
	if newAlleles:
		#~ print (fastapath)
		for newalleleid in newSortedIds:
			towrite+=">"+shortGeneName+"_"+newalleleid+"\n"+localAlleles[newalleleid]+"\n"
		with open(fastapath, 'w') as fp:
			fp.write(towrite)
	towrite=''
	

	if gene in auxBar:
		auxlen=len(auxBar)
		index=auxBar.index(gene)
		print ( "["+"="*index+">"+" "*(auxlen-index)+"] Syncing "+str(int((index/auxlen)*100))+"%")
	
	
	if newAlleles:
		return fastapath

class Result():
	def __init__(self):
		self.val = []

	def update_result(self, val):
		lala=val
		self.val.append(lala)
        
	def get_result(self):
		return self.val
	

def syncSchema(lastSyncServerDate,schemaUri,path2schema,cpu2use):
	
	listFastas2Init=[]
	
	#get list of sequences that are new considering the last date
	params = {}
	params['date'] = lastSyncServerDate
	uri=schemaUri+"/loci"
	r = requests.get(uri, data=params)
	result=r.json()
	#get last element that contains the server time
	servertime=result.pop()
	
	print ("Get NS new alleles since "+str(lastSyncServerDate))

	newDict= defaultdict(list)
	for newAllele in result:
		geneFasta=newAllele['locus_name']['value']
		alleleid=newAllele['allele_id']['value']
		sequenceUri=newAllele['sequence']['value']
		newDict[geneFasta].append([alleleid,sequenceUri])
	
	
	#test down bar
	auxBar=[]
	orderedkeys=sorted(newDict.keys())
	step=int((len(orderedkeys))/10)
	counter=0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter+=step
	
	result = Result()
	
	print ("Adding new alleles to local files")
	#process each gene that has new alleles on ns
	pool = multiprocessing.Pool(cpu2use)
	for gene in orderedkeys:
		#~ lala=process_locus(gene,path2schema,newDict,auxBar)
		p=pool.apply_async(process_locus,args=[gene,path2schema,newDict,auxBar],callback=result.update_result)
	
	pool.close()
	pool.join()
	listFastas2Init=result.get_result()
	
	towrite=''
	syncedGenes=0
	for elem in listFastas2Init:
		if not not elem:
			syncedGenes+=1
			towrite+=elem+"\n"
	with open(os.path.abspath("temp.txt"), 'w') as fp:
		fp.write(towrite)
		

		
	print ("processing the fastas for chewBBACA usability")
	ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'init_schema_4_bbaca.py')
	args = [ScriptPath, '-i', os.path.abspath("temp.txt"), '--cpu', str(cpu2use)]
	proc = subprocess.Popen(args)
	proc.wait()
	os.remove(os.path.abspath("temp.txt"))
	print ("\n################\n" +str(syncedGenes)+" loci had new alleles")
	print ("Schema is now synched")
	
	return servertime				




def main():

	parser = argparse.ArgumentParser(description="Sync a schema from NS")
	parser.add_argument('-p', nargs='?', type=str, help='path of schema', required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

	args = parser.parse_args()

	path2schema = args.p
	cpu2use = args.cpu
	
	
	lastSyncServerDate,schemaUri=get_schema_vars(path2schema)
	
	newserverTime=syncSchema(lastSyncServerDate,schemaUri,path2schema,cpu2use)
	
	schemapath_config=os.path.join(path2schema,".config.txt")
	try:
		with open(schemapath_config,"w") as fp:
			fp.write(newserverTime+"\n"+schemaUri)
	except:
		print ("your schema has no config file")
		raise
	  

if __name__ == "__main__":
    main()
