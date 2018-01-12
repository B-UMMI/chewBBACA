#!/usr/bin/env python3

import requests,json
import argparse
from collections import defaultdict
import multiprocessing

class Result():
	def __init__(self):
		self.val = []

	def update_result(self, val):
		lala=val
		self.val.append(lala)

	def get_result(self):
		return self.val

def get_isol_profile(uri,schema,orderedLoci,auxBar):


	schemaId=(schema.split("/"))[-1]
	uriProfile=uri+"/schemas/"+schemaId
	r = requests.get(uriProfile, timeout=10)
	result=r.json()

	#create a dictionary with each locus as key
	resultsDict= defaultdict(list)
	for locus in orderedLoci:
		resultsDict[locus]=[]

	#fill the dictionary with the alleles
	for locus in result:
		locusName=locus['name']['value']
		try:
			locusAllele=int(locus['id']['value'])
		except:
			locusAllele=0
		resultsDict[locusName].append(locusAllele)

	if uri in auxBar:
		auxlen=len(auxBar)
		index=auxBar.index(uri)
		print ( "["+"="*index+">"+" "*(auxlen-index)+"] Downloading profiles "+str(int((index/auxlen)*100))+"%")

	return [uri,resultsDict]

def get_profiles(species,schema,cpu2use):



	listOfIsolates=species+"/isolates"

	#build dictionary locus_id --> gene_name
	dictIsol={}
	r = requests.get(listOfIsolates)
	result=r.json()


	print ("Number of profiles to Down: "+str(len(result)))
	#~ print (result)
	#~ asdas
	for isolate in result:
		dictIsol[str(isolate['isolates']['value'])]=str(isolate['name']['value'])

	listOfILoci=schema+"/loci"

	#build dictionary locus_id --> gene_name
	dictLoci={}
	r = requests.get(listOfILoci)
	result2=r.json()
	#remove timestamp
	result2.pop()

	print ("Number of Loci on schema: "+str(len(result2)))
	for loci in result2:
		dictLoci[str(loci['name']['value'])]=str(loci['locus']['value'])

	#test down bar
	auxBar=[]
	orderedkeys=sorted(dictIsol.keys())
	step=int((len(orderedkeys))/10)+1
	counter=0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter+=step

	orderedLoci=sorted(dictLoci.keys())

	#get each profile and save results in memory
	pool = multiprocessing.Pool(cpu2use)
	result = Result()
	for uri in orderedkeys:
		p=pool.apply_async(get_isol_profile,args=[uri,schema,orderedLoci,auxBar],callback=result.update_result)

	pool.close()
	pool.join()
	listResults=result.get_result()


	#create final results list of profiles
	listFinalProfile=[["FILE"]+orderedLoci]

	for isolate in listResults:
		aux=[]
		isolUri=isolate[0]
		profile=isolate[1]
		aux.append(dictIsol[isolUri])
		for loci in orderedLoci:
			alleleId=0
			try:
				alleleId=(profile[loci])[0]
			except:
				pass
			aux.append(alleleId)

		listFinalProfile.append(aux)


	with open("profiles.tsv", "w") as f:
		for profile in listFinalProfile:
			#~ print (profile)
			f.write(('\t'.join(map(str, profile)))+"\n")


	print ("Profile is now available at profiles.tsv")

def main():

	parser = argparse.ArgumentParser(description="Download profiles")
	parser.add_argument('--sp', nargs='?', type=str, help='species uri', required=True)
	parser.add_argument('--sc', nargs='?', type=str, help='schema uri', required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

	args = parser.parse_args()

	species = args.sp
	schema = args.sc
	cpu2use = args.cpu

	lala=get_profiles(species,schema,cpu2use)




if __name__ == "__main__":
	main()
