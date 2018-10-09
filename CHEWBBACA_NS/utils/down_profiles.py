#!/usr/bin/env python3

import requests
from collections import defaultdict
import multiprocessing
import csv

class Result():
	def __init__(self):
		self.val = []

	def update_result(self, val):
		lala=val
		self.val.append(lala)

	def get_result(self):
		return self.val

def process_input_profile(input_profile):

	#check if the header loci match the schema loci
	listgenomes=[]
	with open(input_profile) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')

		loci_list = next(reader)
		loci_list.pop(0)

		for row in reader:
			listgenomes.append(row[0])

		loci_list=sorted(loci_list)

	return loci_list,listgenomes

def get_isol_profile(uri,schema,orderedLoci,auxBar):

	#build the uri and get the profile
	schemaId=(schema.split("/"))[-1]
	uriProfile=uri+"/schemas/"+schemaId
	r = requests.get(uriProfile, timeout=30)
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

def get_profiles(species,schema,cpu2use,genomes2Down,inputProfile):

	dictIsol = {}
	listOfILoci=species+"/schemas/"+schema+"/loci"

	#build dictionary {locus_id:gene_name}
	dictLoci={}
	#get list of loci from schema the user provided
	r = requests.get(listOfILoci)
	result2=r.json()
	#remove timestamp
	result2.pop()

	print ("Number of Loci on schema: "+str(len(result2)))
	for loci in result2:
		dictLoci[str(loci['name']['value'])]=str(loci['locus']['value'])



	orderedLoci=sorted(dictLoci.keys())
	genomes_not2Down=[]
	equal_headers=False

	if not inputProfile is None:

		profile_list_loci,genomes_not2Down=process_input_profile(inputProfile)

		equal_headers=set(profile_list_loci)==set(orderedLoci)
		if not equal_headers:
			print("profile given by user doesn't have the same loci, will be ignored")
			genomes_not2Down = []

	#build the uri to get isolates list
	schema=species+"/schemas/"+schema
	listOfIsolates=species+"/isolates"

	#get list of isolates to download
	r = requests.get(listOfIsolates)
	result=r.json()


	print ("Number of profiles existing : "+str(len(result['Isolates'])))

	#parse
	if len(genomes2Down)>0:
		for isolate in result['Isolates']:
			isolateName=str(isolate['name']['value'])
			if isolateName in genomes_not2Down:
				pass
			elif isolateName in genomes2Down:
				dictIsol[str(isolate['isolate']['value'])] = str(isolate['name']['value'])
	else:
		for isolate in result['Isolates']:
			isolateName = str(isolate['name']['value'])
			if isolateName in genomes_not2Down:
				print(isolateName+ " is already on the profile input, will not be downloaded again")
				pass
			else:
				dictIsol[str(isolate['isolate']['value'])]=str(isolate['name']['value'])

	print("Number of profiles to Download : " + str(len(dictIsol.keys())))

	#test download status bar
	auxBar=[]
	orderedkeys=sorted(dictIsol.keys())
	step=int((len(orderedkeys))/10)+1
	counter=0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter+=step

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


	if len(dictIsol.keys()) == 0:
		print("no profiles were downloaded")

	elif equal_headers:
		listFinalProfile.pop(0)
		with open("profiles.tsv", "w") as f:
			with open(inputProfile) as f1:
				for line in f1:
					f.write(line)
			for profile in listFinalProfile:
				# ~ print (profile)
				f.write(('\t'.join(map(str, profile))) + "\n")
		print("Profile is now available at profiles.tsv")

	elif len(dictIsol.keys()) >0:
		with open("profiles.tsv", "w") as f:
			for profile in listFinalProfile:
				#~ print (profile)
				f.write(('\t'.join(map(str, profile)))+"\n")
		print("Profile is now available at profiles.tsv")



def main(species,schema,cpu2use,inputProfile,listGenomes=None):

	server_url=(species.split("/species"))[0]

	r = requests.get(server_url, timeout=5)
	if r.status_code >200:
		print("server returned error with code "+str(r.status_code)+", program will stop since there seems to be a problem contacting the server")
		return
	else:
		print("server is running, will start the program")

	#build list of genomes to download if provided, if not it downloads all profiles which take a while
	try:
		with open(listGenomes) as f:
			content = f.readlines()
		genomes2Down = [x.strip() for x in content]

	except:
		genomes2Down=[]
		pass


	lala=get_profiles(species,schema,cpu2use,genomes2Down,inputProfile)




if __name__ == "__main__":
	main()
