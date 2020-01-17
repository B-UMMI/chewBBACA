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

def get_isol_metadata(isolateURI):

	print(isolateURI)
	r = requests.get(isolateURI, timeout=30)
	result=r.json()
	result=result[0]

	dict_results={'FILE':isolateURI}
	for key in result.keys():
		dict_results[key]=result[key]['value']

	return dict_results

def get_isol_profile(uri,schema,orderedLoci,auxBar):

	#build the uri and get the profile
	schemaId=(schema.split("/"))[-1]
	uriProfile=uri+"/schemas/"+schemaId
	#print(uriProfile)
	r = requests.get(uriProfile, timeout=30)
	result=r.json()

	#create a dictionary with each locus as key
	resultsDict= defaultdict(list)
	for locus in orderedLoci:
		resultsDict[locus]=[0]

	#fill the dictionary with the alleles
	for locus in result:
		locusName=locus['name']['value']
		try:
			locusAllele=int(locus['id']['value'])
			resultsDict[locusName].pop()
			resultsDict[locusName].append(locusAllele)
		except:
			locusAllele=0
			resultsDict[locusName].pop()
			resultsDict[locusName].append(locusAllele)

	if uri in auxBar:
		auxlen=len(auxBar)
		index=auxBar.index(uri)
		print ( "["+"="*index+">"+" "*(auxlen-index)+"] Downloading profiles "+str(int((index/auxlen)*100))+"%")

	return [uri,resultsDict]

def get_profiles(species,schema,cpu2use,genomes2Down,inputProfile,token=None):


	dictIsol = {}
	listOfILoci=species+"/schemas/"+schema+"/loci"

	print("searching for the schema loci...")
	#build dictionary {locus_id:gene_name}
	dictLoci={}
	#get list of loci from schema the user provided
	r = requests.get(listOfILoci)
	result2=r.json()
	result2 = result2["Loci"]
	#serverTime = result2.pop()['date']
	serverTime = r.headers['Server-Date']

	print ("Number of Loci on schema: "+str(len(result2)))
	for loci in result2:
		dictLoci[str(loci['name']['value'])]=str(loci['locus']['value'])



	orderedLoci=sorted(dictLoci.keys())
	genomes_not2Down=[]
	equal_headers=False

	print(inputProfile)

	#if inputProfile:
	if not inputProfile is None:

		print("profile given")

		profile_list_loci,genomes_not2Down=process_input_profile(inputProfile)

		equal_headers=set(profile_list_loci)==set(orderedLoci)
		if not equal_headers:
			print("profile given by user doesn't have the same loci, will be ignored")
			genomes_not2Down = []

	#build the uri to get isolates list
	schema=species+"/schemas/"+schema
	if not token:
		listOfIsolates=species+"/isolates"
	else:
		listOfIsolates = species + "/user/isolates"

	continueDown=True

	while continueDown:

		print("searching for isolates...")
		headersToken={}
		if token:
			headersToken = {'Content-Type': 'application/json', 'Authentication-Token': token}
		params = {}

		#get list of all isolates on server if no list of genomes name provided
		if len(genomes2Down) < 1:

			if len(params.keys()) < 1:
				r = requests.get(listOfIsolates,headers=headersToken)
			else:
				r = requests.get(listOfIsolates,data=params,headers=headersToken)
			print(r.json())
			result=r.json()

			#latestIsolDate=result.pop()
			print(r.headers['Last-Isolate'])
			latestIsolDate =r.headers['Last-Isolate']
			if r.headers['All-Isolates-Returned']=='True':
				continueDown=False
			else:
				params['start'] = latestIsolDate


			print ("Number of profiles existing : "+str(len(result['Isolates'])))

		#if a list of genomes provided, check if genomes on result given by the server
		if len(genomes2Down)>0:

			print("list of genomes was provided")

			for isolateName in genomes2Down:
				params = {}
				params['isolName'] = isolateName
				r2 = requests.get(listOfIsolates,data=params)
				result2 = r2.json()
				#latestIsolDate = result2.pop()
				latestIsolDate = r.headers['Last-Isolate']

				for isolate in result2['Isolates']:
					dictIsol[str(isolate['isolate']['value'])] = str(isolate['name']['value'])
		#add all isolates returned from the server to the list of isolates to get profile
		else:
			for isolate in result['Isolates']:
				isolateName = str(isolate['isolate']['value'])
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

	print("fetching isolates profiles individually...\n")
	if len(orderedkeys)>0:
		for uri in orderedkeys:
			p=pool.apply_async(get_isol_profile,args=[uri,schema,orderedLoci,auxBar],callback=result.update_result)

		pool.close()
		pool.join()
		listResults=result.get_result()

		#get metadata
		print("\nfetching metadata...\n")
		pool = multiprocessing.Pool(cpu2use)
		result2 = Result()

		for uri in orderedkeys:
			p = pool.apply_async(get_isol_metadata, args=[uri], callback=result2.update_result)

		pool.close()
		pool.join()
		listResults2= result2.get_result()
	else:
		print("no genomes to download or fetch metadata")
		return

	#get all possible metadata headers
	allheaders=[]
	for isolate in listResults2:
		for headers in isolate.keys():
			allheaders.append(headers)

	allheaders=sorted(set(allheaders))

	print("\nwriting metadata\n")
	if len(allheaders)>0:
		with open("metadata.tsv", "w") as f:
			allheaders.remove('FILE')
			f.write('FILE\t')
			f.write(('\t'.join(map(str, allheaders))) + "\n")

			for isolate in listResults2:
				towrite=[isolate['FILE']]
				for header in allheaders:
					try:
						towrite.append(isolate[header])
					except:
						towrite.append('NA')
				f.write(('\t'.join(map(str, towrite))) + "\n")

		print("Metadata is now available at metadata.tsv")

	#create final results list of profiles
	listFinalProfile=[["FILE"]+orderedLoci]

	for isolate in listResults:
		aux=[]
		isolUri=isolate[0]
		profile=isolate[1]
		#aux.append(dictIsol[isolUri])
		aux.append(isolUri)
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
		toWrite=''

		with open(inputProfile) as f1:
			for line in f1:
				toWrite+=line
		with open("profiles.tsv", "w") as f:
				#f.write(line)
			f.write(toWrite)
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

	#with open("keys.tsv","w") as f:
	#	f.write("URI\tName"+"\n")
	#	for uri in dictIsol.keys():
	#		f.write(uri+"\t"+dictIsol[uri]+"\n")

def main(species,schema,cpu2use,inputProfile,listGenomes=None,token=None):

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

	print(inputProfile)
	print("laa")
	lala=get_profiles(species,schema,cpu2use,genomes2Down,inputProfile,token)




if __name__ == "__main__":
	main()
