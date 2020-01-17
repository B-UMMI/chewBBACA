#!/usr/bin/env python3

import requests,json
import csv
import multiprocessing
from collections import defaultdict
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import time

try:
    from utils import send_metadata
except ImportError:
    from CHEWBBACA_NS.utils import send_metadata

class Result():
	def __init__(self):
		self.val = []

	def update_result(self, val):
		lala = val
		self.val.append(lala)

	def get_result(self):
		return self.val

def send_post(loci_uri,sequence,token):
	params = {}
	params['sequence'] = sequence
	headers={'Authentication-Token': token}
	url = loci_uri+"/alleles"
	r = requests.post(url, data=params,headers=headers,timeout=60)
	req_code=int(r.status_code)
	allele_url=((r.content).decode("utf-8")).replace('"', '').strip()

	return allele_url,req_code,r.content


def send_sequence(token,sequence,loci_uri,prev_allele_id):

	new_allele_id = "*" + str(prev_allele_id)
	#if no token is provided, send to the server to see if it already exists
	if token == False:
		params = {}
		params['sequence'] = sequence.upper()


		sucess_send = False
		attempts=0
		waitFactor = 4
		while not sucess_send and attempts<5:

			r = requests.get(loci_uri + "/sequences", data=params, timeout=60)
			req_code = int(r.status_code)

			if req_code > 201:
				print("Server returned code " + str(req_code))
				print(loci_uri)
				time.sleep(waitFactor)
				waitFactor=waitFactor*2
				attempts+=1
				#raise
			else:
				result = r.json()
				sucess_send=True
				#if returns a value get the new allele, else the allele is not on the database, return the same id
				try:
					print (result[0]['id']['value'])
					new_allele_id=result[0]['id']['value']
				except:
					new_allele_id="*"+str(prev_allele_id)

	#token was provided, send the allele with token
	else:

		sucess_send = False
		waitFactor = 4
		attempts=0
		while not sucess_send and attempts<5:

			allele_url,reqCode,req_content=send_post(loci_uri,sequence,token)
			#if token is not valid return the new allele as the old one
			if reqCode==401:
				print ("Token is not valid")
				sucess_send=True
				new_allele_id = str(prev_allele_id) + "*"
			elif reqCode==409:
				print("HASH COLLISION!!!1 CONTACT ADMIN")
				print(req_content.decode("utf-8"))
				print("Server returned code " + str(reqCode))
				sucess_send = True

			#something went wrong, retry the request
			elif reqCode>201:
				print ("Server returned code "+str(reqCode))
				print("Retrying in seconds "+str(waitFactor))
				print(loci_uri)
				time.sleep(waitFactor)
				waitFactor = waitFactor * 2
				attempts+=1
			#everything went ok
			else:
				new_allele_id=str(int(allele_url.split("/")[-1]))
				sucess_send = True

	return new_allele_id

#def process_alleles(token, sequence, loci_uri, prev_allele_id,allele_id):

def process_locus(gene,newAlleles,path2Schema,token, loci_uri,auxBar):

	#sort list of new alleles id to start from smaller to larger
	newAlleles.sort()
	max_id_new_allele = int(max(newAlleles))
	gene_full_path = os.path.join(path2Schema, gene)

	results=[]

	# for allele in the locus fasta process it
	for allele in SeqIO.parse(gene_full_path, "fasta", generic_dna):
		sequence = str(allele.seq)

		#if * not in the allele name it's not a local allele only
		if "*" not in str(allele.name):
			continue

		#get allele id, remover all extra from string and only get the int number
		allele_id = int(str(allele.name).split("_")[-1].replace("*", ""))

		#if the allele is on the list of alleles to send from the profile, send it
		if allele_id in newAlleles:
			#orderSequences.append(allele.name)
			#sequencesdict[allele.name] = sequence
			prev_allele_id=allele_id

			sucess_send=False
			waitFactor=4
			attempts=0
			while not sucess_send and attempts<5:

				try:
					new_allele_id = send_sequence(token, sequence, loci_uri, prev_allele_id)
					results.append((prev_allele_id,new_allele_id))
					sucess_send=True

				except Exception as e:

					#something went wrong, will try again after wait facto
					time.sleep(waitFactor)
					waitFactor=waitFactor*2
					attempts+=1
					print(prev_allele_id)
					print("something went wrong running the send_sequence function, retrying in sec"+str(waitFactor))


		#if this allele id is larger than the possible max id allele to send is because there wont be more alleles that will be sent?
		#if max_id_new_allele <= allele_id:
		#	break

	if gene in auxBar:
		auxlen = len(auxBar)
		index = auxBar.index(gene)
		print("[" + "=" * index + ">" + " " * (auxlen - index) + "] Sending alleles " + str(
			int((index / auxlen) * 100)) + "%")


	return [gene,results]


def process_genome(genome,schemaURI,profile,listHeaders,token,auxBar):

	# schema uri is like http://137.205.69.51/app/NS/species/1/schemas/1 remove last 2 and get species uri
	server_ip = schemaURI.split("/")
	server_ip.pop()
	server_ip.pop()
	server_ip = (str.join("/", server_ip)) + "/profiles"

	headers = {'Content-Type': 'application/json', 'Authentication-Token': token}
	event_data = {}
	event_data['profile'] = {genome: profile}
	event_data['headers'] = listHeaders

	# ~ server_return = requests.post(server_ip, headers=headers, json=event_data)

	sucess_send = False
	waitFactor = 4
	attempts=0
	while not sucess_send and attempts<5:

		server_return = requests.post(server_ip, headers=headers, data=json.dumps(event_data))

		if server_return.status_code == 401:
			print("Token is not valid")
			sucess_send = True

		elif server_return.status_code>201:
			time.sleep(waitFactor)
			waitFactor=waitFactor*2
			attempts+=1

		else:
			sucess_send=True

	print(server_return.content.decode("utf-8"))

	if genome in auxBar:
		auxlen = len(auxBar)
		index = auxBar.index(genome)
		print("[" + "=" * index + ">" + " " * (auxlen - index) + "] Sending genomes " + str(
			int((index / auxlen) * 100)) + "%")

	try:
		time.sleep(1)
		genome_new_uri=((server_return.content.decode("utf-8").split("at "))[-1]).strip()
	except:
		genome_new_uri=False

	return [genome,genome_new_uri]

def collect_new_alleles_seq(profileFile,path2Schema, token,schemaURI,cpu2Use,percentMDallowed,metadataFile):

	#read profile to upload and put into a dictionary {genome:[full profile]}
	profileDict2={}
	with open(profileFile) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			profileDict2[row[0]]=row

	#build dictionary with {locus:[alleles per locus]}
	profileDict= defaultdict(list)
	listHeaders=profileDict2.pop("FILE")
	genomes=profileDict2.keys()
	i=0
	while i<len(listHeaders):

		for genome in genomes:
			profileDict[listHeaders[i]].append(profileDict2[genome][i])

		i+=1


	dict_new_alleles= defaultdict(list)


	listGenomes2discard=[]

	#re-read the profile and count the missing data per genome
	with open(profileFile) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		#~ headers = reader.next()
		headers = next(reader)

		for row in reader:

			missing_data_count=0
			missing_data_count+=row.count("LNF")
			missing_data_count+=row.count("NIPH")
			missing_data_count+=row.count("PLOT")
			missing_data_count+=row.count("ASM")
			missing_data_count+=row.count("ALM")

			i=1
			#~ print (int(len(listHeaders)/2))
			#~ print (missing_data_count)

			#if number of missing data is larger than the missing data allowed the genome is discarded
			if int(len(listHeaders)*(percentMDallowed))<=missing_data_count:
				i+=1
				print (row[0]+" discarded from further analysis because missing data count is over "+str(percentMDallowed))
				listGenomes2discard.append(row[0])
				continue

			hasAlleles=False

			#if number of missing data is acceptable, get the alleles that are only local and store them on a dictionary {locus:[new alleles]}
			while i<len(row):

				gene= headers[i]
				if "*" in row[i]:
					new_allele= int(row[i].replace('*',''))
					dict_new_alleles[headers[i]].append(new_allele)

				i+=1


	#get list of loci from schema and build a dictionary {locus_id: gene_name}
	dictgenes={}
	r = requests.get(schemaURI+"/loci")
	result=r.json()

	#remove server time
	#result.pop()

	#get only the locus
	result = result["Loci"]
	serverTime = r.headers['Server-Date']

	for gene in result:
		dictgenes[str(gene['name']['value'])+".fasta"]=str(gene['locus']['value'])


	#read the fasta and process each new allele, send the sequence to the server and replace the attributed allele on the profiles


	# this is for the progress bar
	auxBar = []
	orderedkeys = sorted(dict_new_alleles)
	step = int((len(orderedkeys)) / 10) + 1
	counter = 0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter += step

	# this is for the callback of the sending of alleles
	new_result = Result()

	#send the new alleles per locus
	pool = multiprocessing.Pool(cpu2Use)
	for gene in sorted(dict_new_alleles):
		newAlleles=list(set(dict_new_alleles[gene]))
		loci_uri = dictgenes[gene]
		p = pool.apply_async(process_locus, args=[gene, newAlleles, path2Schema, token, loci_uri,auxBar],
							 callback=new_result.update_result)
	pool.close()
	pool.join()

	#here we build a list of the alleles that were renumbered
	listChangedAlleles = new_result.get_result()
	number_new_alleles=0

	#listChangedAlleles is constituted of [[locus],[[old_id,new_id],[old_id,new_id]]
	for elem in listChangedAlleles:
		try:
			gene = elem[0]
			pass
		except:
			continue
		for elem2 in elem[1]:
			new_allele_id=elem2[1]
			allele_id=elem2[0]

			#if there is a * it's because the token was not provided and the allele was not found on the server, else, replace the old allele for the new
			if "*" in str(new_allele_id):
				continue
			else:
				number_new_alleles += 1
				profileDict[gene] = [new_allele_id if x == "*" + str(allele_id) else x for x in
										profileDict[gene]]

	print(str(number_new_alleles)+" new alleles were uploaded")

	#give the server a minute of rest if a lot of alleles were sent, you already waited a lot so it wont make difference anyway :))
	if number_new_alleles>100000:
		time.sleep(300)

	lists2print=[]
	modifiedProfileDict={}
	i=0

	#create the string of the new profile
	while i<len(profileDict["FILE"]):
		aux=[]
		for header in listHeaders:
			aux.append(profileDict[header][i])
		i+=1
		modifiedProfileDict[aux.pop(0)]=aux
		lists2print.append(aux)

	newProfileStr=('\t'.join(map(str,listHeaders)))+"\n"
	i=0
	for elem in lists2print:
		newProfileStr += str(list(genomes)[i])+"\t"+('\t'.join(map(str, elem)))+"\n"
		i+=1

	#if there is a token, upload the profiles
	if not token == False:

		new_result_genomes = Result()
		# test down bar
		auxBar = []
		orderedkeys = sorted(modifiedProfileDict)
		step = int((len(orderedkeys)) / 10) + 1
		counter = 0
		while counter < len(orderedkeys):
			auxBar.append(orderedkeys[counter])
			counter += step

		pool = multiprocessing.Pool(cpu2Use)
		for genome in sorted(modifiedProfileDict):

			if genome in listGenomes2discard:
				print (genome+" profile has a low locus count, probably a different species")
				continue

			else:
				profile = modifiedProfileDict[genome]
				p = pool.apply_async(process_genome, args=[genome, schemaURI, profile, listHeaders, token,auxBar],callback = new_result_genomes.update_result)

		pool.close()
		pool.join()

		listUploadedGenomes = new_result_genomes.get_result()
		print(str(len(listUploadedGenomes))+" profiles were uploaded")

		#if to send metadata process the metadata info and send it
		if metadataFile:
			dicgenomes2sendMeta={}
			for genome in listUploadedGenomes:
				if genome[1]:
					try:
						dicgenomes2sendMeta[genome[0]]=(genome[1]).replace('"', '').strip()
					except:
						#genome is not on the metadata file
						continue



			print("will send the metadata of the following isolates")
			print(dicgenomes2sendMeta)
			send_metadata.main(metadataFile, cpu2Use,token,str(dicgenomes2sendMeta))



	with open("newProfile.tsv", 'w') as f:
		f.write(newProfileStr)


	return "Done"

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
	percentMDallowed=float(percentMDallowed)


	#get schema uri and date from schema .config file
	lastSyncServerDate,schemaUri=get_schema_vars(pathSchema)

	server_url = (schemaUri.split("/species"))[0]

	r = requests.get(server_url, timeout=5)
	if r.status_code > 200:
		print("server returned error with code " + str(
			r.status_code) + ", program will stop since there seems to be a problem contacting the server")
		return
	else:
		print("server is running, will start the program")

	collect_new_alleles_seq(profileFile,pathSchema,token,schemaUri,cpu2Use,percentMDallowed,metadataFile)



if __name__ == "__main__":
	main()
