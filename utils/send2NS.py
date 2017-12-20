#!/usr/bin/env python3

import requests,json
import csv
import argparse
from collections import defaultdict
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def send_post(loci_uri,sequence,token):
	params = {}
	params['sequence'] = sequence
	headers={'Authentication-Token': token}
	url = loci_uri+"/alleles"
	r = requests.post(url, data=params,headers=headers)
	allele_url=((r.content).decode("utf-8")).replace('"', '').strip()
	
	return allele_url

def read_metadata(metadataFile):
	metadata={}
	with open(metadataFile) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		
		firstrow=next(reader)
		if "FILE" not in firstrow or "ACCESSION" not in firstrow or  "COUNTRY" not in firstrow or  "ST" not in firstrow :
			print (firstrow)
			return "Not correct header"
		
		print (firstrow)
		for row in reader:
			metadata[row[0]]=row[1:]
			
	
	return 	metadata
	

def collect_new_alleles_seq(profileFile,path2Schema, token,schemaURI):
	profileDict2={}
	with open(profileFile) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			profileDict2[row[0]]=row
	
	
	profile=[]
    
    
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
			if int(len(listHeaders)/2)<=missing_data_count:
				i+=1
				print (row[0]+" discarded from further analysis because missing data count is over 50%")
				listGenomes2discard.append(row[0])
				continue

			hasAlleles=False
			while i<len(row):
				
				gene= headers[i]
				if "*" in row[i]:
					new_allele= int(row[i].replace('*',''))
					dict_new_alleles[headers[i]].append(new_allele)
                
				i+=1
                

    
	sequencesdict={}
	orderSequences=[]

	for gene,newAlleles in dict_new_alleles.items():
		newAlleles.sort()
		newAlleles=list(set(newAlleles))
		gene_full_path = os.path.join(path2Schema,gene)
		max_id_new_allele=max(newAlleles)

		#~ print gene_full_path
		for allele in SeqIO.parse(gene_full_path, "fasta", generic_dna):
			sequence=str(allele.seq)
			allele_id=int(str(allele.name).split("_")[-1].replace("*",""))			

			if allele_id in newAlleles:
				orderSequences.append(allele.name)
				sequencesdict[allele.name]=sequence
			if max_id_new_allele<=allele_id:
				break
	
	
	#get dictionary locus_id --> gene_name
	dictgenes={}
	r = requests.get(schemaURI+"/loci")
	result=r.json()
	
	#remove server time
	result.pop()
	for gene in result:
		dictgenes[str(gene['name']['value'])]=str(gene['locus']['value'])
	
	
		
	#read the fasta and process each new allele, send the sequence to the server and replace the attributed allele on the profiles
    
	for allele_name in orderSequences:
		sequence=sequencesdict[allele_name]
        
		#~ sequence=str(allele.seq)
		prev_allele_id=str(allele_name).split("_")[-1]
		allele_id=int(str(allele_name).split("_")[-1].replace("*",""))
		gene=str(allele_name).split("_")[0]
		print ("Sending alleles of: "+gene)

		loci_uri= dictgenes[gene+".fasta"]


		# post the sequence to the respective locus, api return the new allele or the allele if already on db
		#if token is not provided just query the db for the sequence
		try:
			if token == False:
				params = {}
				params['sequence'] = sequence
				r = requests.get(loci_uri+"/sequences",data=params)
				result=r.json()
				try:
					print (result[0]['id']['value'])
					new_allele_id=result[0]['id']['value']
				except:
					new_allele_id=prev_allele_id
			else:
				allele_url=send_post(loci_uri,sequence,token)
				new_allele_id=str(int(allele_url.split("/")[-1]))
		
		except Exception as e:
			# something went wrong,wait and retry
			sleep(5)
			if token == False:
				params = {}
				params['sequence'] = sequence
				r = requests.get(loci_uri+"/sequences",data=params)
				result=r.json()
				try:
					print (result[0]['id']['value'])
					new_allele_id=result[0]['id']['value']
				except:
					new_allele_id=prev_allele_id
			else:
				allele_url=send_post(loci_uri,sequence,token)
				new_allele_id=str(int(allele_url.split("/")[-1]))


		#replace the old id on the profile for the new one
		profileDict[gene+".fasta"]=[new_allele_id if x=="*"+str(allele_id) else x for x in profileDict[gene+".fasta"]]
	
	lists2print=[]
	modifiedProfileDict={}
	i=0
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

	
	if not token == False:
		for genome,profile in modifiedProfileDict.items():
			
			if genome in listGenomes2discard:
				print (genome+" profile has a low locus count, probably a different species")
				continue
		
			# schema uri is like http://137.205.69.51/app/NS/species/1/schemas/1 remove last 2 and get species uri
			server_ip = schemaURI.split("/")
			server_ip.pop()
			server_ip.pop()
			server_ip=(str.join("/",server_ip))+"/profiles"
			
			headers = {'Content-Type': 'application/json','Authentication-Token': token}
			event_data = {}
			event_data['profile']={genome:profile}
			event_data['headers']=listHeaders

			#~ server_return = requests.post(server_ip, headers=headers, json=event_data)
			server_return = requests.post(server_ip, headers=headers, data=json.dumps(event_data))
			print (server_return.content.decode("utf-8"))
	
	with open("newProfile.tsv", 'w') as f:
		f.write(newProfileStr)
	
	
	return "Done"

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

def main():

	parser = argparse.ArgumentParser(description="This synchronizes a local profile and it's new alleles with the NS")
	parser.add_argument('-s', nargs='?', type=str, help='path to schema folder', required=True)
	parser.add_argument('-p', nargs='?', type=str, help='tsv with profile', required=True)
	parser.add_argument('-m', nargs='?', type=str, help='tsv with metadata', required=False, default=False)
	parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)

	args = parser.parse_args()

	profileFile = args.p
	pathSchema = args.s
	token= args.t
	metadataFile= args.m
	
	if not metadataFile==False:
		metadataFile=read_metadata(metadataFile)
	
	if metadataFile == '"Not correct header"':
		print ("Error on metadata file header")
		return
	#~ print (metadata)
	#~ asdasd
	
	lastSyncServerDate,schemaUri=get_schema_vars(pathSchema)
	
	collect_new_alleles_seq(profileFile,pathSchema,token,schemaUri)
	
	  

if __name__ == "__main__":
    main()
