#!/usr/bin/env python3

import requests
from collections import defaultdict
import os
import pickle
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import multiprocessing
import shutil
import sys
import time

try:
	from utils import CommonFastaFunctions
except ImportError:
	from CHEWBBACA_NS.utils import CommonFastaFunctions

try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO


def translateSeq(DNASeq):
	seq=DNASeq
	reversedSeq=False
	tableid=11
	originalSeq=True
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=tableid,cds=True)
	except:
		reversedSeq=True
		originalSeq=False
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			protseq=Seq.translate(myseq, table=tableid,cds=True)

		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				protseq=Seq.translate(myseq, table=tableid,cds=True)
			except:
				reversedSeq=False
				try:
					seq=seq[::-1]
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					protseq=Seq.translate(myseq, table=tableid,cds=True)
				except Exception as e:
					print ("translation error")
					print (e)
					raise
	return str(protseq), seq, originalSeq

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	strDNArevC = ''
	for l in strDNA:

		strDNArevC += basecomplement[l]

	return strDNArevC[::-1]

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

	lastSyncServerDate=lines[0]
	schemaUri=lines[1]

	return lastSyncServerDate,schemaUri

def getBlastScoreRatios(genefile, basepath, doAll, blastPath):

	listAllelesNames = []
	# get bsr for each allele
	for allele in SeqIO.parse(genefile, "fasta", generic_dna):

		listAllelesNames.append(allele.id)

		geneScorePickle = os.path.abspath(genefile) + '_bsr.txt'
		with open(geneScorePickle, 'rb') as f:
			var = pickle.load(f)

	# returning all allele BSR scores and list of alleles for this gene
	return var, listAllelesNames

def reDogetBlastScoreRatios(sequence, alleleI, allelescores2, newGene_Blast_DB_name, picklepath,
							blastPath, listAllelesNames):


	cline = NcbiblastpCommandline(cmd=blastPath, db=newGene_Blast_DB_name, evalue=0.001, outfmt=5,num_threads=1)

	out, err = cline(stdin=sequence)

	allelescore = 0

	psiblast_xml = StringIO(out)
	blast_records = NCBIXML.parse(psiblast_xml)


	matchscore = 0
	for blast_record in blast_records:

		for alignment in blast_record.alignments:

			for match in alignment.hsps:
				matchscore = int(match.score)


	allelescores2[alleleI] = matchscore

	with open(picklepath, 'wb') as f:
		pickle.dump(allelescores2, f)

	return allelescores2, listAllelesNames


def process_locus(gene,path2schema,newAllelesDict,auxBar,blastPath,temppath,bsrTresh):

	bsrTresh = float(bsrTresh)
	maxBsrShort=0.7
	#print(gene)
	geneFile=os.path.join((os.path.abspath(path2schema)),gene)
	shortgeneFile = os.path.join(os.path.dirname(geneFile), "short", os.path.basename(geneFile))
	shortgeneFile = shortgeneFile.replace(".fasta", "_short.fasta")
	geneBaseName = os.path.splitext(os.path.basename(geneFile))[0]

	basepath = os.path.join(temppath, geneBaseName)
	perGeneBasepath = os.path.join(basepath, geneBaseName)

	shortAlleleDict = {}
	localAllelesDict = {}
	listShortAllelesNames = []
	finalFastaString = ''

	alleleShortProt = ''

	for allele in SeqIO.parse(geneFile, "fasta", generic_dna):

		auxid = ((allele.id).split("_"))[-1]
		sequence = str(allele.seq.upper())
		if "*" in auxid:
			localAllelesDict[sequence] = auxid
		else:
			finalFastaString += ">" +geneBaseName+"_"+ str(auxid) + "\n" + str(sequence + "\n")

	dict_local2NS_corr = {}

	for k,v in newAllelesDict.items():

		sequenceuri = v

		sucess_send = False
		waitFactor = 4
		attempts=0
		while not sucess_send and attempts<5:

			r = requests.get(sequenceuri,timeout=30)

			if r.status_code >201:
				print("Server returned code " + str(req_code))
				time.sleep(waitFactor)
				waitFactor = waitFactor * 2
				attempts+=1
			else:
				sucess_send=True

		result = r.json()
		dnaSequence = str(result[0]['sequence']['value'])
		newAllelesDict[k] = dnaSequence

	# check which local alleles are the same as new alleles from NS and keep info on dict
	for alleleIaux in sorted(newAllelesDict):

		newSeq = newAllelesDict[alleleIaux]

		try:
			local_index = localAllelesDict[newSeq]
			dict_local2NS_corr[local_index] = alleleIaux

		except:
			pass

	#print(dict_local2NS_corr)
	for allele in SeqIO.parse(shortgeneFile, "fasta", generic_dna):
		sequence = str(allele.seq.upper())
		auxid = ((allele.id).split("_"))[-1]

		# if there is a short allele that is a local allele to be changed by a NS, start to use the NS index id
		if auxid in dict_local2NS_corr.keys():
			auxid = dict_local2NS_corr[auxid]
		shortAlleleDict[auxid] = sequence
		listShortAllelesNames.append(auxid)
		protseq, Inverted, seq = translateSeq(sequence)
		alleleShortProt += ">" + str(auxid) + "\n" + str(protseq + "\n")

	geneScorePickle = os.path.abspath(shortgeneFile) + '_bsr.txt'

	allelescores, listShortAllelesNames = getBlastScoreRatios(shortgeneFile, basepath, False,
															  blastPath)

	alleleScoresChanged=False
	# replace all local alleles index to be changed by NS on the allelescores dict
	for k, v in dict_local2NS_corr.items():
		try:
			lala = allelescores[k]
			allelescores[v] = allelescores.pop(k)
			auxSeq = newAllelesDict[v]
			localAllelesDict[auxSeq] = v
			alleleScoresChanged=True
		except Exception as e:
			pass

	Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(perGeneBasepath, 1, True, alleleShortProt)

	#calculate the bsrs of the new alleles
	for alleleIaux in sorted(newAllelesDict):


		newSeq = newAllelesDict[alleleIaux]

		protSeq, Inverted, seq = translateSeq(newSeq)

		cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10,
									  max_hsps=10, num_threads=1)
		try:
			out, err = cline(stdin=protSeq)
		except Exception as e:
			print(e)
		psiblast_xml = StringIO(out)
		blast_records = NCBIXML.parse(psiblast_xml)
		bestmatch = [0, 0]
		try:

			# iterate through the blast results
			for blast_record in blast_records:

				for alignment in blast_record.alignments:

					# select the best match
					for match in alignment.hsps:

						# query id comes with query_id, not name of the allele
						alleleMatchid = int((blast_record.query_id.split("_"))[-1])

						# ~ scoreRatio=float(match.score)/float(allelescores[int(alleleMatchid)-1])
						# query_id starts with 1
						alleleMatchid2 = (((listShortAllelesNames[alleleMatchid - 1]).split("_"))[-1])

						try:
							scoreRatio = float(match.score) / float(allelescores[int(alleleMatchid2)])
						except:
							#print(allelescores),alleleMatchid2
							scoreRatio = float(match.score) / float(allelescores[alleleMatchid2])

						# DNAstr = str(currentCDSDict[">" + cdsStrName])
						if scoreRatio == 1 and match.score > bestmatch[0]:
							bestmatch = [match.score, scoreRatio]


						elif (match.score > bestmatch[0] and scoreRatio >= bsrTresh and scoreRatio > bestmatch[1]):
							bestmatch = [match.score, scoreRatio]

		except Exception as e:
			print("some error occurred")
			print(e)
			print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))

		#if float(bestmatch[1]) >= bsrTresh and float(bestmatch[1]) < bsrTresh + 0.1:
		if float(bestmatch[1]) < maxBsrShort:
			# newDNAAlleles2Add2shortFasta +=">" +geneBaseName + "_"+ str(alleleIaux) + "\n" + newSeq + '\n'
			shortAlleleDict[alleleIaux] = newSeq

			alleleShortProt += ">" + str(alleleIaux) + "\n" + newSeq + '\n'
			shortAlleleDict[auxid] = sequence

			geneTransalatedPath2 = os.path.join(basepath, str(
				os.path.basename(shortgeneFile) + '_protein2.fasta'))

			# proteinFastaString += '>' + alleleIaux + '\n' + str(protSeq) + '\n'

			# --- remake blast DB and recalculate the BSR for the locus --- #
			listShortAllelesNames.append(alleleIaux)
			shortAlleleDict[auxid] = sequence

			sequence_2_blast = '>' + alleleIaux + '\n' + str(protSeq)

			Gene_Blast_DB_name2 = CommonFastaFunctions.Create_Blastdb_no_fasta(geneTransalatedPath2, 1,
																			   True, sequence_2_blast)

			Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(perGeneBasepath, 1, True, alleleShortProt)

			#print(alleleIaux+" is a short ###############################################")
			alleleScoresChanged=False
			allelescores, listShortAllelesNames = reDogetBlastScoreRatios(sequence_2_blast,
																		  alleleIaux,
																		  allelescores,
																		  Gene_Blast_DB_name2,
																		  geneScorePickle,
																		  blastPath,
																		  listShortAllelesNames)

		if alleleIaux not in localAllelesDict.values():
			localAllelesDict[newSeq] = alleleIaux

	#print(finalFastaString)

	auxDict = {}
	#remove the * from the local alleles that were on the new alleles synced
	for lala in localAllelesDict.values():
		lala = str(lala)
		new_id = lala.replace("*", "")
		auxDict[lala] = int(new_id)

	newSortedIds = sorted(auxDict, key=auxDict.get)

	#build the fasta string by sorted ids
	for SortedId in newSortedIds:

		for k, v in localAllelesDict.items():
			if v == SortedId:
				finalFastaString += ">" + geneBaseName + "_" + str(v) + "\n" + str(k) + "\n"
				break

	#if new alleles bsr added, create new pickle
	if alleleScoresChanged:
		with open(geneScorePickle, 'wb') as f:
			pickle.dump(allelescores, f)

	#build the fasta string for the short file
	finalShortFastaString = ''
	for k, v in shortAlleleDict.items():
		finalShortFastaString += ">" + geneBaseName + "_" + str(k) + "\n" + str(v) + "\n"

	#print(finalFastaString)
	with open(geneFile, 'w') as fG:
		fG.write(finalFastaString)

	with open(shortgeneFile, 'w') as fG:
		fG.write(finalShortFastaString)

	try:
		shutil.rmtree(basepath)
	except Exception as e:
		print(e)
		pass


	if gene in auxBar:
		auxlen=len(auxBar)
		index=auxBar.index(gene)
		print ( "["+"="*index+">"+" "*(auxlen-index)+"] Syncing "+str(int((index/auxlen)*100))+"%")


	return geneFile


def getNewsequences(newDict,lastSyncServerDate,schemaUri,numberSequences):

	print("Start trying to get new sequences from server")
	#request the new alleles starting on the date given
	params = {}
	params['date'] = lastSyncServerDate
	uri = schemaUri + "/loci"
	r = requests.get(uri, data=params)


	result = r.json()
	result=result["newAlleles"]
	#print(len())
	#print(result["newAlleles"][-1]['date'])

	#remove the server date from the request last result, the request returns the list of alleles sorted meaning the last is the oldest from the list
	#servertime = result.pop()['date']

	#if no Last-Allele on headers is because there are no alleles
	try:
		servertime = r.headers['Last-Allele']
	except:
		servertime = r.headers['Server-Date']
		print("The schema is already up to date at: " + str(lastSyncServerDate))
		return (newDict, True, servertime, numberSequences)

	print("Getting NS new alleles since " + str(lastSyncServerDate))

	#no allele returned, no more new alleles on the db
	#if len(result) < 1:
	#	print("The schema is already up to date at: " + str(lastSyncServerDate))
	#	return (newDict,True,servertime,numberSequences)

	#for each new allele add info to the dictionary
	for newAllele in result:
		geneFasta = (newAllele['locus_name']['value'])+".fasta"
		alleleid = newAllele['allele_id']['value']
		sequenceUri = newAllele['sequence']['value']
		# newDict[geneFasta].append([alleleid,sequenceUri])
		newDict[geneFasta][alleleid] = sequenceUri


	numberSequences +=len(result)
	print(r.headers)
	#if the result is the maximum number of new alleles the server can return, re-do the function, else stop doing the function
	if r.headers['All-Alleles-Returned']=='False':
		return (newDict, False, servertime, numberSequences)
	else:
		return (newDict, True, servertime, numberSequences)

	#if len(result) >= 50000:
	#	return (newDict,False,servertime,numberSequences)


def syncSchema(lastSyncServerDate,schemaUri,path2schema,cpu2use,blastPath,bsrTresh):

	listFastas2Init=[]
	#get list of sequences that are new considering the last date
	allNewSequences=False
	newDict= defaultdict(dict)
	servertime=lastSyncServerDate
	numberSequences=0

	#get all sequences until the number of new sequences is less than 100k, the maximum the server return is 100k
	while not allNewSequences:
		newDict,allNewSequences,servertime,numberSequences = getNewsequences(newDict,servertime,schemaUri,numberSequences)

	print(str(numberSequences) + " new alleles found to be synced locally")

	#the last date is the last allele date from the server
	if numberSequences < 1:
		print("The schema is now up to date at: " + str(lastSyncServerDate))
		return (lastSyncServerDate)


	#for status bar
	auxBar=[]
	orderedkeys=sorted(newDict.keys())
	step=(int((len(orderedkeys))/10))+1
	counter=0
	while counter < len(orderedkeys):
		auxBar.append(orderedkeys[counter])
		counter+=step


	print ("Adding new alleles to local files")
	#process each gene that has new alleles on ns
	pool = multiprocessing.Pool(cpu2use)
	temppath=os.path.join((os.path.abspath(path2schema)), "temp")
	for gene in orderedkeys:
		#~ lala=process_locus(gene,path2schema,newDict,auxBar)
		newDictGene=newDict[gene]
		#temppath = os.path.join((os.path.abspath(path2schema)), gene, "temp")
		p=pool.apply_async(process_locus,args=[gene,path2schema,newDictGene,auxBar,blastPath,temppath,bsrTresh])



	pool.close()
	pool.join()

	try:
		shutil.rmtree(temppath)
	except Exception as e:
		print(e)
		pass

	#listFastas2Init=result.get_result()
	print("processing the fastas for chewBBACA usability")
	# init_schema_4_bbaca.main(os.path.abspath("temp.txt"), cpu2use)
	# os.remove(os.path.abspath("temp.txt"))
	print("\n################\n" + str(len(orderedkeys)) + " loci had new alleles")
	print("Schema is now synched")

	return servertime


def main(path2schema, cpu2use, bsrTresh):

	blastPath = 'blastp'
	lastSyncServerDate, schemaUri=get_schema_vars(path2schema)

	server_url = (schemaUri.split("/species"))[0]

	r = requests.get(server_url, timeout=5)
	if r.status_code > 200:
		print("server returned error with code " + str(
			r.status_code) + ", program will stop since there seems to be a problem contacting the server")
		return
	else:
		print("server is running, will start the program")


	newserverTime=syncSchema(lastSyncServerDate, schemaUri,path2schema,cpu2use,blastPath,bsrTresh)

	if newserverTime==lastSyncServerDate:
		return

	schemapath_config=os.path.join(path2schema, ".config.txt")
	try:
		with open(schemapath_config, "w") as fp:
			fp.write(newserverTime+"\n"+schemaUri)
	except:
		print ("your schema has no config file")
		raise


if __name__ == "__main__":
	main()
