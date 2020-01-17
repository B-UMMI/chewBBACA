#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import pickle
import multiprocessing
import shutil
import sys

try:
	from utils import CommonFastaFunctions
except ImportError:
	from CHEWBBACA_NS.utils import CommonFastaFunctions

#~ from cStringIO import StringIO
try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO

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
	return protseq, seq, originalSeq

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	strDNArevC = ''
	for l in strDNA:

		strDNArevC += basecomplement[l]

	return strDNArevC[::-1]


#def main(geneFile,blastPath,newAllelesDict):
def main():

	geneFile='/home/msilva/Desktop/test_init_schema/schema/ACIBA03208.fasta'
	blastPath='blastp'
	newAllelesDict={'12':'ATGCAAAACATACTTACTTTTATTAAAAATGGTCAAAAATTTACTTTAGATGCAAGCAATATTAAAAGTGTAAGAGATACGGGAAATGGTTTTGAAATCACTCTTAAAAGTGGTGAGATAATTCAAGCTGATTCCTATAGTGTTACTCCTAATACTCAGCCAGTGGAGCAAGTTTCTAGTCAGCAAGTAGAAAATAATTCTGAAGCTGATGAGGAATTTGAAAAAGAAGAAGTTTTAAAAGAAAAGAGTTTATTGGATTCACAAACAATGTTATGGGGTGGTGCTGCGATCGCACTTGGTGGTATTGCGATAGCAGCATCAGGTAGTGGTGGTAGTGATGGTTCAGCACCTTCTGATCAAACTCCTCCAGCAAGTTTAACTAAAATACTGTCTAAAGATGGTAAAGCTGTAAGTGGTTTAACAGAAGCTGGAGCAACTGTTATAGTTGAAAATTCGGCTGGAAAAGTAATCGGTTCTGCCAAAGCAGGAGCAGATGGGACATATTTAATTAATTTAGATAAAGCCTATATTAATGGTGAAATACTAAAAGTTTCTGCTCAAGATACTGCCGGTAATAGTACAGTTAAACAAGAGCTTATCGCTTCTGATATTACAGCTCCTACTTTAACACATGTTATTTCTTCTAACGGTAAAACTATTACAGGTTTGACAGAGGCTAATTCAACTGTCACCGTTAAAGATTCATCTGGAAAAATCATTGGAACTGCAAAGTCTGATAATGATGGAAAATATACAGTCATACTGGATAAAGCTTATTTAAATGGTGAAAATCTTACTATTTCTGCCGAAGATTTAGCTGGTAATAAATCAACAATTCAAACAATTTTGGCTGACGATAAAACTGCGCCAATAGGTTTAACAGTTGCAATTGATACAGCGGGTAAAGTAGTTACAGGTGTAACTGAAGCCAATGCTGTTGTGACAGTTAAGAATGCTGCAGGAATTGTAGTTGGTACAGCTACGGCTGATACCGCAGGTAAATACTCGATTACACTAGATAAAGTATATTTAAATGGTGAAAGCTTAAGCGCTACTGCATCAGATCAATCTGGTAATGCTACAGCACCAAAAACCATTATCGCACCAGATACCACGGCACCATCAAGTTTAACTGCAAGTATTGGTACTGCGGGTAAAGTAGTTACAGGTGTAACTGAAGCCAATGCTGTTGTGACAGTTAAGAATGCTGCAGGAATTGTAGTTGGTACAGCTACGGCTGATACCGTAGGTAAATACTCGATTACACTAGATAAAGTATATTTAAATGGTGAAAGCTTAAATGTAACGGCTGCGGATAAAGCGGGTAATGCTACCGTACCAAAAACAATTGTTGCACCGGATACCACAGCACCATCAAGTCTGACTGCAACGATTGATGCCGCAGGTAAAGCAATTACAGGTGTGACTGAGGCCAATGCAACAGTAACAGTTAAGAACGTGGTAGGAACCATAGTGGGTACAGGTACAGCTGATGCCACAGGTAAATATTCAGTTACCCTAGATAAAATTTATCTAAATGGCGAAAGCTTAAGTACAATAGCGGCAGATAAAGCGGGTAATGCTACCGTACCAAAAATAATAGTTGCACCAGACATTACAGCACCATCAAGTCTGACTGCAACGATTGATGCTGCAGGTAAAGCAATTACAGGTGTGACTGAGGCCAATGCAACAGTAACAGTTAAGGACGTGGTAGGAACCATAGTGGGTACAGGTACAGCCGATGTCACAGGTAAATATTCAGTTACCCTAGATAAAGCTTATCTAAATGGTGAAAGCTTAAATGCTATTGCAGCAGATAAAGTAGGTAATACTACAACACCAAAAACAATTATTGCACCGGACACTACAGCGCCTTCAAGTTTAACTGCAATTATAGATGCTGCAGGTAAAGTGATTACAGGTACGACTGAAGTAGGTGCAAGAGTGACAGTAGAACAAGTAACTGCTGTTTATAAAGAAGTAACTGTTTTAGAAACCCAAACTGTTTTGAATGAGACAGTACAGTCAAATTATTTATCGAAAACATATTCTTTCGAAGTTACAGGTACAAATGCACATGTTTCTTTAAATTTGAGCTCTTCGACTAATTCACTTTCTGGTAGTTATAGCTCAACTCTTTCTGGTGCGAATTTAAATACAAGTCTAACTGGTAGTATCTCTCAAGCTGGAAACGGAAATTATAGTATTGATCTAGCACAAGGCTCGGTATTACCGCCTGGAACCTATACTTTGACAGTAAATTATTCTAGTTCAATAGTCCCCGTTATAAATGTAAATGTAACACAAGAAGTGCCAAAAACTATTTTAGAAGTAGACCATTACGAAACAAAAGTCGTAGGTGCAGCTAATGCAGATGAGGCAGGAAATTATTCCATAACTCTTGATAAAGCTTACTTAAATGGCGAAAGCTTAAGTGTAACGGCAGCAGATCAATCGGATAATAAAACTGAAGTGAAAGAGGTTATTGCACCAGATAGTACGCCACCAATTTTGCATCAACCAACTATTCAAGGTGGATGGACTGAAGGGCAATCAGTACAAGGTACAACTGAAGCAAATGCCACAATAATTGTTAAAAATACTGCAGGTGATGTAATCGGTTCGGCAATAGCTGATGCATCCGGGTATTACAATGTAATTTTGAATACGGTTTATGAAGATGGTGAACTATTAAAAGTAATTGCTGCGGATGCTAAAGGTAATGAAAGTTCAATTAATATAAATACACCTGATATTACAGCACCTATATTAGCTAATTTGTTTAATTATGATGTTAGTACAGATAAAATTATTTTTAATGCGCCAAGTGATAGTTATATTGTTGAGCAAAAAATTGGTGATGCTTGGGTTCAGGTAAATGTTGAGGAAAAATTCGATTGGTTGAATACAGAATTCAGAGTAACTGCTAAAGATCTAGCTGGGAATAGCTCTCAACCACTAACCATTATAATTAATACTGCATCGGGAACTTATAAGCCAACTGATCCTACATTTATTCAAATTATTAAAGGATCTATAGGAAATGATTATCTTTATGGTGGTAATGGGGATGATACATTAGTTTCTAATACAGGCTCTGATTATTTGTATGGTGGCTCAGGTAATGATACGTTGATTTATGGTGGTAATTCCAATGTGTATACGGCTTTACAAGGTCAAGCTGGGAATGATACTTATATTGTAGATAAAGCATTACTTACTTCATCAAGCTCTATTCATATTTTAGATAATGCAGCTGAAGAAAATATCCTTCAATTAAAATCAGTGTCCTCTGGTGATATTTCACTTAAGCAGTCGGATTCACTCATCATAATATCATTCAATGATTCAGCATCGACAATACGTTTTGGTGAGGGGCAATTATCTTCTATTGTATTTGATGATGGAACAGTCTGGAATAAAGCACAAATTGAAGCTAATACTATTGGTAAATTATTAGGTACTGATGCAGCAGATAATCTTCAAGCAGATGCTGAGATTTCAACTATTTATGGTTTGGGAGGCAATGATACGATTCAAGGTGGTGTACAGAATGATTATCTGTATGGTGATGATGGAGACGATACATTAGTTTCTAATACAGGCTCTGATTATTTGTATGGTGGTGCAGGTAATGACACATTGATTTATGGCGGTAATTCCAATGTGTATACGGCTTTACAAGGTCAAGCTGGGAATGACACTTATATTATAGATAAAGCATTACTTACTTCATCAAGCTATATCCATATTTTAGATAATGCAACTGAAGAAAATACCCTTCAATTAAAATCAGTATCCTCTAGTGATATTTCACTTAAGCAATCCGATTCATTAATCATAATATCTTTCAATGATTCAGCATCAACAATACGTTTTGGTAAGGATAATTTATCTTTTATTGTATTTGATGATGGAACTGTTTGGGATAAAGCTCAAATTGAAGCTAATACTATTGGTAAATTATTAGGTACTGATGCAGCAGATAATCTTCAAGCAGATGCTGAGATTTCAACTATTTATAGTTTGGGAGGCAATGATACGATTCAAGGTGGTGTACAGAATGATTATCTGTATGGTGGTGATGGAGATGATACATTAGTTTCTAATACAGGCTCTGATTATTTGTATGGTGGCTCAGGTAATGATACGTTGATTTATGGTGGTAATTCCAATGTGTATACGGCTTTACAAGGTCAAGCTGGGAATGATACTTATATTGTAGATAAAGCATTACTTACTTCATCAAGCTCTATTCATATTTTAGATAATGCAGCTGAAGAAAATATCCTTCAATTAAAATCAGTGTCCTCTGGTGATATTTCACTTAAGCAGTCGGATTCACTCATCATAATATCATTCAATGATTCAGCATCGACAATACGTTTTGGTGAGGGACAATTATCTTCTATTGTATTTGATGATGGAACTGTTTGGGATAAAGCTCAAATTGAACAACATATTGCCGAACCTGTATTTGGCACAACTGGTAACGATGTAATCGAGACAAACATACCAAATCAGAGTTACTCTTATATTCTAGGAGATGGTGCAGATACTGTAGTTTTTAATATTCTCGATAATACAGATAATTTAGGTGGAAATGTTAAGACTGAATGGACAGATTTTAATCTTGTAGAGAATGATAAAATTGATTTTTCTCAACTGCTTATAAATAACAATGGAAACCTTCAAGAATTTATAACGGTAAAAGATACGCAGGCTGGAGTAGTTATGTCTGTAGATAGAGATGGCTCTAATCAATCCACTTATCATTCACAAGAATTAATTCTACTTACTGGTAAGCACTATACGTTAGAGGATTTAATGGCTTCAAATGCATTTATTCATTAA'}
	bsrTresh=0.6

	temppath = os.path.join(os.path.dirname(geneFile), "temp")
	shortgeneFile = os.path.join(os.path.dirname(geneFile), "short", os.path.basename(geneFile))
	shortgeneFile = shortgeneFile.replace(".fasta", "_short.fasta")
	geneBaseName=os.path.splitext(os.path.basename(geneFile))[0]


	basepath = os.path.join(temppath, geneBaseName)
	perGeneBasepath =os.path.join(basepath, geneBaseName)


	shortAlleleDict={}
	localAllelesDict={}
	listShortAllelesNames=[]
	finalFastaString=''

	alleleShortProt=''

	for allele in SeqIO.parse(geneFile, "fasta", generic_dna):

		auxid=((allele.id).split("_"))[-1]
		sequence = str(allele.seq.upper())
		if "*" in auxid:
			localAllelesDict[sequence]=auxid
		else:
			finalFastaString+=">" + str(auxid) + "\n" + str(sequence + "\n")

	dict_local2NS_corr={}
	#check which local alleles are the same as new alleles from NS and keep info on dict
	for alleleIaux in sorted(newAllelesDict):

		newSeq=newAllelesDict[alleleIaux]

		try:
			local_index=localAllelesDict[newSeq]
			dict_local2NS_corr[local_index]=alleleIaux

		except:
			pass


	for allele in SeqIO.parse(shortgeneFile, "fasta", generic_dna):
		sequence = str(allele.seq.upper())
		auxid=((allele.id).split("_"))[-1]

		#if there is a short allele that is a local allele to be changed by a NS, start to use the NS index id
		if auxid in dict_local2NS_corr.keys():
			print("lala")
			auxid=dict_local2NS_corr[auxid]
		shortAlleleDict[auxid]=sequence
		listShortAllelesNames.append(auxid)
		protseq, Inverted, seq = translateSeq(sequence)
		alleleShortProt+=">" + str(auxid) + "\n" + str(protseq + "\n")


	geneScorePickle = os.path.abspath(shortgeneFile) + '_bsr.txt'


	allelescores, listShortAllelesNames = getBlastScoreRatios(shortgeneFile, basepath, False,
																		  blastPath)

	# replace all local alleles index to be changed by NS on the allelescores dict
	print(allelescores)
	for k,v in dict_local2NS_corr.items():
		try:
			lala=allelescores[k]
			allelescores[v] = allelescores.pop(k)
			auxSeq=newAllelesDict[v]
			localAllelesDict[auxSeq]= v
		except Exception as e:
			pass


	Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(perGeneBasepath, 1, True, alleleShortProt)


	for alleleIaux in sorted(newAllelesDict):

		print(alleleIaux)

		newSeq=newAllelesDict[alleleIaux]

		protSeq, Inverted, seq = translateSeq(newSeq)
		protSeq=str(protSeq)

		cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name, evalue=0.001, outfmt=5, max_target_seqs=10,
									  max_hsps=10, num_threads=1)

		out, err = cline(stdin=protSeq)

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
							scoreRatio = float(match.score) / float(allelescores[alleleMatchid2])

						# DNAstr = str(currentCDSDict[">" + cdsStrName])
						print(bestmatch)
						if scoreRatio == 1 and match.score > bestmatch[0]:
							bestmatch = [match.score, scoreRatio]


						elif (match.score > bestmatch[0] and scoreRatio >= bsrTresh and scoreRatio > bestmatch[1]):
							bestmatch = [match.score, scoreRatio]

		except Exception as e:
			print("some error occurred")
			print(e)
			print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))


		if float(bestmatch[1]) >= bsrTresh and float(bestmatch[1]) < bsrTresh + 0.1:
			#newDNAAlleles2Add2shortFasta +=">" +geneBaseName + "_"+ str(alleleIaux) + "\n" + newSeq + '\n'
			shortAlleleDict[alleleIaux]=newSeq

			alleleShortProt+=">" + str(alleleIaux) + "\n" + newSeq + '\n'
			shortAlleleDict[auxid] = sequence

			geneTransalatedPath2 = os.path.join(basepath, str(
				os.path.basename(shortgeneFile) + '_protein2.fasta'))

			#proteinFastaString += '>' + alleleIaux + '\n' + str(protSeq) + '\n'

			# --- remake blast DB and recalculate the BSR for the locus --- #
			listShortAllelesNames.append(alleleIaux)
			shortAlleleDict[auxid] = sequence

			sequence_2_blast = '>' + alleleIaux + '\n' + str(protSeq)

			print(geneTransalatedPath2)
			Gene_Blast_DB_name2 = CommonFastaFunctions.Create_Blastdb_no_fasta(geneTransalatedPath2, 1,
																			   True, sequence_2_blast)

			Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(perGeneBasepath, 1, True, alleleShortProt)

			allelescores, listShortAllelesNames = reDogetBlastScoreRatios(sequence_2_blast,
																						  alleleIaux,
																						  allelescores,
																						  Gene_Blast_DB_name2,
																						  geneScorePickle,
																						  blastPath,
																						  listShortAllelesNames)

		if alleleIaux not in localAllelesDict.values():
			localAllelesDict[newSeq]=alleleIaux

	auxDict={}
	for lala in localAllelesDict.values():
		lala=str(lala)
		new_id=lala.replace("*","")
		auxDict[lala]=int(new_id)

	newSortedIds = sorted(auxDict, key=auxDict.get)


	for SortedId in newSortedIds:

		for k,v in localAllelesDict.items():
			if v==SortedId:
				finalFastaString += ">" + geneBaseName + "_" + str(v) + "\n" + str(k) + "\n"
				break

	finalShortFastaString=''
	for k, v in shortAlleleDict.items():
		finalShortFastaString += ">" + geneBaseName + "_" + str(k) + "\n" + str(v) + "\n"

	with open(geneFile, 'w') as fG:
		fG.write(finalFastaString)

	with open(shortgeneFile, 'w') as fG:
		fG.write(finalShortFastaString)

	try:
		shutil.rmtree(basepath)
	except Exception as e:
		print(e)
		pass







if __name__ == "__main__":
	main()
