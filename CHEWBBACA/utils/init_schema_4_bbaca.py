#!/usr/bin/env python3

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import pickle
import multiprocessing
import shutil
import warnings

try:
	from utils import CommonFastaFunctions
except ImportError:
	from CHEWBBACA.utils import CommonFastaFunctions

#~ from cStringIO import StringIO
try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO


def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return '\n{0}: {1}\n\n'.format('PendingDeprecationWarning',
    							 str(msg))

warnings.formatwarning = custom_formatwarning


def get_Short (gene,auxBar):
	blastPath='blastp'
	genesList=[str(gene)]

	pathtoDir=os.path.join(os.path.dirname(gene),"short")

	try:
		if not os.path.exists(pathtoDir):
			os.makedirs(pathtoDir)
	except Exception as e:
		pass
		#~ print (e)


	for gene in genesList:
		#~ print ("processing " +gene)

		pathtoDir=os.path.join(os.path.dirname(gene),"short")

		shortgene= os.path.join(os.path.dirname(gene),"short",os.path.basename(gene))
		shortgene= shortgene.replace(".fasta","_short.fasta")

		tempgene= os.path.join(os.path.dirname(shortgene),"temp",os.path.basename(gene).replace(".fasta",""))
		tempgeneProt= os.path.join(tempgene,os.path.basename(gene))
		tempgeneProt2= os.path.join(tempgene,os.path.basename(gene))
		tempgeneProt2= tempgeneProt2.replace(".fasta","2.fasta")

		tempgeneProtFasta=''
		tempgeneProt2Fasta=''
		shortfasta=''

		if not os.path.exists(tempgene):
			os.makedirs(tempgene)

		#~ gene_fp2 = HTSeq.FastaReader(gene)

		counter=0
		alleleI=0
		var={}

		geneScorePickle=shortgene+'_bsr.txt'
		selfscores=[]
		fasta_corrected=''
		total_alleles=0
		error_alleles=0
		corrected_alleles=0
		for allele in SeqIO.parse(gene, "fasta"):
			total_alleles+=1
			try:
				translatedSequence,sortedSeq, originalSeq=translateSeq(str(allele.seq.upper()))

				if not originalSeq:
					fasta_corrected+='>'+str(allele.name)+'\n'+str(sortedSeq) + '\n'
					corrected_alleles+=1
				else:
					fasta_corrected+='>'+str(allele.name)+'\n'+str(str(allele.seq.upper())) + '\n'
				#~ alleleI=int(((allele.name).split("_"))[-1])
				alleleI=int(((allele.name).split("_"))[-1])
				if counter<1:

					#add first allele as short and calculate self bsr

					counter+=1

					shortfasta='>'+str(allele.name)+'\n'+str(allele.seq.upper()) + '\n'



					tempgeneProtFasta='>'+str(allele.name)+'\n'+str(translatedSequence) + '\n'

					Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(tempgeneProt, 1, True,tempgeneProtFasta)


					# --- get BLAST score ratio --- #

					cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name, evalue=0.001, outfmt=5,num_threads=1)
					out, err = cline(stdin=tempgeneProtFasta)
					psiblast_xml = StringIO(out)
					blast_records = NCBIXML.parse(psiblast_xml)


					allelescore=[]


					for blast_record in blast_records:

						for alignment in blast_record.alignments:

							for match in alignment.hsps:

								allelescore.append(int(match.score))

					selfbsr=float(allelescore[0])/float(allelescore[0])

					var[alleleI]= allelescore[0]

					selfscores.append(allelescore[0])

				else:

					#calculate selfbsr for each allele

					translatedSequence,sortedSeq, originalSeq=translateSeq(str(allele.seq.upper()))


					tempgeneProt2Fasta='>'+str(allele.name)+'\n'+str(translatedSequence) + '\n'


					Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(tempgeneProt2, 1, True,tempgeneProt2Fasta)


					# --- get BLAST score ratio --- #
					cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name, evalue=0.001, outfmt=5,num_threads=1)
					out, err = cline(stdin=tempgeneProt2Fasta)
					psiblast_xml = StringIO(out)
					blast_records = NCBIXML.parse(psiblast_xml)


					allelescore=[]


					for blast_record in blast_records:

						for alignment in blast_record.alignments:

							for match in alignment.hsps:

								allelescore.append(int(match.score))

					selfscore=allelescore[0]
					selfbsr=float(selfscore)/float(selfscore)


					#calculate bsr for the allele vs all previous alleles


					# --- get BLAST score ratio --- #
					cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name, evalue=0.001, outfmt=5,num_threads=1)
					out, err = cline(stdin=tempgeneProtFasta)
					psiblast_xml = StringIO(out)
					blast_records = NCBIXML.parse(psiblast_xml)



					allelescore=[]
					allelescoreId=[]

					bestbsr=0
					bestscore=0
					for blast_record in blast_records:

						for alignment in blast_record.alignments:

							for match in alignment.hsps:

								alleleMatchid=int((blast_record.query_id.split("_"))[-1])

								bsr=float(match.score)/float(selfscores[int(alleleMatchid-1)])
								if bsr>bestbsr and match.score>bestscore and bsr>=0.6:
									bestbsr=bsr
									bestscore=match.score

					if bestbsr>=0.6 and bestbsr<0.7:

						shortfasta+='>'+str(allele.name)+'\n'+str(allele.seq.upper()) + '\n'

						var[alleleI]= selfscore
						selfscores.append(selfscore)

						tempgeneProtFasta+='>'+str(allele.name)+'\n'+str(translatedSequence) + '\n'




			except Exception as e:
				#~ print ('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
				#print (str(allele.name)+" "+str(e))
				#print ("allele not translatable")
				error_alleles+=1

		#~ print ("processed " +gene)
		
		if error_alleles<total_alleles:

			with open(geneScorePickle,'wb') as f:
				pickle.dump(var, f)

			with open(shortgene,'w') as f:
				f.write(shortfasta)
		
		else:
			print ("ATTENTION!!!111 \n"+str(gene)+" has no correct aleles, the file will be removed!!")
			os.remove(gene)
		
		if corrected_alleles>=1 or error_alleles>=1:
			with open(gene,'w') as f:
				f.write(fasta_corrected)

		#print status bar
		if gene in auxBar:
			auxlen=len(auxBar)
			index=auxBar.index(gene)
			print ( "["+"="*index+">"+" "*(auxlen-index)+"] processed "+str(int((float(index)/auxlen)*100))+"%")



	return	True

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
					#~ print ("translation error")
					#~ print (e)
					raise
	return protseq, seq, originalSeq

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	strDNArevC = ''
	for l in strDNA:

		strDNArevC += basecomplement[l]

	return strDNArevC[::-1]

def check_if_list_or_folder(folder_or_list):
	list_files = []
	# check if given a list of genomes paths or a folder to create schema
	try:
		list_files=[]
		gene_fp = open( folder_or_list, 'r')
		for gene in gene_fp:
			gene = gene.strip()
			list_files.append(gene)
		#~ f = open(folder_or_list, 'r')
		#~ f.close()
		#~ list_files = folder_or_list
	except IOError:

		for gene in os.listdir(folder_or_list):
			if not gene.endswith(".fasta"):
				continue
			try:
				genepath = os.path.join(folder_or_list, gene)
				
				if os.path.isdir(genepath):
					continue
				
				for allele in SeqIO.parse(genepath, "fasta"):
					break
				list_files.append(os.path.abspath(genepath))
			except Exception as e:
				print (e)
				pass

	return list_files

def main(geneFiles,cpu2use):

	warnings.warn('init_schema_4_bbaca will be deprecated in '
				  'a future version. To adapt external schemas '
				  'use the PrepExternalSchema process instead. '
				  '\nType "chewBBACA.py PrepExternalSchema -h" '
				  'to get information about adapting external '
				  'schemas.')

	#~ parser = argparse.ArgumentParser(description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fasta with only the most diverse alleles")
	#~ parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
	#~ parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
	#~
	#~ 
	#~ args = parser.parse_args()
	#~ 
	#~ geneFiles = args.i
	#~ cpu2use = args.cpu

	listGenes=[]
	#~ gene_fp = open( geneFiles, 'r')
	#~ for gene in gene_fp:
	#~ gene = gene.rstrip('\n')
	#~ listGenes.append(gene)
	#~ gene_fp.close()	



	listGenes=check_if_list_or_folder(geneFiles)

	#test down bar
	auxBar=[]
	step=int((len(listGenes))/10)+1
	counter=0
	while counter < len(listGenes):
		auxBar.append(listGenes[counter])
		counter+=step


	tempFolder=''

	pool = multiprocessing.Pool(cpu2use)
	for gene in listGenes:

		pool.apply_async(get_Short,args=[str(gene),auxBar])
		tempFolder= os.path.join(os.path.dirname(gene),"short","temp")
	pool.close()
	pool.join()

	try:
		shutil.rmtree(tempFolder)
	except:
		pass


if __name__ == "__main__":

	main()
