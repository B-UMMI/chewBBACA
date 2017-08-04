#!/usr/bin/env python
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import argparse
import multiprocessing

def check_if_list_or_folder(folder_or_list):
    list_files = []
    # check if given a list of genomes paths or a folder to create schema
    try:
        f = open(folder_or_list, 'r')
        f.close()
        list_files = folder_or_list
    except IOError:

        for gene in os.listdir(folder_or_list):
            try:
                genepath = os.path.join(folder_or_list, gene)
                for allele in SeqIO.parse(genepath, "fasta", generic_dna):
                    break
                list_files.append(os.path.abspath(genepath))
            except Exception as e:
                print e
                pass

    return list_files

def get_Short (gene):

		
	newFasta=''
	for allele in SeqIO.parse(gene, "fasta", generic_dna):		
		try: 
			translatedSequence,sequence=translateSeq(allele.seq)
			newFasta+=">"+allele.id+"\n"+str(sequence)+"\n"
		except Exception as e:
			print e
			pass
			
	with open(gene, "wb") as f:
		f.write(newFasta)

		
	return	True
	
def translateSeq(DNASeq):
	seq=DNASeq
	reversedSeq=False
	tableid=11
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=tableid,cds=True)
	except:
		reversedSeq=True
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
					print "translation error"
					#~ print e
					raise
	return protseq,seq

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]
def main():
			
	parser = argparse.ArgumentParser(description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele")
	parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)

	args = parser.parse_args()
	geneFiles = args.i
	
	
	geneFiles = check_if_list_or_folder(geneFiles)
	if isinstance(geneFiles, list):
		with open("listGenes.txt", "wb") as f:
			for genome in geneFiles:
				f.write(genome + "\n")
		geneFiles = "listGenes.txt"
	
	
	listGenes=[]
	gene_fp = open( geneFiles, 'r')
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		listGenes.append(gene)
	gene_fp.close()	
	
	try:
		os.remove("listGenes.txt")
	except:
		pass
	
	pool = multiprocessing.Pool(2)
	for gene in listGenes:
		
		pool.apply_async(get_Short,args=[str(gene)])
	pool.close()
	pool.join()
	
		
	
	

if __name__ == "__main__":
    main()
