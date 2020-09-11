#!/usr/bin/env python3
import os
import argparse
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	strDNArevC = ''
	for l in strDNA:

		strDNArevC += basecomplement[l]

	return strDNArevC[::-1]

def translateSeq(DNASeq,transTable):
	seq=DNASeq
	tableid=transTable
	reversedSeq=False
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
					print (e)
					raise ValueError(e)
	return protseq,seq,reversedSeq

def curate(geneFile):

	gene2write=''
	#gene_fp2 = HTSeq.FastaReader(geneFile)
	for allele in SeqIO.parse(geneFile, "fasta"):
		sequence = str(allele.seq.upper())
		name = allele.name
		#per gene remove the alleles that are not CDS
		#for allele in gene_fp2:

		# if allele is not multiple of 3 it's useless to try to translate
		if (len(sequence) % 3 != 0):

			pass
		else:
			try:
				protseq,seq,reversedSeq=translateSeq(sequence, 11)
				gene2write+=">"+name+"\n"+sequence+"\n"

			except Exception as err:
				print(err)


	with open(geneFile, "wb") as f:
		f.write(gene2write)
	return True

def main():

	parser = argparse.ArgumentParser(description="This program removes alleles from a set of gene files, given a list of names of the genomes to be removed")
	parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus", required=True)

	args = parser.parse_args()
	genes = args.i
	cpuToUse=args.cpu

	gene_fp = open( genes, 'r')

	pool = multiprocessing.Pool(cpuToUse)
	for gene in gene_fp:
		print (gene)
		gene = gene.rstrip('\n')

		pool.apply_async(curate,args=[str(gene)])

	pool.close()
	pool.join()




if __name__ == "__main__":
	main()
