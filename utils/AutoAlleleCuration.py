#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def main():

    parser = argparse.ArgumentParser(description="This program removes alleles from a set of gene files, given a list of names of the genomes to be removed")
    parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='List of genomes files which to remove', required=True)


    args = parser.parse_args()
    genes = args.i
    genomes = args.g

    gene_fp = open( genes, 'r')
    genomes_fp = open( genomes, 'r')
    listGenomesRemove=[]
    for genomes in genomes_fp:
        listGenomesRemove.append(genomes.rstrip('\n'))



    for gene in gene_fp:
        listSeqTokeep=[]
        gene = gene.rstrip('\n')

        #per gene remove the alleles that belong to the genomes to be removed
        firstallele=False
        for allele in SeqIO.parse(gene, "fasta", generic_dna):


            alleleGenomeName=((str(allele.name)).split("_"))[-1]

            if not alleleGenomeName in listGenomesRemove:
                listSeqTokeep.append(allele)
            else:
                print("removed: "+ str(allele.name))

            alleleGenomeName=((str(allele.name)).split("_"))[-1]

            #if to keep only 1st allele

            """if not firstallele:
                listSeqTokeep.append(allele)
                firstallele=True			
            else:
                print "removed: "+ str(allele.name)
            
            #if to remove the 1st allele
            try:
                int(allele.name)
                print "removed: "+ str(allele.name)
            except:
                listSeqTokeep.append(allele)"""
        my_fasta_file = open( gene, "w" )

        for Seq in listSeqTokeep:

            Seq.write_to_fasta_file( my_fasta_file )
        my_fasta_file.close()



if __name__ == "__main__":
    main()
