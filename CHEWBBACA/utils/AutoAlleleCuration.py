#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO


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
        print(gene)

        #per gene remove the alleles that belong to the genomes to be removed
        firstallele=False
        changed=False
        fasta2Write=''
        for allele in SeqIO.parse(gene, "fasta"):


            alleleGenomeName=((str(allele.name)).split("_"))[-1]

            if not alleleGenomeName in listGenomesRemove:
                fasta2Write+=">"+allele.name+"\n"+allele.seq+"\n"
            else:
                changed=True
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
        
        if changed==True:
            with open(gene,'w') as f:
                f.write(fasta2Write)
        



if __name__ == "__main__":
    main()
