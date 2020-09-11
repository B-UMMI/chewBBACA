#!/usr/bin/env python3
from Bio import SeqIO
import os
import argparse


def get_Short(genesList):
    for gene in genesList:
        # gene = gene.rstrip('\n')
        pathtoDir = os.path.join(os.path.dirname(gene), "short")
        if not os.path.exists(pathtoDir):
            os.makedirs(pathtoDir)
        shortgene = os.path.join(os.path.dirname(gene), "short", os.path.basename(gene))
        shortgene = shortgene.replace(".fasta", "_short.fasta")

        #gene_fp2 = HTSeq.FastaReader(gene)
        for allele in SeqIO.parse(gene, "fasta"):
            fG = open(shortgene, 'w')
            fG.write('>' + str(allele.id) + '\n' + str(allele.seq.upper()) + '\n')
            fG.close()
            break

    return True


def main(geneFiles):
    #~ parser = argparse.ArgumentParser(
        #~ description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele")
    #~ parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
#~ 
    #~ args = parser.parse_args()

    geneFiles = args.i
    listGenes = []
    gene_fp = open(geneFiles, 'r')
    for gene in gene_fp:
        gene = gene.rstrip('\n')
        listGenes.append(gene)
    gene_fp.close()

    get_Short(listGenes)


if __name__ == "__main__":
    main()
