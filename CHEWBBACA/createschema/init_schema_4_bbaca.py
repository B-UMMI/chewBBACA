#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import argparse


def get_Short(genesList):
    for gene in genesList:

        pathtoDir = os.path.join(os.path.dirname(gene), "short")
        if not os.path.exists(pathtoDir):
            os.makedirs(pathtoDir)
        shortgene = os.path.join(os.path.dirname(gene), "short", os.path.basename(gene))
        shortgene = shortgene.replace(".fasta", "_short.fasta")

        for allele in SeqIO.parse(gene, "fasta", generic_dna):
            fG = open(shortgene, 'w')
            fG.write('>' + str(allele.id) + '\n' + str(allele.seq.upper()) + '\n')
            fG.close()
            break

    return True


def main(geneFiles):

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
