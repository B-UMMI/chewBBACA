#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
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
                print(e)
                pass

    return list_files

def sort_fasta (gene,verbose):

    if verbose:
        def verboseprint(*args):
            for arg in args:
                print (arg),
            print()
    else:
        verboseprint = lambda *a: None  # do-nothing function

    newFasta=''
    allelesTranslated=0
    for allele in SeqIO.parse(gene, "fasta", generic_dna):
        try:
            translatedSequence,sequence=translateSeq(allele.seq)
            newFasta+=">"+allele.id+"\n"+str(sequence)+"\n"
            allelesTranslated+=1
            if allelesTranslated==1:
                shortAllele='>' + str(allele.id) + '\n' + str(allele.seq) + '\n'

        except Exception as e:
            verboseprint( (str(e) + " on gene "+gene+" on allele "+str(allele.id)))
            pass

    if allelesTranslated>0:

        pathtoDir = os.path.join(os.path.dirname(gene), "short")
        if not os.path.exists(pathtoDir):
            os.makedirs(pathtoDir)
        shortgene = os.path.join(os.path.dirname(gene), "short", os.path.basename(gene))
        shortgene = shortgene.replace(".fasta", "_short.fasta")
        fG = open(shortgene, 'w')
        fG.write(shortAllele)
        fG.close()

        with open(gene, "w") as f:
            f.write(newFasta)
    else:
        try:
            os.remove(gene)
            print("No alleles represent a CDS. Removed "+ gene)
        except:
            pass


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
                    raise ValueError(e)
    return protseq,seq

def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:

        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def main():
    parser = argparse.ArgumentParser(
        description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele")
    parser.add_argument('-i', nargs='?', type=str, help='path to folder containg the schema fasta files ( alternative a list of fasta files)', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)

    args = parser.parse_args()

    geneFiles = args.i
    cpu2use = args.cpu
    verbose = args.verbose
    
    geneFiles = check_if_list_or_folder(geneFiles)
    if isinstance(geneFiles, list):
        with open("listGenes.txt", "w") as f:
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

    print ("Processing the fastas")
    pool = multiprocessing.Pool(cpu2use)
    for gene in listGenes:

        pool.apply_async(sort_fasta,args=[str(gene),verbose])
    pool.close()
    pool.join()
    
    print ("Done")
    

if __name__ == "__main__":
    main()
