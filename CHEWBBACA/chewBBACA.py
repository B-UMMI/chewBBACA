#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna

try:
	from allelecall import BBACA
	from createschema import PPanGen
	from SchemaEvaluator import ValidateSchema
	from utils import TestGenomeQuality,profile_joiner,init_schema_4_bbaca,uniprot_find,Extract_cgAlleles,RemoveGenes
except ImportError:
	from CHEWBBACA.allelecall import BBACA
	from CHEWBBACA.createschema import PPanGen
	from CHEWBBACA.SchemaEvaluator import ValidateSchema
	from CHEWBBACA.utils import TestGenomeQuality,profile_joiner,init_schema_4_bbaca,uniprot_find,Extract_cgAlleles,RemoveGenes

#~ from allelecall import CommonFastaFunctions,callAlleles_protein3,BBACA


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


def create_schema():

    def msg(name=None):                                                            
        return ''' chewBBACA.py CreateSchema [CreateSchema ...] [-h] -i [I] -o [O] --cpu [CPU] [-b [B]] [--bsr [BSR]] [-t [T]] [-v] [-l [L]]'''

    parser = argparse.ArgumentParser(description="This program creates a schema when provided the genomes",usage=msg())
    parser.add_argument('CreateSchema', nargs='+', help='create a schema')
    parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
    parser.add_argument('-o', nargs='?', type=str, help="Name of the output folder", required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2",
                        required=True)
    parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False, default='blastp')
    parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR similarity", required=False, default=0.6)
    parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False, default=False)
    parser.add_argument('--ptf', nargs='?', type=str, help="provide your own prodigal training file (ptf) path", required=False, default=False)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)
    parser.add_argument('-l', nargs='?', type=int, help="minimum bp locus lenght", required=False, default=201)

    args = parser.parse_args()

    genomeFiles = args.i
    cpuToUse = args.cpu
    outputFile = args.o
    BlastpPath = args.b
    bsr = args.bsr
    chosenTaxon = args.t
    chosenTrainingFile = args.ptf
    verbose = args.verbose
    min_length = args.l

    genomeFiles = check_if_list_or_folder(genomeFiles)
    if isinstance(genomeFiles, list):
        with open("listGenomes2Call.txt", "w") as f:
            for genome in genomeFiles:
                f.write(genome + "\n")
        genomeFiles = "listGenomes2Call.txt"

    
    PPanGen.main(genomeFiles,cpuToUse,outputFile,bsr,BlastpPath,min_length,verbose,chosenTaxon,chosenTrainingFile)

    
    
    try:
        os.remove("listGenomes2Call.txt")
    except:
        pass


def allele_call():

    def msg(name=None):                                                            
        return ''' chewBBACA.py AlleleCall [AlleleCall ...][-h] -i [I] -g [G] -o [O] --cpu [CPU] [-v] [-b [B]][--bsr [BSR]] [--ptf [PTF]] [--fc] [--fr] [--json]
            '''

    parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes when provided a schema",usage=msg())
    parser.add_argument('AlleleCall', nargs='+', help='do allele call')
    parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
    parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2",
                        required=True)
    parser.add_argument("--contained", help=argparse.SUPPRESS, required=False, action="store_true", default=False)
    parser.add_argument("--CDS", help=argparse.SUPPRESS, required=False, action="store_true", default=False)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)
    parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False, default='blastp')
    parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False, default=0.6)
    parser.add_argument('--st', nargs='?', type=float, help="size threshold, default at 0.2 means alleles with size variation of +-20% will be tagged as ASM/ALM", required=False, default=0.2)
    #parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False, default=False)
    parser.add_argument('--ptf', nargs='?', type=str, help="provide the prodigal training file (ptf) path", required=False, default=False)
    parser.add_argument("--fc", help="force continue", required=False, action="store_true", default=False)
    parser.add_argument("--fr", help="force reset", required=False, action="store_true", default=False)
    parser.add_argument("--json", help="report in json file", required=False, action="store_true", default=False)

    args = parser.parse_args()

    genomeFiles = args.i
    genes = args.g
    cpuToUse = args.cpu
    BSRTresh = args.bsr
    sizeTresh = args.st
    verbose = args.verbose
    BlastpPath = args.b
    gOutFile = args.o
    chosenTaxon = False
    chosenTrainingFile = args.ptf
    forceContinue = args.fc
    forceReset = args.fr
    contained = args.contained
    inputCDS = args.CDS
    jsonReport = args.json

    genes2call = check_if_list_or_folder(genes)

    if isinstance(genes2call, list):
        with open("listGenes2Call.txt", "w") as f:
            for genome in genes2call:
                f.write(genome + "\n")
        genes2call = "listGenes2Call.txt"

    
    #try to open as a fasta
    fasta = SeqIO.parse(genomeFiles, "fasta", generic_dna)
    try:
        isFasta=(any(fasta))
    except:
        isFasta=False
    
    #if is a fasta pass as a list of genomes with a single genome, if not check if is a folder or a txt with a list of paths
    if isFasta==True:
        genomes2call=[os.path.abspath(genomeFiles)]
    else:
        genomes2call = check_if_list_or_folder(genomeFiles)
    
    if isinstance(genomes2call, list):
        with open("listGenomes2Call.txt", "w") as f:
            for genome in genomes2call:
                f.write(genome + "\n")
        genomes2call = "listGenomes2Call.txt"

    
    BBACA.main(genomes2call,genes2call,cpuToUse,gOutFile,BSRTresh,BlastpPath,forceContinue,jsonReport,verbose,forceReset,contained,chosenTaxon,chosenTrainingFile,inputCDS,sizeTresh)


    try:
        os.remove("listGenes2Call.txt")
    except:
        pass
    try:
        os.remove("listGenomes2Call.txt")
    except:
        pass

def evaluate_schema():

    def msg(name=None):                                                            
        return ''' chewBBACA.py SchemaEvaluator [SchemaEvaluator ...] [-h] -i [I] [-p] [--log] -l [L] -ta [TA] [-t [T]] [--title [TITLE]] --cpu [CPU] [-s [S]] [--light]
            '''

    parser = argparse.ArgumentParser(
        description="This program analyses a set of gene files, analyzing the alleles CDS and the length of the alleles per gene",usage=msg())
    parser.add_argument('SchemaEvaluator', nargs='+', help='evaluation of a schema')
    parser.add_argument('-i', nargs='?', type=str, help='list genes, directory or .txt file with the full path',
                        required=True)
    parser.add_argument('-p', dest='conserved', action='store_true', help='One bad allele still makes gene conserved',
                        required=False,
                        default=False)
    parser.add_argument('--log', dest='logScale', action='store_true', default=False)
    parser.add_argument('-l', nargs='?', type=str, help='name/location main html file', required=True)
    parser.add_argument('-ta', nargs='?', type=int, help='ncbi translation table', required=False, default=11)
    parser.add_argument('-t', nargs='?', type=float, help='Threshold', required=False, default=0.05)
    parser.add_argument('--title', nargs='?', type=str, help='title on the html', required=False,
                        default="My Analyzed wg/cg MLST Schema - Rate My Schema")
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu to use', required=True)
    parser.add_argument('-s', nargs='?', type=int,
                        help='number of boxplots per page (more than 500 can make the page very slow)', required=False,
                        default=500)
    parser.add_argument('--light', help="skip clustal and mafft run", required=False, action="store_true", default=False)
    
    args = parser.parse_args()
    genes = args.i
    transTable = args.ta
    logScale = args.logScale
    outputpath = args.l
    cpuToUse = args.cpu
    threshold = args.t
    OneBadGeneNotConserved = args.conserved
    splited = args.s
    light = args.light
    title = str(args.title)

    ValidateSchema.main(genes,cpuToUse,outputpath,transTable,threshold,splited,title,logScale,OneBadGeneNotConserved,light)



def test_schema():

    def msg(name=None):                                                            
        return ''' chewBBACA.py TestGenomeQuality [TestGenomeQuality ...] [-h] -i [I] -n [N] -t [T] -s [S] [-o [O]] [-v]
                    '''

    parser = argparse.ArgumentParser(
        description="This program analyze an allele call raw output matrix, returning info on which genomes are responsible for cgMLST loci loss",usage=msg())
    parser.add_argument('TestGenomeQuality', nargs='+', help='test the quality of the genomes on the allele call')
    parser.add_argument('-i', nargs='?', type=str, help='raw allele call matrix file', required=True)
    parser.add_argument('-n', nargs='?', type=int, help='maximum number of iterations', required=True)
    parser.add_argument('-t', nargs='?', type=int, help='maximum threshold of bad calls above 95 percent',
                        required=True)
    parser.add_argument('-s', nargs='?', type=int, help='step between each threshold analysis', required=True)
    parser.add_argument('-o', nargs='?', type=str, help="Folder for the analysis files", required=False, default=".")
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)

    args = parser.parse_args()

    pathOutputfile = args.i
    iterationNumber = args.n
    thresholdBadCalls = args.t
    step = args.s
    out_folder = args.o
    verbose = args.verbose

    TestGenomeQuality.main(pathOutputfile,iterationNumber,thresholdBadCalls,step,out_folder,verbose)



def extract_cgmlst():

    def msg(name=None):                                                            
        return ''' chewBBACA.py ExtractCgMLST [ExtractCgMLST ...] [-h] -i [I] -o [O] [-r [R]] [-g [G]]
                    '''

    parser = argparse.ArgumentParser(description="This program cleans an output file for phyloviz",usage=msg())
    parser.add_argument('ExtractCgMLST', nargs='+', help='clean chewBBACA output')
    parser.add_argument('-i', nargs='?', type=str, help='input file to clean', required=True)
    parser.add_argument('-o', nargs='?', type=str, help='output folder', required=True)
    parser.add_argument('-r', nargs='?', type=str, help='listgenes to remove', required=False, default=False)
    parser.add_argument('-g', nargs='?', type=str, help='listgenomes to remove', required=False, default=False)
    parser.add_argument('-p', nargs='?', type=float, help='maximum presence (e.g 0.95)', required=False, default=1)

    args = parser.parse_args()

    pathOutputfile = args.i
    newfile = args.o
    genes2remove = args.r
    genomes2remove = args.g
    cgMLSTpercent=args.p

    Extract_cgAlleles.main(pathOutputfile,newfile,cgMLSTpercent,genes2remove,genomes2remove)


def remove_genes():

    def msg(name=None):                                                            
        return ''' chewBBACA.py RemoveGenes [RemoveGenes ...][-h] -i [I] -g [G] -o [O] [--inverse]
                    '''

    parser = argparse.ArgumentParser(description="This program removes gens from a tab separated allele profile file",usage=msg())
    parser.add_argument('RemoveGenes', nargs='+', help='remove loci from a chewBBACA profile')
    parser.add_argument('-i', nargs='?', type=str, help='main matrix file from which to remove', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='list of genes to remove', required=True)
    parser.add_argument('-o', nargs='?', type=str, help='output file name', required=True)
    parser.add_argument("--inverse", help="list to remove is actually the one to keep", dest='inverse',
                        action="store_true", default=False)

    args = parser.parse_args()

    args = parser.parse_args()
    mainListFile = args.i
    toRemoveListFile = args.g
    outputfileName = args.o
    inverse = args.inverse
    
    
    RemoveGenes.main(mainListFile,toRemoveListFile,outputfileName,inverse)


def join_profiles():

    def msg(name=None):
        return ''' chewBBACA.py JoinProfiles [RemoveGenes ...][-h] -p1 -p2 -o [O] 
                    '''

    parser = argparse.ArgumentParser(description="This program joins two profiles, returning a single profile file with the common loci",usage=msg())
    parser.add_argument('JoinProfiles', nargs='+', help='join profiles')
    parser.add_argument('-p1', nargs='?', type=str, help='profile 1', required=True)
    parser.add_argument('-p2', nargs='?', type=str, help='profile 2', required=True)
    parser.add_argument('-o', nargs='?', type=str, help='outut file name', required=True)

    args = parser.parse_args()
    profile1 = args.p1
    profile2 = args.p2
    outputFile = args.o

    profile_joiner.main(profile1,profile2,outputFile)


def prep_schema():

    def msg(name=None):                                                            
        return ''' chewBBACA.py PrepExternalSchema [PrepExternalSchema ...][-h] -i [I] --cpu [CPU] [-v]
                    '''

    parser = argparse.ArgumentParser(
        description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele")
    parser.add_argument('PrepExternalSchema', nargs='+', help='prepare a schema for chewbbaca')
    parser.add_argument('-i', nargs='?', type=str, help='path to folder containg the schema fasta files ( alternative a list of fasta files)', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)

    args = parser.parse_args()

    geneFiles = args.i
    cpu2use = args.cpu
    verbose = args.verbose
    
    init_schema_4_bbaca.main(geneFiles,cpu2use)


def find_uniprot():
	
    def msg(name=None):                                                            
        return ''' chewBBACA.py UniprotFinder [UniprotFinder ...][-h] -i [I] -t [T] --cpu [CPU]
                    '''
	
    parser = argparse.ArgumentParser(
        description="This program gets information of each locus created on the schema creation, based on the uniprot database")
    parser.add_argument('UniprotFinder', nargs='+', help='get info about a schema created with chewBBACA')
    parser.add_argument('-i', nargs='?', type=str, help='path to folder containg the schema fasta files ( alternative a list of fasta files)', required=True)
    parser.add_argument('-t', nargs='?', type=str, help='path to proteinID_Genome.tsv file generated', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

    args = parser.parse_args()

    geneFiles = args.i
    tsvFile = args.t
    cpu2use = args.cpu
    
    uniprot_find.main(geneFiles,tsvFile,cpu2use)


def main():
    functions_list = ['CreateSchema', 'AlleleCall', 'SchemaEvaluator', 'TestGenomeQuality', 'ExtractCgMLST','RemoveGenes','PrepExternalSchema','JoinProfiles','UniprotFinder']
    desc_list = ['Create a gene by gene schema based on genomes', 'Perform allele call for target genomes', 'Tool that builds an html output to better navigate/visualize your schema', 'Analyze your allele call output to refine schemas', 'Select a subset of loci without missing data (to be used as PHYLOViZ input)','Remove a provided list of loci from your allele call output','prepare an external schema to be used by chewBBACA','join two profiles in a single profile file','get info about a schema created with chewBBACA']

    version="2.0.10"
    createdBy="Mickael Silva"
    rep="https://github.com/B-UMMI/chewBBACA"
    contact="mickaelsilva@medicina.ulisboa.pt"

    if len(sys.argv)>1 and "version" in sys.argv[1]:
        print (version)
        return
    
    print ("chewBBACA version "+version+" by "+ createdBy+ " at "+ rep+ "\nemail contact: "+ contact)

    try:
        print ("\n")
        if sys.argv[1] == functions_list[0]:
            create_schema()
        elif sys.argv[1] == functions_list[1]:
            allele_call()
        elif sys.argv[1] == functions_list[2]:
            evaluate_schema()
        elif sys.argv[1] == functions_list[3]:
            test_schema()
        elif sys.argv[1] == functions_list[4]:
            extract_cgmlst()
        elif sys.argv[1] == functions_list[5]:
            remove_genes()
        elif sys.argv[1] == functions_list[6]:
            prep_schema()
        elif sys.argv[1] == functions_list[7]:
            join_profiles()
        elif sys.argv[1] == functions_list[8]:
            find_uniprot()
        else:
            print('\n\tUSAGE : chewBBACA.py [module] -h \n')
            print('Select one of the following functions :\n')
            i=0
            while i<len(functions_list):
                print (functions_list[i] +" : "+desc_list[i])
                i+=1
    except Exception as e:
        print (e)
        print('\n\tUSAGE : chewBBACA.py [module] -h \n')
        print('Select one of the following functions :\n')
        i=0
        while i<len(functions_list):
            print (functions_list[i] +" : "+desc_list[i])
            i+=1
if __name__ == "__main__":
    main()
