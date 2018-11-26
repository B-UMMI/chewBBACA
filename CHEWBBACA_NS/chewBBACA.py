#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna

try:
    from allelecall import BBACA
    from SchemaEvaluator import ValidateSchema
    from utils import TestGenomeQuality,profile_joiner,init_schema_4_bbaca,Extract_cgAlleles,RemoveGenes,down_schema,down_profiles,send2NS,sync_schema,send_metadata
except ImportError:
    from CHEWBBACA_NS.allelecall import BBACA
    from CHEWBBACA_NS.SchemaEvaluator import ValidateSchema
    from CHEWBBACA_NS.utils import TestGenomeQuality,profile_joiner,init_schema_4_bbaca,Extract_cgAlleles,RemoveGenes,down_schema,down_profiles,send2NS,sync_schema,send_metadata

import CHEWBBACA_NS


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

                if os.path.isdir(genepath):
                    continue

                for allele in SeqIO.parse(genepath, "fasta", generic_dna):
                    break
                list_files.append(os.path.abspath(genepath))
            except Exception as e:
                print (e)
                pass

    return list_files


def allele_call():

    def msg(name=None):
        return ''' chewBBACA.py AlleleCall [AlleleCall ...][-h] -i [I] -g [G] -o [O] --cpu [CPU] [-v] [-b [B]][--bsr [BSR]] [-t [T]] [--fc] [--fr] [--json]
			'''

    parser = argparse.ArgumentParser(
        description="This program call alleles for a set of genomes when provided a schema", usage=msg())
    parser.add_argument('AlleleCall', nargs='+', help='do allele call')
    parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
    parser.add_argument('-o', nargs='?', type=str, help="prefix for the output directory", required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2",
                        required=True)
    parser.add_argument("--contained", help=argparse.SUPPRESS, required=False, action="store_true", default=False)
    parser.add_argument("--CDS", help=argparse.SUPPRESS, required=False, action="store_true", default=False)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)
    parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False, default='blastp')
    parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False, default=0.6)
    parser.add_argument('--st', nargs='?', type=float,
                        help="size threshold, default at 0.2 means alleles with size variation +-20 percent will be tagged as ASM/ALM",
                        required=False, default=0.2)
    parser.add_argument('--ptf', nargs='?', type=str, help="provide the prodigal training file (ptf) path",
                        required=False, default=False)
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

    # try to open as a fasta
    fasta = SeqIO.parse(genomeFiles, "fasta", generic_dna)
    try:
        isFasta = (any(fasta))
    except:
        isFasta = False

    # if is a fasta pass as a list of genomes with a single genome, if not check if is a folder or a txt with a list of paths
    if isFasta == True:
        genomes2call = [os.path.abspath(genomeFiles)]
    else:
        genomes2call = check_if_list_or_folder(genomeFiles)

    if isinstance(genomes2call, list):
        with open("listGenomes2Call.txt", "w") as f:
            for genome in genomes2call:
                f.write(genome + "\n")
        genomes2call = "listGenomes2Call.txt"

    try:
        chosenTrainingFileLocal = os.path.join(os.path.dirname(CHEWBBACA_NS.__file__), "prodigal_training_files",
                                               chosenTrainingFile)
    except:
        chosenTrainingFileLocal = False
        pass

    if os.path.isfile(chosenTrainingFile):
        pass

    elif os.path.isfile(chosenTrainingFileLocal):
        chosenTrainingFile = chosenTrainingFileLocal
        pass
    elif not chosenTrainingFileLocal:
        pass
    else:
        print(str(chosenTrainingFile) + " file not found")
        return

    BBACA.main(genomes2call, genes2call, cpuToUse, gOutFile, BSRTresh, BlastpPath, forceContinue, jsonReport, verbose,
               forceReset, contained, chosenTaxon, chosenTrainingFile, inputCDS, sizeTresh)

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
        description="This program analyses a set of gene files, analyzing the alleles CDS and the length of the alleles per gene",
        usage=msg())
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
    parser.add_argument('--light', help="skip clustal and mafft run", required=False, action="store_true",
                        default=False)

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

    ValidateSchema.main(genes, cpuToUse, outputpath, transTable, threshold, splited, title, logScale,
                        OneBadGeneNotConserved, light)


def test_schema():

    def msg(name=None):
        return ''' chewBBACA.py TestGenomeQuality [TestGenomeQuality ...] [-h] -i [I] -n [N] -t [T] -s [S] [-o [O]] [-v]
                    '''

    def msg(name=None):
        return ''' chewBBACA.py TestGenomeQuality [TestGenomeQuality ...] [-h] -i [I] -n [N] -t [T] -s [S] [-o [O]] [-v]
                    '''

    parser = argparse.ArgumentParser(
        description="This program analyze an allele call raw output matrix, returning info on which genomes are responsible for cgMLST loci loss",
        usage=msg())
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

    TestGenomeQuality.main(pathOutputfile, iterationNumber, thresholdBadCalls, step, out_folder, verbose)


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

    parser = argparse.ArgumentParser(description="This program removes genes from a tab separated allele profile file",usage=msg())
    parser.add_argument('RemoveGenes', nargs='+', help='remove loci from a chewBBACA profile')
    parser.add_argument('-i', nargs='?', type=str, help='main matrix file from which to remove', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='list of genes to remove', required=True)
    parser.add_argument('-o', nargs='?', type=str, help='output file name', required=True)
    parser.add_argument("--inverse", help="list to remove is actually the one to keep", dest='inverse',
                        action="store_true", default=False)

    args = parser.parse_args()
    mainListFile = args.i
    toRemoveListFile = args.g
    outputfileName = args.o
    inverse = args.inverse

    RemoveGenes.main(mainListFile, toRemoveListFile, outputfileName, inverse)

def join_profiles():

    def msg(name=None):
        return ''' chewBBACA.py JoinProfiles [RemoveGenes ...][-h] -p1 -p2 -o [O] 
                    '''

    parser = argparse.ArgumentParser(description="This program joins two profiles, returning a single profile file with the common loci",usage=msg())
    parser.add_argument('JoinProfiles', nargs='+', help='join profiles')
    parser.add_argument('-p1', nargs='?', type=str, help='profile 1', required=True)
    parser.add_argument('-p2', nargs='?', type=str, help='profile 2', required=True)
    parser.add_argument('-o', nargs='?', type=str, help='output file name', required=True)

    args = parser.parse_args()
    profile1 = args.p1
    profile2 = args.p2
    outputFile = args.o

    profile_joiner.main(profile1, profile2, outputFile)

def prep_schema():

    def msg(name=None):
        return ''' chewBBACA.py PrepExternalSchema [PrepExternalSchema ...][-h] -i [I] --cpu [CPU] [-v]
                    '''

    parser = argparse.ArgumentParser(
        description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele",usage=msg())
    parser.add_argument('PrepExternalSchema', nargs='+', help='prepare a schema for chewbbaca')
    parser.add_argument('-i', nargs='?', type=str, help='path to folder containg the schema fasta files ( alternative a list of fasta files)', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

    args = parser.parse_args()

    geneFiles = args.i
    cpu2use = args.cpu

    init_schema_4_bbaca.main(geneFiles, cpu2use)

def download_schema_NS():

    def msg(name=None):
        return ''' chewBBACA.py DownloadSchema [DownloadSchema ...][-h] -i [I] --cpu [CPU] -p [P]
                    '''

    parser = argparse.ArgumentParser(
        description="This program downloads a schema from NS",usage=msg())
    parser.add_argument('DownloadSchema', nargs='+', help='This program downloads a schema from the nomenclature server to be used, it is advised to use the compressed uri insteado of this program as it may be slow')
    parser.add_argument('-i', nargs='?', type=str, help='uri to NS schema', required=True)
    parser.add_argument('-p', nargs='?', type=str, help='path where to down', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument('--bsr', nargs='?', type=float, help="maximum BSR score to consider as short, close to 1 all alleles will be on short making blast slower", required=False, default=0.6)

    args = parser.parse_args()

    uri2schema = args.i
    path2down = args.p
    cpu2use = args.cpu
    maxBsrShort = args.bsr

    #check if path to down is direcotry
    if not os.path.isdir(path2down):
        print(path2down+" is not a folder")
        return
    if os.listdir(path2down):
        print(path2down+" is not a folder")
        return
    down_schema.main(uri2schema, path2down,cpu2use,maxBsrShort)

def sync_schema_NS():

    def msg(name=None):
        return ''' chewBBACA.py SyncSchema [SyncSchema ...][-h] -i [I] --cpu [CPU] -p [P]
                    '''

    parser = argparse.ArgumentParser(
        description="This program syncs a local schema with NS",usage=msg())
    parser.add_argument('SyncSchema', nargs='+', help='Syncronize a local schema with the NS')
    parser.add_argument('-p', nargs='?', type=str, help='path of schema', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False, default=0.6)

    args = parser.parse_args()

    path2schema = args.p
    cpu2use = args.cpu
    bsrTresh= args.bsr
    try:
        sync_schema.main(path2schema,cpu2use,bsrTresh)
    except Exception as e:
        print (e)

def send_NS():

    def msg(name=None):
        return ''' chewBBACA.py Send2NS [Send2NS ...][-h] -s [S] -t [T] -p [P]
                    '''

    parser = argparse.ArgumentParser(description="Send local profile and respective alleles to NS",usage=msg())
    parser.add_argument('Send2NS', nargs='+', help='send profiles and local alleles to NS')
    parser.add_argument('-s', nargs='?', type=str, help='path to schema folder', required=True)
    parser.add_argument('-p', nargs='?', type=str, help='tsv with profile', required=True)
    parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)
    parser.add_argument('-m', nargs='?', type=str, help='tsv with metadata', required=False, default=False)
    parser.add_argument('--mdr', nargs='?', type=str, help='maximum missing data allowed to fail, default 0.5 (50 percent missing data allowed). 1 == all profiles are uploaded even with 100 percent missing data', required=False, default=0.5)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

    args = parser.parse_args()

    profileFile = args.p
    pathSchema = args.s
    token= args.t
    metadata = args.m
    cpu2use = args.cpu
    percentMDallowed=args.mdr

    send2NS.main(profileFile,pathSchema,token,metadata,percentMDallowed,cpu2use)

def send_meta():

    def msg(name=None):
        return ''' chewBBACA.py SendMetadata [SendMetadata ...][-h] -s [S] -t [T] -p [P]
                    '''

    parser = argparse.ArgumentParser(description="send metadata to isolates on the NS",usage=msg())
    parser.add_argument('SendMetadata', nargs='+', help='send metadata to isolates on the NS')
    parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)
    parser.add_argument('-m', nargs='?', type=str, help='tsv with metadata', required=False, default=False)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

    args = parser.parse_args()

    token= args.t
    metadata = args.m
    cpu2use = args.cpu

    send_metadata.main(metadata,cpu2use,token)


def down_prof():

    def msg(name=None):
        return ''' chewBBACA.py DownloadProfiles [DownloadProfiles ...][-h] --sp [SP] --sc [SC] --cpu [CPU]
                    '''

    parser = argparse.ArgumentParser(
        description="Download profiles from the NS",usage=msg())
    parser.add_argument('DownloadProfiles', nargs='+', help='download profiles from NS')
    parser.add_argument('--sp', nargs='?', type=str, help='species uri', required=True)
    parser.add_argument('--sc', nargs='?', type=str, help='schema id', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument('-r', nargs='?', type=str, help='genomes to down profile', required=False, default=None)
    parser.add_argument('-p', nargs='?', type=str, help='profile with already downloaded profiles for that schema', required=False, default=None)
    parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)


    args = parser.parse_args()

    species = args.sp
    schema = args.sc
    cpu2use = args.cpu
    genomes2Down = args.r
    inputProfile = args.p
    token = args.t

    down_profiles.main(species,schema,cpu2use,inputProfile,genomes2Down,token)


def main():

    functions_list = ['CreateSchema', 'AlleleCall', 'SchemaEvaluator', 'TestGenomeQuality', 'ExtractCgMLST','RemoveGenes','PrepExternalSchema','JoinProfiles',
                      'DownloadSchema',"SyncSchema",'Send2NS','DownloadProfiles','SendMetadata']
    desc_list = ['Create a gene by gene schema based on genomes', 'Perform allele call for target genomes',
                 'Tool that builds an html output to better navigate/visualize your schema', 'Analyze your allele call output to refine schemas',
                 'Select a subset of loci without missing data (to be used as PHYLOViZ input)','Remove a provided list of loci from your allele call output',
                 'prepare an external schema to be used by chewBBACA','join two profiles in a single profile file', 'Download schema from NS',
                 'Syncronize a local schema (downloaded from NS) with NS','Send local profile and respective alleles to NS','Download all profiles of a given species for a given schema','send metadata to isolates on the NS']

    version="1.0.1"
    createdBy="Mickael Silva"
    rep="https://github.com/B-UMMI/chewBBACA/tree/chewie_NS"
    contact="mickaelsilva@medicina.ulisboa.pt"

    if len(sys.argv) > 1 and "version" in sys.argv[1]:
        print(version)
        return


    print ("chewBBACA version "+version+" by "+ createdBy+ " at "+ rep+ "\nemail contact: "+ contact)

    try:
        print ("\n")
        if sys.argv[1] == functions_list[1]:
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
            download_schema_NS()
        elif sys.argv[1] == functions_list[9]:
            sync_schema_NS()
        elif sys.argv[1] == functions_list[10]:
            send_NS()
        elif sys.argv[1] == functions_list[11]:
            down_prof()
        elif sys.argv[1] == functions_list[12]:
            send_meta()
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
        i=1
        while i<len(functions_list):
            print (functions_list[i] +" : "+desc_list[i])
            i+=1

if __name__ == "__main__":
    main()
