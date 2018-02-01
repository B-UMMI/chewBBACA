#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import os
import argparse
import time
import pickle
import shutil
import multiprocessing
import subprocess
try:
	from createschema import CreateSchema,runProdigal
except ImportError:
	from CHEWBBACA.createschema import CreateSchema,runProdigal


def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return "Not found"


def checkGeneStrings(genome1, genome2, newName, basepath, cpu, blastp, createSchemaPath,verbose,bsr):

    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg, end="")
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function

    pathForTemp = os.path.join(basepath, newName)
    if not os.path.exists(pathForTemp):
        os.makedirs(pathForTemp)

    listOfGenomes = [genome1, genome2]

    dictprots = {}
    dictprotsLen = {}
    dictprotsName = {}
    newlistOfCDS = {}
    proteinsEqual = 0
    smallProteins = 0
    protid = 0
    genomeProts = ""
    genomeProtsTrans = ""

    try:

        for genomeFile in listOfGenomes:

            genomename = (os.path.basename(genomeFile)).split(".")
            genomename = genomename[0]

            listOfCDS = {}

            currentCDSDict = {}
            currentGenomeDict = {}
            filepath = os.path.join(basepath, str(os.path.basename(genomeFile)) + "_ORF.txt")
            newfilepath = os.path.join(basepath, str(newName))

            #g_fp = HTSeq.FastaReader(genomeFile)
            for contig in SeqIO.parse(genomeFile, "fasta", generic_dna):
                sequence = str(contig.seq.upper())
                currentGenomeDict[contig.id] = sequence

            # after the first iteration, genomes are already defined by their cds and no longer have a cds dictionary pickle file
            try:
                with open(filepath, 'rb') as f:
                    currentCDSDict = pickle.load(f)
            except:
                for k, v in currentGenomeDict.items():
                    currentCDSDict[k] = [[v]]

            j = 0

            counter = 0
            tsvProtidGenome=""

            for contigTag, value in currentCDSDict.items():

                for protein in value:
                    protid += 1
                    

                    # at first iteration we use the genome file and after a cds only multifasta file
                    try:
                        seq = currentGenomeDict[contigTag][protein[0]:protein[1]].upper()
                        aux=[genomename,contigTag,str(protein[0]),str(protein[1]),str(protid)]
                        tsvProtidGenome+="\n"+'\t'.join(aux)

                    except Exception as e:

                        seq = str(protein[0])

                    try:
                        protseq, orderedSeq = translateSeq(seq)
                        lengthofProt = len(str(protseq))
                    except:
                        verboseprint( str(genome1) + " " + str(genome2))
                        pass
                    # check if any protein with size on dict

                    try:

                        if len(str(protseq)) < 67:
                            smallProteins += 1
                            pass

                        elif dictprotsLen[lengthofProt]:
                            proteinFound = False
                            listproteinsid = dictprotsLen[lengthofProt]
                            for elem in listproteinsid:

                                if protseq == dictprots[elem]:
                                    proteinFound = True
                                    proteinsEqual += 1
                                    break

                            if not proteinFound:
                                dictprotsLen[lengthofProt].append(protid)
                                dictprots[protid] = protseq
                                newlistOfCDS[protid] = orderedSeq
                                try:
                                    protein[1]
                                    idstr = ">" + str(genomename) + "|protein" + str(protid)
                                    idstr2 = ">" + str(genomename) + "|protein" + str(protid)
                                except:
                                    idstr = ">" + str(contigTag)
                                    idstr2 = ">" + str((contigTag.split(" "))[0])
                                genomeProts += idstr + "\n"
                                genomeProtsTrans += idstr2 + "\n"
                                genomeProts += str(orderedSeq) + "\n"
                                genomeProtsTrans += str(protseq) + "\n"
                                dictprotsName[protid] = idstr2
                    except Exception as e:
                        try:    
                            dictprotsLen[lengthofProt] = [protid]
                            dictprots[protid] = protseq
                            newlistOfCDS[protid] = orderedSeq
                            try:
                                protein[1]
                                idstr = ">" + str(genomename) + "|protein" + str(protid)
                                idstr2 = ">" + str(genomename) + "|protein" + str(protid)

                            except:
                                idstr = ">" + str(contigTag)
                                idstr2 = ">" + str((contigTag.split(" "))[0])

                            genomeProts += idstr + "\n"
                            genomeProtsTrans += idstr2 + "\n"
                            genomeProts += str(orderedSeq) + "\n"
                            genomeProtsTrans += str(protseq) + "\n"
                            dictprotsName[protid] = idstr2
                        except:
                            pass

                    else:
                        pass

            listOfCDS = ''
            currentGenomeDict = ''
            currentCDSDict = ''
            with open("proteinID_Genome.tsv", 'a') as f:
                f.write(tsvProtidGenome)

        verboseprint( "Checked equal proteins for: " + str(genome1) + " " + str(genome2))
        verboseprint( "Starting with a total of loci: " + str(protid))
        verboseprint( "equal proteins : " + str(proteinsEqual))
        verboseprint( "small proteins : " + str(smallProteins))

        fastaFile = os.path.join(pathForTemp, newName + ".fasta")
        with open(fastaFile, 'a') as f:
            f.write(genomeProts)

        newlistOfCDS = {}
        genomeProtsTrans = ''
        genomeProts = ''

        # check if any protein is substring of a larger, mantaining the larger one
        # ordering the sequences by length, larger ones first
        auxlist = dictprotsLen.keys()
        auxlist = sorted(auxlist, key=int)
        auxlist = auxlist[::-1]

        finalProtDict = {}
        genomeProtsTrans = ''
        auxprotlist = []
        contained = 0
        finalnumber = 0
        verboseprint( "Looking for contained proteins in : " + str(genome1) + " " + str(genome2))
        counter = 0

        for elem in auxlist:
            counter += 1

            for protid in dictprotsLen[elem]:
                str2 = str(dictprots[protid])

                try:
                    auxprotlist[0]
                    for elem2 in auxprotlist:
                        isContained = False
                        if len(elem2) < len(str2):
                            break
                        else:
                            if str2 in elem2:
                                isContained = True
                                contained += 1
                                break
                    if not isContained:
                        str1 = dictprotsName[protid]
                        genomeProtsTrans += str1 + "\n" + str2 + "\n"
                        finalnumber += 1
                        auxprotlist.append(str2)
                except Exception as e:
                    str1 = dictprotsName[protid]
                    genomeProtsTrans += str1 + "\n" + str2 + "\n"
                    finalnumber += 1
                    auxprotlist.append(str2)

        auxprotlist = []
        dictprots = {}
        dictprotsLen = {}
        dictprotsName = {}

        verboseprint( "number of contained proteins : " + str(contained))
        verboseprint( "total of loci to blast : " + str(finalnumber))

        proteinFile = os.path.join(pathForTemp, newName + "_proteins.fasta")

        with open(proteinFile, 'a') as f:
            f.write(genomeProtsTrans)
        genomeProtsTrans = ''

        # run createschema for the final protogenome
        verboseprint( "running blast will use this number of cpu: " + str(cpu))
        CreateSchema.main(fastaFile,200,cpu,proteinFile,fastaFile,blastp,bsr,verbose)
        #~ proc = subprocess.Popen(
            #~ [createSchemaPath, '-i', fastaFile, '-l', "200", '--cpu', str(cpu), '-p', proteinFile, '-o', fastaFile,
             #~ "-b", blastp], stdout=subprocess.PIPE)
        #~ p_status = proc.wait()
        verboseprint( "finished blast")

        os.remove(proteinFile)


    except Exception as e:
        verboseprint( e )
        return e

    return True


def reverseComplement(strDNA):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def translateSeq(DNASeq):
    seq = DNASeq
    tableid = 11

    #look for ambiguous bases
    try:
        reverseComplement(seq)
    except:
        raise
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        try:
            seq = reverseComplement(seq)
            myseq = Seq(seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)

        except:
            try:
                seq = seq[::-1]
                myseq = Seq(seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
            except:
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                except Exception as e:
                    #print "translation error"
                    #print e
                    raise

    return protseq, seq


def call_proc(cmd):
    p = subprocess.Popen(args=cmd)
    out, err = p.communicate()
    return p


# ================================================ MAIN ================================================ #

def main(genomeFiles,cpuToUse,outputFile,bsr,BlastpPath,min_length,verbose,chosenTaxon,chosenTrainingFile):
    #~ parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes provided a schema")
    #~ parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
    #~ parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
    #~ parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2",
                        #~ required=True)
    #~ parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False, default='blastp')
    #~ parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR similarity", required=False, default=0.6)
    #~ parser.add_argument('-l', nargs='?', type=int, help="minimum bp locus lenght", required=False, default=200)
    #~ parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False, default=False)
    #~ parser.add_argument('--ptf', nargs='?', type=str, help="provide own training file path", required=False, default=False)
    #~ parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        #~ default=False)
#~ 
    #~ args = parser.parse_args()
#~ 
    #~ genomeFiles = args.i
    #~ cpuToUse = args.cpu
    #~ outputFile = args.o
    #~ BlastpPath = args.b
    #~ bsr = args.bsr
    #~ chosenTaxon = args.t
    #~ chosenTrainingFile = args.ptf
    #~ verbose = args.verbose
    #~ min_length = args.l

    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg, end="")
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function


    # avoid user to run the script with all cores available, could impossibilitate any usage when running on a laptop
    if cpuToUse > multiprocessing.cpu_count() - 2:
        print("Warning, you are close to use all your cpus, if you are using a laptop you may be uncapable to perform any action")

    taxonList = {'Campylobacter jejuni': 'trained_campyJejuni.trn',
                 'Acinetobacter baumannii': 'trained_acinetoBaumannii.trn',
                 'Streptococcus agalactiae': 'trained_strepAgalactiae.trn',
                 'Haemophilus influenzae': 'trained_haemoInfluenzae_A.trn',
                 'Yersinia enterocolitica': 'trained_yersiniaEnterocolitica.trn',
                 'Escherichia coli': 'trained_eColi.trn',
                 'Enterococcus faecium': 'trained_enteroFaecium.trn',
                 'Staphylococcus haemolyticus': 'trained_staphHaemolyticus.trn',
                 'Salmonella enterica': 'trained_salmonellaEnterica_enteritidis.trn',
                 'Staphylococcus aureus': 'trained_StaphylococcusAureus.trn',
                 'Streptococcus pneumoniae': 'trained_strepPneumoniae.trn'
                 }
    if isinstance(chosenTaxon, str):
        trainingFolderPAth = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'TrainingFiles4Prodigal'))
        try:
            chosenTaxon = os.path.join(trainingFolderPAth, taxonList[chosenTaxon])

            if os.path.isfile(chosenTaxon):
                print("will use this training file : " + chosenTaxon)
            else:
                print("training file don't exist")
                print(chosenTaxon)
                return "retry"
        except:
            print("Your chosen taxon "+chosenTaxon+" is not attributed, select one from:")
            for elem in taxonList.keys():
                print(elem)
            return "retry"

    if isinstance(chosenTrainingFile, str):
        trainingFolderPAth = os.path.abspath(chosenTrainingFile)
        try:
            chosenTaxon = trainingFolderPAth

            if os.path.isfile(chosenTaxon):
                print ("will use this training file : " + chosenTaxon)
            else:
                print ("training file don't exist "+chosenTaxon)
                return "retry"
        except:
            print ("The training file you provided doesn't exist:")
            print (chosenTaxon)
            return "retry"
    
    scripts_path = os.path.dirname(os.path.realpath(__file__))

    print ("Will use this number of cpus: " + str(cpuToUse))
    print ("Checking all programs are installed")

    print ("Checking Prodigal installed... " + str(which('prodigal')))

    starttime = "\nStarting Script at : " + time.strftime("%H:%M:%S-%d/%m/%Y")
    print (starttime)

    listOfGenomes = []

    fp = open(genomeFiles, 'r')

    for genomeFile in fp:
        genomeFile = genomeFile.rstrip('\n')
        genomeFile = genomeFile.rstrip('\r')
        listOfGenomes.append(genomeFile)

    fp.close()
    listOfGenomes.sort(key=lambda y: y.lower())

    # check if remnant files from previous run exist, prompt user if exists to know if it's his run and want to continue or start a new one

    basepath = os.path.join((os.path.dirname(outputFile)), "temp")
    if not os.path.exists(basepath):
        os.makedirs(basepath)

    # ------------------------------------------------- #
    #           RUN PRODIGAL OVER ALL GENOMES           #
    # ------------------------------------------------- #


    print ("\nStarting Prodigal at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    # Prodigal run on the genomes, one genome per core using n-2 cores (n number of cores)
    print ("chosen taxon :" + str(chosenTaxon))
    pool = multiprocessing.Pool(cpuToUse)
    for genome in listOfGenomes:
        pool.apply_async(runProdigal.main,( str(genome), basepath, str(chosenTaxon)))

    pool.close()
    pool.join()

    print("\nChecking all prodigal processes created the necessary files...")

    listOfORFCreated = []
    for orffile in os.listdir(basepath):
        if orffile.endswith("_ORF.txt"):
            listOfORFCreated.append(orffile)

    if len(listOfGenomes) > len(listOfORFCreated):
        message = "Missing some files from prodigal. " + str(
            (len(listOfGenomes)) - (len(listOfORFCreated))) + " missing files out of " + str(len(listOfGenomes))
        shutil.rmtree(basepath)
        raise ValueError(message)
    else:
        print("All prodigal files necessary were created\n")

    print ("Finishing Prodigal at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    createSchemaPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CreateSchema.py')

    # ---CDS to protein---#

    # translate the genome CDSs, load them into dictionaries and fasta files to be used further ahead
    # listpairs=[]

    pairID = 0
    # while len(processed)<len(toprocess):
    with open("proteinID_Genome.tsv", 'w') as f:
        f.write("Genome\tcontig\tStart\tStop\tprotID")
    while len(listOfGenomes) > 0:

        pair = []
        dictPairs = {}

        for genomeFile in listOfGenomes:
            # toprocess.append(listOfGenomes)
            if len(pair) < 2:
                pair.append(genomeFile)
            else:
                dictPairs[pairID] = pair
                pairID += 1
                pair = []
                pair.append(genomeFile)

        # if total unpair, keep the remainig
        listOfGenomes = []

        if len(pair) == 2:
            dictPairs[pairID] = pair
            pairID += 1

        elif len(pair) > 0:
            listOfGenomes.append(pair[0])

        numberOfPairs = len(dictPairs.items())
        extraCpu = 0
        if numberOfPairs >= cpuToUse:
            pool = multiprocessing.Pool(cpuToUse)
        else:
            pool = multiprocessing.Pool(numberOfPairs)
            extraCpu = cpuToUse - numberOfPairs

        # print dictPairs
        for item in dictPairs.items():
            k = item[0]
            v = item[1]

            newgGenome = "protogenome" + str(k)
            pathFornewgGenome = os.path.join(basepath, newgGenome, newgGenome + ".fasta")
            listOfGenomes.append(pathFornewgGenome)
            extraCpuPerProcess = extraCpu / numberOfPairs
            print( "running analysis for pair : " + str(v[0]) + " " + str(v[1]))
            pool.apply_async(checkGeneStrings,
                             args=[v[0], v[1], newgGenome, basepath, int(extraCpuPerProcess + 1), BlastpPath,
                                   createSchemaPath,verbose,bsr])

        pool.close()
        pool.join()
        verboseprint( "finished running pair analysis")

        if len(listOfGenomes) == 1:

            verboseprint( "___________________\nFinal step : creating the schema")
            lastFile = listOfGenomes.pop()
            CreateSchema.main(lastFile,min_length,cpuToUse,False,outputFile,BlastpPath,bsr,verbose)
            #~ proc = subprocess.Popen(
                #~ [createSchemaPath, '-i', lastFile, '-l', str(min_length), '--cpu', str(cpuToUse), "-b", BlastpPath, "-o",
                 #~ outputFile, "--bsr", str(bsr)])
            #~ p_status = proc.wait()
            verboseprint( "Schema Created sucessfully")

    shutil.rmtree(basepath)

    print (starttime)
    print ("Finished Script at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))


if __name__ == "__main__":
    main()
