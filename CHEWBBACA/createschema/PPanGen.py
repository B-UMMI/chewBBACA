#!/usr/bin/env python3

import os
import time
import pickle
import shutil
import argparse
import subprocess
import datetime as dt
import multiprocessing

from Bio import SeqIO
from Bio.Seq import Seq

try:
    from createschema import CreateSchema
    from utils import runProdigal
except:
    from CHEWBBACA.createschema import CreateSchema
    from CHEWBBACA.utils import runProdigal


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


def checkGeneStrings(genome1, genome2, newName, basepath, cpu, blastp,
                     createSchemaPath, verbose, bsr):

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
            filepath = os.path.join(basepath, '{0}_ORF.txt'.format(os.path.basename(genomeFile)))
            newfilepath = os.path.join(basepath, str(newName))

            for contig in SeqIO.parse(genomeFile, "fasta"):
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
            tsvProtidGenome = ""
            for contigTag, value in currentCDSDict.items():

                for protein in value:
                    protid += 1

                    # at first iteration we use the genome file and after a cds only multifasta file
                    try:
                        seq = currentGenomeDict[contigTag][protein[0]:protein[1]].upper()
                        aux = [genomename, contigTag, str(protein[0]), str(protein[1]), str(protid)]
                        tsvProtidGenome += "\n"+'\t'.join(aux)

                    except Exception as e:

                        seq = str(protein[0])

                    try:
                        protseq, orderedSeq = translateSeq(seq)
                        lengthofProt = len(str(protseq))
                    except:
                        verboseprint(str(genome1) + " " + str(genome2))
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

        verboseprint("Checked equal proteins for: " + str(genome1) + " " + str(genome2))
        verboseprint("Starting with a total of loci: " + str(protid))
        verboseprint("equal proteins : " + str(proteinsEqual))
        verboseprint("small proteins : " + str(smallProteins))

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
        verboseprint("Looking for contained proteins in : " + str(genome1) + " " + str(genome2))
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

        verboseprint("number of contained proteins : " + str(contained))
        verboseprint("total of loci to blast : " + str(finalnumber))

        proteinFile = os.path.join(pathForTemp, newName + "_proteins.fasta")

        with open(proteinFile, 'a') as f:
            f.write(genomeProtsTrans)
        genomeProtsTrans = ''

        # run createschema for the final protogenome
        verboseprint("running blast will use this number of cpu: " + str(cpu))
        CreateSchema.main(fastaFile, 200, cpu, proteinFile, fastaFile, blastp, bsr, verbose)
        verboseprint("finished blast")

        os.remove(proteinFile)

    except Exception as e:
        verboseprint(e)
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

    # look for ambiguous bases
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
                    raise

    return protseq, seq


def call_proc(cmd):
    p = subprocess.Popen(args=cmd)
    out, err = p.communicate()
    return p


def main(genomeFiles, cpuToUse, outputFile, bsr, BlastpPath, min_length,
         verbose, chosenTrainingFile, inputCDS, translation_table, st):

    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg, end="")
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function

    # avoid user to run the script with all cores available, could impossibilitate any usage when running on a laptop
    if cpuToUse > multiprocessing.cpu_count() - 2:
        print('\nWARNING: you provided a --cpu value close to the '
              'maximum number of available CPU cores.\n'
              'This might degrade system performance and lead '
              'to system unresponsiveness.\n')
        time.sleep(2)

    if isinstance(chosenTrainingFile, str):
        trainingFolderPAth = os.path.abspath(chosenTrainingFile)
        try:
            chosenTaxon = trainingFolderPAth

            if os.path.isfile(chosenTaxon):
                print("Prodigal training file: " + chosenTaxon)
            else:
                print("Training file does not exist "+chosenTaxon)
                return "retry"
        except:
            print("The training file you provided does not exist:")
            print(chosenTaxon)
            return "retry"
    else:
        chosenTaxon = False

    scripts_path = os.path.dirname(os.path.realpath(__file__))

    print("Number of CPU cores: " + str(cpuToUse))

    print("\nChecking dependencies...")
    print("Blast installation..." + str(which(str(BlastpPath))))
    print("Prodigal installation..." + str(which('prodigal')))

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%H:%M:%S-%d/%m/%Y')
    print('\nStarted at: {0}\n'.format(start_date_str))

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

    if inputCDS is True:
        CreateSchema.main(listOfGenomes[0], min_length, cpuToUse, False, outputFile, BlastpPath, bsr, verbose)
        shutil.rmtree(basepath)
        return True

    elif inputCDS is False:

        print("Starting Prodigal at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))

        # Prodigal run on the genomes, one genome per core using n-2 cores (n number of cores)
        pool = multiprocessing.Pool(cpuToUse)
        for genome in listOfGenomes:
            pool.apply_async(runProdigal.main, (str(genome), basepath, str(chosenTaxon), translation_table))

        pool.close()
        pool.join()

        print("Finishing Prodigal at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))

        print("\nChecking if Prodigal created all the necessary files...")
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
            print("All files were created.\n")

    createSchemaPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CreateSchema.py')

    # ---CDS to protein---#

    # translate the genome CDSs, load them into dictionaries and fasta files to be used further ahead
    pairID = 0
    with open("proteinID_Genome.tsv", 'w') as f:
        f.write("Genome\tcontig\tStart\tStop\tprotID")
    while len(listOfGenomes) > 0:

        pair = []
        dictPairs = {}

        for genomeFile in listOfGenomes:
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
            print("Running analysis for pair: " + str(v[0]) + " " + str(v[1]))
            pool.apply_async(checkGeneStrings,
                             args=[v[0], v[1], newgGenome, basepath, int(extraCpuPerProcess + 1), BlastpPath,
                                   createSchemaPath, verbose, bsr])

        pool.close()
        pool.join()

        if len(listOfGenomes) == 1:

            verboseprint("___________________\nFinal step : creating the schema")
            lastFile = listOfGenomes.pop()

            CreateSchema.main(lastFile, min_length, cpuToUse, False,
                              outputFile, BlastpPath, bsr, verbose)

            verboseprint("Schema Created sucessfully")

    shutil.rmtree(basepath)

    end_date = dt.datetime.now()
    end_date_str = dt.datetime.strftime(end_date, '%H:%M:%S-%d/%m/%Y')

    delta = end_date - start_date
    minutes, seconds = divmod(delta.total_seconds(), 60)

    print('\nFinished at: {0}'.format(end_date_str))
    print('Elapsed time: {0:.0f}m{1:.0f}s'.format(minutes, seconds))


if __name__ == "__main__":
    main()
