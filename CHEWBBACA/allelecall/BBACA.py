#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import time
import pickle
import shutil
import datetime as dt
import multiprocessing

from Bio import SeqIO
from Bio.Seq import Seq

try:
    from allelecall import callAlleles_protein3
    from utils import (ParalogPrunning,
                       Create_Genome_Blastdb,
                       constants as ct,
                       gene_prediction as gp,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.allelecall import callAlleles_protein3
    from CHEWBBACA.utils import (ParalogPrunning,
                                 Create_Genome_Blastdb,
                                 constants as ct,
                                 gene_prediction as gp,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


def prepGenomes(genomeFile, basepath, verbose, inputCDS):
    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg),
            print()
    else:
        verboseprint = lambda *a: None  # do-nothing function

    listOfCDS = {}
    genomeProts = ""
    currentCDSDict = {}
    currentGenomeDict = {}

    j = 0
    if inputCDS is False:
        # determine file basename
        # splitting too much???
        basename = (os.path.basename(genomeFile)).split('.')[0]
        filepath = os.path.join(basepath, basename + "_ORF.txt")
        with open(filepath, 'rb') as f:
            currentCDSDict = pickle.load(f)

        for contig in SeqIO.parse(genomeFile, "fasta"):
            sequence = str(contig.seq.upper())
            currentGenomeDict[contig.id] = sequence

        for contigTag, value in currentCDSDict.items():

            for protein in value:
                try:
                    seq = currentGenomeDict[contigTag][protein[0]:protein[1]].upper()
                    protseq, inverted, seq = translateSeq(seq, verbose)
                    j += 1
                    if inverted:
                        idstr = ">" + contigTag + "&protein" + str(j) + "&" + str(protein[1]) + "-" + str(protein[0])
                    else:
                        idstr = ">" + contigTag + "&protein" + str(j) + "&" + str(protein[0]) + "-" + str(protein[1])
                    genomeProts += idstr + "\n"
                    listOfCDS[idstr] = seq
                    genomeProts += str(protseq) + "\n"
                except Exception as e:
                    verboseprint((str(e) + " " + str(genomeFile)))
                    pass

    else:
        for contig in SeqIO.parse(genomeFile, "fasta"):
            sequence = str(contig.seq.upper())
            currentGenomeDict[contig.id] = sequence
            try:
                protseq, inverted, seq = translateSeq(sequence, verbose)
                j += 1
                if inverted:
                    idstr = ">" + contig.id + "&protein" + str(j) + "&0-" + str(len(sequence))
                else:
                    idstr = ">" + contig.id + "&protein" + str(j) + "&0-" + str(len(sequence))
                genomeProts += idstr + "\n"
                listOfCDS[idstr] = seq
                genomeProts += str(protseq) + "\n"
            except:
                print(contig.id + " is not translatable to protein, sequence ignored")
                pass

    if j < 2:
        raise ValueError("your genome has something wrong, are you using a genome as a CDS fasta file or vice versa?")

    filepath = os.path.join(basepath, str(os.path.basename(genomeFile)) + "_ORF_Protein.txt")
    with open(filepath, 'wb') as f:
        var = listOfCDS
        pickle.dump(var, f)
    listOfCDS = ''

    filepath = os.path.join(basepath, str(os.path.basename(genomeFile)) + "_Protein.fasta")
    with open(filepath, 'w') as f:
        f.write(genomeProts)
    genomeProts = ''
    var = ''
    currentGenomeDict = ''
    currentCDSDict = ''

    return True


def reverseComplement(strDNA):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def translateSeq(DNASeq, verbose):
    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg),
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function

    seq = DNASeq
    tableid = 11
    inverted = False
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        try:
            seq = reverseComplement(seq)
            myseq = Seq(seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
            inverted = True
        except:
            try:
                seq = seq[::-1]
                myseq = Seq(seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
                inverted = True
            except:
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                    inverted = False
                except Exception as e:
                    verboseprint("translation error")
                    verboseprint(e)
                    raise

    return protseq, inverted, str(seq)


def loci_translation(genesList, listOfGenomes2, verbose):
    gene_fp = open(genesList, 'r')

    genepath = ''
    basepath = ''
    lGenesFiles = []
    argumentsList = []
    noshortgeneFile = []

    for gene in gene_fp:
        gene = gene.rstrip('\n')
        if not gene.endswith(".fasta"):
            continue

        k = 0

        multiple = True
        shortgene = os.path.join(os.path.dirname(gene), "short", os.path.basename(gene))
        shortgene = shortgene.replace(".fasta", "_short.fasta")
        if not os.path.isfile(shortgene):
            noshortgeneFile.append(gene)
            break

        for allele in SeqIO.parse(shortgene, "fasta"):
            sequence = str(allele.seq.upper())
            k += 1
            if len(sequence) % 3 != 0:
                multiple = False
                print("allele " + str(k) + " is not multiple of 3: " + str(gene) + " this gene is to be removed")
                break
            else:
                try:
                    protseq, Inverted, seq = translateSeq(sequence, verbose)
                except:
                    multiple = False
                    print("allele " + str(k) + " is not translating : " + str(gene) + " this gene is to be removed")
                    break
        if multiple:
            lGenesFiles.append(gene)
            genepath = os.path.dirname(gene)
            basepath = os.path.join(genepath, "temp")
            if not os.path.exists(basepath):
                os.makedirs(basepath)
            filepath = os.path.join(basepath, str(os.path.basename(gene)) + "_argList.txt")
            with open(filepath, 'wb') as f:
                var = [gene, listOfGenomes2, genesList]
                pickle.dump(var, f)
            argumentsList.append(filepath)

    gene_fp.close()
    lGenesFiles.sort(key=lambda y: y.lower())
    return genepath, basepath, lGenesFiles, argumentsList, noshortgeneFile


def main(genomes_files, schema_genes, cpu_cores, output_directory,
         blast_score_ratio, blast_path, force_continue, json_report,
         verbose, force_reset, contained, ptf_path, cds_input,
         size_threshold, translation_table, ns, prodigal_mode):

    divideOutput = False
    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%H:%M:%S-%d/%m/%Y')

    listOfGenomes = []
    listOfGenomesBasename = []

    print("Checking if genome files exist...")
    with open(genomes_files, 'r') as fp:
        for genomeFile in fp:

            genomeFile = genomeFile.rstrip('\n')
            genomeFile = genomeFile.rstrip('\r')
            if os.path.isfile(genomeFile):
                listOfGenomes.append(genomeFile)
                listOfGenomesBasename.append(os.path.basename(genomeFile))
            else:
                print("File does not exist, will not be used: " + str(genomeFile))

    if len(listOfGenomes) == 0:
        raise ValueError('ERROR! No usable genome files in ' + str(genomes_files))

    lGenesFiles = []
    print("Checking if gene files exist...")
    with open(schema_genes, 'r') as gene_fp:
        for gene in gene_fp:
            gene = gene.rstrip('\n')
            gene = gene.rstrip('\r')
            if not gene.endswith(".fasta"):
                print(gene)
                continue
            if os.path.isfile(gene):
                lGenesFiles.append(gene)
            else:
                print("File does not exist, will not be used : " + str(gene))

    if len(lGenesFiles) == 0:
        raise ValueError('ERROR! No usable gene files in ' + str(schema_genes))

    # sort the genomes and genes list for an ordered output
    listOfGenomes.sort(key=lambda y: os.path.basename(y).lower())
    listOfGenomesBasename.sort(key=lambda y: y.lower())
    lGenesFiles.sort(key=lambda y: y.lower())

    # create temp folder inside the folder where the first gene is located
    first_gene = lGenesFiles[0]
    genepath = os.path.dirname(first_gene)
    basepath = os.path.join(genepath, "temp")
    testVar = ""
    if os.path.isdir(basepath) and not force_continue and not force_reset:
        testVar = input('\nWe found files from an unfinished run. '
                        'Do you wish to continue from where it stopped?\n'
                        'Answer (Y/yes): ')
    continueRun = False

    if testVar.lower() == "yes" or testVar.lower() == "y" or force_continue:
        if os.path.isdir(basepath):
            continueRun = True

        if force_reset:
            continueRun = True

    if continueRun:
        print("You chose to continue.")
        argumentsList = []
        resultsList = []
        for filee in os.listdir(basepath):
            if filee.endswith("_argList.txt"):
                argumentsList.append(os.path.join(basepath, str(filee)))
            elif filee.endswith("_result.txt"):
                resultsList.append((os.path.join(basepath, str(filee))).replace("_result.txt", "_argList.txt"))

        # remove unfinished directories
        todelFolders = next(os.walk(genepath))[1]
        for foldertodel in todelFolders:
            if "short" in foldertodel or "temp" in foldertodel:
                pass
            else:
                shutil.rmtree(os.path.join(genepath, str(foldertodel)))

    if not continueRun:

        try:
            # user decided not to continue although there was a run that had been stopped before finishing, files will be removed for a new run
            if os.path.isdir(basepath):
                print("You chose to start a new allele call.")
                print("Removing existing files...")
                shutil.rmtree(basepath)

            # translate the loci
            genepath, basepath, lGenesFiles, argumentsList, noShort = loci_translation(schema_genes, listOfGenomes, verbose)

            # starting a fresh allele call process, going to run prodigal for all genomes on the list

            if len(argumentsList) == 0:
                shutil.rmtree(basepath)
                raise ValueError('ERROR! AT LEAST ONE GENE FILE MUST CONTAIN AN ALLELE REPRESENTING A CDS')

            if len(noShort) > 0:
                shutil.rmtree(basepath)
                raise ValueError('ERROR! These loci have no short gene file: ' + str(noShort))

            # ------------------------------------------------- #
            #           RUN PRODIGAL OVER ALL GENOMES           #
            # ------------------------------------------------- #

            # if the input is a draft genome
            if cds_input is False:

                print("\nStarting Prodigal at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))
                print('Training file: {0}'.format(ptf_path))

                # divide input genomes into equal number of sublists for
                # maximum process progress resolution
                prodigal_inputs = im.divide_list_into_n_chunks(listOfGenomes,
                                                               len(listOfGenomes))
                # add common arguments to all sublists
                common_args = [basepath, ptf_path, translation_table,
                               prodigal_mode, gp.main]
                prodigal_inputs = [i+common_args for i in prodigal_inputs]

                # run Prodigal to predict genes
                prodigal_results = mo.map_async_parallelizer(prodigal_inputs,
                                                             mo.function_helper,
                                                             cpu_cores,
                                                             show_progress=True)

                # determine if Prodigal predicted genes for all genomes
                failed_outdir = os.path.dirname(genepath)
                failed, failed_file = gp.check_prodigal_results(prodigal_results,
                                                                failed_outdir)

                if len(failed) > 0:
                    print('\n\nFailed to predict genes for {0} genomes.'.format(len(failed)))
                    print('Info for failed cases stored in: {0}'.format(failed_file))
                    shutil.rmtree(basepath)
                    sys.exit('Please remove invalid inputs and try again.')

                print("\nFinishing Prodigal at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))

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
                    print("All files were created.")

            # ---CDS to protein---#

            # translate the genome CDSs, load them into dictionaries and fasta files to be used further ahead
            print("\nTranslating genomes...")
            pool = multiprocessing.Pool(cpu_cores)
            for genomeFile in listOfGenomes:
                pool.apply_async(prepGenomes, args=[str(genomeFile), basepath, verbose, cds_input])
            pool.close()
            pool.join()

            print('Creating Blast databases for all genomes...\n')
            # creation of the Databases for each genome, one genome per core using n cores
            pool = multiprocessing.Pool(cpu_cores)
            makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)
            for genomeFile in listOfGenomes:
                filepath = os.path.join(basepath, str(os.path.basename(genomeFile)) + "_Protein.fasta")
                os.makedirs(os.path.join(basepath, str(os.path.basename(genomeFile))))

                pool.apply_async(Create_Genome_Blastdb.main,
                                 (makeblastdb_path, filepath, os.path.join(basepath,
                                  str(os.path.basename(genomeFile))), str(os.path.basename(genomeFile)), False))

            pool.close()
            pool.join()

        except Exception as e:
            exc_type, exc_obj, tb = sys.exc_info()
            lineno = tb.tb_lineno
            print(lineno)

            if json_report:
                runReport = {'finalStatus': 'error : ' + str(e)}
                with open('reportStatus.txt', 'w') as outfile:
                    json.dump(runReport, outfile)
            else:
                print(e)
            raise ValueError(e)
    else:

        # user decided to continue the stopped run, finding out which genes were left to run
        argumentsList = list(set(argumentsList) - set(resultsList))
        argumentsList = sorted(argumentsList)

    print
    print("Starting Allele Calling at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    # Run the allele call, one gene per core using n cores
    pool = multiprocessing.Pool(cpu_cores)
    for argList in argumentsList:
        pool.apply_async(callAlleles_protein3.main,
                         (str(argList), basepath, blast_path, str(verbose), blast_score_ratio, size_threshold, ns))

    pool.close()
    pool.join()

    print("\nFinished Allele Calling at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    print("\nWrapping up the results...")

    output = []
    output2 = []
    for gene in lGenesFiles:
        filepath = os.path.join(basepath, os.path.basename(gene) + "_result.txt")
        filepath2 = os.path.join(basepath, os.path.basename(gene) + "_result2.txt")
        with open(filepath, 'rb') as f:
            var = pickle.load(f)
            output.append(var)
        with open(filepath2, 'rb') as f:
            var = pickle.load(f)
            output2.append(var)

    numberOfLoci = len(output[0][1])
    print('{0}\n {1} genomes used for {2} loci'.format('#'*50, numberOfLoci, len(output)))
    numberexactmatches = 0
    for gene in output:
        for gAllele in gene[1]:
            try:
                allele_id = int(gAllele) if '*' not in gAllele else int(gAllele.replace('*', ''))
                numberexactmatches += 1
            except:
                pass

    print('\n Used a BSR of: {0}'.format(blast_score_ratio))
    print('\n {0} exact matches found out of {1}'.format(numberexactmatches, len(output[0][1])*len(lGenesFiles)))
    exact_percentage = float((numberexactmatches*100) / (numberOfLoci*len(lGenesFiles)))
    exact_percentage = '{:.2f}'.format(exact_percentage)
    print('\n {0} percent of exact matches\n{1}'.format(exact_percentage, '#'*50))
    print('\nWriting output files...\n')

    # wrapping up the results into a matrix and save in a file
    try:
        phylovout = []
        phylovout2 = []
        genesnames = []
        statistics = []

        for gene in lGenesFiles:
            genename = gene.split("/")
            genename = genename[len(genename) - 1]
            genesnames.append(genename)
        gene = 0

        containedOutpWrite = ''
        for geneOut in output:
            genome = 0
            alleleschema = []
            while genome < numberOfLoci:

                genename = (geneOut[1][genome]).split("_")
                if len(genename) != 1:
                    alleleschema.append(genename[1])
                else:
                    alleleschema.append(genename[0])
                genome += 1

            phylovout.append(alleleschema)
            if len(geneOut[0]) > 0:
                containedOutpWrite += lGenesFiles[gene] + "\n"
            for contained2 in geneOut[0]:
                containedOutpWrite += contained2[0] + "\t" + contained2[1] + "-->" + contained2[2] + "\n"
            gene += 1

        for geneOut in output2:
            gene = 0
            alleleschema = []
            while gene < len(output2[0]):
                genename = (geneOut[gene])

                alleleschema.append(genename)
                gene += 1
            phylovout2.append(alleleschema)

        genome = 0

        genesHeader = "FILE" + "\t" + ('\t'.join(map(str, genesnames)))
        finalphylovinput = genesHeader
        finalphylovinput2 = genesHeader

        if divideOutput:
            allelesDict = {}
            contigDict = {}
            statsDict = {}

        while genome < len(listOfGenomes):
            auxList = []
            currentGenome = listOfGenomesBasename[genome]
            statsaux = [0] * 7  # EXC INF LNF PLOT NIPH ALM ASM
            finalphylovinput += "\n" + currentGenome + "\t"
            for gene in phylovout:

                val = str(gene[genome])
                auxList.append(val)
                if "INF" in val:
                    statsaux[1] += 1
                elif "LNF" in val:
                    statsaux[2] += 1
                elif "PLOT" in val:
                    statsaux[3] += 1
                elif "NIPH" in val:
                    statsaux[4] += 1
                elif "ALM" in val:
                    statsaux[5] += 1
                elif "ASM" in val:
                    statsaux[6] += 1
                else:
                    statsaux[0] += 1

            if divideOutput:
                allelesDict[currentGenome] = ('\t'.join(map(str, auxList)))
            finalphylovinput += ('\t'.join(map(str, auxList)))
            # append genome profile to list with new profiles to append
            # to schema profiles master file
            genome += 1
            statistics.append(statsaux)

        #########################
        # This removes INF from the matrix
        # finalphylovinput = finalphylovinput.replace("INF-","")

        genome = 0
        while genome < len(listOfGenomes):
            auxList = []
            currentGenome = listOfGenomesBasename[genome]
            finalphylovinput2 += "\n" + currentGenome + "\t"
            for gene in phylovout2:
                val = str(gene[genome])
                auxList.append(val)

            if divideOutput:
                contigDict[currentGenome] = ('\t'.join(map(str, auxList)))
            finalphylovinput2 += ('\t'.join(map(str, auxList)))
            genome += 1

        statsHeader = 'Genome\tEXC\tINF\tLNF\tPLOT\tNIPH\tALM\tASM'
        statswrite = statsHeader
        # create formatted stats for stdout
        width = max([len(g) for g in listOfGenomesBasename]) + 2
        stdout_Header = ('{:<{width}}  {:^5}  {:^5}  {:^5}  {:^5}  {:^5}  {:^5}  {:^5}'
                         ''.format('Genome', 'EXC', 'INF', 'LNF', 'PLOT', 'NIPH', 'ALM', 'ASM', width=width))
        stdout_line = '-'*len(stdout_Header)
        stdout_stats = [stdout_line, stdout_Header, stdout_line]
        genome = 0
        while genome < len(listOfGenomes):
            auxList = []
            currentGenome = listOfGenomesBasename[genome]
            statswrite += "\n" + currentGenome + "\t"
            for k in statistics[genome]:
                auxList.append(str(k))

            if divideOutput:
                statsDict[currentGenome] = ('\t'.join(map(str, auxList)))
            statswrite += ('\t'.join(map(str, auxList)))
            stdout_stats.append('{:<{width}}  {:^5}  {:^5}  {:^5}  {:^5}  {:^5}  {:^5}  {:^5}'
                                ''.format(currentGenome, *auxList, width=width))
            genome += 1

        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        outputfolder = os.path.join(output_directory, "results_" + str(time.strftime("%Y%m%dT%H%M%S")))
        os.makedirs(outputfolder)

        stdout_stats.append(stdout_line)
        stdout_stats = '\n'.join(stdout_stats)
        print(stdout_stats)

        if json_report:
            runReport = {'finalStatus': 'success'}
            with open(os.path.join(outputfolder, "reportStatus.json"), 'w') as outfile:
                json.dump(runReport, outfile)

            aux = []
            runReport = {}
            for allelename in ((finalphylovinput.splitlines()[0]).split('\t'))[1:]:
                aux.append(allelename)
            runReport['header'] = aux

            for line in (finalphylovinput.splitlines())[1:]:
                aux2 = line.split('\t')
                genome = aux2[0]
                runReport[genome] = aux2[1:]

            with open(os.path.join(outputfolder, "results_alleles.json"), 'w') as outfile:
                json.dump(runReport, outfile)

            aux = []
            runReport = {}
            for allelename in ((statswrite.splitlines()[0]).split('\t'))[1:]:
                aux.append(allelename)
            runReport['header'] = aux

            for line in (statswrite.splitlines())[1:]:
                aux2 = line.split('\t')
                genome = aux2[0]
                runReport[genome] = aux2[1:]

            with open(os.path.join(outputfolder, "results_statistics.json"), 'w') as outfile:
                json.dump(runReport, outfile)

        elif not divideOutput:
            with open(os.path.join(outputfolder, "results_alleles.tsv"), 'w') as f:
                f.write(finalphylovinput)

            with open(os.path.join(outputfolder, "results_statistics.tsv"), 'w') as f:
                f.write(str(statswrite))

            with open(os.path.join(outputfolder, "results_contigsInfo.tsv"), 'w') as f:
                f.write(str(finalphylovinput2))
            if contained:
                with open(os.path.join(outputfolder, "results_contained.txt"), 'w') as f:
                    f.write(str(containedOutpWrite))
            with open(os.path.join(outputfolder, "logging_info.txt"), 'w') as f:
                f.write('Started Script at: {0}'.format(start_date_str))
                f.write("\nFinished Script at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))
                f.write("\nNumber of genomes: " + str(len(listOfGenomes)))
                f.write("\nNumber of loci: " + str(len(lGenesFiles)))
                f.write("\nUsed this number of CPU cores: " + str(cpu_cores))
                f.write("\nUsed a bsr of: " + str(blast_score_ratio))

            print('\nChecking the existence of paralog genes...')

            ParalogPrunning.main(os.path.join(outputfolder, "results_contigsInfo.tsv"), outputfolder)

        else:
            for genome in listOfGenomesBasename:
                currentGenome = os.path.splitext(genome)[0]
                perGenomeFolder = os.path.join(outputfolder, currentGenome)
                os.makedirs(perGenomeFolder)
                with open(os.path.join(perGenomeFolder, currentGenome + "_statistics.txt"), 'w') as f:
                    f.write(statsHeader + "\n")
                    f.write(genome)
                    f.write(statsDict[genome])
                with open(os.path.join(perGenomeFolder, currentGenome + "_contigsInfo.txt"), 'w') as f:
                    f.write(genesHeader + "\n")
                    f.write(genome)
                    f.write(contigDict[genome])
                with open(os.path.join(perGenomeFolder, currentGenome + "_alleles.txt"), 'w') as f:
                    f.write(genesHeader + "\n")
                    f.write(genome)
                    f.write(allelesDict[genome])

    except Exception as e:
        exc_type, exc_obj, tb = sys.exc_info()
        lineno = tb.tb_lineno
        print(lineno, e)
        if json_report:
            runReport = {'finalStatus': 'error : ' + str(e) + ' at line: ' + str(lineno)}
            with open(os.path.join(outputfolder, 'reportStatus.json'), 'w') as outfile:
                json.dump(runReport, outfile)

    if len(os.listdir(outputfolder)) > 0:
        # delete all temp files
        shutil.rmtree(basepath)
    else:
        print('Could not create some output files. Will not delete temp files.')


if __name__ == "__main__":
    main()
