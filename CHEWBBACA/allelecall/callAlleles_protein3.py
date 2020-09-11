#!/usr/bin/env python3
from io import StringIO
from Bio import SeqIO
import sys
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
from collections import Counter
import os
from Bio.Blast import NCBIXML
try:
    from utils import CommonFastaFunctions
except:
    from CHEWBBACA.utils import CommonFastaFunctions
import time
import pickle
import shutil
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


def getBlastScoreRatios(genefile, basepath, doAll, verbose, blastPath):
    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg),
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function

    allelescores = []
    alleleProt = ''
    alleleAllProt = ''
    alleleList = []
    alleleI = 0
    alleleIlist = []
    listAllelesNames = []
    # calculate bsr for each allele
    for allele in SeqIO.parse(genefile, "fasta"):

        # usually first allele name is just > 1 and after that it has > gene_id_genome
        aux = allele.id.split("_")
        if len(aux) < 2:
            alleleI = str(aux[0])
        else:
            alleleI = str(aux[-1])

        # try to translate the allele
        alleleIlist.append(alleleI)
        alleleList.append(str(allele.seq.upper()))
        listAllelesNames.append(allele.id)

        translatedSequence, x = translateSeq(str(allele.seq.upper()))

        if translatedSequence == '':
            print("cannot translate allele on bsr calculation")
            pass

        # calculate BSR for the allele
        else:
            alleleProt = ">" + str(alleleI) + "\n" + str(translatedSequence + "\n")
            alleleAllProt += ">" + str(alleleI) + "\n" + str(translatedSequence + "\n")
            proteinfastaPath = os.path.join(basepath, str(os.path.basename(genefile) + '_protein2.fasta'))

            # new db for each allele to blast it against himself
            Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb_no_fasta(proteinfastaPath, 1, True, alleleProt)

            # if bsr hasn't been calculated, do the BLAST
            if doAll:

                verboseprint("Starting Blast alleles at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

                # --- get BLAST score ratio --- #
                cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name,
                                              evalue=0.001, outfmt=5, num_threads=1)
                out, err = cline(stdin=alleleProt)
                psiblast_xml = StringIO(out)
                blast_records = NCBIXML.parse(psiblast_xml)

                allelescore = 0

                verboseprint("Blasted alleles at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

                for blast_record in blast_records:

                    for alignment in blast_record.alignments:

                        for match in alignment.hsps:
                            allelescores.append(int(match.score))

                geneScorePickle = os.path.abspath(genefile) + '_bsr.txt'
                verboseprint("________")
                var = dict(zip(alleleIlist, allelescores))
                with open(geneScorePickle, 'wb') as f:
                    pickle.dump(var, f)

            # bsr had already been calculated, load it to memory
            else:
                geneScorePickle = os.path.abspath(genefile) + '_bsr.txt'
                with open(geneScorePickle, 'rb') as f:
                    var = pickle.load(f)
                # needs to convert dictionaries that have integer keys
                # to string keys
                var = {str(k): v for k, v in var.items()}

    proteinfastaPath = os.path.join(basepath, str(os.path.basename(genefile) + '_protein.fasta'))
    with open(proteinfastaPath, "w") as f:
        f.write(alleleAllProt)

    # returning all allele BSR scores and list of alleles for this gene
    return var, alleleList, listAllelesNames


def reDogetBlastScoreRatios(sequence, basepath, alleleI, allelescores2, newGene_Blast_DB_name, alleleList2, picklepath,
                            verbose, blastPath, listAllelesNames):
    if verbose:
        def verboseprint(*args):
            for arg in args:
                print(arg),
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function

    alleleProt = ''

    verboseprint("Starting Blast of new alleles to calculate BSR at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    cline = NcbiblastpCommandline(cmd=blastPath, db=newGene_Blast_DB_name, evalue=0.001, outfmt=5, num_threads=1)

    out, err = cline(stdin=sequence)

    allelescore = 0
    psiblast_xml = StringIO(out)
    blast_records = NCBIXML.parse(psiblast_xml)

    verboseprint("Blasted alleles at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    found = False
    matchscore = 0
    for blast_record in blast_records:

        for alignment in blast_record.alignments:

            for match in alignment.hsps:
                matchscore = int(match.score)

    allelescores2[alleleI] = matchscore
    with open(picklepath, 'wb') as f:
        pickle.dump(allelescores2, f)

    return allelescores2, alleleList2, listAllelesNames


def reverseComplement(strDNA):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def translateSeq(DNASeq):
    seq = DNASeq
    reversedSeq = False
    tableid = 11
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        reversedSeq = True
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
                reversedSeq = False
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                except Exception as e:
                    print("translation error")
                    print(e)
                    protseq = ""
    return protseq, seq


# ======================================================== #
#            Allele calling and classification             #
# ======================================================== #
def main(input_file, temppath, blastPath, verbose, bsrTresh, sizeTresh, ns):

    if verbose == 'True':
        verbose = True
    else:
        verbose = False

    argumentList = []
    with open(input_file, 'rb') as f:
        argumentList = pickle.load(f)

    if verbose:
        def verboseprint(*args):

            for arg in args:
                print(arg),
            print
    else:
        verboseprint = lambda *a: None

    geneFile = argumentList[0]

    verboseprint("Using gene: " + str(geneFile))
    shortgeneFile = os.path.join(os.path.dirname(argumentList[0]), "short", os.path.basename(argumentList[0]))
    shortgeneFile = shortgeneFile.replace(".fasta", "_short.fasta")
    genomesList = argumentList[1]
    genesList = argumentList[2]

    newListgenes = []
    with open(genesList, 'r') as gene_fp:
        for gene in gene_fp:
            gene = gene.rstrip('\n')
            gene = gene.rstrip('\r')
            newListgenes.append(gene)

    statusbar = float(newListgenes.index(str(geneFile))) / len(newListgenes)
    locusnumber = (newListgenes.index(str(geneFile)))
    totalocusnumber = len(newListgenes)
    basepath = os.path.join(temppath, os.path.splitext(geneFile)[0])
    newDNAAlleles2Add2Fasta = ''
    newDNAAlleles2Add2shortFasta = ''
    proteinFastaString = ''

    print('\rProcessing ' + os.path.basename(geneFile) + ". Start " + time.strftime("%H:%M:%S-%d/%m/%Y") + " Locus " + str(
          locusnumber) + " of " + str(totalocusnumber) + ". Done " + str(int(statusbar * 100)) + "%.", end="")

    if not os.path.exists(basepath):
        os.makedirs(basepath)

    fullAlleleList = []
    fullAlleleNameList = []
    alleleI = 0
    # get full list of alleles from main gene file and last allele number id
    for allele in SeqIO.parse(geneFile, "fasta"):
        aux = allele.id.split("_")
        if len(aux) < 2:
            alleleI = aux[0]
        else:
            alleleI = aux[-1]

        fullAlleleList.append(str(allele.seq.upper()))
        fullAlleleNameList.append(allele.id)

    resultsList = []
    i = 0
    perfectMatchIdAllele = []
    perfectMatchIdAllele2 = []
    allelescores = {}
    listShortAllelesNames = []

    verboseprint("Getting BSR at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    geneScorePickle = os.path.abspath(shortgeneFile) + '_bsr.txt'

    # check if bsr as arealdy been calculated and recalculate it if necessary

    if os.path.isfile(geneScorePickle):
        allelescores, alleleList, listShortAllelesNames = getBlastScoreRatios(shortgeneFile, basepath, False, verbose,
                                                                              blastPath)

    else:
        allelescores, alleleList, listShortAllelesNames = getBlastScoreRatios(shortgeneFile, basepath, True, verbose,
                                                                              blastPath)

    with open(os.path.join(basepath, str(os.path.basename(shortgeneFile) + '_protein.fasta')), 'r') as myfile:
        proteinFastaString = myfile.read()
        proteinFastaString += "\n"

    verboseprint("Finished BSR at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

    verboseprint("starting allele call blast at: " + time.strftime("%H:%M:%S-%d/%m/%Y"))
    for genomeFile in genomesList:
        verboseprint(genomeFile)
        bestmatch = [0, 0, False, '',
                     0]  # score, score ratio, perfectmatch, key name of the DNA sequence string, allele ID
        currentGenomeDict = {}
        currentCDSDict = {}

        # load the CDS from the genome to a dictionary
        filepath = os.path.join(temppath, str(os.path.basename(genomeFile)) + "_ORF_Protein.txt")

        with open(filepath, 'rb') as f:
            currentCDSDict = pickle.load(f)

        try:
            intersection = set(fullAlleleList).intersection(currentCDSDict.values())
            intersection = list(intersection)
            if len(intersection) > 1:
                perfectMatchIdAllele.append('NIPHEM')
                perfectMatchIdAllele2.append('NIPHEM')
                verboseprint(os.path.basename(genomeFile) + " has " + str(
                    len(intersection)) + " multiple exact match : " + os.path.basename(
                    geneFile) + " MULTIPLE ALLELES as EXACT MATCH")
                raise ValueError("MULTIPLE ALLELES as EXACT MATCH")

            elif len(intersection) == 1:
                alleleStr = intersection[0]
                # it doenst return both keys with equal values
                # ~ elem=currentCDSDict.keys()[currentCDSDict.values().index(alleleStr)]

                elem = [key for key, value in currentCDSDict.items() if value == alleleStr]
                if len(elem) > 1:
                    perfectMatchIdAllele.append('NIPHEM')
                    perfectMatchIdAllele2.append('NIPHEM')
                    verboseprint(os.path.basename(genomeFile) + " has " + str(
                        len(intersection)) + " multiple exact match : " + os.path.basename(
                        geneFile) + " MULTIPLE ALLELES as EXACT MATCH")
                    raise ValueError("MULTIPLE ALLELES as EXACT MATCH")

                contigname = elem[0].split("&")
                matchLocation = contigname[2]
                # starting CDS base need to be +1
                matchLocation = matchLocation.split("-")
                matchLocation = [int(matchLocation[0]) + 1, int(matchLocation[1])]
                contigname = (contigname[0]).replace(">", "")
                alleleName = ''
                alleleMatchid = 0

                alleleName = fullAlleleNameList[fullAlleleList.index(alleleStr)]
                alleleMatchid = (alleleName.split("_"))[-1]

                try:
                    containedInfo = (alleleName.split("_"))[1]
                except:
                    containedInfo = ''

                perfectMatchIdAllele.append(alleleMatchid)

                if matchLocation[0] > matchLocation[1]:
                    perfectMatchIdAllele2.append(
                        str(contigname) + "&" + str(matchLocation[0]) + "-" + str(matchLocation[1]) + "&" + "-")
                else:

                    perfectMatchIdAllele2.append(
                        str(contigname) + "&" + str(matchLocation[0]) + "-" + str(matchLocation[1]) + "&" + "+")

                # check if atributed allele is contained or contains
                if containedInfo == "CD":
                    resultsList.append([(os.path.basename(genomeFile)), str(alleleMatchid), containedInfo.rstrip()])
                elif containedInfo == "CS":
                    resultsList.append([(os.path.basename(genomeFile)), str(alleleMatchid), containedInfo.rstrip()])
                else:
                    pass

                raise ValueError("EQUAL")
        except Exception as e:
            # ~ exc_type, exc_obj, exc_tb = sys.exc_info()
            # ~ fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            # ~ print(exc_tb.tb_lineno)
            # ~ print e
            continue

        else:
            verboseprint("Blasting alleles on genome at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

            blast_out_file = os.path.join(basepath, "blastdbs/" + os.path.basename(geneFile) + '_List.xml')

            Gene_Blast_DB_name = os.path.join(temppath, str(os.path.basename(genomeFile)) + "/" + str(
                os.path.basename(genomeFile)) + "_db")

            cline = NcbiblastpCommandline(cmd=blastPath, db=Gene_Blast_DB_name, evalue=0.001,
                                          outfmt=5, max_target_seqs=10, max_hsps=10, num_threads=1)

            out, err = cline(stdin=proteinFastaString)

            psiblast_xml = StringIO(out)
            blast_records = NCBIXML.parse(psiblast_xml)

            verboseprint("Blasted alleles on genome at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

            alleleSizes = []
            for allele in fullAlleleList:
                alleleSizes.append(len(allele))

            biggestSizeAllele = max(alleleSizes)

            # get mode allele size
            moda = max(set(alleleSizes), key=alleleSizes.count)
            contador = Counter(alleleSizes).most_common()

            # if most common allele size appears 1 time, get first allele size
            if (contador[0])[1] == 1:
                moda = alleleSizes[0]

            try:

                # iterate through the blast results
                for blast_record in blast_records:

                    locationcontigs = []

                    for alignment in blast_record.alignments:

                        # select the best match
                        for match in alignment.hsps:

                            # query id comes with query_id, not name of the allele
                            # the query will always be a representative allele
                            # we get the index of the representative sequence
                            alleleMatchid = int((blast_record.query_id.split("_"))[-1])

                            # query_id starts with 1 and we have to subtract 1 to get index 0
                            # we get the identifier of the representative (including '*') with the int index
                            alleleMatchid2 = (((listShortAllelesNames[alleleMatchid - 1]).split("_"))[-1])

                            scoreRatio = float(match.score) / float(allelescores[alleleMatchid2])

                            cdsStrName = (alignment.title.split(" "))[1]

                            AlleleDNAstr = alleleList[int(alleleMatchid) - 1]
                            verboseprint("BSR : " + str(scoreRatio))

                            if scoreRatio >= bsrTresh:
                                locationcontigs.append(cdsStrName)

                            # select the best match from BLAST results
                            if scoreRatio == 1 and match.score > bestmatch[0]:
                                bestmatch = [match.score, scoreRatio, False, cdsStrName, int(alleleMatchid), match,
                                             len(AlleleDNAstr)]

                            elif (match.score > bestmatch[0] and scoreRatio >= bsrTresh and scoreRatio > bestmatch[1] and bestmatch[2] is False):
                                bestmatch = [match.score, scoreRatio, False,
                                             cdsStrName, int(alleleMatchid), match,
                                             len(AlleleDNAstr)]

                verboseprint("Classifying the match at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

                # if no best match was found it's a Locus Not Found

                # check for ambiguious bases
                if not bestmatch[0] == 0:
                    alleleStr = currentCDSDict[">" + bestmatch[3]]
                    listFoundAmbiguities = []
                    listambiguousBases = ['K', 'M', 'R', 'Y', 'S', 'W', 'B', 'V', 'H', 'D', 'X', 'N', '-', '.']
                    listFoundAmbiguities = [e for e in listambiguousBases if e in alleleStr]

                if bestmatch[0] == 0 or len(listFoundAmbiguities) > 0:

                    ###################
                    # LOCUS NOT FOUND #
                    ###################
                    if bestmatch[0] == 0:
                        perfectMatchIdAllele.append('LNF')
                        perfectMatchIdAllele2.append('LNF')
                        verboseprint("Locus not found, no matches \n")
                    else:

                        perfectMatchIdAllele.append('LNF')
                        perfectMatchIdAllele2.append('LNF')
                        verboseprint("Locus has strange base \n")

                # if more than one BSR >0.6 in two different CDSs it's a Non Paralog Locus
                elif len(list(set(locationcontigs))) > 1:
                    verboseprint("NIPH", "")
                    perfectMatchIdAllele.append('NIPH')
                    perfectMatchIdAllele2.append('NIPH')
                    for elem in locationcontigs:
                        verboseprint(elem)

                # if match with BSR >0.6 and not equal DNA sequences
                else:
                    # load the contig info of the genome to a dictionary
                    for contig in SeqIO.parse(genomeFile, "fasta"):
                        currentGenomeDict[contig.id] = len(str(contig.seq.upper()))

                    match = bestmatch[5]
                    geneLen = bestmatch[6]
                    alleleStr = currentCDSDict[">" + bestmatch[3]]
                    contigname = bestmatch[3]

                    contigname = contigname.split("&")
                    matchLocation = contigname[2]
                    matchLocation = matchLocation.split("-")
                    matchLocation = [int(matchLocation[0]) + 1, matchLocation[1]]
                    contigname = contigname[0]

                    bestMatchContigLen = currentGenomeDict[contigname]

                    protSeq, alleleStr = translateSeq(alleleStr)
                    # get extra space to the right and left between the allele and match and check if it's still inside the contig

                    rightmatchAllele = geneLen - ((int(match.query_end) + 1) * 3)
                    leftmatchAllele = ((int(match.query_start) - 1) * 3)

                    Reversed = False
                    # ~ if Reversed swap left and right contig extra
                    Reversed = False
                    if int(matchLocation[1]) < int(matchLocation[0]):
                        rightmatchContig = bestMatchContigLen - int(matchLocation[0])
                        leftmatchContig = int(matchLocation[1])
                        aux = rightmatchAllele
                        rightmatchAllele = leftmatchAllele
                        leftmatchAllele = aux
                        Reversed = True

                    else:
                        rightmatchContig = bestMatchContigLen - int(matchLocation[1])
                        leftmatchContig = int(matchLocation[0])

                    ###########################
                    # LOCUS ON THE CONTIG TIP #
                    ###########################

                    # check if contig is smaller than the matched allele
                    if leftmatchContig < leftmatchAllele and rightmatchContig < rightmatchAllele:

                        perfectMatchIdAllele.append('LOTSC')
                        perfectMatchIdAllele2.append('LOTSC')

                        verboseprint(match, contigname, geneFile, leftmatchAllele, rightmatchAllele,
                                     "Locus is bigger than the contig \n")

                    elif leftmatchContig < leftmatchAllele:

                        perfectMatchIdAllele.append('PLOT3')
                        perfectMatchIdAllele2.append('PLOT3')

                        verboseprint(match, contigname, geneFile, leftmatchAllele, rightmatchAllele,
                                     "Locus is on the 3' tip of the contig \n")

                    elif rightmatchContig < rightmatchAllele:

                        perfectMatchIdAllele.append('PLOT5')
                        perfectMatchIdAllele2.append('PLOT5')

                        verboseprint(match, contigname, geneFile, leftmatchAllele, rightmatchAllele,
                                     "Locus is on the 5' tip of the contig \n")

                    elif sizeTresh is not None and (float(len(alleleStr)) > moda + (moda * sizeTresh)):

                        verboseprint("Locus is larger than mode", moda, alleleStr)

                        perfectMatchIdAllele.append('ALM')
                        perfectMatchIdAllele2.append('ALM')

                    elif sizeTresh is not None and (float(len(alleleStr)) < moda - (moda * sizeTresh)):

                        verboseprint("Locus is smaller than mode", moda, alleleStr)

                        perfectMatchIdAllele.append('ASM')
                        perfectMatchIdAllele2.append('ASM')

                    else:
                        #######################
                        # ADD INFERRED ALLELE #		# a new allele
                        #######################

                        wasContained = False
                        tagAuxC = 'S'
                        if ns is True:
                            alleleIaux = "*"+str(int(alleleI.replace("*", ""))+1)
                            alleleI = alleleIaux
                        else:
                            alleleI = str(int(alleleI)+1)
                            alleleIaux = alleleI

                        for alleleaux in fullAlleleList:
                            if alleleStr in alleleaux:
                                alleleName = fullAlleleNameList[fullAlleleList.index(alleleaux)]
                                alleleMatchid = (alleleName.split("_"))[-1]
                                tagAuxC = 'CD' + alleleMatchid.rstrip()
                                resultsList.append([(os.path.basename(genomeFile)), str(alleleIaux), tagAuxC])
                                break
                            elif alleleaux in alleleStr:
                                alleleName = fullAlleleNameList[fullAlleleList.index(alleleaux)]
                                alleleMatchid = (alleleName.split("_"))[-1]
                                tagAuxC = 'CS' + alleleMatchid.rstrip()
                                resultsList.append([(os.path.basename(genomeFile)), str(alleleIaux), tagAuxC])
                                break

                        if not wasContained:
                            tagAux = 'INF'

                            perfectMatchIdAllele.append(tagAux + "-" + str(alleleIaux))

                            if not Reversed:
                                perfectMatchIdAllele2.append(str(contigname) + "&" + str(matchLocation[0]) + "-" + str(
                                    matchLocation[1]) + "&" + "+")
                            else:
                                perfectMatchIdAllele2.append(str(contigname) + "&" + str(matchLocation[0]) + "-" + str(
                                    matchLocation[1]) + "&" + "-")

                            verboseprint(
                                "New allele! Adding allele " + tagAux + str(alleleIaux) + " to the database\n")

                            # --- add the new allele to the gene fasta --- #
                            appendAllele = '>' + str((((os.path.basename(geneFile)).split("."))[0]).replace("_",
                                                                                                            "-")) + "_" + tagAuxC + "_" + (
                                               str(os.path.basename(genomeFile))).replace("_", "-") + "_" + time.strftime("%d/%m/%YT%H:%M:%S") + '_' + str(
                                alleleIaux)

                            newDNAAlleles2Add2Fasta += appendAllele + "\n" + alleleStr + '\n'

                            fullAlleleList.append(alleleStr)
                            fullAlleleNameList.append(appendAllele)

                            if float(bestmatch[1]) >= bsrTresh and float(bestmatch[1]) < bsrTresh + 0.1:

                                newDNAAlleles2Add2shortFasta += appendAllele + "\n" + alleleStr + '\n'

                                geneTransalatedPath2 = os.path.join(basepath, str(
                                    os.path.basename(shortgeneFile) + '_protein2.fasta'))
                                geneTransalatedPath = os.path.join(basepath, str(
                                    os.path.basename(shortgeneFile) + '_protein.fasta'))

                                proteinFastaString += '>' + alleleIaux + '\n' + str(protSeq) + '\n'

                                match = bestmatch[5]

                                # --- remake blast DB and recalculate the BSR for the locus --- #
                                alleleList.append(alleleStr)
                                listShortAllelesNames.append(appendAllele)

                                sequence_2_blast = '>' + alleleIaux + '\n' + str(protSeq)
                                Gene_Blast_DB_name2 = CommonFastaFunctions.Create_Blastdb_no_fasta(geneTransalatedPath2, 1, True, sequence_2_blast)

                                verboseprint("Re-calculating BSR at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))
                                allelescores, alleleList, listShortAllelesNames = reDogetBlastScoreRatios(sequence_2_blast,
                                                                                                          basepath,
                                                                                                          alleleIaux,
                                                                                                          allelescores,
                                                                                                          Gene_Blast_DB_name2,
                                                                                                          alleleList,
                                                                                                          geneScorePickle,
                                                                                                          verbose,
                                                                                                          blastPath,
                                                                                                          listShortAllelesNames)
                                verboseprint("Done Re-calculating BSR at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))

            except Exception as e:
                print("some error occurred")
                print(e)
                print('Error on line {}'.format(sys.exc_info()[-1].tb_lineno))
                perfectMatchIdAllele2.append("ERROR")
                perfectMatchIdAllele.append("ERROR")

    # add new alleles to the locus fasta file
    if len(newDNAAlleles2Add2Fasta) > 5:
        with open(geneFile, 'a') as fG:
            fG.write(newDNAAlleles2Add2Fasta)
    if len(newDNAAlleles2Add2shortFasta) > 5:
        with open(shortgeneFile, 'a') as fG:
            fG.write(newDNAAlleles2Add2shortFasta)

    final = (resultsList, perfectMatchIdAllele)
    verboseprint("Finished allele calling at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))
    filepath = os.path.join(temppath, os.path.basename(geneFile) + "_result.txt")
    filepath2 = os.path.join(temppath, os.path.basename(geneFile) + "_result2.txt")
    with open(filepath, 'wb') as f:
        pickle.dump(final, f)
    with open(filepath2, 'wb') as f:
        pickle.dump(perfectMatchIdAllele2, f)
    shutil.rmtree(basepath)
    return True


if __name__ == "__main__":
    main()
