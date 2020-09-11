#!/usr/bin/env python3

import os
import time
import pickle
import shutil
import argparse
import collections

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline

try:
    from createschema import init_schema_4_bbaca
    from utils import CommonFastaFunctions
except:
    from CHEWBBACA.createschema import init_schema_4_bbaca
    from CHEWBBACA.utils import CommonFastaFunctions


def which(program):

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


def reverseComplement(strDNA):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def translateSeq(DNASeq, genename):
    seq = DNASeq
    tableid = 11
    inverted = False

    # look for ambiguous base
    try:
        reverseComplement(seq)
    except:
        protseq = ""

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
                    protseq = ""
    return protseq, seq, inverted


def main(genes, sizethresh, cpuToUse, proteinFIlePath, outputFIlePath,
         BlastpPath, bsr, verbose):

    if verbose:
        def verboseprint(*args):

            for arg in args:
                print(arg),
            print
    else:
        verboseprint = lambda *a: None

    # translate to protein and create new file
    abspath = os.path.abspath(genes)
    filename = os.path.basename(genes)
    abspath = abspath.replace(filename, '')
    proteinfile = os.path.join(abspath, 'proteins.fasta')

    geneDict = {}
    protDict = {}
    orderedprotDict = collections.OrderedDict()
    alreadyIn = []
    totalgenes = 0
    repeatedgenes = 0
    smallgenes = 0
    nottranslatable = 0

    verboseprint("Checking translatability of the loci:\n")

    if not proteinFIlePath:
        with open(proteinfile, "w") as f:

            for gene in SeqIO.parse(genes, "fasta"):
                dnaseq = str(gene.seq.upper())
                protseq, seq, y = translateSeq(dnaseq, gene.id)
                totalgenes += 1
                if len(protseq) > 1:

                    if str(protseq) in alreadyIn:
                        repeatedgenes += 1

                    elif len(str(seq)) < sizethresh:
                        smallgenes += 1

                    else:
                        alreadyIn.append(str(protseq))
                        protname = ">" + str(gene.id) + "\n"

                        f.write(protname + str(protseq) + "\n")
                        protDict[protname] = str(protseq)
                        geneDict[str(gene.name)] = dnaseq
                else:
                    nottranslatable += 1
                    continue

            verboseprint(str(nottranslatable) + " not translatable out of " + str(totalgenes))

            verboseprint("\nChecking if repeated protein sequences:\n")

            orderedprotList = []
            orderedprotList = sorted(protDict.items(), key=lambda x: len(x[1]), reverse=True)

            i = 0
            while i < len(orderedprotList):
                elem = orderedprotList[i]
                orderedprotDict[elem[0]] = elem[1]
                i += 1

        verboseprint(str(repeatedgenes) + " repeated loci out of " + str(totalgenes))
        verboseprint(str(smallgenes) + " loci out of " + str(totalgenes) + " smaller than " + str(sizethresh) + "bp")
        verboseprint("\nprotein file created\n")

        # first step -  remove genes contained in other genes or 100% equal genes
        # list of results - the output of the function
        resultsList = []

        auxDict = {}
        g = 0
        j = 0

        verboseprint("Checking if protein sequences are contained in others...")

        # for each gene from all the annotated genes - starting with an empty dictionary, only add a new gene if the "to be added gene" is not contained or equal to a gene already added to the dictionary
        auxprot = []

        for elem in orderedprotDict.items():

            contained = False

            prot = str(elem[1])
            if any(prot in x for x in auxprot):
                g += 1
                contained = True

            else:
                auxDict[elem[1]] = elem[0]
                auxprot.append(str(elem[1]))

            j += 1
        verboseprint(str(g)+" loci are contained in other genes\n")

        # overwrite the original file, obtaining a new file with unique genes
        with open(proteinfile, "w") as f:
            allsequences = ''
            for k, v in auxDict.items():
                allsequences += v + k + "\n"
            f.write(allsequences)

    else:
        proteinfile = proteinFIlePath
        totalgenes = 0
        smallgenes = 0
        proteinfile = proteinFIlePath
        for gene in SeqIO.parse(genes, "fasta"):
            dnaseq = str(gene.seq.upper())

            protname = ">" + str(gene.id) + "\n"
            geneDict[str(gene.name)] = dnaseq

    verboseprint("Starting Blast")

    geneFile = os.path.abspath(proteinfile)
    Gene_Blast_DB_name = CommonFastaFunctions.Create_Blastdb(geneFile, 1, True)

    geneF = os.path.splitext(geneFile)[0]
    blast_out_file = geneF + '.xml'
    # ------------------------------ RUNNING BLAST ------------------------------ #
    if cpuToUse:
        cline = NcbiblastpCommandline(cmd=BlastpPath, query=geneFile, db=Gene_Blast_DB_name, evalue=0.001,
                                      out=blast_out_file, outfmt=5, num_threads=int(cpuToUse))
    else:
        cline = NcbiblastpCommandline(cmd=BlastpPath, query=geneFile, db=Gene_Blast_DB_name, evalue=0.001,
                                      out=blast_out_file, outfmt=5, num_threads=1)
    blast_records = CommonFastaFunctions.runBlastParser(cline, blast_out_file)
    verboseprint("Finished blast")

    toRemove = []
    genesToKeep = []
    log = ["removed\tcause\texplanation"]
    for blast_record in blast_records:

        allelename = blast_record.query
        allelename = allelename.split(" ")
        allelename = allelename[0]
        alleleLength = len(geneDict[allelename])

        try:
            # if gene A is not on the toRemove list yet, add to genesToKeep list
            if str(blast_record.query) not in toRemove:
                genesToKeep.append(blast_record.query)

                i = 0
                # if first alignement is not against self, gene B is bigger than gene A and very simillar - remove gene A from genesToKeep and add gene B instead
                if not str(blast_record.query) == str((blast_record.alignments[0]).hit_def):
                    genesToKeep.remove(str(blast_record.query))
                    toRemove.append(str(blast_record.query))
                    log.append(str(blast_record.query) + "\t" + str(
                        (blast_record.alignments[0]).hit_def) + "\t" + "2 is first best match")

                    # if gene B is not on the toRemove list, add to genesToKeep list
                    if str((blast_record.alignments[0]).hit_def) not in toRemove:
                        genesToKeep.append(str((blast_record.alignments[0]).hit_def))

                    raise Exception

                selfblastscore = (((blast_record.alignments[0]).hsps)[0]).score

                while i < len(blast_record.alignments):
                    align = blast_record.alignments[i]

                    match = (align.hsps)[0]
                    scoreRatio = float(match.score) / float(selfblastscore)

                    alleleLength2 = len(geneDict[str(align.hit_def)])

                    # if good match and gene B not in toremove list
                    if (scoreRatio > bsr and not str(align.hit_def) == str(blast_record.query) and str(
                            align.hit_def) not in toRemove):

                        # if gene B is bigger than gene A, keep bigger gene B
                        if alleleLength2 > alleleLength:
                            genesToKeep.append(str(align.hit_def))
                            genesToKeep.remove(str(blast_record.query))
                            toRemove.append(str(blast_record.query))
                            log.append(str(blast_record.query) + "\t" + str(
                                align.hit_def) + "\t" + "2 is bigger and bsr >" + str(bsr))

                            raise Exception
                        # else add gene B to toremove list
                        elif str(align.hit_def) in genesToKeep:
                            genesToKeep.remove(str(align.hit_def))
                            toRemove.append(str(align.hit_def))
                            log.append(str(align.hit_def) + "\t" + str(
                                blast_record.query) + "\t" + "2 is bigger and bsr >" + str(bsr))

                    i += 1

            # else gene A is on toRemove list, add all similar genes (not in genesToKeep) list to the toRemove list
            else:

                i = 0
                selfblastscore = 0
                for align in blast_record.alignments:
                    if not (str(align.hit_def) == str(blast_record.query)):
                        selfblastscore = ((align.hsps)[0]).score
                        # print "gene "+str(align.hit_def)+" is larger than gene "+str(blast_record.query)
                        raise Exception

                while i < len(blast_record.alignments):
                    align = blast_record.alignments[i]
                    match = (align.hsps)[0]
                    scoreRatio = float(match.score) / float(selfblastscore)

                    if align.hit_def not in genesToKeep and not str(align.hit_def) == str(
                            blast_record.query) and scoreRatio > bsr:
                        toRemove.append(align.hit_def)
                        log.append(str(align.hit_def) + "\t" + str(
                            blast_record.query) + "\t" + "2 was on the removed list and bsr >" + str(bsr))
                    else:
                        pass

                    i += 1

        except Exception as e:
            # print e
            pass

    genesToKeep = list(set(genesToKeep))
    toRemove = list(set(toRemove))
    s = set(toRemove)
    notcommonToKeep = [x for x in genesToKeep if x not in s]

    pathfiles = os.path.dirname(geneFile)
    pathfiles = pathfiles + "/"
    listfiles = []

    removedparalogs = 0
    removedsize = 0
    totalgenes = 0
    rest = 0
    concatenatedFile = ''
    schema_folder_path = os.path.join(pathfiles, 'schema_seed')

    if not os.path.exists(schema_folder_path) and not proteinFIlePath and not outputFIlePath:
        os.makedirs(schema_folder_path)
    elif not proteinFIlePath and outputFIlePath:
        os.makedirs(outputFIlePath)

    for contig in SeqIO.parse(genes, "fasta"):
        totalgenes += 1
        name2 = contig.id

        if name2 not in toRemove and name2 in genesToKeep:
            if int(len(contig.seq)) > sizethresh:
                namefile = contig.name
                namefile = namefile.replace("|", "_")
                namefile = namefile.replace("_", "-")
                namefile = namefile.replace("(", "")
                namefile = namefile.replace(")", "")
                namefile = namefile.replace("'", "")
                namefile = namefile.replace("\"", "")
                namefile = namefile.replace(":", "")

                if not proteinFIlePath and not outputFIlePath:
                    newFile = os.path.join(schema_folder_path, namefile + ".fasta")
                    listfiles.append(newFile)
                    with open(newFile, "w") as f:
                        f.write(">" + namefile + "_1\n" + str(contig.seq).upper() + "\n")
                elif not proteinFIlePath and outputFIlePath:
                    newFile = os.path.join(outputFIlePath, namefile + ".fasta")
                    listfiles.append(newFile)
                    with open(newFile, "w") as f:
                        f.write(">" + namefile + "_1\n" + str(contig.seq).upper() + "\n")
                else:
                    concatenatedFile += ">" + contig.id + " \n" + str(contig.seq.upper()) + "\n"

                rest += 1

            else:
                removedsize += 1
        else:

            removedparalogs += 1

    if proteinFIlePath and outputFIlePath:
        with open(outputFIlePath, "w") as f:
            f.write(concatenatedFile)
    elif not proteinFIlePath and outputFIlePath:
        init_schema_4_bbaca.get_Short(listfiles)
        verboseprint("\nRemoved "+str(removedparalogs)+" with a high similarity (BSR>"+str(bsr)+")")
        print("\nCreated schema with "+str(rest)+" loci.")
        os.remove(proteinfile)

    # create short folder
    else:
        init_schema_4_bbaca.get_Short(listfiles)
        verboseprint("\nRemoved "+str(removedparalogs)+" with a high similarity (BSR>"+str(bsr)+")")
        print("\nCreated schema with "+str(rest)+" loci.")
        os.remove(proteinfile)

    shutil.rmtree(os.path.join(pathfiles, 'blastdbs'))

    os.remove(blast_out_file)


if __name__ == "__main__":
    main()
