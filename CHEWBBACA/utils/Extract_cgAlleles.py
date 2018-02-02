#!/usr/bin/env python3

import csv
import numpy as np
from numpy import array
import argparse
import os


def presAbs(d3, listgenomesRemove,outputfile,cgPercent):
    d2c = np.copy(d3)

    geneslist = d2c[:1, :]
    genomeslist = d2c[:, :1]
    genomeslist = (genomeslist.tolist())

    geneslistaux = []
    for genome in genomeslist:
        geneslistaux.append(genome[0])

    # remove genomes
    listidstoremove = []
    for genome in geneslistaux:
        if genome in listgenomesRemove:
            print ("removed genome : "+genome)
            rowid = geneslistaux.index(genome)
            listidstoremove.append(rowid)
            # geneslistaux.pop(rowid)

    listidstoremove = sorted(listidstoremove, reverse=True)
    for idtoremove in listidstoremove:
        d3 = np.delete(d3, idtoremove, 0)
        d2c = np.delete(d2c, idtoremove, 0)

    print ("building the presence and abscence matrix...")
    row = 1
    row2Del = []
    while row < d2c.shape[0]:
        column = 1
        while column < d2c.shape[1]:
            try:

                aux = int(d2c[row, column])
                if aux > 0:
                    d2c[row, column] = 1
                else:
                    d2c[row, column] = 0
                    row2Del.append(int(column))
            except:
                try:
                    aux = str((d2c[row, column])).replace("INF-", "")
                    aux = int(aux)
                    d2c[row, column] = 1
                except Exception as e:
                    d2c[row, column] = 0
                    row2Del.append(int(column))

            column += 1

        row += 1


    if cgPercent <float(1):
        column = 1
        row2Del=[]
        total=int(d2c.shape[0])-1
        while column < d2c.shape[1]:
            L=d2c[:,column][1:].astype(np.int)
            present=np.count_nonzero(L==1)
            percentPresence=(float(present)/float(total))
            if percentPresence <cgPercent:
                row2Del.append(int(column))
            column += 1
    print ("presence and abscence matrix built")

    d2d = d2c.tolist()

    with open(os.path.join(outputfile,"Presence_Abscence.tsv"), "w") as f:

        writer = csv.writer(f, delimiter='	')
        writer.writerows(d2d)

    row2Del = list(set(row2Del))

    return d2c, d3, row2Del


def clean(inputfile, outputfile, totaldeletedgenes, rangeFloat, toremovegenes, toremovegenomes,cgPercent):
    # open the raw file to be clean

    with open(inputfile) as f:
        reader = csv.reader(f, delimiter="\t")
        d = list(reader)

    originald2 = array(d)

    # get presence abscence matrix
    d2, originald2, del2CG = presAbs(originald2, toremovegenomes,outputfile,cgPercent)

    genomeslist = d2[1:, :1]
    geneslist = (d2[:1, 1:])[0]

    originald2 = originald2.T
    d2 = d2.T
    rowid = (d2.shape[0]) - 1
    deleted = 0
    abscenceMatrix = True
    balldel = 0
    cgMLST = 0

    numbergenomes = len(genomeslist)
    pontuationmatrix = [0] * numbergenomes
    lostgenesList = []

    # clean the original matrix, using the information on the presence/abscence matrix
    print ("processing the matrix")
    while rowid > 0:
        columnid = 1
        genomeindex = 0
        print (str(rowid) + "/" + str(d2.shape[0]))
        # ~ print type((d2[rowid,1:])[0]
        if rowid in del2CG:
            originald2 = np.delete(originald2, rowid, 0)
            d2 = np.delete(d2, rowid, 0)
            totaldeletedgenes += 1
            deleted += 1
            # ~ rowid-=1
            columnid = 1

        elif geneslist[rowid - 1] in toremovegenes:

            originald2 = np.delete(originald2, rowid, 0)
            # ~ d2=np.delete(d2, rowid, 0)
            totaldeletedgenes += 1
            deleted += 1
            # ~ rowid-=1
            columnid = 1
        # ~ break
        else:
            cgMLST += 1

        rowid -= 1

    #remove INF and other missing data tags from the profile
    list2Replace=['LNF','PLOT3','PLOT5','ASM','ALM','NIPHEM','NIPH','LOTSC']
    rowid=0

    for row in originald2:
        auxrow=[]
        for elem in row:
            if elem in list2Replace:
                elem=0
            elif "INF-" in elem:
                elem=elem.replace('INF-', '')
            auxrow.append(elem)


        originald2[rowid]=auxrow
        rowid+=1

    originald2 = originald2.T

    #count number of missing data per genome
    rowid=1
    missingDataCount=[["FILE","number of missing data","percentage"]]
    while rowid<originald2.shape[0]:
        mdCount=originald2[rowid].tolist().count("0")
        missingDataCount.append( [originald2[rowid][0],mdCount,float("{0:.2f}".format(float(mdCount)/int(originald2.shape[1])*100))])
        rowid+=1

    geneslist = (originald2[:1, 1:])[0]
    originald2 = originald2.tolist()

    # write the output file

    with open(os.path.join(outputfile,"cgMLST.tsv"), "w") as f:
        writer = csv.writer(f, delimiter='	')
        writer.writerows(originald2)
    with open(os.path.join(outputfile,"cgMLSTschema.txt"), "w") as f:
        for gene in geneslist:
            f.write(gene+"\n")

    with open(os.path.join(outputfile,"mdata_stats.tsv"), "w") as f:
        for stats in missingDataCount:
            f.write(('\t'.join(map(str, stats)))+"\n")

    #~ statswrite += ('\t'.join(map(str, auxList)))

    print ("deleted : %s loci" % totaldeletedgenes)
    print ("total loci remaining : " + str(cgMLST))


def main(pathOutputfile,newfile,percent,genesToRemoveFile,genomesToRemoveFile):

    if not os.path.exists(newfile):
        os.makedirs(newfile)

    genesToRemove = []
    genomesToRemove = []

    if not genesToRemoveFile==False:
        with open(genesToRemoveFile, "r") as fp:

            for geneFile in fp:
                geneFile = geneFile.rstrip('\n')
                geneFile = geneFile.rstrip('\r')
                geneFile = (geneFile.split('\t'))[0]

                genesToRemove.append(geneFile)


    if not genomesToRemoveFile==False:
        with open(genomesToRemoveFile, "r") as fp:

            for genomeFile in fp:
                genomeFile = genomeFile.rstrip('\n')
                genomeFile = genomeFile.rstrip('\r')

                genomesToRemove.append(genomeFile)



    clean(pathOutputfile, newfile, 0, 0.2, genesToRemove, genomesToRemove,percent)


if __name__ == "__main__":
    main()
