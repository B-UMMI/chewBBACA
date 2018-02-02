#!/usr/bin/env python3
import csv
import numpy as np
from numpy import array
import argparse
import operator
from collections import Counter
import os

def main(contigsfile,out_folder):

    #~ parser = argparse.ArgumentParser(description="Check if locus are being represented more than once")
    #~ parser.add_argument('-i', nargs='?', type=str, help='contig info file', required=True)
    #~ parser.add_argument('-o', nargs='?', type=str, help="Folder for the analysis files", required=False,default=".")
    #~
    #~
    #~ args = parser.parse_args()
    #~ contigsfile = args.i
    #~ out_folder = args.o

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    with open(contigsfile) as f:
        reader = csv.reader(f, delimiter="\t")
        d = list(reader)

    d2 = array(d)

    genelist= d2[:1,:]
    genelist=genelist.tolist()[0]
    genomeslist= d2[:,:1]
    genomeslist=genomeslist.tolist()

    i=1
    paralogs=[]
    pontuationDict={}
    pontuationDict2={}

    while i<len(genomeslist):

        genomeSchema2=d2[i]
        genomeSchema=genomeSchema2.tolist()

        i+=1

        #per each line(genome) get the contigs+positions that are repeated
        #remove the repeated elements that are not contigs+positions

        equalelem=[x for x, y in Counter(genomeSchema).items() if y > 1]

        if "LNF" in equalelem:
            equalelem.remove("LNF")
        if "LOT3" in equalelem:
            equalelem.remove("LOT3")
        if "LOT5" in equalelem:
            equalelem.remove("LOT5")
        if "LOTSC" in equalelem:
            equalelem.remove("LOTSC")
        if "PLOT3" in equalelem:
            equalelem.remove("PLOT3")
        if "PLOT5" in equalelem:
            equalelem.remove("PLOT5")
        if "PLOTSC" in equalelem:
            equalelem.remove("PLOTSC")
        if "PLOT" in equalelem:
            equalelem.remove("PLOT")
        if "NIPH" in equalelem:
            equalelem.remove("NIPH")
        if "NIPHEM" in equalelem:
            equalelem.remove("NIPHEM")
        if "ALM" in equalelem:
            equalelem.remove("ALM")
        if "ASM" in equalelem:
            equalelem.remove("ASM")
        if "allele incomplete" in equalelem:
            equalelem.remove("allele incomplete")
        if "undefined" in equalelem:
            equalelem.remove("ASM")
        if "small match" in equalelem:
            equalelem.remove("small match")


        #give a +1(paralog) point to a locus per each reapeated contig+position

        for elem in equalelem:

            ii = np.where(genomeSchema2 == elem)[0]

            for index in ii:
                paralogs.append(genelist[index])


                if genelist[index] in pontuationDict:
                    pontuationDict[genelist[index]]+=1
                else:
                    pontuationDict[genelist[index]]=1
        j=0

        #give a +1(problem) point to a locus per each reapeated contig+position
        for elem in genomeSchema:

            if elem == "LNF" or "LOT" in elem or "PLOT" in elem or elem == "NIPH" or elem == "NIPHEM" or elem == "ALM" or elem == "ASM" or elem == "allele incomplete" or elem == "undefined" or elem == "small match":
                if genelist[j] in pontuationDict2:
                    pontuationDict2[genelist[j]]+=1
                else:
                    pontuationDict2[genelist[j]]=1
            j+=1
    ordered = sorted(pontuationDict.items(), key=operator.itemgetter(1))

    paralogs=set(paralogs)
    print("Detected number of paralog loci: "+ str(len(paralogs)))


    #write file with a overrepresented locus per line, the number of times the locus is overrepresented, problems and total of overrepresentation+problems

    with open(os.path.join(out_folder,"RepeatedLoci.txt"), "w") as f:
        f.write("gene\tPC\tNDC\n")
        for k,v in ordered:
            try:
                troubledLocus=str(pontuationDict2[k])
            except:
                troubledLocus="0"
                pass
            f.write( (k+ "\t"+str( v)+"\t"+ troubledLocus+"\n"))


if __name__ == "__main__":
    main()
