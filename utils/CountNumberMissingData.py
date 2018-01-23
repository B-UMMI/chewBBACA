#!/usr/bin/env python3
import csv
from numpy import array
import argparse


def main():

    parser = argparse.ArgumentParser(description="Check if locus are being represented more than once")
    parser.add_argument('-i', nargs='?', type=str, help='raw file with allele call', required=True)


    args = parser.parse_args()
    contigsfile = args.i


    with open(contigsfile) as f:
        reader = csv.reader(f, delimiter="\t")
        d = list(reader)

    d2 = array(d)

    genelist= d2[:1,:]
    genelist=genelist.tolist()[0]
    genomeslist= d2[:,:1]
    genomeslist=genomeslist.tolist()
    genomeslist2=[]

    i=1
    pontuationDict={}


    while i<len(genomeslist):

        genomeSchema2=d2[i]
        genomeSchema=genomeSchema2.tolist()
        genomeSchema.pop(0)
        genomeSchema= "\t".join(genomeSchema)

        numberLOT=genomeSchema.count("LOT")
        numberLNF=genomeSchema.count("LNF")
        numberNIPH=genomeSchema.count("NIPH")
        numberASM=genomeSchema.count("ASM")
        numberALM=genomeSchema.count("ALM")

        numberMissingData=numberLOT+numberLNF+numberNIPH+numberASM+numberALM+0
        #if numberMissingData>0:
        pontuationDict[(genomeslist[i])[0]]=numberMissingData
        genomeslist2.append((genomeslist[i])[0])
        i+=1

    for k,v in pontuationDict.items():
        print(k,v)


if __name__ == "__main__":
    main()
