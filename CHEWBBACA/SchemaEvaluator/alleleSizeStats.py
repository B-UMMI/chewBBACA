#!/usr/bin/env python3
import os
import argparse
from Bio import SeqIO
from operator import itemgetter
import numpy
import math

#~ def main():
#~ 
    #~ parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
    #~ parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
    #~ parser.add_argument('-t', nargs='?', type=float, help='threshold', required=False)
    #~ parser.add_argument('-c', nargs='?', type=bool, help='All alleles must be conserved', required=False)
    #~ parser.add_argument('-r', nargs='?', type=bool, help='Return values', required=False)
#~ 
    #~ args = parser.parse_args()
    #~ genes = args.i
#~ 
#~ 
    #~ try:
        #~ threshold=float(args.t)
    #~ except:
        #~ threshold=0.05
        #~ pass
    #~ try:
        #~ OneNotConserved=bool(args.c)
    #~ except:
        #~ OneNotConserved=False
        #~ pass
#~ 
    #~ try:
        #~ ReturnValues=bool(args.r)
    #~ except:
        #~ ReturnValues=False
        #~ pass
#~ 
    #~ return getStats(genes,threshold,OneNotConserved,ReturnValues)



def main(genes,threshold,OneNotConserved,ReturnValues,logScale,outputpath,split_thresh):

    gene_fp = open( genes, 'r')

    statistics=[]

    genenumber=0
    conservedlengthgenes=[]
    notconservedlengthgenes=[]
    genesWoneAllele=[]
    z=0
    print("Genes with only 1 allele:")
    modaStats=[]
    allsizes=[]
    tabStats=[]
    allAllelesStats=[]
    allNumberAlleles=[]
    allNumberAllelesMean=[]
    allNumberAllelesMedian=[]
    genesList=[]

    htmlgenespath=os.path.join(outputpath,"genes_html/")
    relpath=os.path.relpath(htmlgenespath,outputpath)

    #get stats for each locus file
    for gene in gene_fp:

        gene = gene.rstrip('\n')

        genesList.append(gene)
        #gene_fp2 = HTSeq.FastaReader(gene)
        maxsize=0
        minsize=99999999999
        sizesum=0
        allelenumber=0
        aux=[0,0]
        sizes=[]


        #per gene get all sizes, minimin size, maximum size, media and mode
        for allele in SeqIO.parse(gene, "fasta"):

            allelesize=len(str(allele.seq))
            sizes.append(allelesize)
            if 	allelesize>=maxsize or allelesize<=minsize:
                if allelesize>maxsize:
                    maxsize=allelesize
                if allelesize<minsize:
                    minsize=allelesize
            else:
                aux.append(allelesize)
            sizesum+=allelesize
            allelenumber+=1

        allAllelesStats.append(allelenumber)

        aux[0]=	minsize
        #mean=float(sizesum)/allelenumber

        aux[1]=	maxsize

        sizesnpList=numpy.array(sizes)
        tabStats.append([gene,numpy.amin(sizesnpList),numpy.amax(sizesnpList),numpy.mean(sizesnpList),numpy.std(sizesnpList)])

        i=0

        moda=max(set(sizes), key=sizes.count)
        median=numpy.median(numpy.array(sizes))
        mean=numpy.mean(numpy.array(sizes))

        modaStats.append(moda)


        #get ratio between number of alleles outside conserved threshold
        for size in sizes:

            if (not float(size)> moda*(1+threshold) and not float(size)< moda*(1-threshold)):
                i+=1
        rate = i/float(allelenumber)

        #check if the gene is conserved considering the threshold and the -p parameter
        #get locus with only 1 allele
        if not OneNotConserved and (rate>=1 or (len(sizes)-i)<2) :
            if allelenumber==1:
                z+=1
                genesWoneAllele.append(os.path.join(relpath,(os.path.basename(gene)).replace(".fasta",".html")))

                print(os.path.basename(gene))
            conservedlengthgenes.append(os.path.basename(gene))
        elif OneNotConserved and rate>=1  :
            #print i,allelenumber
            if allelenumber==1:
                z+=1
                genesWoneAllele.append(os.path.join(relpath,(os.path.basename(gene)).replace(".fasta",".html")))

                print(os.path.basename(gene))
            conservedlengthgenes.append(os.path.basename(gene))
        else:
            #notconservedlengthgenes.append(os.path.basename(gene))
            notconservedlengthgenes.append(os.path.join(relpath,(os.path.basename(gene)).replace(".fasta",".html")))



        sizes.append(gene)
        sizes.append(sizes[0])
        sizes.append(median)



        allsizes.append(sizes)
        if logScale:
            allNumberAlleles.append([gene,math.log10(moda),math.log10(allelenumber)])
            allNumberAllelesMean.append([gene,math.log10(mean),math.log10(allelenumber)])
            allNumberAllelesMedian.append([gene,math.log10(median),math.log10(allelenumber)])
        else:
            allNumberAlleles.append([gene,moda,allelenumber])
            allNumberAllelesMean.append([gene,mean,allelenumber])
            allNumberAllelesMedian.append([gene,median,allelenumber])

    print("\n"+str(z)+ " genes with only one allele\n")

    #order genes by median
    sortbymedia=sorted(allsizes, key=itemgetter(-1))
    #sortbymedia.reverse()

    #order genes by number of alleles
    sortbyNumberAlleles=sorted(allNumberAlleles, key=itemgetter(1))
    sortbyNumberAllelesMean=sorted(allNumberAllelesMean, key=itemgetter(1))
    sortbyNumberAllelesMedian=sorted(allNumberAllelesMedian, key=itemgetter(1))

    #sortbyNumberAlleles.reverse()





    orderedlistgene=[]

    for elem in sortbymedia:
        elem.pop(-1)
        elem.pop(-1)
        orderedlistgene.append(elem.pop(-1))


    orderedlistgene2=[]
    sortbyNumberAllelesx=[]
    for elem in sortbyNumberAlleles:
        #orderedlistgene2.append(os.path.basename(elem.pop(0)))
        orderedlistgene2.append(os.path.join(relpath,(os.path.basename(elem.pop(0))).replace(".fasta",".html")))
        sortbyNumberAllelesx.append(elem.pop(0))

    orderedlistgene3=[]
    sortbyNumberAllelesMeanx=[]
    for elem in sortbyNumberAllelesMean:
        #orderedlistgene3.append(os.path.basename(elem.pop(0)))
        orderedlistgene3.append(os.path.join(relpath,(os.path.basename(elem.pop(0))).replace(".fasta",".html")))
        sortbyNumberAllelesMeanx.append(elem.pop(0))

    orderedlistgene4=[]
    sortbyNumberAllelesMedianx=[]
    for elem in sortbyNumberAllelesMedian:
        #orderedlistgene4.append(os.path.basename(elem.pop(0)))
        orderedlistgene4.append(os.path.join(relpath,(os.path.basename(elem.pop(0))).replace(".fasta",".html")))
        sortbyNumberAllelesMedianx.append(elem.pop(0))



    print(str(len(conservedlengthgenes)) +" conserved genes")
    print(str(len(notconservedlengthgenes)) +" not conserved genes")


    print("\nBuilding box plot...")

    if not ReturnValues	:
        with open ("notconserved.txt","w") as f:
            for gene in notconservedlengthgenes:
                f.write(str(gene)+"\n")

    else:

        print("Creating the allele number plot")

        orderedlistgene2_basename=[]
        for elem in orderedlistgene2:
            orderedlistgene2_basename.append(os.path.basename(elem))

        orderedlistgene2_html=[]
        for elem in orderedlistgene2:
            orderedlistgene2_html.append(os.path.join(relpath,(os.path.basename(elem)).replace(".fasta",".html")))


        sortbyNumberAlleles=[item for sublist in sortbyNumberAlleles for item in sublist]
        sortbyNumberAllelesMean=[item for sublist in sortbyNumberAllelesMean for item in sublist]
        sortbyNumberAllelesMedian=[item for sublist in sortbyNumberAllelesMedian for item in sublist]

        #numberallelesplotMedianhtml=[[sortbyNumberAllelesx,sortbyNumberAlleles,orderedlistgene2],[sortbyNumberAllelesMeanx,sortbyNumberAllelesMean,orderedlistgene3],[sortbyNumberAllelesMedianx,sortbyNumberAllelesMedian,orderedlistgene4]]
        numberallelesplotMedianhtml=[[sortbyNumberAllelesx,sortbyNumberAlleles,orderedlistgene2],[sortbyNumberAllelesMeanx,sortbyNumberAllelesMean,orderedlistgene3],[sortbyNumberAllelesMedianx,sortbyNumberAllelesMedian,orderedlistgene4]]


        print("Creating the allele size histogram")

        with open('locus_stats.tsv', "w") as f:
            f.write("Locus\tMode_value\tnumber_alleles\n")
            for elem2 in list(zip(genesList, modaStats,allAllelesStats)):
                f.write(os.path.basename(elem2[0])+"\t"+str(elem2[1])+"\t"+str(elem2[2])+"\n")

        histplothtml=modaStats

        finalList=[]
        finalNameList=[]
        finalLinkList=[]
        subList=[]
        nameSublist=[]
        linkSublist=[]
        i=0
        j=0
        os.path.join(relpath,(os.path.basename(elem)).replace(".fasta",".html"))
        while j<len(sortbymedia):

            if i >= split_thresh:
                #print len(subList)
                finalList.append(subList)
                subList=[]
                finalNameList.append(nameSublist)
                nameSublist=[]
                finalLinkList.append(linkSublist)
                linkSublist=[]
                i=0
                j-=1
            else:
                subList.append(sortbymedia[j])
                nameSublist.append(os.path.basename(orderedlistgene[j]))
                linkSublist.append(os.path.join(relpath,(os.path.basename(orderedlistgene[j])).replace(".fasta",".html")))
                i+=1
            j+=1

        if len(subList)>0:
            finalList.append(subList)
            finalNameList.append(nameSublist)
            finalLinkList.append(linkSublist)

        orderedlistgene=finalNameList
        boxplothtml=finalList

        return notconservedlengthgenes,len(genesList),genesWoneAllele,boxplothtml,histplothtml,numberallelesplotMedianhtml,orderedlistgene,j,finalLinkList,allAllelesStats

if __name__ == "__main__":
    main()
