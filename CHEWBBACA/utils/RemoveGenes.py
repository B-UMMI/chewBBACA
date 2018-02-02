#!/usr/bin/env python3

import csv
import argparse


def main(mainListFile,toRemoveListFile,outputfileName,inverse):
    #~ parser = argparse.ArgumentParser(description="This program removes gens from a tab separated allele profile file")
    #~ parser.add_argument('-i', nargs='?', type=str, help='main matrix file from which to remove', required=True)
    #~ parser.add_argument('-g', nargs='?', type=str, help='list of genes to remove', required=True)
    #~ parser.add_argument('-o', nargs='?', type=str, help='output file name', required=True)
    #~ parser.add_argument("--inverse", help="list to remove is actually the one to keep", dest='inverse',
    #~ action="store_true", default=False)
    #~
    #~ args = parser.parse_args()
    #~ mainListFile = args.i
    #~ toRemoveListFile = args.g
    #~ outputfileName = args.o
    #~ inverse = args.inverse

    if inverse:
        FilesToRemove = ['File', 'FILE', 'file']
    else:
        FilesToRemove = []
    with open(toRemoveListFile) as f:
        for File in f:
            File = File.rstrip('\n')
            File = File.rstrip('\r')
            File = (File.split('\t'))[0]
            FilesToRemove.append(File)
    # print FilesToRemove

    with open(mainListFile, 'r') as tsvin, open(outputfileName + ".tsv", "w") as csvout:
        tsvin = csv.reader(tsvin, delimiter='\t')

        listindextoremove = []
        for firstline in tsvin:
            for gene in firstline:
                if gene in FilesToRemove and not inverse:
                    listindextoremove.append(firstline.index(gene))
                elif gene not in FilesToRemove and inverse:
                    listindextoremove.append(firstline.index(gene))

            for elem in reversed(listindextoremove):
                del firstline[elem]
            csvout.write(('\t'.join(firstline)) + "\n")
            break

        for line in tsvin:
            for elem in reversed(listindextoremove):
                del line[elem]
            csvout.write(('\t'.join(line)) + "\n")


if __name__ == "__main__":
    main()
