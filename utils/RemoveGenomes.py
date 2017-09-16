#!/usr/bin/env python3

import argparse


def main():
    parser = argparse.ArgumentParser(
        description="This program removes genomes from a tab separated allele profile file")
    parser.add_argument('-i', nargs='?', type=str, help='main list file', required=True)
    parser.add_argument('-l', nargs='?', type=str, help='to remove list file', required=True)
    parser.add_argument('-o', nargs='?', type=str, help='output file', required=True)

    args = parser.parse_args()
    mainListFile = args.i
    toRemoveListFile = args.l
    outputfile = args.o

    FilesToRemove = []
    with open(toRemoveListFile) as f:
        for File in f:
            FilesToRemove.append(File.replace("\n", ''))
    print(FilesToRemove)
    towrite = ''

    with open(mainListFile) as f:
        for line in f:
            # print (line.split("\t"))[0]
            if (line.split("\t"))[0] not in FilesToRemove:
                towrite += line

    with open(outputfile, 'w') as f:
        f.write(towrite)


if __name__ == "__main__":
    main()
