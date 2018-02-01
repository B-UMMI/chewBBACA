#!/usr/bin/env python3
import sys
import os


def main(questionDB,directory,genomeFile,nucleotide):

    name = directory + "/" + genomeFile + "_db"

    if nucleotide==True:
        os.system("makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log")
    else:
        os.system("makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log")

    return True


if __name__ == "__main__":
    main()
