#!/usr/bin/env python3
import sys
import os


def main(questionDB,directory,genomeFile):

    name = directory + "/" + genomeFile + "_db"

    try:
        nucleotide = sys.argv[4]
        os.system("makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log")
    except:
        os.system("makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log")

    return True


if __name__ == "__main__":
    main()
