#!/usr/bin/env python3

from Bio.Blast import NCBIXML
import os


def ensure_dir(f):
    if not os.path.isdir(f):
        # print "No BLAST db dir was found. Creating it..."
        os.makedirs(f)


# overwrite has to be 1 or 0 (True or False)
def Create_Blastdb(questionDB, overwrite, dbtypeProt):
    base = os.path.basename(questionDB)
    dirname = os.path.dirname(questionDB)
    isProt = dbtypeProt

    if len(dirname) == 0:
        dirname = '.'
    basename = os.path.splitext(base)[0]
    ensure_dir(dirname + "/blastdbs")
    name = dirname + "/blastdbs/" + basename + "_db"

    if not os.path.isfile(name + ".nin") and not os.path.isfile(name + ".nhr") and not os.path.isfile(name + ".nsq"):

        if not isProt:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log")
        else:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log")

    elif overwrite:
        if not isProt:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log")
        else:
            os.system(
                "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log")

    else:
        print("BLAST DB files found. Using existing DBs..")

    return name


def Create_Blastdb_no_fasta(questionDB, overwrite, dbtypeProt, sequence):

    base = os.path.basename(questionDB)
    dirname = os.path.dirname(questionDB)
    isProt = dbtypeProt

    if len(dirname) == 0:
        dirname = '.'
    basename = os.path.splitext(base)[0]
    ensure_dir(dirname + "/blastdbs")
    name = dirname + "/blastdbs/" + basename + "_db"

    if not os.path.isfile(name + ".nin") and not os.path.isfile(name + ".nhr") and not os.path.isfile(name + ".nsq"):

        if not isProt:
            os.system('echo "'+sequence+'" | makeblastdb -in - -title titulo -out ' + name + ' -dbtype nucl -logfile ' + name + '_blast.log')
        else:
            os.system('echo "'+sequence+'" | makeblastdb -in - -title titulo -out ' + name + ' -dbtype prot -logfile ' + name + '_blast.log')


    elif overwrite:
        if not isProt:
            os.system('echo "'+sequence+'" | makeblastdb -in - -title titulo -out ' + name + ' -dbtype nucl -logfile ' + name + '_blast.log')
        else:
            os.system('echo "'+sequence+'" | makeblastdb -in - -title titulo -out ' + name + ' -dbtype prot -logfile ' + name + '_blast.log')

    else:
        print("BLAST DB files found. Using existing DBs..")

    return name


def runBlast(cline, bOutFile, locus_sbjct):
    os.system(str(cline))
    rec = open(bOutFile)
    blast_record = NCBIXML.read(rec)

    if os.path.isfile(locus_sbjct):
        os.remove(locus_sbjct)
    os.remove(bOutFile)

    return blast_record


def runBlastParser(cline, bOutFile):
    os.system(str(cline))
    rec = open(bOutFile)
    blast_records = NCBIXML.parse(rec)

    #	if os.path.isfile(locus_sbjct):
    #		os.remove(locus_sbjct)

    # os.remove(bOutFile)

    return blast_records
