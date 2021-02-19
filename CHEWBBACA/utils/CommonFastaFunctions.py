#!/usr/bin/env python3


import os

from Bio.Blast import NCBIXML


def ensure_dir(f):
    if not os.path.isdir(f):
        # print "No BLAST db dir was found. Creating it..."
        os.makedirs(f)


# overwrite has to be 1 or 0 (True or False)
def Create_Blastdb(makeblastdb_path, questionDB, overwrite, dbtypeProt):
    base = os.path.basename(questionDB)
    dirname = os.path.dirname(questionDB)
    isProt = dbtypeProt

    if len(dirname) == 0:
        dirname = '.'
    basename = os.path.splitext(base)[0]
    ensure_dir(dirname + "/blastdbs")
    name = dirname + "/blastdbs/" + basename + "_db"

    cmd_template = '{0} -in {1} -out {2} -dbtype {3} -logfile {2}_blast.log'

    if not os.path.isfile(name + ".nin") and not os.path.isfile(name + ".nhr") and not os.path.isfile(name + ".nsq"):

        if not isProt:
            os.system(cmd_template.format(makeblastdb_path, questionDB, name, 'nucl'))
        else:
            os.system(cmd_template.format(makeblastdb_path, questionDB, name, 'prot'))

    elif overwrite:
        if not isProt:
            os.system(cmd_template.format(makeblastdb_path, questionDB, name, 'nucl'))
        else:
            os.system(cmd_template.format(makeblastdb_path, questionDB, name, 'prot'))

    else:
        print("BLAST DB files found. Using existing DBs..")

    return name


def Create_Blastdb_no_fasta(makeblastdb_path, questionDB, overwrite, dbtypeProt, sequence):

    base = os.path.basename(questionDB)
    dirname = os.path.dirname(questionDB)
    isProt = dbtypeProt

    if len(dirname) == 0:
        dirname = '.'
    basename = os.path.splitext(base)[0]
    ensure_dir(dirname + "/blastdbs")
    name = dirname + "/blastdbs/" + basename + "_db"

    cmd_template = 'echo "{0}" | {1} -in - -title titulo -out {2} -dbtype {3} -logfile {2}_blast.log'

    if not os.path.isfile(name + ".nin") and not os.path.isfile(name + ".nhr") and not os.path.isfile(name + ".nsq"):

        if not isProt:
            os.system(cmd_template.format(sequence, makeblastdb_path, name, 'nucl'))
        else:
            os.system(cmd_template.format(sequence, makeblastdb_path, name, 'prot'))

    elif overwrite:
        if not isProt:
            os.system(cmd_template.format(sequence, makeblastdb_path, name, 'nucl'))
        else:
            os.system(cmd_template.format(sequence, makeblastdb_path, name, 'prot'))

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

    return blast_records
