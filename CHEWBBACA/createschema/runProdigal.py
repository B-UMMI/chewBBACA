#!/usr/bin/env python3
import os
import pickle
import subprocess


def main(input_file,tempPath,choosenTaxon,translation_table):

    contigsFasta = input_file

    basepath = tempPath

    # ------------ #
    # RUN PRODIGAL #
    # ------------ #
    # prodigal_path='prodigal'

    if choosenTaxon == "False":

        proc = subprocess.Popen(
            ['prodigal', '-i', contigsFasta, '-c', '-m', '-g', str(translation_table), '-p', 'single', '-f', 'sco', '-q'],
            stdout=subprocess.PIPE)
    else:
        proc = subprocess.Popen(
            ['prodigal', '-i', contigsFasta, '-c', '-m', '-g', str(translation_table), '-p', 'single', '-f', 'sco', '-q', '-t',
             choosenTaxon], stdout=subprocess.PIPE)

    cdsDict = {}
    tempList = []
    line = ' '
    while line != '':

        # when it finds a contig tag
        if "seqhdr" in line:
            # add contig to cdsDict and start new entry

            if len(tempList) > 0:

                # --- brute force parsing of the contig tag - better solution is advisable --- #

                i = 0
                for l in contigTag:
                    if l == ' ':
                        break
                    i += 1
                contigTag = contigTag[:i]

                cdsDict[contigTag.replace("\r", "")] = tempList
                tempList = []

            contigTag = line.split('"')[-2]

        # when it finds a line with cds indexes
        elif line[0] == '>':

            # parsing
            cdsL = line.split('_')

            # --- each element of this list is a pair of indices - the start and the end of a CDS --- #

            tempList.append([int(cdsL[1]) - 1, int(
                cdsL[2])])  # start index correction needed because prodigal indexes start in 1 instead of 0

        # reads the stdout from 'prodigal'
        line = proc.stdout.readline().decode("utf-8")

    # ADD LAST
    if len(tempList) > 0:

        # --- brute force parsing of the contig tag - better solution is advisable --- #

        i = 0
        for l in contigTag:
            if l == ' ':
                break
            i += 1
        contigTag = contigTag[:i]

        cdsDict[contigTag.replace("\r", "")] = tempList

    filepath = os.path.join(basepath, str(os.path.basename(contigsFasta)) + "_ORF.txt")
    with open(filepath, 'wb') as f:
        var = cdsDict
        pickle.dump(var, f)

    print("done prodigal run on:" + str(os.path.basename(contigsFasta)))

    return True


if __name__ == "__main__":
    main()
