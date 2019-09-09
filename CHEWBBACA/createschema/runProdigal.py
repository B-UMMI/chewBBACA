#!/usr/bin/env python3
import sys
import os
import subprocess
import pickle

def main(input_file,tempPath,choosenTaxon):

    contigsFasta = input_file
    basepath = tempPath


    # ------------ #
    # RUN PRODIGAL #
    # ------------ #


    if choosenTaxon == "False":
    
        proc = subprocess.Popen(
            ['prodigal', '-i', contigsFasta, '-c', '-m', '-g', '11', '-p', 'single', '-f', 'sco', '-q'],
            stdout=subprocess.PIPE)
    else:
        proc = subprocess.Popen(
            ['prodigal', '-i', contigsFasta, '-c', '-m', '-g', '11', '-p', 'single', '-f', 'sco', '-q', '-t',
                choosenTaxon], stdout=subprocess.PIPE)
    
    cdsDict = {}
    tempList = []
    # Reads the stdout from Prodigal
    prodigal_out = proc.stdout.readlines()

    # Parse 'Prodigal's output 
    for line in prodigal_out:
        line_decoded = line.decode("utf-8")
        
        if "seqhdr" in line_decoded:
            seqid = line_decoded.split('"')[1].split()[0]

        # Obtain the start and end positions of the CDSs
        elif ">" in line_decoded:
            cdsL = line_decoded.split("_")
            
            # Start index correction needed because Prodigal indexes start in 1 instead of 0
            start_position = int(cdsL[1]) - 1
            end_position = int(cdsL[2])
            
            # Strand of the CDS. 1 if sense (+), else 0 (-) (antisense)
            strand = 1 if cdsL[-1].strip() == "+" else 0
            tempList = [start_position, end_position, strand]
        
            if len(tempList) > 0:
                if seqid not in cdsDict:
                    # Add the sequence ID as the key and the list of lists of the CDSs' positions as the value
                    cdsDict[seqid] = [tempList]
                else:
                    cdsDict[seqid] += [tempList]
                    tempList = []
    
    # Add the cdsDict to a file
    file_basename = os.path.basename(contigsFasta).split('.')[0]
    filepath = os.path.join(basepath, file_basename + "_ORF.txt")
    with open(filepath, 'wb') as f:
        var = cdsDict
        pickle.dump(var, f)
    
    print("done prodigal run on: " + str(os.path.basename(contigsFasta)))

    return True


if __name__ == "__main__":
    main()
