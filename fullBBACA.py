#!/usr/bin/env python
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import sys
import os
import argparse
import time
import multiprocessing
import subprocess


def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return "Not found"


def call_proc(cmd):
    p = subprocess.Popen(args=cmd)
    out, err = p.communicate()
    return p


def check_if_list_or_folder(folder_or_list):
    list_files = []
    # check if given a list of genomes paths or a folder to create schema
    try:
        f = open(folder_or_list, 'r')
        f.close()
        list_files = folder_or_list
    except IOError:

        for gene in os.listdir(folder_or_list):
            try:
                genepath = os.path.join(folder_or_list, gene)
                for allele in SeqIO.parse(genepath, "fasta", generic_dna):
                    break
                list_files.append(os.path.abspath(genepath))
            except Exception as e:
                print e
                pass

    return list_files


# ================================================ MAIN ================================================ #

def main():
    parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes provided a schema")
    genes_group = parser.add_mutually_exclusive_group(required=True)
    genes_group.add_argument("--cs", nargs='?', type=str, help='folder of genomes to use for schema creation',
                             default=False)
    genes_group.add_argument("--genes", nargs='?', type=str, help='list of target genes fasta file', default=False)
    parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False, default=0.6)
    parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2",
                        required=True)
    parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False, default='blastp')
    parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False, default=False)
    parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
    parser.add_argument('--genomes', nargs='?', type=str, help='List of genome files (list of fasta files)',
                        required=True)
    parser.add_argument("--auto", help="full automatic pipeline", required=False, action="store_true", default=False)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", dest='verbose', action="store_true",
                        default=False)
    args = parser.parse_args()

    genomes2CreateSchema = args.cs
    genes2call = args.genes
    BSRTresh = args.bsr
    cpuToUse = args.cpu
    BlastpPath = args.b
    chosenTaxon = args.t
    gOutFile = os.path.abspath(args.o)
    genomes2call = args.genomes
    auto = args.auto
    verbose = args.verbose

    # avoid user to run the script with all cores available, could impossibilitate any usage when running on a laptop
    if cpuToUse > multiprocessing.cpu_count() - 2:
        print "Warning, you are close to use all your cpus, if you are using a laptop you may be uncapable to perform any action"

    taxonList = {'Campylobacter_Jejuni': 'trained_campyJejuni.trn',
                 'Acinetobacter_Baumannii': 'trained_acinetoBaumannii.trn',
                 'Streptococcus_Agalactiae': 'trained_strepAgalactiae.trn',
                 'Haemophilus_Influenzae': 'trained_haemoInfluenzae_A.trn',
                 'Yersinia_Enterocolitica': 'trained_yersiniaEnterocolitica.trn',
                 'Escherichia_Coli': 'trained_eColi.trn',
                 'Enterococcus_Faecium': 'trained_enteroFaecium.trn',
                 'Staphylococcus_Haemolyticus': 'trained_staphHaemolyticus.trn',
                 'Salmonella_Enterica_enteritidis': 'trained_salmonellaEnterica_enteritidis.trn'
                 }
    if isinstance(chosenTaxon, basestring):
        trainingFolderPAth = os.path.abspath(os.path.join(os.path.dirname(__file__), 'TrainingFiles4Prodigal'))
        try:
            os.path.join(trainingFolderPAth, taxonList[chosenTaxon])

        except:
            print "Your chosen taxon is not attributed, select one from:"
            for elem in taxonList.keys():
                print elem
            return "retry"

    print BlastpPath

    print ("Will use this number of cpus: " + str(cpuToUse))
    print ("Checking all programs are installed")

    print ("Checking Blast installed... " + str(which(str(BlastpPath))))
    print ("Checking Prodigal installed... " + str(which('prodigal')))

    # check version of Blast

    proc = subprocess.Popen([BlastpPath, '-version'], stdout=subprocess.PIPE)
    line = proc.stdout.readline()
    if not "blastp: 2.5." in str(line) and not "blastp: 2.6." in str(line):
        print "your blast version is " + str(line)
        print "update your blast to 2.5.0 or above, will exit program"
        sys.exit()
    else:
        print "blast version is up to date, the program will continue"

    main_starttime = "\nStarting full chewBBACA at : " + time.strftime("%H:%M:%S-%d/%m/%Y")
    print (main_starttime)

    if genes2call:

        # check if given a list of genes paths or a folder with genes fastas
        genes2call = check_if_list_or_folder(genes2call)

        if isinstance(genes2call, list):
            with open("listGenes2Call.txt", "wb") as f:
                for genome in genes2call:
                    f.write(genome + "\n")
            genes2call = "listGenes2Call.txt"

    # if user want to create schema
    str_genomes_4_sc = ''
    if genomes2CreateSchema:

        schemaSeedDir = os.path.normpath("schema_seed")

        # check if given a list of genomes paths or a folder to create schema

        genomes2CreateSchema = check_if_list_or_folder(genomes2CreateSchema)

        if isinstance(genomes2CreateSchema, list):
            with open("listGenomes4Schema.txt", "wb") as f:
                for genome in genomes2CreateSchema:
                    f.write(genome + "\n")
                    str_genomes_4_sc += genome + "\n"
            genomes2CreateSchema = "listGenomes4Schema.txt"

        # create the schema
        print("Starting the schema creation script")

        ppanScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'createschema/PPanGen.py')

        args = [ppanScriptPath, '-i', genomes2CreateSchema, '--cpu', str(cpuToUse), "-t", chosenTaxon, "-o",
                schemaSeedDir,
                "--bsr", str(BSRTresh), "-b", "blastp"]

        if verbose:
            args.append('-v')

        proc = subprocess.Popen(args)
        proc.wait()

        # create list with genes to call and call bbaca
        genes2call = "listGenes2Call.txt"

        with open(genes2call, "wb") as f:
            for gene in os.listdir(schemaSeedDir):
                try:
                    genepath = os.path.join(schemaSeedDir, gene)
                    #gene_fp2 = HTSeq.FastaReader(genepath)
                    for allele in SeqIO.parse(genepath, "fasta", generic_dna):
                        break
                    f.write(os.path.abspath(genepath) + "\n")
                except Exception as e:
                    print e
                    pass

    # else user provided the genes, check if its a folder or a list of files
    else:

        genes2call = check_if_list_or_folder(genes2call)
        if isinstance(genes2call, list):
            with open("listGenes2Call.txt", "wb") as f:
                for genome in genes2call:
                    f.write(genome + "\n")
            genes2call = "listGenes2Call.txt"

        # check if short folder exists on the genes folder, if not, create it

        with open(genes2call, 'r') as f:
            first_gene = f.readline()

        short_folder = os.path.join(os.path.dirname(os.path.abspath(first_gene)), "short")
        if not os.path.isdir(short_folder):
            ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils/init_schema_4_bbaca.py')
            proc = subprocess.Popen([ScriptPath, '-i', genes2call])
            proc.wait()

    # create the list of genomes to use for the call, including the use on the schema creation

    genomes2call = check_if_list_or_folder(genomes2call)

    if isinstance(genomes2call, list):
        with open("listGenomes2Call.txt", "wb") as f:
            for genome in genomes2call:
                f.write(genome + "\n")
        genomes2call = "listGenomes2Call.txt"

    with open("listGenomes2Call.txt") as f:
        str_genomes_4_sc += f.read()

    with open("listGenomes2Call.txt", "wb") as f:
        f.write(str_genomes_4_sc)

    # run allele call on the genomes
    print("Starting the allele call script")
    ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'allelecall/BBACA.py')
    args = [ScriptPath, '-i', genomes2call, '-g', genes2call, '--cpu', str(cpuToUse), "-t",
            chosenTaxon, "-o", "results", "--bsr", str(BSRTresh)]
    if auto:
        args.append('--fc')
    if verbose:
        args.append('-v')

    proc = subprocess.Popen(args)
    proc.wait()

    # select the last results folder created
    last_result_folder = max([os.path.join("results", d) for d in os.listdir("results")], key=os.path.getmtime)

    last_result_folder = os.path.abspath(last_result_folder)

    # detect and remove paralogs

    print("Detecting the paralog loci and removing them")

    if not os.path.exists(gOutFile):
        os.makedirs(gOutFile)

    ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils/ParalogPrunning.py')
    proc = subprocess.Popen(
        [ScriptPath, '-i', os.path.join(last_result_folder, "results_contigsInfo.tsv"), '-o', gOutFile])
    proc.wait()

    ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils/RemoveGenes.py')
    proc = subprocess.Popen([ScriptPath, '-i', os.path.join(last_result_folder, "results_alleles.tsv"), '-o',
                             os.path.join(gOutFile, "results_alleles_no_paralog"), '-g',
                             os.path.join(gOutFile, "RepeatedLoci.txt")])
    proc.wait()

    # testqualitygenomes, show plot and ask for threshold to use for phyloviz, set default to 20
    print("Starting the test quality genomes script")
    max_thresh = 200
    step = 5
    max_iteration = 13
    ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils/TestGenomeQuality.py')
    proc = subprocess.Popen(
        [ScriptPath, '-i', os.path.join(gOutFile, "results_alleles_no_paralog.tsv"), '-n', str(max_iteration), '-t',
         str(max_thresh), '-s', str(step), '-o', gOutFile])
    proc.wait()

    # get the genomes to remove at specific threshold, 10 if number of genomes is low
    num_lines = sum(1 for line in open(genomes2call))
    if num_lines < 50:
        thresh = 10
    else:
        thresh = 20

    if not auto:
        list_possible_thresh = range(0, max_thresh, step)
        testVar = raw_input(
            "(wait for a browser page to auto-open a plot) which threshold you want to use? (value on X axis) : ")

        try:
            thresh = str(int(testVar))
            if int(testVar) not in list_possible_thresh:
                print ("seems the threshold you selected was not found, the program will use " + str(thresh))
        except:
            print "the value you use is incorrect : " + testVar + " , the program will use " + str(thresh)

    list_genomes2remove = ''
    file_genomres2remove = os.path.abspath("listGenomes2Del.txt")
    for line in open(os.path.join(gOutFile, "removedGenomes.txt")):
        token = "using a threshold of " + str(thresh) + " at"
        if line.startswith(token):
            list_genomes2remove = line.split('\t')
            del list_genomes2remove[0]
            break

    with open(file_genomres2remove, "wb") as f:
        for genome2del in list_genomes2remove:
            f.write(genome2del + "\n")

    # extract the cgMLST using the genomes removed at specific threshold
    print("Extracting the cgMLST after removing the genomes you selected on the previous step")
    ScriptPath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils/Extract_cgAlleles.py')
    proc = subprocess.Popen([ScriptPath, '-i', os.path.join(gOutFile, "results_alleles_no_paralog.tsv"),
                             '-o', gOutFile, '-g', 'listGenomes2Del.txt'])
    proc.wait()

    # send data to phyloviz online ????




    # run shema evaluator

    # remove the list files
    os.remove(genomes2call)
    os.remove(genes2call)
    os.remove(file_genomres2remove)

    print (main_starttime)
    print ("Finished full chewBBACA at : " + time.strftime("%H:%M:%S-%d/%m/%Y"))


if __name__ == "__main__":
    main()
