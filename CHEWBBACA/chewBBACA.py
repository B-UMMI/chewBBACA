#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHORS

    Mickael Silva
    github: @

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import sys
import platform
import argparse

from allelecall import BBACA
from createschema import PPanGen
from SchemaEvaluator import ValidateSchema
from PrepExternalSchema import PrepExternalSchema
from utils import (TestGenomeQuality, profile_joiner,
                   uniprot_find, Extract_cgAlleles,
                   RemoveGenes, auxiliary_functions as aux)

from CHEWBBACA_NS import (down_schema, load_schema,
                          sync_schema, down_profiles,
                          send2NS, send_metadata,
                          stats_requests)

import CHEWBBACA


current_version = '2.1.0'


def create_schema():

    def msg(name=None):
        return '''chewBBACA.py CreateSchema [CreateSchema ...] [-h]
                 -i [I] -o [O] --cpu [CPU] [-b [B]] [--bsr [BSR]]
                 [--ptf [PTF]] [-v] [-l [L]]'''

    parser = argparse.ArgumentParser(description='Creates a wgMLST '
                                                 'schema based on a '
                                                 'set of input genomes.',
                                     usage=msg())

    parser.add_argument('CreateSchema', nargs='+',
                        help='')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory with the genomes FASTA '
                             'files or to a file with a list of paths to '
                             'the FASTA files, one per line.')

    parser.add_argument('-o', nargs='?', type=str, required=True,
                        dest='schema_dir',
                        help='Output folder where the schema will be created.')

    parser.add_argument('--cpu', nargs='?', type=int, required=True,
                        dest='cpu_num',
                        help='Number of CPU cores/threads that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores/threads).')

    parser.add_argument('-b', nargs='?', type=str, required=False,
                        default='blastp', dest='blastp_path',
                        help='Path to the BLASTp executables.')

    parser.add_argument('--CDS', required=False, action='store_true',
                        default=False, dest='cds_input',
                        help='Input is a FASTA file with one representative '
                             'sequence per gene in the schema.')

    parser.add_argument('--bsr', nargs='?', type=float, required=False,
                        default=0.6, dest='bsr_value',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--ptf', nargs='?', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('-v', '--verbose', action='store_true',
                        default=False, dest='verbose',
                        help='Increased output verbosity during execution.')

    parser.add_argument('-l', nargs='?', type=int, required=False,
                        default=201, dest='min_seq_len',
                        help='Minimum sequence length accepted for a '
                             'coding sequence to be included in the schema.')

    args = parser.parse_args()

    input_files = args.input_files
    schema_dir = args.schema_dir
    cpu_num = args.cpu_num
    blastp_path = args.blastp_path
    cds_input = args.cds_input
    bsr = args.bsr_value
    ptf_path = args.ptf_path
    verbose = args.verbose
    min_length = args.min_seq_len

    if cds_input is True:
        input_files = [os.path.abspath(input_files)]
    else:
        input_files = aux.check_input_type(input_files, 'listGenomes2Call.txt')

    # determine if user provided training file
    # if not, try to find species training file included in chewBBACA
    ptf_path = aux.check_ptf(ptf_path, CHEWBBACA.__file__)

    if ptf_path is False:
        print('Could not find a valid Prodigal training file.')
        return 1

    # start CreateSchema process
    PPanGen.main(input_files, cpu_num, schema_dir,
                 bsr, blastp_path, min_length,
                 verbose, ptf_path, cds_input)

    # remove temporary file with paths
    # to genome files
    if os.path.isfile(input_files):
        os.remove(input_files)


def allele_call():

    def msg(name=None):
        return '''chewBBACA.py AlleleCall [AlleleCall ...][-h]
                 -i [I] -g [G] -o [O] --cpu [CPU] [-v] [-b [B]]
                 [--bsr [BSR]] [--st] [--ptf [PTF]] [--fc] '
                 [--fr] [--json]'''

    parser = argparse.ArgumentParser(description='Performs allele calling '
                                                 'on a set of input genomes.',
                                     usage=msg())

    parser.add_argument('AlleleCall', nargs='+', help='')
    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory with the genomes FASTA '
                             'files or to a file with a list of paths to '
                             'the FASTA files, one per line.')
    parser.add_argument('-g', nargs='?', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the schema directory with the'
                             ' genes FASTA files.')
    parser.add_argument('-o', nargs='?', type=str, required=True,
                        dest='output_dir',
                        help='Output directory where the allele '
                             'calling results will be stored.')
    parser.add_argument('--cpu', nargs='?', type=int, required=True,
                        dest='cpu_num',
                        help='Number of CPU cores/threads that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores/threads).')
    parser.add_argument('--contained', action='store_true', required=False,
                        default=False, dest='contained',
                        help=argparse.SUPPRESS)
    parser.add_argument('--CDS', action='store_true', required=False,
                        default=False, dest='cds_input',
                        help=argparse.SUPPRESS)
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        dest='verbose',
                        help='Increased output verbosity during execution.')
    parser.add_argument('-b', nargs='?', type=str, required=False,
                        default='blastp', dest='blastp_path',
                        help='Path to the BLASTp executables.')
    parser.add_argument('--bsr', nargs='?', type=float, required=False,
                        default=0.6, dest='bsr_value',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')
    parser.add_argument('--st', nargs='?', type=float, required=False,
                        default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent will be classified as ASM/ALM')
    parser.add_argument('--ptf', nargs='?', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')
    parser.add_argument('--fc', action='store_true', required=False,
                        default=False, dest='force_continue',
                        help='Continue allele call process that '
                             'was interrupted.')
    parser.add_argument('--fr', action='store_true', required=False,
                        default=False, dest='force_reset',
                        help='Force process reset even if there '
                             'are temporary files from a previous '
                             'process that was interrupted.')
    parser.add_argument('--json', action='store_true', required=False,
                        default=False, dest='json_report',
                        help='Output report in JSON format.')

    args = parser.parse_args()

    input_files = args.input_files
    schema_dir = args.schema_dir
    cpu_num = args.cpu_num
    bsr = args.bsr_value
    size_threshold = args.size_threshold
    verbose = args.verbose
    blastp_path = args.blastp_path
    output_dir = args.output_dir
    ptf_path = args.ptf_path
    force_continue = args.force_continue
    force_reset = args.force_reset
    contained = args.contained
    cds_input = args.cds_input
    json_report = args.json_report
    chosen_taxon = False

    schema_genes = aux.check_input_type(schema_dir, 'listGenes2Call.txt')

    ptf_path = aux.check_ptf(ptf_path, CHEWBBACA.__file__)

    if ptf_path is False:
        print('Could not find a valid Prodigal training file.')
        return 1

    # if is a fasta pass as a list of genomes with a single genome,
    # if not check if is a folder or a txt with a list of paths
    genomes_files = aux.check_input_type(input_files, 'listGenomes2Call.txt')

    BBACA.main(genomes_files, schema_genes, cpu_num,
               output_dir, bsr, blastp_path,
               force_continue, json_report, verbose,
               force_reset, contained, chosen_taxon,
               ptf_path, cds_input, size_threshold)

    # check parameters values in config file and alter if needed
    config_file = os.path.join(schema_dir, '.schema_config')
    with open(config_file, 'rb') as pf:
        params = pickle.load(pf)

    # check bsr
    if bsr not in params['bsr']:
        params['bsr'].append(bsr)
    # check training file
    if ptf_path not in params['prodigal_training_file']:
        params['prodigal_training_file'].append(ptf_path)
    # check chewie version
    if current_version not in params['chewBBACA_version']:
        params['chewBBACA_version'].append(current_version)

    # save updated schema config file
    with open(config_file, 'wb') as pf:
        pickle.dump(params, pf)    

    # remove temporary files with paths to genomes
    # and schema files files
    if os.path.isfile(schema_genes) is True:
        os.remove(schema_genes)
    if os.path.isfile(genomes_files):
        os.remove(genomes_files)


def evaluate_schema():

    def msg(name=None):
        return '''chewBBACA.py SchemaEvaluator [SchemaEvaluator ...] [-h]
                 -i [I] [-p] [--log] -l [L] -ta [TA] [-t [T]]
                 [--title [TITLE]] --cpu [CPU] [-s [S]] [--light]'''

    parser = argparse.ArgumentParser(description='This program analyses a '
                                                 'set of gene files, '
                                                 'analyzing the alleles '
                                                 'CDS and the length of '
                                                 'the alleles per gene',
                                     usage=msg())

    parser.add_argument('SchemaEvaluator', nargs='+',
                        help='evaluation of a schema')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='',
                        help='list genes, directory or .txt file '
                             'with the full path')

    parser.add_argument('-p', action='store_true', required=False,
                        default=False, dest='conserved',
                        help='One bad allele still makes gene conserved.')

    parser.add_argument('--log', action='store_true', default=False,
                        dest='logScale',
                        help='')

    parser.add_argument('-l', nargs='?', type=str, required=True,
                        dest='',
                        help='name/location main html file')

    parser.add_argument('-ta', nargs='?', type=int, required=False,
                        default=11, dest='',
                        help='ncbi translation table')

    parser.add_argument('-t', nargs='?', type=float, required=False,
                        default=0.05, dest='',
                        help='Threshold')

    parser.add_argument('--title', nargs='?', type=str, required=False,
                        default='My Analyzed wg/cg MLST Schema - Rate My Schema',
                        dest='',
                        help='title on the html')

    parser.add_argument('--cpu', nargs='?', type=int, required=True,
                        dest='',
                        help='number of cpu to use')

    parser.add_argument('-s', nargs='?', type=int, required=False,
                        default=500, dest='',
                        help='number of boxplots per page (more than '
                        '500 can make the page very slow).')

    parser.add_argument('--light', action='store_true', required=False,
                        default=False, dest='',
                        help='skip clustal and mafft run')

    args = parser.parse_args()
    genes = args.i
    transTable = args.ta
    logScale = args.logScale
    outputpath = args.l
    cpuToUse = args.cpu
    threshold = args.t
    OneBadGeneNotConserved = args.conserved
    splited = args.s
    light = args.light
    title = str(args.title)

    ValidateSchema.main(genes, cpuToUse, outputpath,
                        transTable, threshold, splited,
                        title, logScale, OneBadGeneNotConserved,
                        light)


def test_schema():

    def msg(name=None):
        return '''chewBBACA.py TestGenomeQuality [TestGenomeQuality ...] [-h]
                 -i [I] -n [N] -t [T] -s [S] [-o [O]] [-v]'''

    parser = argparse.ArgumentParser(description='This program analyzes an '
                                                 'allele call raw output '
                                                 'matrix, returning info on '
                                                 'which genomes are '
                                                 'responsible for cgMLST '
                                                 'loci loss',
                                     usage=msg())

    parser.add_argument('TestGenomeQuality', nargs='+',
                        help='test the quality of the genomes on the '
                             'allele call')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='',
                        help='raw allele call matrix file')

    parser.add_argument('-n', nargs='?', type=int, required=True,
                        dest='',
                        help='maximum number of iterations')

    parser.add_argument('-t', nargs='?', type=int, required=True,
                        dest='',
                        help='maximum threshold of bad calls above 95 percent')

    parser.add_argument('-s', nargs='?', type=int, required=True,
                        dest='',
                        help='step between each threshold analysis')

    parser.add_argument('-o', nargs='?', type=str, required=False,
                        default=".", dest='',
                        help="Folder for the analysis files")

    parser.add_argument("-v", "--verbose", action="store_true", default=False,
                        dest='verbose',
                        help="increase output verbosity")

    args = parser.parse_args()

    pathOutputfile = args.i
    iterationNumber = args.n
    thresholdBadCalls = args.t
    step = args.s
    out_folder = args.o
    verbose = args.verbose

    TestGenomeQuality.main(pathOutputfile, iterationNumber,
                           thresholdBadCalls, step,
                           out_folder, verbose)


def extract_cgmlst():

    def msg(name=None):
        return '''chewBBACA.py ExtractCgMLST [ExtractCgMLST ...] [-h]
                  -i [I] -o [O] [-r [R]] [-g [G]]'''

    parser = argparse.ArgumentParser(description='This program cleans an '
                                                 'output file for phyloviz',
                                     usage=msg())

    parser.add_argument('ExtractCgMLST', nargs='+',
                        help='clean chewBBACA output')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='',
                        help='input file to clean')

    parser.add_argument('-o', nargs='?', type=str, required=True,
                        dest='',
                        help='output folder')

    parser.add_argument('-r', nargs='?', type=str, required=False,
                        default=False, dest='',
                        help='listgenes to remove')

    parser.add_argument('-g', nargs='?', type=str, required=False,
                        default=False, dest='',
                        help='listgenomes to remove')

    parser.add_argument('-p', nargs='?', type=float, required=False,
                        default=1, dest='',
                        help='maximum presence (e.g 0.95)')

    args = parser.parse_args()

    pathOutputfile = args.i
    newfile = args.o
    genes2remove = args.r
    genomes2remove = args.g
    cgMLSTpercent = args.p

    Extract_cgAlleles.main(pathOutputfile, newfile,
                           cgMLSTpercent, genes2remove,
                           genomes2remove)


def remove_genes():

    def msg(name=None):
        return '''chewBBACA.py RemoveGenes [RemoveGenes ...][-h]
                  -i [I] -g [G] -o [O] [--inverse]'''

    parser = argparse.ArgumentParser(description='This program removes genes '
                                                 'from a tab separated allele '
                                                 'profile file',
                                     usage=msg())

    parser.add_argument('RemoveGenes', nargs='+',
                        help='remove loci from a chewBBACA profile')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='',
                        help='main matrix file from which to remove')

    parser.add_argument('-g', nargs='?', type=str, required=True,
                        help='list of genes to remove')

    parser.add_argument('-o', nargs='?', type=str, required=True,
                        help='output file name')

    parser.add_argument("--inverse", action="store_true", default=False,
                        dest='inverse',
                        help="list to remove is actually the one to keep")

    args = parser.parse_args()
    mainListFile = args.i
    toRemoveListFile = args.g
    outputfileName = args.o
    inverse = args.inverse

    RemoveGenes.main(mainListFile, toRemoveListFile, outputfileName, inverse)


def join_profiles():

    def msg(name=None):
        return '''chewBBACA.py JoinProfiles [RemoveGenes ...][-h]
                  -p1 -p2 -o [O]'''

    parser = argparse.ArgumentParser(description='This program joins two '
                                                 'profiles, returning a '
                                                 'single profile file with '
                                                 'the common loci',
                                     usage=msg())

    parser.add_argument('JoinProfiles', nargs='+',
                        help='join profiles')

    parser.add_argument('-p1', nargs='?', type=str, required=True,
                        help='profile 1')

    parser.add_argument('-p2', nargs='?', type=str, required=True,
                        help='profile 2')

    parser.add_argument('-o', nargs='?', type=str, required=True,
                        help='output file name')

    args = parser.parse_args()
    profile1 = args.p1
    profile2 = args.p2
    outputFile = args.o

    profile_joiner.main(profile1, profile2, outputFile)


def prep_schema():

    def msg(name=None):
        return '''chewBBACA.py PrepExternalSchema [PrepExternalSchema ...][-h]
                 -i [I] --cpu [CPU] [-v]'''

    parser = argparse.ArgumentParser(description='This script enables the '
                                                 'adaptation of external '
                                                 'schemas so that the loci '
                                                 'and alleles present in '
                                                 'those schemas can be used '
                                                 'with chewBBACA. During '
                                                 'the process, alleles that '
                                                 'do not correspond to a '
                                                 'complete CDS or that cannot '
                                                 'be translated are discarded '
                                                 'from the final schema. One '
                                                 'or more alleles of each '
                                                 'gene/locus will be chosen '
                                                 'as representatives and '
                                                 'included in the "short" '
                                                 'directory.')

    parser.add_argument('PrepExternalSchema', nargs='+',
                        help='Adapt an external schema to be used with '
                        'chewBBACA.')

    parser.add_argument('-i', type=str, required=True, dest='input_files',
                        help='Path to the folder containing the fasta files, '
                             'one fasta file per gene/locus (alternatively, '
                             'a file with a list of paths can be given).')

    parser.add_argument('-o', type=str, required=True, dest='output_directory',
                        help='The directory where the output files will be '
                        'saved (will create the directory if it does not '
                        'exist).')

    parser.add_argument('--cpu', type=int, required=False, default=1,
                        dest='core_count',
                        help='The number of CPU cores to use (default=1).')

    parser.add_argument('--bsr', type=float, required=False, default=0.6,
                        dest='blast_score_ratio',
                        help='The BLAST Score Ratio value that will be '
                        'used to adapt the external schema (default=0.6).')

    parser.add_argument('--len', type=int, required=False, default=0,
                        dest='minimum_length',
                        help='Minimum sequence length accepted. Sequences with'
                        ' a length value smaller than the value passed to this'
                        ' argument will be discarded (default=0).')

    parser.add_argument('--tbl', type=int, required=False, default=11,
                        dest='translation_table',
                        help='Genetic code to use for CDS translation.'
                        ' (default=11, for Bacteria and Archaea)')

    args = parser.parse_args()

    gene_files = args.input_files
    output_dir = args.output_directory
    processes = args.core_count
    bsr = args.blast_score_ratio
    min_len = args.minimum_length
    trans_table = args.translation_table

    PrepExternalSchema.main(gene_files, output_dir, processes,
                            bsr, min_len, trans_table)


def find_uniprot():

    def msg(name=None):
        return '''chewBBACA.py UniprotFinder [UniprotFinder ...][-h]
                  -i [I] -t [T] --cpu [CPU]'''

    parser = argparse.ArgumentParser(description='This program gets '
                                                 'information of each '
                                                 'locus created on the '
                                                 'schema creation, based '
                                                 'on the uniprot database')

    parser.add_argument('UniprotFinder', nargs='+',
                        help='get info about a schema created with chewBBACA')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='',
                        help='path to folder containg the schema fasta '
                             'files ( alternative a list of fasta files)')

    parser.add_argument('-t', nargs='?', type=str, required=True,
                        dest='',
                        help='path to proteinID_Genome.tsv file generated')

    parser.add_argument('--cpu', nargs='?', type=int, required=False,
                        default=1, dest='',
                        help='number of cpu')

    args = parser.parse_args()

    geneFiles = args.i
    tsvFile = args.t
    cpu2use = args.cpu

    uniprot_find.main(geneFiles, tsvFile, cpu2use)


def download_schema_NS():

    def msg(name=None):
        return '''chewBBACA.py Download_NS_Schema [Download_NS_Schema ...][-h]
                 -ns_schema [NS_SCHEMA] -ns_species [NS_SPECIES]
                 -down_folder [DOWN_FOLDER] --cpu [CPU] --ns_url [NS_URL]'''

    parser = argparse.ArgumentParser(description='This program downloads '
                                                 'a schema from the NS.',
                                     usage=msg())

    parser.add_argument('Download_NS_Schema', nargs='+',
                        help='This program downloads a schema from '
                             'the NS.')

    parser.add_argument('-ns_schema', nargs='?', type=str, required=True,
                        dest='ns_schema',
                        help='The URI, integer identifier or description of '
                             'the schema to download from the NS.')

    parser.add_argument('-ns_species', nargs='?', type=str, required=True,
                        dest='ns_species',
                        help='The integer identifier or name of the species '
                             'that the schema is associated to in the NS.')

    parser.add_argument('-down_folder', type=str, required=True,
                        dest='down_folder',
                        help='Output folder to which the schema will '
                             'be saved.')

    parser.add_argument('--cpu', nargs='?', type=int, required=False,
                        default=1, dest='cpu_num',
                        help='Number of CPU cores/threads that will '
                             'be passed to the PrepExternalSchema to '
                             'determine representatives and create the '
                             'final schema.')

    parser.add_argument('--ns_url', type=str, required=False,
                        default='http://127.0.0.1:5000/NS/api/',
                        dest='ns_url',
                        help='The base URL for the Nomenclature Server.')

    args = parser.parse_args()

    ns_schema = args.ns_schema
    ns_species = args.ns_species
    down_folder = args.down_folder
    cpu_num = args.cpu_num
    ns_url = args.ns_url

    down_schema.main(ns_schema, ns_species, down_folder,
                     cpu_num, ns_url)


def load_schema_NS():

    def msg(name=None):
        return '''chewBBACA.py Load_NS_Schema [Load_NS_Schema ...][-h]
                 -schema [SCHEMA] -ns_species [NS_SPECIES]
                 -schema_desc [SCHEMA_DESC] -loci_prefix [LOCI_PREFIX]
                 --thr [THREADS] --ns_url [NS_URL]
                 --continue_up [CONTINUE_UP]'''

    parser = argparse.ArgumentParser(description='This program uploads '
                                                 'a schema to the NS.',
                                     usage=msg())

    parser.add_argument('Load_NS_Schema', nargs='+',
                        help='This program loads a schema to '
                             'the NS.')

    parser.add_argument('-schema', nargs='?', type=str, required=True,
                        dest='schema',
                        help='Path to the directory with the local schema '
                             'files.')

    parser.add_argument('-ns_species', nargs='?', type=str, required=True,
                        dest='ns_species',
                        help='The integer identifier or name of the species '
                             'that the schema will be associated to in '
                             'the NS.')

    parser.add_argument('-schema_desc', type=str, required=True,
                        dest='schema_desc',
                        help='A brief and meaningful description that '
                             'should help understand the type and content '
                             'of the schema.')

    parser.add_argument('-loci_prefix', nargs='?', type=str, required=True,
                        dest='loci_prefix',
                        help='Prefix included in the name of each locus of '
                             'the schema.')

    parser.add_argument('--thr', nargs='?', type=int, required=False,
                        default=20, dest='threads',
                        help='Number of threads to use to upload the alleles '
                             'of the schema.')

    parser.add_argument('--ns_url', type=str, required=False,
                        default='http://127.0.0.1:5000/NS/api/',
                        dest='ns_url',
                        help='The base URL for the Nomenclature Server.')

    parser.add_argument('--continue_up', required=False, default='no',
                        dest='continue_up', choices=['no', 'yes'],
                        help='If the process should check if the schema '
                             'upload was interrupted and try to finish it.')

    args = parser.parse_args()

    schema = args.schema
    ns_species = args.ns_species
    schema_desc = args.schema_desc
    loci_prefix = args.loci_prefix
    thr = args.threads
    ns_url = args.ns_url
    continue_up = args.continue_up

    load_schema.main(schema, ns_species, schema_desc,
                     loci_prefix, thr, ns_url,
                     continue_up)


def sync_schema_NS():

    def msg(name=None):
        return '''chewBBACA.py Sync_NS_Schema [Sync_NS_Schema ...][-h]
                  -schema_dir [SCHEMA_DIR] --cpu [CPU] --ns_url [NS_URL]
                  --submit [SUBMIT]'''

    parser = argparse.ArgumentParser(description='This program syncs a local '
                                                 'schema with NS',
                                     usage=msg())

    parser.add_argument('Sync_NS_Schema', nargs='+',
                        help='Synchronize a local schema, previously '
                             'downloaded from the NS, with its latest '
                             'version in the NS.')

    parser.add_argument('-schema_dir', nargs='?', type=str, required=True,
                        dest='schema_dir',
                        help='Path to the directory with the local schema '
                             'files.')

    parser.add_argument('--cpu', nargs='?', type=int, required=False,
                        default=1, dest='cpu_num',
                        help='Number of CPU cores/threads that will '
                             'be passed to the PrepExternalSchema to '
                             'determine representatives for the '
                             'updated version of the schema.')

    parser.add_argument('--ns_url', nargs='?', type=str, required=False,
                        default='http://127.0.0.1:5000/NS/api/',
                        dest='ns_url',
                        help='The base URL for the Nomenclature Server.')

    parser.add_argument('--submit', nargs='?', type=str, required=False,
                        default='no', dest='submit',
                        help='If the local alleles that are not in the NS '
                             'should be uploaded to update the NS schema.')

    args = parser.parse_args()

    schema_dir = args.schema_dir
    cpu_num = args.cpu_num
    ns_url = args.ns_url
    submit = args.submit

    sync_schema.main(schema_dir, cpu_num,
                     ns_url, submit)


def send_NS():

    def msg(name=None):
        return ''' chewBBACA.py Send2NS [Send2NS ...][-h] -s [S] -t [T] -p [P]
                    '''

    parser = argparse.ArgumentParser(description="Send local profile and respective alleles to NS",usage=msg())
    parser.add_argument('Send2NS', nargs='+', help='send profiles and local alleles to NS')
    parser.add_argument('-s', nargs='?', type=str, help='path to schema folder', required=True)
    parser.add_argument('-p', nargs='?', type=str, help='tsv with profile', required=True)
    parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)
    parser.add_argument('-m', nargs='?', type=str, help='tsv with metadata', required=False, default=False)
    parser.add_argument('--mdr', nargs='?', type=str, help='maximum missing data allowed to fail, default 0.5 (50 percent missing data allowed). 1 == all profiles are uploaded even with 100 percent missing data', required=False, default=0.5)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

    args = parser.parse_args()

    profileFile = args.p
    pathSchema = args.s
    token= args.t
    metadata = args.m
    cpu2use = args.cpu
    percentMDallowed=args.mdr

    send2NS.main(profileFile,pathSchema,token,metadata,percentMDallowed,cpu2use)


def send_meta():

    def msg(name=None):
        return ''' chewBBACA.py SendMetadata [SendMetadata ...][-h] -s [S] -t [T] -p [P]
                    '''

    parser = argparse.ArgumentParser(description="send metadata to isolates on the NS",usage=msg())
    parser.add_argument('SendMetadata', nargs='+', help='send metadata to isolates on the NS')
    parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)
    parser.add_argument('-m', nargs='?', type=str, help='tsv with metadata', required=False, default=False)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)

    args = parser.parse_args()

    token= args.t
    metadata = args.m
    cpu2use = args.cpu

    send_metadata.main(metadata,cpu2use,token)


def down_prof():

    def msg(name=None):
        return ''' chewBBACA.py DownloadProfiles [DownloadProfiles ...][-h] --sp [SP] --sc [SC] --cpu [CPU]
                    '''

    parser = argparse.ArgumentParser(
        description="Download profiles from the NS",usage=msg())
    parser.add_argument('DownloadProfiles', nargs='+', help='download profiles from NS')
    parser.add_argument('--sp', nargs='?', type=str, help='species uri', required=True)
    parser.add_argument('--sc', nargs='?', type=str, help='schema id', required=True)
    parser.add_argument('--cpu', nargs='?', type=int, help='number of cpu', required=False, default=1)
    parser.add_argument('-r', nargs='?', type=str, help='genomes to down profile', required=False, default=None)
    parser.add_argument('-p', nargs='?', type=str, help='profile with already downloaded profiles for that schema', required=False, default=None)
    parser.add_argument('-t', nargs='?', type=str, help='private token', required=False, default=False)

    args = parser.parse_args()

    species = args.sp
    schema = args.sc
    cpu2use = args.cpu
    genomes2Down = args.r
    inputProfile = args.p
    token = args.t

    down_profiles.main(species,schema,cpu2use,inputProfile,genomes2Down,token)


def ns_stats():

    def msg(name=None):
        return '''chewBBACA.py NSStats [NSStats ...][-h] -m [MODE]
                  --url [URL] --taxon [TAXON]'''

    parser = argparse.ArgumentParser(description='',
                                     usage=msg())

    parser.add_argument('NSStats', nargs='+', help='')

    parser.add_argument('-m', type=str, required=True, dest='stats_mode',
                        help='')

    parser.add_argument('--url', type=str, required=False, dest='base_url',
                        default='http://127.0.0.1:5000/NS/api/',
                        help='')

    parser.add_argument('--taxon', type=str, required=False, dest='taxon',
                        default='0', help='')

    args = parser.parse_args()

    mode = args.stats_mode
    ns_url = args.base_url
    species = args.taxon

    stats_requests.main(mode, ns_url, species)


def main():

    functions_info = {'CreateSchema': ['Create a gene by gene schema based on '
                                       'genomes',
                                       create_schema],
                      'AlleleCall': ['Perform allele call for target genomes',
                                     allele_call],
                      'SchemaEvaluator': ['Tool that builds an html output '
                                          'to better navigate/visualize '
                                          'your schema',
                                          evaluate_schema],
                      'TestGenomeQuality': ['Analyze your allele call output '
                                            'to refine schemas',
                                            test_schema],
                      'ExtractCgMLST': ['Select a subset of loci without '
                                        'missing data (to be used as '
                                        'PHYLOViZ input)',
                                        extract_cgmlst],
                      'RemoveGenes': ['Remove a provided list of loci from '
                                      'your allele call output',
                                      remove_genes],
                      'PrepExternalSchema': ['Adapt an external schema to be '
                                             'used with chewBBACA.',
                                             prep_schema],
                      'JoinProfiles': ['join two profiles in a single profile '
                                       'file',
                                       join_profiles],
                      'UniprotFinder': ['get info about a schema created with '
                                        'chewBBACA',
                                        find_uniprot],
                      'DownloadSchema': ['Download schema from NS',
                                         download_schema_NS],
                      'LoadSchema': ['Upload a schema to the NS',
                                     load_schema_NS],
                      'SyncSchema': ['Syncronize a local schema (downloaded '
                                     'from NS) with NS',
                                     sync_schema_NS],
                      'Send2NS': ['Send local profile and respective alleles '
                                  'to NS',
                                  send_NS],
                      'DownloadProfiles': ['Download all profiles of a given '
                                           'species for a given schema',
                                           down_prof],
                      'SendMetadata': ['send metadata to isolates on the NS',
                                       send_meta],
                      'NSStats': ['',
                                  ns_stats]}

    version = '2.1.0'
    authors = 'Mickael Silva, Pedro Cerqueira, Rafael Mamede'
    repository = 'https://github.com/B-UMMI/chewBBACA'
    contacts = ('rmamede@medicina.ulisboa.pt, '
                'pedro.cerqueira@medicina.ulisboa.pt')

    # Check python version, if fail, exit with message
    try:
        python_version = platform.python_version()
        assert tuple(map(int, python_version.split('.'))) >= (3, 4, 0)

    except AssertionError:
        print('Python version found: {} '.format(platform.python_version()))
        print('Please use version Python >= 3.4')
        sys.exit(0)

    # print help if no command is passed
    if len(sys.argv) == 1:
        print('\n\tUSAGE : chewBBACA.py [module] -h \n')
        print('Select one of the following functions :\n')
        for f in functions_info:
            print('{0}: {1}'.format(f, functions_info[f][0]))
        sys.exit(0)

    if len(sys.argv) > 1 and 'version' in sys.argv[1]:
        print(version)
        return

    print('chewBBACA version: {0}'.format(version))
    print('authors: {0}'.format(authors))
    print('github: {0}'.format(repository))
    print('contacts: {0}'.format(contacts))

    process = sys.argv[1]
    if process in functions_info:
        functions_info[process][1]()
    else:
        print('\n\tUSAGE : chewBBACA.py [module] -h \n')
        print('Select one of the following functions:\n')
        for f in functions_info:
            print('{0}: {1}'.format(f, functions_info[f][0]))


if __name__ == "__main__":
    main()
