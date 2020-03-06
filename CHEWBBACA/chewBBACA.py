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
import pickle
import shutil
import hashlib
import platform
import argparse

from allelecall import BBACA
from createschema import PPanGen
from SchemaEvaluator import ValidateSchema
from PrepExternalSchema import PrepExternalSchema
from utils import (TestGenomeQuality, profile_joiner,
                   uniprot_find, Extract_cgAlleles,
                   RemoveGenes, sqlite_functions as sq,
                   auxiliary_functions as aux,
                   constants as cnts)

from CHEWBBACA_NS import (down_schema, load_schema,
                          sync_schema, down_profiles,
                          send2NS, send_metadata,
                          stats_requests)

import CHEWBBACA


current_version = '2.1.0'


# custom functions to validate arguments type and value
def bsr_type(arg, min_value=cnts.BSR_MIN, max_value=cnts.BSR_MAX):
    """
    """

    try:
        farg = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid float value: {0}'.format(arg))

    if farg < min_value or farg > max_value:
        raise argparse.ArgumentTypeError('value must be > {0} '
                                         'and <= {1}.'.format(min_value, max_value))

    return farg


def minimum_sequence_length_type(arg, min_value=cnts.MSL_MIN, max_value=cnts.MSL_MAX):
    """
    """

    try:
        iarg = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid int value: {0}'.format(arg))

    if iarg < min_value or iarg > max_value:
        raise argparse.ArgumentTypeError('value must be > {0} '
                                         'and <= {1}.'.format(min_value, max_value))

    return iarg


def size_threshold_type(arg, min_value=cnts.ST_MIN, max_value=cnts.ST_MAX):
    """
    """

    try:
        farg = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid float value: {0}'.format(arg))

    if farg < min_value:
        raise argparse.ArgumentTypeError('value must be > {0} '
                                         'and <= {1}.'.format(min_value, max_value))

    return farg


def translation_table_type(arg, genetic_codes=cnts.GENETIC_CODES):
    """
    """

    try:
        iarg = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError('invalid int value: {0}'.format(arg))

    gc_list = list(genetic_codes.keys())
    if arg not in gc_list:
        # format available genetic codes into list
        lines = []
        for k, v in genetic_codes.items():
            new_line = '\t{0}: {1}'.format(k, v)
            lines.append(new_line)

        gc_table = '\n{0}\n'.format('\n'.join(lines))

        raise argparse.ArgumentTypeError('value must correspond to '
                                         'one of the accepted genetic '
                                         'codes\n\nAccepted genetic codes:\n{0}'.format(gc_table))

    return arg


def binary_file_hash(binary_file):
    """
    """

    with open(binary_file, 'rb') as bf:
        file_hash = hashlib.blake2b()
        file_text = bf.read()
        file_hash.update(file_text)
        file_hash = file_hash.hexdigest()

    return file_hash


def write_schema_config(blast_score_ratio, ptf_hash,
                        translation_table, minimum_sequence_length,
                        chewie_version, size_threshold, output_directory):
    """
    """

    params = {}
    params['bsr'] = [blast_score_ratio]
    params['prodigal_training_file'] = [ptf_hash]
    params['translation_table'] = [translation_table]
    params['minimum_locus_length'] = [minimum_sequence_length]
    params['chewBBACA_version'] = [chewie_version]
    params['size_threshold'] = [size_threshold]

    config_file = os.path.join(output_directory, '.schema_config')
    with open(config_file, 'wb') as cf:
        pickle.dump(params, cf)

    return [os.path.isfile(config_file), config_file]


def write_gene_list(schema_dir):
    """
    """

    schema_files = [file for file in os.listdir(schema_dir) if '.fasta' in file]
    schema_list_file = os.path.join(schema_dir, '.genes_list')
    with open(schema_list_file, 'wb') as sl:
        pickle.dump(schema_files, sl)

    return [os.path.isfile(schema_list_file), schema_list_file]


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

    # prog is the name of the program 'chewBBACA.py'
    def __init__(self, prog, indent_increment=2, max_help_position=56, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    # override split lines method
    def _split_lines(self, text, width):
        lines = super()._split_lines(text, width) + ['']
        return lines


def create_schema():

    def msg(name=None):
        # simple command to create schema from genomes
        simple_cmd = ('chewBBACA.py CreateSchema -i <input_files> '
                                                '-o <output_directory> '
                                                '--ptf <ptf_path>')
        # command to create schema from genomes with non-default parameters
        params_cmd = ('chewBBACA.py CreateSchema -i <input_files> '
                                                '-o <output_directory> '
                                                '--ptf <ptf_path>\n'
                                                '\t\t\t    --cpu <cpu_cores> '
                                                '--bsr <blast_score_ratio> '
                                                '--l <minimum_length>\n'
                                                '\t\t\t    --t <translation_table> '
                                                '--st <size_threshold>')
        # command to create schema from single FASTA
        cds_cmd = ('chewBBACA.py CreateSchema -i <input_file> '
                                             '-o <output_directory> '
                                             '--ptf <ptf_path> '
                                             '--CDS')

        usage_msg = ('\nCreate schema from input genomes:\n  {0}\n'
                     '\nCreate schema from input genomes with non-default parameters:\n  {1}\n'
                     '\nCreate schema from single FASTA file:\n  {2}'.format(simple_cmd, params_cmd, cds_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='CreateSchema',
                                     description='Creates a wgMLST '
                                                 'schema based on a '
                                                 'set of input genomes.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('CreateSchema', nargs='+',
                        help='')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='Output directory where the schema will be created.')

    parser.add_argument('-ptf', type=str, required=True,
                        dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('--bsr', type=bsr_type,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', type=minimum_sequence_length_type,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length accepted for a '
                             'coding sequence to be included in the schema.')

    parser.add_argument('--t', type=translation_table_type,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences.')

    parser.add_argument('--st', type=size_threshold_type,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent will be classified as ASM/ALM.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores).')

    parser.add_argument('--b', type=str, required=False,
                        default='blastp', dest='blastp_path',
                        help='Path to the BLASTp executables.')

    parser.add_argument('--CDS', required=False, action='store_true',
                        dest='cds_input',
                        help='Input is a FASTA file with one representative '
                             'sequence per gene in the schema.')

    parser.add_argument('--v', required=False, action='store_true',
                        dest='verbose',
                        help='Increased output verbosity during execution.')

    args = parser.parse_args()

    input_files = args.input_files
    output_directory = args.output_directory
    ptf_path = args.ptf_path
    blast_score_ratio = args.blast_score_ratio
    minimum_length = args.minimum_length
    translation_table = args.translation_table
    size_threshold = args.size_threshold
    cpu_cores = args.cpu_cores
    blastp_path = args.blastp_path
    cds_input = args.cds_input
    verbose = args.verbose

    # check if ptf exists
    if os.path.isfile(ptf_path) is False:
        message = ('Cannot find specified Prodigal training file.\nPlease provide a '
                   'valid training file.\n\nYou can create a training '
                   'file for a species of interest with the following command:\n  '
                   'prodigal -i <reference_genome> -t <training_file.trn> -p single\n\n'
                   'It is strongly advised to provide a high-quality and closed genome '
                   'for the training process.')
        sys.exit(message)

    if cds_input is True:
        input_files = [os.path.abspath(input_files)]
    else:
        input_files = aux.check_input_type(input_files, 'listGenomes2Call.txt')

    # start CreateSchema process
    PPanGen.main(input_files, cpu_cores, output_directory,
                 blast_score_ratio, blastp_path, minimum_length,
                 verbose, ptf_path, cds_input,
                 translation_table, size_threshold)

    # copy training file to schema directory
    shutil.copy(ptf_path, output_directory)

    # determine PTF checksum
    ptf_hash = binary_file_hash(ptf_path)

    # write schema config file
    schema_config = write_schema_config(blast_score_ratio, ptf_hash,
                                        translation_table, minimum_length,
                                        current_version, size_threshold,
                                        output_directory)

    # create hidden file with genes/loci list
    genes_list_file = write_gene_list(output_directory)

    # remove temporary file with paths
    # to genome files
    if os.path.isfile(input_files):
        os.remove(input_files)


def allele_call():

    def msg(name=None):
        # simple command to perform AlleleCall with schema deafult parameters
        simple_cmd = ('chewBBACA.py AlleleCall -i <input_files> '
                                              '-g <schema_directory> '
                                              '-o <output_directory> ')
        # command to perform AlleleCall with non-default parameters
        params_cmd = ('chewBBACA.py AlleleCall -i <input_files> '
                                              '-g <schema_directory> '
                                              '-o <output_directory> '
                                              '--ptf <ptf_path>\n'
                                              '\t\t\t  --cpu <cpu_cores> '
                                              '--bsr <blast_score_ratio> '
                                              '--l <minimum_length>\n'
                                              '\t\t\t  --t <translation_table> '
                                              '--st <size_threshold>')
        # command to perform AlleleCall with single Fasta file
        # cds_cmd = ('chewBBACA.py AlleleCall -i <input_file> '
        #                                    '-o <output_directory> '
        #                                    '--ptf <ptf_path> '
        #                                    '--CDS')

        usage_msg = ('\nPerform AlleleCall with schema default parameters:\n  {0}\n'
                     '\nPerform AlleleCall with non-default parameters:\n  {1}\n'.format(simple_cmd, params_cmd))
                     #'\nPerform AlleleCall with single FASTA file that contains coding sequences:\n  {2}'.format(simple_cmd, params_cmd, cds_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='AlleleCall',
                                     description='Performs allele calling to determine the '
                                                 'allelic profiles of a set of input genomes. '
                                                 'The process identifies new alleles, assigns '
                                                 'an integer identifier to those alleles and '
                                                 'adds them to the schema.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter,
                                     epilog='It is strongly advised to perform AlleleCall '
                                            'with the default schema parameters to ensure '
                                            'more consistent results.')

    parser.add_argument('AlleleCall', nargs='+', help='')

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory with the genomes FASTA '
                             'files or to a file with a list of paths to '
                             'the FASTA files, one per line.')

    parser.add_argument('-g', type=str, required=True,
                        dest='schema_directory',
                        help='Path to the schema directory with the'
                             ' genes FASTA files.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='Output directory where the allele '
                             'calling results will be stored.')

    parser.add_argument('--ptf', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file. '
                             'Default is to get training file from '
                             'schema directory.')

    parser.add_argument('--bsr', type=bsr_type, required=False,
                        default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--t', type=translation_table_type, required=False,
                        default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences '
                             '(default=11).')

    parser.add_argument('--st', type=size_threshold_type, required=False,
                        default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent will be classified as ASM/ALM')

    parser.add_argument('--cpu', type=int, required=False, default=1,
                        dest='cpu_cores',
                        help='Number of CPU cores/threads that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores/threads).')

    parser.add_argument('--b', type=str, required=False,
                        default='blastp', dest='blastp_path',
                        help='Path to the BLASTp executables.')

    parser.add_argument('--contained', action='store_true', required=False,
                        default=False, dest='contained',
                        help=argparse.SUPPRESS)

    parser.add_argument('--CDS', action='store_true', required=False,
                        default=False, dest='cds_input',
                        help=argparse.SUPPRESS)

    parser.add_argument('--json', action='store_true', required=False,
                        dest='json_report',
                        help='Output report in JSON format.')

    parser.add_argument('--fc', action='store_true', required=False,
                        dest='force_continue',
                        help='Continue allele call process that '
                             'was interrupted.')

    parser.add_argument('--fr', action='store_true', required=False,
                        dest='force_reset',
                        help='Force process reset even if there '
                             'are temporary files from a previous '
                             'process that was interrupted.')

    parser.add_argument('--v', action='store_true',
                        dest='verbose',
                        help='Increased output verbosity during execution.')

    args = parser.parse_args()

    input_files = args.input_files
    schema_directory = args.schema_directory
    output_directory = args.output_directory
    ptf_path = args.ptf_path
    blast_score_ratio = args.blast_score_ratio
    translation_table = args.translation_table
    size_threshold = args.size_threshold
    cpu_cores = args.cpu_cores
    blastp_path = args.blastp_path
    contained = args.contained
    cds_input = args.cds_input
    json_report = args.json_report
    force_continue = args.force_continue
    force_reset = args.force_reset
    verbose = args.verbose
    chosen_taxon = False

    schema_genes = aux.check_input_type(schema_directory, 'listGenes2Call.txt')

    # check parameters values in config file and alter if needed
    config_file = os.path.join(schema_directory, '.schema_config')
    with open(config_file, 'rb') as pf:
        params = pickle.load(pf)

    unmatch_params = {}
    # check bsr
    if blast_score_ratio not in params['bsr']:
        unmatch_params['bsr'] = blast_score_ratio
    # check chewie version
    if current_version not in params['chewBBACA_version']:
        unmatch_params['chewBBACA_version'] = current_version
    # check size threshold value
    if size_threshold not in params['size_threshold']:
        unmatch_params['size_threshold'] = size_threshold
    if translation_table not in params['translation_table']:
        unmatch_params['translation_table'] = translation_table

    if len(unmatch_params) > 0:
        print('Provided arguments values differ from arguments values used for schema creation:')
        print('\n'.join(map(str, list(unmatch_params.values()))))
        params_answer = input('Continuing might invalidate the schema. Continue?\n')
        if params_answer.lower() not in ['y', 'yes']:
            sys.exit(0) 
        else:
            for p in unmatch_params:
                params[p].append(unmatch_params[p])

    # default is to get the training file in schema directory
    print(ptf_path)
    if ptf_path is False:
        for file in os.listdir(schema_directory):
            if file.endswith('.trn'):
                ptf_path = os.path.join(schema_directory, file)
        if os.path.isfile(ptf_path) is False:
            sys.exit('There is no valid training file in schema directory.')
    # if user provides a training file
    else:
        if os.path.isfile(ptf_path) is False:
            sys.exit('Provided Prodigal training file does not exist.')

    # determine PTF checksum
    ptf_hash = binary_file_hash(ptf_path)

    if ptf_hash not in params['prodigal_training_file']:
        ptf_num = len(params['prodigal_training_file'])
        if ptf_num == 1:
            print('Prodigal training file is not the one used to create the schema.')
            ptf_answer = input('Using this training file will invalidate the schema.\nContinue process?\n')
        if ptf_num > 1:
            print('Prodigal training file is not any of the {0} used in previous runs.'.format(ptf_num))
            ptf_answer = input('Continue?\n')

        if ptf_answer.lower() not in ['y', 'yes']:
            sys.exit(0)
        else:
            params['prodigal_training_file'].append(ptf_hash)
            unmatch_params['prodigal_training_file'] = ptf_hash
    print(ptf_path)

    # save updated schema config file
    if len(unmatch_params) > 0:
        with open(config_file, 'wb') as cf:
            pickle.dump(params, cf)

    # if is a fasta pass as a list of genomes with a single genome,
    # if not check if is a folder or a txt with a list of paths
    genomes_files = aux.check_input_type(input_files, 'listGenomes2Call.txt')

    BBACA.main(genomes_files, schema_genes, cpu_cores,
               output_directory, blast_score_ratio,
               blastp_path, force_continue, json_report,
               verbose, force_reset, contained, chosen_taxon,
               ptf_path, cds_input, size_threshold,
               translation_table)

    # add profiles to SQLite database
    results_folder = os.path.join(output_directory, [f for f in os.listdir(output_directory) if 'results' in f][0])
    results_matrix = os.path.join(results_folder, 'results_alleles.tsv')
    insert_date = results_folder.split('_')[-1]

    # verify that database directory exists
    database_directory = os.path.join(schema_directory, 'profiles_database')
    # create if it does not exist
    if os.path.isdir(database_directory) is False:
        os.mkdir(database_directory)

    # also need to check for database file
    database_file = os.path.join(database_directory, 'profiles.db')
    if os.path.isfile(database_file) is False:
        sq.create_database_structure(database_file)
        # insert loci list into loci table
        sq.insert_loci(database_file, results_matrix)
    
    # insert whole matrix
    a = sq.insert_allelecall_matrix(results_matrix, database_file, insert_date)

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

        # simple command to adapt external schema with default arguments values
        simple_cmd = ('chewBBACA.py PrepExternalSchema -i <input_files> '
                                                      '-o <output_directory> '
                                                      '-ptf <ptf_path> ')

        # command to adapt external schema with non-default arguments values
        params_cmd = ('chewBBACA.py PrepExternalSchema -i <input_files> '
                                                      '-o <output_directory> '
                                                      '-ptf <ptf_path> ')

        usage_msg = ('\nAdapt external schema (one FASTA file per schema gene):\n{0}\n'
                     '\nAdapt external schema with non-default parameters:\n{1}\n'.format(simple_cmd, params_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='PrepExternalSchema',
                                     description='This script enables the '
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
                                                 'directory.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

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

    parser.add_argument('-ptf', type=str, required=True,
                        dest='ptf_path',
                        help='Path to the Prodigal training file. '
                             'Default is to get training file in '
                             'schema directory.')

    parser.add_argument('--bsr', type=bsr_type,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='The BLAST Score Ratio value that will be '
                        'used to adapt the external schema.')

    parser.add_argument('--l', type=minimum_sequence_length_type,
                        required=False, default=0, dest='minimum_length',
                        help='Minimum sequence length accepted. Sequences with'
                        ' a length value smaller than the value passed to this'
                        ' argument will be discarded.')

    parser.add_argument('--t', type=translation_table_type,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code to use for CDS translation.')

    parser.add_argument('--st', type=size_threshold_type,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent will be classified as ASM/ALM.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='The number of CPU cores to use.')

    args = parser.parse_args()

    input_files = args.input_files
    output_directory = args.output_directory
    ptf_path = args.ptf_path
    blast_score_ratio = args.blast_score_ratio
    minimum_length = args.minimum_length
    translation_table = args.translation_table
    size_threshold = args.size_threshold
    cpu_cores = args.cpu_cores

    # check if ptf exists
    if os.path.isfile(ptf_path) is False:
        message = ('Cannot find specified Prodigal training file.\nPlease provide a '
                   'valid training file.\n\nYou can create a training '
                   'file for a species of interest with the following command:\n  '
                   'prodigal -i <reference_genome> -t <training_file.trn> -p single\n\n'
                   'It is strongly advised to provide a high-quality and closed genome '
                   'for the training process.')
        sys.exit(message)

    PrepExternalSchema.main(input_files, output_directory, cpu_cores,
                            blast_score_ratio, minimum_length,
                            translation_table, ptf_path,
                            size_threshold)

    # copy training file to schema directory
    shutil.copy(ptf_path, output_directory)

    # determine PTF checksum
    ptf_hash = binary_file_hash(ptf_path)

    # write schema config file
    schema_config = write_schema_config(blast_score_ratio, ptf_hash,
                                        translation_table, minimum_length,
                                        current_version, size_threshold,
                                        output_directory)

    # create hidden file with genes/loci list
    genes_list_file = write_gene_list(output_directory)


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
    wiki = 'https://github.com/B-UMMI/chewBBACA/wiki'
    tutorial = 'https://github.com/B-UMMI/chewBBACA_tutorial'
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
        print('\n\tUSAGE: chewBBACA.py [module] -h \n')
        print('Select one of the following functions :\n')
        for f in functions_info:
            print('{0}: {1}'.format(f, functions_info[f][0]))
        sys.exit(0)

    if len(sys.argv) > 1 and 'version' in sys.argv[1]:
        print(version)
        return

    print('\nchewBBACA version: {0}'.format(version))
    print('Authors: {0}'.format(authors))
    print('Github: {0}'.format(repository))
    print('Wiki: {0}'.format(wiki))
    print('Tutorial: {0}'.format(tutorial))
    print('Contacts: {0}\n'.format(contacts))

    process = sys.argv[1]
    if process in functions_info:
        functions_info[process][1]()
    else:
        print('\n\tUSAGE: chewBBACA.py [module] -h \n')
        print('Select one of the following functions:\n')
        for f in functions_info:
            print('{0}: {1}'.format(f, functions_info[f][0]))


if __name__ == "__main__":
    main()
