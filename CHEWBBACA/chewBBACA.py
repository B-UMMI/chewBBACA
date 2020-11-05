#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This is the main script of the chewBBACA suite.

"""


import os
import sys
import time
import pickle
import shutil
import datetime
import platform
import argparse

try:
    from __init__ import __version__
    from allelecall import BBACA
    from createschema import PPanGen
    from SchemaEvaluator import ValidateSchema
    from PrepExternalSchema import PrepExternalSchema
    from utils import (TestGenomeQuality, profile_joiner,
                       uniprot_find, Extract_cgAlleles,
                       RemoveGenes, sqlite_functions as sq,
                       auxiliary_functions as aux,
                       constants as cnst,
                       parameters_validation as pv)

    from utils.parameters_validation import ModifiedHelpFormatter

    from CHEWBBACA_NS import (down_schema, load_schema,
                              sync_schema, stats_requests)
except:
    from CHEWBBACA import __version__
    from CHEWBBACA.allelecall import BBACA
    from CHEWBBACA.createschema import PPanGen
    from CHEWBBACA.SchemaEvaluator import ValidateSchema
    from CHEWBBACA.PrepExternalSchema import PrepExternalSchema
    from CHEWBBACA.utils import (TestGenomeQuality, profile_joiner,
                                 uniprot_find, Extract_cgAlleles,
                                 RemoveGenes, sqlite_functions as sq,
                                 auxiliary_functions as aux,
                                 constants as cnst,
                                 parameters_validation as pv)

    from CHEWBBACA.utils.parameters_validation import ModifiedHelpFormatter

    from CHEWBBACA.CHEWBBACA_NS import (down_schema, load_schema,
                                        sync_schema, stats_requests)


version = __version__


def create_schema():

    def msg(name=None):
        # simple command to create schema from genomes
        simple_cmd = ('chewBBACA.py CreateSchema -i <input_files> '
                                                '-o <output_directory> '
                                                '-ptf <ptf_path>')
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

    parser.add_argument('--ptf', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('--bsr', type=pv.bsr_type,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', type=pv.minimum_sequence_length_type,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length accepted for a '
                             'coding sequence to be included in the schema.')

    parser.add_argument('--t', type=pv.translation_table_type,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences.')

    parser.add_argument('--st', type=pv.size_threshold_type,
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

    header = 'chewBBACA - CreateSchema'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

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
    if ptf_path is not False:
        ptf_val = aux.check_ptf(ptf_path)
        if ptf_val[0] is False:
            sys.exit(ptf_val[1])

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
    if ptf_path is not False:
        shutil.copy(ptf_path, output_directory)
        # determine PTF checksum
        ptf_hash = aux.hash_file(ptf_path, 'rb')
    else:
        ptf_hash = ''

    # write schema config file
    schema_config = aux.write_schema_config(blast_score_ratio, ptf_hash,
                                            translation_table, minimum_length,
                                            version, size_threshold,
                                            output_directory)

    # create hidden file with genes/loci list
    genes_list_file = aux.write_gene_list(output_directory)

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

    parser.add_argument('--bsr', type=pv.bsr_type, required=False,
                        default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', type=pv.minimum_sequence_length_type,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length accepted for a '
                             'coding sequence to be included in the schema.')

    parser.add_argument('--t', type=pv.translation_table_type, required=False,
                        default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences '
                             '(default=11).')

    parser.add_argument('--st', type=pv.size_threshold_type, required=False,
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

    parser.add_argument('--db', required=False, action='store_false',
                        dest='store_profiles',
                        help='If the profiles in the output matrix '
                             'should be stored in the local SQLite '
                             'database.')

    parser.add_argument('--v', required=False, action='store_true',
                        dest='verbose',
                        help='Increased output verbosity during execution.')

    args = parser.parse_args()

    header = 'chewBBACA - AlleleCall'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

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
    store_profiles = args.store_profiles
    verbose = args.verbose
    minimum_length = args.minimum_length

    timeout = 30

    # change schema's directory value if user provided file with list of genes
    genes_list = schema_directory
    if os.path.isfile(schema_directory) is True:
        with open(genes_list, 'r') as f:
            schema_directory = f.readline().strip()
            schema_directory = os.path.dirname(schema_directory)

    # check parameters values in config file and alter if needed
    config_file = os.path.join(schema_directory, '.schema_config')

    # legacy schemas do not have config file, create one if user wants to continue
    if os.path.isfile(config_file) is False:
        if force_continue is False:
            prompt = ('It seems that your schema was created with chewBBACA 2.1.0 or lower.\n'
                      'It is highly recommended that you run the PrepExternalSchema '
                      'process to guarantee full compatibility with the new chewBBACA '
                      'version.\nIf you wish to continue, the AlleleCall process will convert '
                      'the schema to v{0}, but will not determine if schema structure respects '
                      'configuration values.\nDo you wish to proceed?\n'.format(version))
            proceed = aux.input_timeout(prompt, timeout)
        else:
            proceed = 'yes'

        if proceed.lower() not in ['y', 'yes']:
            sys.exit('Exited.')
        else:
            if ptf_path is not False:
                print('\nAdding Prodigal training file to schema...')
                ptf_val = aux.check_ptf(ptf_path)
                if ptf_val[0] is False:
                    sys.exit(ptf_val[1])
                # copy training file to schema directory
                shutil.copy(ptf_path, schema_directory)
                print('Created {0}'.format(os.path.join(schema_directory,
                                                        os.path.basename(ptf_path))))

            # determine PTF checksum
            if ptf_path is not False:
                ptf_hash = aux.hash_file(ptf_path, 'rb')
            else:
                ptf_hash = ''

            print('\nCreating file with schema configs...')
            # write schema config file
            schema_config = aux.write_schema_config(blast_score_ratio, ptf_hash,
                                                    translation_table, minimum_length,
                                                    version, size_threshold,
                                                    schema_directory)
            print('Created {0}'.format(os.path.join(schema_directory,
                                                    '.schema_config')))

            print('\nCreating file with list of genes...')
            # create hidden file with genes/loci list
            genes_list_file = aux.write_gene_list(schema_directory)
            print('Created {0}\n'.format(os.path.join(schema_directory,
                                                      '.genes_list')))

    # read schema configs
    with open(config_file, 'rb') as pf:
        schema_params = pickle.load(pf)

    # run parameters values
    run_params = {'bsr': blast_score_ratio,
                  'chewBBACA_version': version,
                  'prodigal_training_file': ptf_path,
                  'translation_table': translation_table,
                  'size_threshold': size_threshold}

    # mismatched schema and run parameters values
    unmatch_params = {k: v for k, v in run_params.items()
                      if v not in schema_params[k] and k not in ['prodigal_training_file', 'chewBBACA_version']}

    if len(unmatch_params) > 0:
        print('Provided arguments values differ from arguments '
              'values used for schema creation:\n')
        params_diffs = [[p, ':'.join(map(str, schema_params[p])), str(unmatch_params[p])] for p in unmatch_params]
        params_diffs_text = ['{:^20} {:^20} {:^10}'.format('Argument', 'Schema', 'Provided')]
        params_diffs_text += ['{:^20} {:^20} {:^10}'.format(p[0], p[1], p[2]) for p in params_diffs]
        print('\n'.join(params_diffs_text))
        if force_continue is False:
            prompt = ('\nContinuing might lead to results not consistent with '
                      'previous runs.\nProviding parameters values that differ '
                      'from the ones used for schema creation will also '
                      'invalidate the schema for uploading and synchronization '
                      'with the Chewie-NS.\nContinue?\n')
            params_answer = aux.input_timeout(prompt, timeout)
        else:
            params_answer = 'yes'

        if params_answer.lower() not in ['y', 'yes']:
            sys.exit('Exited.')
        else:
            for p in unmatch_params:
                schema_params[p].append(unmatch_params[p])

    # default is to get the training file in schema directory
    if run_params['prodigal_training_file'] is False:
        # deal with multiple training files
        schema_ptfs = [file for file in os.listdir(schema_directory) if file.endswith('.trn')]
        if len(schema_ptfs) > 1:
            sys.exit('Found more than one Prodigal training file in schema directory.\n'
                     'Please maintain only the training file used in the schema creation process.')
        elif len(schema_ptfs) == 1:
            ptf_path = os.path.join(schema_directory, schema_ptfs[0])
        elif len(schema_ptfs) == 0:
            print('There is no Prodigal training file in schema\'s directory.')
    # if user provides a training file
    else:
        ptf_val = aux.check_ptf(ptf_path)
        if ptf_val[0] is False:
            sys.exit(ptf_val[1])

    # determine PTF checksum
    if ptf_path is not False:
        ptf_hash = aux.hash_file(ptf_path, 'rb')
    else:
        ptf_hash = ''

    if ptf_hash not in schema_params['prodigal_training_file']:
        ptf_num = len(schema_params['prodigal_training_file'])
        if force_continue is False:
            if ptf_num == 1:
                print('Prodigal training file is not the one used to create the schema.')
                prompt = ('Using this training file might lead to results not '
                          'consistent with previous runs and invalidate the '
                          'schema for usage with the NS.\nContinue process?\n')
                ptf_answer = aux.input_timeout(prompt, timeout)
            if ptf_num > 1:
                print('Prodigal training file is not any of the {0} used in previous runs.'.format(ptf_num))
                prompt = ('Continue?\n')
                ptf_answer = aux.input_timeout(prompt, timeout)
        else:
            ptf_answer = 'yes'

        if ptf_answer.lower() not in ['y', 'yes']:
            sys.exit('Exited.')
        else:
            schema_params['prodigal_training_file'].append(ptf_hash)
            unmatch_params['prodigal_training_file'] = ptf_hash

    # save updated schema config file
    if len(unmatch_params) > 0:
        with open(config_file, 'wb') as cf:
            pickle.dump(schema_params, cf)

    # if is a fasta pass as a list of genomes with a single genome,
    # if not check if is a folder or a txt with a list of paths
    schema_genes = aux.check_input_type(genes_list, 'listGenes2Call.txt')
    genomes_files = aux.check_input_type(input_files, 'listGenomes2Call.txt')

    # determine if schema was downloaded from the Chewie-NS
    ns_config = os.path.join(schema_directory, '.ns_config')
    ns = os.path.isfile(ns_config)

    BBACA.main(genomes_files, schema_genes, cpu_cores,
               output_directory, blast_score_ratio,
               blastp_path, force_continue, json_report,
               verbose, force_reset, contained,
               ptf_path, cds_input, size_threshold,
               translation_table, ns)

    if store_profiles is True:
        # add profiles to SQLite database
        # parent results folder might have several results folders
        results_folders = [os.path.join(output_directory, file)
                           for file in os.listdir(output_directory) if 'results' in file]
        # create datetime objects and sort to get latest
        insert_dates = [(file, datetime.datetime.strptime(file.split('_')[-1], '%Y%m%dT%H%M%S'))
                        for file in results_folders]
        sorted_insert_dates = sorted(insert_dates, key=lambda x: x[1], reverse=True)
        results_matrix = os.path.join(sorted_insert_dates[0][0], 'results_alleles.tsv')
        insert_date = sorted_insert_dates[0][0].split('_')[-1]

        # verify that database directory exists
        database_directory = os.path.join(schema_directory, 'profiles_database')
        # create if it does not exist
        if os.path.isdir(database_directory) is False:
            os.mkdir(database_directory)

        # also need to check for database file
        database_file = os.path.join(database_directory, 'profiles.db')
        if os.path.isfile(database_file) is False:
            print('\nCreating SQLite database to store profiles...', end='')
            try:
                sq.create_database(database_file)
                # insert loci list into loci table
                total_loci = sq.insert_loci(database_file, results_matrix)
                print('done.')
                print('Inserted {0} loci into database.'.format(total_loci))
            except Exception:
                print('WARNING: Could not create database file. Will not store profiles.')

        # insert whole matrix
        if os.path.isfile(database_file) is not False:
            print('\nSending allelic profiles to SQLite database...', end='')
            try:
                total_profiles = sq.insert_allelecall_matrix(results_matrix, database_file, insert_date)
                print('done.')
                print('Inserted {0} profiles ({1} total, {2} total unique).'.format(total_profiles[0], total_profiles[1], total_profiles[2]))
            except Exception as e:
                print(e)
                print('WARNING: Could not store profiles in local database.')

    # remove temporary files with paths to genomes
    # and schema files files
    if os.path.isfile(schema_genes) is True:
        os.remove(schema_genes)
    if os.path.isfile(genomes_files):
        os.remove(genomes_files)


def evaluate_schema():

    def msg(name=None):
        # simple command to evaluate schema or set of loci
        simple_cmd = ('chewBBACA.py SchemaEvaluator -i <input_files> '
                                                   '-l <output_file> '
                                                   '--cpu <cpu_cores>')

        usage_msg = ('\nEvaluate schema with default parameters:\n  {0}\n'.format(simple_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='SchemaEvaluator',
                                     description='Evaluate the number of alelles and allele size '
                                                 'variation for the loci in a schema or for a set '
                                                 'of selected loci. Provide information about '
                                                 'problematic alleles per locus and individual pages '
                                                 'for each locus with a plot with allele size, a Neighbor '
                                                 'Joining tree based on a multiple sequence alignment (MSA) '
                                                 'and a visualization of the MSA.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('SchemaEvaluator', nargs='+',
                        help='Evaluates a set of loci.')

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='Path to the schema\'s directory or path to a file containing '
                             'the paths to the FASTA files of the loci that will be evaluated, '
                             'one per line.')

    parser.add_argument('-l', type=str, required=True,
                        dest='output_file',
                        help='Path to the output HTML file.')

    parser.add_argument('-p', action='store_true', required=False,
                        default=False, dest='conserved',
                        help='If all alleles must be within the threshold for the '
                             'locus to be considered as having low length variability.')

    parser.add_argument('--log', action='store_true', default=False,
                        dest='log_scale',
                        help='Apply log scale transformation to the yaxis '
                             'of the plot with the number of alleles per locus.')

    parser.add_argument('-ta', type=int, required=False,
                        default=11, dest='translation_table',
                        help='Genetic code used to translate coding '
                             'sequences.')

    parser.add_argument('-t', type=float, required=False,
                        default=0.05, dest='threshold',
                        help='Allele size variation threshold. If an allele has '
                             'a size within the interval of the locus mode -/+ '
                             'the threshold, it will be considered a conserved '
                             'allele.')

    parser.add_argument('--title', type=str, required=False,
                        default='My Analyzed wg/cg MLST Schema - Rate My Schema',
                        dest='title',
                        help='Title displayed on the html page.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores to use to run the process.')

    parser.add_argument('-s', type=int, required=False,
                        default=500, dest='split_range',
                        help='Number of boxplots displayed in the plot area (more than '
                             '500 can lead to performance issues).')

    parser.add_argument('--light', action='store_true', required=False,
                        default=False, dest='light_mode',
                        help='Skip clustal and mafft.')

    args = parser.parse_args()

    header = 'chewBBACA - SchemaEvaluator'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    input_files = args.input_files
    output_file = args.output_file
    log_scale = args.log_scale
    translation_table = args.translation_table
    cpu_cores = args.cpu_cores
    threshold = args.threshold
    conserved = args.conserved
    split_range = args.split_range
    light_mode = args.light_mode
    title = str(args.title)

    ValidateSchema.main(input_files, cpu_cores, output_file,
                        translation_table, threshold, split_range,
                        title, log_scale, conserved,
                        light_mode)


def test_schema():

    def msg(name=None):
        # simple command to evaluate genome quality
        simple_cmd = ('chewBBACA.py TestGenomeQuality -i <input_file> '
                                                   '-n <max_iteration> '
                                                   '-t <max_threshold>'
                                                   '-s <step>')

        usage_msg = ('\nEvaluate genome quality with default parameters:\n  {0}\n'.format(simple_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(description='This process evaluates the quality of genomes '
                                                 'based on the results of the AlleleCall process.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('TestGenomeQuality', nargs='+',
                        help='Evaluate the quality of input genomes based '
                             'on allele calling results.')

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='Path to file with a matrix of allelic profiles.')

    parser.add_argument('-n', type=int, required=True,
                        dest='max_iteration',
                        help='Maximum number of iterations.')

    parser.add_argument('-t', type=int, required=True,
                        dest='max_threshold',
                        help='Maximum threshold of bad calls above 95 percent.')

    parser.add_argument('-s', type=int, required=True,
                        dest='step',
                        help='Step between each threshold analysis.')

    parser.add_argument('-o', type=str, required=False,
                        default='.', dest='output_directory',
                        help='Path to the output directory that will '
                             'store output files')

    parser.add_argument('-v', '--verbose', action='store_true',
                        default=False, dest='verbose',
                        help='Increase stdout verbosity.')

    args = parser.parse_args()

    header = 'chewBBACA - TestGenomeQuality'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    input_file = args.input_file
    max_iteration = args.max_iteration
    max_threshold = args.max_threshold
    step = args.step
    output_directory = args.output_directory
    verbose = args.verbose

    TestGenomeQuality.main(input_file, max_iteration,
                           max_threshold, step,
                           output_directory, verbose)


def extract_cgmlst():

    def msg(name=None):
        # simple command to determine loci that constitute cgMLST
        simple_cmd = ('  chewBBACA.py ExtractCgMLST -i <input_file> '
                                                   '-o <output_directory> ')

        # command to determine cgMLST with custom threshold
        threshold_cmd = ('  chewBBACA.py ExtractCgMLST -i <input_file> '
                                                   '-o <output_directory> '
                                                   '\n\t\t\t     --p <threshold>')

        # command to get information about a single schema
        remove_cmd = ('  chewBBACA.py ExtractCgMLST -i <input_file> '
                                                   '-o <output_directory> '
                                                   '\n\t\t\t     --r <genes2remove> '
                                                   '--g <genomes2remove>')

        usage_msg = ('\nDetermine cgMLST:\n{0}\n'
                     '\nDetermine cgMLST based on non-default threshold:\n{1}\n'
                     '\nRemove genes and genomes from matrix:\n{2}\n'
                     ''.format(simple_cmd, threshold_cmd, remove_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='ExtractCgMLST',
                                     description='Determines the set of '
                                                 'loci that constitute the '
                                                 'core genome based on a '
                                                 'threshold.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('ExtractCgMLST', nargs='+',
                        help='Determines the set of '
                             'loci that constitute the '
                             'core genome based on a '
                             'threshold.')

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='Path to input file containing a matrix with '
                             'allelic profiles.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='Path to the directory where the process '
                             'will store output files.')

    parser.add_argument('--p', '-p', type=float, required=False,
                        default=1, dest='threshold',
                        help='Genes that constitute the core genome '
                             'must be in a proportion of genomes that is '
                             'at least equal to this value.')

    parser.add_argument('--r', '-r', type=str, required=False,
                        default=False, dest='genes2remove',
                        help='Path to file with a list of genes/columns to '
                             'remove from the matrix (one gene identifier '
                             'per line).')

    parser.add_argument('--g', '-g', type=str, required=False,
                        default=False, dest='genomes2remove',
                        help='Path to file with a list of genomes/rows to '
                             'remove from the matrix (one genome identifier '
                             'per line).')

    args = parser.parse_args()

    header = 'chewBBACA - ExtractCgMLST'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    input_file = args.input_file
    output_directory = args.output_directory
    threshold = args.threshold
    genes2remove = args.genes2remove
    genomes2remove = args.genomes2remove

    Extract_cgAlleles.main(input_file, output_directory, threshold,
                           genes2remove, genomes2remove)


def remove_genes():

    def msg(name=None):

        # simple command to remove a set of genes from a matrix with allelic profiles
        simple_cmd = ('  chewBBACA.py RemoveGenes -i <input_file> '
                                                   '-g <genes_list> '
                                                   '-o <output_file>')

        usage_msg = ('\nRemove a set of genes from a matrix with allelic profiles:\n{0}\n'.format(simple_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(description='Remove loci from a matrix with allelic profiles.',
                                     usage=msg())

    parser.add_argument('RemoveGenes', nargs='+',
                        help='Remove loci from a matrix with allelic profiles.')

    parser.add_argument('-i', type=str, required=True,
                        dest='input_file',
                        help='TSV file that contains a matrix with allelic profiles '
                             'determined by the AlleleCall process.')

    parser.add_argument('-g', type=str, required=True,
                        dest='genes_list',
                        help='File with the list of genes to remove.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_file',
                        help='Path to the output file that will be created with the '
                             'new matrix.')

    parser.add_argument('--inverse', action='store_true', default=False,
                        dest='inverse',
                        help='List of genes that is provided is the list of genes to '
                             'keep and all other genes should be removed.')

    args = parser.parse_args()

    header = 'chewBBACA - RemoveGenes'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    input_file = args.input_file
    genes_list = args.genes_list
    output_file = args.output_file
    inverse = args.inverse

    RemoveGenes.main(input_file, genes_list, output_file, inverse)


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

    header = 'chewBBACA - JoinProfiles'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    profile1 = args.p1
    profile2 = args.p2
    outputFile = args.o

    profile_joiner.main(profile1, profile2, outputFile)


def prep_schema():

    def msg(name=None):

        # simple command to adapt external schema with default arguments values
        simple_cmd = ('  chewBBACA.py PrepExternalSchema -i <input_files> '
                                                      '-o <output_directory> '
                                                      '--ptf <ptf_path> ')

        # command to adapt external schema with non-default arguments values
        params_cmd = ('  chewBBACA.py PrepExternalSchema -i <input_files> '
                                                      '-o <output_directory> '
                                                      '--ptf <ptf_path>\n'
                                                      '\t\t\t\t  --cpu <cpu_cores> '
                                                      '--bsr <blast_score_ratio> '
                                                      '--l <minimum_length>\n'
                                                      '\t\t\t\t  --t <translation_table> '
                                                      '--st <size_threshold>')

        usage_msg = ('\nAdapt external schema (one FASTA file per schema gene):\n\n{0}\n'
                     '\nAdapt external schema with non-default parameters:\n\n{1}\n'.format(simple_cmd, params_cmd))

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

    parser.add_argument('--ptf', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file that '
                             'will be associated with the adapted schema.')

    parser.add_argument('--bsr', type=pv.bsr_type,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='The BLAST Score Ratio value that will be '
                             'used to adapt the external schema (default=0.6).')

    parser.add_argument('--l', type=pv.minimum_sequence_length_type,
                        required=False, default=0, dest='minimum_length',
                        help='Minimum sequence length accepted. Sequences with'
                             ' a length value smaller than the value passed to this'
                             ' argument will be discarded (default=0).')

    parser.add_argument('--t', type=pv.translation_table_type,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code to use for CDS translation.'
                             ' (default=11, for Bacteria and Archaea)')

    parser.add_argument('--st', type=pv.size_threshold_type,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent when compared to the representative '
                             'will not be included in the final schema.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='The number of CPU cores to use (default=1).')

    args = parser.parse_args()

    header = 'chewBBACA - PrepExternalSchema'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    input_files = args.input_files
    output_directory = args.output_directory
    ptf_path = args.ptf_path
    blast_score_ratio = args.blast_score_ratio
    minimum_length = args.minimum_length
    translation_table = args.translation_table
    size_threshold = args.size_threshold
    cpu_cores = args.cpu_cores

    # check if ptf exists
    if ptf_path is not False:
        ptf_val = aux.check_ptf(ptf_path)
        if ptf_val[0] is False:
            sys.exit(ptf_val[1])

    PrepExternalSchema.main(input_files, output_directory, cpu_cores,
                            blast_score_ratio, minimum_length,
                            translation_table, ptf_path,
                            size_threshold)

    # copy training file to schema directory
    if ptf_path is not False:
        shutil.copy(ptf_path, output_directory)

    # determine PTF checksum
    if ptf_path is not False:
        ptf_hash = aux.hash_file(ptf_path, 'rb')
    else:
        ptf_hash = ''

    # write schema config file
    schema_config = aux.write_schema_config(blast_score_ratio, ptf_hash,
                                            translation_table, minimum_length,
                                            version, size_threshold,
                                            output_directory)

    # create hidden file with genes/loci list
    genes_list_file = aux.write_gene_list(output_directory)


def find_uniprot():

    def msg(name=None):

        # simple command to determine annotations for the loci in a schema
        simple_cmd = ('  chewBBACA.py UniprotFinder -i <input_files> '
                                                   '-t <protein_table> '
                                                   '--cpu <cpu_cores>')

        usage_msg = ('\nFind annotations for loci in a schema:\n\n{0}\n'.format(simple_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='UniprotFinder',
                                     description='This process determines loci annotations based '
                                                 'on exact matches found in the UniProt database.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('UniprotFinder', nargs='+',
                        help='Determine annotations for loci in a schema.')

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='Path to the schema\'s directory or to a file with '
                             'a list of paths to loci FASTA files, one per line.')

    parser.add_argument('-t', type=str, required=True,
                        dest='protein_table',
                        help='Path to the "proteinID_Genome.tsv" file created by '
                             'the CreateSchema process.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='The number of CPU cores to use during the process.')

    args = parser.parse_args()

    header = 'chewBBACA - UniprotFinder'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    input_files = args.input_files
    protein_table = args.protein_table
    cpu_cores = args.cpu_cores

    uniprot_find.main(input_files, protein_table, cpu_cores)


def download_schema():

    def msg(name=None):
        # simple command to download a schema from the NS
        simple_cmd = ('  chewBBACA.py DownloadSchema -sp <species_id> '
                                                  '-sc <schema_id> '
                                                  '-o <download_folder> ')

        # command to download a schema from the NS with non-default arguments values
        params_cmd = ('  chewBBACA.py DownloadSchema -sp <species_id> '
                                                  '-sc <schema_id> '
                                                  '-o <download_folder>\n'
                                                  '\t\t\t      --cpu <cpu_cores> '
                                                  '--ns <nomenclature_server_url> ')

        usage_msg = ('\nDownload schema:\n{0}\n'
                     '\nDownload schema with non-default parameters:\n{1}\n'.format(simple_cmd, params_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='DownloadSchema',
                                     description='This program downloads '
                                                 'a schema from the NS.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('DownloadSchema', nargs='+',
                        help='This program downloads a schema from '
                             'the NS.')

    parser.add_argument('-sp', type=str, required=True,
                        dest='species_id',
                        help='The integer identifier or name of the species '
                             'that the schema is associated to in the NS.')

    parser.add_argument('-sc', type=str, required=True,
                        dest='schema_id',
                        help='The URI, integer identifier or name of '
                             'the schema to download from the NS.')

    parser.add_argument('-o', type=str, required=True,
                        dest='download_folder',
                        help='Output folder to which the schema will '
                             'be saved.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores that will '
                             'be passed to the PrepExternalSchema process to '
                             'determine representatives and create the '
                             'final schema.')

    parser.add_argument('--ns', type=pv.validate_ns_url, required=False,
                        default='main',
                        dest='nomenclature_server',
                        help='The base URL for the Nomenclature Server. '
                             'The default value, "main", will establish a '
                             'connection to "https://chewbbaca.online/", '
                             '"tutorial" to "https://tutorial.chewbbaca.online/" '
                             'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
                             'Users may also provide the IP address to other '
                             'Chewie-NS instances.')

    parser.add_argument('--d', type=str, required=False,
                        default=None,
                        dest='date',
                        help='Download schema with state from specified date. '
                             'Must be in the format "Y-m-dTH:M:S".')

    parser.add_argument('--latest', required=False, action='store_true',
                        dest='latest',
                        help='If the compressed version that is available is '
                             'not the latest, downloads all loci and constructs '
                             'schema locally.')

    args = parser.parse_args()

    header = 'chewBBACA - DownloadSchema'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    ns_species = args.species_id
    ns_schema = args.schema_id
    download_folder = args.download_folder
    cpu_cores = args.cpu_cores
    nomenclature_server = args.nomenclature_server
    date = args.date
    latest = args.latest

    down_schema.main(ns_species, ns_schema, download_folder,
                     cpu_cores, nomenclature_server, date,
                     latest)


def upload_schema():

    def msg(name=None):
        # simple command to load a schema to the NS
        simple_cmd = ('  chewBBACA.py LoadSchema -i <schema_directory> '
                                              '-sp <species_id> '
                                              '-sn <schema_name>\n'
                                              '\t\t\t  -lp <loci_prefix> ')

        # command to load a schema to the NS with non-default arguments values
        params_cmd = ('  chewBBACA.py LoadSchema -i <schema_directory> '
                                              '-sp <species_id> '
                                              '-sn <schema_name>\n'
                                              '\t\t\t  -lp <loci_prefix> '
                                              '--thr <threads> '
                                              '--ns <nomenclature_server_url>')

        # command to continue schema upload that was interrupted or aborted
        continue_cmd = ('  chewBBACA.py LoadSchema -i <schema_directory> '
                                                '-sp <species_id> '
                                                '-sn <schema_name>\n'
                                                '\t\t\t  --continue_up')

        usage_msg = ('\nLoad schema:\n{0}\n'
                     '\nLoad schema with non-default parameters:\n{1}\n'
                     '\nContinue schema upload that was interrupted or aborted:\n{2}\n'.format(simple_cmd, params_cmd, continue_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='LoadSchema',
                                     description='This program uploads '
                                                 'a schema to the NS.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('LoadSchema', nargs='+',
                        help='This program loads a schema to '
                             'the NS.')

    parser.add_argument('-i', type=str, required=True,
                        dest='schema_directory',
                        help='Path to the directory of the schema to upload.')

    parser.add_argument('-sp', type=str, required=True,
                        dest='species_id',
                        help='The integer identifier or name of the species '
                             'that the schema will be associated to in '
                             'the NS.')

    parser.add_argument('-sn', type=str, required=True,
                        dest='schema_name',
                        help='A brief and meaningful name that '
                             'should help understand the type and content '
                             'of the schema.')

    parser.add_argument('-lp', type=str, required=True,
                        dest='loci_prefix',
                        help='Prefix included in the name of each locus of '
                             'the schema.')

    parser.add_argument('--df', type=str, required=False,
                        dest='description_file', default=None,
                        help='Path to a text file with a description '
                             'about the schema. Markdown syntax is supported '
                             'in order to offer greater customizability of '
                             'the rendered description in the Frontend. '
                             'Will default to the schema\'s name if the user '
                             'does not provide a valid path for a file.')

    parser.add_argument('--a', type=str, required=False,
                        dest='annotations', default=None,
                        help='Path to a TSV file with loci annotations. '
                             'The first column has loci identifiers '
                             '(w/o .fasta extension), the second has user '
                             'annotations and the third has custom '
                             'annotations.')

    parser.add_argument('--cpu', type=int, required=False,
                        dest='cpu_cores', default=1,
                        help='Number of CPU cores that will '
                             'be used in the Schema Pre-processing step.')

    parser.add_argument('--thr', type=int, required=False,
                        default=20, dest='threads',
                        help='Number of threads to use to search for '
                             'annotations on UniProt')

    parser.add_argument('--ns', type=pv.validate_ns_url, required=False,
                        default='main',
                        dest='nomenclature_server',
                        help='The base URL for the Nomenclature Server. '
                             'The default value, "main", will establish a '
                             'connection to "https://chewbbaca.online/", '
                             '"tutorial" to "https://tutorial.chewbbaca.online/" '
                             'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
                             'Users may also provide the IP address to other '
                             'Chewie-NS instances.')

    parser.add_argument('--continue_up', required=False, action='store_true',
                        dest='continue_up',
                        help='If the process should check if the schema '
                             'upload was interrupted and try to finish it.')

    args = parser.parse_args()

    header = 'chewBBACA - LoadSchema'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    schema_directory = args.schema_directory
    species_id = args.species_id
    schema_name = args.schema_name
    loci_prefix = args.loci_prefix
    description_file = args.description_file
    annotations = args.annotations
    cpu_cores = args.cpu_cores
    threads = args.threads
    nomenclature_server = args.nomenclature_server
    continue_up = args.continue_up

    load_schema.main(schema_directory, species_id, schema_name,
                     loci_prefix, description_file, annotations,
                     cpu_cores, threads, nomenclature_server,
                     continue_up)


def synchronize_schema():

    def msg(name=None):
        # simple command to synchronize a schema with its NS version
        simple_cmd = ('  chewBBACA.py SyncSchema -sc <schema_directory> ')

        # command to synchronize a schema with its NS version with non-default arguments values
        params_cmd = ('  chewBBACA.py SyncSchema -sc <schema_directory> '
                                                '--cpu <cpu_cores> '
                                                '--ns <nomenclature_server_url>')

        # command to submit novel local alleles
        submit_cmd = ('  chewBBACA.py SyncSchema -sc <schema_directory> --submit')

        usage_msg = ('\nSync schema:\n{0}\n'
                     '\nSync schema with non-default parameters:\n{1}\n'
                     '\nSync schema and send novel local alleles to the NS:\n{2}\n'.format(simple_cmd, params_cmd, submit_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='SyncSchema',
                                     description='This program syncs a local '
                                                 'schema with NS',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('SyncSchema', nargs='+',
                        help='Synchronize a local schema, previously '
                             'downloaded from the NS, with its latest '
                             'version in the NS.')

    parser.add_argument('-sc', type=str, required=True,
                        dest='schema_directory',
                        help='Path to the directory with the schema to be'
                             'synced.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores that will '
                             'be used to determine new representatives '
                             'if the process downloads new alleles from '
                             'the Chewie-NS.')

    parser.add_argument('--ns', type=pv.validate_ns_url, required=False,
                        default=None,
                        dest='nomenclature_server',
                        help='The base URL for the Nomenclature Server. '
                             'The default option will get the base URL from the '
                             'schema\'s URI. It is also possible to specify other '
                             'options that are available in chewBBACA\'s configs, '
                             'such as: "main" will establish a connection to '
                             '"https://chewbbaca.online/", "tutorial" to '
                             '"https://tutorial.chewbbaca.online/" and "local" '
                             'to "http://127.0.0.1:5000/NS/api/" (localhost). '
                             'Users may also provide the IP address to other '
                             'Chewie-NS instances.')

    parser.add_argument('--submit', required=False,
                        action='store_true', dest='submit',
                        help='If the process should identify new alleles '
                             'in the local schema and send them to the '
                             'NS. (only users with permissons level of '
                             'Contributor can submit new alleles).')

    args = parser.parse_args()

    header = 'chewBBACA - SyncSchema'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    schema_directory = args.schema_directory
    cpu_cores = args.cpu_cores
    nomenclature_server = args.nomenclature_server
    submit = args.submit

    sync_schema.main(schema_directory, cpu_cores,
                     nomenclature_server, submit)


def ns_stats():

    def msg(name=None):
        # simple command to list species and totals
        simple_cmd = ('  chewBBACA.py NSStats -m species ')

        # command to list all schemas for a species
        schemas_cmd = ('  chewBBACA.py NSStats -m schemas --sp <species_id> ')

        # command to get information about a single schema
        schema_cmd = ('  chewBBACA.py NSStats -m schemas --sp <species_id> '
                                             '--sc <schema_id>')

        usage_msg = ('\nList species and totals:\n{0}\n'
                     '\nList all schemas for a species and associated information:\n{1}\n'
                     '\nGet information about a particular schema:\n{2}\n'
                     ''.format(simple_cmd, schemas_cmd, schema_cmd))

        return usage_msg

    parser = argparse.ArgumentParser(prog='NSStats',
                                     description='Retrieve basic information '
                                                 'about the species and schemas in '
                                                 'the Chewie-NS.',
                                     usage=msg(),
                                     formatter_class=ModifiedHelpFormatter)

    parser.add_argument('NSStats', nargs='+',
                        help='')

    parser.add_argument('-m', type=str, required=True,
                        dest='stats_mode', choices=['species', 'schemas'],
                        help='The process can retrieve the list of species '
                             '("species" option) in the Chewie-NS or the '
                             'list of schemas for a species '
                             '("schemas" option).')

    parser.add_argument('--ns', type=pv.validate_ns_url, required=False,
                        default='main',
                        dest='nomenclature_server',
                        help='The base URL for the Nomenclature Server. '
                             'The default value, "main", will establish a '
                             'connection to "https://chewbbaca.online/", '
                             '"tutorial" to "https://tutorial.chewbbaca.online/" '
                             'and "local" to "http://127.0.0.1:5000/NS/api/" (localhost). '
                             'Users may also provide the IP address to other '
                             'Chewie-NS instances.')

    parser.add_argument('--sp', type=str, required=False,
                        dest='species_id', default=None,
                        help='The integer identifier of a '
                             'species in the Chewie-NS.')

    parser.add_argument('--sc', type=str, required=False,
                        dest='schema_id', default=None,
                        help='The integer identifier of a schema in '
                             'the Chewie-NS.')

    args = parser.parse_args()

    header = 'chewBBACA - NSStats'
    hf = '='*(len(header)+4)
    print('{0}\n  {1}\n{0}'.format(hf, header, hf))

    mode = args.stats_mode
    nomenclature_server = args.nomenclature_server
    species_id = args.species_id
    schema_id = args.schema_id

    stats_requests.main(mode, nomenclature_server, species_id, schema_id)


def main():

    functions_info = {'CreateSchema': ['Create a gene-by-gene schema based on '
                                       'a set of input genomes.',
                                       create_schema],
                      'AlleleCall': ['Determine the allelic profiles of a set of '
                                     'input genomes based on a schema.',
                                     allele_call],
                      'SchemaEvaluator': ['Tool that builds an html output '
                                          'to better navigate/visualize '
                                          'your schema.',
                                          evaluate_schema],
                      'TestGenomeQuality': ['Analyze your allele call output '
                                            'to refine schemas.',
                                            test_schema],
                      'ExtractCgMLST': ['Determines the set of '
                                        'loci that constitute the '
                                        'core genome based on a '
                                        'threshold.',
                                        extract_cgmlst],
                      'RemoveGenes': ['Remove a provided list of loci from '
                                      'your allele call output.',
                                      remove_genes],
                      'PrepExternalSchema': ['Adapt an external schema to be '
                                             'used with chewBBACA.',
                                             prep_schema],
                      'JoinProfiles': ['Join two profiles in a single profile '
                                       'file.',
                                       join_profiles],
                      'UniprotFinder': ['Retrieve annotations for loci in a schema.',
                                        find_uniprot],
                      'DownloadSchema': ['Download a schema from the Chewie-NS.',
                                         download_schema],
                      'LoadSchema': ['Upload a schema to the Chewie-NS.',
                                     upload_schema],
                      'SyncSchema': ['Synchronize a schema with its remote version '
                                     'in the Chewie-NS.',
                                     synchronize_schema],
                      'NSStats': ['Retrieve basic information about the species '
                                  'and schemas in the Chewie-NS.',
                                  ns_stats]}

    authors = 'Mickael Silva, Pedro Cerqueira, Rafael Mamede'
    repository = 'https://github.com/B-UMMI/chewBBACA'
    wiki = 'https://github.com/B-UMMI/chewBBACA/wiki'
    tutorial = 'https://github.com/B-UMMI/chewBBACA_tutorial'
    contacts = 'imm-bioinfo@medicina.ulisboa.pt'

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
