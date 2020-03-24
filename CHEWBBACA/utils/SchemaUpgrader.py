#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

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
import argparse

import auxiliary_functions as aux


def main(schema_directory, bsr, prodigal_training_file, translation_table,
         minimum_locus_length, size_threshold, chewBBACA_version):

    params = {'bsr': bsr,
              'prodigal_training_file': prodigal_training_file,
              'translation_table': translation_table,
              'minimum_locus_length': minimum_locus_length,
              'chewBBACA_version': chewBBACA_version,
              'size_threshold': size_threshold}

    print('\nYou have provided the following parameters values:')

    # display parameter values provided by the user
    params_lines = []
    params_header = '\n{:<25} {:<100}'.format('Parameter', 'Value')
    params_lines.append(params_header)
    params_lines += ['{:<25} {:<100}'.format(k, v) for k, v in params.items()]
    params_table = '\n'.join(params_lines)

    print(params_table)

    proc = input('\nDo you wish to proceed with provided parameters values?\n')

    if proc.lower() not in ['y', 'yes']:
        sys.exit('Exited.')

    # copy Prodigal training file to schema directory
    ptf_name = os.path.basename(params['prodigal_training_file'])
    ptf_path = os.path.join(schema_directory, ptf_name)
    if os.path.isfile(ptf_path) is False:
        shutil.copy(params['prodigal_training_file'], ptf_path)

    # create file with schema configs
    # determine PTF checksum
    ptf_hash = aux.binary_file_hash(ptf_path)
    params['prodigal_training_file'] = ptf_hash
    config_file = os.path.join(schema_directory, '.schema_config')
    with open(config_file, 'wb') as cf:
        pickle.dump(params, cf)

    # create hidden file with genes/loci list
    genes_list_file = aux.write_gene_list(schema_directory)


def parse_arguments():


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='schema_directory',
                        help='')

    parser.add_argument('--bsr', type=float, required=True,
                        dest='blast_score_ratio',
                        help='')

    parser.add_argument('--ptf', type=str, required=True,
                        dest='prodigal_training_file',
                        help='')

    parser.add_argument('--tt', type=int, required=True,
                        dest='translation_table',
                        help='')

    parser.add_argument('--ml', type=int, required=True,
                        dest='minimum_locus_length',
                        help='')

    parser.add_argument('--st', type=float, required=True,
                        dest='size_threshold',
                        help='')

    parser.add_argument('--cv', type=str, required=False,
                        dest='chewBBACA_version',
                        default='2.1.0',
                        help='')

    args = parser.parse_args()

    return [args.schema_directory, args.blast_score_ratio,
            args.prodigal_training_file, args.translation_table,
            args.minimum_locus_length, args.size_threshold,
            args.chewBBACA_version]


if __name__ == '__main__':

    arguments = parse_arguments()

    main(arguments[0], arguments[1], arguments[2],
         arguments[3], arguments[4], arguments[5],
         arguments[6])
