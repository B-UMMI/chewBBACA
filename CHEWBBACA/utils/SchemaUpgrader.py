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


import argparse


def main(schema_directory, bsr, chewBBACA_version, prodigal_training_file,
         translation_table, minimum_locus_length, size_threshold, word_size,
         cluster_sim, representative_filter, intraCluster_filter):






def parse_arguments():


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True, dest='schema_directory',
                        help='')

    parser.add_argument('--bsr', type=float, required=True, dest='blast_score_ratio',
                        help='')

    parser.add_argument('--ptf', type=str, required=True, dest='prodigal_training_file',
                        help='')

    parser.add_argument('--tt', type=int, required=True, dest='translation_table',
                        help='')

    parser.add_argument('--ml', type=int, required=True, dest='minimum_locus_length',
                        help='')

    parser.add_argument('--st', type=float, required=True, dest='size_threshold',
                        help='')
    
    parser.add_argument('--cv', type=str, required=False, dest='chewBBACA_version',
                        default='2.1.0',
                        help='')

    args = parser.parse_args()

    return [args.stats_mode, args.base_url, args.taxon]


