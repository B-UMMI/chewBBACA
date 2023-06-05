#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module removes a set of loci from a TSV file with
results from the AlleleCall process.

Code documentation
------------------
"""


import csv
import argparse

import pandas as pd

try:
    from utils import file_operations as fo
except ModuleNotFoundError:
    from CHEWBBACA.utils import file_operations as fo


def main(input_file, genes_list, output_file, inverse):

    # read genes list
    with open(genes_list, 'r') as infile:
        genes_list = list(csv.reader(infile, delimiter='\t'))
        genes_list = [g[0] for g in genes_list[1:]]

    # get list of loci in allele call results
    loci = fo.read_lines(input_file, strip=True, num_lines=1)
    loci = loci[0].split('\t')[1:]
    print('Total loci: {0}'.format(len(loci)))

    if inverse is True:
        columns_to_keep = [g for g in loci if g in genes_list]
    else:
        columns_to_keep = [g for g in loci if g not in genes_list]

    columns_to_remove = (len(loci)) - len(columns_to_keep)
    print('Loci to remove: {0}'.format(columns_to_remove))

    # include first column with sample ids
    columns_to_keep = ['FILE'] + columns_to_keep
    df = pd.read_csv(input_file, usecols=columns_to_keep,
                     sep='\t', dtype=str)

    # save dataframe to file
    df.to_csv(output_file, header=True, sep='\t', index=False)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_file',
                        help='TSV file that contains a matrix with '
                             'allelic profiles determined by the '
                             'AlleleCall process.')

    parser.add_argument('-gl', '--genes-list', type=str,
                        required=True, dest='genes_list',
                        help='File with the list of genes to '
                             'remove, one identifier per line.')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='Path to the output file.')

    parser.add_argument('--inverse', action='store_true',
                        required=False, dest='inverse',
                        help='List of genes that is provided '
                             'is the list of genes to keep and '
                             'all other genes should be removed.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
