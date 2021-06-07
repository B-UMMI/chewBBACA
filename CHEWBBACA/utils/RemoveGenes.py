#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module 

Expected input
--------------

The process expects the following variables whether through command
line execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_file`` : 

    - e.g.: ``/home/user/results_alleles.tsv``

- ``-gl``, ``genes_list`` : 

    - e.g.: ``/home/user/genes_list.txt``

- ``-o``, ``output_file`` : 

    - e.g.: ``/home/user/selected_results.tsv``

- ``--inverse`` : 

Code documentation
------------------
"""


import csv
import argparse
import pandas as pd


def get_headers(files, delimiter='\t'):
    """ Gets the headers (first line) from a set of
        files.

        Parameters
        ----------
        files : list
            List with paths to files.

        Returns
        -------
        headers : list
            List with the first line in each file
            (each header is a sublist of elements
            separated based on the delimiter).
    """

    headers = []
    for file in files:
        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter=delimiter)
            header = next(reader)
            headers.append(header)

    return headers


#input_file = '/home/rfm/Lab_Analyses/GAS_PrepExternalSchema/wgMLST_schema/spyogenes_schema/AlleleCall_final_rounds/05_06_2021/spyogenes_ena_results/results_alleles_ena1.tsv'
#genes_list = '/home/rfm/Lab_Analyses/GAS_PrepExternalSchema/wgMLST_schema/spyogenes_schema/AlleleCall_final_rounds/05_06_2021/spyogenes_ena_results/rm_genes.txt'
#output_file = '/home/rfm/Lab_Analyses/GAS_PrepExternalSchema/wgMLST_schema/spyogenes_schema/AlleleCall_final_rounds/05_06_2021/spyogenes_ena_results/rm_matrix.tsv'
#inverse = 'False'
def main(input_file, genes_list, output_file, inverse):

    # read genes list
    with open(genes_list, 'r') as infile:
        genes_list = list(csv.reader(infile, delimiter='\t'))
        genes_list = [g[0] for g in genes_list]

    # get list of loci in allele call results
    loci = get_headers([input_file])[0]
    print('Total loci: {0}'.format(len(loci)-1))

    if inverse is True:
        columns_to_keep = [g for g in loci if g in genes_list]
    else:
        columns_to_keep = [g for g in loci if g not in genes_list]

    columns_to_remove = (len(loci)-1) - len(columns_to_keep)
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
                        help='')

    parser.add_argument('-gl', '--genes-list', type=str,
                        required=True, dest='genes_list',
                        help='')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='')

    parser.add_argument('--inverse', action='store_true',
                        required=False, dest='inverse',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
