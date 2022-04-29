#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module enables

Expected input
--------------


Code documentation
------------------
"""


import os
import sys
import zlib
import shutil
import hashlib
import argparse

import pandas as pd
from Bio import SeqIO

try:
    from utils import (file_operations as fo,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (file_operations as fo,
                                 multiprocessing_operations as mo)


def hash_column(column, locus_file, hashing_function):
    """
    """

    # read Fasta file with locus alleles
    locus_alleles = {(rec.id).split('_')[-1]: str(rec.seq)
                     for rec in SeqIO.parse(locus_file, 'fasta')}
    
    hashed_alleles = {}
    for seqid, seq in locus_alleles.items():
        # hash function does not accept string object, encode to get bytes object
        hashed_seq = hashing_function(seq.encode())
        # bitwise operation to convert crc32 and adler32 hashes to unsigned
        # integer and ensure the computed value is the same for Python 2 & 3
        if isinstance(hashed_seq, int):
            hashed_seq &= 0xffffffff
        hashed_alleles[seqid] = hashed_seq

    # faster than replace or map with update to avoid adding NaN
    hashed_column = column.apply(lambda x: hashed_alleles.get(x, x))

    return hashed_column


def hash_profiles(profiles_table, loci_ids, loci_files, hashing_function,
                  nrows, skiprows, include_header, output_directory):
    """
    """

    current_rows = pd.read_csv(profiles_table, delimiter='\t', dtype=str,
                               skiprows=skiprows, nrows=nrows)

    current_samples = current_rows['FILE']
    hashed_profiles = [current_samples]
    for locus in loci_ids:
        locus_column = current_rows[locus]
        hashed_column = hash_column(locus_column, loci_files[locus], hashing_function)
        hashed_profiles.append(hashed_column)

    hashed_df = pd.concat(hashed_profiles, axis=1)
    start = skiprows.stop
    stop = skiprows.stop+nrows
    output_file = fo.join_paths(output_directory, ['df_{0}-{1}.tsv'.format(start, stop)])
    hashed_df.to_csv(output_file, sep='\t', index=False, header=include_header)

    return output_file
    

# profiles_table = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/small_dataset.tsv'
# schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/Streptococcus_pyogenes_wgMLST_schema/'
# output_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/'
# hash_type = 'crc32'
# threads = 2
# nrows = 1000
def main(profiles_table, schema_directory, output_directory, hash_type, threads, nrows):

    # get hash function
    try:
        hashing_function = getattr(hashlib, hash_type)
    except Exception:
        hashing_function = getattr(zlib, hash_type)
    except Exception:
        sys.exit('{0} hash function is not available in hashlib or zlib.'.format(hash_function))

    # get loci identifiers
    with open(profiles_table, 'r') as infile:
        header = infile.readlines(1)[0].split()
        loci_ids = header[1:]

    loci_files = {}
    for locus in loci_ids:
        locus_file = fo.join_paths(schema_directory, [locus])
        if locus_file.endswith('.fasta') is False:
            locus_file += '.fasta'
        loci_files[locus] = locus_file

    sample_ids = pd.read_csv(profiles_table, delimiter='\t', dtype=str, usecols=['FILE'])

    # create multiprocessing inputs
    multi_inputs = []
    include_header = True
    # divide and process by row chunks
    for i in range(0, len(sample_ids), nrows):
        multi_inputs.append([profiles_table, loci_ids, loci_files,
                             hashing_function, nrows, range(1, i+1),
                             include_header, output_directory, hash_profiles])
        include_header = False

    output_files = mo.map_async_parallelizer(multi_inputs, mo.function_helper, threads)

    # concatenate all files
    output_file = fo.join_paths(output_directory, ['final_df.tsv'])
    fo.concatenate_files(output_files, output_file)

    # delete intermediate dataframes
    fo.remove_files(output_files)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--profiles', type=str, required=True,
                        dest='profiles_table',
                        help='')

    parser.add_argument('-s', '--schema-directory', type=str, required=True,
                        dest='schema_directory',
                        help='')

    parser.add_argument('-o', '--output-directory', type=str, required=True,
                        dest='output_directory',
                        help='')

    parser.add_argument('-hf', '--hash-type', type=str, required=False,
                        dest='hash_type',
                        help='')

    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=1, dest='threads',
                        help='')

    parser.add_argument('-n', '--nrows', type=int, required=False,
                        default=1000, dest='nrows',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
