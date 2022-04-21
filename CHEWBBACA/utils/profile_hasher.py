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
import hashlib
import argparse

import pandas as pd
from Bio import SeqIO


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


def main(profiles_table, schema_directory, output_file, hash_type):

    # get hash function
    try:
        hashing_function = getattr(hashlib, hash_type)
    except Exception:
        hashing_function = getattr(zlib, hash_type)
    except Exception:
        sys.exit('{0} hash function is not available in hashlib or zlib.'.format(hash_function))

    with open(profiles_table, 'r') as infile:
        loci_ids = infile.readlines(1)[0].split()[1:]

    # get column with sample identifiers
    sample_ids = pd.read_csv(profiles_table, usecols=['FILE'], delimiter='\t', dtype=str)['FILE']
    hashed_profiles = [sample_ids]
    # process by column chunks to avoid high memory usage with huge files
    # read row chunks instead so that it is not needed to have the full hashed matrix in memory?
    for i in range(0, len(loci_ids), 250):
        current_loci = loci_ids[i:i+250]
        loci_columns = pd.read_csv(profiles_table, usecols=current_loci, delimiter='\t', dtype=str)
        for locus in current_loci:
            locus_file = os.path.join(schema_directory, locus)
            if locus_file.endswith('.fasta') is False:
                locus_file += '.fasta'
            
            locus_column = loci_columns[locus]
            hashed_column = hash_column(locus_column, locus_file, hashing_function)
            hashed_profiles.append(hashed_column)
        
        print('Processed {0}'.format(i+250))

    hashed_df = pd.concat(hashed_profiles, axis=1)
    hashed_df.to_csv(output_file, sep='\t', index=False)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--profiles', type=str, required=True,
                        dest='profiles_table',
                        help='')

    parser.add_argument('-s', '--schema-directory', type=str, required=True,
                        dest='schema_directory',
                        help='')

    parser.add_argument('-o', '--output-file', type=str, required=True,
                        dest='output_file',
                        help='')

    parser.add_argument('-hf', '--hash-type', type=str, required=False,
                        dest='hash_type',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
