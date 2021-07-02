#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module joins allele calling results from different
runs. It can concatenate files with allelic profiles
for the same set of loci or create a new file with the
allelic profiles for the loci that were common between
all input files.

Expected input
--------------

The process expects the following variables whether through command
line execution or invocation of the :py:func:`main` function:

- ``-p``, ``profiles`` : Path to files containing the results from
  the AlleleCall process (results_alleles.tsv).

    - e.g.: ``/home/user/profile1.tsv /home/user/profile2.tsv``

- ``-o``, ``output_file`` : Path to the output file.

    - e.g.: ``/home/user/joined_profiles.tsv``

- ``--common`` : Create file with profiles for the set of
  common loci.

Code documentation
------------------
"""


import os
import sys
import argparse

import pandas as pd

try:
    from utils import file_operations as fo
except:
    from CHEWBBACA.utils import file_operations as fo


def concatenate_profiles(files, loci_list, output_file):
    """ Concatenates TSV files with allele calling results
        for a common set of loci.

        Parameters
        ----------
        files : list
            List with the paths to the TSV files with
            allele calling results.
        loci_list : list
            List with the identifiers of common loci.
        output_file : str
            Path to the output file.

        Returns
        -------
        total_profiles : int
            Number of profiles written to the output file.
    """

    total_profiles = 0
    for f in files:
        # create dataframe with columns for loci in list
        profiles = pd.read_csv(f, sep='\t',
                               usecols=loci_list, dtype=str)
        # reorder columns to prevent concatenating results
        # with same loci but with different loci/column order
        profiles = profiles[loci_list]
        total_profiles += len(profiles)
        # save dataframe to file
        profiles.to_csv(output_file, mode='a',
                        header=not os.path.exists(output_file),
                        sep='\t', index=False)

    return total_profiles


def main(profiles, output_file, common):

    if len(profiles) == 1:
        sys.exit('Provided a single file. Nothing to do.')

    headers = fo.get_headers(profiles)

    if common is False:
        # check if headers are equal
        if all([set(headers[0]) == set(h) for h in headers[1:]]) is True:
            print('Profiles have {0} loci.'.format(len(headers[0])-1))

            total_profiles = concatenate_profiles(profiles,
                                                  headers[0],
                                                  output_file)
        else:
            sys.exit('Files have different sets of loci. Please provide '
                     'files with the results for the same set of loci or '
                     'provide the "--common" parameter to create a file '
                     'with the results for the set of common loci.')
    else:
        # determine set of common loci
        common_loci = headers[0]
        for h in headers[1:]:
            # determine common loci ordered based on headers in first file
            latest_common = [l for l in common_loci if l in h]
            # update set of common loci so far
            common_loci = latest_common

        if len(latest_common) <= 1:
            sys.exit('Profiles do not have loci in common.')

        print('Profiles have {0} loci in common.'
              ''.format(len(latest_common)-1))
        total_profiles = concatenate_profiles(profiles,
                                              latest_common,
                                              output_file)

    print('Joined {0} files with a total of {1} '
          'profiles.'.format(len(profiles), total_profiles))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-p', '--profiles', nargs='+', type=str,
                        required=True, dest='profiles',
                        help='Path to files containing the results from '
                             'the AlleleCall process (results_alleles.tsv).')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='Path to the output file.')

    parser.add_argument('--common', action='store_true',
                        required=False, dest='common',
                        help='Create file with profiles for '
                             'the set of common loci.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
