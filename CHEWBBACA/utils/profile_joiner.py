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


import csv
import sys
import shutil
import argparse


def count_lines(file):
    """ Counts the number of lines in a file.

        Parameters
        ----------
        file : str
            Path to a file.

        Returns
        -------
        num_lines : int
            Total number of lines in the file.
    """

    with open(file, 'rb') as infile:
        num_lines = sum(1 for line in infile)

    return num_lines


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


def concatenate_profiles(profiles, output_file):
    """ Concatenates TSV files with allele calling results
        for the same set of loci.

        Parameters
        ----------
        profiles : list
            List with the paths to the TSV files with
            allele calling results.
        output_file : str
            Path to the output file.

        Returns
        -------
        total_profiles : int
            Number of profiles written to the output file.
    """

    total_profiles = 0
    with open(output_file, 'w') as outfile:
        for i, file in enumerate(profiles):
            with open(file, 'r') as infile:
                # skip header if not first file
                if i != 0:
                    infile.readline()
                # copy file content to output file
                shutil.copyfileobj(infile, outfile)
                # count number of profiles (excluding header)
                total_profiles += count_lines(file) - 1

    return total_profiles


def concatenate_common_loci(profiles, common_loci, output_file):
    """ Creates a file with the allele calling results for a set
        of loci common to several profiles.

        Parameters
        ----------
        profiles : list
            List with the paths to the TSV files with
            allele calling results.
        common_loci : list
            List with loci identifiers that are common
            between all profiles.
        output_file : str
            Path to the output file.

        Returns
        -------
        total_profiles : int
            Number of profiles written to the output file.
    """

    total_profiles = 0
    with open(output_file, 'w') as outfile:
        # write header
        outfile.write('\t'.join(common_loci)+'\n')
        for file in profiles:
            with open(file, 'r') as infile:
                lines = list(csv.reader(infile, delimiter='\t'))
                header = lines[0]
                profiles = lines[1:]
                total_profiles += len(profiles)
                # get loci indexes
                loci_indexes = [header.index(l) for l in common_loci]
                # get all profiles
                common_profiles = [[p[l] for l in loci_indexes]
                                   for p in profiles]
                # write profiles to file
                tsv_profiles = ['\t'.join(p) for p in common_profiles]
                text_profiles = '\n'.join(tsv_profiles)
                outfile.write(text_profiles+'\n')

    return total_profiles


def main(profiles, output_file, common):

    if len(profiles) == 1:
        sys.exit('Provided a single file. Nothing to do.')

    headers = get_headers(profiles)

    if common is False:
        # check if headers are equal
        if all([set(headers[0]) == set(h) for h in headers[1:]]) is True:
            print('Profiles have {0} loci.'.format(len(headers[0])-1))
            total_profiles = concatenate_profiles(profiles, output_file)
        else:
            sys.exit('Files have different sets of loci. Please provide '
                     'files with the results for the same set of loci or '
                     'provide the "--common" parameter to create a file '
                     'with the results for the set of common loci.')
    else:
        # determine set of common loci
        common_loci = set(headers[0])
        for h in headers[1:]:
            # determine common loci ordered
            latest_common = [l for l in h if l in common_loci]
            # update set of common loci so far
            common_loci = set(latest_common)

        if len(latest_common) <= 1:
            sys.exit('Profiles do not have loci in common.')

        print('Profiles have {0} loci in common.'
              ''.format(len(latest_common)-1))
        total_profiles = concatenate_common_loci(profiles,
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
