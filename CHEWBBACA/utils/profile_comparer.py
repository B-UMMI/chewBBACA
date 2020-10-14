#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 19:08:47 2020

@author: rfm
"""


import os
import csv
import pickle
import hashlib
import argparse

from Bio import SeqIO


def pickle_loader(input_file):
    """ Use the Pickle module to de-serialize an object.

        Parameters
        ----------
        input_file : str
            Path to file with byte stream to be
            de-serialized.

        Returns
        -------
        data : type
            Variable that refers to the de-serialized
            object.
    """

    with open(input_file, 'rb') as pi:
        data = pickle.load(pi)

    return data


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Args:
            fasta_path (str): full path to the FASTA file.

        Returns:
            dictionary that has sequences ids as keys and DNA
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {(rec.id).split('_')[-1]: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def sequence_kmerizer(sequence, k_value, offset=1, position=False):
    """
    """

    if position is False:
        kmers = [sequence[i:i+k_value] for i in range(0, len(sequence)-k_value+1, offset)]
    elif position is True:
        kmers = [(sequence[i:i+k_value], i) for i in range(0, len(sequence)-k_value+1, offset)]

    return kmers


def mask_missing(alleles_sets):
    """
    """

    current_sets_miss = []
    for s in alleles_sets:
        try:
            new_set = [int(j) for j in s[0]]
            current_sets_miss.append([new_set, s[1]])
        except:
            current_sets_miss.append(None)

    return current_sets_miss


profiles_path = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/Annotations/profiles_test_masked.tsv'
schema_path = '/home/rfm/Desktop/rfm/Lab_Analyses/Bacgentrack_data/sagalactiae_schema'
k_value = 3
output_file = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/Annotations/similarity_matrix'



def matrix_hashes(profiles_path, schema_path, output_file):

    # import profiles
    with open(profiles_path, 'r') as p:
        lines = list(csv.reader(p, delimiter='\t'))

    loci = lines[0][1:]
    samples = [l[0] for l in lines[1:]]
    profiles = [l[1:] for l in lines[1:]]

    # divide profiles into sets of k alleles
    allelic_sets = [sequence_kmerizer(p, k_value, position=True) for p in profiles]

    # mask sets with missing data
    masked_sets = [mask_missing(s) for s in allelic_sets]

    # groups sets from same allelic sets
    # vary window to group several sets and decide which are lexicographically smaller???
    grouped_sets = []
    for i in range(len(masked_sets[0])):
        grouped_sets.append([s[i] for s in masked_sets])

    # determine_profiles hashes
    count = 0
    hashes = []
    for g in grouped_sets:
        # get identifiers of loci in set
        current_loci = loci[count:count+k_value]
        count += 1

        # get fasta records for loci set
        loci_seqs = []
        for l in current_loci:
            file_path = os.path.join(schema_path, l)
            loci_seqs.append(import_sequences(file_path))

        # get allele sequences
        alleles = []
        for allele_set in g:
            if allele_set is not None:
                alleles.append([loci_seqs[i][str(a)] for i, a in enumerate(allele_set[0])])
            else:
                alleles.append(None)

        # determine hashes
        alleles_hashes = []
        for allele_set in alleles:
            if allele_set is not None:
                concat = ''.join(allele_set)
                set_hash = hashlib.sha256(concat.encode('utf-8')).hexdigest()
                alleles_hashes.append(set_hash)
            else:
                alleles_hashes.append(None)

        hashes.append(alleles_hashes)

    # regroup hashes from same profiles
    profiles_hashes = []
    for s in range(len(samples)):
        current_profile = [hashes[i][s] for i in range(len(hashes))]
        profiles_hashes.append([samples[s], current_profile])

    # save results
    with open(output_file, 'wb') as outfile:
        pickle.dump(profiles_hashes, outfile)

    print('Determined hashes for allele sets of {0} '
          'profiles.'.format(len(profiles_hashes)))


profiles1 = '/home/rfm/Desktop/rfm/Lab_Analyses/GAS_PrepExternalSchema/Annotations/sim_file'
profiles2 = profiles1

def sims_matrix(profiles1, profiles2):
    """
    """

    profiles1_hashes = pickle_loader(profiles1)
    profiles2_hashes = pickle_loader(profiles2)

    hashes_results = []
    # need to remove NoneType from hashes list
    for p in profiles1_hashes:
        # determine
        sims = [len(set(p[1]).intersection(k[1]))/len(set(p[1])) for k in profiles2_hashes]
        sims = [round(s, 2) for s in sims]
        sims2 = [len(set(p[1]).intersection(k[1]))/len(set(k[1])) for k in profiles2_hashes]
        sims2 = [round(s, 2) for s in sims2]
        final_sims = list(zip(sims, sims2))
        hashes_results.append([p[0], final_sims])

    header = ['FILE'] + [n[0] for n in profiles2_hashes]
    header = '\t'.join(header)
    outlines = [header] + ['{0}\t{1}'.format(k[0], '\t'.join([str(n) for n in k[1]])) for k in hashes_results]

    with open(output_file, 'w') as out:
        out.write('\n'.join(outlines))





def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', type=str, required=True,
                        choices=['hash', 'comp'],
                        dest='mode',
                        help='')

    parser.add_argument('-p', nargs='?', type=str, required=True,
                        dest='profiles_path',
                        help='')

    parser.add_argument('-q', nargs='?', type=str, required=False,
                        dest='comp_path',
                        help='')

    parser.add_argument('-s', type=str, required=True,
                        dest='schema_path')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_file')

    args = parser.parse_args()

    return [args.mode, args.profiles_path, args.comp_path,
            args.schema_path, args.output_file]


if __name__ == '__main__':

    args = parse_arguments()

    if args[0] == 'hash':
        matrix_hashes(args[1], args[3], args[4])
    elif args[0] == 'comp':
        sims_matrix(args[1], args[2], args[4])
