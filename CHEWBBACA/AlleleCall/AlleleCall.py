#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables ...

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_files`` : Path to the directory that contains the input
  FASTA files. Alternatively, a single file with a list of paths to FASTA
  files, one per line.

    - e.g.: ``/home/user/genomes``

- ``-o``, ``output_directory`` : Output directory where the process will
  store intermediate files and create the schema's directory.

    - e.g.: ``/home/user/schemas/new_schema``

- ``--ptf``, ``ptf_path`` : Path to the Prodigal training file.

    - e.g.: ``/home/user/training_files/species.trn``

- ``--bsr``, ``blast_score_ratio`` : BLAST Score Ratio value.

    - e.g.: ``0.6``

- ``--l``, ``minimum_length`` : Minimum sequence length. Coding sequences
  shorter than this value are excluded.

    - e.g.: ``201``

- ``--t``, ``translation_table`` : Genetic code used to predict genes and
  to translate coding sequences.

    - e.g.: ``11``

- ``--st``, ``size_threshold`` : CDS size variation threshold. Added to the
  schema's config file and used to identify alleles with a length value that
  deviates from the locus length mode during the allele calling process.

    - e.g.: ``0.2``

- ``--w``, ``word_size`` : word size used to generate k-mers during the
  clustering step.

    - e.g.: ``5``

- ``--ws``, ``window_size`` : window size value. Number of consecutive
  k-mers included in each window to determine a minimizer.

    - e.g.: ``5``

- ``--cs``, ``clustering_sim`` : clustering similarity threshold. Minimum
  decimal proportion of shared distinct minimizers for a sequence to be
  added to a cluster.

    - e.g.: ``0.2``

- ``--rf``, ``representative_filter`` : representative similarity threshold.
  Clustered sequences are excluded if they share this proportion of distinct
  minimizers with the cluster representative.

    - e.g.: ``0.9``

- ``--if``, ``intra_filter`` : intra-cluster similarity threshold. Clustered
  sequences are excluded if they share this proportion of distinct minimizers
  with another clustered sequence of equal or greater length.

    - e.g.: ``0.9``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores used to run the process.

    - e.g.: ``4``

- ``--b``, ``blast_path`` : Path to the BLAST executables.

    - e.g.: ``/home/software/blast``

- ``--pm``, ``prodigal_mode`` : Prodigal running mode.

    - e.g.: ``single``

- ``--CDS``, ``cds_input`` : If provided, input is a single or several FASTA
  files with coding sequences (skips gene prediction and CDS extraction).

    - e.g.: ``/home/user/coding_sequences_files``

- ``--no-cleanup``, ``no_cleanup`` : If provided, intermediate files
  generated during process execution are not removed at the end.

Code documentation
------------------
"""


import os
import sys
import math
import argparse

from Bio import SeqIO

try:
    from utils import (constants as ct,
                       blast_wrapper as bw,
                       core_functions as cf,
                       gene_prediction as gp,
                       file_operations as fo,
                       fasta_operations as fao,
                       sequence_clustering as sc,
                       sequence_manipulation as sm,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 blast_wrapper as bw,
                                 core_functions as cf,
                                 gene_prediction as gp,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_clustering as sc,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


def exact_matches(fasta_file, hashes):
    """
    """

    # import fasta records
    records = fao.import_sequences(fasta_file)

    # translate records
    translated_records = set([str(sm.translate_sequence(seq, 11))
                              for seq in records.values()])

    records_hashes = [im.hash_sequence(seq)
                      for seq in translated_records]

    matches = [h for h in hashes if h in records_hashes]

    return matches


input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids320.txt'
output_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/test_allelecall'
ptf_path = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/Streptococcus_agalactiae.trn'
blast_score_ratio = 0.6
minimum_length = 201
translation_table = 11
size_threshold = 0.2
word_size = 5
window_size = 5
clustering_sim = 0.2
representative_filter = 0.9
intra_filter = 0.9
cpu_cores = 6
blast_path = '/home/rfm/Software/anaconda3/envs/ns/bin'
prodigal_mode = 'single'
cds_input = False
schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed'
def allele_calling(input_files, schema_directory, output_directory, ptf_path,
                   blast_score_ratio, minimum_length, translation_table,
                   size_threshold, word_size, window_size, clustering_sim,
                   representative_filter, intra_filter, cpu_cores, blast_path,
                   prodigal_mode, cds_input):
    """
    """

    # define directory for temporary files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # read file with paths to input files
    fasta_files = fo.read_lines(input_files, strip=True)

    # sort paths to FASTA files
    fasta_files = im.sort_data(fasta_files, sort_key=lambda x: x.lower())

    if cds_input is False:

        print('Number of inputs: {0}'.format(len(fasta_files)))

        # gene prediction step
        gp_results = cf.predict_genes(fasta_files, ptf_path,
                                      translation_table, prodigal_mode,
                                      cpu_cores, temp_directory,
                                      output_directory)

        fasta_files, prodigal_path = gp_results

        # CDS extraction step
        cds_files = cf.extract_genes(fasta_files, prodigal_path,
                                     cpu_cores, temp_directory,
                                     output_directory)
    else:
        cds_files = fasta_files
        print('Number of inputs: {0}'.format(len(cds_files)))

    # create directory to store files from pre-process steps
    preprocess_dir = fo.join_paths(temp_directory, ['3_cds_preprocess'])
    fo.create_directory(preprocess_dir)

    # DNA sequences deduplication step
    distinct_dna_template = 'distinct_cds_{0}.fasta'
    distinct_seqs_file, dna_shelve = cf.exclude_duplicates(cds_files, preprocess_dir, cpu_cores,
                                                           distinct_dna_template)

    # sequence translation step
    ts_results = cf.translate_sequences(list(valid_seqids.keys()), distinct_seqs_file,
                                        preprocess_dir, translation_table,
                                        minimum_length, cpu_cores)

    dna_file, protein_file, ut_seqids, ut_lines = ts_results

    excluded += [valid_seqids[seqid] for seqid in ut_seqids]
    valid_dna_seqids = {k: v
                        for k, v in valid_seqids.items()
                        if k not in ut_seqids}

    # write info about invalid alleles to file
    invalid_alleles_file = fo.join_paths(output_directory,
                                         ['invalid_cds.txt'])
    invalid_alleles = im.join_list(ut_lines+ss_lines, '\n')
    fo.write_to_file(invalid_alleles, invalid_alleles_file, 'w', '\n')
    print('Info about untranslatable and small sequences '
          'stored in {0}'.format(invalid_alleles_file))

    # protein sequences deduplication step
    distinct_prot_template = 'distinct_prots_{0}.fasta'
    ds_results = cf.exclude_duplicates([protein_file], preprocess_dir, 1,
                                       distinct_prot_template, True)

    distinct_protein_seqs, distinct_prots_file = ds_results

    distinct_protein_seqids = [v[0] for k, v in distinct_protein_seqs.items()]
    distinct_protein_seqids.sort(key=lambda y: y.lower())

    # write protein FASTA file
    qc_protein_file = os.path.join(preprocess_dir, 'filtered_proteins.fasta')
    indexed_protein_file = SeqIO.index(protein_file, 'fasta')
    fao.get_sequences_by_id(indexed_protein_file, distinct_protein_seqids,
                            qc_protein_file)

    # write DNA FASTA file
    qc_dna_file = os.path.join(preprocess_dir, 'filtered_dna.fasta')
    indexed_dna_file = SeqIO.index(dna_file, 'fasta')
    fao.get_sequences_by_id(indexed_dna_file, distinct_protein_seqids,
                            qc_dna_file)

    print('\nKept {0} sequences after filtering the initial '
          'sequences.'.format(len(distinct_protein_seqids)))

    # determine protein exact matches to schema seqs
    schema_loci = [file for file in os.listdir(schema_directory) if '.fasta' in file]
    schema_loci = [fo.join_paths(schema_directory, [file])
                   for file in schema_loci]

    # use protein sequences hashes to search for exact matches in loci
    prot_hashes = list(distinct_protein_seqs.keys())

    # find exact matches
    prot_matches = {}
    for locus in schema_loci:
        matches = exact_matches(locus, prot_hashes)
        if len(matches) > 0:
            prot_matches[locus] = matches

    # get identifiers for alleles that code for proteins
    prot_matches_ids = {}
    for locus, v in prot_matches.items():
        matches_ids = []
        for h in v:
            matches_ids += distinct_protein_seqs[h]
        prot_matches_ids[locus] = matches_ids

    # determine the allele identifier for exact matches
    # for each locus with exact matches
    # import DNA sequences from locus
    # use index for file with distinct DNA sequences to get DNA sequences
    # check if they are in locus sequences
    # assign ID if they are
    # new allele if they are not
    

    # remove identifiers from exact matches and perform clustering



def main(input_files, schema_directory, output_directory, ptf_path,
         blast_score_ratio, minimum_length, translation_table,
         size_threshold, word_size, window_size, clustering_sim,
         representative_filter, intra_filter, cpu_cores, blast_path,
         cds_input, prodigal_mode, no_cleanup):

    print('Prodigal training file: {0}'.format(ptf_path))
    print('CPU cores: {0}'.format(cpu_cores))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum sequence length: {0}'.format(minimum_length))
    print('Size threshold: {0}'.format(size_threshold))
    print('Word size: {0}'.format(word_size))
    print('Window size: {0}'.format(window_size))
    print('Clustering similarity: {0}'.format(clustering_sim))
    print('Representative filter: {0}'.format(representative_filter))
    print('Intra-cluster filter: {0}'.format(intra_filter))

    results = allele_calling(input_files, schema_directory, output_directory,
                             ptf_path, blast_score_ratio, minimum_length,
                             translation_table, size_threshold, word_size,
                             window_size, clustering_sim, representative_filter,
                             intra_filter, cpu_cores, blast_path,
                             prodigal_mode, cds_input)

    # remove temporary files
    if no_cleanup is False:
        fo.delete_directory(results[1])


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', nargs='?', type=str,
                        required=True, dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line.')

    parser.add_argument('-g', '--schema-directory', type=str,
                        required=True, dest='schema_directory',
                        help='')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Output directory where the process will store '
                             'intermediate files and create the schema\'s '
                             'directory.')

    parser.add_argument('--ptf', '--training-file', type=str,
                        required=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('--bsr', '--blast-score-ratio', type=float,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', '--minimum-length', type=int,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length value. Coding sequences '
                             'shorter than this value are excluded.')

    parser.add_argument('--t', '--translation-table', type=int,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences.')

    parser.add_argument('--st', '--size-threshold', type=float,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. Added to the '
                             'schema\'s config file and used to identify '
                             'alleles with a length value that deviates from '
                             'the locus length mode during the allele calling '
                             'process.')

    parser.add_argument('--cpu', '--cpu-cores', type=int,
                        required=False, default=1, dest='cpu_cores',
                        help='Number of CPU cores that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores).')

    parser.add_argument('--b', '--blast-path', type=str,
                        required=False, default='', dest='blast_path',
                        help='Path to the BLAST executables.')

    parser.add_argument('--pm', '--prodigal-mode', required=False,
                        choices=['single', 'meta'],
                        default='single', dest='prodigal_mode',
                        help='Prodigal running mode.')

    parser.add_argument('--CDS', required=False, action='store_true',
                        dest='cds_input',
                        help='If provided, input is a single or several FASTA '
                             'files with coding sequences.')

    parser.add_argument('--no-cleanup', required=False, action='store_true',
                        dest='no_cleanup',
                        help='If provided, intermediate files generated '
                             'during process execution are not removed at '
                             'the end.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main()