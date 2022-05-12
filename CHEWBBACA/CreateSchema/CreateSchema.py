#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the creation of a whole genome multi locus sequence
typing (wgMLST) schema seed.

The process starts by predicting and extracting genes from all input
files/assemblies (skipped if a single or several FASTA files with coding
sequences are passed as input and the ``--CDS`` flag is included in the
command). This is followed by the identification and removal of repeated
and small sequences from the set of coding sequences extracted from all
files/assemblies. 'Distinct' DNA sequences are translated and the resulting
protein sequences go through another deduplication step (repeated protein
sequences result from DNA sequences that code for the same protein). Protein
sequences are sorted and clustered based on the percentage of shared distinct
minimizers. Two filters are applied to each cluster to identify and exclude
sequences that are highly similar to the cluster's representative or to other
sequences in the cluster (the threshold is set at >=90% shared minimizers).
BLASTp is used to determine the similarity of the remaining sequences in each
cluster, if the clusters are not singletons (have a single protein sequence),
and exclude sequences based on high BLAST Score Ratio (BSR) values. A final
step uses BLASTp to determine the similarity of this set of sequences and
determine the 'distinct' set of loci that constitute the schema seed. The
schema seed contains 1 FASTA file per distinct locus and a single
representative allele for each locus.

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

- ``--n``, ``schema_name`` : Name given to the folder that will store the
  schema files.

    - e.g.: ``my_schema``

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
                       file_operations as fo,
                       fasta_operations as fao,
                       sequence_manipulation as sm,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 blast_wrapper as bw,
                                 core_functions as cf,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


def create_schema_structure(schema_seed_fasta, output_directory,
                            temp_directory, schema_name):
    """ Creates the schema seed directory with one FASTA file per
        distinct locus and the `short` directory with the FASTA files
        used to save the representative sequences.

        Parameters
        ----------
        schema_seed_fasta : str
            Path to the FASTA file that contains the sequences that
            constitute the schema seed. Each FASTA record in the file
            is a representative sequence chosen for a locus.
        output_directory : str
            Path to the main output directory of the process.
        temp_directory : str
            Path to the directory where the FASTA file with the
            records for the schema seed will be saved to.
        schema_name : str
            Name for the schema's directory.

        Returns
        -------
        schema_files : list
            List with the paths to the FASTA files in the schema seed.
    """

    # add allele identifier to all sequences
    schema_records = ['>{0}\n{1}'.format(im.replace_multiple_characters(rec.id, ct.CHAR_REPLACEMENTS) + '_1', str(rec.seq))
                      for rec in SeqIO.parse(schema_seed_fasta, 'fasta')]

    final_records = os.path.join(temp_directory, 'schema_seed_final.fasta')
    fo.write_lines(schema_records, final_records)

    schema_dir = fo.join_paths(output_directory, [schema_name])
    fo.create_directory(schema_dir)

    # create directory and schema files
    filenames = [record.id[:-2] for record in SeqIO.parse(final_records, 'fasta')]
    filenames = [im.replace_multiple_characters(f, ct.CHAR_REPLACEMENTS) for f in filenames]
    file_paths = (fo.join_paths(schema_dir, ['{0}.fasta'.format(f)]) for f in filenames)
    schema_files = fao.split_fasta(final_records, file_paths, 1)
    fo.create_short(schema_files, schema_dir)

    return schema_files


def create_schema_seed(input_files, output_directory, schema_name, ptf_path,
                       blast_score_ratio, minimum_length, translation_table,
                       size_threshold, word_size, window_size, clustering_sim,
                       representative_filter, intra_filter, cpu_cores, blast_path,
                       prodigal_mode, cds_input):
    """ Creates a schema seed based on a set of FASTA files that
        contain genome assemblies or coding sequences.
    """

    # define directory for temporary files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # read file with paths to input files
    fasta_files = fo.read_lines(input_files, strip=True)

    # sort paths to FASTA files
    fasta_files = im.sort_iterable(fasta_files, sort_key=lambda x: x.lower())

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
    ds_results = cf.exclude_duplicates(cds_files, preprocess_dir, cpu_cores,
                                       distinct_dna_template, False)

    schema_seqids, distinct_seqs_file = ds_results

    # determine small sequences step
    ss_results = cf.exclude_small(distinct_seqs_file, minimum_length)

    small_seqids, ss_lines = ss_results

    # exclude seqids of small sequences
    schema_seqids = list(set(schema_seqids) - set(small_seqids))
    schema_seqids = im.sort_iterable(schema_seqids, sort_key=lambda x: x.lower())

    # sequence translation step
    ts_results = cf.translate_sequences(schema_seqids, distinct_seqs_file,
                                        preprocess_dir, translation_table,
                                        minimum_length, cpu_cores)

    dna_file, protein_file, ut_seqids, ut_lines = ts_results

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
                                       distinct_prot_template)

    distinct_protein_seqs, distinct_prots_file = ds_results

    distinct_protein_seqs.sort(key=lambda y: y.lower())

    # write protein FASTA file
    qc_protein_file = os.path.join(preprocess_dir, 'filtered_proteins.fasta')
    indexed_protein_file = SeqIO.index(protein_file, 'fasta')
    fao.get_sequences_by_id(indexed_protein_file, distinct_protein_seqs,
                            qc_protein_file)

    # write DNA FASTA file
    qc_dna_file = os.path.join(preprocess_dir, 'filtered_dna.fasta')
    indexed_dna_file = SeqIO.index(dna_file, 'fasta')
    fao.get_sequences_by_id(indexed_dna_file, distinct_protein_seqs,
                            qc_dna_file)

    print('\nKept {0} sequences after filtering the initial '
          'sequences.'.format(len(distinct_protein_seqs)))

    # protein clustering step
    # read protein sequences
    proteins = fao.import_sequences(qc_protein_file)

    # create directory to store clustering data
    clustering_dir = fo.join_paths(temp_directory, ['4_clustering'])
    fo.create_directory(clustering_dir)

    cs_results = cf.cluster_sequences(proteins, word_size, window_size,
                                      clustering_sim, None, True,
                                      1, 1, clustering_dir, cpu_cores,
                                      'clusters', True, False)

    # clustering pruning step
    cp_results = cf.cluster_representative_filter(cs_results,
                                                  representative_filter,
                                                  clustering_dir,
                                                  'repfilter')

    clusters, excluded_seqids = cp_results

    # remove excluded seqids
    schema_seqids = list(set(distinct_protein_seqs) - excluded_seqids)

    # intra cluster pruner step
    cip_results = cf.cluster_intra_filter(clusters, proteins,
                                          word_size, intra_filter,
                                          clustering_dir,
                                          'intrafilter')

    clusters, intra_excluded = cip_results

    # remove excluded seqids - we get set of sequences from clusters
    # plus singletons
    schema_seqids = list(set(schema_seqids) - set(intra_excluded))

    # BLASTp clusters step
    blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
    makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)

    if len(clusters) > 0:
        blasting_dir = fo.join_paths(clustering_dir, ['cluster_blaster'])
        fo.create_directory(blasting_dir)

        blast_results, ids_dict = cf.blast_clusters(clusters, proteins,
                                                    blasting_dir, blastp_path,
                                                    makeblastdb_path, cpu_cores,
                                                    'blast')

        blast_files = im.flatten_list(blast_results)

        # compute and exclude based on BSR
        blast_excluded_alleles = [sm.apply_bsr(fo.read_tabular(file),
                                               indexed_dna_file,
                                               blast_score_ratio,
                                               ids_dict)
                                  for file in blast_files]

        # merge bsr results
        blast_excluded_alleles = im.flatten_list(blast_excluded_alleles)

        blast_excluded_alleles = [ids_dict[seqid] for seqid in blast_excluded_alleles]
        schema_seqids = list(set(schema_seqids) - set(blast_excluded_alleles))
        print('\n\nRemoved {0} sequences based on high BSR value with '
              'other sequences.'.format(len(set(blast_excluded_alleles))))

    # perform final BLAST to identify similar sequences that do not
    # share many/any kmers
    print('Total of {0} sequences to compare in final BLAST.'.format(len(schema_seqids)))

    # sort seqids before final BLASTp to ensure consistent results
    schema_seqids = im.sort_iterable(schema_seqids, sort_key=lambda x: x.lower())

    # create directory for final BLASTp
    final_blast_dir = fo.join_paths(temp_directory, ['5_final_blast'])
    fo.create_directory(final_blast_dir)

    # create FASTA file with remaining sequences
    beta_file = os.path.join(final_blast_dir, 'pre_schema_seed.fasta')
    fao.get_sequences_by_id(proteins, schema_seqids, beta_file)

    # change sequence identifiers to avoid BLAST error
    # related to sequence header ength limite
    integer_seqids = os.path.join(final_blast_dir, 'pre_schema_seed_int.fasta')
    ids_dict2 = fao.integer_headers(beta_file, integer_seqids)

    # create BLASTp database
    blast_db = fo.join_paths(final_blast_dir, ['pre_schema_seed_int'])
    db_stderr = bw.make_blast_db(makeblastdb_path, integer_seqids, blast_db, 'prot')

    if len(db_stderr) > 0:
        sys.exit(db_stderr)

    # divide FASTA file into groups of 100 sequences to reduce
    # execution time for large sequence sets
    file_num = math.ceil(len(schema_seqids)/100)
    filenames = ['split{0}'.format(i+1) for i in range(0, file_num)]
    file_paths = (fo.join_paths(final_blast_dir, [f]) for f in filenames)
    splitted_fastas = fao.split_fasta(integer_seqids, filenames, 100)

    blast_outputs = ['{0}/{1}_blast_out.tsv'.format(final_blast_dir,
                                                    fo.file_basename(file, suffix=False))
                     for file in splitted_fastas]

    # add common arguments to all sublists
    blast_inputs = [[blastp_path, blast_db, file,
                     blast_outputs[i], 1, 1, bw.run_blast]
                    for i, file in enumerate(splitted_fastas)]

    print('Performing final BLASTp...')
    blast_stderr = mo.map_async_parallelizer(blast_inputs,
                                             mo.function_helper,
                                             cpu_cores,
                                             show_progress=True)

    blast_stderr = im.flatten_list(blast_stderr)
    if len(blast_stderr) > 0:
        sys.exit(blast_stderr)

    # concatenate files with BLASTp results
    blast_output = fo.join_paths(final_blast_dir, ['blast_out_concat.tsv'])
    blast_output = fo.concatenate_files(blast_outputs, blast_output)

    final_excluded = sm.apply_bsr(fo.read_tabular(blast_output),
                                  indexed_dna_file,
                                  blast_score_ratio,
                                  ids_dict2)
    final_excluded = [ids_dict2[seqid] for seqid in final_excluded]

    schema_seqids = list(set(schema_seqids) - set(final_excluded))

    print('\nRemoved {0} sequences that were highly similar '
          'to other sequences.'.format(len(final_excluded)))

    output_schema = os.path.join(final_blast_dir, 'schema_seed.fasta')

    # create file with the schema representative sequences
    fao.get_sequences_by_id(indexed_dna_file, schema_seqids, output_schema)

    schema_files = create_schema_structure(output_schema, output_directory,
                                           final_blast_dir, schema_name)

    return [schema_files, temp_directory]


def main(input_files, output_directory, schema_name, ptf_path,
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

    results = create_schema_seed(input_files, output_directory, schema_name,
                                 ptf_path, blast_score_ratio, minimum_length,
                                 translation_table, size_threshold, word_size,
                                 window_size, clustering_sim, representative_filter,
                                 intra_filter, cpu_cores, blast_path,
                                 prodigal_mode, cds_input)

    # remove temporary files
    if no_cleanup is False:
        exists = fo.delete_directory(results[1])
        if exists is True:
            print('Could not delete intermediate files located in {0}'.format(results[1]))

    # print message about schema that was created
    print('Created schema seed with {0} loci.'.format(len(results[0])))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', nargs='?', type=str,
                        required=True, dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Output directory where the process will store '
                             'intermediate files and create the schema\'s '
                             'directory.')

    parser.add_argument('--n', '--schema-name', type=str,
                        required=False, default='schema_seed',
                        dest='schema_name',
                        help='Name given to the folder that will store '
                             'the schema files.')

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


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
