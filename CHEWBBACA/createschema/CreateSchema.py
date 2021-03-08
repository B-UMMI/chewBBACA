#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the creation of a whole genome multi locus sequence
typing (wgMLST) schema seed.

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

Code documentation
------------------
"""


import os
import sys
import argparse

from Bio import SeqIO

try:
    from utils import (constants as ct,
                       blast_wrapper as bw,
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
                                 gene_prediction as gp,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_clustering as sc,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


def predict_genes(fasta_files, ptf_path, translation_table,
                  prodigal_mode, cpu_cores, temp_directory,
                  output_directory):
    """ Runs Prodigal to predict coding sequences in FASTA
        files with genomic sequence.

        Parameters
        ----------
        fasta_files : list
            List with paths to FASTA files with genomic
            sequences.
        ptf_path : str
            Path to the Prodigal training file. Should
            be False if a training file is not provided.
        translation_table : int
            Genetic code used to predict and translate
            coding sequences.
        prodigal_mode : str
            Prodigal execution mode.
        cpu_cores : int
            Number of process that will run Prodigal in
            parallel.
        temp_directory : str
            Path to the directory where temporary files
            should be stored in.
        output_directory : str
            Path to the main directory of the process.

        Returns
        -------
        A list with the following elements:
            fasta_files : list
                Input list without the paths to the files
                that Prodigal could not predict genes for.
            prodigal_path : str
                Path to the directory with the files with
                Prodigal results.

        Raises
        ------
        SystemExit
            If Prodigal could not predict genes for any of
            the input files.
    """

    # create directory to store files with Prodigal results
    prodigal_path = fo.join_paths(temp_directory, ['1_gene_prediction'])
    fo.create_directory(prodigal_path)

    # run Prodigal to determine CDSs for all input genomes
    print('\nPredicting gene sequences...\n')
    # divide input genomes into equal number of sublists for
    # maximum process progress resolution
    prodigal_inputs = im.divide_list_into_n_chunks(fasta_files,
                                                   len(fasta_files))
    # add common arguments to all sublists
    common_args = [prodigal_path, ptf_path, translation_table,
                   prodigal_mode, gp.main]
    prodigal_inputs = [i+common_args for i in prodigal_inputs]

    # run Prodigal to predict genes
    prodigal_results = mo.map_async_parallelizer(prodigal_inputs,
                                                 mo.function_helper,
                                                 cpu_cores,
                                                 show_progress=True)

    # determine if Prodigal predicted genes for all genomes
    failed, failed_file = gp.check_prodigal_results(prodigal_results,
                                                    output_directory)

    if len(failed) > 0:
        print('\nFailed to predict genes for {0} genomes.'.format(len(failed)))
        print('Info for failed cases stored in: {0}'.format(failed_file))

    # remove failed genomes from paths
    for f in failed:
        fasta_files.remove(f[0])

    if len(fasta_files) == 0:
        sys.exit('\nCould not predict gene sequences from any '
                 'of the input files.\nPlease provide input files '
                 'in the accepted FASTA format.')

    return [fasta_files, prodigal_path]


def extract_genes(fasta_files, prodigal_path, cpu_cores,
                  temp_directory, output_directory):
    """ Extracts coding sequences from FASTA files with genomic
        sequences and saves coding sequences and info about coding
        sequences to files.

        Parameters
        ----------
        fasta_files : list
            List with paths to FASTA files with genomic sequences.
        prodigal_path : str
            Path to the directory with the files with Prodigal
            results.
        cpu_cores : int
            Number of process that will extract coding sequences
            in parallel.
        temp_directory : str
            Path to the directory where temporary files should
            be stored in.
        output_directory : str
            Path to the main directory of the process.

        Returns
        -------
        cds_files : list
            List with paths to the FASTA files that contain
            extracted coding sequences.
    """

    # divide inputs into at least 20 sublists for 5% process
    # progress resolution
    num_chunks = 20 if cpu_cores < 20 else cpu_cores
    extractor_inputs = im.divide_list_into_n_chunks(fasta_files, num_chunks)
    # add common arguments and unique index/identifier
    cds_extraction_path = os.path.join(temp_directory,
                                       '2_cds_extraction')
    fo.create_directory(cds_extraction_path)
    extractor_inputs = [[extractor_inputs[i-1], prodigal_path,
                         cds_extraction_path, i, gp.cds_batch_extractor]
                        for i in range(1, len(extractor_inputs)+1)]

    # extract coding sequences
    print('\n\nExtracting coding sequences...\n')
    extracted_cdss = mo.map_async_parallelizer(extractor_inputs,
                                               mo.function_helper,
                                               cpu_cores,
                                               show_progress=True)

    total_extracted = sum([f[2] for f in extracted_cdss])
    print('\n\nExtracted a total of {0} coding sequences from {1} '
          'genomes.'.format(total_extracted, len(fasta_files)))

    # create full table file
    table_files = [f[0] for f in extracted_cdss]
    table_file = fo.join_paths(output_directory, ['cds_info.tsv'])
    table_header = 'Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n'
    fo.concatenate_files(table_files, table_file, table_header)
    fo.remove_files(table_files)

    cds_files = [f[1] for f in extracted_cdss]

    return cds_files


def exclude_duplicates(fasta_files, temp_directory, cpu_cores,
                       outfile_template):
    """ Identifies duplicated sequences in FASTA files and
        selects a distinct set of sequences.

        Parameters
        ----------
        fasta_files : list
            List with paths to FASTA files.
        temp_directory : str
            Path to the directory where new files will be
            created.
        cpu_cores : int
            Number of deduplication processes to run in
            parallel.
        outfile_template : str
            Template for the name of output files.

        Returns
        -------
        A list with the following elements:
            distinct_seqids : list
                List with the sequence identifiers of distinct
                sequences.
            distinct_seqs : str
                Path to the FASTA file with distinct sequences.
    """

    dedup_inputs = [[file,
                     fo.join_paths(temp_directory,
                                   [outfile_template.format(i+1)]),
                     sm.determine_distinct]
                    for i, file in enumerate(fasta_files)]

    # determine distinct sequences (keeps 1 seqid per sequence)
    print('\nRemoving duplicated sequences...')
    dedup_results = mo.map_async_parallelizer(dedup_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=False)

    # determine number of duplicated sequences
    repeated = sum([d[0] for d in dedup_results])

    # one last round after concatenating files
    if len(dedup_inputs) > 1:
        dedup_files = [f[1] for f in dedup_inputs]
        cds_file = fo.join_paths(temp_directory, ['distinct_seqs_concat.fasta'])
        cds_file = fo.concatenate_files(dedup_files, cds_file)
        distinct_seqs = fo.join_paths(temp_directory, ['distinct_seqs.fasta'])
        dedup_results = sm.determine_distinct(cds_file, distinct_seqs)
        repeated_seqs, distinct_seqids = dedup_results

        repeated += repeated_seqs

        print('Removed {0} repeated sequences.'.format(repeated))

        return [distinct_seqids, distinct_seqs]
    else:
        print('Removed {0} repeated sequences.'.format(repeated))

        return [dedup_results[0][1], dedup_inputs[0][1]]


def exclude_small(fasta_file, minimum_length):
    """ Identifies sequences smaller that a specified length
        value.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file.
        minimum_length : int
            Sequences with a length value below this value are
            considered small.

        Returns
        -------
        small_seqids : list
            List with the sequence identifiers of small
            sequences.
    """

    # determine small sequences and keep their seqids
    print('\nRemoving sequences smaller than {0}...'.format(minimum_length))
    small_seqids = sm.determine_small(fasta_file, minimum_length)

    ss_lines = ['{0}: smaller than {1} chars'.format(seqid, minimum_length)
                for seqid in small_seqids]

    print('Removed {0} sequences shorter than '
          '{1}.'.format(len(small_seqids), minimum_length))

    return [small_seqids, ss_lines]


def translate_sequences(sequence_ids, sequences_file, temp_directory,
                        translation_table, minimum_length, cpu_cores):
    """ Translates DNA sequences, returns information about
        sequences that are untranslatable and saves
        translatable DNA sequences and proteins to FASTA
        files.

        Parameters
        ----------
        sequence_ids : list
            List with the identifiers of the sequences that
            should be translated.
        sequences_file : str
            Path to the FASTA files with the DNA sequences.
        temp_directory : str
            Path to the directory where new files will be
            created.
        translation_table : int
            Genetic code used to predict and translate
            coding sequences.
        minimum_length : int
            Sequences with a length value below this value are
            considered small and are excluded.
        cpu_cores : int
            Number of translation processes to run in
            parallel.

        Returns
        -------
        A list with the following elements:
            dna_file : str
                Path to the FASTA with all DNA sequences
                that could be translated.
            protein_file : str
                Path to the FASTA file with the translated
                sequences.
            untrans_seqids : list
                List with the identifiers of the sequences
                that could not be translated.
            untrans_lines : list
                List with one string per untranslatable
                sequence. Each string has a sequence
                identifier and a small description that
                indicates why the sequence could not be
                translated.
    """

    # translate distinct DNA sequences
    print('\nTranslating {0} DNA sequences...'.format(len(sequence_ids)))

    # divide inputs into sublists
    translation_inputs = im.divide_list_into_n_chunks(sequence_ids, cpu_cores)

    # create paths to files with DNA sequences of translatable sequences
    dna_template = fo.join_paths(temp_directory, ['dna_{0}.fasta'])
    dna_files = [dna_template.format(i+1)
                 for i in range(len(translation_inputs))]

    # create paths to files with protein sequences
    protein_template = fo.join_paths(temp_directory, ['protein_{0}.fasta'])
    protein_files = [protein_template.format(i+1)
                     for i in range(len(translation_inputs))]

    # add common args to sublists
    common_args = [sequences_file, translation_table, minimum_length]
    translation_inputs = [[translation_inputs[i], *common_args, dna_files[i],
                           protein_files[i], sm.translate_coding_sequences]
                          for i in range(0, len(dna_files))]

    # translate sequences
    translation_results = mo.map_async_parallelizer(translation_inputs,
                                                    mo.function_helper,
                                                    cpu_cores,
                                                    show_progress=False)

    # concatenate files
    dna_file = fo.join_paths(temp_directory, ['dna.fasta'])
    dna_file = fo.concatenate_files(dna_files, dna_file)
    protein_file = fo.join_paths(temp_directory, ['protein.fasta'])
    protein_file = fo.concatenate_files(protein_files, protein_file)

    # determine sequences that could not be translated
    untrans_lines = []
    untrans_seqids = []
    for res in translation_results:
        if len(res[0]) > 0:
            untrans_lines.extend(['{0}: {1}'.format(r[0], r[1])
                                  for r in res[0]])
            untrans_seqids.extend([r[0] for r in res[0]])

    print('Removed {0} DNA sequences that could not be '
          'translated.'.format(len(untrans_seqids)))

    return [dna_file, protein_file, untrans_seqids, untrans_lines]


def cluster_sequences(sequences, word_size, window_size, clustering_sim,
                      representatives, grow_clusters, kmer_offset,
                      seq_num_cluster, temp_directory, cpu_cores,
                      file_prefix, divide, position):
    """ Clusters sequences based on the proportion of shared kmers.

        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        word_size : int
            Value k for the kmer size.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.
        representatives : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values.
        grow_clusters : bool
            If it is allowed to create new clusters.
        kmer_offset : int
            Value to indicate offset of consecutive kmers.
        minimizer : bool
            If clustering should be based on shared minimizers.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.
        temp_directory : str
            Path to the directory where the clustering results
            will be saved to.
        cpu_cores : int
            Number of clustering processes to run in parallel.
        file_prefix : str
            A prefix to include in the names of created files.

        Returns
        -------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
    """

    # sort sequences by length
    sorted_seqs = {k: v for k, v in im.sort_data(sequences.items(),
                                                 sort_key=lambda x: len(x[1]),
                                                 reverse=True)}

    if divide is True:
        # divide sequences into sublists
        # do not divide based on number of available cores as it may
        # lead to different results with different number of cores
        cluster_inputs = im.split_iterable(sorted_seqs,
                                           int(len(sorted_seqs)/40+10))
    else:
        cluster_inputs = [sorted_seqs]

    common_args = [word_size, window_size, clustering_sim,
                   representatives, grow_clusters, kmer_offset,
                   position, seq_num_cluster,
                   sc.clusterer]
    cluster_inputs = [[c, *common_args] for c in cluster_inputs]

    # cluster proteins in parallel
    print('\nClustering sequences...')
    clustering_results = mo.map_async_parallelizer(cluster_inputs,
                                                   mo.function_helper,
                                                   cpu_cores)

    # merge clusters
    clusters = [d[0] for d in clustering_results]
    clusters = im.merge_dictionaries(clusters)
    rep_sequences = [d[1] for d in clustering_results]
    rep_sequences = im.merge_dictionaries(rep_sequences)

    # perform clustering with representatives
    if len(cluster_inputs) > 1:
        # cluster representatives
        rep_clusters = sc.clusterer(rep_sequences, word_size,
                                    window_size, clustering_sim,
                                    representatives, grow_clusters,
                                    kmer_offset, position,
                                    seq_num_cluster)

        merged_clusters = {}
        for k, v in rep_clusters[0].items():
            # merge clusters whose representatives are similar
            for n in v:
                # representatives from other clusters are added with
                # similarity score against new representative
                # clustered sequences from other clusters are added
                # with similarity score against the old representative
                add_seqids = [n] + [s for s in clusters[n[0]] if s[0] != n[0]]
                merged_clusters.setdefault(k, []).extend(add_seqids)

        clusters = merged_clusters

    print('Clustered {0} sequences into {1} '
          'clusters'.format(len(sorted_seqs), len(clusters)))

    # sort clusters
    clusters = {k: v for k, v in im.sort_data(clusters.items())}

    # write file with clustering results
    clusters_out = os.path.join(temp_directory,
                                '{0}.txt'.format(file_prefix))
    sc.write_clusters(clusters, clusters_out)

    return clusters


def cluster_representative_filter(clusters, representative_filter,
                                  output_directory, file_prefix):
    """ Excludes sequences from clusters based on the proportion
        of shared kmers with the representative. After removing
        highly similar sequences, excludes clusters that are
        singletons (only contain the representative).

        Parameters
        ----------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
        representative_filter : float
            Similarity threshold value. Sequences with
            equal or greater similarity value with the cluster's
            representative are excluded from clusters.
        output_directory : str
            Path to the directory where the clustering results
            will be saved to.
        file_prefix : str
            A prefix to include in the names of created files.

        Returns
        -------
        A list with the following elements:
            pruned_clusters : dict
                Clusters without the sequences that were highly
                similar to the cluster's representative and without
                the clusters that were singletons (only contained
                the representative).
            excluded_seqids : list
                List with the sequence identifiers of the sequences
                that were excluded from the clusters.
    """

    # remove sequences that are very similar to representatives
    pruning_results = sc.representative_pruner(clusters,
                                               representative_filter)

    pruned_clusters, excluded_seqids = pruning_results

    # get identifiers of excluded sequences
    # determine set because same seqids could be in several clusters
    excluded_seqids = set([e[0] for e in excluded_seqids])

    print('Removed {0} sequences based on high similarity with '
          'cluster representative.'.format(len(excluded_seqids)))

    # remove excluded seqids from clusters without high representative
    # similarity
    pruned_clusters = {k: [e for e in v if e[0] not in excluded_seqids]
                       for k, v in pruned_clusters.items()}

    # write file with pruning results
    pruned_out = os.path.join(output_directory,
                              '{0}_clusters.txt'.format(file_prefix))
    sc.write_clusters(pruned_clusters, pruned_out)

    # identify singletons and exclude those clusters
    singletons = im.select_clusters(pruned_clusters, 0)
    print('Found {0} singletons.'.format(len(singletons)))

    pruned_clusters = im.remove_entries(pruned_clusters, singletons)

    # determine number of sequences that still need to be evaluated
    # +1 to include representative
    clustered_sequences = sum([len(v)+1 for k, v in pruned_clusters.items()])
    print('Remaining sequences after representative and singleton '
          'pruning: {0}'.format(clustered_sequences))

    return [pruned_clusters, excluded_seqids]


def cluster_intra_filter(clusters, sequences, word_size,
                         intra_filter, output_directory,
                         file_prefix):
    """ Determines similarity between clustered sequences and
        excludes sequences that are highly similar to other clustered
        sequences.

        Parameters
        ----------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        word_size : int
            Value k for the kmer size.
        intra_filter : float
            Similarity threshold value. Sequences with
            equal or greater similarity value with other
            clustered sequences are excluded from clusters.
        output_directory : str
            Path to the directory where the clustering results
            will be saved to.
        file_prefix : str
            A prefix to include in the names of created files.

        Returns
        -------
        A list with the following elements:
            clusters : dict
                Clusters without the sequences that were highly
                similar to other clustered sequences.
            intra_excluded : list
                List with the identifiers of the sequences that
                were excluded.
    """

    # identify clusters with more than 1 sequence
    intra_clusters = {k: v for k, v in clusters.items() if len(v) > 1}

    excluded_seqids, excluded_sims = sc.intra_cluster_sim(intra_clusters,
                                                          sequences,
                                                          word_size,
                                                          intra_filter)

    intra_excluded = [v for k, v in excluded_seqids.items()]
    intra_excluded = im.flatten_list(intra_excluded)
    # get identifiers of excluded sequences
    # determine set because same seqids could be in several clusters
    intra_excluded = set(intra_excluded)
    print('Removed {0} sequences.'.format(len(intra_excluded)))

    # remove excluded seqids from clusters without high intra-similarity
    pruned_clusters = {k: [e for e in v if e[0] not in intra_excluded]
                       for k, v in clusters.items()}

    # write excluded to file
    intrasim_out = os.path.join(output_directory,
                                '{0}_excluded.txt'.format(file_prefix))
    sc.write_clusters(excluded_sims, intrasim_out)
    # write clusters to file
    intrasim_out = os.path.join(output_directory,
                                '{0}_clusters.txt'.format(file_prefix))
    sc.write_clusters(pruned_clusters, intrasim_out)

    # add key because it is representative identifier
    clustered_sequences = sum([len(v)+1 for k, v in pruned_clusters.items()])
    print('Remaining sequences after intra cluster pruning: '
          '{0}'.format(clustered_sequences))

    return [pruned_clusters, intra_excluded]


def blast_clusters(clusters, sequences, output_directory,
                   blastp_path, makeblastdb_path, cpu_cores,
                   file_prefix):
    """ Performs all-against-all comparisons between sequences in
        the same cluster through alignments with BLAST.

        Parameters
        ----------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        output_directory : str
            Path to the directory where the clustering results
            will be saved to.
        blastp_path : str
            Path to the BLAST executables.
        cpu_cores : int
            Number of BLAST processes to run in parallel.
        file_prefix : str
            A prefix to include in the names of created files.

        Returns
        -------
        A list with the following elements:
            blast_results : list
                List with paths to the files with BLAST results
                (one file per input).
            ids_dict : dict
                Dictionary that maps sequence identifiers to
                shorter and unique integer identifiers used
                to avoid errors during BLAST related with
                sequence headers/identifiers that exceed the
                length limit allowed by BLAST.
    """

    print('Clusters to BLAST: {0}'.format(len(clusters)))

    # create FASTA file with sequences in clusters
    clustered_seqs_file = fo.join_paths(output_directory,
                                         ['{0}_clustered_proteins.fasta'.format(file_prefix)])
    clustered_sequences = [[k]+[e[0] for e in v] for k, v in clusters.items()]
    clustered_sequences = im.flatten_list(clustered_sequences)
    # do not include duplicate identifiers
    clustered_sequences = list(set(clustered_sequences))
    fao.get_sequences_by_id(sequences, clustered_sequences,
                            clustered_seqs_file)

    # create FASTA file with replaced headers to avoid header
    # length limitation in BLAST
    integer_clusters = fo.join_paths(output_directory,
                                      ['{0}_clustered_proteins_int.fasta'.format(file_prefix)])
    ids_dict = fao.integer_headers(clustered_seqs_file, integer_clusters)

    # create BLAST DB
    blast_db = fo.join_paths(output_directory,
                             ['{0}_clustered_proteins_int'.format(file_prefix)])
    db_stderr = bw.make_blast_db(makeblastdb_path, integer_clusters,
                                 blast_db, 'prot')
    if len(db_stderr) > 0:
        sys.exit(db_stderr)

    blast_results_dir = os.path.join(output_directory,
                                     '{0}_results'.format(file_prefix))
    os.mkdir(blast_results_dir)

    # create files with replaced sequence identifiers per cluster
    seqids_to_blast = sc.blast_inputs(clusters, blast_results_dir, ids_dict)

    # distribute clusters per available cores
    process_num = 20 if cpu_cores <= 20 else cpu_cores
    splitted_seqids = mo.split_genes_by_core(seqids_to_blast,
                                             process_num,
                                             'seqcount')

    common_args = [integer_clusters, blast_results_dir, blastp_path,
                   blast_db, sc.cluster_blaster]

    splitted_seqids = [[s, *common_args] for s in splitted_seqids]

    # create the FASTA files with the protein sequences before BLAST?
    print('BLASTing protein sequences in each cluster...\n')

    # BLAST each sequences in a cluster against every sequence in that cluster
    blast_results = mo.map_async_parallelizer(splitted_seqids,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    print('\n\nFinished BLASTp. Determining schema representatives...')

    return [blast_results, ids_dict]


def create_schema_structure(input_files, output_directory, temp_directory,
                            schema_name):
    """
    """

    # add allele identifier to all sequences
    schema_records = ['>{0}\n{1}'.format(im.replace_multiple_characters(rec.id, ct.CHAR_REPLACEMENTS) + '_1', str(rec.seq))
                      for rec in SeqIO.parse(input_files, 'fasta')]

    final_records = os.path.join(temp_directory, 'schema_seed_final.fasta')
    fo.write_lines(schema_records, final_records)

    schema_dir = fo.join_paths(output_directory, [schema_name])
    fo.create_directory(schema_dir)

    # create directory and schema files
    filenames = (record.id[:-2] for record in SeqIO.parse(final_records, 'fasta'))
    schema_files = fao.split_fasta(final_records, schema_dir, 1, filenames)
    fo.create_short(schema_files, schema_dir)

    return schema_files


def create_schema_seed(input_files, output_directory, schema_name, ptf_path,
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
        gp_results = predict_genes(fasta_files, ptf_path,
                                   translation_table, prodigal_mode,
                                   cpu_cores, temp_directory,
                                   output_directory)

        fasta_files, prodigal_path = gp_results

        # CDS extraction step
        cds_files = extract_genes(fasta_files, prodigal_path,
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
    ds_results = exclude_duplicates(cds_files, preprocess_dir, cpu_cores,
                                    distinct_dna_template)

    schema_seqids, distinct_seqs_file = ds_results

    # determine small sequences step
    ss_results = exclude_small(distinct_seqs_file, minimum_length)

    small_seqids, ss_lines = ss_results

    # exclude seqids of small sequences
    schema_seqids = list(set(schema_seqids) - set(small_seqids))
    schema_seqids = im.sort_data(schema_seqids, sort_key=lambda x: x.lower())

    # sequence translation step
    ts_results = translate_sequences(schema_seqids, distinct_seqs_file,
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
    ds_results = exclude_duplicates([protein_file], preprocess_dir, 1,
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

    cs_results = cluster_sequences(proteins, word_size, window_size,
                                   clustering_sim, None, True,
                                   1, 1, clustering_dir, cpu_cores,
                                   'clusters', True, False)

    # clustering pruning step
    cp_results = cluster_representative_filter(cs_results,
                                               representative_filter,
                                               clustering_dir,
                                               'repfilter')

    clusters, excluded_seqids = cp_results

    # remove excluded seqids
    schema_seqids = list(set(distinct_protein_seqs) - excluded_seqids)

    # intra cluster pruner step
    cip_results = cluster_intra_filter(clusters, proteins,
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

    blasting_dir = fo.join_paths(clustering_dir, ['cluster_blaster'])
    fo.create_directory(blasting_dir)

    blast_results, ids_dict = blast_clusters(clusters, proteins,
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
    #############################
    ### write list of excluded alleles per cluster to file!!!!
    blast_excluded_alleles = [ids_dict[seqid] for seqid in blast_excluded_alleles]
    schema_seqids = list(set(schema_seqids) - set(blast_excluded_alleles))

    # perform final BLAST to identify similar sequences that do not
    # share many/any kmers
    print('Total of {0} sequences to compare in final BLAST.'.format(len(schema_seqids)))

    # sort seqids before final BLASTp to ensure consistent results
    schema_seqids = im.sort_data(schema_seqids, sort_key=lambda x: x.lower())

    # create directory for final BLASTp
    final_blast_dir = fo.join_paths(temp_directory, ['5_final_blast'])
    fo.create_directory(final_blast_dir)

    beta_file = os.path.join(final_blast_dir, 'pre_schema_seed.fasta')
    fao.get_sequences_by_id(proteins, schema_seqids, beta_file)

    integer_seqids = os.path.join(final_blast_dir, 'pre_schema_seed_int.fasta')
    ids_dict2 = fao.integer_headers(beta_file, integer_seqids)

    blast_db = fo.join_paths(final_blast_dir, ['pre_schema_seed_int'])
    db_stderr = bw.make_blast_db(makeblastdb_path, integer_seqids, blast_db, 'prot')

    if len(db_stderr) > 0:
        sys.exit(db_stderr)

    blast_output = '{0}/{1}_blast_out.tsv'.format(final_blast_dir,
                                                  'pre_schema_seed')
    blast_stderr = bw.run_blast(blastp_path, blast_db, integer_seqids,
                                blast_output, 1, cpu_cores)

    if len(blast_stderr) > 0:
        sys.exit(blast_stderr)

    final_excluded = sm.apply_bsr(fo.read_tabular(blast_output),
                                  indexed_dna_file,
                                  blast_score_ratio,
                                  ids_dict2)
    final_excluded = [ids_dict2[seqid] for seqid in final_excluded]

    schema_seqids = list(set(schema_seqids) - set(final_excluded))

    print('Removed {0} loci that were similar to other loci '
          'in the schema.'.format(len(final_excluded)))

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
    print('Intra filter: {0}'.format(intra_filter))

    results = create_schema_seed(input_files, output_directory, schema_name,
                                 ptf_path, blast_score_ratio, minimum_length,
                                 translation_table, size_threshold, word_size,
                                 window_size, clustering_sim, representative_filter,
                                 intra_filter, cpu_cores, blast_path,
                                 prodigal_mode, cds_input)

    # remove temporary files
    if no_cleanup is False:
        fo.delete_directory(results[1])

    # print message about schema that was created
    print('Created schema seed with {0} loci.'.format(len(results[0])))


def parse_arguments():


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='Output directory where the process will store '
                             'intermediate files and create the schema\'s '
                             'directory.')

    parser.add_argument('--n', type=str, required=False,
                        default='schema_seed', dest='schema_name',
                        help='Name to give to folder that will store '
                             'the schema files.')

    parser.add_argument('--ptf', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('--bsr', type=float,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', type=int,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length accepted for a '
                             'coding sequence to be included in the schema.')

    parser.add_argument('--t', type=int,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences.')

    parser.add_argument('--st', type=float,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent will be classified as ASM/ALM.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores).')

    parser.add_argument('--b', type='blastp', required=False,
                        default='', dest='blast_path',
                        help='Path to the BLAST executables.')

    parser.add_argument('--CDS', required=False, action='store_true',
                        dest='cds_input',
                        help='Input is a FASTA file with one representative '
                             'sequence per gene in the schema.')

    parser.add_argument('--pm', required=False, choices=['single', 'meta'],
                        default='single', dest='prodigal_mode',
                        help='Prodigal running mode.')

    parser.add_argument('--no_cleanup', required=False,
                        action='store_true', dest='no_cleanup',
                        help='Delete intermediate files at the end.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
