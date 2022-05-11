#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to

Code documentation
------------------
"""


import os
import sys

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
                  prodigal_mode, cpu_cores, output_directory,
                  parent_directory):
    """ Runs Prodigal to predict coding sequences from FASTA
        files with genomic sequences.

    Parameters
    ----------
    fasta_files : list
        List of paths to FASTA files with genomic
        sequences.
    ptf_path : str
        Path to the Prodigal training file. Should
        be NoneType if a training file is not provided.
    translation_table : int
        Genetic code used to predict and translate
        coding sequences.
    prodigal_mode : str
        Prodigal execution mode.
    cpu_cores : int
        Number of processes that will run Prodigal in
        parallel.
    output_directory : str
        Path to the directory where output files
        with Prodigal's results will be stored in.
    parent_directory : str
        Path to the main output directory of the process.

    Returns
    -------
    failed_info : list
        List that contains a list with the stderr for the
        cases that Prodigal failed to predict genes for
        and the path to the file with information about
        failed cases. Returns NoneType if gene prediction
        succeeded for all inputs.
    """

    # divide input genomes into equal number of sublists for
    # maximum process progress resolution
    prodigal_inputs = im.divide_list_into_n_chunks(fasta_files,
                                                   len(fasta_files))

    common_args = [output_directory, ptf_path,
                   translation_table, prodigal_mode]

    # add common arguments to all sublists
    prodigal_inputs = im.multiprocessing_inputs(prodigal_inputs,
                                                common_args,
                                                gp.main)

    # run Prodigal to predict genes
    prodigal_results = mo.map_async_parallelizer(prodigal_inputs,
                                                 mo.function_helper,
                                                 cpu_cores,
                                                 show_progress=True)

    # determine if Prodigal predicted genes for all genomes
    failed_info = gp.check_prodigal_results(prodigal_results,
                                              parent_directory)

    return failed_info


def extract_genes(fasta_files, prodigal_path, cpu_cores,
                  temp_directory, parent_directory):
    """ Extracts coding sequences from FASTA files with genomic
        sequences and saves coding sequences and info about coding
        sequences.

    Parameters
    ----------
    fasta_files : list
        List of paths to FASTA files with genomic sequences.
    prodigal_path : str
        Path to the directory with the files with Prodigal
        results.
    cpu_cores : int
        Number of processes that will extract coding sequences
        in parallel.
    temp_directory : str
        Path to the directory where FASTA files with extracted
        coding sequences will be stored in.
    parent_directory : str
        Path to the main output directory of the process.

    Returns
    -------
    cds_files : list
        List with paths to the FASTA files that contain
        extracted coding sequences.
    total_extracted : int
        Total number of coding sequences extracted from the
        input fasta files.
    """

    # divide inputs into 15 sublists for ~7% process
    # progress resolution
    num_chunks = 15
    extractor_inputs = im.divide_list_into_n_chunks(fasta_files, num_chunks)

    # add common arguments and unique index/identifier
    output_index = [i+1 for i in range(len(extractor_inputs))]
    extractor_inputs = im.aggregate_iterables([extractor_inputs, output_index])
    extractor_inputs = im.multiprocessing_inputs(extractor_inputs,
                                                 [prodigal_path, temp_directory],
                                                 gp.cds_batch_extractor)

    # extract coding sequences
    extracted_cdss = mo.map_async_parallelizer(extractor_inputs,
                                               mo.function_helper,
                                               cpu_cores,
                                               show_progress=True)

    total_extracted = sum([f[2] for f in extracted_cdss])

    # create full table file
    table_files = [f[0] for f in extracted_cdss]
    table_file = fo.join_paths(parent_directory, ['cds_info.tsv'])
    fo.concatenate_files(table_files, table_file, ct.CDS_TABLE_HEADER)
    fo.remove_files(table_files)

    cds_files = [f[1] for f in extracted_cdss]

    return [cds_files, total_extracted]


def exclude_duplicates(fasta_files, temp_directory, cpu_cores,
                       outfile_template, ids_map, ids=False, polyline=False):
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

    # create groups of inputs for multiprocessing
    output_files = [fo.join_paths(temp_directory,
                                  [outfile_template.format(i+1)])
                    for i, file in enumerate(fasta_files)]

    inputs = im.aggregate_iterables([fasta_files, output_files])

    dedup_inputs = im.multiprocessing_inputs(inputs,
                                             [ids_map, ids],
                                             sm.determine_distinct)

    # determine distinct sequences (keeps 1 seqid per sequence)
    dedup_results = mo.map_async_parallelizer(dedup_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=False)

    distinct_seqids = []
    if polyline is True:
        # merge results
        repeated = 0
        merged_results = {}
        for p in dedup_results:
            r = fo.pickle_loader(p)
            for k, v in r.items():
                if k in merged_results:
                    stored_ids = im.polyline_decoding(merged_results[k])
                    merged_results[k] = im.polyline_encoding([stored_ids[0]]+sorted(stored_ids[1:]+v[1:]))
                else:
                    distinct_seqids.append(v[0])
                    merged_results[k] = im.polyline_encoding([int(v[0].split('protein')[-1])]+v[1:])
                repeated += len(v) -1

        # determine number of duplicated sequences
        # minus 2 so we do not count the seqid and genome integer
        # identifier for the representative record
        repeated = repeated - len(merged_results)
    else:
        # merge results
        # only using during the protein deduplication and first ID comes duplicated...
        repeated = 0
        merged_results = {}
        for p in dedup_results:
            r = fo.pickle_loader(p)
            for k, v in r.items():
                integer_list = []
                # ignore first entry because it's equal to the second entry
                for i in v[1:]:
                    split = i.split('-protein')
                    integer_list.extend([ids_map[split[0]], int(split[1])])

                if k in merged_results:
                    stored_ids = im.polyline_decoding(merged_results[k])
                    merged_results[k] = im.polyline_encoding([stored_ids+integer_list])
                else:
                    merged_results[k] = im.polyline_encoding(integer_list)
                    distinct_seqids.append(v[0])

                # exclude first and second that match sequence that will represent
                repeated += len(v) - 2


    # concatenate Fasta files from parallel processes
    dedup_files = [f[1] for f in dedup_inputs]
    cds_file = fo.join_paths(temp_directory, ['distinct_seqs_concat.fasta'])
    cds_file = fo.concatenate_files(dedup_files, cds_file)

    # create index for concatenated Fasta
    cds_index = SeqIO.index(cds_file, 'fasta')

    # define filename for file with distinct sequences
    distinct_seqs = fo.join_paths(temp_directory,
                                  [outfile_template.format('d')])
    # get the representative record for each distinct sequence
    fao.get_sequences_by_id(cds_index, distinct_seqids,
                            distinct_seqs, 20000)

    return [merged_results, distinct_seqs, repeated]


def exclude_small(fasta_file, minimum_length, variation=0):
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
        ss_lines : list
            List with one string per small sequence. Each string
            represents an exception message for a sequence that
            is small.
    """

    # determine small sequences and keep their seqids
    print('\nRemoving sequences smaller than {0} '
          'nucleotides...'.format(minimum_length), end='')
    small_seqids = sm.determine_small(fasta_file, minimum_length, variation)

    ss_lines = ['{0}: smaller than {1} chars'.format(seqid, minimum_length)
                for seqid in small_seqids]

    print('removed {0} sequences.'.format(len(small_seqids)))

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
        Path to the FASTA file with the DNA sequences.
    temp_directory : str
        Path to the directory where new files will be
        created.
    translation_table : int
        Genetic code used to translate coding sequences.
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
    translation_inputs = im.aggregate_iterables([translation_inputs,
                                                 dna_files,
                                                 protein_files])
    translation_inputs = im.multiprocessing_inputs(translation_inputs,
                                                   common_args,
                                                   sm.translate_coding_sequences)

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

    return [dna_file, protein_file, untrans_seqids, untrans_lines]


def cluster_sequences(sequences, word_size, window_size, clustering_sim,
                      representatives, grow_clusters, kmer_offset,
                      seq_num_cluster, temp_directory, cpu_cores,
                      file_prefix, divide, position):
    """ Clusters sequences based on the proportion of shared minimizers.

        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        word_size : int
            Value k for the k-mer size.
        window_size : int
            Value for the window size/number of consecutive
            k-mers per window.
        clustering_sim : float
            Similarity threshold to add a sequence to
            a cluster.
        representatives : dict
            Dictionary with k-mers as keys and a list with
            identifiers of sequences that contain that k-mer
            as values.
        grow_clusters : bool
            If it is allowed to create new clusters.
        kmer_offset : int
            Value to indicate offset of consecutive kmers.
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
        divide : bool
            If input sequences should be divided into smaller
            groups that can be processed in parallel.
        position : bool
            True if the start position for each k-mer should be saved.
            False otherwise.

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
    sorted_seqs = {k: v for k, v in im.sort_iterable(sequences.items(),
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
    print('\nClustering sequences based on the proportion '
          'of shared distinct minimizers...')
    clustering_results = mo.map_async_parallelizer(cluster_inputs,
                                                   mo.function_helper,
                                                   cpu_cores)

    # merge clusters
    clusters = [d[0] for d in clustering_results]
    clusters = im.merge_dictionaries(clusters)
    rep_sequences = [d[1] for d in clustering_results]
    rep_sequences = im.merge_dictionaries(rep_sequences)

    # perform clustering with representatives
    # this step does not run for AlleleCall!
    if len(cluster_inputs) > 1 and any(len(r) > 0 for r in rep_sequences) is True:
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
                # with similarity score against their representative
                add_seqids = [n] + [s for s in clusters[n[0]] if s[0] != n[0]]
                merged_clusters.setdefault(k, []).extend(add_seqids)

        clusters = merged_clusters

    print('Clustered {0} sequences into {1} '
          'clusters.'.format(len(sorted_seqs), len(clusters)))

    # sort clusters
    clusters = {k: v for k, v in im.sort_iterable(clusters.items())}

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
          'the cluster representative.'.format(len(excluded_seqids)))

    # remove excluded seqids from clusters without high representative
    # similarity
    pruned_clusters = {k: [e for e in v if e[0] not in excluded_seqids]
                       for k, v in pruned_clusters.items()}

    # write file with pruning results
    pruned_out = os.path.join(output_directory,
                              '{0}_clusters.txt'.format(file_prefix))
    sc.write_clusters(pruned_clusters, pruned_out)

    # identify singletons and exclude those clusters
    singletons = im.select_keys(pruned_clusters, 0)
    print('Identified and removed {0} singletons.'.format(len(singletons)))

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
            Value k for the k-mer size.
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
    print('Removed {0} sequences based on high similarity with '
          'other clustered sequences.'.format(len(intra_excluded)))

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
    print('Remaining sequences after intra-cluster pruning: '
          '{0}'.format(clustered_sequences))

    return [pruned_clusters, intra_excluded]


def blast_clusters(clusters, sequences, output_directory,
                   blastp_path, makeblastdb_path, cpu_cores,
                   file_prefix, only_rep=False):
    """ Uses BLAST to align sequences in the same clusters.

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
            Path to the `BLASTp` executable.
        makeblastdb_path : str
            Path to the `makeblastdb` executable.
        cpu_cores : int
            Number of BLASTp processes to run in parallel.
        file_prefix : str
            A prefix to include in the names of created files.

        Returns
        -------
        A list with the following elements:
            blast_results : list
                List with paths to the files with BLASTp results
                (one file per cluster).
            ids_dict : dict
                Dictionary that maps sequence identifiers to
                shorter and unique integer identifiers used
                to avoid errors during BLAST execution related with
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
    seqids_to_blast = sc.blast_inputs(clusters, blast_results_dir, ids_dict,
                                      only_rep)

    # distribute clusters per available cores
    process_num = 20 if cpu_cores <= 20 else cpu_cores
    splitted_seqids = mo.split_genes_by_core(seqids_to_blast,
                                             process_num,
                                             'seqcount')
    
    common_args = [integer_clusters, blast_results_dir, blastp_path,
                   blast_db, only_rep, sc.cluster_blaster]

    splitted_seqids = [[s, *common_args] for s in splitted_seqids]

    print('BLASTing protein sequences in each cluster...\n')

    # BLAST each sequences in a cluster against every sequence in that cluster
    blast_results = mo.map_async_parallelizer(splitted_seqids,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    return [blast_results, ids_dict]


def compute_bsr(subject_score, query_score):
    """ Computes the BLAST Score Ratio for an alignment
        between two sequences.

    Parameters
    ----------
    subject_score : float
        Alignment raw score computed by BLAST.
    query_score : float
        Raw score computed by BLAST for the
        self-alignment of the query sequence.

    Returns
    -------
    bsr : float
        BLAST Score Ratio for the alignment.
    """

    bsr = subject_score / query_score

    return bsr
