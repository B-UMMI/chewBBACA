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
except ModuleNotFoundError:
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
                  prodigal_mode, cpu_cores, output_directory):
    """Execute Prodigal to predict coding sequences from Fasta files.

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
    failed = {line[0]: line[1]
              for line in prodigal_results
              if line[1] == 0
              or isinstance(line[1], str) is True}

    return failed


def extract_genes(fasta_files, prodigal_path, cpu_cores,
                  temp_directory):
    """Extract coding sequences from FASTA files.

    Extract coding sequences from genomic Fasta files and
    saves coding sequences and info about coding sequences.

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
    table_file = fo.join_paths(temp_directory, [ct.CDS_COORDINATES_BASENAME])
    fo.concatenate_files(table_files, table_file, ct.CDS_TABLE_HEADER)
    fo.remove_files(table_files)

    cds_files = [f[1] for f in extracted_cdss]

    return [cds_files, total_extracted, table_file]


def exclude_duplicates(fasta_files, temp_directory, cpu_cores,
                       outfile_template, ids_map, protein=False,
                       only_seqids=False):
    """Identify duplicated sequences and select distinct set of sequences.

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
                                  [outfile_template.format(str(i+1)+'.fasta')])
                    for i, file in enumerate(fasta_files)]

    inputs = im.aggregate_iterables([fasta_files, output_files])

    dedup_inputs = im.multiprocessing_inputs(inputs,
                                             [ids_map[0]],
                                             sm.determine_distinct)

    # determine distinct sequences
    dedup_results = mo.map_async_parallelizer(dedup_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=False)

    repeated = 0
    merged_results = {}
    distinct_seqids = []
    # merge results
    for p in dedup_results:
        r = fo.pickle_loader(p)
        for k, v in r.items():
            if k in merged_results:
                stored_ids = im.polyline_decoding(merged_results[k])
                if protein is False:
                    rest = [v[i] for i in range(1, len(v), 2)]
                    merged_results[k] = im.polyline_encoding(stored_ids+rest)
                else:
                    merged_results[k] = im.polyline_encoding(stored_ids+v)
                repeated += (len(v)/2)
            else:
                seqid = '{0}-protein{1}'.format(ids_map[1][v[1]], v[0])
                distinct_seqids.append(seqid)
                if protein is False:
                    rep = v[0:2]
                    rest = [v[i] for i in range(3, len(v), 2)]
                    merged_results[k] = im.polyline_encoding(rep+rest)
                else:
                    merged_results[k] = im.polyline_encoding(v)

                repeated += (len(v)/2) - 1

    # save table with deduplicated records
    hash_table_file = fo.join_paths(temp_directory, [outfile_template.format('merged.hashtable')])
    fo.pickle_dumper(merged_results, hash_table_file)

    # remove intermediate deduplication tables
    fo.remove_files(dedup_results)

    # concatenate Fasta files from parallel processes
    dedup_files = [f[1] for f in dedup_inputs]
    cds_file = fo.join_paths(temp_directory, [outfile_template.format('concat.fasta')])
    cds_file = fo.concatenate_files(dedup_files, cds_file)

    # create index for concatenated Fasta
    cds_index = fao.index_fasta(cds_file)

    # define filename for file with distinct sequences
    distinct_seqs = fo.join_paths(temp_directory,
                                  [outfile_template.format('merged.fasta')])
    # get the representative record for each distinct sequence
    fao.get_sequences_by_id(cds_index, distinct_seqids,
                            distinct_seqs, 50000)

    fo.remove_files(dedup_files+[cds_file])

    if only_seqids is False:
        return [merged_results, distinct_seqs, repeated]
    else:
        return [distinct_seqids, distinct_seqs, repeated]


def exclude_small(fasta_file, minimum_length, variation=0):
    """Identify sequences smaller that a specified length value.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.
    minimum_length : int
        Sequences with a length value below this value are
        considered small.
    variation : float
        Accept sequences with length variation of up to
        minus (`minimum_length`*`variation`).

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
    small_seqids = sm.determine_small(fasta_file, minimum_length, variation)

    ss_lines = ['{0}: smaller than {1} chars'.format(seqid, minimum_length)
                for seqid in small_seqids]

    return [small_seqids, ss_lines]


def translate_sequences(sequence_ids, sequences_file, temp_directory,
                        translation_table, minimum_length, cpu_cores):
    """Translate DNA sequences.

    Returns information about sequences that are untranslatable
    and saves translatable DNA sequences and proteins to FASTA files.

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

    # create paths to files with protein sequences
    protein_template = fo.join_paths(temp_directory, ['translated_cds_{0}.fasta'])
    protein_files = [protein_template.format(i+1)
                     for i in range(len(translation_inputs))]

    # add common args to sublists
    common_args = [sequences_file, translation_table, minimum_length]
    translation_inputs = im.aggregate_iterables([translation_inputs,
                                                 protein_files])
    translation_inputs = im.multiprocessing_inputs(translation_inputs,
                                                   common_args,
                                                   sm.translate_coding_sequences)

    # translate sequences
    translation_results = mo.map_async_parallelizer(translation_inputs,
                                                    mo.function_helper,
                                                    cpu_cores,
                                                    show_progress=True)

    # concatenate files
    protein_file = fo.join_paths(temp_directory, ['translated_cds_concat.fasta'])
    protein_file = fo.concatenate_files(protein_files, protein_file)
    fo.remove_files(protein_files)

    # determine sequences that could not be translated
    untrans_lines = []
    untrans_seqids = []
    for res in translation_results:
        if len(res[0]) > 0:
            untrans_lines.extend(['{0}: {1}'.format(r[0], r[1])
                                  for r in res[0]])
            untrans_seqids.extend([r[0] for r in res[0]])

    return [protein_file, untrans_seqids, untrans_lines]


def cluster_sequences(sequences, word_size, window_size, clustering_sim,
                      representatives, grow_clusters, kmer_offset,
                      seq_num_cluster, temp_directory, cpu_cores,
                      divide, position):
    """Cluster sequences based on the proportion of shared minimizers.

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

    if divide is not None:
        # divide sequences into sublists
        # do not divide based on number of available cores as it may
        # lead to different results with different number of cores
        # divide into clusters with fixed number of sequences
        cluster_inputs = im.split_iterable(sorted_seqs, divide)
    else:
        cluster_inputs = [sorted_seqs]

    common_args = [word_size, window_size, clustering_sim,
                   representatives, grow_clusters, kmer_offset,
                   position, seq_num_cluster,
                   sc.clusterer]
    cluster_inputs = [[c, *common_args] for c in cluster_inputs]

    # cluster proteins in parallel
    clustering_results = mo.map_async_parallelizer(cluster_inputs,
                                                   mo.function_helper,
                                                   cpu_cores,
                                                   show_progress=True)

    # merge clusters
    clusters = [d[0] for d in clustering_results]
    clusters = im.merge_dictionaries(clusters)
    rep_sequences = [d[1] for d in clustering_results]
    rep_sequences = im.merge_dictionaries(rep_sequences)

    # perform clustering with representatives
    # this step does not run for AlleleCall
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

    # sort clusters
    clusters = {k: v for k, v in im.sort_iterable(clusters.items())}

    # write file with clustering results
    clusters_out = os.path.join(temp_directory, 'clusters.txt')
    sc.write_clusters(clusters, clusters_out)

    return clusters


def cluster_representative_filter(clusters, representative_filter,
                                  output_directory):
    """Exclude sequences highly similar to cluster representatives.

    After removing highly similar sequences, excludes clusters
    that are singletons (only contain the representative).

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

    # remove excluded seqids from clusters without high representative
    # similarity
    pruned_clusters = {k: [e for e in v if e[0] not in excluded_seqids]
                       for k, v in pruned_clusters.items()}

    # write file with pruning results
    pruned_out = os.path.join(output_directory, 'clusters.txt')
    sc.write_clusters(pruned_clusters, pruned_out)

    # identify singletons and exclude those clusters
    singletons = im.select_keys(pruned_clusters, 0)

    pruned_clusters = im.prune_dictionary(pruned_clusters, singletons)

    # determine number of sequences that still need to be evaluated
    # +1 to include representative
    clustered_sequences = sum([len(v)+1 for k, v in pruned_clusters.items()]) + len(singletons)

    # write list of excluded seqids to file
    excluded_outfile = os.path.join(output_directory, 'excluded.txt')
    fo.write_lines(excluded_seqids, excluded_outfile)

    return [pruned_clusters, excluded_seqids, singletons, clustered_sequences]


def cluster_intra_filter(clusters, sequences, word_size,
                         intra_filter, output_directory):
    """Exclude sequences that are highly similar to other clustered sequences.

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

    # remove excluded seqids from clusters without high intra-similarity
    pruned_clusters = {k: [e for e in v if e[0] not in intra_excluded]
                       for k, v in clusters.items()}

    # write excluded to file
    intrasim_out = os.path.join(output_directory, 'excluded.txt')
    sc.write_clusters(excluded_sims, intrasim_out)
    # write clusters to file
    intrasim_out = os.path.join(output_directory, 'clusters.txt')
    sc.write_clusters(pruned_clusters, intrasim_out)

    return [pruned_clusters, intra_excluded]


def blast_clusters(clusters, sequences, output_directory,
                   blastp_path, makeblastdb_path, cpu_cores,
                   only_rep=False):
    """Use BLAST to align sequences in the same clusters.

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
    # create FASTA file with sequences in clusters
    clustered_seqs_file = fo.join_paths(output_directory,
                                        ['clustered_sequences.fasta'])
    clustered_sequences = [[k]+[e[0] for e in v] for k, v in clusters.items()]
    clustered_sequences = im.flatten_list(clustered_sequences)
    # do not include duplicate identifiers
    clustered_sequences = list(set(clustered_sequences))
    fao.get_sequences_by_id(sequences, clustered_sequences,
                            clustered_seqs_file)

    # create FASTA file with replaced headers to avoid header
    # length limitation in BLAST
    integer_clusters = fo.join_paths(output_directory,
                                     ['clustered_sequences_integer_headers.fasta'])
    ids_dict = fao.integer_headers(clustered_seqs_file, integer_clusters)

    # create BLAST DB
    blast_db = fo.join_paths(output_directory,
                             ['clustered_sequences'])
    db_stderr = bw.make_blast_db(makeblastdb_path, integer_clusters,
                                 blast_db, 'prot')
    if len(db_stderr) > 0:
        sys.exit(db_stderr)

    blast_results_dir = os.path.join(output_directory,
                                     'BLAST_results')
    fo.create_directory(blast_results_dir)

    # create files with replaced sequence identifiers per cluster
    seqids_to_blast = sc.blast_seqids(clusters, blast_results_dir, ids_dict,
                                      only_rep)

    # distribute clusters per available cores
    process_num = 20 if cpu_cores <= 20 else cpu_cores
    splitted_seqids = mo.distribute_loci(seqids_to_blast, process_num, 'seqcount')

    common_args = [integer_clusters, blast_results_dir, blastp_path,
                   blast_db, only_rep, sc.cluster_blaster]

    splitted_seqids = [[s, *common_args] for s in splitted_seqids]

    # BLAST sequences in a cluster against every sequence in that cluster
    blast_results = mo.map_async_parallelizer(splitted_seqids,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    return [blast_results, ids_dict]


def compute_bsr(subject_score, query_score):
    """Compute the BLAST Score Ratio for an alignment between two sequences.

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


def determine_self_scores(fasta_file, output_directory, makeblastdb_path,
                          blast_path, db_type, blast_threads):
    """Compute the self-alignment raw score for sequences in a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file that contains the sequences.
    output_directory : str
        Path to the working directory.
    makeblastdb_path : str
        Path to the 'makeblastdb' executable.
    blast_path : str
        Path to the BLASTp/n executable.
    db_type : str
        Type of the database, nucleotide (nuc) or
        protein (prot).
    blast_threads : int
        Number of threads/cores used to run BLAST.

    Returns
    -------
    self_scores : dict
        Dictionary with sequence identifiers as keys and
        tuples with the sequence length and the raw score
        of the self-alignment as values.
    """
    # change identifiers to shorten and avoid BLAST error
    # related with sequence header length
    integer_fasta = fasta_file.replace('.fasta', '_integer_headers.fasta')
    ids_map = fao.integer_headers(fasta_file, integer_fasta, start=1, limit=5000)

    blast_db = fo.join_paths(output_directory,
                             [fo.file_basename(integer_fasta, False)])
    # will not work if file contains duplicates
    db_stderr = bw.make_blast_db(makeblastdb_path, integer_fasta,
                                 blast_db, db_type)

    if len(db_stderr) > 0:
        return db_stderr

    # split Fasta file to BLAST short sequences (<30aa) separately
    # only possible to have alleles <30aa with non-default schemas
    above_outfile, below_outfile = fao.split_seqlength(integer_fasta,
                                                       output_directory,
                                                       ct.BLAST_TASK_THRESHOLD['blastp'])

    if above_outfile is not None:
        # divide FASTA file into groups of 100 sequences to reduce
        # execution time for large sequence sets
        splitted_fastas = fao.split_seqcount(above_outfile[0], output_directory, 100)

        # create TXT with list of sequence identifiers
        seqids_files = []
        for f in splitted_fastas:
            seqids = list(f[1])
            seqids_file = fo.join_paths(output_directory, [fo.file_basename(f[0], False)])
            fo.write_lines(seqids, seqids_file)
            seqids_files.append(seqids_file)
    # this should not happen or be very rare, but just in case
    else:
        splitted_fastas = []

    # create directory to store results from final BLASTp
    final_blastp_dir = fo.join_paths(output_directory, ['BLAST_results'])
    fo.create_directory(final_blastp_dir)
    blast_outputs = ['{0}/{1}_blast_out.tsv'.format(final_blastp_dir,
                                                    fo.file_basename(file[0], False))
                     for file in splitted_fastas]

    # add common arguments to all sublists
    blast_inputs = [[blast_path, blast_db, file[0],
                     blast_outputs[i], 1, 1, seqids_files[i], 'blastp', None, ct.IGNORE_RAISED, bw.run_blast]
                    for i, file in enumerate(splitted_fastas)]

    # add file with short sequences
    if below_outfile is not None:
        seqids = list(below_outfile[1])
        seqids_file = fo.join_paths(output_directory, [fo.file_basename(below_outfile[0], False)])
        fo.write_lines(seqids, seqids_file)
        below_blastout = '{0}/{1}_blast_out.tsv'.format(final_blastp_dir,
                                                        fo.file_basename(below_outfile[0], False))
        blast_outputs.append(below_blastout)
        blast_inputs.append([blast_path, blast_db, below_outfile[0],
                             below_blastout, 1, 1, seqids_file,
                             'blastp-short', None, ct.IGNORE_RAISED, bw.run_blast])

    blast_stderr = mo.map_async_parallelizer(blast_inputs,
                                             mo.function_helper,
                                             blast_threads,
                                             show_progress=False)

    blast_stderr = im.flatten_list(blast_stderr)
    if len(blast_stderr) > 0:
        sys.exit(blast_stderr)

    # concatenate files with BLASTp results
    blast_output = fo.join_paths(final_blastp_dir, ['blast_out_concat.tsv'])
    blast_output = fo.concatenate_files(blast_outputs, blast_output)

    current_results = fo.read_tabular(blast_output)
    # get raw score and sequence length
    # multiply by 3 to get DNA sequence length and add 3 to count stop codon
    self_results = [line for line in current_results if line[0] == line[4]]
    self_scores = {}
    for line in self_results:
        # multiply by 3 to get DNA sequence length
        # add 3 to count stop codon
        dna_length = (int(line[3])*3)+3
        self_scores[ids_map[line[0]]] = (dna_length, float(line[6]))

    return self_scores
