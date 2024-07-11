#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains core functions used by chewBBACA's
modules.

Code documentation
------------------
"""


import os

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
    if ptf_path is not None:
        # Read training file to create GeneFinder object
        training_data = gp.read_training_file(ptf_path)
        # Create GeneFinder object based on training data
        gene_finder = gp.create_gene_finder(training_data, True, True, False)
    elif ptf_path is None and prodigal_mode == 'meta':
        # Create GeneFinder object to run in meta mode
        gene_finder = gp.create_gene_finder(None, True, True, True)
    else:
        gene_finder = None

    common_args = [output_directory, gene_finder, translation_table]

    # Divide inputs into equal number of sublists for maximum process
    # progress resolution
    pyrodigal_inputs = im.divide_list_into_n_chunks(list(fasta_files.items()),
                                                    len(fasta_files))

    # Add common arguments to all sublists
    pyrodigal_inputs = im.multiprocessing_inputs(pyrodigal_inputs,
                                                 common_args,
                                                 gp.predict_genome_genes)

    # Run Pyrodigal to predict genes
    # Need to use ThreadPool. Pyrodigal might hang when using Pool
    pyrodigal_results = mo.map_async_parallelizer(pyrodigal_inputs,
                                                  mo.function_helper,
                                                  cpu_cores,
                                                  show_progress=True,
                                                  pool_type='threadpool')

    # Determine if Pyrodigal predicted genes for all genomes
    failed = {line[0][0]: line[1]
              for line in pyrodigal_results
              if line[1] == 0
              or isinstance(line[1], str) is True}

    cds_counts = {line[0][1]: line[1]
                  for line in pyrodigal_results
                  if isinstance(line[1], int) is True}

    total_cds = sum([line[1]
                     for line in pyrodigal_results
                     if isinstance(line[1], int) is True])

    cds_fastas = [line[2] for line in pyrodigal_results if line[2] is not None]

    cds_hashes = {line[0][1]: line[3] for line in pyrodigal_results if line[3] is not None}

    # Merge dictionaries with info about CDSs closes to contig tips
    close_to_tip = [line[4] for line in pyrodigal_results if len(line[4]) > 0]
    close_to_tip = im.merge_dictionaries(close_to_tip)

    return [failed, total_cds, cds_fastas, cds_hashes, cds_counts, close_to_tip]


def merge_dna_dedup(dedup_files, ids_map):
    """Select distinct DNA sequences based on sequence deduplication results.

    Parameters
    ----------
    dedup_files : list
        List with the paths to the files with the results
        from sequence deduplication.
    ids_map : dict
        List with two dictionaries mapping input
        identifiers to integer identifiers and vice
        versa.

    Return
    ------
    merged_results : dict
        Dictionary with distinct sequence hashes as keys
        and lists with the protid and input integer identifier
        of the sequence chosen as representative followed by the
        list of other inputs that contain the sequence
        (e.g. 02bf9c4f695ce7...: [3567, 1, 2, 3, 10] ).
        The lists are encoded with polyline encoding.
    total_duplicated : dict
        Total number of times a sequence was duplicated.
    representative_seqids : list
        List with the sequence identifiers of the
        representative sequences.
    """
    merged_results = {}
    total_duplicated = 0
    representative_seqids = []
    for file in dedup_files:
        results = fo.pickle_loader(file)
        # Sequence hash and list with protid:inputID pairs
        for hashid, ids in results.items():
            stored_ids = merged_results.get(hashid, [])
            # Already have a representative seqid
            if len(stored_ids) > 0:
                stored_ids = im.polyline_decoding(stored_ids)
                # Get input integer identifiers
                new_ids = [ids[i] for i in range(1, len(ids), 2)]
                total_duplicated += (len(ids)/2)
            # New distinct sequence
            else:
                # Choose first sequence seqid as representative
                rep_seqid = '{0}-protein{1}'.format(ids_map[1][ids[1]], ids[0])
                representative_seqids.append(rep_seqid)
                new_ids = ids[0:2] + [ids[i] for i in range(3, len(ids), 2)]
                # Do not count representative as duplicate
                total_duplicated += (len(ids)/2) - 1

            merged_results[hashid] = im.polyline_encoding(stored_ids+new_ids)

    return [merged_results, total_duplicated, representative_seqids]


def merge_protein_dedup(dedup_files, ids_map):
    """Select distinct proteins based on sequence deduplication results.

    Parameters
    ----------
    dedup_files : list
        List with the paths to the files with the results
        from sequence deduplication.
    ids_map : dict
        List with two dictionaries mapping input
        identifiers to integer identifiers and vice
        versa.

    Return
    ------
    merged_results : dict
        Dictionary with distinct sequence hashes as keys
        and lists with protid and input integer identifiers
        pairs. The lists are encoded with polyline encoding.
    total_duplicated : dict
        Total number of times a sequence was duplicated.
    representative_seqids : list
        List with the sequence identifiers of the
        representative sequences.
    """
    merged_results = {}
    total_duplicated = 0
    representative_seqids = []
    for file in dedup_files:
        results = fo.pickle_loader(file)
        # Sequence hash and list with protid:inputID pairs
        for hashid, ids in results.items():
            stored_ids = merged_results.get(hashid, [])
            new_ids = ids
            # Already have a representative seqid
            if len(stored_ids) > 0:
                stored_ids = im.polyline_decoding(stored_ids)
                total_duplicated += (len(ids)/2)
            else:
                rep_seqid = '{0}-protein{1}'.format(ids_map[1][ids[1]], ids[0])
                representative_seqids.append(rep_seqid)
                total_duplicated += (len(ids)/2) - 1
            
            merged_results[hashid] = im.polyline_encoding(stored_ids+new_ids)

    return [merged_results, total_duplicated, representative_seqids]


def merge_get_seqids(dedup_files, ids_map):
    """Select distinct sequences based on sequence deduplication results.

    Parameters
    ----------
    dedup_files : list
        List with the paths to the files with the results
        from sequence deduplication.
    ids_map : dict
        List with two dictionaries mapping input
        identifiers to integer identifiers and vice
        versa.

    Return
    ------
    total_duplicated : dict
        Total number of times a sequence was duplicated.
    representative_seqids : list
        List with the sequence identifiers of the
        representative sequences.
    """
    total_duplicated = 0
    representative_seqids = []
    representative_hashes = set()
    for file in dedup_files:
        results = fo.pickle_loader(file)
        # Sequence hash and list with protid:inputID pairs
        for hashid, ids in results.items():
            if hashid not in representative_hashes:
                rep_seqid = '{0}-protein{1}'.format(ids_map[1][ids[1]], ids[0])
                representative_seqids.append(rep_seqid)
                representative_hashes.add(hashid)
                total_duplicated += (len(ids)/2) - 1
            else:
                total_duplicated += (len(ids)/2)

    return [None, total_duplicated, representative_seqids]


def exclude_duplicates(fasta_files, temp_directory, cpu_cores,
                       ids_map, protein=False, only_seqids=False):
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
    ids_map : list
        List with two dictionaries mapping input
        identifiers to integer identifiers and vice
        versa.
    protein : bool
        If sequences to deduplicate are proteins.
    only_seqids : bool
        Return only the list of distinct seqids.
        Do not return the hash table with encoded results.

    Returns
    -------
    A list with the following elements:
        distinct_seqids : list
            List with the sequence identifiers of distinct
            sequences.
        distinct_seqs : str
            Path to the FASTA file with distinct sequences.
    """
    # Define output files paths
    output_files = [fo.join_paths(temp_directory, [f'distinct_{str(i+1)}.fasta'])
                    for i, file in enumerate(fasta_files)]
    # Create input lists to pass to multiprocessing
    dedup_inputs = im.aggregate_iterables([fasta_files, output_files])
    dedup_inputs = im.multiprocessing_inputs(dedup_inputs, [ids_map[0]],
                                             sm.determine_distinct)

    # Determine set of distinct sequences per file
    dedup_results = mo.map_async_parallelizer(dedup_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=False)

    if only_seqids is True:
        merge_results = merge_get_seqids(dedup_results, ids_map)
    else:
        if protein is False:
            merge_results = merge_dna_dedup(dedup_results, ids_map)  
        else:
            merge_results = merge_protein_dedup(dedup_results, ids_map)

    hash_table, duplicated_count, representative_seqids = merge_results

    if hash_table is not None:
        hash_table_file = fo.join_paths(temp_directory, ['distinct.hashtable'])
        fo.pickle_dumper(hash_table, hash_table_file)

    # Concatenate Fasta files from parallel processes
    dedup_outfiles = [files[1] for files in dedup_inputs]
    merge_fasta = fo.join_paths(temp_directory, ['concat.fasta'])
    merge_fasta = fo.concatenate_files(dedup_outfiles, merge_fasta)

    # Create index for concatenated Fasta and get representative sequences
    merge_fasta_index = fao.index_fasta(merge_fasta)
    representative_fasta = fo.join_paths(temp_directory, ['distinct.fasta'])
    fao.get_sequences_by_id(merge_fasta_index, representative_seqids,
                            representative_fasta, 50000)

    # Remove intermediate files
    fo.remove_files(dedup_outfiles+dedup_results+[merge_fasta])

    return [hash_table, representative_seqids,
            representative_fasta, duplicated_count]


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
    # Divide inputs into sublists
    translation_inputs = im.divide_list_into_n_chunks(sequence_ids, cpu_cores)

    # Create paths to files with protein sequences
    protein_template = fo.join_paths(temp_directory, ['translated_{0}.fasta'])
    protein_files = [protein_template.format(i+1)
                     for i in range(len(translation_inputs))]

    # Add common args to sublists
    common_args = [sequences_file, translation_table, minimum_length]
    translation_inputs = im.aggregate_iterables([translation_inputs,
                                                 protein_files])
    translation_inputs = im.multiprocessing_inputs(translation_inputs,
                                                   common_args,
                                                   sm.translate_coding_sequences)

    # Translate sequences
    translation_results = mo.map_async_parallelizer(translation_inputs,
                                                    mo.function_helper,
                                                    cpu_cores,
                                                    show_progress=True)

    # Concatenate files
    protein_file = fo.join_paths(temp_directory, ['translated.fasta'])
    protein_file = fo.concatenate_files(protein_files, protein_file)
    fo.remove_files(protein_files)

    # Determine sequences that could not be translated
    untrans_lines = []
    untrans_seqids = []
    for res in translation_results:
        if len(res[0]) > 0:
            untrans_lines.extend([f'{r[0]}: {r[1]}' for r in res[0]])
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
    divide : int
        Create smaller groups of sequences that can
        be processed in parallel.
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
    # Sort sequences by length
    sorted_seqs = {k: v for k, v in im.sort_iterable(sequences.items(),
                                                     sort_key=lambda x: len(x[1]),
                                                     reverse=True)}

    if divide is not None:
        # Divide input sequences into smaller groups
        cluster_inputs = im.split_iterable(sorted_seqs, divide)
    else:
        cluster_inputs = [sorted_seqs]

    common_args = [word_size, window_size, clustering_sim,
                   representatives, grow_clusters, kmer_offset,
                   position, seq_num_cluster,
                   sc.clusterer]
    cluster_inputs = [[c, *common_args] for c in cluster_inputs]

    # Cluster proteins
    clustering_results = mo.map_async_parallelizer(cluster_inputs,
                                                   mo.function_helper,
                                                   cpu_cores,
                                                   show_progress=True)

    # Merge clusters
    clusters = [d[0] for d in clustering_results]
    clusters = im.merge_dictionaries(clusters)
    rep_sequences = [d[1] for d in clustering_results]
    rep_sequences = im.merge_dictionaries(rep_sequences)

    # Perform clustering with representatives
    # This step does not run for AlleleCall as we do
    # not need to compare schema representatives
    # It runs for CreateSchema because the representatives
    # from each cluster group were not compared
    # No need to run this step if input sequences were not divided into smaller groups
    # AlleleCall does not select new representatives, so this does not run
    if len(cluster_inputs) > 1 and len(rep_sequences) > 0:
        # Cluster representatives
        rep_clusters = sc.clusterer(rep_sequences, word_size,
                                    window_size, clustering_sim,
                                    representatives, grow_clusters,
                                    kmer_offset, position,
                                    seq_num_cluster)

        merged_clusters = {}
        for k, v in rep_clusters[0].items():
            # Merge clusters whose representatives are similar
            for n in v:
                # Representatives from other clusters are added with
                # similarity score against new representative
                # Clustered sequences from other clusters are added
                # with similarity score against their representative
                add_seqids = [n] + [s for s in clusters[n[0]] if s[0] != n[0]]
                merged_clusters.setdefault(k, []).extend(add_seqids)

        clusters = merged_clusters

    # Sort clusters
    clusters = {k: v for k, v in im.sort_iterable(clusters.items())}

    # Write file with clustering results
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
                   blastdb_aliastool_path, only_rep=False):
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
        Path to the FASTA file that contains the sequences in
        the clusters.
    output_directory : str
        Path to the directory where the clustering results
        will be saved to.
    blastp_path : str
        Path to the `BLASTp` executable.
    makeblastdb_path : str
        Path to the `makeblastdb` executable.
    cpu_cores : int
        Number of BLASTp processes to run in parallel.
	blastdb_aliastool_path : str
        Path to the blastalias_tool executable to convert
		seqid files to binary format.
    only_rep

    Returns
    -------
    blast_results : list
        List with paths to the files with BLASTp results
        (one file per cluster).
    """
    # Create directory to store BLASTp database
    blast_db_dir = fo.join_paths(output_directory, ['BLASTp_db'])
    fo.create_directory(blast_db_dir)
    # Create BLAST DB
    blast_db = fo.join_paths(blast_db_dir, ['distinct_proteins'])
    db_std = bw.make_blast_db(makeblastdb_path, sequences, blast_db, 'prot')

    blastp_results_dir = os.path.join(output_directory, 'BLASTp_outfiles')
    fo.create_directory(blastp_results_dir)

    # Create TXT files with the list of sequences per cluster
    seqids_to_blast = sc.blast_seqids(clusters, blastp_results_dir, only_rep)

    # Distribute clusters per available cores
    process_num = 20 if cpu_cores <= 20 else cpu_cores
    splitted_seqids = mo.distribute_loci(seqids_to_blast, process_num, 'seqcount')
    common_args = [sequences, blastp_results_dir, blastp_path,
                   blast_db, blastdb_aliastool_path, only_rep, sc.cluster_blaster]
    splitted_seqids = [[s, *common_args] for s in splitted_seqids]

    # BLAST sequences in a cluster against every sequence in that cluster
    blast_results = mo.map_async_parallelizer(splitted_seqids,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    return [blast_results, blastp_results_dir]


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
                          blast_path, db_type, blast_threads, blastdb_aliastool_path):
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
	blastdb_aliastool_path : str
        Path to the blastalias_tool executable to convert
		seqid files to binary format.

    Returns
    -------
    self_scores : dict
        Dictionary with sequence identifiers as keys and
        tuples with the sequence length and the raw score
        of the self-alignment as values.
    """
    blast_db_dir = fo.join_paths(output_directory, ['BLASTp_db'])
    fo.create_directory(blast_db_dir)
    blast_db = fo.join_paths(blast_db_dir, [fo.file_basename(fasta_file, False)])
    # Will not work if file contains duplicates
    db_std = bw.make_blast_db(makeblastdb_path, fasta_file, blast_db, db_type)

    # Split Fasta file to BLAST short sequences (<30aa) separately
    # only possible to have alleles <30aa with non-default schemas
    above_outfile, below_outfile = fao.split_seqlength(fasta_file,
                                                       output_directory,
                                                       ct.BLAST_TASK_THRESHOLD['blastp'])

    total_seqids = above_outfile[0] + below_outfile[0]
    all_seqids = above_outfile[1] + below_outfile[1]

    if above_outfile[0] > 0:
        above_outdir = fo.join_paths(output_directory, ['above'])
        fo.create_directory(above_outdir)
        # Divide FASTA file into groups of 100 sequences to reduce
        # execution time for large sequence sets
        split_fastas = fao.split_seqcount(above_outfile[2], above_outdir, 100)

        # Create TXT with list of sequence identifiers
        seqids_files = []
        for f in split_fastas:
            seqids = list(f[1])
            seqids_file = fo.join_paths(above_outdir, [fo.file_basename(f[0], False)])
            fo.write_lines(seqids, seqids_file)
            seqids_files.append(seqids_file)

		binary_seqid_files = []
		for file in seqids_files:
			binary_file = f'{file}.bin'
			blast_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
													file,
													binary_file)
			binary_seqid_files.append(binary_file)
		seqids_files = binary_seqid_files
    # This should not happen or be very rare, but just in case
    else:
        split_fastas = []

    # Create directory to store results from final BLASTp
    final_blastp_dir = fo.join_paths(output_directory, ['BLASTp_results'])
    fo.create_directory(final_blastp_dir)
    blast_outputs = ['{0}/{1}_blast_out.tsv'.format(final_blastp_dir,
                                                    fo.file_basename(file[0], False))
                     for file in split_fastas]

    # Add common arguments to all sublists
    blast_inputs = [[blast_path, blast_db, file[0],
                     blast_outputs[i], 1, 1,
                     seqids_files[i], 'blastp', None,
                     None, bw.run_blast]
                    for i, file in enumerate(split_fastas)]

    # Add file with short sequences
    if below_outfile[0] > 0:
        seqids = below_outfile[1]
        seqids_file = fo.join_paths(output_directory, [fo.file_basename(below_outfile[2], False)])
        fo.write_lines(seqids, seqids_file)

		binary_file = f'{seqids_file}.bin'
		blast_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
												seqids_file,
												binary_file)
		seqids_file = binary_file

        below_blastout = '{0}/{1}_blastout.tsv'.format(final_blastp_dir,
                                                        fo.file_basename(below_outfile[2], False))
        blast_outputs.append(below_blastout)
        blast_inputs.append([blast_path, blast_db, below_outfile[2],
                             below_blastout, 1, 1,
                             seqids_file, 'blastp-short', None,
                             None, bw.run_blast])

    blast_results = mo.map_async_parallelizer(blast_inputs,
                                              mo.function_helper,
                                              blast_threads,
                                              show_progress=False)

    # Concatenate files with BLASTp results
    blast_output = fo.join_paths(final_blastp_dir, ['blastout_concat.tsv'])
    blast_output = fo.concatenate_files(blast_outputs, blast_output)

    current_results = fo.read_tabular(blast_output)
    # Get raw score and sequence length
    # multiply by 3 to get DNA sequence length and add 3 to count stop codon
    self_results = [line for line in current_results if line[0] == line[4]]
    self_scores = {}
    for line in self_results:
        # multiply by 3 to get DNA sequence length
        # add 3 to count stop codon
        dna_length = (int(line[3])*3)+3
        self_scores[line[0]] = (dna_length, float(line[6]))

    # Determine if we got the score for all representatives
    if len(self_scores) < total_seqids:
        missing = [seqid
                   for seqid in all_seqids
                   if seqid not in [line[0] for line in self_results]]
        # Index FASTA file
        concat_reps_index = fao.index_fasta(fasta_file)
        # Get all representatives that do not have a self-score
        for seqid in missing:
            current_rep = concat_reps_index[seqid]
            record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [current_rep.id, str(current_rep.seq)])
            rep_file = fo.join_paths(output_directory, [f'{current_rep.id}_solo.fasta'])
            fo.write_lines([record], rep_file)
            # Create file with representative seqid to only compare against self
            id_file = fo.join_paths(output_directory, [f'{current_rep.id}_ids.txt'])
            fo.write_lines([current_rep.id], id_file)

			binary_file = f'{id_file}.bin'
			blast_std = bw.run_blastdb_aliastool(blastdb_aliastool_path,
													id_file,
													binary_file)
			id_file = binary_file
            rep_blastout = fo.join_paths(output_directory, [f'{current_rep.id}_blastout.tsv'])
            # Cannot get self-alignemnt for some sequences if composition-based stats is enabled
            blast_std = bw.run_blast(blast_path, blast_db, rep_file,
                                     rep_blastout, 1, 1,
                                     id_file, 'blastp', None, 0)
            rep_results = fo.read_tabular(rep_blastout)
            if len(rep_results) > 0:
                dna_length = (int(rep_results[0][3])*3)+3
                self_scores[rep_results[0][0]] = (dna_length, float(rep_results[0][6]))
            else:
                print('Could not determine the self-alignment raw '
                      f'score for {rep_results[0][0]}')

    return self_scores


def write_coordinates_file(coordinates_file, output_file):
    """Write genome CDS coordinates to a TSV file.

    Parameters
    ----------
    coordinates_file : str
        Path to the pickle file that contains data about
        the CDSs coordinates.
    output_file : str
        Path to the output TSV file.
    """
    data = fo.pickle_loader(coordinates_file)
    lines = [coords for h, coords in data[0].items()]
    lines = im.flatten_list(lines)
    lines = ['\t'.join(line) for line in lines]
    fo.write_lines(lines, output_file)
