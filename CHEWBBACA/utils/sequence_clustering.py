#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to sequence clustering
based on k-mers.

Code documentation
------------------
"""


import os
from collections import Counter

from Bio import SeqIO

try:
    from utils import (file_operations as fo,
                       iterables_manipulation as im,
                       blast_wrapper as bw,
                       fasta_operations as fao,
                       constants as ct)
except:
    from CHEWBBACA.utils import (file_operations as fo,
                                 iterables_manipulation as im,
                                 blast_wrapper as bw,
                                 fasta_operations as fao,
                                 constants as ct)


def select_representatives(kmers, reps_groups, clustering_sim):
    """ Determines the set of clusters that a sequence
        can be added to based on the decimal proportion
        of shared distinct kmers.

        Parameters
        ----------
        kmers : list or set
            Set of kmers determined by decomposing a single
            sequence.
        reps_groups : dict
            Dictionary with kmers as keys and sequence
            identifiers of sequences that contain that
            kmer as values.
        clustering_sim : float
            Sequences are added to clusters if they
            share a minimum decimal proportion of
            distinct kmers with a cluster representative.

        Returns
        -------
        selected_reps : list
            List with a tuple per cluster/representative
            that the sequence can be added to. Each tuple
            has the identifier of the cluster representative
            and the decimal proportion of shared distinct
            kmers.
    """

    current_reps = [reps_groups[k] for k in kmers if k in reps_groups]
    current_reps = im.flatten_list(current_reps)

    # count number of kmer hits per representative
    counts = Counter(current_reps)
    selected_reps = [(k, v/len(kmers))
                     for k, v in counts.items()
                     if v/len(kmers) >= clustering_sim]

    # sort by identifier and then by similarity to always get same order
    selected_reps = sorted(selected_reps, key=lambda x: x[0])
    selected_reps = sorted(selected_reps, key=lambda x: x[1], reverse=True)

    return selected_reps


def pick_excluded(query, match, query_seq, match_seq):
    """ Determines which sequence should be excluded
        based on sequence length.

        Parameters
        ----------
        query : str
            Query sequence identifier.
        match : tup
            Tuple with a sequence identifier and
            the decimal proportion of shared
            distinct kmers with the query sequence.
        query_seq : str
            Query sequence.
        match_seq : str
            Match sequence.

        Returns
        -------
        picked : tup
            Tuple with the sequence identifier of the
            sequence that should be excluded, the sequence
            identifier of the sequence that matched with
            the excluded sequence and the decimal proportion
            of shared distinct kmers between both sequences.
    """

    # if query sequence is longer, keep it and exclude candidate
    if len(query_seq) >= len(match_seq):
        picked = (match[0], query, match[1])
    # otherwise, exclude query sequence
    elif len(match_seq) > len(query_seq):
        picked = (query, match[0], match[1])

    return picked


def intra_cluster_sim(clusters, sequences, word_size, intra_filter):
    """ Determines the percentage of shared kmers/minimizers
        between sequences in the same cluster and excludes
        sequences that are similar to other sequences in the
        cluster.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        sequences : dict
            Dictionary with sequence identifiers as keys
            and sequences as values.
        word_size : int
            Value k for the kmer size.
        intra_filter : float
            Similarity threshold value. If two sequences in
            the same cluster have a similarity value equal
            or greater to this value, the shorter sequence
            will be excluded from the cluster.

        Returns
        -------
        excluded_seqids : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with the sequence identifiers of sequences
            that were excluded from the cluster as values.
        excluded_sims : dict
            A dictionary with the identifiers of sequences
            that are cluster representatives as keys
            and a list with tuples as values. Each tuple
            contains the sequence identifier and the
            similarity value for a match that led to
            an exclusion.
    """

    excluded_seqids = {}
    excluded_sims = {}
    for representative, clustered in clusters.items():
        # get identifiers of sequences in the cluster
        clustered_ids = [c[0] for c in clustered]
        clustered_seqs = {seqid: sequences[seqid] for seqid in clustered_ids}

        # create kmer index
        kmers_mapping, cluster_kmers = im.kmer_index(clustered_seqs, word_size)

        excluded = []
        similarities = []
        # for each sequence in the cluster
        for seqid, kmers in cluster_kmers.items():
            if seqid not in excluded:
                query_kmers = kmers

                # select sequences with same kmers
                sims = select_representatives(query_kmers,
                                              kmers_mapping,
                                              intra_filter)

                if len(sims) > 1:
                    # exclude current sequence
                    candidates = [s for s in sims if s[0] != seqid]
                    for c in candidates:
                        picked = pick_excluded(seqid, c,
                                               clustered_seqs[seqid],
                                               clustered_seqs[c[0]])
                        similarities.append(picked)
                        excluded.append(picked[0])

        # convert into set first to remove possible duplicates
        excluded_seqids[representative] = list(set(excluded))
        excluded_sims[representative] = similarities

    return [excluded_seqids, excluded_sims]


def minimizer_clustering(sorted_sequences, word_size, window_size, position,
                         offset, clusters, reps_sequences, reps_groups,
                         seq_num_cluster, clustering_sim):
    """ Cluster sequences based on the decimal proportion of
        shared distinct minimizers.

        Parameters
        ----------
        sorted_sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values. Sorted by decreasing sequence
            length.
        word_size : int
            Value k for the kmer size.
        window_size : int
            Window size used to determine minimizers.
        position : bool
            If minimizer sequence position should be stored.
        offset : int
            Value to indicate offset of consecutive kmers.
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the decimal proportion of shared
            distinct minimizers and the length of the clustered
            sequence. This dictionary should be empty at the
            start of the clustering step during the CreateSchema
            process. For the AlleleCall process, the dictionary
            should contain the identifiers of the loci
            representatives.
        reps_sequences : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            their sequences as values. This dictionary should
            be empty at the start of the clustering step
            during the CreateSchema process. For the AlleleCall
            process, the dictionary should contain the
            identifiers of the loci representatives as keys
            and their sequences as values.
        reps_groups : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values. This dictionary should be empty at the
            start of the clustering step during the CreateSchema
            process. For the AlleleCall process, the dictionary
            should contain all kmers contained in the set of
            the schema's representatives sequences.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.

        Returns
        -------
        A list with the following elements:
            clusters : dict
                Dictionary with the identifiers of sequences
                that are clusters representatives as keys and
                a list with tuples as values. Each tuple has
                the identifier of a sequence that was added to
                the cluster, the decimal proportion of shared
                distinct minimizers and the length of the clustered
                sequence.
            reps_sequences : dict
                Dictionary with the identifiers of sequences
                that are cluster representatives as keys and
                their sequences as values.
            reps_groups : dict
                Dictionary with kmers as keys and a list with
                identifiers of sequences that contain that kmer
                as values.
    """

    for protid, protein in sorted_sequences.items():

        minimizers = im.determine_minimizers(protein, window_size,
                                             word_size, offset=offset,
                                             position=position)

        distinct_minimizers = set(minimizers)

        selected_reps = select_representatives(distinct_minimizers,
                                               reps_groups,
                                               clustering_sim)

        top = (len(selected_reps)
               if len(selected_reps) < seq_num_cluster
               else seq_num_cluster)

        # sort to get most similar at index 0
        if len(selected_reps) > 0:
            for i in range(0, top):
                clusters[selected_reps[i][0]].append((protid,
                                                      selected_reps[i][1],
                                                      len(protein),
                                                      len(minimizers),
                                                      len(distinct_minimizers)))
        else:
            for k in distinct_minimizers:
                reps_groups.setdefault(k, []).append(protid)

            clusters[protid] = [(protid, 1.0, len(protein),
                                len(minimizers), len(distinct_minimizers))]
            reps_sequences[protid] = protein

    return [clusters, reps_sequences, reps_groups]


def clusterer(sorted_sequences, word_size, window_size,
              clustering_sim, representatives, grow,
              offset, position, seq_num_cluster):
    """ Cluster sequences based on the decimal proportion of
        shared distinct minimizers.

        Parameters
        ----------
        sorted_sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values. Sorted by decreasing sequence
            length.
        word_size : int
            Value k for the kmer size.
        window_size : int
            Window size used to determine minimizers.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.
        representatives : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values.
        grow : bool
            If it is allowed to create new clusters.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If minimizer sequence position should be stored.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.

        Returns
        -------
        A list with the following elements:
            clusters : dict
                Dictionary with the identifiers of sequences
                that are cluster representatives as keys and
                a list with tuples as values. Each tuple has
                the identifier of a sequence that was added to
                the cluster, the decimal proportion of shared
                distinct minimizers and the length of the clustered
                sequence.
            reps_sequences : dict
                Dictionary with the identifiers of sequences
                that are cluster representatives as keys and
                their sequences as values.
    """

    clusters = {}
    reps_sequences = {}
    if representatives is None:
        reps_groups = {}
    else:
        reps_groups = representatives
        clusters_ids = set(im.flatten_list(list(reps_groups.values())))
        clusters = {rep: [] for rep in clusters_ids}

    cluster_results = minimizer_clustering(sorted_sequences, word_size,
                                           window_size, position,
                                           offset, clusters,
                                           reps_sequences, reps_groups,
                                           seq_num_cluster, clustering_sim)

    return cluster_results[0:2]


def write_clusters(clusters, outfile):
    """ Writes information about clusters to file.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the decimal proportion of shared
            distinct kmers/minimizers and the length of the
            clustered sequence.
        outfile : str
            Path to the file that will be created to save
            information about clusters.
    """

    cluster_lines = []
    for rep, seqids in clusters.items():
        current_cluster = []
        current_cluster.append('>{0}'.format(rep))
        clustered = [', '.join(['{}']*len(s)).format(*s)
                     for s in seqids]
        current_cluster.extend(clustered)
        cluster_lines.append(current_cluster)

    # sort by number of lines to get clusters with more sequences first
    cluster_lines = im.sort_data(cluster_lines,
                                 sort_key=lambda x: len(x), reverse=True)
    cluster_lines = im.flatten_list(cluster_lines)
    cluster_text = im.join_list(cluster_lines, '\n')

    fo.write_to_file(cluster_text, outfile, 'w', '\n')


def representative_pruner(clusters, sim_cutoff):
    """ Removes sequences from clusters based on a similarity
        threshold.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        sim_cutoff : float
            Similarity threshold value. Sequences with
            equal or greater similarity value are excluded
            from clusters.

        Returns
        -------
        A list with the following elements:
            pruned_clusters : dict
                Input dictionary without values for the
                sequences that had similarity values
                equal or greater that defined threshold.
            excluded : list
                List with a list per cluster. Each list
                has the sequences that were excluded from
                a cluster.
    """

    excluded = []
    pruned_clusters = {}
    for rep, seqids in clusters.items():
        pruned_clusters[rep] = [seqid
                                for seqid in seqids
                                if seqid[1] < sim_cutoff]
        excluded.extend([seqid
                         for seqid in seqids
                         if seqid[1] >= sim_cutoff and seqid[0] != rep])

    return [pruned_clusters, excluded]


def cluster_blaster(seqids, sequences, output_directory,
                    blast_path, blastdb_path):
    """ Aligns sequences in the same cluster with BLAST.

        Parameters
        ----------
        seqids : list
            List with cluster identifiers.
        sequences : str
            Path to the FASTA file with the protein sequences
            in all clusters.
        output_directory : str
            Path to the directory where FASTA files with
            the sequences in each cluster and files with
            BLAST results will be written to.
        blast_path : str
            Path to BLAST executables
        blastdb_path : str
            Path to a BLAST database.

        Returns
        -------
        out_files : list
            List with the paths to the files with the BLAST
            results for each cluster.
    """

    indexed_fasta = SeqIO.index(sequences, 'fasta')

    out_files = []
    for cluster in seqids:

        cluster_id = cluster
        ids_file = os.path.join(output_directory,
                                '{0}_ids.txt'.format(cluster_id))

        with open(ids_file, 'r') as clstr:
            cluster_ids = [l.strip() for l in clstr.readlines()]

        fasta_file = os.path.join(output_directory,
                                  '{0}_protein.fasta'.format(cluster_id))
        # create file with protein sequences
        fao.get_sequences_by_id(indexed_fasta, cluster_ids, fasta_file)

        blast_output = os.path.join(output_directory,
                                    '{0}_blast_out.tsv'.format(cluster_id))

        # Use subprocess to capture errors and warnings
        stderr = bw.run_blast(blast_path, blastdb_path, fasta_file,
                              blast_output, 1, 1, ids_file,
                              ignore=ct.IGNORE_RAISED)

        if len(stderr) > 0:
            raise ValueError('\n'.join(stderr))

        out_files.append(blast_output)

    return out_files


def blast_inputs(clusters, output_directory, ids_dict):
    """ Creates files with the identifiers of the sequences
        in each cluster.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of cluster
            representatives as keys and a list with tuples
            as values (each tuple has the identifier of a
            sequence that is in the cluster, the decimal
            proportion of shared minimizers and the length
            of that sequence).
        output_directory : str
            Path to the directory where files with identifiers
            will be created.
        ids_dict : dict
            Dictionary that maps sequence identifiers to
            shorter and unique identifiers that will be
            saved in the files and used as sequence
            identifiers during BLAST to avoid errors
            related with sequence headers/identifiers
            that exceed length limit allowed by BLAST.

        Returns
        -------
        ids_to_blast : list
            List with the identifiers of all clusters.
    """

    rev_ids = {v: k for k, v in ids_dict.items()}

    ids_to_blast = []
    for rep in clusters:

        cluster_file = os.path.join(output_directory,
                                    '{0}_ids.txt'.format(rev_ids[rep]))
        cluster_ids = [rev_ids[rep]] + [rev_ids[seqid[0]] for seqid in clusters[rep]]
        cluster_lines = im.join_list(cluster_ids, '\n')
        fo.write_to_file(cluster_lines, cluster_file, 'w', '')
        ids_to_blast.append((rev_ids[rep], len(cluster_ids)))

    return ids_to_blast
