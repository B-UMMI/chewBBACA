#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

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


def intra_cluster_sim(clusters, sequences, word_size, intra_filter):
    """ Determines the percentage of shared kmers/minimizers
        between sequences in the same cluster and excludes
        sequences from a clusters sequences that are similar
        to other sequences in the cluster.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
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
            or greater to this value, one of the sequences
            will be excluded from the cluster.

        Returns
        -------
        excluded_dict : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            lists as values. Each list has two elements:
            a list with the identifiers of the sequences
            that were excluded from the clusters and a
            dictionary with sequences identifiers as keys
            and tuples with sequence identifiers and the
            similarity value for each match that led to
            an eclusion.
    """

    excluded_seqids = {}
    excluded_sims = {}
    for representative, clustered in clusters.items():
        # get identifiers of sequences in the cluster
        clustered_ids = [c[0] for c in clustered]
        clustered_seqs = {seqid: sequences[seqid] for seqid in clustered_ids}

        # get all kmers per sequence
        kmers_mapping = {}
        cluster_kmers = {}
        for seqid, seq in clustered_seqs.items():
            minimizers = im.determine_minimizers(seq, word_size, word_size, position=False)
            kmers = set(minimizers)

            # dict with sequence indentifiers and kmers
            cluster_kmers[seqid] = kmers

            # create dict with kmers as keys and list
            # of sequences with given kmers as values
            for kmer in kmers:
                kmers_mapping.setdefault(kmer, []).append(seqid)

        excluded = []
        similarities = []
        # for each sequence in the cluster
        for seqid, kmers in cluster_kmers.items():
            if seqid not in excluded:
                query_kmers = kmers
                # determine sequences that also have the same kmers
                current_reps = [kmers_mapping[kmer] for kmer in query_kmers if kmer in kmers_mapping]
                current_reps = im.flatten_list(current_reps)

                # count number of common kmers with other sequences
                counts = Counter(current_reps)
                # determine sequences that are equal or above a similarty threshold
                current_reps = [(s, v/len(kmers)) for s, v in counts.items() if v/len(kmers) >= intra_filter]
                # sort to get most similar first
                sims = sorted(current_reps, key=lambda x: x[1], reverse=True)

                if len(sims) > 1:
                    # exclude current sequence
                    candidates = [s for s in sims if s[0] != seqid]
                    for c in candidates:
                        # if query sequence is longer, keep it and exclude candidate
                        if len(clustered_seqs[seqid]) >= len(clustered_seqs[c[0]]):
                            similarities.append((c[0], seqid, c[1]))
                            excluded.append(c[0])
                        # otherwise, exclude query sequence
                        elif len(clustered_seqs[c[0]]) > len(clustered_seqs[seqid]):
                            similarities.append(([seqid, c[0], c[1]]))
                            excluded.append(seqid)

        # convert into set first to remove possible duplicates
        excluded_seqids[representative] = list(set(excluded))
        excluded_sims[representative] = similarities

    return [excluded_seqids, excluded_sims]


def select_minimizer_cluster(kmers, reps_groups, clustering_sim):
    """
    """

    current_reps = [reps_groups[k] for k in kmers if k in reps_groups]
    current_reps = im.flatten_list(current_reps)

    # count number of kmer hits per representative
    counts = Counter(current_reps)
    selected_reps = [(k, v/len(kmers))
                     for k, v in counts.items()
                     if v/len(kmers) >= clustering_sim]

    selected_reps = sorted(selected_reps, key=lambda x: x[0])
    selected_reps = sorted(selected_reps, key=lambda x: x[1], reverse=True)

    return selected_reps


def minimizer_clustering(sorted_sequences, word_size, window_size, position,
                         offset, clusters, reps_sequences, reps_groups,
                         seq_num_cluster, clustering_sim):
    """
    """

    for protid, protein in sorted_sequences.items():

        minimizers = im.determine_minimizers(protein, window_size,
                                             word_size, offset=offset,
                                             position=position)
        kmers = list(set(minimizers))

        selected_reps = select_minimizer_cluster(kmers, reps_groups,
                                                 clustering_sim)

        top = (len(selected_reps)
               if len(selected_reps) < seq_num_cluster
               else seq_num_cluster)

        # sort to get most similar at index 0
        if len(selected_reps) > 0:
            for i in range(0, top):
                clusters[selected_reps[i][0]].append((protid,
                                                      selected_reps[i][1],
                                                      len(protein)))
        else:
            for k in kmers:
                reps_groups.setdefault(k, []).append(protid)

            clusters[protid] = [(protid, 1.0, len(protein))]
            reps_sequences[protid] = protein

    return [clusters, reps_sequences, reps_groups]


##sorted_sequences = {'GCA_000007265-protein1917': 'MNTNIPSQFLPFFIPLILLQVILIIIALLKLRKITKTNYLSKPVWVLIILFVNLLGPIAFLSLEGKNA',
##                    'GCA_000007265-protein1969': 'MNTEIPNQLLPFLVPLIMLQGALIIISLVKLSKLKMTKHLSKPVWFLVIIFLNIIGPIAFLILEGNNE'}
#proteins_file = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/sagalactiae_schema/temp/postCluster.fasta'
#proteins = {rec.id: str(rec.seq) for rec in SeqIO.parse(proteins_file, 'fasta')}
#sorted_sequences = proteins
#word_size = 5
#window_size = 5
#position = True
#clustering_sim = 0.2
#clustering_mode = 'greedy'
#representatives = None
#grow = True
#offset = 1
#minimizer = True
#seq_num_cluster = 1
#temp_directory = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/sagalactiae_schema/temp'
#cpu_cores = 6
#clustering_type = 'minstrobes'
#num_strobes = 1
##clustering_type = 'minimizers'


def cluster_sequences(sorted_sequences, word_size, window_size,
                      clustering_sim, representatives, grow,
                      offset, position, seq_num_cluster):
    """ Cluster sequences based on shared percentage of kmers/minimizers.

        Parameters
        ----------
        sorted_sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values. Sorted by decreasing sequence
            length.
        word_size : int
            Value k for the kmer size.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.
        mode : str
            Clustering mode.
        representatives : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values.
        grow : bool
            If it is allowed to create new clusters.
        offset : int
            Value to indicate offset of consecutive kmers.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.

        Returns
        -------
        A list with the following elements:
            clusters : dict
                Dictionary with the identifiers of sequences
                that are clusters representatives as keys and
                a list with tuples as values. Each tuple has
                the identifier of a sequence that was added to
                the cluster, the percentage of shared
                kmers/minimizers and the length of the clustered
                sequence.
            reps_sequences : dict
                Dictionary with the identifiers of sequences
                that are clusters representatives as keys and
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
    """ Writes clusters to file.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        outfile : str
            Path to the file that will be
            created to save clusters.
    """

    cluster_lines = []
    for rep, seqids in clusters.items():
        current_cluster = []
        current_cluster.append('>{0}'.format(rep))
        clustered = ['\t{0}, {1}, {2}'.format(s[0], s[1], s[2])
                     for s in seqids]
        current_cluster.extend(clustered)
        cluster_lines.append(current_cluster)

    # sort by number of lines to get clusters with more sequences first
    cluster_lines = im.sort_data(cluster_lines, sort_key=lambda x: len(x), reverse=True)
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
        inputs : list
            List with clusters identifiers, the path to
            BLAST executables, the path to a BLAST database,
            the path to the directory that contains files with
            the identifiers of the sequences in each cluster
            and where FASTA files with the sequences in each
            cluster and files with BLAST results will be
            written to and the path to the FASTA file with
            the protein sequences in all clusters.

        Returns
        -------
        out_files : list
            List with the paths to the files with BLAST
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
            sequence that is in the cluster, the percentage
            of shared minimizers and the length os that
            sequence).
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
