#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 23:57:49 2020

@author: rfm
"""


import os
import itertools
from collections import Counter
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Args:
            strDNA (str): string representing a DNA sequence.

        Returns:
            revC_dna (str): the reverse complement of the DNA sequence, without
            lowercase letters.

        Example:
            >>> reverse_complement('ATCGgcaNn')
            'NNTGCCGAT'
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A',
                       'n': 'N', 'N': 'N'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    bases = [base_complement[base] for base in bases]

    complement_strand = ''.join(bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def reverse_str(string):
    """ Reverse character order in input string.

        Args:
            string (str): string to be reversed.

        Returns:
            revstr (str): reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.

        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Tries to translate an input DNA sequence in specified orientation
        and stores exceptions when the input sequence cannot be translated.

        Args:
            dna_sequence (str): string representing a DNA sequence.
            method (str): a string specifying the way the sequence will
            be oriented to attempt translation.
            table_id (int): translation table identifier.

        Returns:
            List with following elements if translation is successful:
                protseq (str): string representing the translated DNA sequence.
                myseq (str): string representing the DNA sequence in the
                orientation used to translate it.
            Otherwise, returns string derived from captured exception.
    """

    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def check_str_alphabet(string, alphabet):
    """ Determine if a string only has characters from specified
        alphabet.

        Args:
            string (str): input string.
            alphabet (str): string that has all characters from desired
            alphabet.

        Returns:
            "True" if sequence only has characters from specified
            alphabet and string "ambiguous or invalid characters" if
            it any of its characters is not in the alphabet.
    """

    valid_chars = alphabet
    if all(n in valid_chars for n in string) is True:
        return True
    else:
        return 'ambiguous or invalid characters'


def check_str_multiple(string, number):
    """ Determine if length of input string is multiple of
        a specified number.

        Args:
            string (str): input string.
            number (int): integer that will be used to check if sequence
            length is multiple of.

        Returns:
            "True" if the length of the sequence is a multiple of the
            specified number and "sequence length is not a multiple of number"
            if condition is not satisfied.
    """

    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def translate_dna(dna_sequence, table_id, min_len):
    """ Checks if sequence is valid and attempts to translate it,
        calling several functions to ensure that the sequence only has
        'ACTG', is multiple of 3 and that it can be translated in any of 4
        different orientations. Stores exceptions so that it is possible to
        understand the sequence could not be translated.

        Args:
            dna_sequence (str):
            table_id (int):

        Returns:
            If the sequence can be translated,
            a list with following elements:
                sequence (list): a list with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
                coding_strand (str): the strand orientation that had could be
                translated.
            Otherwise:
                exception_str (str): a string containing the exceptions that
                determined that the sequence could not be translated.
    """

    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the string is DNA, without ambiguous bases
    valid_dna = check_str_alphabet(original_seq, 'ACTG')
    if valid_dna is not True:
        return valid_dna

    # check if sequence is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return valid_length

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        sequence, exception_collector = retranslate(original_seq,
                                                    translating_methods[i],
                                                    table_id, strands[i],
                                                    exception_collector)

        i += 1
        if i == len(strands) or isinstance(sequence, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(sequence, list):
        return [sequence, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str


def retranslate(sequence, method, table_id, strands, exception_collector):
    """ Sends sequence for translation and collects exceptions when
        the sequence cannot be translated.

        Args:
            sequence (str): string representing a DNA sequence.
            method (str): a string specifying the sequence orientation
            that should be used to attempt translation.
            table_id (int): translation table identifier.
            strands (list): list with 4 different orientations that can
            be checked.
            exception_collector (list): list used to store all exceptions
            arising from translation attempts.

        Returns:
            A list with following elements, if the sequence can be translated:
                translated_seq (list): a list with the protein sequence and
                with the DNA sequence in the orientation used for translation.
                exception_collector (list): a list with the exceptions that are
                captured when the sequence could not be translated.
            Otherwise:
                translated_seq (str): a string with the exception/reason why
                the sequence could not be translated.
                exception_collector (list): list with all exception that have
                been captured during translation attempts of the current
                sequence.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def sequence_kmerizer(sequence, k_value):
    """
    """

    kmers = [sequence[i:i+k_value] for i in range(0, len(sequence)-k_value+1)]

    return kmers


def cluster_sequences(sorted_sequences, word_size,
                      clustering_sim, mode):
    """
    """

    clusters = {}
    reps_groups = {}
    for prot in sorted_sequences:
        protid = prot[0]
        protein = prot[1]
        kmers = sequence_kmerizer(protein, word_size)
        if len(clusters) == 0:
            clusters[protid] = [(protid, 1.0)]
            for k in kmers:
                # initiallize as set to avoid duplication of
                # ids per kmer that will lead to overestimation
                # of similarity
                reps_groups[k] = set([protid])
        else:
            current_reps = []
            for k in kmers:
                if k in reps_groups:
                    current_reps.extend(list(reps_groups[k]))

            # count number of kmer hits per representative
            # 0.20
            counts = Counter(current_reps)
            current_reps = [(k, v/len(kmers)) for k, v in counts.items() if v/len(kmers) >= clustering_sim]

            sims = sorted(current_reps, key=lambda x: x[1], reverse=True)

            if len(sims) > 0:

                if mode == 'greedy':
                    clusters[sims[0][0]].append((protid, sims[0][1]))
                elif mode == 'full':
                    for s in sims:
                        clusters[s[0]].append((protid, s[1]))
            else:
                for k in kmers:
                    if k in reps_groups:
                        reps_groups[k].add(protid)
                    else:
                        reps_groups[k] = set([protid])
                clusters[protid] = [(protid, 1.0)]

    return clusters


def cluster_prunner(clusters, sim_cutoff):
    """
    """

    excluded = []
    prunned_clusters = {}
    for rep, seqids in clusters.items():
        # this removes the representative from the cluster
        prunned_clusters[rep] = [seqid for seqid in seqids if seqid[1] > sim_cutoff]
        excluded.extend([seqid for seqid in seqids if seqid[1] <= sim_cutoff and seqid[0] != rep])

    return [prunned_clusters, excluded]


def get_sequences_by_id(sequences_index, seqids, out_file):
    """
    """

    selected = []
    seq_limit = 5000
    for i, seqid in enumerate(seqids):
        identifier = sequences_index[seqid].id
        header = '>{0}'.format(identifier)
        sequence = str(sequences_index[seqid].seq)
        selected.append(header)
        selected.append(sequence)

        if len(selected) // 2 == seq_limit or i+1 == len(seqids):

            lines = join_list(selected, '\n')
            write_to_file(lines, out_file, 'a', '\n')

            selected = []


def write_to_file(text, output_file, write_mode, end_char):
    """
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def join_list(lst, link):
    """
    """

    return link.join(lst)


def blast_inputs(clusters, output_directory):
    """
    """

    ids_to_blast = []
    for i in clusters:

        cluster_file = os.path.join(output_directory,
                                    '{0}_ids.txt'.format(i))
        cluster_ids = [i] + [seq[0] for seq in clusters[i]]
        cluster_lines = join_list(cluster_ids, '\n')
        write_to_file(cluster_lines, cluster_file, 'w', '')
        ids_to_blast.append(i)

    return ids_to_blast


def split_blast_inputs_by_core(blast_inputs, threads, blast_files_dir):
    """
    """

    splitted_ids = [[] for cpu in range(threads)]
    splitted_values = [[] for cpu in range(threads)]
    cluster_sums = [0] * threads
    i = 0
    for cluster in blast_inputs:
        cluster_file = os.path.join(blast_files_dir,
                                    '{0}_ids.txt'.format(cluster))
        with open(cluster_file, 'r') as infile:
            cluster_seqs = [line.strip() for line in infile.readlines()]
        splitted_values[i].append(len(cluster_seqs))
        splitted_ids[i].append(cluster)
        cluster_sums[i] += len(cluster_seqs)
        i = cluster_sums.index(min(cluster_sums))

    return [s for s in splitted_ids if len(s) > 0]


def cluster_blaster(inputs):
    """
    """

    blast_inputs = inputs[0:-4]
    blastp_path = inputs[-4]
    output_directory = inputs[-2]
    proteins_file = inputs[-1]
    blast_db = inputs[-3]

    indexed_fasta = SeqIO.index(proteins_file, 'fasta')

    for cluster in blast_inputs:

        cluster_id = cluster
        ids_file = os.path.join(output_directory,
                                '{0}_ids.txt'.format(cluster_id))

        with open(ids_file, 'r') as clstr:
            cluster_ids = [l.strip() for l in clstr.readlines()]

        fasta_file = os.path.join(output_directory,
                                  '{0}_protein.fasta'.format(cluster_id))
        # create file with protein sequences
        get_sequences_by_id(indexed_fasta, cluster_ids, fasta_file)

        blast_output = os.path.join(output_directory,
                                    '{0}_blast_out.tsv'.format(cluster_id))
        blast_command = ('{0} -db {1} -query {2} -out {3} '
                         '-outfmt "6 qseqid sseqid score" '
                         '-max_hsps 1 -num_threads {4} -evalue '
                         '0.001 -seqidlist {5}'.format(blastp_path,
                                                       blast_db,
                                                       fasta_file,
                                                       blast_output,
                                                       1, ids_file))

        os.system(blast_command)


# schema_short_dir = ''
# # get schema gene list
# schema_genes = [os.path.join(schema_short_dir, file) for file in os.listdir(schema_short_dir) if '.fasta' in file]

# # import all sequences
# protein_file = 'all_proteins.fasta'
# with open(protein_file, 'a') as pf:
#     proteins = []
#     for locus in schema_genes:
#         for record in SeqIO.parse(locus, 'fasta'):
#             dna_seq = str(record.seq)
#             prot = translate_dna(dna_seq, 11, 0)
#             prot = str(prot[0][0])
#             protid = record.id
#             proteins.append('>{0}\n{1}'.format(protid, prot))

#     proteins_text = '\n'.join(proteins)
#     pf.write(proteins_text)


# prots = {}
# for record in SeqIO.parse(protein_file, 'fasta'):
#     seqid = record.id
#     prot_seq = str(record.seq)
#     prots[seqid] = prot_seq

# # sort proteins by length and alphabetically
# sorted_prots = sorted(list(prots.items()),
#                       key=lambda x: (-len(x[1]), x[0]))

# clusters = cluster_sequences(sorted_prots, 4, 0.20, 'greedy')

# # remove singletons
# non_singletons = {k: v for k, v in clusters.items() if len(v) > 1}

# # remove clustered sequences that do not reach a certain threshold
# prunned_clusters = cluster_prunner(non_singletons, 0.70)[0]

# # remove clusters that only have representative
# final_clusters = {k: v for k, v in prunned_clusters.items() if len(v) > 1}

# # we have to BLAST the sequences in the remaining clusters

# # write all proteins to a file
# indexed_prots = SeqIO.index(protein_file, 'fasta')

# seqids = [[e[0] for e in v] for v in list(final_clusters.values())]
# seqids = list(itertools.chain.from_iterable(seqids))

# get_sequences_by_id(indexed_prots, seqids, 'test_paralog_out')

# makedb_cmd = 'makeblastdb -in {0} -parse_seqids -dbtype prot >/dev/null 2>&1'.format('test_paralog_out')
# os.system(makedb_cmd)

# # BLAST necessary sequences against only the sequences from the same cluster
# blast_db = os.path.join('test_paralog_out')

# blast_results_dir = os.path.join('blast_results')
# os.mkdir(blast_results_dir)

# seqids_to_blast = blast_inputs(final_clusters, blast_results_dir)

# # distribute clusters per available cores, try to group inputs into
# # even groups in terms of number of clusters and sum of number of sequences
# # per inputs group
# splitted_seqids = split_blast_inputs_by_core(seqids_to_blast,
#                                              6,
#                                              blast_results_dir)

# for s in range(len(splitted_seqids)):
#     splitted_seqids[s].append('blastp')
#     splitted_seqids[s].append(blast_db)
#     splitted_seqids[s].append(blast_results_dir)
#     splitted_seqids[s].append('test_paralog_out')


# p = Pool(processes = 6)
# r = p.map_async(cluster_blaster, splitted_seqids)
# r.wait()
