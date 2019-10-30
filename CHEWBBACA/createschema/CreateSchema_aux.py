#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION





"""


import os
import csv
import pickle
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

import new_utils
import init_schema_4_bbaca


def join_list(lst, link):
    """
    """

    return link.join(lst)


def genome_id(file_path):
    """ Get file basename without extension suffix.

        Args:
            file_path (str): the path to the file.
            Full path or relative path.

        Returns:
            file (str): the file basename without the
            file extension suffix.
    """

    file = os.path.basename(file_path)

    file = file.split('.')[0]

    return file


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Args:
            fasta_path (str): full path to the FASTA file.

        Returns:
            dictionary that has sequences ids as keys and DNA
            sequences as values.
    """

    seqs_dict = {}
    # use BioPython to read FASTA file and get each contig sequence
    for record in SeqIO.parse(fasta_path, 'fasta', generic_dna):
        # seq object has to be converted to string
        sequence = str(record.seq.upper())
        seqid = record.id

        # add contig id as key and DNA sequence as value
        seqs_dict[seqid] = sequence

    return seqs_dict


def extract_coding_sequences(reading_frames, contigs, starting_id):
    """ Extracts CDSs from contigs based on the start
        and stop codon positions determined by Prodigal.

        Args:
            reading_frames (str): full path to the ORF file derived
            from Prodigal.
            contigs (dict): a dictionary with contigs ids as keys
            and contigs sequences as values.
            starting_id (int): integer identifier to give to the
            first CDS extracted and that will be incremented to
            serve as identifier for subsequent CDSs.

        Returns:
            coding_sequences (dict): dictionary with CDSs ids as keys and
            CDSs DNA sequences as values.
            coding_sequences_info (list): list with a sublist for each
            extracted CDS. Sublists have information about the extracted CDS
            (identifier of the contig where the CDS was found, start position
            in the contig, stop position in the contig, sequence identifier
            attributed to that CDS and the strand that coded for that CDS.)
    """

    # load binary file with a list of lists for each contig
    # each sublist has a start codon and stop codon positions
    # in the contig
    with open(reading_frames, 'rb') as orf_file:
        rfs = pickle.load(orf_file)

    seqid = starting_id
    coding_sequences = {}
    coding_sequences_info = []
    for contig_id, frames in rfs.items():
        # for each start and stop codon in the contig
        for coding_sequence in frames:
            start_codon = coding_sequence[0]
            stop_codon = coding_sequence[1]
            strand = coding_sequence[2]
            # extract CDS sequence
            cds_sequence = contigs[contig_id][start_codon:stop_codon].upper()
            # check coding strand to change sequence orientation
            # if needed
            if strand == 0:
                cds_sequence = reverse_complement(cds_sequence)

            # store CDS with unique id
            coding_sequences[seqid] = cds_sequence

            # store CDS information
            coding_sequences_info.append([contig_id, str(start_codon),
                                          str(stop_codon), str(seqid),
                                          str(strand)])

            # increment seqid
            seqid += 1

    return [coding_sequences, coding_sequences_info]


def reverse_str(string):
    """ Reverse character order in input string.

        Args:
            string (str): string to be reversed.

        Returns:
            revstr (str): reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Args:
            dna_sequence (str): string representing a DNA sequence.

        Returns:
            reverse_complement_strand (str): the reverse complement
            of the DNA sequence, without lowercase letters.
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    complement_bases = []
    for base in bases:
        if base in base_complement:
            complement_bases.append(base_complement[base])
        else:
            complement_bases.append(base.upper())

    complement_strand = join_list(complement_bases, '')

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


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
    # sequences must be a complete and valid CDS
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


def translate_dna(dna_sequence, table_id):
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
        exception_str = join_list(exception_collector, ',')
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


def write_to_file(text, output_file, write_mode, end_char):
    """
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def translate_coding_sequences(sequences_file, valid_seqs, dna_valid_file,
                               protein_valid_file, table_id):
    """ Translates CDSs into protein sequences.

        Args:
            cdss (dict): a dictionary with CDSs ids as keys and CDSs DNA
            sequences as values.

        Returns:
            prots (dict): a dictionary with CDSs/proteins ids as keys and
            protein sequences as values.
            trans_state (dict): a dictionary with the CDSs/proteins ids as
            keys and the DNA strand that coded for those proteins.
            ambiguous (dict): a dictionary with CDSs/proteins ids as keys
            and CDSs DNA sequences that had ambiguous bases as values.
            untranslatable (dict): a dictionary with CDSs/proteins ids as
            keys and CDSs DNA sequences that could not be translated.
    """

    # define limit of records to keep in memory
    dna_lines = []
    total_seqs = 0
    prot_lines = []
    line_limit = 5000
    invalid_alleles = []
    dna_seqs_file = dna_valid_file
    protein_seqs_file = protein_valid_file

    cds_index = SeqIO.index(sequences_file, 'fasta')

    for i, seqid in enumerate(valid_seqs):
        sequence = str(cds_index.get(seqid).seq)

        translation = translate_dna(sequence, table_id)

        if isinstance(translation, list):
            dna_lines.append('>{0}'.format(seqid))
            dna_lines.append(translation[0][1])
            prot_lines.append('>{0}'.format(seqid))
            prot_lines.append(str(translation[0][0]))
            total_seqs += 1
        # if returned value is a string, translation failed and
        # string contains exceptions
        elif isinstance(translation, str):
            invalid_alleles.append([seqid, translation])

        if len(dna_lines)//2 == line_limit or i+1 == len(valid_seqs):

            dna_lines = join_list(dna_lines, '\n')
            write_to_file(dna_lines, dna_seqs_file, 'a', '\n')
            dna_lines = []

            prot_lines = join_list(prot_lines, '\n')
            write_to_file(prot_lines, protein_seqs_file, 'a', '\n')
            prot_lines = []

    return [invalid_alleles, total_seqs]


def determine_repeated(sequences_file, repeated_output, unique_fasta):
    """
    """

    total = 0
    seqs_dict = {}
    for record in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        sequence = str(record.seq.upper())
        seqid = record.id

        # store each sequence and first seqid found
        # with that sequence
        if sequence not in seqs_dict:
            seqs_dict[sequence] = seqid
        elif sequence in seqs_dict:
            total += 1

    # get sequences with more than one seqid
    out_limit = 5000
    out_seqs = []
    for seq, seqid in seqs_dict.items():
        header = '>{0}'.format(seqid)
        sequence = seq
        out_seqs.append(header)
        out_seqs.append(sequence)
        if len(out_seqs)/2 == out_limit:
            out_seqs = join_list(out_seqs, '\n')
            write_to_file(out_seqs, unique_fasta, 'a', '\n')
            out_seqs = []

    out_seqs = join_list(out_seqs, '\n')
    write_to_file(out_seqs, unique_fasta, 'a', '\n')
    unique_seqids = list(seqs_dict.values())

    return [total, unique_seqids]


def determine_small(sequences_file, minimum_length):
    """ Find protein sequences that are shorter than desired length.

        Args:
            prots (dict): a dictionary with protein ids as keys and protein
            sequences as values.
            min_len (int): Proteins with a number of amino acids lower than
            this value are considered small.

        Returns:
            small_proteins (dict): a dictionary with the ids of small proteins
            as keys and their amino acid sequence as values.
    """

    small_seqs = []
    valid_seqs = []
    for record in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        #sequence = sequences_index.get(seqid).seq
        sequence = str(record.seq)
        seqid = record.id

        if len(sequence) < minimum_length:
            small_seqs.append(seqid)
        else:
            valid_seqs.append(seqid)

    return [small_seqs, valid_seqs]


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


def create_fasta_lines(sequences, genome_id):
    """

        Args:

        Returns:

        Example:
    """

    lines = []
    for seqid in sequences:
        header = '>' + genome_id + '-protein' + str(seqid)
        sequence = sequences[seqid]

        lines.append(header)
        lines.append(sequence)

    return lines


def write_fasta(fasta_lines, output_file):
    """
    """

    joined_lines = join_list(fasta_lines, '\n')

    write_to_file(joined_lines, output_file, 'a', '\n')


#file_name = protein_table 
#cds_info = genome_info[1]

def write_protein_table(file_name, genome_id, cds_info):
    """
    """

    table_lines = [[genome_id] + protein_info
                   for protein_info in cds_info]
    table_lines = [join_list(line, '\t') for line in table_lines]
    table_text = join_list(table_lines, '\n')
    write_to_file(table_text, file_name, 'a', '\n')


def read_blast_tabular(blast_tabular_file):
    """ Read a file with BLAST results in tabular format

        Args:
            blast_tabular_file (str): path to output file of BLAST.

        Returns:
            blasting_results (list): a list with a sublist per line
            in the input file.
    """

    with open(blast_tabular_file, 'r') as blastout:
        blasting_results = []
        reader = csv.reader(blastout, delimiter='\t')
        for row in reader:
            blasting_results.append(row)

    return blasting_results


def get_schema_reps(dna_index, schema_seqids, protein_file, output_file):
    """
    """

    schema_lines = []
    for seqid in schema_seqids:
        header = '>{0}'.format(seqid)
        sequence = str(dna_index[seqid].seq)

        schema_lines.append(header)
        schema_lines.append(sequence)

    schema_text = join_list(schema_lines, '\n')
    write_to_file(schema_text, output_file, 'w', '')


def build_schema(schema_file, output_path):
    """
    """

    total_genes = 0
    schema_files = []
    for record in SeqIO.parse(schema_file, 'fasta', IUPAC.unambiguous_dna):
        file_name = record.name
        file_name = new_utils.replace_multiple_characters(file_name)

        new_file = '{0}{1}'.format(file_name, '.fasta')
        new_file_path = os.path.join(output_path, new_file)
        schema_files.append(new_file_path)

        header = '>{0}_1'.format(file_name)
        sequence = str(record.seq).upper()
        file_text = join_list([header, sequence], '\n')
        write_to_file(file_text, new_file_path, 'w', '\n')

        total_genes += 1

    init_schema_4_bbaca.get_Short(schema_files)
    print('\nTotal of {0} loci that constitute the schema.'.format(total_genes))


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

    return splitted_ids


def divide_list_into_n_chunks(list_to_divide, n):
    """
    """

    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    return sublists


def apply_bsr(inputs):
    """
    """

    blast_file = inputs[0]
    blast_results = read_blast_tabular(blast_file)
    self_scores = {r[0]: r[2] for r in blast_results if r[0] == r[1]}
    # do not include self-scores lines, no need to evaluate those hits
    blast_results = [r for r in blast_results if r[0] != r[1]]

    fasta_file = inputs[1]
    lengths = {}
    for k in self_scores:
        record = fasta_file.get(k)
        sequence = str(record.seq)
        lengths[k] = len(sequence)

    bsr = inputs[2]

    excluded_alleles = []
    for res in blast_results:

        query = res[0]
        hit = res[1]
        score = res[2]

        if query not in excluded_alleles:
            #print(res)
            # try to apply BSR strategy
            try:
                self_blast_score = self_scores[query]

                query_length = lengths[query]
                hit_length = lengths[hit]
                blast_score_ratio = float(score) / float(self_blast_score)

                # BSR has to be greater than threshold, just as in the original function
                if blast_score_ratio >= bsr and hit not in excluded_alleles:

                    if hit_length > query_length and query not in excluded_alleles:
                        excluded_alleles.append(query)
                        #print(query)

                    elif hit_length <= query_length:
                        excluded_alleles.append(hit)
                        #print(hit)
            # it might not work because there is no self score for
            # some sequences due to low complexity regions...
            except Exception:
                excluded_alleles.append(query)

    return excluded_alleles


def sequence_kmerizer(sequence, k_value):
    """
    """

    kmers = [sequence[i:i+k_value] for i in range(0, len(sequence)-k_value+1)]

    return kmers


def cluster_sequences(sorted_sequences, word_filter, filtering_sim,
                      word_size, clustering_sim, mode):
    """
    """

    clusters = {}
    reps_groups = {}
    for prot in sorted_sequences:
        protid = prot[0]
        protein = prot[1]
        if len(clusters) == 0:
            clusters[protid] = [(protid, 1.0)]
            kmers = sequence_kmerizer(protein, word_size)
            longmers = sequence_kmerizer(protein, word_filter)
            for k in longmers:
                # initiallize as set to avoid duplication of
                # ids per kmer that will lead to overestimation
                # of similarity
                reps_groups[k] = set([protid])
        else:
            kmers = sequence_kmerizer(protein, word_size)
            longmers = sequence_kmerizer(protein, word_filter)
            current_reps = []
            for k in longmers:
                if k in reps_groups:
                    current_reps.extend(list(reps_groups[k]))

            # count number of kmer hits per representative
            # 0.20
            counts = Counter(current_reps)
            current_reps = [(k, v/len(kmers)) for k, v in counts.items() if v/len(kmers) >= clustering_sim]

#            intersections = []
#            for r in current_reps:
#                shared_kmers = len(set.intersection(set(reps[r]), set(kmers)))
#                intersections.append((r, shared_kmers/len(kmers)))

#            sims = [s for s in intersections if s[1] >= clustering_sim]
            # sort to get most similar at index 0
            sims = sorted(current_reps, key=lambda x: x[1], reverse=True)

            if len(sims) > 0:

                if mode == 'greedy':
                    clusters[sims[0][0]].append((protid, sims[0][1]))
                elif mode == 'full':
                    for s in sims:
                        clusters[s[0]].append((protid, s[1]))
            else:
                for k in longmers:
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
        prunned_clusters[rep] = [seqid for seqid in seqids if seqid[1] < sim_cutoff]
        excluded.extend([seqid for seqid in seqids if seqid[1] >= sim_cutoff and seqid[0] != rep])

    return [prunned_clusters, excluded]


def determine_singletons(clusters):
    """
    """

    singletons = {k: v for k, v in clusters.items() if len(v) == 0}

    return singletons


def remove_clusters(clusters, cluster_ids):
    """
    """

    new_clusters = {k: v for k, v in clusters.items() if k not in cluster_ids}

    return new_clusters


def intra_cluster_sim(clusters, protein_file):
    """
    """

    excluded_dict = {}
    for k, v in clusters.items():
        cluster_ids = v

        # get sequences
        cluster_sequences = {}
        for seqid in cluster_ids:
            cluster_sequences[seqid[0]] = str(protein_file[seqid[0]].seq)

        # get all kmers per sequence
        kmers_mapping = {}
        cluster_kmers = {}
        for seqid, prot in cluster_sequences.items():
            prot_kmers = sequence_kmerizer(prot, 4)
            cluster_kmers[seqid] = prot_kmers
            
            for kmer in prot_kmers:
                if kmer in kmers_mapping:
                    kmers_mapping[kmer].add(seqid)
                else:
                    kmers_mapping[kmer] = set([seqid])

        sims_cases = {}
        excluded = []
        for seqid, kmers in cluster_kmers.items():
            if seqid not in excluded:
                query_kmers = kmers
                current_reps = []
                for kmer in query_kmers:
                    if kmer in kmers_mapping:
                        current_reps.extend(list(kmers_mapping[kmer]))
                
                counts = Counter(current_reps)
                current_reps = [(s, v/len(kmers)) for s, v in counts.items() if v/len(kmers) >= 0.8]
                
                sims = sorted(current_reps, key=lambda x: x[1], reverse=True)
                
                if len(sims) > 1:
                    candidates = [s for s in sims if s[0] != seqid]
                    for c in candidates:
                    
                        if len(query_kmers) >= len(cluster_kmers[c[0]]):
                            sims_cases[c[0]] = (seqid, c[1])
                            excluded.append(c[0])
                        elif len(cluster_kmers[c[0]]) > len(query_kmers):
                            sims_cases[seqid] = (c[0], c[1])
                            excluded.append(seqid)

        excluded_dict[k] = [list(set(excluded)), sims_cases]

    return excluded_dict
