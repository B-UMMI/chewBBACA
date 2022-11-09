#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with biological sequence
manipulation.

Code documentation
------------------
"""


from Bio.Seq import Seq
from collections import Counter

try:
    from utils import (constants as ct,
                       file_operations as fo,
                       iterables_manipulation as im,
                       fasta_operations as fao)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 iterables_manipulation as im,
                                 fasta_operations as fao)


def translate_sequence(dna_str, table_id):
    """Translate a DNA sequence using the BioPython package.

    Parameters
    ----------
    dna_str : str
        String representing a DNA sequence.
    table_id : int
        Translation table identifier.

    Returns
    -------
    protseq : Bio.Seq.Seq
        Protein sequence created by translating the
        input DNA sequence.
    """
    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """Attempt to translate a DNA sequence in specified orientation.

    Attempts to translate an input DNA sequence in specified
    orientation and stores exceptions when the input sequence
    cannot be translated.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.
    method : str
        Sequence orientation to attempt translation.
    table_id : int
        Translation table identifier.

    Returns
    -------
    If the sequence can be translated:
        protseq : Bio.Seq.Seq
            Translated DNA sequence.
        myseq : str
            String representing the DNA sequence in the
            orientation used to translate it.
    Otherwise, returns string with the description of the
    exception that was raised.
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
            myseq = im.reverse_complement(myseq, ct.DNA_BASES)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = im.reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = im.reverse_str(myseq)
            myseq = im.reverse_complement(myseq, ct.DNA_BASES)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def translate_dna(dna_sequence, table_id, min_len):
    """Determine if a DNA sequence is valid and translate it.

    Checks if sequence is valid and attempts to translate
    it, calling several functions to ensure that the sequence
    only has 'ACTG', is multiple of 3 and that it can be
    translated in any of 4 different orientations. Stores
    exceptions so that it is possible to understand why the
    sequence could not be translated.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.
    table_id : int
        Translation table identifier.
    min_len : int
        Minimum sequence length. Sequences shorter
        than this value are not translated.

    Returns
    -------
    If the sequence can be translated:
        sequence : list
            List with two elemets, the protein sequence
            and the DNA sequence in the correct orientation.
        coding_strand : str
            The sequence orientation that codes for the
            protein.
    Otherwise:
        exception_str : str
            A string containing the exceptions that
            explain why the the sequence could not be
            translated.
    """
    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the sequence has ambiguous bases
    valid_dna = im.check_str_alphabet(original_seq, ct.DNA_BASES)
    if valid_dna is not True:
        return 'ambiguous or invalid characters'

    # check if sequence size is multiple of three
    valid_length = im.check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return 'sequence length is not a multiple of 3'

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        translated_seq = translate_dna_aux(original_seq, translating_methods[i], table_id)
        if not isinstance(translated_seq, list):
            exception_collector.append('{0}({1})'.format(strands[i],
                                                         translated_seq.args[0]))

        i += 1
        if i == len(strands) or isinstance(translated_seq, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(translated_seq, list):
        return [translated_seq, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str


def determine_duplicated_seqs(sequences):
    """Create mapping between sequences and sequence identifiers.

    Parameters
    ----------
    sequences : dict
        Dictionary with sequence identifiers as keys and
        sequences as values.

    Returns
    -------
    equal_seqs : dict
        Dictionary with sequences as keys and sequence
        identifiers that are associated with each
        sequence as values.
    """
    equal_seqs = {}
    for seqid, seq in sequences.items():
        # if protein sequence was already added as key
        if seq in equal_seqs:
            # append new protid
            equal_seqs[seq].append(seqid)
        # else add new protein sequence as key and protid
        # as value
        else:
            equal_seqs[seq] = [seqid]

    return equal_seqs


def determine_longest(seqids, sequences):
    """Find the longest sequence in a set of sequences.

    Parameters
    ----------
    seqids : list
        List with sequence identifiers.
    sequences : dict
        Dictionary with sequence identifiers as keys
        and sequences as values.

    Returns
    -------
    chosen : str
        Sequence identifier of the longest sequence.
    """
    seqids_tups = [(seqid, sequences[seqid]) for seqid in seqids]
    sorted_tups = sorted(seqids_tups, key=lambda x: len(x[1]), reverse=True)
    chosen = sorted_tups[0][0]

    return chosen


def determine_mode(values):
    """Compute the mode based on a list of integers.

    Parameters
    ----------
    values : list
        List with integer values.

    Returns
    -------
    modes : list
        The most frequent integer values.
    """
    # determine frequency of each length value
    counts = Counter(values)

    # order by most common first
    most_common = counts.most_common()

    # get most common
    modes = [most_common[0][0]]

    # determine if there are more length values that are just as common
    modes += [m[0] for m in most_common[1:] if m[1] == most_common[0][1]]

    return modes


def mode_filter(sequences, size_threshold):
    """Find sequences with length that deviates from the mode value.

    Determines the mode from a set of input sequences and
    identifies sequences that have a length value smaller
    or greater than the mode based on a threshold.

    Parameters
    ----------
    sequences : dict
        Dictionary with sequence identifiers as keys and
        sequences as values.
    size_threshold : float
        Sequences with +/- this value * mode will be
        reported as above or below the mode.

    Returns
    -------
    A list with the following variables:
        modes : list
            List with mode values determined based on the
            length of input sequences.
        alm : list
            List with the sequence identifiers of the
            sequences that are above the mode value by
            mode*size_threshold.
        asm : list
            List with the sequence identifiers of the
            sequences that are below the mode value by
            mode*size_threshold.
        seqs_lengths : dict
            Dictionary with sequence identifiers as keys
            and sequence lengths as values.
    """
    # determine length value of all sequences
    seqs_lengths = {seqid: len(seq) for seqid, seq in sequences.items()}

    # determine mode/s
    modes = determine_mode(list(seqs_lengths.values()))

    # determine top and bot length value limits
    max_mode = max(modes)
    top_limit = max_mode + (max_mode*size_threshold)
    min_mode = min(modes)
    bot_limit = min_mode - (min_mode*size_threshold)

    # determine sequences that are below or above limits
    alm = [seqid for seqid, length in seqs_lengths.items()
           if length > top_limit]
    asm = [seqid for seqid, length in seqs_lengths.items()
           if length < bot_limit]

    return [modes, alm, asm, seqs_lengths]


def get_seqs_dicts(fasta_path, gene_id, table_id, min_len, size_threshold):
    """Translate DNA sequences in a FASTA file.

    Translates the set of alleles from a gene. Identifies
    sequences that cannot be translated according to the
    criteria enforced by chewBBACA.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file with DNA sequences.
    gene_id : str
        Gene identifier.
    table_id : int
        Translation table identifier.
    min_len : int
        Minimum sequence length. Sequences shorter
        than this value are not translated.
    size_threshold : float
        Sequences with +/- this value * mode will be
        reported as above or below the mode.

    Returns
    -------
    List with following elements:
        dna_seqs : dict
            Dictionary with sequence identifiers as keys
            and DNA sequences as values.
        prot_seqs : dict
            Dictionary with protein identifiers as keys
            and Protein sequences as values. Keys are
            consecutive integers to enable alignment
            with BLASTp without getting exceptions due
            to long sequence identifiers.
        invalid : list
            List with sequence identifiers of alleles
            that are not valid because they could not be
            translated.
        seqids_map : dict
            Dictionary with consecutive integers as keys
            and original allele identifiers as values.
        total_seqs : int
            Total number of sequences that was processed
            (including invalid alleles).
    """
    sequences = fao.import_sequences(fasta_path)

    # translate sequences
    translated_seqs = {k: translate_dna(v, table_id, min_len)
                       for k, v in sequences.items()}

    # add locus identifier to headers
    # some headers might only have the allele identifier
    seqids = list(translated_seqs.keys())
    new_seqids = {}
    for i in seqids:
        new_id = '{0}_{1}'.format(gene_id, i.split('_')[-1])
        new_seqids[i] = new_id

    # switch ids
    sequences = {new_seqids[k]: v
                 for k, v in translated_seqs.items()}

    valid = {k: v
             for k, v in sequences.items()
             if isinstance(v, list) is True}
    invalid = [[k, v]
               for k, v in sequences.items()
               if isinstance(v, list) is False]

    seqid = 1
    seqids_map = {}
    dna_seqs = {}
    prot_seqs = {}
    for k, v in valid.items():
        seqids_map[str(seqid)] = k
        dna_seqs[k] = v[0][1]
        prot_seqs[str(seqid)] = str(v[0][0])
        seqid += 1

    if size_threshold is not None and len(prot_seqs) > 0:
        # remove alleles based on length mode and size threshold
        modes, alm, asm, alleles_lengths = mode_filter(dna_seqs, size_threshold)
        excluded = set(asm + alm)

        dna_seqs = {seqid: seq
                    for seqid, seq in dna_seqs.items()
                    if seqid not in excluded}
        prot_seqs = {seqid: seq
                     for seqid, seq in prot_seqs.items()
                     if seqids_map[seqid] not in excluded}

        modes_concat = ':'.join(map(str, modes))
        st_percentage = int(size_threshold*100)
        invalid += [[s, 'allele greater than {0}% locus length mode '
                     '({1}>{2})'.format(st_percentage,
                                        alleles_lengths[s],
                                        modes_concat)]
                    for s in alm]
        invalid += [[s, 'allele smaller than {0}% locus length mode '
                     '({1}<{2})'.format(st_percentage,
                                        alleles_lengths[s],
                                        modes_concat)]
                    for s in asm]

    total_seqs = len(translated_seqs)

    return [dna_seqs, prot_seqs, invalid, seqids_map, total_seqs]


def translate_coding_sequences(seqids, protein_file, sequences_file,
                               translation_table, minimum_length):
    """Translate coding sequences.

    Parameters
    ----------
    seqids : list
        List with the sequence identifiers of the sequences
        to be translated.
    protein_file : str
        Path to a file to save protein sequences.
    sequences_file : str
        Path to the FASTA file that contains the DNA sequences.
    translation_table : int
        Translation table identifier.
    minimum_length : int
        The minimum sequence length value.

    Returns
    -------
    A list with following elements:
        invalid_alleles : list
            List with one sublist per invalid allele.
            Each sublist contains a sequence identifer
            and the exception message returned after
            attempting translation.
        total_seqs : int
            Total number of DNA sequences that were
            translated.
    """
    # define limit of records to keep in memory
    total_seqs = 0
    prot_lines = []
    line_limit = 10000
    invalid_alleles = []
    cds_index = fao.index_fasta(sequences_file)
    for i, seqid in enumerate(seqids):
        sequence = str(cds_index.get(seqid).seq)

        translation = translate_dna(sequence, translation_table, minimum_length)
        if isinstance(translation, list):
            prot_lines.append('>{0}'.format(seqid))
            prot_lines.append(str(translation[0][0]))
            total_seqs += 1
        # if returned value is a string, translation failed and
        # string contains exceptions
        elif isinstance(translation, str):
            invalid_alleles.append([seqid, translation])

        if len(prot_lines)//2 == line_limit or i+1 == len(seqids):
            prot_lines = im.join_list(prot_lines, '\n')
            fo.write_to_file(prot_lines, protein_file, 'a', '\n')
            prot_lines = []

    return [invalid_alleles, total_seqs]


def determine_distinct(sequences_file, unique_fasta, map_ids):
    """Identify duplicated sequences in a FASTA file.

    Parameters
    ----------
    sequences_file : str
        Path to a FASTA file.
    unique_fasta : str
        Path to a FASTA file that will be created to
        store distinct sequences.
    map_ids : dict
        Dictionary with mapping between genome string
        identifiers and genome integer identifiers.

    Returns
    -------
    pickle_out : str
        Pickled file that contains a dictionary with sequence
        hashes as keys and a list with pairs of protein
        identifiers and genome identifiers as values
    """
    out_seqs = []
    duplicates = {}
    exausted = False
    # limit of 10000 Fasta records in memory
    out_limit = 10000
    seq_generator = fao.sequence_generator(sequences_file)
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            # seq object has to be converted to string
            seqid = record.id
            sequence = str(record.seq.upper())

            # use digest() instead of hexdigest() to reduce memory usage?
            seq_hash = im.hash_sequence(sequence)

            # add unseen sequence to Fasta file with distinct sequences
            if seq_hash not in duplicates:
                recout = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [seqid, sequence])
                out_seqs.append(recout)

            # add CDS hash as key
            # add genome integer identifier and protein identifier to values list
            # genome identifier and protein identifier can be used to fetch sequences
            genome_id, protid = seqid.split('-protein')
            genome_id = map_ids[genome_id]
            duplicates.setdefault(seq_hash, []).extend([int(protid), int(genome_id)])
        else:
            exausted = True

        # write Fasta records to file
        if len(out_seqs) == out_limit or exausted is True:
            if len(out_seqs) > 0:
                out_seqs = im.join_list(out_seqs, '\n')
                fo.write_to_file(out_seqs, unique_fasta, 'a', '\n')
                # reset list to avoid writing same records multiple times
                out_seqs = []

    # save dictionary with genome integer identifiers per distinct sequence
    # to pickle and only return file path to avoid keeping all dicts from
    # parallel processes in memory
    pickle_out = unique_fasta.replace('.fasta', '.duplicates')
    fo.pickle_dumper(duplicates, pickle_out)

    return pickle_out


def determine_small(sequences_file, minimum_length, variation=0):
    """Find protein sequences that are shorter than a specified length.

    Parameters
    ----------
    sequences_file : str
        Path to a FASTA file.
    minimum_length : int
        Sequences with a length value below this value
        are considered small.
    variation : float
        Accept sequences with length variation of up to
        minus (`minimum_length`*`variation`).

    Returns
    -------
    small_seqids : list
        List with the identifiers of small sequences.
    """
    variation = minimum_length - (minimum_length*variation)
    small_seqids = []
    seq_generator = fao.sequence_generator(sequences_file)
    for record in seq_generator:
        # seq object has to be converted to string
        sequence = str(record.seq)
        seqid = record.id

        if len(sequence) < variation:
            small_seqids.append(seqid)

    return small_seqids


def apply_bsr(blast_results, fasta_file, bsr, ids_dict):
    """Find similar sequences based on the BLAST Score Ratio.

    Parameters
    ----------
    blast_results : list
        List with the path to a file with BLAST
        results in tabular format.
    fasta_file : str
        Path to a FASTA file that contains the
        sequences that were aligned.
    bsr : float
        The BSR value to use as threshold
    ids_dict : dict
        Dictionary with the mapping between
        sequence identifiers used for BLAST and
        the original sequence identifiers.

    Returns
    -------
    excluded_alleles : list
        List with the identifiers of the sequences
        that were highly similar to other sequences.
    """
    self_scores = {r[0]: r[-1] for r in blast_results if r[0] == r[4]}
    # do not include self-scores lines, no need to evaluate those hits
    blast_results = [r for r in blast_results if r[0] != r[4]]

    lengths = {}
    for k in self_scores:
        record = fasta_file.get(ids_dict[k])
        sequence = str(record.seq)
        lengths[k] = len(sequence)

    excluded_alleles = []
    for res in blast_results:

        query = res[0]
        hit = res[4]
        score = res[-1]

        if query not in excluded_alleles:
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

                    elif hit_length <= query_length:
                        excluded_alleles.append(hit)
            # it might not work because there is no self score for
            # some sequences due to low complexity regions...
            except Exception:
                excluded_alleles.append(query)

    return excluded_alleles
