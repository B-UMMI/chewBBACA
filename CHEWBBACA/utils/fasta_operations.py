#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to work with FASTA files.

Code documentation
------------------
"""


from Bio import SeqIO
from Bio.SeqIO import FastaIO

try:
    from utils import (file_operations as fo,
                       iterables_manipulation as im,
                       sequence_manipulation as sm,
                       constants as ct)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (file_operations as fo,
                                 iterables_manipulation as im,
                                 sequence_manipulation as sm,
                                 constants as ct)


def sequence_generator(input_file):
    """Create a SeqRecord iterator.

    Parameters
    ----------
    input_file : str
        Path to a Fasta file.

    Returns
    -------
    records : Bio.SeqIO.FastaIO.FastaIterator
        SeqRecord iterator.
    """
    records = SeqIO.parse(input_file, 'fasta')

    return records


def index_fasta(fasta_file):
    """Create index to retrieve data from a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.

    Returns
    -------
    fasta_index : Bio.File._IndexedSeqFileDict
        FASTA file index.
    """
    fasta_index = SeqIO.index(fasta_file, 'fasta')

    return fasta_index


def import_sequences(input_file):
    """Import sequences from a FASTA file.

    Parameters
    ----------
    input_file : str
        Path to a FASTA file.

    Returns
    -------
    seqs_dict : dict
        Dictionary with sequence identifiers as keys and
        sequences as values.
    """
    records = sequence_generator(input_file)
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def count_sequences(fasta_file):
    """Count the number of sequences in a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.

    Returns
    -------
    total_seqs : int
        Number of sequences in the FASTA file.
    """
    with open(fasta_file, 'r') as infile:
        sequence_headers = [line for line in infile if line.startswith('>')]
        total_seqs = len(sequence_headers)

    return total_seqs


def write_records(records, output_file):
    """Write FASTA records (BioPython SeqRecord) to a file.

    Parameters
    ----------
    records : list
        List with BioPython SeqRecord objects.
    output_file : str
        Path to the output file.
    """
    with open(output_file, 'w') as output_handle:
        fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(records)


def integer_headers(input_fasta, output_fasta, start=1, limit=50000):
    """Switch sequence headers in Fasta file by integer values.

    Parameters
    ----------
    input_fasta : str
        Path to a FASTA file.
    output_fasta : str
        Path to the output file with modified headers.
    start : int
        Integer value for the first identifier.
    limit : int
        Maximum number of FASTA records to keep in
        memory.

    Returns
    -------
    ids_map : dict
        Dictionary with mapping between integer and original
        headers.
    """
    seqs = []
    ids_map = {}
    exausted = False
    seq_generator = sequence_generator(input_fasta)
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            new_id = 'seq_{0}'.format(start)
            ids_map[new_id] = record.id
            sequence = str(record.seq)
            new_rec = fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
                                       [new_id, sequence])
            seqs.append(new_rec)
            start += 1
        elif record is None:
            exausted = True

        if len(seqs) == limit or exausted is True:
            fo.write_lines(seqs, output_fasta, write_mode='a')
            seqs = []

    return ids_map


def fasta_str_record(record_template, record_data):
    """Create the string representation of a FASTA record.

    Parameters
    ----------
    record_template : str
        String template to construct the FASTA record.
    record_data : list
        List with the elements to add to the string.

    Returns
    -------
    record : str
        String representation of the FASTA record.
    """
    record = record_template.format(*record_data)

    return record


def fasta_lines(template, records_data):
    """Create a list with FASTA records.

    Parameters
    ----------
    template : str
        String template to construct the FASTA record.
    records_data : list
        A list with one sublist per FASTA record.
        Each sublist contains the elements to insert
        inside the template placeholders.

    Returns
    -------
    seqs_lines : list
        A list with strings representing FASTA records.
    """
    seqs_lines = [fasta_str_record(template, arg) for arg in records_data]

    return seqs_lines


def validate_fasta(file_path):
    """Check if a file is a FASTA file.

    Parameters
    ----------
    file_path : str
        Path to the file.

    Returns
    -------
    True if file has valid FASTA format,
    False otherwise.
    """
    records = SeqIO.parse(file_path, 'fasta')

    # returns False if it was not a FASTA file
    return any(records)


def filter_non_fasta(files):
    """Select FASTA files from a list with file paths.

    Parameters
    ----------
    files : list
        A list that contains file paths.

    Returns
    -------
    fasta_files : list
        List that contains paths to FASTA files.
    """
    fasta_files = [file for file in files if validate_fasta(file) is True]

    return fasta_files


def sequence_lengths(fasta_file, hashed=False):
    """Determine length of sequences in a FASTA file.

    Read Fasta file and create dictionary with mapping
    between sequence identifiers and sequence lengths.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.
    hashed : bool
        If False, sequence identifiers are extracted from
        sequence headers. If True, sequence hashes will be
        used as keys.

    Returns
    -------
    lengths : dict
        Dictionary with sequence identifiers as keys and
        sequence lengths as values.
    """
    records = sequence_generator(fasta_file)
    if hashed is False:
        lengths = {rec.id: len(rec.seq) for rec in records}
    else:
        lengths = {im.hash_sequence(str(rec.seq)):
                   len(rec.seq) for rec in records}

    return lengths


def get_sequences_by_id(sequences, seqids, output_file, limit=50000):
    """Retrieve sequences from indexed FASTA file or dictionary.

    Retrieves sequences based on sequence identifiers from a FASTA
    file indexed with SeqIO.index or from a dictionary with sequence
    identifiers as keys and sequences as values.

    Parameters
    ----------
    sequences : dict or Bio.File._IndexedSeqFileDict
        Dictionary with seqids as keys and sequences as values or
        a FASTA file index created with BioPython.
    seqids : list
        List with the identifiers of the sequences that should be
        retrieved.
    output_file : str
        Path to the FASTA file to which selected sequences will
        be saved.
    limit : int
        Maximum number of sequences that will be kept in memory
        at a time (to avoid keeping huge datasets in memory).

    Returns
    -------
    total_selected : int
        Total number of records written to the output file.
    """
    # using a generator
    if type(sequences) == dict:
        seqs = ((seqid, sequences[seqid]) for seqid in seqids)
    else:
        seqs = ((seqid, str(sequences[seqid].seq)) for seqid in seqids)

    records = []
    total_selected = 0
    exhausted = False
    while exhausted is False:
        record = next(seqs, None)
        if record is not None:
            record = fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
                                      [record[0], record[1]])
            records.append(record)
        else:
            exhausted = True

        # write records when it reaches the maximum number of records to
        # keep in memory or there are no records left to fetch
        if len(records) == limit or exhausted is True:
            fo.write_lines(records, output_file, write_mode='a')
            total_selected += len(records)
            records = []

    return total_selected


def split_seqcount(fasta_path, output_directory, max_seqs):
    """Split a FASTA file based on a maximum number of sequences per file.

    Parameters
    ----------
    fasta_path : str
        Path to a FASTA file.
    output_directory : str
        Path to the output directory.
    max_seqs : int
        Split FASTA file into files with a maximum number
        of sequences equal to this value.

    Returns
    -------
    splitted_files : list
        List with paths to the new files that were
        created by splitting the input FASTA file.
    """
    file_count = 1
    exhausted = False
    current_recs = []
    splitted_files = []
    record_generator = sequence_generator(fasta_path)
    while exhausted is False:
        record = next(record_generator, None)
        if record is not None:
            current_recs.append(record)
        else:
            exhausted = True

        if len(current_recs) == max_seqs or exhausted is True:
            if len(current_recs) > 0:
                file_path = fo.join_paths(output_directory,
                                          ['seqcount{0}.fasta'.format(file_count)])
                seqids = (rec.id for rec in current_recs)
                splitted_files.append([file_path, seqids])
                write_records(current_recs, file_path)
                current_recs = []
                file_count += 1

    return splitted_files


def split_seqlength(fasta_path, output_directory, length_cutoff):
    """Split a FASTA file based on a sequence length threshold.

    Parameters
    ----------
    fasta_path : str
        Path to a FASTA file.
    output_directory : str
        Path to the output directory.
    length_cutoff : int
        Sequence length threshold used to split the FASTA file.

    Returns
    -------
    List with two tuples: the first contains the path to a FASTA
    file with sequences equal or above the length cutoff and the
    list of sequence identifiers of those sequences or is None if
    there were no sequences above the cutoff. The second tuple
    contains the same data for the sequences below the cutoff or
    is None if there were no sequences below the cutoff.
    """
    length_values = sequence_lengths(fasta_path)
    below_cutoff = [seqid for seqid, length in length_values.items()
                    if length < length_cutoff]
    above_cutoff = list(set(length_values) - set(below_cutoff))

    fasta_index = index_fasta(fasta_path)
    if len(above_cutoff) > 0:
        above_outfile = fo.join_paths(output_directory, ['above_cutoff.fasta'])
        above_count = get_sequences_by_id(fasta_index, above_cutoff, above_outfile)
        above_seqids = (seqid for seqid in above_cutoff)
        above_data = [above_outfile, above_seqids]
    else:
        above_data = None

    # file has sequences shorter than cutoff value
    if len(below_cutoff) > 0:
        below_outfile = fo.join_paths(output_directory, ['below_cutoff.fasta'])
        below_count = get_sequences_by_id(fasta_index, below_cutoff, below_outfile)
        below_seqids = (seqid for seqid in below_cutoff)
        below_data = [below_outfile, below_seqids]
    else:
        below_data = None

    return [above_data, below_data]


def fasta_stats(fasta_file):
    """Determine the number of sequences in a FASTA file and length stats.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file.

    Returns
    -------
    fasta_file : str
        Path to the FASTA file.
    total_seqs: int
        Total number of records in the FASTA file.
    mean_length: float
        Mean sequence length.
    """
    seq_lengths = sequence_lengths(fasta_file)
    min_length = min(seq_lengths.values())
    max_length = max(seq_lengths.values())
    mean_length = sum(seq_lengths.values())/len(seq_lengths)
    total_seqs = len(seq_lengths)

    return [fasta_file, total_seqs, min_length, max_length, mean_length]


def translate_fasta(input_fasta, output_directory, translation_table):
    """Translate DNA sequences in a FASTA file.

    Parameters
    ----------
    input_fasta : str
        Path to the FASTA file that contains the DNA sequences
        to translate.
    output_directory : str
        Path to the output directory where the FASTA file with
        protein sequences will be written to.
    translation_table : int
        Genetic code used to translate DNA sequences.

    Returns
    -------
    input_fasta : str
        Path to the input FASTA file.
    protein_file : str
        Path to the FASTA file that contains the translated
        sequences.
    Also returns the number of sequences that were translated
    successfully.
    """
    records = sequence_generator(input_fasta)
    translated_records = [[rec.id,
                           str(sm.translate_dna(str(rec.seq),
                                                translation_table, 0)[0][0])]
                          for rec in records]
    translated_lines = fasta_lines(ct.FASTA_RECORD_TEMPLATE,
                                   translated_records)

    basename = fo.file_basename(input_fasta, True).replace('.fasta', '_protein.fasta')
    protein_file = fo.join_paths(output_directory, [basename])

    fo.write_lines(translated_lines, protein_file)

    return [input_fasta, protein_file, len(translated_records)]
