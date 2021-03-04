#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to work with FASTA files and
FASTA records.

Code documentation
------------------
"""


import os
import hashlib

from Bio import SeqIO
from Bio.SeqIO import FastaIO

try:
    from utils import (file_operations as fo,
                       iterables_manipulation as im,
                       constants as ct)
except:
    from CHEWBBACA.utils import (file_operations as fo,
                                 iterables_manipulation as im,
                                 constants as ct)


def count_sequences(fasta_file):
    """ Counts the number of sequences in a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file.

        Returns
        -------
        total_seqs : int
            Number of sequences in the input FASTA file.
    """

    records = SeqIO.parse(fasta_file, 'fasta')
    total_seqs = len(list(records))

    return total_seqs


def write_records(records, output_file):
    """ Writes FASTA records (BioPython SeqRecord) to a file.

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


def integer_headers(input_fasta, output_fasta, start=1, limit=5000):
    """ Switches FASTA records headers in a file by integer
        values.

        Parameters
        ----------
        input_fasta : str
            Path to the a FASTA file.
        output_fasta : str
            Path to the output file with modified headers.
        start : int
            Integer value of first identifier.
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
    seq_generator = SeqIO.parse(input_fasta, 'fasta')
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            new_id = 'seq_{0}'.format(start)
            ids_map[new_id] = record.id
            sequence = str(record.seq)
            new_rec = '>{0}\n{1}'.format(new_id, sequence)
            seqs.append(new_rec)
            start += 1
        elif record is None:
            exausted = True

        if len(seqs) == limit or exausted is True:
            fo.write_lines(seqs, output_fasta)
            seqs = []

    return ids_map


def create_fasta_lines(sequences, prefix):
    """ Creates FASTA records in string format.

        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys
            and sequences as values.
        prefix : str
            Prefix to include in sequences headers.

        Returns
        -------
        lines : list
            List with Fasta records in string format.
    """

    template = '>{0}-protein{1}\n{2}'

    lines = [template.format(prefix, seqid, sequence)
             for seqid, sequence in sequences.items()]

    return lines


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.

        Returns
        -------
        seqs_dict : dict
            Dictionary that has sequences ids as keys and
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def is_fasta(file_path):
    """ Checks if a file is a FASTA file.

        Parameters
        ----------
        file_path : str
            Path to the file.

        Returns
        -------
        True if file has valid FASTA format,
        False otherwise.
    """

    with open(file_path, 'r') as handle:
        try:
            fasta = SeqIO.parse(handle, 'fasta')
        except:
            fasta = [False]

        return any(fasta)


def filter_non_fasta(files):
    """ Creates a new list of files names/paths that only contains
        FASTA files.

        Parameters
        ----------
        files : list
            A list with files names/paths.

        Returns
        -------
        fasta_files : list
            List with files names/paths that have FASTA
            format.
    """

    fasta_files = [file for file in files if is_fasta(file) is True]

    return fasta_files


def fasta_lines(identifiers, sequences_dictionary):
    """ Creates list with line elements for a FASTA file based
        on the sequence identifiers passed.

        Parameters
        ----------
        identifiers : list
            A list with the identifiers of sequences that will be
            included in the list.
        sequences_dictionary : dict
            A dictionary with sequence identifiers as keys and
            sequences as values.

        Returns
        -------
        seqs_lines : list
            A list with strings representing the header of
            the sequence and the sequence for each of the specified
            sequence identifiers.
    """

    seqs_lines = ['>{0}\n{1}\n'.format(seqid, sequences_dictionary[seqid])
                  for seqid in identifiers]

    return seqs_lines


def sequences_lengths(fasta_file):
    """ Determines the length of all sequences in a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file with sequences.

        Returns
        -------
        lengths : dict
            Dictionary with the `fasta_file` basename as key and
            a nested dictionary with sequences SHA256 hashes as
            keys and sequences lengths as values.
    """

    basename = os.path.basename(fasta_file)
    lengths = {basename: {hashlib.sha256(str(rec.seq).encode('utf-8')).hexdigest(): len(rec.seq)
                          for rec in SeqIO.parse(fasta_file, 'fasta')}}

    return lengths


def fasta_str_record(seqid, sequence):
    """ Creates the string representation of a FASTA record.

        Parameters
        ----------
        seqid : str
            Sequence identifier to include in the header.
        sequence : str
            String representing DNA or Protein sequence.

        Returns
        -------
        record : str
            String representation of the FASTA record.
    """

    record = '>{0}\n{1}'.format(seqid, sequence)

    return record


def get_sequences_by_id(sequences, seqids, out_file, limit=5000):
    """ Retrieves sequences from an indexed FASTA file.

        Parameters
        ----------
        sequences : dict or Bio.File._IndexedSeqFileDict
            Dictionary with seqids as keys and sequences
            as values or a Fasta file index created with
            BioPython.
        seqids : list
            List with the identifiers of the sequences
            that should be retrieved.
        out_file : str
            Path to the FASTA file to which selected
            sequences will be saved.
        limit : int
            Maximum number of sequences that will be
            kept in memory at a time (to avoid keeping
            huge datasets in memory).

        Returns
        -------
        Creates a file with the sequences that have the
        identifiers in the input list.
    """

    if type(sequences) == dict:
        seqs = [(seqid, sequences[seqid]) for seqid in seqids]
    else:
        seqs = [(seqid, str(sequences[seqid].seq)) for seqid in seqids]

    records = []
    for seq in seqs:
        record = fasta_str_record(seq[0], seq[1])
        records.append(record)

        if len(records) == limit or seq[0] == seqids[-1]:
            lines = im.join_list(records, '\n')
            fo.write_to_file(lines, out_file, 'a', '\n')
            records = []


def split_fasta(fasta_path, output_path, num_seqs, filenames):
    """ Splits a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.
        output_path : str
            Path to the output directory where new FASTA
            files will be created.
        num_seqs : int
            Split FASTA file into files with this number
            of sequences.
        filenames : list
            List with names to attribute to new files.

        Returns
        -------
        splitted_files : list
            List with paths to the new files that were
            created by splitting the input FASTA file.
    """

    splitted_files = []
    current_recs = []
    records = [rec for rec in SeqIO.parse(fasta_path, 'fasta')]
    for record in records:
        current_recs.append(record)
        if len(current_recs) == num_seqs or record.id == records[-1].id:
            file_name = filenames.__next__()
            file_name = im.replace_multiple_characters(file_name, ct.CHAR_REPLACEMENTS)

            new_file = fo.join_paths(output_path,
                                     ['{0}{1}'.format(file_name, '.fasta')])

            splitted_files.append(new_file)

            write_records(current_recs, new_file)

            current_recs = []

    return splitted_files


def gene_seqs_info(fasta_file):
    """ Determines the total number of sequences and the mean
        length of sequences in a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file.

        Returns
        -------
        stats : list
            A list with the path to the FASTA file, the
            total number of records in that file and
            the sequence mean length for the sequences
            in the file.
    """

    seq_generator = SeqIO.parse(fasta_file, 'fasta')
    seqs_lengths = [len(seq) for seq in seq_generator]
    mean_length = sum(seqs_lengths)/len(seqs_lengths)
    total_seqs = len(seqs_lengths)
    stats = [fasta_file, total_seqs, mean_length]

    return stats
