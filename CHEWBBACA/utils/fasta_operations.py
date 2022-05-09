#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to work with FASTA files.

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
                       sequence_manipulation as sm,
                       blast_wrapper as bw,
                       constants as ct)
except:
    from CHEWBBACA.utils import (file_operations as fo,
                                 iterables_manipulation as im,
                                 sequence_manipulation as sm,
                                 blast_wrapper as bw,
                                 constants as ct)


def sequence_generator(input_file):
    """ Creates a SeqRecord iterator.

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
    """ Creates index to retrieve data from a FASTA file.

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
    """ Imports sequences from a FASTA file.

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
    """ Counts the number of sequences in a FASTA file.

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
        total_seqs = len([l for l in infile if l.startswith('>')])

    return total_seqs


def sequence_lengths(fasta_file):
    """ Read Fasta file and create dictionary
        with mapping between sequence identifiers
        and sequence lengths.

    Parameters
    ----------
    fasta_file : str
        Path to a Fasta file.

    Returns
    -------
    lengths : dict
        Dictionary with sequence identifiers
        as keys and sequence lengths as values.
    """

    sequences = sequence_generator(fasta_file)
    lengths = {rec.id: len(rec.seq) for rec in sequences}

    return lengths


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
            new_rec = '>{0}\n{1}'.format(new_id, sequence)
            seqs.append(new_rec)
            start += 1
        elif record is None:
            exausted = True

        if len(seqs) == limit or exausted is True:
            fo.write_lines(seqs, output_fasta, write_mode='a')
            seqs = []

    return ids_map


def fasta_str_record(record_template, record_data):
    """ Creates the string representation of a FASTA record.

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
    """ Creates a list with FASTA records.

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

    seqs_lines = [fasta_str_record(template, *arg) for arg in records_data]

    return seqs_lines


def validate_fasta(file_path):
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

    records = SeqIO.parse(file_path, 'fasta')

    # returns False is it was not a FASTA file
    return any(records)


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


def sequences_lengths_hash(fasta_file):
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


def get_sequences_by_id(sequences, seqids, output_file, limit=50000):
    """ Retrieves sequences from a FASTA file indexed with
        SeqIO.index or from a dictionary with sequence
        identifiers as keys and sequences as values.

    Parameters
    ----------
    sequences : dict or Bio.File._IndexedSeqFileDict
        Dictionary with seqids as keys and sequences
        as values or a Fasta file index created with
        BioPython.
    seqids : list
        List with the identifiers of the sequences
        that should be retrieved.
    output_file : str
        Path to the FASTA file to which selected
        sequences will be saved.
    limit : int
        Maximum number of sequences that will be
        kept in memory at a time (to avoid keeping
        huge datasets in memory).

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
            record = fasta_str_record(record[0], record[1])
            records.append(record)
        else:
            exhausted = True

        # write records when it reaches the maximum number fo records to
        # keep in memory or there are no records left to fetch
        if len(records) == limit or exhausted is True:
            fo.write_lines(records, output_file, write_mode='a')
            total_selected += len(records)
            records = []

    return total_selected


def exclude_sequences_by_id(sequences, identifiers, output_file):
    """ Creates a FASTA file with the sequences whose
        identifiers are not in the input list.

    Parameters
    ----------
    sequences : dict or Bio.File._IndexedSeqFileDict
        Dictionary with seqids as keys and sequences
        as values or a Fasta file index created with
        BioPython.
    identifiers : list
        List with the sequence identifiers that should
        not be included in the output Fasta file.
    output_file : str
        Path to output file created to store selected
        sequences.

    Returns
    -------
    total_selected : int
        Total number of sequences written to the output
        file.
    """

    selected_ids = (rec for rec in sequences if rec not in identifiers)
    total_selected = get_sequences_by_id(sequences, selected_ids, output_file)

    return total_selected


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
    filenames : gen
        Generator with names to attribute to new files.

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


def get_self_scores(fasta_file, output_directory, blast_threads,
                    blastp_path, makeblastdb_path):
    """ Aligns a set of sequences against itself to determine
        the raw score of the self-alignment.

    Parameters
    ----------
    fasta_file : str
        Path to a FASTA file with protein sequences.
    output_directory : str
        Path to the directory where intermediate files
        will be created.
    blast_threads : int
        Number of threads for BLASTp execution.
    blastp_path : str
        Path to the BLASTp executable.
    makeblastdb_path : str
        Path to the makeblastdb executable.

    Returns
    -------
    self_lines_ids : dict
        Dictionary with sequences identifiers as keys
        and the BLASTp raw score from self-alignment.
    """

    basename = fo.file_basename(fasta_file, suffix=False)

    integer_seqids = fo.join_paths(output_directory,
                                   ['{0}_int.fasta'.format(basename)])
    ids_dict = integer_headers(fasta_file, integer_seqids)

    blastdb = fo.join_paths(output_directory, ['{0}_db'.format(basename)])
    stderr = bw.make_blast_db(makeblastdb_path, integer_seqids,
                              blastdb, 'prot')

    blastout = fo.join_paths(output_directory, ['self_blastout.tsv'])
    self_results = bw.run_blast(blastp_path, blastdb, integer_seqids,
                                blastout, threads=blast_threads,
                                max_targets=1)

    self_lines = fo.read_tabular(blastout)
    self_lines_ids = {ids_dict[l[0]]: l[-1] for l in self_lines}

    return self_lines_ids


def translate_fasta(input_fasta, output_directory, translation_table):
    """ Translates DNA sequences in a FASTA file.

    Parameters
    ----------
    input_fasta : str
        Path to the FASTA file that contains the DNA
        sequences to translate.
    output_directory : str
        Path to the output directory where the FASTA file
        with protein sequences will be writen to.
    translation_table : int
        Genetic code used to translate DNA sequences.

    Returns
    -------
    protein_file : str
        Path to the FASTA file with translated sequences.
    """

    records = import_sequences(input_fasta)
    translated_records = [[seqid, str(sm.translate_dna(seq, translation_table, 0)[0][0])]
                          for seqid, seq in records.items()]
    translated_lines = fasta_lines(ct.FASTA_RECORD_TEMPLATE, translated_records)

    basename = fo.file_basename(input_fasta, True).replace('.fasta', '_protein.fasta')
    protein_file = fo.join_paths(output_directory, [basename])

    fo.write_lines(translated_lines, protein_file)

    return [input_fasta, protein_file]


def determine_self_scores(work_directory, fasta_file, makeblastdb_path,
                          blast_path, db_type, blast_threads):
    """ Determines the self-alignment raw score for the
        sequences in a FASTA file.

    Parameters
    ----------
    work_directory : str
        Path to the working directory.
    fasta_file : str
        Path to the FASTA file that contains the sequences.
    makeblastdb_path : str
        Path to the 'maskeblastdb' executable.
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
        the sequence length and self-alignment raw scores
        as values.
    """

    # change identifiers to shorten and avoid BLAST error related with sequence header length
    output_fasta = fo.join_paths(work_directory, [fo.file_basename(fasta_file, False)+'_intids.fasta'])
    ids_map = integer_headers(fasta_file, output_fasta, start=1, limit=5000)

    blast_db = fo.join_paths(work_directory, [fo.file_basename(output_fasta, False)])
    # will not work if file contains duplicates
    db_stderr = bw.make_blast_db(makeblastdb_path, output_fasta, blast_db, db_type)

    if len(db_stderr) > 0:
        return db_stderr

    output_blast = fo.join_paths(work_directory, [fo.file_basename(fasta_file, False)+'_blastout.tsv'])
    blastp_stderr = bw.run_blast(blast_path, blast_db, output_fasta,
                                 output_blast, threads=blast_threads)

    current_results = fo.read_tabular(output_blast)
    # get raw score and sequence length
    # multiply by 3 to get DNA sequence length and add 3 to count stop codon
    self_scores = {ids_map[l[0]]: ((int(l[3])*3)+3, float(l[6]))
                   for l in current_results
                   if l[0] == l[4]}

    return self_scores
