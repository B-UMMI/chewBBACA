#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with the execution
of the BLAST software (https://www.ncbi.nlm.nih.gov/books/NBK279690/).

Code documentation
------------------
"""


import subprocess

try:
    from utils import (constants as ct,
                       iterables_manipulation as im)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 iterables_manipulation as im)


def make_blast_db(makeblastdb_path, input_fasta, output_path, db_type,
                  ignore=None):
    """Create a BLAST database.

    Parameters
    ----------
    makeblastdb_path : str
        Path to the 'makeblastdb' executable.
    input_fasta : str
        Path to the FASTA file that contains the sequences that
        will be added to the BLAST database.
    output_path : str
        Path to the directory where the database files will be
        created. Database files will have the same basename as
        the `input_fasta`.
    db_type : str
        Type of the database, nucleotide (nuc) or protein (prot).
    ignore : list or NoneType
        List with BLAST warnings that should be ignored.

    Returns
    -------
    stderr : list
        A list with the warnings and errors raised by BLAST.
    """
    # use '-parse-seqids' to be able to retrieve/align sequences by identifier
    blastdb_cmd = [makeblastdb_path, '-in', input_fasta,
                   '-out', output_path, '-parse_seqids',
                   '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = makedb_cmd.stderr.readlines()

    # ignore errors/warnings provided to `ignore`
    if len(stderr) > 0:
        stderr = im.decode_str(stderr, 'utf8')
        if ignore is not None:
            stderr = im.filter_list(stderr, ignore)

    return stderr


def determine_blast_task(sequences, blast_type='blastp'):
    """Determine the type of BLAST task to execute.

    It is necessary to define the BLAST task if any of the
    sequences to align is shorter that 50 base pairs for
    BLASTn or 30 amino acids for BLASTp.

    Parameters
    ----------
    sequences : list
        List that contains strings representing DNA or
        protein sequences.
    blast_type : str
        Used to define the type of application, 'blastn'
        or 'blastp'.

    Returns
    -------
    blast_task : str
        A string that indicates the type of BLAST task to
        execute based on the minimum sequence size.

    Notes
    -----
    More information about the task option at:
        https://www.ncbi.nlm.nih.gov/books/NBK569839/
    """
    # get sequence length threshold for BLAST application
    length_threshold = ct.BLAST_TASK_THRESHOLD[blast_type]
    sequence_lengths = [len(p) for p in sequences]
    minimum_length = min(sequence_lengths)
    if minimum_length < length_threshold:
        blast_task = '{0}-short'.format(blast_type)
    else:
        blast_task = blast_type

    return blast_task


def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None, blast_task=None,
              max_targets=None, ignore=None):
    """Execute BLAST to align sequences against a BLAST database.

    Parameters
    ----------
    blast_path : str
        Path to the BLAST application executable.
    blast_db : str
        Path to the BLAST database.
    fasta_file : str
        Path to the FASTA file with sequences to align against
        the BLAST database.
    blast_output : str
        Path to the file that will be created to store the
        results.
    max_hsps : int
        Maximum number of High Scoring Pairs per pair of aligned
        sequences.
    threads : int
        Number of threads/cores used to run BLAST.
    ids_file : str
        Path to a file with sequence identifiers, one per line.
        Sequences will only be aligned to the sequences in the
        BLAST database that match any of the identifiers in this
        file.
    blast_task : str
        Type of BLAST task.
    max_targets : int
        Maximum number of target/subject sequences to align
        against.
    ignore : list or None
        List with BLAST warnings that should be ignored.

    Returns
    -------
    stderr : list
        A list with the warnings and errors raised by BLAST.
    """
    # do not retrieve hits with high probability of occuring by change (-evalue=0.001)
    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', ct.BLAST_DEFAULT_OUTFMT,
                  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
                  '-evalue', '0.001']

    # add file with list of sequence identifiers to align against
    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])
    # add type of BLASTp or BLASTn task
    if blast_task is not None:
        blast_args.extend(['-task', blast_task])
    # add maximum number of target sequences to align against
    if max_targets is not None:
        blast_args.extend(['-max_target_seqs', str(max_targets)])

    blast_proc = subprocess.Popen(blast_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()

    # ignore errors/warnings provided to `ignore`
    if len(stderr) > 0:
        stderr = im.decode_str(stderr, 'utf8')
        if ignore is not None:
            stderr = im.filter_list(stderr, ignore)

    return stderr
