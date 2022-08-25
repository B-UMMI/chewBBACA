#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with gene prediction
with Prodigal and extraction of predicted coding sequences
from FASTA files.

Code documentation
------------------
"""


import os
import subprocess

try:
    from utils import (constants as ct,
                       file_operations as fo,
                       fasta_operations as fao,
                       iterables_manipulation as im)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 iterables_manipulation as im)


def extract_genome_cds(reading_frames, contigs, starting_id):
    """Extract coding sequence from contigs.

    Extracts coding sequences from FASTA files based on the
    start and stop coordinates predicted by Prodigal.

    Parameters
    ----------
    reading_frames : str
        Path to the pickled file with the coordinates to
        extract the CDSs predicted by Prodigal.
    contigs : dict
        Dictionary with contig identifiers as keys and
        contig sequences as values.
    starting_id : int
        Integer identifier attributed to the first CDS
        and that will be incremented to serve as identifier
        for subsequent CDSs.

    Returns
    -------
    coding_sequences : dict
        Dictionary with coding sequences ids as keys and
        coding sequences as values.
    coding_sequences_info : list
        List with a sublist for each extracted CDS. Sublists
        have information about the extracted CDS (identifier
        of the contig the CDS was identified in, start position
        in the contig, stop position in the contig, sequence
        identifier attributed to that CDS and the strand that
        coded for that CDS).
    """
    seqid = starting_id
    coding_sequences = {}
    coding_sequences_info = []
    for contig_id, frames in reading_frames.items():
        sequence = contigs[contig_id]
        # for each start and stop codon in the contig
        for cds in frames:
            # start position is 0-based, stop position is upper-bound exclusive
            start_pos = cds[0]
            stop_pos = cds[1]
            strand = cds[2]
            # extract CDS sequence
            cds_sequence = im.extract_single_cds(sequence, *cds).upper()
            seq_hash = im.hash_sequence(cds_sequence)

            # store CDS with unique id
            coding_sequences[seqid] = cds_sequence

            # store CDS information
            coding_sequences_info.append([contig_id, str(start_pos),
                                          str(stop_pos), str(seqid),
                                          str(strand), seq_hash])

            # increment seqid
            seqid += 1

    return [coding_sequences, coding_sequences_info]


def write_coordinates_tsv(cds_info, genome_id, output_file):
    """Write a TSV file with coding sequence coordinates.

    Parameters
    ----------
    cds_info : list
        List with information about each coding sequence
        identified in the genome (contig identifier,
        CDS start position, CDS stop position, CDS
        identifier and CDS coding strand).
    genome_id : str
        Identifier of the genome to add to first field
        of every new line.
    output_file : str
        Path to the output file to which info will
        be saved.
    """
    # write TSV file
    table_lines = [[genome_id] + protein_info[:-1]
                   for protein_info in cds_info]
    table_lines = [im.join_list(line, '\t') for line in table_lines]
    table_text = im.join_list(table_lines, '\n')
    fo.write_to_file(table_text, output_file, 'a', '\n')


def write_coordinates_pickle(cds_info, contig_lengths, output_file):
    """Save coordinates for coding sequences predicted by Prodigal.

    Parameters
    ----------
    cds_info : list
        List with information about each coding sequence
        identified in the genome (contig identifier,
        CDS start position, CDS stop position, CDS
        identifier and CDS coding strand).
    contig_lengths : dict
        Dictionary with contig identifiers as keys and
        contig lengths as values.
    output_file : str
        Path to the output file to which info will
        be saved.
    """
    # write pickle with CDS hash to CDS info dictionary
    pickle_data = [{}, contig_lengths]
    # create dictionary to map CDS hash to CDS location
    for p in cds_info:
        # make sure to store CDS duplicated in the genome
        pickle_data[0].setdefault(p[-1], []).append(p[:-1])

    fo.pickle_dumper(pickle_data, output_file)


def save_extracted_cds(genome, identifier, orf_file, protein_table, cds_file):
    """Extract coding sequences from a FASTA file.

    Extracts coding sequences from a FASTA file based on
    Prodigal's predictions. Writes coding sequences to a
    FASTA file and information about coding sequences to
    a TSV file.

    Parameters
    ----------
    genome : str
        Path to the FASTA file with the FASTA sequences for
        a genome.
    identifier : str
        Genome identifier to add to FASTA records headers
        and to the first field in the TSV file.
    orf_file : str
        Path to the file with Prodigal results.
    protein_table : str
        Path to the TSV file to which coding sequences
        information will be written.
    cds_file : str
        Path to the FASTA file to which coding sequences
        will be written.

    Returns
    -------
    total_cds : int
        Total number of coding sequences extracted from
        the genome.
    """
    # import contigs for current genome/assembly
    contigs = fao.import_sequences(genome)
    # determine contig lengths
    contig_lengths = {k: len(v) for k, v in contigs.items()}
    # extract coding sequences from contigs
    reading_frames = fo.pickle_loader(orf_file)
    genome_info = extract_genome_cds(reading_frames, contigs, 1)
    # save coding sequences to file
    # create records and write them to file
    cds_data = [[identifier, k, v] for k, v in genome_info[0].items()]
    cds_lines = fao.fasta_lines(ct.FASTA_CDS_TEMPLATE, cds_data)
    fo.write_lines(cds_lines, cds_file, write_mode='a')

    write_coordinates_tsv(genome_info[1], identifier, protein_table)
    pickle_out = os.path.join(os.path.dirname(protein_table), identifier+'.cds_hash')
    write_coordinates_pickle(genome_info[1], contig_lengths, pickle_out)

    total_cds = len(genome_info[0])

    return total_cds


def cds_batch_extractor(genomes, index, prodigal_path, temp_directory):
    """Extract coding sequences from a set of genomes.

    Parameters
    ----------
    input_data : list
        List with a set of paths for FASTA files with
        genomic sequences, followed by the path to the
        directory with files with Prodigal resutls, the
        path to the temporary directory for all files and
        directories that will be read and written and
        an index/identifier to add to the output files
        with coding sequences and coding sequences info.

    Returns
    -------
    A list with the following elements:
        protein_table : str
            Path to the TSV file to which coding sequence
            coordinates was written.
        cds_file : str
            Path to the FASTA file to which coding sequences
            were written.
        batch_total : int
            Total number of coding sequences extracted from
            the set of input genomes.
    """
    protein_table = fo.join_paths(temp_directory,
                                  ['cds_coordinates_{0}.tsv'.format(index)])

    cds_file = fo.join_paths(temp_directory,
                             ['coding_sequences_{0}.fasta'.format(index)])

    batch_total = 0
    for g in genomes:
        # determine Prodigal ORF file path for current genome
        identifier = fo.file_basename(g, False)
        orf_file_path = fo.join_paths(prodigal_path,
                                      ['{0}.cds_coordinates'.format(identifier)])
        total = save_extracted_cds(g, identifier, orf_file_path,
                                   protein_table, cds_file)
        batch_total += total

    return [protein_table, cds_file, batch_total]


def run_prodigal(input_file, translation_table, mode, ptf_path):
    """Execute Prodigal.

    Parameters
    ----------
    input_file : str
        Path to a FASTA file.
    translation_table : int
        Genetic code.
    mode : str
        Prodigal execution mode ('single' is the default,
        'meta' should be used to predict genes from smaller
        contigs).
    ptf_path : str or None
        Path to the training file.

    Returns
    -------
    stdout : bytes
        Prodigal's stdout.
    stderr : bytes
        Prodigal's stderr.
    """
    if ptf_path is not None:
        proc = subprocess.Popen(['prodigal', '-i', input_file, '-c',
                                 '-m', '-g', str(translation_table), '-p',
                                 mode, '-f', 'sco', '-q', '-t', ptf_path],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    elif ptf_path is None:
        proc = subprocess.Popen(['prodigal', '-i', input_file, '-c',
                                 '-m', '-g', str(translation_table), '-p',
                                 mode, '-f', 'sco', '-q'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

    # Read the stdout from Prodigal
    stdout = proc.stdout.readlines()
    stderr = proc.stderr.readlines()

    return [stdout, stderr]


def main(input_file, output_dir, ptf_path, translation_table, mode):

    stdout, stderr = run_prodigal(input_file, translation_table, mode, ptf_path)

    # this has to be changed to include the full file name!
    genome_basename = fo.file_basename(input_file, False)

    if len(stderr) > 0:
        stderr = [line.decode('utf-8').strip() for line in stderr]
        stderr = [line for line in stderr if line != '']
        error = ' '.join(stderr)
        return [input_file, error]

    # Parse output
    lines = [line.decode('utf-8').strip() for line in stdout]

    # determine contigs headers indexes
    contigs_headers = [line for line in lines if 'seqhdr' in line]
    contigs_ids = [line.split('"')[1].split()[0] for line in contigs_headers]
    contigs_idx = [lines.index(line) for line in contigs_headers] + [len(lines)]

    # get CDSs' positions for each contig
    contigs_pos = {contigs_ids[i]: lines[contigs_idx[i]+1:contigs_idx[i+1]]
                   for i in range(len(contigs_ids))}

    # exclude contigs without coding sequences
    contigs_pos = {k: v[1:] for k, v in contigs_pos.items() if len(v) > 1}

    # +/1 for sense, -/0 for antisense
    strand_trans = {'+': 1, '-': 0}

    # split and convert list elements
    contigs_pos = {k: [p.split('_')[1:] for p in v]
                   for k, v in contigs_pos.items()}
    contigs_pos = {k: [[int(p[0])-1, int(p[1]), strand_trans[p[2]]]
                   for p in v] for k, v in contigs_pos.items()}

    total_contigs = {k: len(v) for k, v in contigs_pos.items()}
    total_genome = sum(total_contigs.values())

    if total_genome > 0:
        # save positions in file
        filepath = os.path.join(output_dir, genome_basename + '.cds_coordinates')
        fo.pickle_dumper(contigs_pos, filepath)

    status = [input_file, total_genome]

    return status
