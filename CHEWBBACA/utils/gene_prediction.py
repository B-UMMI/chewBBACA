#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with gene prediction
with Pyrodigal.

Code documentation
------------------
"""


import pyrodigal

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


def create_gene_finder(training_data, closed, mask, meta):
    """Create a Pyrodigal GeneFinder object.

    Parameters
    ----------
    training_data : pyrodigal.TrainingInfo
        A training info instance used to predict genes in single
        mode.
    closed : bool
        True to prevent prediction of partial genes at edges of
        sequences, False otherwise.
    meta: bool
        True to run Prodigal in `meta` mode (uses pre-trained
        profiles).

    Returns
    -------
    gene_finder : pyrodigal.GeneFinder
        A GeneFinder object configured based on provided arguments.
    """
    gene_finder = pyrodigal.GeneFinder(training_info=training_data,
                                       closed=closed,
                                       mask=mask,
                                       meta=meta)

    return gene_finder


def train_gene_finder(gene_finder, sequences, translation_table):
    """Train a Pyrodigal GeneFinder object based on a set of sequences.

    Parameters
    ----------
    gene_finder : pyrodigal.GeneFinder
        A GeneFinder object.
    sequences : list
        Sequences used to train the GeneFinder (list
        of bytes objects).
    translation_table : int
        Translation table to use.

    Return
    ------
    gene_finder : pyrodigal.GeneFinder
        A GeneFinder object configured based on provided arguments.
    """
    gene_finder.train(*sequences, translation_table=translation_table)

    return gene_finder


def read_training_file(training_file):
    """Load training info for Pyrodigal from Prodigal training file.

    Parameters
    ----------
    training_file : str
        Path to Prodigal training file.

    Returns
    -------
    training_data : pyrodigal.TrainingInfo
        The deserialized training info.
    """
    with open(training_file, 'rb') as infile:
        training_data = pyrodigal.TrainingInfo.load(infile)

    return training_data


def get_gene_info(contig_id, genome_id, protid, genes):
    """Get genes information from a pyrodigal.Genes object.

    Parameters
    ----------
    contig_id : str
        The unique identifier of the sequence/contig.
    genome_id : str
        The unique identifier of the genome/file.
    protid : int
        The integer identifier to attriute to the first gene.
    genes : pyrodigal.Genes
        The list of genes predicted by Prodigal.

    Returns
    -------
    gene_info : list
        List with one sublist per gene predicted. Each sublist
        includes the sequence SHA256 hash, the DNA sequence, the
        genome identifier, the contig identifier, the start position
        in the sequence, the end position, the integer identifier and
        the strand the gene was identified in.
    protid : int
        The integer identifier to attribute to the first gene
        in the next sequence/contig.
    """
    gene_info = []
    for gene in genes:
        sequence = gene.sequence()
        sequence_hash = im.hash_sequence(sequence)
        gene_info.append([sequence_hash, sequence, genome_id, contig_id,
                          str(gene.begin), str(gene.end), str(protid),
                          str(gene.strand)])
        protid += 1

    return gene_info, protid


def write_gene_fasta(gene_info, output_file):
    """Write a FASTA file based on the results returned by `get_gene_info`.

    Parameters
    ----------
    gene_info : list
        List with the data for the genes returned by `get_gene_info`.
    output_file : str
        Path to the output FASTA file.
    """
    fasta_sequences = []
    for gene in gene_info:
        fasta_str = ct.FASTA_CDS_TEMPLATE.format(gene[2], gene[6], gene[1])
        fasta_sequences.append(fasta_str)
    fo.write_lines(fasta_sequences, output_file)


def write_coordinates_pickle(gene_info, contig_sizes, output_file):
    """Write gene coordinates to a pickle file.

    Parameters
    ----------
    gene_info : list
        List with the data for the genes returned by `get_gene_info`.
    contig_sizes : dict
        Dictionary with contig/sequence identifiers as keys and
        contig/sequence size as values.
    output_file : str
    Path to the output file.
    """
    gene_coordinates = {}
    for gene in gene_info:
        gene_coordinates.setdefault(gene[0], []).append(gene[2:])
    fo.pickle_dumper([gene_coordinates, contig_sizes], output_file)


def predict_genome_genes(input_file, output_directory, gene_finder,
                         translation_table):
    """Predict genes for sequences in a FASTA file.

    Parameters
    ----------
    input_file : str
        Path to the FASTA file.
    output_directory : str
        Path to the output_directory to store files with
        the results.
    gene_finder : pyrodigal.GeneFinder
        A GeneFinder object.
    translation_table : int
        Translation table used to configure the GeneFinder
        (None type if the GeneFinder does not need to be
         configured).

    Returns
    -------
    input_file : str
        Path to the input FASTA file.
    total_genome : int
        Total number of genes predicted.
    fasta_outfile : str
        Path to the output FASTA file that contains the
        predited gene sequences.
    coordinates_outfile : str
        Path to the output pickle file that contains the gene
        coordinates and contig size data.
    """
    # Get genome unique identifier
    genome_basename = input_file[1]
    records = fao.sequence_generator(input_file[0])
    records = {rec.id: bytes(rec.seq) for rec in records}
    contig_sizes = {recid: len(sequence)
                    for recid, sequence in records.items()}

    # Train based on input sequences
    # Only train if object does not contain training info
    # and if it won't run in meta mode
    if gene_finder.training_info is None and gene_finder.meta is False:
        gene_finder = train_gene_finder(gene_finder,
                                        records.values(),
                                        translation_table)

    # Predict genes for all input contigs
    contig_genes = {}
    for recid, sequence in records.items():
        genes = gene_finder.find_genes(sequence)
        contig_genes[recid] = genes

    # Extract data from Gene objects
    protid = 1
    gene_info = []
    for recid, genes in contig_genes.items():
        data = get_gene_info(recid, genome_basename, protid, genes)
        gene_info.extend(data[0])
        protid = data[1]

    total_genome = len(gene_info)
    fasta_outfile = None
    coordinates_outfile = None
    if total_genome > 0:
        # Create FASTA file with DNA sequences
        fasta_outfile = fo.join_paths(output_directory,
                                      [f'{genome_basename}.fasta'])
        write_gene_fasta(gene_info, fasta_outfile)

        # Save gene coordinates and contig sizes to pickle
        coordinates_outfile = fo.join_paths(output_directory,
                                            [f'{genome_basename}.coordinates'])
        write_coordinates_pickle(gene_info, contig_sizes, coordinates_outfile)

    return [input_file, total_genome, fasta_outfile, coordinates_outfile]
