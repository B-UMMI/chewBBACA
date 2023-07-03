#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with gene prediction
with Prodigal.

Code documentation
------------------
"""


import os
import subprocess

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


def create_orf_finder(training_data, closed, mask, meta):
    """
    """

    orf_finder = pyrodigal.OrfFinder(training_info=training_data, closed=closed, mask=mask, meta=meta)

    return orf_finder


def train_orf_finder(orf_finder, sequences, translation_table):
    """
    """
    orf_finder.train(sequences, translation_table=translation_table)

    return orf_finder


def read_training_file(training_file):
    """
    """
    with open(training_file, 'rb') as infile:
        training_data = pyrodigal.TrainingInfo.load(infile)

    return training_data


def predict_genes(orf_finder, sequence):
    """
    """
    genes = orf_finder.find_genes(sequence)

    return genes


def get_gene_info(contig_id, genome_id, protid, genes):
    """
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
    """
    """
    fasta_sequences = []
    for gene in gene_info:
        fasta_str = ct.FASTA_CDS_TEMPLATE.format(gene[2], gene[6], gene[1])
        fasta_sequences.append(fasta_str)
    fo.write_lines(fasta_sequences, output_file)


def write_coordinates_pickle(gene_info, contig_sizes, output_file):
    """
    """
    gene_coordinates = {}
    for gene in gene_info:
        gene_coordinates.setdefault(gene[0], []).append(gene[2:])
    fo.pickle_dumper([gene_coordinates, contig_sizes], output_file)


def predict_genome_genes(input_file, output_directory, orf_finder, translation_table):

    # Get genome unique identifier
    genome_basename = input_file[1]
    records = fao.sequence_generator(input_file[0])
    records = {rec.id: bytes(rec.seq) for rec in records}
    contig_sizes = {recid: len(sequence) for recid, sequence in records.items()}

    # Train based on input sequences
    # Only train if object does not contain training info and if it won't run in meta mode
    if orf_finder.training_info is None and orf_finder.meta is False:
        orf_finder = train_orf_finder(orf_finder, *records.values(), translation_table)

    # Predict genes for all input contigs
    contig_genes = {}
    for recid, sequence in records.items():
        genes = orf_finder.find_genes(sequence)
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
        fasta_outfile = fo.join_paths(output_directory, [f'{genome_basename}.fasta'])
        write_gene_fasta(gene_info, fasta_outfile)

        # Save gene coordinates and contig sizes to pickle
        coordinates_outfile = fo.join_paths(output_directory, [f'{genome_basename}.coordinates'])
        write_coordinates_pickle(gene_info, contig_sizes, coordinates_outfile)

    return [input_file, total_genome, fasta_outfile, coordinates_outfile]
