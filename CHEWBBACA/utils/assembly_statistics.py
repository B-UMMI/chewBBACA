#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script computes assembly statistics.
"""


import math
import statistics
from pathlib import Path

from Bio.SeqUtils import gc_fraction

try:
	from utils import (fasta_operations as fao,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (fasta_operations as fao,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


def calc_n50(contig_sizes):
    """Calculate the assembly N50.

    Parameters
    ----------
    contig_sizes : list
        List containing the size of the contigs.

    Returns
    -------
    cl : int
        Assembly N50.
    """
    # Sort the contig sizes in descending order
    contig_sizes.sort(reverse=True)

    # Calculate n50
    total_bp = sum(contig_sizes)
    n50_threshold = total_bp * 0.5
    for cs in contig_sizes:
        total_bp -= cs
        if total_bp <= n50_threshold:
            # Return the value of the largest contig equal to or below n50 threshold
            return cs


def determine_assembly_stats(input_file):
    """Analyse a set of genome assemblies.

    Parameters
    ----------
    input_file : list
        Path to a genome assembly in FASTA format.

    Returns
    -------
    stats : dict
        Dictionary with the input basename as key and a list containing 
        the sample identifier, number of contigs, average contig size, 
        N50, total assembly length, GC content, and total missing data 
        as values.
    """
    # Get the sample basename/ID from the file
    basename = Path(input_file).stem

    # Import FASTA records
    records = fao.import_sequences(input_file)

    # Calculate the weighted GC content
    records_gc_content = [gc_fraction(seq) for seqid, seq in records.items()]
    records_lengths = [len(seq) for seqid, seq in records.items()]
    total_length = sum(records_lengths)
    weights = [length/total_length for length in records_lengths]
    gc_content = statistics.fmean(records_gc_content, weights)
    gc_content = round(gc_content, 2)

    # Get the total number of contigs in the assembly file
    nr_contigs = len(records)

    # Calculate the average contig size
    avg_size = statistics.mean(records_lengths)
    avg_size = round(avg_size, 2)

    # Calculate the total assembly length
    total_length = sum(records_lengths)

    # Calculate the N50
    n50 = calc_n50(records_lengths)

    # Determine missing data
    missing_data = sum([seq.count('N') for seqid, seq in records.items()])

    stats = {basename: [nr_contigs, avg_size, n50, total_length, gc_content, missing_data]}

    return stats


def compute_assembly_stats(input_files, cpu_cores):

    # Sort inputs
    input_files = sorted(input_files)
    # Create inputs for multiprocessing
    stats_inputs = im.divide_list_into_n_chunks(input_files, len(input_files))
    common_args = []
    stats_inputs = im.multiprocessing_inputs(stats_inputs, common_args, determine_assembly_stats)
    stats_results = mo.map_async_parallelizer(stats_inputs, mo.function_helper, cpu_cores, show_progress=True, pool_type='threadpool')
    # Merge into a single dictionary
    assembly_statistics = im.merge_dictionaries(stats_results)

    return assembly_statistics
