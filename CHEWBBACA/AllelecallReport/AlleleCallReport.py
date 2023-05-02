#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


Code documentation
------------------
""" 


import os
import math
import json
import shutil
import pickle
import statistics
from collections import Counter

import pandas as pd
from Bio.Align.Applications import MafftCommandline

try:
    from utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo)


def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)


def count_lines(file):
    """
    """
    with open(file, 'rb') as infile:
        c_generator = _count_generator(infile.raw.read)
        # count each \n
        count = sum(buffer.count(b'\n') for buffer in c_generator)

    return count


input_files = '/home/rmamede/Desktop/test_chewbbaca320/spneumo_results'
schema_directory = '/home/rmamede/Desktop/test_chewbbaca320/spneumo_schema/schema_seed'
output_directory = '/home/rmamede/Desktop/AlleleCallReport_test'
cpu_cores = 6
annotations = '/home/rmamede/Desktop/test_chewbbaca320/spneumo_schema/spneumo_annotations/schema_seed_annotations.tsv'
def main(input_files, schema_directory, output_directory, annotations,
         cpu_cores):

    # Create temp directory
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # Read "results_statistics.tsv" to get class counts per sample
    sample_statistics_file = fo.join_paths(input_files,
                                           ['results_statistics.tsv'])
    sample_counts = pd.read_csv(sample_statistics_file, delimiter='\t')

    # Get total number of samples
    total_samples = len(sample_counts)
    # Get sample identifiers
    sample_ids = sample_counts['FILE'].tolist()

    # Data for stacked bar plot with class counts per sample
    sample_data = []
    for c in ct.ALLELECALL_CLASSIFICATIONS:
        sample_data.append(sample_counts[c].tolist())

    # Read "loci_summary_stats.tsv" to get class counts per locus
    loci_statistics_file = fo.join_paths(input_files,
                                    ['loci_summary_stats.tsv'])
    loci_counts = pd.read_csv(loci_statistics_file, delimiter='\t')

    # Get total number of loci
    total_loci = len(loci_counts)
    # Get loci identifiers
    loci_ids = loci_counts['Locus'].tolist()

    # Data for stacked bar plot with class counts per locus
    loci_data = []
    for c in ct.ALLELECALL_CLASSIFICATIONS:
        loci_data.append(loci_counts[c].tolist())

### Sample Stats

    sample_stats = {}
    cds_coordinates_file = fo.join_paths(input_files, ['cds_coordinates.tsv'])
    with open(cds_coordinates_file, 'r') as infile:
        header = infile.__next__()
        sample_lines = []
        for line in infile:
            current_line = line.strip().split('\t')
            if len(sample_lines) == 0:
                sample_lines.append(current_line)
            else:
                if current_line[0] == sample_lines[0][0]:
                    sample_lines.append(current_line)
                else:
                    sample_file = fo.join_paths(temp_directory, [f'{sample_lines[0][0]}'])
                    sample_outlines = [im.join_list(line, '\t') for line in sample_lines]
                    fo.write_lines(sample_outlines, sample_file)
                    sample_lines = [current_line]

            sample_stats.setdefault(current_line[0], [[], 0])
            sample_stats[current_line[0]][1] += 1
            if current_line[1] not in sample_stats[current_line[0]][0]:
                sample_stats[current_line[0]][0].append(current_line[1])
        
        sample_file = fo.join_paths(temp_directory, [f'{sample_lines[0][0]}'])
        sample_outlines = [im.join_list(line, '\t') for line in sample_lines]
        fo.write_lines(sample_outlines, sample_file)

    # Get classification sum per sample
    # Can I get the exact number of classified CDSs???
    for i, sample in enumerate(sample_stats):
        sample_line = sample_counts.iloc[i]
        valid = sum(sample_line[1:3])
        invalid = sum(sample_line[3:11])
        total_cds = sample_stats[sample][1]
        classified_loci = valid + invalid
        classified_proportion = round(classified_loci/total_cds, 3)
        classified_proportion_schema = round(classified_loci/total_loci, 3)
        sample_stats[sample].append(int(classified_loci))
        sample_stats[sample].append(float(classified_proportion))
        sample_stats[sample].append(float(classified_proportion_schema))
        sample_stats[sample].append(int(valid))
        sample_stats[sample].append(int(invalid))

    sample_stats_columns = ['Sample', 'Total Contigs', 'Total CDSs',
                            'Classified Loci', 'Proportion of Classified CDSs',
                            'Proportion of Classified Loci',
                            'Valid Classes', 'Invalid Classes']
    
    sample_stats_rows = []
    for k, v in sample_stats.items():
        sample_stats_rows.append([k, len(v[0]), *v[1:]])

### Loci Stats

    loci_stats = []
    for i in range(len(loci_counts)):
        locus_line = loci_counts.iloc[i]
        locus_id = locus_line['Locus']
        total_cds = int(locus_line['Total_CDS'])
        valid = sum(locus_line[1:3])
        invalid = sum(locus_line[3:11])
        sample_proportion = round(valid/total_samples, 3)
        loci_stats.append([locus_id, int(total_cds), int(valid),
                           int(invalid), float(sample_proportion)])

    loci_stats_columns = ['Locus', 'Total CDSs', 'Valid Classes',
                          'Invalid Classes', 'Proportion Samples']

### Loci Annotations

    # Read loci annotations from TSV file
    annotation_values = [[], []]
    if annotations is not None:
        annotation_lines = fo.read_tabular(annotations)
        # Add file columns
        annotation_values[0].extend(annotation_lines[0])
        # Add sublists with lines
        # Only keep lines for loci in the schema
        loci_annotations = [line for line in annotation_lines[1:]
                            if line[0] in loci_ids]
        annotation_values[1].extend(loci_annotations)





    # Add total CDSs, total classified CDSs, total per class to summary table
    # Count number of CDSs
    cds_coordinates_file = fo.join_paths(input_files, ['cds_coordinates.tsv'])
    total_cds = count_lines(cds_coordinates_file) - 1

    # Count total number of classified CDSs per class
    loci_sums = loci_counts[loci_counts.columns[1:]].sum(axis=0)
    loci_sums = loci_sums.tolist()

    summary_columns = ('Total Samples\tTotal Loci\tTotal CDSs\tTotal '
                       'Classified\tEXC\tINF\tPLOT3\tPLOT5\tLOTSC\tNIPH\t'
                       'NIPHEM\tALM\tASM\tPAMA\tLNF')
    summary_columns = summary_columns.split('\t')

    summary_rows = [total_samples, total_loci, total_cds, loci_sums[-1]]
    summary_rows.extend(loci_sums[:-1])

    report_data = {"summaryData": [{"columns": summary_columns},
                                   {"rows": [summary_rows]}],
                   "sample_ids": sample_ids,
                   "sample_data": sample_data,
                   "loci_ids": loci_ids,
                   "loci_data": loci_data,
                   "sample_stats": [{"columns": sample_stats_columns},
                                    {"rows": sample_stats_rows}],
                   "loci_stats": [{"columns": loci_stats_columns},
                                  {"rows": loci_stats}],
                   "annotations": [{"columns": annotation_values[0]},
                                   {"rows": annotation_values[1]}],
                   }

    # Write Schema Report HTML file
    report_html = ct.ALLELECALL_REPORT_HTML
    report_html = report_html.format(json.dumps(report_data))
    report_html_file = fo.join_paths(output_directory, ['main_report.html'])
    fo.write_to_file(report_html, report_html_file, 'w', '\n')

    # Create loci reports
    # if loci_reports:

    # Create sample reports
    # if sample_reports:


# Add info alert about valid (EXC+INF) and invalid (PLOT3+PLOT5+NIPH+NIPHEM+ASM+ALM)
