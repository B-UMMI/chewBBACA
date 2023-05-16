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
import subprocess
import statistics
import pandas as pd
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
        multiprocessing_operations as mo,
        distance_matrix as dm)
    from ExtractCgMLST import Extract_cgAlleles
except ModuleNotFoundError:
    from CHEWBBACA.utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo,
        distance_matrix as dm)
    from CHEWBBACA.ExtractCgMLST import Extract_cgAlleles


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


# locus = loci_ids[0]
# schema_directory = '/home/rmamede/Desktop/test_chewbbaca320/spneumo_schema/schema_seed'
# allelic_profiles = allelic_profiles_file
# output_directory = fasta_dir
def compute_locus_statistics(locus, schema_directory, allelic_profiles, output_directory):
    """
    """
    # Read locus column
    df = pd.read_csv(allelic_profiles, usecols=['FILE', locus], delimiter='\t', dtype=str)
    locus_column = df[locus].tolist()
    sample_ids = df['FILE'].tolist()

    # Read locus FASTA file
    locus_fasta = fo.join_paths(schema_directory, [f'{locus}.fasta'])
    alleles = fao.import_sequences(locus_fasta)
    alleles = {k.split('_')[-1]: v for k, v in alleles.items()}

    sequences = []
    for i, allele in enumerate(locus_column):
        clean_class = allele.split('INF-')[-1]
        if clean_class in alleles:
            current_sample = sample_ids[i]
            seq = alleles[clean_class]
            record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [current_sample, seq])
            sequences.append(record)

    fasta_file = fo.join_paths(output_directory, [f'{locus}.fasta'])
    fo.write_lines(sequences, fasta_file)

    return fasta_file


def call_mafft(genefile):
    """Call MAFFT to generate an alignment.

    Parameters
    ----------
    genefile : str
        Path to a FASTA file with the sequences to align.

    Returns
    -------
    Path to the file with the computed MSA if successful, False otherwise.
    """
    try:
        mafft_cline = MafftCommandline(input=genefile,
                                       adjustdirection=True,
                                       treeout=True,
                                       thread=1,
                                       retree=1,
                                       maxiterate=0,
                                       )
        stdout, stderr = mafft_cline()
        path_to_save = genefile.replace(".fasta", "_aligned.fasta")
        with open(path_to_save, "w") as handle:
            handle.write(stdout)

        return path_to_save
    except Exception as e:
        print(e)
        return False


def run_fasttree(alignment_file, tree_file):
    """
    """
    proc = subprocess.Popen(['FastTree', '-fastest', '-nosupport',
                             '-noml', '-out', tree_file, alignment_file],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # Read the stdout from FastTree
    stdout = proc.stdout.readlines()
    stderr = proc.stderr.readlines()

    return [stdout, stderr]


input_files = '/home/rmamede/Desktop/CAML_2023/presentation/lynskey_results'
schema_directory = '/home/rmamede/Desktop/CAML_2023/presentation/spyogenes_schema_chewieNS'
output_directory = '/home/rmamede/Desktop/CAML_2023/presentation/AlleleCallReport_test'
cpu_cores = 6
annotations = '/home/rmamede/Desktop/CAML_2023/presentation/spyogenes_annotations/annotations.tsv'
loci_reports = True
translation_table = 11
def main(input_files, schema_directory, output_directory, annotations,
         cpu_cores):

    # Create temp directory
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # Read "results_statistics.tsv" to get class counts per sample
    sample_statistics_file = fo.join_paths(input_files,
                                           ['results_statistics.tsv'])
    sample_counts = pd.read_csv(sample_statistics_file, delimiter='\t')
    # Sort based on decreasing number of EXC
    sample_counts = sample_counts.sort_values(by=['EXC'], ascending=False)

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
    # Sort based on decreasing number of EXC
    loci_counts = loci_counts.sort_values(by=['EXC'], ascending=False)

    # Get total number of loci
    total_loci = len(loci_counts)
    # Get loci identifiers
    loci_ids = loci_counts['Locus'].tolist()

    # Data for stacked bar plot with class counts per locus
    loci_data = []
    for c in ct.ALLELECALL_CLASSIFICATIONS:
        loci_data.append(loci_counts[c].tolist())

### Sample Stats

    cds_coordinates_dir = fo.join_paths(temp_directory, ['coordinates'])
    fo.create_directory(cds_coordinates_dir)

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
                    sample_file = fo.join_paths(cds_coordinates_dir, [f'{sample_lines[0][0]}'])
                    sample_outlines = [im.join_list(line, '\t') for line in sample_lines]
                    fo.write_lines(sample_outlines, sample_file)
                    sample_lines = [current_line]

            sample_stats.setdefault(current_line[0], [[], 0])
            sample_stats[current_line[0]][1] += 1
            if current_line[1] not in sample_stats[current_line[0]][0]:
                sample_stats[current_line[0]][0].append(current_line[1])
        
        sample_file = fo.join_paths(cds_coordinates_dir, [f'{sample_lines[0][0]}'])
        sample_outlines = [im.join_list(line, '\t') for line in sample_lines]
        fo.write_lines(sample_outlines, sample_file)

    # Get classification sum per sample
    # Can I get the exact number of classified CDSs???
    for i, sample in enumerate(sample_stats):
        sample_line = sample_counts.iloc[i]
        valid = sum(sample_line[1:3])
        invalid = sum(sample_line[3:11])
        total_cds = sample_stats[sample][1]
        identified_loci = valid + invalid
        ## Can I get the exact number of classified CDSs? I think this only counts 1 CDS per NIPH, NIPHEM, ...
        classified_proportion = round(identified_loci/total_cds, 3)
        identified_proportion = round(identified_loci/total_loci, 3)
        sample_stats[sample].append(float(classified_proportion))
        sample_stats[sample].append(int(identified_loci))
        sample_stats[sample].append(float(identified_proportion))
        sample_stats[sample].append(int(valid))
        sample_stats[sample].append(int(invalid))

    sample_stats_columns = ['Sample', 'Total Contigs', 'Total CDSs',
                            'Proportion of Classified CDSs', 'Identified Loci',
                            'Proportion of Identified Loci',
                            'Valid Classifications', 'Invalid Classifications']
    
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

    # Sort loci based on decreasing sample presence
    loci_stats = sorted(loci_stats, key=lambda x: x[4], reverse=True)

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

    summary_columns = ['Total Samples', 'Total Loci', 'Total CDSs',
                       'Total CDSs Classified', 'EXC', 'INF',
                       'PLOT3', 'PLOT5', 'LOTSC', 'NIPH',
                       'NIPHEM', 'ALM', 'ASM', 'PAMA', 'LNF']

    summary_rows = [total_samples, total_loci, total_cds, loci_sums[-1]]
    summary_rows.extend(loci_sums[:-1])

    # TSV file with allelic profiles
    allelic_profiles_file = fo.join_paths(input_files, ['results_alleles.tsv'])

### Presence absence matrix
    # import matrix with allelic profiles
    matrix = pd.read_csv(allelic_profiles_file, header=0, index_col=0,
                         sep='\t', low_memory=False)

    # mask missing data
    masked_matrix = matrix.apply(im.replace_chars)

    # build presence/absence matrix
    pa_matrix, pa_outfile = Extract_cgAlleles.presAbs(masked_matrix, temp_directory)

    # Distance matrix



    # Create loci reports
    if loci_reports:
        # Get allele identifiers for samples and create FASTA files
        # Create temporary directory to store FASTA files
        fasta_dir = fo.join_paths(temp_directory, ['fasta_files'])
        fo.create_directory(fasta_dir)

        inputs = im.divide_list_into_n_chunks(loci_ids, len(loci_ids))

        common_args = [schema_directory, allelic_profiles_file, fasta_dir]

        # Add common arguments to all sublists
        inputs = im.multiprocessing_inputs(inputs,
                                           common_args,
                                           compute_locus_statistics)

        # compute loci statistics
        print('Computing loci statistics...')
        results = mo.map_async_parallelizer(inputs,
                                            mo.function_helper,
                                            cpu_cores,
                                            show_progress=True)

        # Translate FASTA files
        translation_inputs = im.divide_list_into_n_chunks(results, len(results))

        common_args = [fasta_dir, translation_table]

        # Add common arguments to all sublists
        inputs = im.multiprocessing_inputs(translation_inputs,
                                           common_args,
                                           fao.translate_fasta)

        results = mo.map_async_parallelizer(inputs,
                                            mo.function_helper,
                                            cpu_cores,
                                            show_progress=True)

        protein_files = [r[1] for r in results]

        # Align sequences with MAFFT
        alignmnet_inputs = im.divide_list_into_n_chunks(protein_files, len(protein_files))
        common_args = []

        # Add common arguments to all sublists
        inputs = im.multiprocessing_inputs(alignmnet_inputs,
                                           common_args,
                                           call_mafft)

        results = mo.map_async_parallelizer(inputs,
                                            mo.function_helper,
                                            cpu_cores,
                                            show_progress=True)
        mafft_files = [file for file in results if file is not False]

        # Determine core genome at 100%
        cgMLST_genes, cgMLST_counts = Extract_cgAlleles.compute_cgMLST(pa_matrix,
                                                                       sample_ids,
                                                                       1,
                                                                       len(sample_ids))

        # Concatenate alignment files from core
        cgMLST_genes = cgMLST_genes.tolist()
        cgMLST_alignment_files = {locus: fo.join_paths(fasta_dir, [f'{locus}_protein_aligned.fasta'])
                                  for locus in cgMLST_genes}
        sample_alignment_files = []
        # Not dealing well with '*' in allele ids
        for sample in sample_ids:
            print(sample)
            alignment = ''
            for locus, file in cgMLST_alignment_files.items():
                # Open alignmnet file
                current_file = file
                seqs = fao.import_sequences(current_file)
                try:
                    alignment += seqs[sample]
                except Exception as e:
                    print(f'Could not get {sample} allele for locus {locus}.')
            # Save cgMLST alignmnet for sample
            sample_cgMLST_alignment = fo.join_paths(fasta_dir, [f'{sample}_cgMLST_alignment.fasta'])
            align_record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [sample, alignment])
            fo.write_lines([align_record], sample_cgMLST_alignment)
            sample_alignment_files.append(sample_cgMLST_alignment)

        # Concatenate all cgMLST alignmnet records
        full_alignment = fo.join_paths(fasta_dir, ['cgMLST_alignment.fasta'])
        fo.concatenate_files(sample_alignment_files, full_alignment)
        
        # Compute NJ tree with FastTree
        out_tree = fo.join_paths(fasta_dir, ['cgMLST.tree'])
        run_fasttree(full_alignment, out_tree)

    # Create sample reports
    # if sample_reports:

    tree_data = fo.read_file(out_tree)        

    # Sort Presence-Absence matrix based on decreasing loci presence
    sorted_loci = [x[0] for x in loci_stats]
    pa_matrix = pa_matrix[sorted_loci]
    pa_lines = pa_matrix.values.tolist()

    # Compute distance matrix
    dm_output_dir = fo.join_paths(temp_directory, ['distance_matrix'])
    dm.main(allelic_profiles_file, dm_output_dir, cpu_cores, True)

    # Import distance matrix and sort
    dm_file = '/home/rmamede/Desktop/CAML_2023/presentation/AlleleCallReport_test/temp/distance_matrix/results_alleles_allelic_differences_symmetric.tsv'
    distance_m = matrix = pd.read_csv(dm_file, header=0, index_col=0,
                                      sep='\t', low_memory=False)
    # Sort based on M1UK first
    lynskey_ids = fo.read_lines('/home/rmamede/Desktop/CAML_2023/presentation/lynskey_sorted.txt')
    distance_m = distance_m[lynskey_ids]
    distance_m = distance_m.reindex(lynskey_ids)
    dm_lines = distance_m.values.tolist()

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
                   "presence_absence": pa_lines,
                   "distance_matrix": dm_lines,
                   "cgMLST_tree": tree_data,
                   }

    # Write Schema Report HTML file
    report_html = ct.ALLELECALL_REPORT_HTML
    report_html = report_html.format(json.dumps(report_data))
    report_html_file = fo.join_paths(output_directory, ['main_report.html'])
    fo.write_to_file(report_html, report_html_file, 'w', '\n')        


# Add info alert about valid (EXC+INF) and invalid (PLOT3+PLOT5+NIPH+NIPHEM+ASM+ALM)
