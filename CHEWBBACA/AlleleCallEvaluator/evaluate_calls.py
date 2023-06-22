#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module generates an interactive HTML report for a schema with
statistics about the allele size per locus and an analysis to
identify alleles that are not complete coding sequences or that
deviate from the allele mode size and minimum length. There is
also the option to include loci annotations provided in a TSV
file and to create a detailed HTML report for each locus that
can include a protein MSA and a NJ Tree.

Code documentation
------------------
"""


import json
import subprocess
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
    from ExtractCgMLST import determine_cgmlst
except ModuleNotFoundError:
    from CHEWBBACA.utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo,
        distance_matrix as dm)
    from CHEWBBACA.ExtractCgMLST import determine_cgmlst


# locus = loci_ids[0]
# schema_directory = '/home/rmamede/Desktop/test_chewbbaca320/spneumo_schema/schema_seed'
# allelic_profiles = allelic_profiles_file
# output_directory = fasta_dir
def create_locus_fasta(locus, schema_directory, allelic_profiles, output_directory):
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


input_files = '/home/rmamede/Desktop/AlleleCallReport_test/lynskey_results'
schema_directory = '/home/rmamede/Desktop/AlleleCallReport_test/spyogenes_schema_chewieNS'
output_directory = '/home/rmamede/Desktop/AlleleCallReport_test/AlleleCallReport_test'
cpu_cores = 6
annotations = '/home/rmamede/Desktop/AlleleCallReport_test/spyogenes_annotations/annotations.tsv'
loci_reports = True
translation_table = 11
def main(input_files, schema_directory, output_directory, annotations,
         cpu_cores, light, no_pa, no_dm, no_tree):

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

    # Compute data for Sample Stats table
    cds_coordinates_dir = fo.join_paths(temp_directory, ['coordinates'])
    fo.create_directory(cds_coordinates_dir)

    sample_stats = {sid: [0, 0, []] for sid in sample_ids}
    cds_coordinates_file = fo.join_paths(input_files, ['cds_coordinates.tsv'])
    with open(cds_coordinates_file, 'r') as infile:
        # Skip header line
        header = infile.__next__()
        for line in infile:
            current_line = line.strip().split('\t')
            if current_line[1] not in sample_stats[current_line[0]][2]:
                sample_stats[current_line[0]][2].append(current_line[1])
                sample_stats[current_line[0]][0] += 1
            sample_stats[current_line[0]][1] += 1

    # Get classification sum per sample
    dataset_total_cds = 0
    for i, sample in enumerate(sample_stats):
        sample_line = sample_counts.iloc[i].tolist()
        # Only counts 1 CDS per NIPH/NIPHEM
        valid = sum(sample_line[1:3])
        invalid = sum(sample_line[3:11])
        total_cds = sample_stats[sample][1]
        dataset_total_cds += sample_stats[sample][1]
        identified_loci = valid + invalid
        classified_proportion = round(identified_loci/total_cds, 3)
        identified_proportion = round(identified_loci/total_loci, 3)
        sample_stats[sample].append(float(classified_proportion))
        sample_stats[sample].append(int(identified_loci))
        sample_stats[sample].append(float(identified_proportion))
        sample_stats[sample].append(int(valid))
        sample_stats[sample].append(int(invalid))
        sample_stats[sample].extend(list(map(int, sample_line[1:])))

    sample_stats_columns = ['Sample', 'Total Contigs', 'Total CDSs',
                            'Proportion of Classified CDSs', 'Identified Loci',
                            'Proportion of Identified Loci',
                            'Valid Classifications', 'Invalid Classifications',
                            ] + sample_counts.columns.tolist()[1:]

    sample_stats_rows = []
    for k, v in sample_stats.items():
        sample_stats_rows.append([k, v[0], v[1], *v[3:]])

    # Compute data for the Loci Stats table
    loci_stats = []
    for i in range(len(loci_counts)):
        locus_line = loci_counts.iloc[i]
        locus_id = locus_line['Locus']
        total_cds = int(locus_line['Total_CDS'])
        valid = sum(locus_line[1:3])
        invalid = sum(locus_line[3:11])
        sample_proportion = round(valid/total_samples, 3)
        loci_stats.append([locus_id, int(total_cds), int(valid),
                           int(invalid), float(sample_proportion),
                           *list(map(int, locus_line[1:12]))])

    # Sort loci based on decreasing sample presence
    loci_stats = sorted(loci_stats, key=lambda x: x[4], reverse=True)

    loci_stats_columns = ['Locus', 'Total CDSs', 'Valid Classes',
                          'Invalid Classes', 'Proportion Samples',
                          ] + loci_counts.columns.tolist()[1:-1]

    # Add Loci annotations
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

    # Data for Summary Table
    # Count total number of classified CDSs per class
    loci_sums = loci_counts[loci_counts.columns[1:]].sum(axis=0)
    loci_sums = loci_sums.tolist()

    summary_columns = ['Total Samples', 'Total Loci', 'Total CDSs',
                       'Total CDSs Classified', 'EXC', 'INF',
                       'PLOT3', 'PLOT5', 'LOTSC', 'NIPH',
                       'NIPHEM', 'ALM', 'ASM', 'PAMA', 'LNF']

    summary_rows = [total_samples, total_loci, dataset_total_cds,
                    loci_sums[-1], *loci_sums[:-1]]

    # Define path to TSV file that contains allelic profiles
    allelic_profiles_file = fo.join_paths(input_files, ['results_alleles.tsv'])
    # Import matrix with allelic profiles
    profiles_matrix = pd.read_csv(allelic_profiles_file, header=0, index_col=0,
                                  sep='\t', low_memory=False)

    # Mask missing data
    masked_profiles = profiles_matrix.apply(im.replace_chars)

    # Compute Presence-Absence matrix
    pa_matrix, pa_outfile = determine_cgmlst.presAbs(masked_profiles, temp_directory)
    # Sort Presence-Absence matrix based on decreasing loci presence
    sorted_loci = [x[0] for x in loci_stats]
    pa_matrix = pa_matrix[sorted_loci]
    pa_lines = pa_matrix.values.tolist()

    # Compute distance matrix
    dm_output_dir = fo.join_paths(temp_directory, ['distance_matrix'])
    dm_file = dm.main(allelic_profiles_file, dm_output_dir, cpu_cores, True)
    # Import distance matrix
    distance_m = pd.read_csv(dm_file[0], header=0, index_col=0,
                             sep='\t', low_memory=False)
    dm_lines = distance_m.values.tolist()

    # Create FASTA files with alleles identified in samples
    # Create temporary directory to store FASTA files
    fasta_dir = fo.join_paths(temp_directory, ['fasta_files'])
    fo.create_directory(fasta_dir)

    inputs = im.divide_list_into_n_chunks(loci_ids, len(loci_ids))

    common_args = [schema_directory, allelic_profiles_file, fasta_dir]

    # Add common arguments to all sublists
    inputs = im.multiprocessing_inputs(inputs,
                                       common_args,
                                       create_locus_fasta)

    # Create FASTA files with identified alleles
    print('Creating FASTA files with identified alleles...')
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
    mafft_files = {fo.file_basename(file, False).split('_protein_aligned')[0]: file
                   for file in results if file is not False}

    # Determine core genome at 100%
    cgMLST_genes, cgMLST_counts = determine_cgmlst.compute_cgMLST(pa_matrix,
                                                                  sample_ids,
                                                                  1,
                                                                  len(sample_ids))
    cgMLST_genes = cgMLST_genes.tolist()

    ### Concatenate all alignment files and index with BioPython to get sequences faster?

    # Concatenate alignment files from core
    sample_alignment_files = []
    # Not dealing well with '*' in allele ids
    for sample in sample_ids:
        print(sample)
        alignment = ''
        for locus in cgMLST_genes:
            # Open alignmnet file
            current_file = mafft_files[locus]
            seqs = fao.import_sequences(current_file)
            try:
                alignment += seqs[sample]
            except Exception as e:
                print(f'Could not get {sample} allele for locus {locus}.')
        # Save cgMLST alignmnet for sample
        cgMLST_alignment_outfile = fo.join_paths(fasta_dir, [f'{sample}_cgMLST_alignment.fasta'])
        align_record = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE, [sample, alignment])
        fo.write_lines([align_record], cgMLST_alignment_outfile)
        sample_alignment_files.append(cgMLST_alignment_outfile)

    # Concatenate all cgMLST alignmnet records
    full_alignment = fo.join_paths(fasta_dir, ['cgMLST_alignment.fasta'])
    fo.concatenate_files(sample_alignment_files, full_alignment)

    # Compute NJ tree with FastTree
    out_tree = fo.join_paths(fasta_dir, ['cgMLST.tree'])
    run_fasttree(full_alignment, out_tree)
    tree_data = fo.read_file(out_tree)

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

    # Write report HTML file
    report_html = ct.ALLELECALL_REPORT_HTML
    report_html = report_html.format(json.dumps(report_data))
    report_html_file = fo.join_paths(output_directory, ['main_report.html'])
    fo.write_to_file(report_html, report_html_file, 'w', '\n')
