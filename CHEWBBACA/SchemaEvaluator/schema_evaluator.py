#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module generates an interactive report that allows the user to explore
the diversity (number of alleles) at each locus, the variation of allele
sizes per locus and the presence of alleles that are not CDSs
(when evaluating schemas called by other algorithms).

Code documentation
------------------
"""


import os
import json
import shutil
import pickle
import statistics
from collections import Counter

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


def call_mafft(genefile):
    """Call MAFFT to generate an alignment.

    Parameters
    ----------
    genefile : str
        A string with the name/path for
        the FASTA file.

    Returns
    -------
    bool
        True if sucessful, False otherwise.
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


def compute_locus_statistics(locus, translation_table, length_threshold):
    """
    """
    allele_lengths = fao.sequence_lengths(locus)
    # sort based on sequence length
    allele_lengths = {x[0]: x[1]
                      for x in sorted(allele_lengths.items(),
                                      key=lambda item: item[1])}
    lengths = list(allele_lengths.values())
    seqids = list(allele_lengths.keys())
    allele_ids = [seqid.split('_')[-1] for seqid in seqids]

    # number of alleles
    nr_alleles = len(lengths)

    # minimum and maximum values
    max_length = max(lengths)
    min_length = min(lengths)

    # Summary statistics
    median_length = round(statistics.median(lengths))
    mean_length = round(sum(lengths) / len(lengths))
    mode_length = sm.determine_mode(lengths)[0]

    # standard deviation
    if nr_alleles > 1:
        locus_sd = statistics.stdev(lengths)
    else:
        locus_sd = 0.0

    # q1 and q3
    if nr_alleles > 1:
        half = int(nr_alleles // 2)
        q1 = statistics.median(lengths[:half])
        q3 = statistics.median(lengths[-half:])
    else:
        q1 = lengths[0]
        q3 = lengths[0]

    # Conserved alleles
    # get ratio between number of alleles outside conserved threshold
    alleles_within_threshold = 0
    top_threshold = mode_length*(1+length_threshold)
    bot_threshold = mode_length*(1-length_threshold)
    for size in lengths:
        if not size > top_threshold and not size < bot_threshold:
            alleles_within_threshold += 1

    ratio = alleles_within_threshold / nr_alleles

    conserved = True if ratio >= 1 else False

    translated = [sm.translate_dna(str(record.seq), translation_table, length_threshold)
                  for record in fao.sequence_generator(locus)]

    stopC = translated.count('Extra in frame stop codon found')
    notStart = translated.count('is not a start codon')+translated.count('is not a stop codon')
    notMultiple = translated.count('sequence length is not a multiple of 3')
    shorter = translated.count('sequence shorter than')
    validCDS = nr_alleles - sum([stopC, notStart, notMultiple, shorter])

    results = [fo.file_basename(locus, False),
               nr_alleles,
               max_length,
               min_length,
               f'{min_length}-{max_length}',
               median_length,
               mean_length,
               mode_length,
               locus_sd,
               q1,
               q3,
               lengths,
               allele_ids,
               conserved,
               [fo.file_basename(locus, False),
                nr_alleles,
                validCDS,
                notMultiple,
                notStart,
                stopC,
                shorter]
               ]

    return results


def locus_report(locus_file, locus_data, translation_dir, html_dir,
                 translation_table, minimum_length, light, add_sequences):
    """
    """
    locus = locus_data[0]
    allele_lengths = locus_data[11]
    allele_ids = locus_data[12]

    # determine counts per distinct allele size
    counts = list(Counter(allele_lengths).items())
    sorted_counts = sorted(counts, key=lambda x: x[0])
    counts_data = list(zip(*sorted_counts))

    locus_rows = [locus,
                  locus_data[1],
                  locus_data[14][2],
                  locus_data[14][3],
                  locus_data[14][4],
                  locus_data[14][5],
                  locus_data[14][6],
                  locus_data[4],
                  locus_data[5],
                  locus_data[7]]

    # translate alleles
    _, protein_file, _ = fao.translate_fasta(locus_file,
                                             translation_dir,
                                             translation_table)

    dna_sequences = {"sequences": []}
    protein_sequences = {"sequences": []}
    if add_sequences is True:
        protein_records = fao.sequence_generator(protein_file)
        for record in protein_records:
            protein_sequences["sequences"].append({"name": (record.id).split('_')[-1],
                                                   "sequence": str(record.seq)})

        dna_records = fao.sequence_generator(locus_file)
        for record in dna_records:
            dna_sequences["sequences"].append({"name": (record.id).split('_')[-1],
                                               "sequence": str(record.seq)})

    phylo_data = {"phylo_data": []}
    msa_data = {"sequences": []}
    if light is False:
        if locus_data[1] > 1:
            alignment_file = call_mafft(protein_file)
            # get MSA data
            with open(alignment_file, 'r') as infile:
                alignmnet_text = infile.read()
                msa_data['sequences'] = alignmnet_text.replace(f'{locus}_', '')

            # get Tree data
            # get the phylocanvas data
            tree_file = alignment_file.replace('_aligned.fasta', '.fasta.tree')
            with open(tree_file, 'r') as phylo:
                phylo_data = phylo.read()

            for i in range(1, locus_data[1]+1):
                phylo_data = phylo_data.replace(f'{i}_{locus}_', '')
            phylo_data = phylo_data.replace('\n', '')

            phylo_data = {"phylo_data": phylo_data}

    locus_columns = ct.LOCUS_COLUMNS
    locus_columns[-4] = locus_columns[-4].format(minimum_length)
    # need to include '.' at start to work properly when referencing local files
    locus_data = {"summaryData": [{"columns": locus_columns},
                                  {"rows": [locus_rows]}],
                  "lengths": allele_lengths,
                  "ids": allele_ids,
                  "counts": [list(counts_data[0]), list(counts_data[1])],
                  "phylo": phylo_data,
                  "msa": msa_data,
                  "dna": dna_sequences,
                  "protein": protein_sequences}

    locus_html = ct.LOCUS_REPORT_HTML.format(json.dumps(locus_data))

    locus_html_file = fo.join_paths(html_dir, [f'{locus}.html'])
    fo.write_to_file(locus_html, locus_html_file, 'w', '\n')

    return locus_html_file


# schema_directory = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_schema/schema_seed'
# output_directory = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/schema_evaluation'
# genes_list = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_loci.txt'
# annotations = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_schema/schema_annotations/schema_seed_annotations.tsv'
# translation_table = 11
# threshold = 0.05
# minimum_length = 0
# conserved = False
# cpu_cores = 4
# loci_reports = True
# light = False
# add_sequences = True
def main(schema_directory, output_directory, genes_list, annotations,
         translation_table, threshold, minimum_length, conserved,
         cpu_cores, loci_reports, light, add_sequences):

    # create directory to store intermediate files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    schema_files = fo.read_lines(genes_list)
    # sort based on locus name
    schema_files = sorted(schema_files)

    # check if the schema was created with chewBBACA
    config_file = os.path.join(schema_directory, ".schema_config")
    if os.path.exists(config_file):
        # get the schema configs
        with open(config_file, "rb") as cf:
            chewie_config = pickle.load(cf)
        print("The schema was created with chewBBACA {0}.".format(
              chewie_config["chewBBACA_version"][0]))

    # Check minimum length value
    if minimum_length is None:
        minimum_length = chewie_config.get("minimum_locus_length", [0])[0]

    if threshold is None:
        threshold = chewie_config.get("size_threshold", [0])[0]

    length_threshold = minimum_length - (minimum_length*threshold)

    # Calculate the summary statistics and other information about each locus.
    inputs = im.divide_list_into_n_chunks(schema_files, len(schema_files))

    common_args = [translation_table, length_threshold]

    # add common arguments to all sublists
    inputs = im.multiprocessing_inputs(inputs, common_args, compute_locus_statistics)

    # compute statistics
    results = mo.map_async_parallelizer(inputs,
                                        mo.function_helper,
                                        cpu_cores,
                                        show_progress=True)

    # group values for the same statistic
    data = list(zip(*results))

    analysis_columns = ct.LOCI_ANALYSIS_COLUMNS
    analysis_columns[-1] = analysis_columns[-1].format(minimum_length)

    not_conserved_message = ('Locus size is considered not conserved if >1 '
                             f'allele are outside the mode +/- {threshold} '
                             'size. Loci with only 1 allele outside the '
                             'threshold are considered conserved.')

    annotation_values = []
    if annotations is not None:
        annotation_lines = fo.read_tabular(annotations)
        annotation_values.append(annotation_lines[0])
        annotation_values.append(annotation_lines[1:])

    column_data = []
    row_data = []
    # check if it is a chewBBACA schema
    if len(chewie_config) > 0:
        translation_table_config = chewie_config["translation_table"][0]
        minimum_length_config = chewie_config["minimum_locus_length"][0]

        column_data.extend(ct.SCHEMA_SUMMARY_TABLE_HEADERS_CHEWIE)
        row_data.extend([chewie_config["chewBBACA_version"][0], chewie_config["bsr"][0]])

    # build the total data dictionary
    column_data.extend(ct.SCHEMA_SUMMARY_TABLE_HEADERS)
    column_data[-1] = column_data[-1].format(minimum_length)

    notMultiple_sum = sum([l[3] for l in data[-1]])
    stopC_sum = sum([l[5] for l in data[-1]])
    notStart_sum = sum([l[4] for l in data[-1]])
    shorter_sum = sum([l[6] for l in data[-1]])
    invalid_sum = sum([notMultiple_sum, stopC_sum,
                       notStart_sum, shorter_sum])
    valid_sum = sum(data[1]) - invalid_sum
    row_data.extend([len(data[0]),
                     sum(data[1]),
                     valid_sum,
                     invalid_sum,
                     notMultiple_sum,
                     notStart_sum,
                     stopC_sum,
                     shorter_sum])

    schema_data = {"summaryData": [{"columns": column_data},
                                   {"rows": [row_data]}],
                   "loci": list(data[0]),
                   "total_alleles": list(data[1]),
                   "max": list(data[2]),
                   "min": list(data[3]),
                   "median": list(data[5]),
                   "mode": list(data[7]),
                   "q1": list(data[9]),
                   "q3": list(data[10]),
                   "annotations": [{"columns": annotation_values[0]},
                                   {"rows": annotation_values[1]}],
                   "analysis": [{"columns": analysis_columns},
                                {"rows": list(data[-1])}],
                   "lociReports": 1 if loci_reports else 0}

    # Write HTML file
    schema_html = ct.SCHEMA_REPORT_HTML.format(json.dumps(schema_data))

    schema_html_file = fo.join_paths(output_directory, ['schema_report.html'])
    fo.write_to_file(schema_html, schema_html_file, 'w', '\n')

    if loci_reports is True:
        translation_dir = fo.join_paths(temp_directory, ['translated_loci'])
        fo.create_directory(translation_dir)
        html_dir = fo.join_paths(output_directory, ['loci_reports'])
        fo.create_directory(html_dir)
        schema_files = {fo.file_basename(file, False): file for file in schema_files}
        loci_data = results

        inputs = [[schema_files[d[0]], d] for d in loci_data]

        common_args = [translation_dir, html_dir, translation_table,
                       minimum_length, light, add_sequences]

        # add common arguments to all sublists
        inputs = im.multiprocessing_inputs(inputs, common_args, locus_report)

        # compute statistics
        loci_htmls = mo.map_async_parallelizer(inputs,
                                               mo.function_helper,
                                               cpu_cores,
                                               show_progress=True)

    # Copy the JS bundle files to the respective directories
    # JS bundle used by schema report
    script_path = os.path.dirname(os.path.abspath(__file__))
    shutil.copy(os.path.join(script_path, "SchemaEvaluator",
                             "resources", "schema_bundle.js"),
                output_directory)
    if loci_reports is True:
        # JS bundle used by loci reports
        shutil.copy(os.path.join(script_path, "SchemaEvaluator",
                                 "resources", "loci_bundle.js"),
                    html_dir)

    fo.delete_directory(temp_directory)
