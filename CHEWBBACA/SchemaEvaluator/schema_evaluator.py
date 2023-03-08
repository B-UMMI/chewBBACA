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
import math
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


def compute_locus_statistics(locus, translation_table, minimum_length, size_threshold):
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
    alleles_below_threshold = 0
    alleles_above_threshold = 0
    top_threshold = math.floor(mode_length*(1+size_threshold))
    bot_threshold = math.ceil(mode_length*(1-size_threshold))
    for size in lengths:
        if size < bot_threshold:
            alleles_below_threshold += 1
        elif size > top_threshold:
            alleles_above_threshold += 1

    translated = {(record.id).split('_')[-1]: sm.translate_dna(str(record.seq), translation_table, minimum_length)
                  for record in fao.sequence_generator(locus)}

    exceptions = {k: v for k, v in translated.items() if type(v) == str}

    exceptions_values = list(exceptions.values())
    exceptions_lines = [[k, v] for k, v in exceptions.items()]

    stopC = exceptions_values.count('Extra in frame stop codon found')
    notStart = exceptions_values.count('is not a start codon')+exceptions_values.count('is not a stop codon')
    notMultiple = exceptions_values.count('sequence length is not a multiple of 3')
    shorter = exceptions_values.count('sequence shorter than')
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
               [fo.file_basename(locus, False),
                nr_alleles,
                validCDS,
                notMultiple,
                notStart,
                stopC,
                shorter,
                alleles_below_threshold,
                alleles_above_threshold],
               bot_threshold,
               top_threshold,
               exceptions_lines
               ]

    return results


def locus_report(locus_file, locus_data, annotation_columns,
                 annotation_values, translation_dir, html_dir,
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
                  locus_data[13][2],
                  locus_data[13][3],
                  locus_data[13][4],
                  locus_data[13][5],
                  locus_data[13][6],
                  locus_data[4],
                  locus_data[5],
                  locus_data[7],
                  locus_data[13][7],
                  locus_data[13][8]]

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

    locus_columns = ct.LOCUS_COLUMNS.format(minimum_length,
                                            locus_data[14],
                                            locus_data[15])
    locus_columns = locus_columns.split('\t')
    # need to include '.' at start to work properly when referencing local files
    locus_html_data = {"summaryData": [{"columns": locus_columns},
                                       {"rows": [locus_rows]}],
                       "annotations": [{"columns": annotation_columns},
                                       {"rows": [annotation_values]}],
                       "lengths": allele_lengths,
                       "ids": allele_ids,
                       "counts": [list(counts_data[0]), list(counts_data[1])],
                       "phylo": phylo_data,
                       "msa": msa_data,
                       "dna": dna_sequences,
                       "protein": protein_sequences,
                       "botThreshold": locus_data[14],
                       "topThreshold": locus_data[15],
                       "invalidAlleles": locus_data[16]}

    locus_html = ct.LOCUS_REPORT_HTML
    locus_html = locus_html.format(json.dumps(locus_html_data))

    locus_html_file = fo.join_paths(html_dir, [f'{locus}.html'])
    fo.write_to_file(locus_html, locus_html_file, 'w', '\n')

    return locus_html_file


# schema_directory = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_schema/schema_seed'
# output_directory = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/schema_evaluation'
# genes_list = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_loci.txt'
# annotations = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_schema/schema_annotations/schema_seed_annotations.tsv'
# translation_table = 11
# size_threshold = 0.2
# minimum_length = 0
# cpu_cores = 4
# loci_reports = True
# light = False
# add_sequences = True
def main(schema_directory, output_directory, genes_list, annotations,
         translation_table, size_threshold, minimum_length,
         cpu_cores, loci_reports, light, add_sequences):

    # create directory to store intermediate files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    schema_files = fo.read_lines(genes_list)
    # sort based on locus name
    schema_files = sorted(schema_files)

    # check if the schema was created with chewBBACA
    config_file = os.path.join(schema_directory, ".schema_config")
    chewie_config = {}
    creation_message = ('Did not find information about config values used '
                        'to create the schema.')
    if os.path.exists(config_file):
        # get the schema configs
        with open(config_file, "rb") as cf:
            chewie_config = pickle.load(cf)
        print("The schema was created with chewBBACA {0}.".format(
              chewie_config["chewBBACA_version"][0]))
        creation_message = ('Schema created with chewBBACA '
                            f'v{chewie_config.get("chewBBACA_version")[0]}, '
                            'BLAST Score Ratio of '
                            f'{chewie_config.get("bsr")[0]}, '
                            'minimum length of '
                            f'{chewie_config.get("minimum_locus_length")[0]},'
                            ' size threshold of '
                            f'{chewie_config.get("size_threshold")[0]}, '
                            'and a translation table of '
                            f'{chewie_config.get("translation_table")[0]}.')
        if chewie_config.get("prodigal_training_file")[0] is not None:
            creation_message += (' A Prodigal training file was used for '
                                 'gene prediction.')
        else:
            creation_message += (' No Prodigal training file used for gene '
                                 'prediction.')

    if minimum_length is None:
        minimum_length = chewie_config.get("minimum_locus_length",
                                           [ct.MINIMUM_LENGTH_DEFAULT])[0]

    if size_threshold is None:
        size_threshold = chewie_config.get("size_threshold",
                                           [ct.SIZE_THRESHOLD_DEFAULT])[0]

    if translation_table is None:
        translation_table = chewie_config.get("translation_table",
                                              [ct.GENETIC_CODES_DEFAULT])[0]

    evaluation_message = ('Schema evaluated with minimum length of '
                          f'{minimum_length}, size threshold of '
                          f'{size_threshold}, and translation table of '
                          f'{translation_table}.')

    # Calculate the summary statistics and other information about each locus.
    inputs = im.divide_list_into_n_chunks(schema_files, len(schema_files))

    common_args = [translation_table, minimum_length, size_threshold]

    # add common arguments to all sublists
    inputs = im.multiprocessing_inputs(inputs, common_args, compute_locus_statistics)

    # compute statistics
    results = mo.map_async_parallelizer(inputs,
                                        mo.function_helper,
                                        cpu_cores,
                                        show_progress=True)

    # group values for the same statistic
    data = list(zip(*results))

    analysis_columns = ct.LOCI_ANALYSIS_COLUMNS.format(minimum_length)
    analysis_columns = analysis_columns.split('\t')

    annotation_values = [[], []]
    if annotations is not None:
        annotation_lines = fo.read_tabular(annotations)
        annotation_values[0].extend(annotation_lines[0])
        annotation_values[1].extend(annotation_lines[1:])

    # build the total data dictionary
    column_data = ct.SCHEMA_SUMMARY_TABLE_HEADERS.format(minimum_length)
    column_data = column_data.split('\t')

    notMultiple_sum = sum([l[3] for l in data[13]])
    stopC_sum = sum([l[5] for l in data[13]])
    notStart_sum = sum([l[4] for l in data[13]])
    shorter_sum = sum([l[6] for l in data[13]])
    below_sum = sum([l[7] for l in data[13]])
    above_sum = sum([l[8] for l in data[13]])
    invalid_sum = sum([notMultiple_sum, stopC_sum,
                       notStart_sum, shorter_sum])
    valid_sum = sum(data[1]) - invalid_sum
    row_data = [len(data[0]),
                sum(data[1]),
                valid_sum,
                invalid_sum,
                notMultiple_sum,
                notStart_sum,
                stopC_sum,
                shorter_sum,
                below_sum,
                above_sum]

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
                                {"rows": list(data[13])}],
                   "lociReports": 1 if loci_reports else 0,
                   "evaluationConfig": evaluation_message,
                   "creationConfig": creation_message}

    # Write HTML file
    schema_html = ct.SCHEMA_REPORT_HTML
    schema_html = schema_html.format(json.dumps(schema_data))

    schema_html_file = fo.join_paths(output_directory, ['schema_report.html'])
    fo.write_to_file(schema_html, schema_html_file, 'w', '\n')

    if loci_reports is True:
        # create mapping between locus ID and annotations
        if annotations is not None:
            annotations_dict = {a[0]: a for a in annotation_values[1]}
        else:
            annotations_dict = {}

        translation_dir = fo.join_paths(temp_directory, ['translated_loci'])
        fo.create_directory(translation_dir)
        html_dir = fo.join_paths(output_directory, ['loci_reports'])
        fo.create_directory(html_dir)
        schema_files = {fo.file_basename(file, False): file for file in schema_files}
        loci_data = results

        inputs = [[schema_files[d[0]], d,
                   annotation_values[0],
                   annotations_dict.get(d[0], [])] for d in loci_data]

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
