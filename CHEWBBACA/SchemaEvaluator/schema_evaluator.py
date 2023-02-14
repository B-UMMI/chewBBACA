#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module generates an interactive report that allows the user to explore
the diversity (number of alleles) at each locus, the variation of allele
sizes per locus and the presence of alleles that are not CDSs
(when evaluating schemas called by other algorithms).

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_files`` : Path to the folder containing the fasta files,
  one fasta file per gene/locus (alternatively, a file with a list of paths
  can be given).

    - e.g.: ``/home/user/schemas/schema_dir``

- ``-l``, ``output_directory`` : The directory where the output files will
  be saved (will create the directory if it does not exist).

    - e.g.: ``/home/user/schemaReport``

- ``-ta``, ``translation_table`` : Genetic code to use for CDS
  translation (default=11, for Bacteria and Archaea).

    - e.g.: ``11``

- ``--cpu``, ``cpu_cores`` : The number of CPU cores to use (default=1).

    - e.g.: ``4``

- ``--light``, ``light_mode`` : Skip clustal and mafft (default=False)

Code documentation
------------------
"""


import os
import sys
import json
import pickle
import statistics

from Bio.Align.Applications import MafftCommandline

try:
    from utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        multiprocessing_operations as mo)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (
        constants as ct,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        multiprocessing_operations as mo)


schema_directory = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_schema/schema_seed'
output_directory = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/schema_evaluation'
annotations = '/home/rmamede/Desktop/Brucella_Mostafa/chewbbaca3/test_schema/schema_annotations/schema_seed_annotations.tsv'
translation_table = 11
threshold = 0.05
minimum_length = 0
conserved = False
cpu_cores = 4
light = False
def main(schema_directory, output_directory, annotations, translation_table,
         threshold, minimum_length, conserved, cpu_cores, light):
    """Create a file with pre-computed data for the Plotly charts.

    Parameters
    ----------
    schema_dir : list
        A list with names/paths for FASTA files.
    translation_table: int
        The translation table to be used.
    output_path : str
        The directory where the output files will
        be saved.
    annotations : str
        Path to the output file of the UniprotFinder
        module
    cpu_to_use: int
        Number of CPU cores to use for multiprocessing.
    minimum_length: int
        Minimum sequence length accepted in nt.
    size_threshold: int
        CDS size variation threshold.
    chewie_schema: bool
        Identifies the schema as a chewBBACA created schema.
    show_progress: bool
        Shows a progress bar for multiprocessing.

    Returns
    -------
    out_path : str
        The directory where the output files will
        be saved.
    """
    # create directory to store intermediate files
    fo.create_directory(output_directory)
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    schema_files = fo.listdir_fullpath(schema_directory)
    schema_files = fo.filter_files(schema_files, ct.FASTA_SUFFIXES)
    # only keep files whose content is typical of a FASTA file
    schema_files = fao.filter_non_fasta(schema_files)
    # sort based on locus name
    schema_files = sorted(schema_files)

    if len(schema_files) == 0:
        sys.exit('The schema directory is empty. Please check your '
                 'path. Exiting...')

    # get chewBBACA schema config values
    config_file = os.path.join(schema_directory, '.schema_config')
    if os.path.isfile(config_file):
        with open(config_file, "rb") as cf:
            chewie_config = pickle.load(cf)
    else:
        chewie_config = {}

    # Check minimum length value
    if minimum_length is None:
        minimum_length = chewie_config.get("minimum_locus_length", [0])[0]

    if threshold is None:
        threshold = chewie_config.get("size_threshold", [0])[0]

    length_threshold = minimum_length - (minimum_length*threshold)

    # Calculate the summary statistics and other information about each locus.
    # can add multiprocessing
    loci = []
    total_alleles = []
    max_values = []
    min_values = []
    size_ranges = []
    median_values = []
    mean_values = []
    mode_values = []
    sd_values = []
    q1_values = []
    q3_values = []
    all_lengths = []
    all_ids = []
    conserved_loci = []
    analysis_values = []
    notMultiple_values = []
    stopC_values = []
    notStart_values = []
    shorter_values = []
    valid_cds_values = []
    for locus in schema_files:
        allele_lengths = fao.sequence_lengths(locus)
        # sort based on sequence length
        allele_lengths = {x[0]: x[1] 
                          for x in sorted(allele_lengths.items(),
                                          key = lambda item: item[1])}
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
        top_threshold = mode_length*(1+threshold)
        bot_threshold = mode_length*(1-threshold)
        for size in lengths:
            if not size > top_threshold and not size < bot_threshold:
                alleles_within_threshold += 1

        ratio = alleles_within_threshold / nr_alleles

        loci.append(fo.file_basename(locus, False))
        total_alleles.append(nr_alleles)
        max_values.append(max_length)
        min_values.append(min_length)
        size_ranges.append(f'{min_length}-{max_length}')
        median_values.append(median_length)
        mean_values.append(mean_length)
        mode_values.append(mode_length)
        sd_values.append(locus_sd)
        q1_values.append(q1)
        q3_values.append(q3)
        all_lengths.append(lengths)
        all_ids.append(allele_ids)

        if ratio >= 1:
            conserved_loci.append(True)
        else:
            conserved_loci.append(False)

        translated = [sm.translate_dna(str(record.seq), translation_table, length_threshold)
                      for record in fao.sequence_generator(locus)]

        stopC = translated.count('Extra in frame stop codon found')
        notStart = translated.count('is not a start codon')+translated.count('is not a stop codon')
        notMultiple = translated.count('sequence length is not a multiple of 3')
        shorter = translated.count('sequence shorter than')
        notMultiple_values.append(notMultiple)
        stopC_values.append(stopC)
        notStart_values.append(notStart)
        shorter_values.append(shorter)
        validCDS = nr_alleles - sum([stopC, notStart, notMultiple, shorter])
        valid_cds_values.append(validCDS)

        analysis_values.append([fo.file_basename(locus, False),
                                notMultiple, notStart,
                                stopC, shorter, validCDS])

    analysis_columns = ["Locus", "not Multiple 3", "no Start Codon",
                        "inner Stop Codon", "Shorter than", "Valid CDSs"]
    analysis_data = [analysis_columns, analysis_values]

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

        column_data.extend(["chewBBACA version", "BLAST Score Ratio"])
        row_data.extend([chewie_config["chewBBACA_version"][0], chewie_config["bsr"][0]])

    # build the total data dictionary
    column_data.extend(["Total loci", "Total alleles", "Total alleles mult3",
                        "Total alleles stopC", "Total alleles notStart",
                        "Total alleles shorter", "Total invalid alleles"])
    
    notMultiple_sum = sum(notMultiple_values)
    stopC_sum = sum(stopC_values)
    notStart_sum = sum(notStart_values)
    shorter_sum = sum(shorter_values)
    row_data.extend([len(loci), sum(total_alleles),
                     notMultiple_sum, stopC_sum,
                     notStart_sum, shorter_sum,
                     sum([notMultiple_sum, stopC_sum,
                          notStart_sum, shorter_sum])])

    schema_data = {"summaryData": [{"columns": column_data},
                                   {"rows": [row_data]}],
                   "loci": loci,
                   "total_alleles": total_alleles,
                   "max": max_values,
                   "min": min_values,
                   "median": median_values,
                   "mode": mode_values,
                   "q1": q1_values,
                   "q3": q3_values,
                   "annotations": [{"columns": annotation_values[0]},
                                   {"rows": annotation_values[1]}],
                   "analysis": [{"columns": analysis_data[0]},
                                {"rows": analysis_data[1]}]}

    # Write HTML file
    schema_html = """<!DOCTYPE html>
    <html lang="en">
        <head>
            <meta charset="UTF-8" />
            <meta name="viewport" content="width=device-width, initial-scale=1.0" />
            <title>Schema Evaluator - React Edition</title>
        </head>
        <body style="background-color: #f6f6f6">
            <noscript> You need to enable JavaScript to run this app. </noscript>
            <div id="root"></div>
            <script> preComputedData = {0} </script>
            <script src="./main.js"></script>
        </body>
    </html>""".format(json.dumps(schema_data))

    schema_html_file = fo.join_paths(output_directory, ['schema_report.html'])
    fo.write_to_file(schema_html, schema_html_file, 'w', '\n')

    # data for loci reports
    # create directory to store translated schema
    translation_dir = fo.join_paths(temp_directory, ['translated_loci'])
    fo.create_directory(translation_dir)
    loci_htmls_dir = fo.join_paths(output_directory, ['loci_reports'])
    fo.create_directory(loci_htmls_dir)
    loci_htmls = []
    locus_columns = ["Locus", "Number of Alleles", "Alleles not multiple of 3",
                     "Alleles w/ >1 stop codon", "Alleles wo/ Start/Stop codon",
                     "Alleles shorter than 201", "Size Range (bp)",
                     "Alleles Median", "Alleles Mode"]
    for i, locus in enumerate(loci):
        allele_lengths = all_lengths[i]
        allele_ids = all_ids[i]

        locus_rows = [locus, total_alleles[i], notMultiple_values[i],
                      stopC_values[i], notStart_values[i], shorter_values[i],
                      size_ranges[i], median_values[i], mode_values[i]]

        if light is False:
            if total_alleles[i] > 1:
                # translate alleles
                _, protein_file, _ = fao.translate_fasta(schema_files[i],
                                                         translation_dir,
                                                         translation_table)
                alignment_file = call_mafft(protein_file)
                # get MSA data
                msa_data = {"sequences": []}
                msa_records = fao.sequence_generator(alignment_file)
                for record in msa_records:
                    msa_data["sequences"].append({"name": (record.id).split('_')[-1], "sequence": str(record.seq)})

                # get Tree data
                # get the phylocanvas data
                tree_file = alignment_file.replace('_aligned.fasta', '.fasta.tree')
                with open(tree_file, 'r') as phylo:
                    phylo_data = phylo.read()
                    
                phylo_data = {"phylo_data": phylo_data}

                # get sequences for Sequence Logo
                with open(protein_file, 'r') as infile:
                    protein_sequences = infile.read()
            else:
                msa_data = "undefined"
                phylo_data = "undefined"
                protein_sequences = "undefined"

        locus_data = {"summaryData": [{"columns": locus_columns},
                                      {"rows": [locus_rows]}],
                       "lengths": allele_lengths,
                       "ids": allele_ids,
                       "msa": msa_data,
                       "phylo": phylo_data,
                       "logo": protein_sequences}

        locus_html = """<!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="UTF-8" />
                <meta name="viewport" content="width=device-width, initial-scale=1.0" />
                <title>Schema Evaluator - Individual Analysis</title>
            </head>
            <body style="background-color: #f6f6f6">
                <noscript> You need to enable JavaScript to run this app. </noscript>
                <div id="root"></div>
                <script> preComputedDataInd = {0} </script>
                <script src="./main_ind.js"></script>
            </body>
        </html>""".format(json.dumps(locus_data))

        locus_html_file = fo.join_paths(loci_htmls_dir, [f'{locus}.html'])
        fo.write_to_file(locus_html, locus_html_file, 'w', '\n')
        loci_htmls.append(locus_html_file)

    return [schema_html, loci_htmls]


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
        mafft_cline = MafftCommandline(
            input=genefile,
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
