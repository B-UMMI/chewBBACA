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


import os
import math
import json
import shutil
import pickle
import statistics
from collections import Counter

try:
    from utils import (
        constants as ct,
        mafft_wrapper as mw,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (
        constants as ct,
        mafft_wrapper as mw,
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo)


def compute_quartiles(values):
    """Compute the q1, median and q3.

    Parameters
    ----------
    values : list
        List with integers and/or floats.

    Returns
    -------
    A list with the q1, median and q3 computed
    based on the input list.
    """
    median = round(statistics.median(values))
    # q1 and q3
    if len(values) > 1:
        half = int(len(values) // 2)
        q1 = statistics.median(values[:half])
        q3 = statistics.median(values[-half:])
    else:
        q1 = values[0]
        q3 = values[0]

    return [q1, median, q3]


def compute_mean(values, decimal_place=None):
    """Compute the mean.

    Parameters
    ----------
    values : list
        List with integers and/or floats.
    decimal_place : int
        Rounding precision.

    Returns
    -------
    mean : int or float
        The mean computed based on the input list.
    """
    mean = round(sum(values) / len(values), decimal_place)

    return mean


def compute_sd(values):
    """Compute the standard deviation.

    Parameters
    ----------
    values : list
        List with integers and/or floats.

    Returns
    -------
    sd : float
        The standard deviation computed based on the input list.
    """
    # standard deviation
    if len(values) > 1:
        sd = statistics.stdev(values)
    else:
        sd = 0.0

    return sd


def outside_threshold(values, threshold, limit):
    """Determine if values are outside a threshold.

    Parameters
    ----------
    values : list
        List with integers and/or floats.
    threshold : int or float
        Threshold value.
    limit : str
        Check if the values are above ('top') or below
        ('bot') the threshold.

    Returns
    -------
    outside_indices : list
        List with the indices of the values that were
        outside the threshold.
    """
    outside_indices = []
    for i, v in enumerate(values):
        if limit == 'top':
            if v > threshold:
                outside_indices.append((i, v))
        elif limit == 'bot':
            if v < threshold:
                outside_indices.append((i, v))

    return outside_indices


def count_translation_exception_categories(exceptions):
    """Count the number of ocurrences for each translation exception.

    Parameters
    ----------
    exceptions : list
        List with the strings corresponding to the exceptions captured
        during sequence translation.

    Returns
    -------
    inframe_stop : int
        Count for in-frame stop codons.
    no_start : int
        Count for sequences that do not start with a start codon.
    no_stop : int
        Count for sequences that do not end with a stop codon.
    incomplete : int
        Count for sequences with a length value that is not a
        multiple of 3.
    ambiguos : int
        Count for sequences that include ambiguous bases.
    """
    inframe_stop = sum([1 for exc in exceptions
                        if ct.TRANSLATION_EXCEPTIONS[0] in exc.split(',')[0]])
    no_start = sum([1 for exc in exceptions
                    if ct.TRANSLATION_EXCEPTIONS[1] in exc.split(',')[0]])
    no_stop = sum([1 for exc in exceptions
                   if ct.TRANSLATION_EXCEPTIONS[2] in exc.split(',')[0]])
    incomplete = sum([1 for exc in exceptions
                      if ct.TRANSLATION_EXCEPTIONS[3] in exc.split(',')[0]])
    ambiguous = sum([1 for exc in exceptions
                     if ct.TRANSLATION_EXCEPTIONS[4] in exc.split(',')[0]])

    return [inframe_stop, no_start, no_stop, incomplete, ambiguous]


def reformat_translation_exceptions(exceptions, allele_lengths):
    """Reformat translation exceptions.

    Parameters
    ----------
    exceptions : list
        List with one sublist for each exception. Each sublist
        includes an allele identifier and the exception string.
    allele_lengths : dict
        Dictionary mapping the allele identifier to the allele
        sequence length.

    Returns
    -------
    formatted_exceptions : dict
        Dictionary mapping the allele identifiers to a list
        with the allele identifier, the exception category of
        the first exception that was captured and the list of
        all exceptions captured for each allele.
    """
    formatted_exceptions = {}
    for exception in exceptions:
        allele_id = exception[0]
        sequence_length = allele_lengths[int(allele_id)]
        raw_exception = exception[1]
        if raw_exception == 'sequence length is not a multiple of 3':
            exception_category = 'Incomplete ORF'
            exception_description = ('Sequence length is not a multiple '
                                     f'of 3 ({sequence_length}bp);')
        if raw_exception == 'ambiguous or invalid characters':
            exception_category = 'Ambiguous Bases'
            exception_description = ('Sequence contains ambiguous bases;')
        if 'codon' in raw_exception:
            # split exception string
            split_exception = raw_exception.split(',')
            clean_exceptions = []
            for e in split_exception:
                estr = e.replace('(', ':')
                estr = estr.replace(')', ';')
                estr = estr.replace('.', '')
                clean_exceptions.append(estr)
            if 'Extra' in clean_exceptions[0]:
                exception_category = 'In-frame Stop Codon'
            else:
                exception_category = 'Missing Start/Stop Codon'
            exception_description = ' '.join(clean_exceptions)

        formatted_exceptions.setdefault(allele_id,
                                        [allele_id,
                                         exception_category]).append(exception_description)

    return formatted_exceptions


def get_alleleID(allele_seqid):
    """
    """
    if '_*' not in allele_seqid:
        allele_id = int(allele_seqid.split('_')[-1])
    else:
        allele_id = int(allele_seqid.split('_*')[-1])

    return allele_id


def compute_locus_statistics(locus, translation_table, minimum_length,
                             size_threshold, translation_dir):
    """Compute sequence length statistics for a locus.

    Parameters
    ----------
    locus : str
        Path to the locus FASTA file.
    translation_table : int
        Translation table to use for sequence translation.
    minimum_length : int
        Minimum sequence length value used to exclude sequences.
    size_threshold : float
        Sequence size variation threshold value used to compute bot
        and top sequence length thresholds.
    translation_dir : str
        Path to the directory where the FASTA files with the
        translated alleles will be saved to.

    Returns
    -------
    results : list
        A list with the length statistics and the sequence translation
        exceptions captured for each sequence in the FASTA file.
    """
    locus_id = fo.file_basename(locus, False)
    sequence_lengths = fao.sequence_lengths(locus)
    ns_alleles = [get_alleleID(k) for k in sequence_lengths if '*' in k]
    ns_alleles = sorted(ns_alleles)
    # sort based on sequence length
    # Determine if any allele identifiers include "*"
    # "*" is added to novel alleles in schemas downloaded from Chewie-NS
    allele_lengths = {}
    for x in sorted(sequence_lengths.items(), key=lambda item: item[1]):
        allele_id = get_alleleID(x[0])
        allele_lengths[allele_id] = x[1]

    lengths = list(allele_lengths.values())
    allele_ids = list(allele_lengths.keys())

    # Determine missing allele ids
    missing_ids = im.find_missing(allele_ids)

    # Get total number of alleles
    nr_alleles = len(lengths)

    # Summary statistics
    max_length = max(lengths)
    min_length = min(lengths)
    size_range = f'{min_length}-{max_length}'
    mean_length = compute_mean(lengths, 0)
    locus_sd = compute_sd(lengths)
    q1, median, q3 = compute_quartiles(lengths)
    mode_length = sm.determine_mode(lengths)[0]

    # Get index of alleles above size threshold
    top_threshold = math.floor(mode_length*(1+size_threshold))
    above_threshold = outside_threshold(lengths, top_threshold, 'top')
    # Get ids based on indices
    above_threshold = [(allele_ids[i[0]], i[1]) for i in above_threshold]
    # Get index of alleles below size threshold
    bot_threshold = math.ceil(mode_length*(1-size_threshold))
    below_threshold = outside_threshold(lengths, bot_threshold, 'bot')
    # Get ids based on indices
    below_threshold = [(allele_ids[i[0]], i[1]) for i in below_threshold]

    # Translate alleles and capture translation exceptions
    _, protein_file, _, exceptions = fao.translate_fasta(locus,
                                                         translation_dir,
                                                         translation_table)
    # If some sequence headers include "*", reformat FASTA file to remove "*"
    # Reformatting protein files uses less disk space than reformatting DNA
    if len(ns_alleles) > 0:
        records = fao.sequence_generator(protein_file)
        records = [[(rec.id).replace('*', ''), str(rec.seq)] for rec in records]
        output_lines = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, records)
        fo.remove_files([protein_file])
        fo.write_lines(output_lines, protein_file)

    # Determine distinct proteins
    translated_alleles = fao.import_sequences(protein_file)
    distinct_proteins = sm.determine_duplicated_seqs(translated_alleles)
    distinct_records = [[(v[0]).split('_')[-1], k]
                        for k, v in distinct_proteins.items()]
    distinct_lines = fao.fasta_lines(ct.FASTA_RECORD_TEMPLATE, distinct_records)
    distinct_file = fo.join_paths(translation_dir, [f'{locus_id}_distinct.fasta'])
    fo.write_lines(distinct_lines, distinct_file)
    distinct_ids = [[get_alleleID(i) for i in v]
                    for k, v in distinct_proteins.items()]
    distinct_ids = [[v[0], v] for v in distinct_ids]
    # Sort in order of decreasing length
    distinct_ids = sorted(distinct_ids, key=lambda x: len(x[1]), reverse=True)
    distinct_ids = sorted(distinct_ids, key=lambda x: x[0])

    exceptions = {str(get_alleleID(exc[0])): exc[1] for exc in exceptions}
    exceptions_values = list(exceptions.values())
    exceptions_lines = [[k, v] for k, v in exceptions.items()]
    # Count number of ocurrences for each translation exception
    exception_counts = count_translation_exception_categories(exceptions_values)

    # Determine list of valid allele ids
    valid_ids = [str(i) for i in allele_ids if str(i) not in exceptions]

    # Do not count shorter alleles as invalid alleles
    # Only count as invalid alleles that cannot be translated
    invalidCDS = sum(exception_counts)
    validCDS = nr_alleles - invalidCDS
    validated_proportion = round(validCDS / nr_alleles, 3)

    # determine sequences shorter than minimum length
    short_ids = [(allele_ids[i], v)
                 for i, v in enumerate(lengths)
                 if v < minimum_length]

    # format exception strings for Exceptions Table Component in loci reports
    formatted_exceptions = reformat_translation_exceptions(exceptions_lines,
                                                           allele_lengths)

    # Add info about short alleles
    short_category = f'Alleles < {minimum_length}bp'
    short_description = 'Sequence shorter than {0}bp ({1}bp);'
    for i in short_ids:
        if str(i[0]) in formatted_exceptions:
            formatted_exceptions[str(i[0])].append(short_description.format(minimum_length, i[1]))
        else:
            formatted_exceptions[str(i[0])] = [str(i[0]),
                                               short_category,
                                               short_description.format(minimum_length, i[1])]

    # Add info about alleles above threshold
    above_category = 'Alleles above size threshold'
    above_description = 'Sequence above size threshold ({0}bp)'
    for i in above_threshold:
        if str(i[0]) in formatted_exceptions:
            formatted_exceptions[str(i[0])].append(above_description.format(i[1]))
        else:
            formatted_exceptions[str(i[0])] = [str(i[0]),
                                               above_category,
                                               above_description.format(i[1])]

    # Add info about alleles below threshold
    below_category = 'Alleles below size threshold'
    below_description = 'Sequence below size threshold ({0}bp)'
    for i in below_threshold:
        if str(i[0]) in formatted_exceptions:
            formatted_exceptions[str(i[0])].append(below_description.format(i[1]))
        else:
            formatted_exceptions[str(i[0])] = [str(i[0]),
                                               below_category,
                                               below_description.format(i[1])]

    # Join exceptions per allele
    formatted_exceptions = [[v[0], v[1], ' '.join(v[2:])]
                            for k, v in formatted_exceptions.items()]

    # Unpack exception counts
    inframe_stop, no_start, no_stop, incomplete, ambiguous = exception_counts

    results = [locus_id, nr_alleles, max_length, min_length,
               size_range, median, mean_length, mode_length,
               locus_sd, q1, q3, lengths, allele_ids,
               [locus_id, nr_alleles, validCDS, invalidCDS,
                validated_proportion, len(distinct_ids), incomplete, ambiguous,
                no_start+no_stop, inframe_stop, len(short_ids),
                len(below_threshold), len(above_threshold), len(missing_ids)],
               bot_threshold, top_threshold, formatted_exceptions,
               protein_file, distinct_file, missing_ids, valid_ids,
               ns_alleles, distinct_ids]

    return results


def locus_report(locus_file, locus_data, annotation_columns,
                 annotation_values, html_dir, minimum_length,
                 light, add_sequences):
    """Create the HTML report for a locus.

    Parameters
    ----------
    locus_file : str
        Path to the locus FASTA file.
    locus_data : list
        List with the locus data returned by `compute_locus_statistics`.
    annotation_columns : list
        List with the column headers for the annotations table.
    annotation_values : list
        List with the locus annotations.
    html_dir : str
        Path to the directory where the loci HTML reports will
        be saved to.
    minimum_length : int
        Minimum sequence length.
    light : bool
        True to compute and add the MSA and NJ to the report.
        False otherwise.
    add_sequences : bool
        True to add the allele DNA sequences and the translated
        allele sequences to the report.

    Returns
    -------
    locus_html_file : str
        Path to the locus report HTML file.
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
                  locus_data[13][7],
                  locus_data[13][8],
                  locus_data[13][9],
                  locus_data[13][10],
                  locus_data[4],
                  locus_data[5],
                  locus_data[7],
                  locus_data[13][11],
                  locus_data[13][12],
                  locus_data[13][13]]

    dna_sequences = {"sequences": []}
    protein_sequences = {"sequences": []}
    # Include DNA and Protein sequences if --add-sequences was provided
    if add_sequences is True:
        protein_records = fao.sequence_generator(locus_data[17])
        for record in protein_records:
            allele_id = get_alleleID(record.id)
            protein_sequences["sequences"].append({"name": allele_id,
                                                   "sequence": str(record.seq)})
        dna_records = fao.sequence_generator(locus_file)
        for record in dna_records:
            allele_id = get_alleleID(record.id)
            dna_sequences["sequences"].append({"name": allele_id,
                                               "sequence": str(record.seq)})

    # Get data for MSA and NJ Tree if --light flag was not provided
    phylo_data = {"phylo_data": []}
    msa_data = {"sequences": []}
    if light is False:
        if locus_data[13][2] > 1 and locus_data[13][5] > 1:
            output_directory = os.path.dirname(locus_data[18])
            alignment_file = mw.call_mafft(locus_data[18], output_directory)
            # Get MSA data
            alignment_text = fo.read_file(alignment_file)
            msa_data['sequences'] = alignment_text

            # Get Tree data
            # Get the phylocanvas data
            tree_file = alignment_file.replace('_aligned.fasta', '.fasta.tree')
            phylo_data = fo.read_file(tree_file)

            # Start by substituting greatest value to avoid substituting
            # smaller values contained in greater values
            for i in range(locus_data[1], 0, -1):
                phylo_data = phylo_data.replace(f'{i}_', '')

            phylo_data = phylo_data.replace('\n', '')
            phylo_data = {"phylo_data": phylo_data}

    locus_columns = ct.LOCUS_COLUMNS.format(minimum_length,
                                            locus_data[14],
                                            locus_data[15])
    locus_columns = locus_columns.split('\t')

    # Gather data in dictionary to store in HTML
    locus_html_data = {"summaryData": [{"columns": locus_columns},
                                       {"rows": [locus_rows]}],
                       "annotations": [{"columns": annotation_columns},
                                       {"rows": [annotation_values]}],
                       "lengths": allele_lengths,
                       "ids": allele_ids,
                       "counts": [list(counts_data[0]), list(counts_data[1])],
                       "phylo": phylo_data,
                       "validIDs": locus_data[20],
                       "msa": msa_data,
                       "dna": dna_sequences,
                       "protein": protein_sequences,
                       "botThreshold": locus_data[14],
                       "topThreshold": locus_data[15],
                       "invalidAlleles": [{"columns": ct.INVALID_ALLELES_COLUMNS},
                                          {"rows": locus_data[16]}],
                       "nsAlleles": locus_data[21],
                       "distinctAlleles": [{"columns": ct.DISTINCT_ALLELES_COLUMNS},
                                           {"rows": locus_data[22]}],
                       }

    # Add data to HTML string and save HTML string into file
    locus_html = ct.LOCUS_REPORT_HTML
    locus_html = locus_html.format(json.dumps(locus_html_data))

    locus_html_file = fo.join_paths(html_dir, [f'{locus}.html'])
    fo.write_to_file(locus_html, locus_html_file, 'w', '\n')

    return locus_html_file


def main(schema_directory, output_directory, genes_list, annotations,
         translation_table, size_threshold, minimum_length,
         cpu_cores, loci_reports, light, add_sequences):

    # Create directory to store intermediate files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    schema_files = fo.read_lines(genes_list)
    # Sort based on locus name
    schema_files = sorted(schema_files)
    # Map basename to full path
    loci_basenames = {fo.file_basename(file, False): file
                      for file in schema_files}

    # Check if the schema was created with chewBBACA
    config_file = os.path.join(schema_directory, ".schema_config")
    if os.path.exists(config_file):
        # get the schema configs
        with open(config_file, "rb") as cf:
            chewie_config = pickle.load(cf)
        chewie_version = (chewie_config["chewBBACA_version"][0]).replace('chewBBACA ', '')
        print("The schema was created with chewBBACA v{0}.".format(chewie_version))
        # Message displayed in the Schema Creation Alert
        creation_message = ('Schema created with chewBBACA '
                            f'v{chewie_version}, '
                            'BLAST Score Ratio of '
                            f'{chewie_config.get("bsr")[0]}, '
                            'minimum length of '
                            f'{chewie_config.get("minimum_locus_length")[0]},'
                            ' size threshold of '
                            f'{chewie_config.get("size_threshold")[0]}, '
                            'and translation table '
                            f'{chewie_config.get("translation_table")[0]}.')
        # check if schema has a Prodigal training file
        if chewie_config.get("prodigal_training_file")[0] is not None:
            creation_message += (' A Prodigal training file was used for '
                                 'gene prediction.')
        else:
            creation_message += (' No Prodigal training file used for gene '
                                 'prediction.')
    else:
        chewie_config = {}
        # Message displayed in the Schema Creation Alert
        creation_message = ('Did not find information about config values '
                            'used to create the schema.')

    # Use schema config values if user did not provide values
    # Use default values defined in constants if there are no chewie config
    if minimum_length is None:
        minimum_length = chewie_config.get("minimum_locus_length",
                                           [ct.MINIMUM_LENGTH_DEFAULT])[0]

    if size_threshold is None:
        size_threshold = chewie_config.get("size_threshold",
                                           [ct.SIZE_THRESHOLD_DEFAULT])[0]

    if translation_table is None:
        translation_table = chewie_config.get("translation_table",
                                              [ct.GENETIC_CODES_DEFAULT])[0]

    # Message displayed in the Schema Evaluation Alert
    evaluation_message = ('Schema evaluated with minimum length of '
                          f'{minimum_length}, size threshold of '
                          f'{size_threshold}, and translation table '
                          f'{translation_table}.')

    # Create temporary directory to store translated alleles
    translation_dir = fo.join_paths(temp_directory, ['translated_loci'])
    fo.create_directory(translation_dir)

    # Calculate the summary statistics and other information about each locus.
    inputs = im.divide_list_into_n_chunks(schema_files, len(schema_files))

    common_args = [translation_table, minimum_length,
                   size_threshold, translation_dir]

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

    # group values for the same statistic
    data = list(zip(*results))

    # Read loci annotations from TSV file
    annotation_values = [[], []]
    if annotations is not None:
        annotation_lines = fo.read_tabular(annotations)
        # Add file columns
        annotation_values[0].extend(annotation_lines[0])
        # Add sublists with lines
        # Only keep lines for loci in the schema
        loci_annotations = [line for line in annotation_lines[1:]
                            if line[0] in loci_basenames]
        annotation_values[1].extend(loci_annotations)

    print(f'\nProvided annotations for {len(annotation_values[1])} '
          'loci in the schema.')

    # Columns for the Summary Data Table in the Schema Report
    summary_columns = ct.SCHEMA_SUMMARY_TABLE_HEADERS.format(minimum_length)
    summary_columns = summary_columns.split('\t')

    # Get total number of alleles with a length value not multiple of 3
    notMultiple_sum = sum([subdata[6] for subdata in data[13]])
    # Get total number of alleles with an in-frame stop codon
    stopC_sum = sum([subdata[9] for subdata in data[13]])
    # Get total number of alleles with no start or stop codon
    notStart_sum = sum([subdata[8] for subdata in data[13]])
    # Get total number of alleles shorter than the minimum length value
    shorter_sum = sum([subdata[10] for subdata in data[13]])
    # Get totla number of alleles with ambiguous bases
    ambiguous_sum = sum([subdata[7] for subdata in data[13]])
    # Get total number of alleles below or above the sequence length thresholds
    below_sum = sum([subdata[11] for subdata in data[13]])
    above_sum = sum([subdata[12] for subdata in data[13]])
    # Get total number of valid and invalid alleles
    invalid_sum = sum([subdata[3] for subdata in data[13]])
    valid_sum = sum(data[1]) - invalid_sum
    # Create list with row values for Summary Data Table in the Schema Report
    summary_rows = [len(data[0]),
                    sum(data[1]),
                    valid_sum,
                    invalid_sum,
                    notMultiple_sum,
                    ambiguous_sum,
                    notStart_sum,
                    stopC_sum,
                    shorter_sum,
                    below_sum,
                    above_sum]

    # Columns for the Allele Analysis Table component
    analysis_columns = ct.LOCI_ANALYSIS_COLUMNS.format(minimum_length)
    analysis_columns = analysis_columns.split('\t')
    analysis_rows = list(data[13])

    # Message displayed if schema was downloaded from Chewie-NS
    ns_schema = []
    ns_config = os.path.join(schema_directory, '.ns_config')
    if os.path.exists(ns_config):
        total_ns = sum([len(subdata) for subdata in data[21]])
        ns_schema.append('This schema was downloaded from Chewie-NS and '
                         f'contains {total_ns} alleles added since the last '
                         'synchronization with the remote schema.')

    # Data in the Schema Report HTML
    schema_data = {"summaryData": [{"columns": summary_columns},
                                   {"rows": [summary_rows]}],
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
                                {"rows": analysis_rows}],
                   "lociReports": 1 if loci_reports else 0,
                   "evaluationConfig": evaluation_message,
                   "creationConfig": creation_message,
                   "nsSchema": ns_schema}

    # Write Schema Report HTML file
    schema_html = ct.SCHEMA_REPORT_HTML
    schema_html = schema_html.format(json.dumps(schema_data))
    schema_html_file = fo.join_paths(output_directory, [ct.SCHEMA_REPORT_BASENAME])
    fo.write_to_file(schema_html, schema_html_file, 'w', '\n')

    # Compute data for loci reports
    if loci_reports is True:
        # create mapping between locus ID and annotations
        if annotations is not None:
            annotations_dict = {a[0]: a for a in annotation_values[1]}
        else:
            annotations_dict = {}

        # Create directory to store loci HTML reports
        html_dir = fo.join_paths(output_directory, ['loci_reports'])
        fo.create_directory(html_dir)

        # Previous data includes values per locus
        loci_data = results

        inputs = [[loci_basenames[d[0]], d,
                   annotation_values[0],
                   annotations_dict.get(d[0], [])] for d in loci_data]

        common_args = [html_dir, minimum_length, light, add_sequences]

        # add common arguments to all sublists
        inputs = im.multiprocessing_inputs(inputs, common_args, locus_report)

        # Create loci reports
        # Compute MSA and NJ tree if --light is not True
        # Add DNA and Protein sequences if --add-sequences is provided
        print('Creating loci reports...')
        loci_htmls = mo.map_async_parallelizer(inputs,
                                               mo.function_helper,
                                               cpu_cores,
                                               show_progress=True)

    # Copy the JS bundle files to the respective directories
    # JS bundle used by schema report
    script_path = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(script_path)
    # When chewBBACA is installed
    try:
        shutil.copy(fo.join_paths(script_path,
                                  ['report_template_components',
                                   'src',
                                   'bundles',
                                   'SchemaEvaluator',
                                   'schema_report',
                                   'report_bundle.js']),
                    output_directory)
    # For development
    except Exception as e:
        shutil.copy(fo.join_paths(parent_dir,
                                  ['report_template_components',
                                   'src',
                                   'bundles',
                                   'SchemaEvaluator',
                                   'schema_report',
                                   'report_bundle.js']),
                    output_directory)
    if loci_reports is True:
        # JS bundle used by loci reports
        try:
            shutil.copy(fo.join_paths(script_path,
                                      ['report_template_components',
                                       'src',
                                       'bundles',
                                       'SchemaEvaluator',
                                       'loci_reports',
                                       'report_bundle.js']),
                        html_dir)
        except Exception as e:
            shutil.copy(fo.join_paths(parent_dir,
                                      ['report_template_components',
                                       'src',
                                       'bundles',
                                       'SchemaEvaluator',
                                       'loci_reports',
                                       'report_bundle.js']),
                        html_dir)

    # Delete all temporary files
    fo.delete_directory(temp_directory)

    print(f'\nResults available in {output_directory}.')
