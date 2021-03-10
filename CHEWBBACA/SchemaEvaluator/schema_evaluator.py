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
import csv
import sys
import json
import pickle
import itertools
import statistics
import multiprocessing
from operator import itemgetter
from collections import Counter


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline

try:
    from utils import (
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo,
    )
except:
    from CHEWBBACA.utils import (
        file_operations as fo,
        fasta_operations as fao,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo,
    )


# Schema Evaluator Auxiliary Functions
def gene_seqs_info_schema_evaluator(gene, threshold, conserved):
    """Determines the total number of alleles and the mean length
    of allele sequences per gene.

    Parameters
    ----------
    gene : string
        a string with names/paths for FASTA files.

    Returns
    -------
    genes_info : list
        a list with a sublist for each input
        gene file. Each sublist contains a gene identifier, the
        total number of alleles for that gene and the mean length
        of allele sequences for that gene.
    """

    seq_generator = SeqIO.parse(gene, "fasta")
    alleles_lengths = [len(allele) for allele in seq_generator]
    alleles_lengths.sort()

    if len(alleles_lengths) == 0:
        sys.exit("At least one file is empty or it doesn't exist. Exiting...")

    # number of alleles
    nr_alleles = len(alleles_lengths)

    # minimum and maximum values
    max_length = max(alleles_lengths)
    min_length = min(alleles_lengths)

    # Summary statistics
    median_length = round(statistics.median(alleles_lengths))
    mean_length = round(sum(alleles_lengths) / len(alleles_lengths))
    mode_length = Counter(alleles_lengths).most_common()[0][0]

    # Conserved alleles
    # get ratio between number of alleles outside conserved threshold
    alleles_within_threshold = 0
    for size in alleles_lengths:
        if not float(size) > mode_length * (1 + threshold) and not float(
            size
        ) < mode_length * (1 - threshold):
            alleles_within_threshold += 1

    ratio = alleles_within_threshold / float(nr_alleles)

    if not conserved and (
        ratio >= 1 or len(alleles_lengths) - alleles_within_threshold < 2
    ):
        if nr_alleles == 1:
            genes_info = {
                "gene": gene,
                "nr_alleles": nr_alleles,
                "min_length": min_length,
                "max_length": max_length,
                "median_length": median_length,
                "mean_length": mean_length,
                "mode_length": mode_length,
                "alleles_lengths": alleles_lengths,
                "conserved": True,
                "one_allele_only": True,
            }
        else:
            genes_info = {
                "gene": gene,
                "nr_alleles": nr_alleles,
                "min_length": min_length,
                "max_length": max_length,
                "median_length": median_length,
                "mean_length": mean_length,
                "mode_length": mode_length,
                "alleles_lengths": alleles_lengths,
                "conserved": True,
                "one_allele_only": False,
            }
    elif conserved and ratio >= 1:
        if nr_alleles == 1:
            genes_info = {
                "gene": gene,
                "nr_alleles": nr_alleles,
                "min_length": min_length,
                "max_length": max_length,
                "median_length": median_length,
                "mean_length": mean_length,
                "mode_length": mode_length,
                "alleles_lengths": alleles_lengths,
                "conserved": True,
                "one_allele_only": True,
            }
        else:
            genes_info = {
                "gene": gene,
                "nr_alleles": nr_alleles,
                "min_length": min_length,
                "max_length": max_length,
                "median_length": median_length,
                "mean_length": mean_length,
                "mode_length": mode_length,
                "alleles_lengths": alleles_lengths,
                "conserved": True,
                "one_allele_only": False,
            }
    else:
        genes_info = {
            "gene": gene,
            "nr_alleles": nr_alleles,
            "min_length": min_length,
            "max_length": max_length,
            "median_length": median_length,
            "mean_length": mean_length,
            "mode_length": mode_length,
            "alleles_lengths": alleles_lengths,
            "conserved": False,
            "one_allele_only": False,
        }

    return genes_info


def gene_seqs_info_individual_schema_evaluator(gene):
    """Determines the total number of alleles and
    the mean length of allele sequences for each locus.

    Parameters
    ----------
    gene : string
        a string with names/paths for FASTA files.

    Returns
    -------
    genes_info : list
        a list with a sublist for each input
        gene file. Each sublist contains a gene identifier, the
        total number of alleles for that gene and the mean length
        of allele sequences for that gene.
    """

    seq_generator = SeqIO.parse(gene, "fasta")

    # Allele IDs and lengths
    allele_ids = []
    alleles_lengths = []
    for allele in seq_generator:
        allele_ids.append(allele.id)
        alleles_lengths.append(len(allele))

    alleles_lengths.sort()

    # number of alleles
    nr_alleles = len(alleles_lengths)

    # minimum and maximum values
    max_length = max(alleles_lengths)
    min_length = min(alleles_lengths)

    # size range
    size_range = "{0}-{1}".format(min_length, max_length)

    # Summary statistics
    median_length = round(statistics.median(alleles_lengths))
    mode_length = Counter(alleles_lengths).most_common()[0][0]

    genes_info = {
        "gene": gene,
        "nr_alleles": nr_alleles,
        "size_range": size_range,
        "median_length": median_length,
        "allele_ids": allele_ids,
        "alleles_lengths": alleles_lengths,
        "mode_length": mode_length,
    }

    return genes_info


def gene_seqs_info_boxplot(schema_dir):
    """Determines boxplot statistics.

    Parameters
    -----------
    schema_dir : list
        a list with names/paths for FASTA files.

    Returns
    -------
    json_to_file : dict
        a dict with a subdict for each input
        gene file. Each subdict contains information about each
        gene such as, mode of the allele sizes, number of alleles
        and summary statistics (min, max, median, mode and quartiles).
    """

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    json_to_file = {"boxplot_data": []}

    for g in schema_files:
        seq_generator = SeqIO.parse(g, "fasta")
        alleles_lengths = [len(allele) for allele in seq_generator]
        alleles_lengths.sort()

        # locus name
        locus_name = os.path.split(g)[1]

        # number of alleles
        nr_alleles = len(alleles_lengths)

        # minimum and maximum values
        loci_max = max(alleles_lengths)
        loci_min = min(alleles_lengths)

        # standard deviation
        if nr_alleles > 1:
            locus_sd = statistics.stdev(alleles_lengths)
        else:
            locus_sd = 0.0

        # median
        median_length = round(statistics.median(alleles_lengths))

        # mean
        mean_length = round(sum(alleles_lengths) / nr_alleles)

        # q1 and q3
        if nr_alleles > 1:
            half = int(nr_alleles // 2)
            q1 = statistics.median(alleles_lengths[:half])
            q3 = statistics.median(alleles_lengths[-half:])
        else:
            q1 = alleles_lengths[0]
            q3 = alleles_lengths[0]

        json_to_file["boxplot_data"].append(
            {
                "locus_name": locus_name,
                "nr_alleles": nr_alleles,
                "max": loci_max,
                "min": loci_min,
                "sd": locus_sd,
                "median": median_length,
                "mean": mean_length,
                "q1": q1,
                "q3": q3,
            }
        )

    # Sort data by locus_name
    for k in json_to_file:
        json_to_file[k] = sorted(json_to_file[k], key=itemgetter("locus_name"))

    return json_to_file


# Functions that obtain the data for panel E
def create_cds_df(
    schema_file, minimum_length, minimum_length_to_translate, translation_table
):
    """Detects alleles that aren't CDSs.

    Parameters
    ----------
    schema_dir : list
        a list with names/paths for FASTA files.
    minimum_length: int
        Minimum sequence length accepted in nt.
    translation_table: int
        the translation table to be used.

    Returns
    -------
    res_sorted : dict
        a dict obtained from a dataframe containing the
        number of non-CDSs detected. It will be used to
        populate a table in the report.
    """
    res = {"stats": []}

    gene_res = {"Gene": os.path.split(schema_file)[1]}

    gene_res["Number of alleles"] = fao.count_sequences(schema_file)

    stopC = 0
    notStart = 0
    notMultiple = 0
    shorter = 0
    CDS = 0
    allele_ids = []

    for allele in SeqIO.parse(schema_file, "fasta"):

        # FASTA headers examples: >allele_1 or >1_2
        if "_" in allele.id:
            allele_ids.append(int(allele.id.split("_")[-1]))
        # FASTA header example: >1
        else:
            allele_ids.append(int(allele.id))

        ola = sm.translate_dna(
            str(allele.seq), translation_table, minimum_length_to_translate
        )

        if "sequence length is not a multiple of 3" in ola:
            notMultiple += 1
        elif "Extra in frame stop codon found" in ola:
            stopC += 1
        elif "is not a start codon" in ola:
            notStart += 1
        elif "is not a stop codon" in ola:
            notStart += 1
        elif "sequence shorter than" in ola:
            shorter += 1
        else:
            CDS += 1

    if len(im.find_missing(allele_ids)) > 0:
        missing_allele_ids = im.find_missing(allele_ids)
    else:
        missing_allele_ids = ["None"]

    gene_res["Alleles not multiple of 3"] = notMultiple
    gene_res["Alleles w/ >1 stop codons"] = stopC
    gene_res["Alleles wo/ Start/Stop Codon"] = notStart
    gene_res["Alleles shorter than {0} nucleotides".format(minimum_length)] = shorter
    gene_res["Total Invalid Alleles"] = notMultiple + stopC + notStart + shorter
    gene_res["Missing Allele IDs"] = missing_allele_ids
    gene_res["CDS"] = CDS

    res["stats"].append(gene_res)

    res_sorted = sorted(res["stats"], key=itemgetter("Number of alleles"))

    return res_sorted


def create_pre_computed_data(
    schema_dir,
    translation_table,
    output_path,
    annotations,
    cpu_to_use,
    minimum_length,
    size_threshold,
    threshold,
    conserved,
    chewie_schema=False,
    show_progress=False,
):
    """Creates a file with pre-computed data for
    the Schema Evaluator plotly charts.

    Parameters
    ----------
    schema_dir : list
        a list with names/paths for FASTA files.
    translation_table: int
        the translation table to be used.
    output_path : str
        the directory where the output files will
        be saved.
    annotations : str
        path to the output file of the UniprotFinder
        module
    cpu_to_use: int
        number of CPU cores to use for multiprocessing.
    minimum_length: int
        minimum sequence length accepted in nt.
    size_threshold: int
        CDS size variation threshold.
    chewie_schema: bool
        identifies the schema as a chewBBACA created schema.
    show_progress: bool
        shows a progress bar for multiprocessing.

    Returns
    -------
    out_path : str
        the directory where the output files will
        be saved.

    """

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    if len(schema_files) < 1:
        sys.exit("The schema directory is empty. Please check your path. Exiting...")

    # Check if files are empty
    empty_files_1 = [fo.is_file_empty(f) for f in schema_files]
    empty_files_2 = [fo.is_file_empty_3(f) for f in schema_files]

    if True in empty_files_1:
        sys.exit("At least one file is empty or it doesn't exist. Exiting...")
    elif True in empty_files_2:
        sys.exit("At least one file is empty or it doesn't exist. Exiting...")

    out_path = os.path.join(output_path, "SchemaEvaluator_pre_computed_data")
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # Check minimum length value
    if minimum_length is None:

        minimum_length = 0
        minimum_length_to_translate = minimum_length - (minimum_length * size_threshold)

        if chewie_schema:
            # read config file to get chewBBACA parameters
            config_file = os.path.join(schema_dir, ".schema_config")
            with open(config_file, "rb") as cf:
                chewie_schema_configs = pickle.load(cf)

            minimum_length = chewie_schema_configs["minimum_locus_length"][0]
            minimum_length_to_translate = minimum_length - (
                minimum_length * chewie_schema_configs["size_threshold"][0]
            )

    minimum_length_to_translate = minimum_length - (
        minimum_length * size_threshold
    )  # set the minimum length value for translation

    if not os.listdir(out_path):
        # Calculate the summary statistics and other information about each locus.
        print("\nCalculating summary statistics...\n")

        results = []

        pool_main = multiprocessing.Pool(processes=cpu_to_use)

        rawr_main = pool_main.starmap_async(
            gene_seqs_info_schema_evaluator,
            zip(schema_files, itertools.repeat(threshold), itertools.repeat(conserved)),
            chunksize=1,
            callback=results.extend,
        )

        if show_progress is True:
            completed = False
            while completed is False:
                completed = mo.progress_bar(rawr_main, len(schema_files))

        rawr_main.wait()

        # Calculate the summary statistics for each individual locus
        print("\nCalculating individual summary statistics...\n")

        results_individual = []

        pool_ind = multiprocessing.Pool(processes=cpu_to_use)

        rawr_ind = pool_ind.map_async(
            gene_seqs_info_individual_schema_evaluator,
            schema_files,
            chunksize=1,
            callback=results_individual.extend,
        )

        if show_progress is True:
            completed = False
            while completed is False:
                completed = mo.progress_bar(rawr_ind, len(schema_files))

        rawr_ind.wait()

        pre_computed_data = {
            "mode": [],
            "total_alleles": [],
            "scatter_data": [],
        }

        pre_computed_data_individual = {}

        # Get data for each locus
        for res_ind in results_individual:

            pre_computed_data_individual[os.path.split(res_ind["gene"])[1]] = {
                "nr_alleles": res_ind["nr_alleles"],
                "size_range": res_ind["size_range"],
                "alleles_median": res_ind["median_length"],
                "locus_ids": res_ind["allele_ids"],
                "allele_sizes": res_ind["alleles_lengths"],
                "alleles_mode": res_ind["mode_length"],
            }

        # Get the data for panels A-C.

        total_number_of_loci = 0
        total_number_of_alleles = 0
        not_conserved = []
        one_allele_only = []

        for res in results:

            total_number_of_loci += 1

            total_number_of_alleles += res["nr_alleles"]

            # Get the mode for each locus.
            pre_computed_data["mode"].append(
                {
                    "locus_name": os.path.split(res["gene"])[1],
                    "alleles_mode": res["mode_length"],
                }
            )

            # Get the number of alleles for each locus.
            pre_computed_data["total_alleles"].append(
                {
                    "locus_name": os.path.split(res["gene"])[1],
                    "nr_alleles": res["nr_alleles"],
                }
            )

            # Get summary statistics (min, max, mean, mode and median) for each locus.
            pre_computed_data["scatter_data"].append(
                {
                    "locus_name": os.path.split(res["gene"])[1],
                    "locus_id": os.path.split(res["gene"])[1],
                    "nr_alleles": res["nr_alleles"],
                    "alleles_mean": res["mean_length"],
                    "alleles_median": res["median_length"],
                    "alleles_min": res["min_length"],
                    "alleles_max": res["max_length"],
                    "alleles_mode": res["mode_length"],
                }
            )

            # get the not conserved loci
            if not res["conserved"]:
                not_conserved.append({"gene": os.path.split(res["gene"])[1]})

            # get loci with only 1 allele
            if res["one_allele_only"]:
                one_allele_only.append({"gene": os.path.split(res["gene"])[1]})

        if len(not_conserved) == 0:
            not_conserved = "undefined"

        if len(one_allele_only) == 0:
            one_allele_only = "undefined"

        not_conserved_message = '"Locus size is considered not conserved if >1 allele are outside the mode +/- {0} size. Loci with only 1 allele outside the threshold are considered conserved."'.format(
            threshold
        )

        # sort pre_computed_data by locus_name
        for k in pre_computed_data:
            pre_computed_data[k] = sorted(
                pre_computed_data[k], key=itemgetter("locus_name")
            )

        # Get data for panel D
        print("\nGenerating data to populate a boxplot...\n")
        boxplot_data = gene_seqs_info_boxplot(schema_dir)

        # Get data for panel E
        print("\nAnalysing CDSs...\n")

        pool_cds = multiprocessing.Pool(processes=cpu_to_use)

        cds_multi = []

        rawr_cds = pool_cds.starmap_async(
            create_cds_df,
            zip(
                schema_files,
                itertools.repeat(minimum_length),
                itertools.repeat(minimum_length_to_translate),
                itertools.repeat(translation_table),
            ),
            chunksize=1,
            callback=cds_multi.extend,
        )

        if show_progress is True:
            completed = False
            while completed is False:
                completed = mo.progress_bar(rawr_cds, len(schema_files))

        rawr_cds.wait()

        # flatten the CDS data output list
        flat_multi_out = im.flatten_list(cds_multi)

        # sort data for the CDS table
        data_ind = sorted(
            flat_multi_out,
            key=itemgetter(
                "Alleles not multiple of 3",
                "Alleles w/ >1 stop codons",
                "Alleles wo/ Start/Stop Codon",
            ),
            reverse=True,
        )

        # sort data for CDS scatterplot
        hist_data_sort = sorted(flat_multi_out, key=itemgetter("Number of alleles"))

        # get total invalid alleles
        total_invalid_alleles = 0

        for k in data_ind:
            total_invalid_alleles += (
                k["Alleles not multiple of 3"]
                + k["Alleles w/ >1 stop codons"]
                + k["Alleles wo/ Start/Stop Codon"]
                + k["Alleles shorter than {0} nucleotides".format(minimum_length)]
            )

        # calculate total invalid alleles per class
        total_alleles_mult3 = sum(k1["Alleles not multiple of 3"] for k1 in data_ind)
        total_alleles_stopC = sum(k2["Alleles w/ >1 stop codons"] for k2 in data_ind)
        total_alleles_notStart = sum(
            k3["Alleles wo/ Start/Stop Codon"] for k3 in data_ind
        )
        total_alleles_shorter = sum(
            k3["Alleles shorter than {0} nucleotides".format(minimum_length)]
            for k3 in data_ind
        )

        # organize data for CDS scatterplot
        hist_data = {}

        hist_data["genes"] = [g["Gene"] for g in hist_data_sort]
        hist_data["total_alleles"] = [
            float(na["Number of alleles"]) for na in hist_data_sort
        ]
        hist_data["mult3"] = [
            float(mult3["Alleles not multiple of 3"]) for mult3 in hist_data_sort
        ]
        hist_data["stopC"] = [
            float(stopC["Alleles w/ >1 stop codons"]) for stopC in hist_data_sort
        ]
        hist_data["notStart"] = [
            float(notStart["Alleles wo/ Start/Stop Codon"])
            for notStart in hist_data_sort
        ]
        hist_data["shorter"] = [
            float(s["Alleles shorter than {0} nucleotides".format(minimum_length)])
            for s in hist_data_sort
        ]
        hist_data["CDS_Alleles"] = [float(cds["CDS"]) for cds in hist_data_sort]

        # check if the user provided annotations
        if annotations is None:
            uniprot_finder_missing_keys = [
                "genome",
                "contig",
                "start",
                "stop",
                "coding_strand",
                "name",
                "url",
            ]
            for d in data_ind:
                d.update(dict.fromkeys(uniprot_finder_missing_keys, "Not provided"))

        else:
            with open(annotations, "r") as a:
                annotations_reader = csv.reader(a, delimiter="\t")
                # skip header
                next(annotations_reader, None)
                annotations_data = {
                    rows[0]: [
                        rows[1],
                        rows[2],
                        rows[3],
                        rows[4],
                        rows[5],
                        rows[6],
                        rows[7],
                        rows[8],
                    ]
                    for rows in annotations_reader
                }

            for d in data_ind:
                try:
                    d["Gene"] in annotations_data
                    d.update(
                        {
                            "genome": annotations_data[d["Gene"]][0],
                            "contig": annotations_data[d["Gene"]][1],
                            "start": annotations_data[d["Gene"]][2],
                            "stop": annotations_data[d["Gene"]][3],
                            "coding_strand": "sense"
                            if annotations_data[d["Gene"]][5] == "1"
                            else "antisense",
                            "name": annotations_data[d["Gene"]][6],
                            "url": annotations_data[d["Gene"]][7],
                        }
                    )
                except KeyError:
                    uniprot_finder_missing_keys = [
                        "genome",
                        "contig",
                        "start",
                        "stop",
                        "coding_strand",
                        "name",
                        "url",
                    ]
                    d.update(dict.fromkeys(uniprot_finder_missing_keys, "-"))

        # check if it is a chewBBACA schema
        if chewie_schema:
            # read config file to get chewBBACA parameters
            config_file = os.path.join(schema_dir, ".schema_config")
            with open(config_file, "rb") as cf:
                chewie_schema_configs = pickle.load(cf)

            translation_table_config = chewie_schema_configs["translation_table"][0]
            minimum_length_config = chewie_schema_configs["minimum_locus_length"][0]

            # build the total data dictionary
            total_data = [
                {
                    "chewBBACA_version": chewie_schema_configs["chewBBACA_version"][0],
                    "bsr": chewie_schema_configs["bsr"][0],
                    "total_loci": total_number_of_loci,
                    "total_alleles": total_number_of_alleles,
                    "total_alleles_mult3": total_alleles_mult3,
                    "total_alleles_stopC": total_alleles_stopC,
                    "total_alleles_notStart": total_alleles_notStart,
                    "total_alleles_shorter": total_alleles_shorter,
                    "total_invalid_alleles": total_invalid_alleles,
                }
            ]

            # check if config parameters are the same as the user input
            if (
                minimum_length_config == minimum_length
                and translation_table_config == translation_table
            ):
                message = '"Schema created and evaluated with minimum length of {0} and translation table {1}."'.format(
                    minimum_length, translation_table
                )
            else:
                message = '"Schema created with minimum length of {0} and translation table {1} and evaluated with minimum length of {2} and translation table {3}."'.format(
                    minimum_length_config,
                    translation_table_config,
                    minimum_length,
                    translation_table,
                )
        else:
            # build the total data dictionary
            total_data = [
                {
                    "total_loci": total_number_of_loci,
                    "total_alleles": total_number_of_alleles,
                    "total_alleles_mult3": total_alleles_mult3,
                    "total_alleles_stopC": total_alleles_stopC,
                    "total_alleles_notStart": total_alleles_notStart,
                    "total_alleles_shorter": total_alleles_shorter,
                    "total_invalid_alleles": total_invalid_alleles,
                }
            ]
            message = '"Schema evaluated with minimum length of {0} and translation table {1}."'.format(
                minimum_length, translation_table
            )

        # Write HTML file
        print("\nWriting main report HTML file...\n")
        html_template_global = """
        <!DOCTYPE html>
        <html lang="en">
            <head>
                <meta charset="UTF-8" />
                <meta name="viewport" content="width=device-width, initial-scale=1.0" />
                <title>Schema Evaluator - React Edition</title>
            </head>
            <body style="background-color: #f6f6f6">
                <noscript> You need to enable JavaScript to run this app. </noscript>
                <div id="root"></div>
                <script> const _preComputedData = {0} </script>
                <script> const _preComputedDataInd = {1} </script>
                <script> const _preComputedDataBoxplot = {2} </script>
                <script> const _cdsDf = {3} </script>
                <script> const _cdsScatter = {4} </script>
                <script> const _totalData = {5} </script>
                <script> const _notConserved = {6} </script>
                <script> const _oneAlleleOnly = {7} </script>
                <script> const _message = {8} </script>
                <script> const _notConservedMessage = {9} </script>
                <script src="./main.js"></script>
            </body>
        </html>
        """.format(
            json.dumps(pre_computed_data),
            json.dumps(pre_computed_data_individual, sort_keys=True),
            json.dumps(boxplot_data),
            json.dumps(data_ind, sort_keys=True),
            json.dumps(hist_data, sort_keys=True),
            json.dumps(total_data),
            json.dumps(not_conserved),
            json.dumps(one_allele_only),
            message,
            not_conserved_message,
        )

        html_file_path = os.path.join(out_path, "schema_evaluator_report.html")

        with open(html_file_path, "w") as html_fh:
            html_fh.write(html_template_global)

        # Write the file to the pre_computed_data directory.
        pre_computed_data_path = os.path.join(out_path, "pre_computed_data.json")
        with open(pre_computed_data_path, "w") as out:
            json.dump(pre_computed_data, out)

        # Write the locus individual file to the pre_computed_data directory.
        pre_computed_data_ind_path = os.path.join(
            out_path, "pre_computed_data_ind.json"
        )
        with open(pre_computed_data_ind_path, "w") as out_ind:
            json.dump(pre_computed_data_individual, out_ind, sort_keys=True)

        # Write the boxplot pre_computed_data
        pre_computed_data_boxplot_path = os.path.join(
            out_path, "pre_computed_data_boxplot.json"
        )
        with open(pre_computed_data_boxplot_path, "w") as box_outfile:
            json.dump(boxplot_data, box_outfile)

        # Write the CDS Analysis (Panel E) data files
        cds_df_path = os.path.join(out_path, "cds_df.json")

        cds_scatter_path = os.path.join(out_path, "cds_scatter.json")

        with open(cds_df_path, "w") as cds_df_json:
            json.dump(data_ind, cds_df_json)

        with open(cds_scatter_path, "w") as cds_scatter_json:
            json.dump(hist_data, cds_scatter_json)

        return out_path
    else:
        print("Files have already been created. Moving on to the report...\n")
        return out_path


def make_protein_record(nuc_record, record_id):
    """Returns a new SeqRecord with the
    translated sequence (default table).

    Parameters
    ----------
    nuc_record : str
        protein sequence.
    record_id: str
        record id.

    Returns
    -------
    SeqRecord
        a SeqRecord object with the
        translated record.

    """
    return SeqRecord(seq=nuc_record, id=record_id, description="")


def create_protein_files(
    schema_dir,
    output_path,
    cpu_to_use,
    minimum_length,
    size_threshold,
    translation_table,
    chewie_schema=False,
    show_progress=False,
):
    """Generates FASTA files with the protein
    sequence of the schema loci.

    Parameters
    ----------
    schema_dir : list
        a list with names/paths for FASTA files.
    output_path : str
        the directory where the output files will
        be saved.
    cpu_to_use: int
        number of CPU cores to use for
        multiprocessing.
    minimum_length: int
        minimum sequence length accepted in nt.
    translation_table: int
        the translation table to be used.
    chewie_schema: bool
        identifies the schema as a chewBBACA created schema.
    show_progress: bool
        shows a progress bar for multiprocessing.

    Returns
    -------
    out_path : str
        the directory where the output files will
        be saved.
    """

    out_path = os.path.join(output_path, "prot_files")
    exception_path = os.path.join(out_path, "exceptions")

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    if not os.path.exists(exception_path):
        os.mkdir(exception_path)

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    # Check minimum length value
    if minimum_length is None:

        minimum_length = 0
        minimum_length_to_translate = minimum_length - (minimum_length * size_threshold)

        if chewie_schema:
            # read config file to get chewBBACA parameters
            config_file = os.path.join(schema_dir, ".schema_config")
            with open(config_file, "rb") as cf:
                chewie_schema_configs = pickle.load(cf)

            minimum_length = chewie_schema_configs["minimum_locus_length"][0]
            minimum_length_to_translate = minimum_length - (
                minimum_length * chewie_schema_configs["size_threshold"][0]
            )

    # set the minimum length value for translation
    minimum_length_to_translate = minimum_length - (
        minimum_length * size_threshold
    )

    print("\nTranslating....\n")

    pool_prot = multiprocessing.Pool(processes=cpu_to_use)

    rawr_prot = pool_prot.starmap_async(
        generate_protein_files,
        zip(
            schema_files,
            itertools.repeat(output_path),
            itertools.repeat(minimum_length_to_translate),
            itertools.repeat(translation_table),
        ),
        chunksize=1,
    )

    if show_progress is True:
        completed = False
        while completed is False:
            completed = mo.progress_bar(rawr_prot, len(schema_files))

    rawr_prot.wait()

    return out_path


def generate_protein_files(fasta, output_path, minimum_length, translation_table):
    """Generates FASTA files with the protein
    sequence of the schema loci.

    Parameters
    ----------
    schema_dir : list
        a list with names/paths for FASTA files.
    output_path : str
        the directory where the output files will
        be saved.
    minimum_length: int
        minimum sequence length accepted in nt.
    translation_table: int
        the translation table to be used.

    Returns
    -------
    None
    """

    out_path = os.path.join(output_path, "prot_files")
    exception_path = os.path.join(out_path, "exceptions")

    file_name_split = os.path.split(fasta)[1]
    prot_file_name = file_name_split.replace(".fasta", "_prot.fasta")
    exc_file_name = file_name_split.replace(".fasta", "_exceptions.json")

    out_file = os.path.join(out_path, prot_file_name)
    exc_file = os.path.join(exception_path, exc_file_name)

    # to be SeqRecord list
    proteins = []
    exceptions = []

    for allele in SeqIO.parse(fasta, "fasta"):
        prot = sm.translate_dna(str(allele.seq), translation_table, minimum_length)

        if isinstance(prot, list):
            tets = make_protein_record(prot[0][0], allele.id.split("_")[-1])
            proteins.append(tets)
        elif isinstance(prot, str):
            if "sense" in prot:
                prot2 = prot.split(",")[0]
            else:
                prot2 = prot
            exc = {"allele": allele.id, "exception": prot2}
            exceptions.append(exc)

    SeqIO.write(proteins, out_file, "fasta")

    with open(exc_file, "w") as ef:
        json.dump(exceptions, ef)


def call_mafft(genefile):
    """Calls MAFFT to generate an alignment.

    Parameters
    ----------
    genefile : str
        a string with the name/path for
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
        path_to_save = genefile.replace("_prot.fasta", "_aligned.fasta")
        with open(path_to_save, "w") as handle:
            handle.write(stdout)
        return True

    except Exception as e:
        print(e)
        return False


def run_mafft(protein_file_path, cpu_to_use, show_progress=False):
    """Run MAFFT with multprocessing and saves the output.

    Parameters
    ----------
    protein_file_path : str
        a string with the name/path for
        the protein FASTA file.
    cpu_to_use : int
        the number of cpu to use for
        multiprocessing.
    show_progress : bool
        If a progress bar should be
        displayed.

    Returns
    -------
    None.
    """

    print("\nRunning MAFFT...\n")

    protein_files = [
        os.path.join(protein_file_path, file)
        for file in os.listdir(protein_file_path)
        if "_prot.fasta" in file
    ]

    pool = multiprocessing.Pool(cpu_to_use)

    rawr = pool.map_async(call_mafft, protein_files, chunksize=1)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = mo.progress_bar(rawr, len(protein_files))

    rawr.wait()


def write_individual_html(
    input_files,
    pre_computed_data_path,
    protein_file_path,
    output_path,
    minimum_length,
    chewie_schema=False,
):
    """Writes HTML files for each locus.

    Parameters
    ----------
    input_files : str
        a string with the name/path for
        the schema files.
    pre_computed_data_path : str
        a string with the name/path for
        the pre-computed data files.
    protein_file_path : str
        a string with the name/path for
        the protein FASTA file.
    output_path : str
        the directory where the output
        files will be saved.
    minimum_length: int
        minimum sequence length accepted in nt.
    chewie_schema: bool
        identifies the schema as a chewBBACA created schema.

    Returns
    -------
    None.
    """

    print("\nWriting individual report HTML files...\n")

    out_path = os.path.join(output_path, "html_files")

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    schema_files = [
        os.path.splitext(f)[0] for f in os.listdir(input_files) if ".fasta" in f
    ]

    pre_computed_data_file = os.path.join(
        pre_computed_data_path, "pre_computed_data_ind.json"
    )

    cds_df_path = os.path.join(
        output_path, "SchemaEvaluator_pre_computed_data", "cds_df.json"
    )

    exceptions_path = os.path.join(
        output_path, "SchemaEvaluator_pre_computed_data", "prot_files", "exceptions"
    )

    # Read the pre_computed data file
    with open(pre_computed_data_file, "r") as pre_comp_file:
        pre_computed_data_individual = json.load(pre_comp_file)

    with open(cds_df_path, "r") as cds_file:
        cds_json_data = json.load(cds_file)

    for sf in schema_files:

        sf_fasta = "{0}.fasta".format(sf)

        # Get the precomputed data for tables and plots
        pre_computed_data_individual_sf = {
            "locus_name": str(sf),
            "data": pre_computed_data_individual[sf_fasta],
        }

        # Get CDS data for table
        cds_ind_data = [e for e in cds_json_data if sf_fasta == e["Gene"]][0]

        # Read the exceptions file
        exceptions_filename_path = os.path.join(
            exceptions_path, "{0}_exceptions.json".format(sf)
        )
        with open(exceptions_filename_path, "r") as ef:
            exc_data = json.load(ef)

        # get the msa data
        msa_file_path = os.path.join(protein_file_path, "{0}_aligned.fasta".format(sf))

        msa_data = {"sequences": []}

        msa_seq_gen = SeqIO.parse(msa_file_path, "fasta")

        allele_ids = [allele.id for allele in msa_seq_gen]

        if len(allele_ids) > 1:
            for allele in SeqIO.parse(msa_file_path, "fasta"):
                msa_data["sequences"].append(
                    {"name": allele.id, "sequence": str(allele.seq)}
                )
        else:
            msa_data = "undefined"

        # get the phylocanvas data
        phylo_file_path = os.path.join(
            protein_file_path, "{0}_prot.fasta.tree".format(sf)
        )
        if os.path.exists(phylo_file_path):
            with open(phylo_file_path, "r") as phylo:
                phylo_data = phylo.read()

            phylo_data_json = {"phylo_data": phylo_data}
        else:
            phylo_data_json = "undefined"

        if minimum_length is None:

            minimum_length = 0
            if chewie_schema:
                # read config file to get chewBBACA parameters
                config_file = os.path.join(input_files, ".schema_config")
                with open(config_file, "rb") as cf:
                    chewie_schema_configs = pickle.load(cf)

                minimum_length = chewie_schema_configs["minimum_locus_length"][0]

        html_template_individual = """
            <!DOCTYPE html>
            <html lang="en">
                <head>
                    <meta charset="UTF-8" />
                    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
                    <title>Schema Evaluator - Individual Analysis</title>
                </head>
                <body style="background-color: #f6f6f6">
                    <noscript> You need to enable JavaScript to run this app. </noscript>
                    <div id="root"></div>
                    <script> const _preComputedDataInd = {0} </script>
                    <script> const _exceptions = {1} </script>
                    <script> const _cdsDf = {2} </script>
                    <script> const _msaData = {3} </script>
                    <script> const _phyloData = {4} </script>
                    <script> const _minLen = {5} </script>
                    <script src="./main_ind.js"></script>
                </body>
            </html>
            """.format(
            json.dumps(pre_computed_data_individual_sf, sort_keys=True),
            json.dumps(exc_data, sort_keys=True),
            json.dumps(cds_ind_data, sort_keys=True),
            json.dumps(msa_data, sort_keys=True),
            json.dumps(phylo_data_json, sort_keys=True),
            json.dumps(int(minimum_length), sort_keys=True),
        )

        html_file_path = os.path.join(out_path, "{0}_individual_report.html".format(sf))

        with open(html_file_path, "w") as html_fh:
            html_fh.write(html_template_individual)

        # if chewie_schema:
        #     # read config file to get chewBBACA parameters
        #     config_file = os.path.join(input_files, ".schema_config")
        #     with open(config_file, "rb") as cf:
        #         chewie_schema_configs = pickle.load(cf)

        #     minimum_length_config = chewie_schema_configs["minimum_locus_length"][0]

        #     html_template_individual = """
        #     <!DOCTYPE html>
        #     <html lang="en">
        #         <head>
        #             <meta charset="UTF-8" />
        #             <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        #             <title>Schema Evaluator - Individual Analysis</title>
        #         </head>
        #         <body style="background-color: #f6f6f6">
        #             <noscript> You need to enable JavaScript to run this app. </noscript>
        #             <div id="root"></div>
        #             <script> const _preComputedDataInd = {0} </script>
        #             <script> const _exceptions = {1} </script>
        #             <script> const _cdsDf = {2} </script>
        #             <script> const _msaData = {3} </script>
        #             <script> const _phyloData = {4} </script>
        #             <script> const _minLen = {5} </script>
        #             <script src="./main_ind.js"></script>
        #         </body>
        #     </html>
        #     """.format(
        #         json.dumps(pre_computed_data_individual_sf, sort_keys=True),
        #         json.dumps(exc_data, sort_keys=True),
        #         json.dumps(cds_ind_data, sort_keys=True),
        #         json.dumps(msa_data, sort_keys=True),
        #         json.dumps(phylo_data_json, sort_keys=True),
        #         json.dumps(int(minimum_length_config), sort_keys=True),
        #     )

        #     html_file_path = os.path.join(
        #         out_path, "{0}_individual_report.html".format(sf)
        #     )

        #     with open(html_file_path, "w") as html_fh:
        #         html_fh.write(html_template_individual)

        # else:

        # html_template_individual = """
        # <!DOCTYPE html>
        # <html lang="en">
        #     <head>
        #         <meta charset="UTF-8" />
        #         <meta name="viewport" content="width=device-width, initial-scale=1.0" />
        #         <title>Schema Evaluator - Individual Analysis</title>
        #     </head>
        #     <body style="background-color: #f6f6f6">
        #         <noscript> You need to enable JavaScript to run this app. </noscript>
        #         <div id="root"></div>
        #         <script> const _preComputedDataInd = {0} </script>
        #         <script> const _exceptions = {1} </script>
        #         <script> const _cdsDf = {2} </script>
        #         <script> const _msaData = {3} </script>
        #         <script> const _phyloData = {4} </script>
        #         <script> const _minLen = {5} </script>
        #         <script src="./main_ind.js"></script>
        #     </body>
        # </html>
        # """.format(
        #     json.dumps(pre_computed_data_individual_sf, sort_keys=True),
        #     json.dumps(exc_data, sort_keys=True),
        #     json.dumps(cds_ind_data, sort_keys=True),
        #     json.dumps(msa_data, sort_keys=True),
        #     json.dumps(phylo_data_json, sort_keys=True),
        #     json.dumps(int(minimum_length), sort_keys=True),
        # )

        # html_file_path = os.path.join(
        #     out_path, "{0}_individual_report.html".format(sf)
        # )

        # with open(html_file_path, "w") as html_fh:
        #     html_fh.write(html_template_individual)
