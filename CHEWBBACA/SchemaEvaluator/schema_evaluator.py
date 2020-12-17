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
import statistics
import pandas as pd
import multiprocessing
from operator import itemgetter
from collections import Counter


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalwCommandline

try:
    from utils import auxiliary_functions as aux
except:
    from CHEWBBACA.utils import auxiliary_functions as aux


# Schema Evaluator Auxiliary Functions

def gene_seqs_info_schema_evaluator(gene):
    """ Determines the total number of alleles and the mean length
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

    genes_info = [
        gene,
        nr_alleles,
        min_length,
        max_length,
        median_length,
        mean_length,
        mode_length,
        alleles_lengths
    ]

    return genes_info


def gene_seqs_info_individual_schema_evaluator(gene):
    """ Determines the total number of alleles and 
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
    size_range = f"{min_length}-{max_length}"

    # Summary statistics
    median_length = round(statistics.median(alleles_lengths))
    mode_length = Counter(alleles_lengths).most_common()[0][0]

    genes_info = [
        gene,
        nr_alleles,
        size_range,
        median_length,
        allele_ids,
        alleles_lengths,
        mode_length,
    ]

    return genes_info


def gene_seqs_info_boxplot(schema_dir):
    """ Determines boxplot statistics.

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

        json_to_file["boxplot_data"].append({
            "locus_name": locus_name,
            "nr_alleles": nr_alleles,
            "max": loci_max,
            "min": loci_min,
            "sd": locus_sd,
            "median": median_length,
            "mean": mean_length,
            "q1": q1,
            "q3": q3
        })

    # Sort data by locus_name
    for k in json_to_file:
        json_to_file[k] = sorted(json_to_file[k], key=itemgetter("locus_name"))

    return json_to_file


# Functions that obtain the data for panel E
def create_cds_df(schema_dir, translation_table):
    """ Detects alleles that aren't CDSs.

        Parameters
        ----------
        schema_dir : list 
            a list with names/paths for FASTA files.
        translation_table: int
            the translation table to be used.

        Returns
        -------
        data_index_records_dict : dict
            a dict obtained from a dataframe containing the
            number of non-CDSs detected. It will be used to
            populate a table in the report.
        hist_data : dict
            a dict obtained from a dataframe containing the
            number of non-CDSs detected. It will be used to
            populate a chart in the report.

    """

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    # schema_files.sort()

    res = {"stats": []}

    for f in schema_files:

        gene_res = {"Gene": os.path.split(f)[1]}

        gene_res["Number of alleles"] = aux.count_sequences(f)

        stopC = 0
        notStart = 0
        notMultiple = 0
        CDS = 0

        for allele in SeqIO.parse(f, "fasta"):

            ola = aux.translate_dna(str(allele.seq), translation_table, 201)

            if "sequence length is not a multiple of 3" in ola:
                notMultiple += 1
            elif "Extra in frame stop codon found" in ola:
                stopC += 1
            elif "is not a start codon" in ola:
                notStart += 1
            elif "is not a stop codon" in ola:
                notStart += 1
            else:
                CDS += 1

        gene_res["Alleles not multiple of 3"] = notMultiple
        gene_res["Alleles w/ >1 stop codons"] = stopC
        gene_res["Alleles wo/ Start/Stop Codon"] = notStart
        gene_res["CDS"] = CDS

        res["stats"].append(gene_res)

        res_sorted = sorted(res["stats"], key=itemgetter("Number of alleles"))

        hist_data = {}

        hist_data["genes"] = [g["Gene"] for g in res_sorted]
        hist_data["total_alleles"] = [
            float(na["Number of alleles"]) for na in res_sorted]
        hist_data["mult3"] = [
            float(mult3["Alleles not multiple of 3"]) for mult3 in res_sorted]
        hist_data["stopC"] = [
            float(stopC["Alleles w/ >1 stop codons"]) for stopC in res_sorted]
        hist_data["notStart"] = [
            float(notStart["Alleles wo/ Start/Stop Codon"]) for notStart in res_sorted]
        hist_data["CDS_Alleles"] = [float(cds["CDS"]) for cds in res_sorted]

    res_sorted_reverse = sorted(res["stats"], key=itemgetter(
        "Number of alleles"), reverse=True)

    data_index = pd.DataFrame.from_dict(res_sorted_reverse)

    data_index = data_index.sort_values(
        ["Alleles not multiple of 3", "Alleles w/ >1 stop codons",
            "Alleles wo/ Start/Stop Codon"],
        ascending=[False, False, False]
    )

    data_index_records_dict = data_index.to_dict(orient="records")

    return data_index_records_dict, hist_data


def create_pre_computed_data(schema_dir, translation_table, output_path):
    """ Creates a file with pre-computed data for 
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

    # schema_files.sort()

    if len(schema_files) < 1:
        sys.exit("The schema directory is empty. Please check your path. Exiting...")

    empty_files_1 = [aux.is_file_empty(f) for f in schema_files]

    if True in empty_files_1:
        sys.exit("At least one file is empty or it doesn't exist. Exiting...")

    out_path = os.path.join(output_path, "SchemaEvaluator_pre_computed_data")
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if not os.listdir(out_path):
        # Calculate the summary statistics and other information about each locus.
        results = [gene_seqs_info_schema_evaluator(i) for i in schema_files]

        # Calculate the summary statistics for each individual locus
        results_individual = [
            gene_seqs_info_individual_schema_evaluator(i) for i in schema_files]

        pre_computed_data = {
            "mode": [],
            "total_alleles": [],
            "scatter_data": [],
        }

        pre_computed_data_individual = {}

        # Get data for each locus
        for res_ind in results_individual:

            pre_computed_data_individual[os.path.split(res_ind[0])[1]] = {
                "nr_alleles": res_ind[1],
                "size_range": res_ind[2],
                "alleles_median": res_ind[3],
                "locus_ids": res_ind[4],
                "allele_sizes": res_ind[5],
                "alleles_mode": res_ind[6],
            }

        # Get the data for panels A-C.
        for res in results:

            # Get the mode for each locus.
            pre_computed_data["mode"].append(
                {"locus_name": os.path.split(
                    res[0])[1], "alleles_mode": res[-1]}
            )

            # Get the number of alleles for each locus.
            pre_computed_data["total_alleles"].append(
                {"locus_name": os.path.split(res[0])[1], "nr_alleles": res[1]}
            )

            # Get summary statistics (min, max, mean, mode and median) for each locus.
            pre_computed_data["scatter_data"].append(
                {
                    "locus_name": os.path.split(res[0])[1],
                    "locus_id": os.path.split(res[0])[1],
                    "nr_alleles": res[1],
                    "alleles_mean": res[5],
                    "alleles_median": res[4],
                    "alleles_min": res[2],
                    "alleles_max": res[3],
                    "alleles_mode": res[6],
                }
            )

        # sort pre_computed_data by locus_name
        for k in pre_computed_data:
            pre_computed_data[k] = sorted(
                pre_computed_data[k], key=itemgetter("locus_name"))

        # Get data for panel D
        boxplot_data = gene_seqs_info_boxplot(schema_dir)

        # Get data for panel E
        data_ind, hist_data = create_cds_df(schema_dir, translation_table)

        # Write HTML file
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
                <script src="./main.js"></script>
            </body>
        </html>
        """.format(
            json.dumps(pre_computed_data),
            json.dumps(pre_computed_data_individual, sort_keys=True),
            json.dumps(boxplot_data),
            json.dumps(data_ind, sort_keys=True),
            json.dumps(hist_data, sort_keys=True)
        )

        html_file_path = os.path.join(
            out_path, "schema_evaluator_report.html"
        )

        with open(html_file_path, "w") as html_fh:
            html_fh.write(html_template_global)

        # Write the file to the pre_computed_data directory.
        pre_computed_data_path = os.path.join(
            out_path, "pre_computed_data.json")
        with open(pre_computed_data_path, "w") as out:
            json.dump(pre_computed_data, out)

        # Write the locus individual file to the pre_computed_data directory.
        pre_computed_data_ind_path = os.path.join(
            out_path, "pre_computed_data_ind.json")
        with open(pre_computed_data_ind_path, "w") as out_ind:
            json.dump(pre_computed_data_individual, out_ind, sort_keys=True)

        # Write the boxplot pre_computed_data
        pre_computed_data_boxplot_path = os.path.join(
            out_path, "pre_computed_data_boxplot.json")
        with open(pre_computed_data_boxplot_path, "w") as box_outfile:
            json.dump(boxplot_data, box_outfile)

        # Write the CDS Analysis (Panel E) data files
        cds_df_path = os.path.join(
            out_path, "cds_df.json"
        )

        cds_scatter_path = os.path.join(
            out_path, "cds_scatter.json"
        )

        with open(cds_df_path, "w") as cds_df_json:
            json.dump(data_ind, cds_df_json)

        with open(cds_scatter_path, "w") as cds_scatter_json:
            json.dump(hist_data, cds_scatter_json)

        return out_path
    else:
        print("Files have already been created. Moving on to the report...\n")
        return out_path


def make_protein_record(nuc_record, record_id):
    """ Returns a new SeqRecord with the 
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
    return SeqRecord(
        seq=nuc_record,
        id="trans_" + record_id,
    )


def create_protein_files(schema_dir, output_path):
    """ Generates FASTA files with the protein
        sequence of the schema loci.

        Parameters
        ----------
        schema_dir : list 
            a list with names/paths for FASTA files.
        output_path : str
            the directory where the output files will
            be saved.

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

    # if not os.listdir(out_path):
    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    print("Translating....\n")

    for f in schema_files:

        # print(f)
        file_name_split = os.path.split(f)[1]
        prot_file_name = file_name_split.replace(".fasta", "_prot.fasta")
        exc_file_name = file_name_split.replace(".fasta", "_exceptions.json")

        out_file = os.path.join(out_path, prot_file_name)
        exc_file = os.path.join(exception_path, exc_file_name)

        # to be SeqRecord list
        proteins = []
        exceptions = []

        for allele in SeqIO.parse(f, "fasta"):
            prot = aux.translate_dna(str(allele.seq), 11, 201)

            if isinstance(prot, list):
                tets = make_protein_record(prot[0][0], allele.id)
                proteins.append(tets)
            elif isinstance(prot, str):
                # exc = [allele.id, ola]
                if "sense" in prot:
                    prot2 = prot.split(",")[0]
                else:
                    prot2 = prot
                exc = {
                    "allele": allele.id,
                    "exception": prot2
                }
                exceptions.append(exc)

        SeqIO.write(proteins, out_file, "fasta")

        with open(exc_file, "w") as ef:
            json.dump(exceptions, ef)

    print("Done!")

    return out_path


def call_mafft(genefile):
    """ Calls MAFFT to generate an alignment.

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
        mafft_cline = MafftCommandline(input=genefile, adjustdirection=True)
        stdout, stderr = mafft_cline()
        path_to_save = genefile.replace("_prot.fasta", "_aligned.fasta")
        with open(path_to_save, "w") as handle:
            handle.write(stdout)
        return True

    except Exception as e:
        print(e)
        return False


def run_mafft(protein_file_path, cpu_to_use, show_progress=False):
    """ Run MAFFT with multprocessing and saves the output.

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

    print("Running MAFFT...\n")

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
            completed = aux.progress_bar(rawr, len(protein_files))

    rawr.wait()


def call_clustalw(genefile):
    """ Call ClustalW to generate a Neighbour Joining tree.

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
        clw_cline = ClustalwCommandline(
            "clustalw", infile=genefile, tree=True)
        clw_cline()

        return True

    except Exception as e:
        print(e)
        return False


def run_clustalw(protein_file_path, cpu_to_use, show_progress=False):
    """ Run ClustalW with multiprocessing and save output.

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

    print("\nRunning Clustal...\n")

    protein_files = [
        os.path.join(protein_file_path, file)
        for file in os.listdir(protein_file_path)
        if "_aligned.fasta" in file
    ]

    pool = multiprocessing.Pool(cpu_to_use)

    rawr = pool.map_async(call_clustalw, protein_files, chunksize=1)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = aux.progress_bar(rawr, len(protein_files))

    rawr.wait()


def write_individual_html(input_files, pre_computed_data_path, protein_file_path, output_path):
    """ Writes HTML files for each locus.

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

        Returns
        -------
        None.
    """

    print("\nWriting HTML files...")

    out_path = os.path.join(output_path, "html_files")

    if not os.path.exists(out_path):
        os.mkdir(out_path)

    schema_files = [
        os.path.splitext(f)[0]
        for f in os.listdir(input_files)
        if ".fasta" in f
    ]

    pre_computed_data_file = os.path.join(
        pre_computed_data_path, "pre_computed_data_ind.json")

    cds_df_path = os.path.join(
        output_path, "SchemaEvaluator_pre_computed_data",  "cds_df.json"
    )

    exceptions_path = os.path.join(
        output_path, "SchemaEvaluator_pre_computed_data",  "prot_files", "exceptions"
    )

    # Read the pre_computed data file
    with open(pre_computed_data_file, "r") as pre_comp_file:
        pre_computed_data_individual = json.load(pre_comp_file)

    with open(cds_df_path, "r") as cds_file:
        cds_json_data = json.load(cds_file)

    for sf in schema_files:

        # Get the precomputed data for tables and plots
        pre_computed_data_individual_sf = {"locus_name": str(
            sf), "data": pre_computed_data_individual[f"{sf}.fasta"]}

        # Get CDS data for table
        cds_ind_data = [e for e in cds_json_data if sf in e["Gene"]][0]

        # print(json.dumps(cds_ind_data, sort_keys=True))

        # Read the exceptions file
        exceptions_filename_path = os.path.join(
            exceptions_path, f"{sf}_exceptions.json")
        with open(exceptions_filename_path, "r") as ef:
            exc_data = json.load(ef)

        # get the msa data
        msa_file_path = os.path.join(protein_file_path, f"{sf}_aligned.fasta")

        msa_data = {"sequences": []}

        for allele in SeqIO.parse(msa_file_path, "fasta"):
            msa_data["sequences"].append(
                {"name": allele.id, "sequence": str(allele.seq)})

        # get the phylocanvas data
        phylo_file_path = os.path.join(protein_file_path, f"{sf}_aligned.ph")
        # print(phylo_file_path)
        if os.path.exists(phylo_file_path):
            with open(phylo_file_path, "r") as phylo:
                phylo_data = phylo.read()

            phylo_data_json = {"phylo_data": phylo_data}
        else:
            phylo_data_json = []

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
                <script src="./main_ind.js"></script>
            </body>
        </html>
        """.format(
            json.dumps(pre_computed_data_individual_sf, sort_keys=True),
            json.dumps(exc_data, sort_keys=True),
            json.dumps(cds_ind_data, sort_keys=True),
            json.dumps(msa_data, sort_keys=True),
            json.dumps(phylo_data_json, sort_keys=True)
        )

        html_file_path = os.path.join(
            out_path, f"{sf}_individual_report.html"
        )

        with open(html_file_path, "w") as html_fh:
            html_fh.write(html_template_individual)
