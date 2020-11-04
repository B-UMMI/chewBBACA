import os
import json
import statistics
import pandas as pd
from operator import itemgetter
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from CHEWBBACA.utils import auxiliary_functions as aux


# Schema Evaluator Auxiliary Functions

def gene_seqs_info_schema_evaluator(gene):
    """ Determines the total number of alleles and the mean length
        of allele sequences per gene.

        Args:
            genes_list (list): a list with names/paths for FASTA
            files.
        Returns:
            genes_info (list): a list with a sublist for each input
            gene file. Each sublist contains a gene identifier, the
            total number of alleles for that gene and the mean length
            of allele sequences for that gene.
    """

    seq_generator = SeqIO.parse(gene, "fasta")
    alleles_lengths = [len(allele) for allele in seq_generator]
    alleles_lengths.sort()

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
    """
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
    """Determines boxplot statistics.

    Args:
        genes_list (list): a list with names/paths for FASTA
        files.
    Returns:
        genes_info (dict): a dict with a subdict for each input
        gene file. Each subdict contains information about each
        gene such as, mode of the allele sizes, number of alleles
        and summary statistics (min, max, median and mode).
    """

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    loci_list = []
    loci_min = []
    loci_q1 = []
    loci_median = []
    loci_q3 = []
    loci_max = []
    loci_mean = []
    loci_sd = []
    alleles_counts = []

    for g in schema_files:
        seq_generator = SeqIO.parse(g, "fasta")
        alleles_lengths = [len(allele) for allele in seq_generator]
        alleles_lengths.sort()

        # locus name
        locus_name = os.path.split(g)[1]
        loci_list.append(locus_name)

        # number of alleles
        nr_alleles = len(alleles_lengths)
        alleles_counts.append(nr_alleles)

        # minimum and maximum values
        loci_max.append(max(alleles_lengths))
        loci_min.append(min(alleles_lengths))
        # standard deviation
        if nr_alleles > 1:
            locus_sd = statistics.stdev(alleles_lengths)
        else:
            locus_sd = 0.0
        loci_sd.append(locus_sd)

        # # median
        median_length = round(statistics.median(alleles_lengths))
        loci_median.append(median_length)
        # # mean
        mean_length = round(sum(alleles_lengths) / nr_alleles)
        loci_mean.append(mean_length)
        # q1 and q3
        if nr_alleles > 1:
            half = int(nr_alleles // 2)
            q1 = statistics.median(alleles_lengths[:half])
            q3 = statistics.median(alleles_lengths[-half:])
        else:
            q1 = alleles_lengths[0]
            q3 = alleles_lengths[0]
        loci_q1.append(q1)
        loci_q3.append(q3)

    json_to_file = {
        "loci": loci_list,
        "min": loci_min,
        "q1": loci_q1,
        "median": loci_median,
        "q3": loci_q3,
        "max": loci_max,
        "mean": loci_mean,
        "sd": loci_sd,
        "nr_alleles": alleles_counts,
    }

    return json_to_file


def create_pre_computed_data(schema_dir, output_path):
    """Creates a file with pre-computed data for the Schema Evaluator plotly charts."""

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    out_path = os.path.join(output_path, "SchemaEvaluator_pre_computed_data")
    if not os.path.exists(out_path):
        os.mkdir(out_path)

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
            {"locus_name": os.path.split(res[0])[1], "alleles_mode": res[-1]}
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

    # Get data for panel D
    boxplot_data = gene_seqs_info_boxplot(schema_dir)

    # Write the file to the (future) pre_computed_data directory.
    pre_computed_data_path = os.path.join(out_path, "pre_computed_data.json")
    with open(pre_computed_data_path, "w") as out:
        json.dump(pre_computed_data, out)

    # Write the locus individual file to the pre_computed_data directory.
    pre_computed_data_ind_path = os.path.join(
        out_path, "pre_computed_data_ind.json")
    with open(pre_computed_data_ind_path, "w") as out_ind:
        json.dump(pre_computed_data_individual, out_ind)

    # Write the boxplot pre_computed_data
    pre_computed_data_boxplot_path = os.path.join(
        out_path, "pre_computed_data_boxplot.json")
    with open(pre_computed_data_boxplot_path, "w") as box_outfile:
        json.dump(boxplot_data, box_outfile)

    return out_path


def create_cds_df(schema_dir):
    """


    Parameters
    ----------
    schema_dir : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    schema_files = [
        os.path.join(schema_dir, file)
        for file in os.listdir(schema_dir)
        if ".fasta" in file
    ]

    res = {"stats": []}

    for f in schema_files:

        gene_res = {"Gene": os.path.split(f)[1]}

        gene_res["Number of alleles"] = aux.count_sequences(f)

        stopC = 0
        notStart = 0
        notMultiple = 0
        CDS = 0

        for allele in SeqIO.parse(f, "fasta"):

            ola = aux.translate_dna(str(allele.seq), 11, 201)

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

    return data_index, hist_data


def create_locus_ind_df(locus_ind_dict, option):
    """


    Parameters
    ----------
    locus_ind_dict : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    option_dict = {option: locus_ind_dict[option]}

    df_to_drop_cols = pd.DataFrame.from_dict(option_dict, orient="index")

    df_to_rename = df_to_drop_cols.drop(["locus_ids", "allele_sizes"], axis=1)

    columns_to_rename = {
        "nr_alleles": "Total Number of Alleles",
        "size_range": "Size Range (bp)",
        "alleles_median": "Median Size (bp)",
        "alleles_mode": "Alleles Mode"}

    df_to_return = df_to_rename.rename(columns=columns_to_rename)

    return df_to_return


def make_protein_record(nuc_record, record_id):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(
        seq=nuc_record,
        id="trans_" + record_id,
    )


def create_protein_files(schema_dir, output_path):
    """


    Parameters
    ----------
    loci_list : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

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

    print("Translating....")

    for f in schema_files:

        # print(f)
        file_name_split = os.path.split(f)[1]
        prot_file_name = file_name_split.replace(".fasta", "_prot.fasta")
        exc_file_name = file_name_split.replace(".fasta", "_exceptions.txt")

        out_file = os.path.join(out_path, prot_file_name)
        exc_file = os.path.join(exception_path, exc_file_name)

        # to be SeqRecord list
        proteins = []
        exceptions = []

        for allele in SeqIO.parse(f, "fasta"):
            ola = aux.translate_dna(str(allele.seq), 11, 201)

            if isinstance(ola, list):
                tets = make_protein_record(ola[0][0], allele.id)
                proteins.append(tets)
            elif isinstance(ola, str):
                exc = [allele.id, ola]
                exceptions.append(exc)


        SeqIO.write(proteins, out_file, "fasta")

        with open(exc_file, "w") as ef:
            ef.writelines('\n'.join([': '.join(map(str, ex))
                                     for ex in exceptions]))


    print("Done!")
