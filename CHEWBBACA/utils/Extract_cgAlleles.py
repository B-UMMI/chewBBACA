#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module determines the set of genes in the core genome based on
a matrix with allelic profiles and a threshold that defines the
proportion of genomes a gene must be present in to be included in
the core genome.

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_file`` : Path to input file containing a matrix with
  allelic profiles.

    - e.g.: ``/home/user/chewie/results/matrix``

- ``-o``, ``output_directory`` : Path to the directory where the process
  will store output files.

    - e.g.: ``/home/user/chewie/results/output_directory``

- ``--p``, ``threshold`` : Genes that constitute the core genome must be
  in a proportion of genomes that is at least equal to this value.

    - e.g.: ``0.95``

- ``--r``, ``genes2remove`` : Path to file with a list of genes/columns to
  remove from the matrix (one gene identifier per line).

    - e.g.: ``home/user/results/genes.txt``

- ``--g``, ``genomes2remove`` : Path to file with a list of genomes/rows
  to remove from the matrix (one genome identifier per line).

Code documentation
------------------
"""


import os
import numpy as np
import pandas as pd


def replace_chars(column):
    """ Replaces all non-numeric characters
        in a column with allele identifiers.

        Parameters
        ----------
        column : pandas.core.series.Series
            Pandas dataframe column.

        Returns
        -------
        replace_missing : pandas.core.series.Series
            Input column with cells that only contain
            numeric characters.
    """

    replace_inf = column.replace(to_replace='INF-',
                                 value='', regex=True)
    # replace '*' in novel alleles from schemas in Chewie-NS
    # double escape due to warning that might become error in future
    replace_inf = replace_inf.replace(to_replace='\\*',
                                      value='', regex=True)
    replace_missing = replace_inf.replace(to_replace='\\D+.*',
                                          value='0', regex=True)

    return replace_missing


def binarize_matrix(column):
    """ Converts a Pandas dataframe column values
        into numeric values.

        Parameters
        ----------
        column : pandas.core.series.Series
            Pandas dataframe column.

        Returns
        -------
        Numpy array corresponding to the input column
        with numeric values equal to 1 for cells that
        had valid allele identifiers and equal to 0
        for cells that had missing data.
    """

    replace_missing = replace_chars(column)
    coln = pd.to_numeric(replace_missing)

    return np.int64(coln > 0)


def above_threshold(column, column_length, threshold):
    """ Determines if a gene is present in a proportion
        of genomes equal or greater than a threshold.

        Parameters
        ----------
        column : pandas.core.series.Series
            Pandas dataframe column.
        column_length : int
            Number of cells/genes in the column.
        threshold : float
            Core genome determination threshold.

        Returns
        -------
        bool
            A boolean, True if gene is equal or above
            threshold, False otherwise.
    """

    return (np.sum(column) / column_length) >= threshold


def remove_genomes(matrix, genomesToRemove):
    """ Removes rows from a Pandas dataframe if the
        index identifier matches the identifier of
        a genome to remove.

        Parameters
        ----------
        matrix : pandas.core.frame.DataFrame
            Pandas dataframe with allelic profiles.
            Each row has the allelic profile of a genome
            and each column has the allele identifiers
            determined for a single gene.
        genomesToRemove : list
            List with the set of genomes to remove.

        Returns
        -------
        pruned_matrix : pandas.core.frame.DataFrame
            Input dataframe without the rows whose
            index matched an identifier of a genome
            to remove.
    """

    to_remove_bool = matrix.index.isin(genomesToRemove)
    pruned_matrix = matrix.loc[~ to_remove_bool]

    # remove genomes
    for genome in matrix.index:
        if genome in genomesToRemove:
            print('Removed genome: {0}'.format(genome))

    return pruned_matrix


def remove_genes(matrix, presence_absence, genesToRemove, threshold):
    """ Determines genes that are in a proportion of genomes
        above or below a threshold.

        Parameters
        ----------
        matrix : pandas.core.frame.DataFrame
            Pandas dataframe with allelic profiles.
            Each row has the allelic profile of a genome
            and each column has the allele identifiers
            determined for a single gene.
        presence_absence : pandas.core.frame.DataFrame
            Pandas dataframe with numeric values equal to
            1 for the cells that had valid allele identifiers
            and equal to 0 for missing data.
        genesToRemove : list
            List with a set of genes to exclude from the core
            genome.
        threshold : float
            Core genome determination threshold.

        Returns
        -------
        pruned_matrix : pandas.core.frame.DataFrame
            Input dataframe without the columns whose
            headers matched an identifier of a gene
            to remove or that was below the threshold.
        genes_to_delete : set
            Set with identifiers of genes that were not included
            in the core genome.
    """

    # determine genes at or above threshold
    pa_rows, _ = presence_absence.shape
    is_above_threshold = presence_absence.apply(above_threshold,
                                                args=(pa_rows, threshold,))
    below_threshold = matrix.columns[~ is_above_threshold]

    genes_to_delete = set(genesToRemove).union(set(below_threshold))

    # get columns with genes to keep
    to_keep = ~matrix.columns.isin(genes_to_delete)
    pruned_matrix = matrix.loc[:, to_keep]

    return [pruned_matrix, genes_to_delete]


def presAbs(matrix, output_directory):
    """ Creates a presence absence matrix.

        Parameters
        ----------
        matrix : pandas.core.frame.DataFrame
            Pandas dataframe with allelic profiles.
            Each row has the allelic profile of a genome
            and each column has the allele identifiers
            determined for a single gene.
        output_directory : str
            Path to the directory where the TSV file with
            the presence absence matrix will be stored.

        Returns
        -------
        presence_absence : pandas.core.frame.DataFrame
            Pandas dataframe with numeric values equal to
            1 for the cells that had valid allele identifiers
            and equal to 0 for missing data.
    """

    presence_absence = matrix.apply(binarize_matrix)

    pa_path = os.path.join(output_directory, 'Presence_Absence.tsv')
    presence_absence.to_csv(pa_path, sep='\t')

    return presence_absence


def missing_data_table(presence_absence):
    """ Determines missing data per genome.

        Parameters
        ----------
        presence_absence : pandas.core.frame.DataFrame
            Pandas dataframe with numeric values equal to
            1 for the cells that have valid allele identifiers
            and equal to 0 for missing data.

        Returns
        -------
        missing_data_df : pandas.core.frame.DataFrame
            Dataframe with number of missing genes and
            percentage of missing genes per genome.
    """

    _, n_genes = presence_absence.shape
    genes_present = presence_absence.apply(np.count_nonzero, axis=1)

    missing_data = {'FILE': presence_absence.index,
                    'missing': n_genes - genes_present,
                    'percentage': 1 - (genes_present / n_genes)}

    missing_data_df = pd.DataFrame(missing_data,
                                   columns=['FILE',
                                            'missing',
                                            'percentage'])

    return missing_data_df


def determine_cgMLST(input_file, output_directory, genesToRemove,
                     genomesToRemove, threshold):
    """ Determines the cgMLST based on an input matrix of allelic
        profiles.

        Parameters
        ----------
        input_file : str
            Path a TSV file with allelic profiles for a set
            of genomes.
        output_directory : str
            Path to the directory where the process will
            store output files.
        genesToRemove : list
            List with a set of genes to remove from the
            analysis.
        genomesToRemove : list
            List with a set of genomes to remove from the
            analysis.
        threshold : float
            Core genome determination threshold.

        Returns
        -------
        List with the paths to three files:

        - Path a TSV file with the cgMLST matrix.
        - Path to a TXT file with the list of genes that
          constitute the core genome.
        - Path to a TSV file with the information about
          missing data per genome.
    """

    # import matrix with allelic profiles
    matrix = pd.read_csv(input_file, header=0, index_col=0,
                         sep='\t', low_memory=False)

    # remove genomes
    genome_pruned = remove_genomes(matrix, genomesToRemove)

    # build presence/absence matrix
    print('\nBuilding presence and absence matrix...', end='')
    presence_absence = presAbs(genome_pruned, output_directory)
    print('done.')

    # remove genes
    print('Determining genes in the core genome...', end='')
    gene_pruned, removed_genes = remove_genes(genome_pruned, presence_absence,
                                              genesToRemove, threshold)
    print('done.')

    # count number of missing data per genome
    print('Determining missing data per genome...', end='')
    missing_data_df = missing_data_table(presence_absence)
    print('done.')

    # write cgMLST matrix
    cgmlst_path = os.path.join(output_directory, 'cgMLST.tsv')
    gene_pruned.to_csv(cgmlst_path, sep='\t')

    # write genes in cgMLST to file
    loci_path = os.path.join(output_directory, 'cgMLSTschema.txt')
    pd.Series(list(gene_pruned.columns.values)).to_csv(loci_path,
                                                       index=False,
                                                       header=False)

    # write data with missing data stats
    mdata_path = os.path.join(output_directory, 'mdata_stats.tsv')
    missing_data_df.to_csv(mdata_path, sep='\t', index=False)

    retained = len(gene_pruned.columns)
    print('\nCore genome composed of {0}/{1} genes.'
          ''.format(retained, len(matrix.columns)))

    return [cgmlst_path, loci_path, mdata_path]


def main(input_file, output_directory, threshold,
         genes2remove, genomes2remove):

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    genesToRemove = []
    if genes2remove:
        with open(genes2remove, 'r') as gr:
            genesToRemove = gr.read().splitlines()
    genomesToRemove = []
    if genomes2remove:
        with open(genomes2remove, 'r') as gr:
            genomesToRemove = gr.read().splitlines()

    determine_cgMLST(input_file, output_directory, genesToRemove,
                     genomesToRemove, threshold)


if __name__ == "__main__":

    main()
