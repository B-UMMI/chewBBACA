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
     'allelic profiles.

    - e.g.: ``/home/user/chewie/results/matrix``

- ``-o``, ``output_directory`` : Path to the directory where the process
  will store output files.

    - e.g.: ``/home/user/chewie/results/output_directory``

- ``--t``, ``threshold`` : Genes that constitute the core genome must be
  in a proportion of genomes that is at least equal to this value.

    - e.g.: ``0.95``

- ``--s``, ``step`` : Number of genomes added to the cgMLST computation
  at each step.

    - e.g.: ``5``

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
from plotly.offline import plot
import plotly.graph_objects as go

try:
    from utils import (constants as ct,
                       file_operations as fo,
                       iterables_manipulation as im)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 iterables_manipulation as im)


def binarize_matrix(column):
    """Convert a Pandas dataframe column values into numeric values.

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
    coln = pd.to_numeric(column)

    return np.int64(coln > 0)


def remove_genomes(matrix, genomesToRemove):
    """Remove rows from an allele calling matrix.

    Remove rows from a Pandas dataframe if the
    index identifier matches the identifier of
    a genome to remove.

    Parameters
    ----------
    matrix : pandas.core.frame.DataFrame
        Pandas dataframe with allelic profiles.
        Each row has the allelic profile of a genome
        and each column has the allele identifiers
        determined for a locus.
    genomesToRemove : list
        List of genomes to remove.

    Returns
    -------
    pruned_matrix : pandas.core.frame.DataFrame
        Input dataframe without the rows whose
        index matched an identifier of a genome
        to remove.
    """
    # determine row indexes that match any genome to remove
    to_remove_bool = matrix.index.isin(genomesToRemove)
    # create new matrix without rows that matched any genome to remove
    pruned_matrix = matrix.loc[~ to_remove_bool]

    return pruned_matrix


def remove_genes(matrix, genesToRemove):
    """Remove columns from an allele calling matrix.

    Parameters
    ----------
    matrix : pandas.core.frame.DataFrame
        Pandas dataframe with allelic profiles.
        Each row has the allelic profile of a genome
        and each column has the allele identifiers
        determined for a locus.
    genesToRemove : list
        List of loci to remove.
    """
    matrix = matrix[matrix.columns[~matrix.columns.isin(genesToRemove)]]

    return matrix


def above_threshold(column, column_length, threshold):
    """Determine if gene presence is equal or above a threshold.

    Parameters
    ----------
    column : pandas.core.series.Series
        Pandas dataframe column with presence (1) and absence (0)
        values for a locus in a set of genomes.
    column_length : int
        Number of genomes in the dataset.
    threshold : float
        Core genome determination threshold.

    Returns
    -------
    bool
        True if gene is equal or above threshold, False otherwise.
    """
    return (np.sum(column) / column_length) >= threshold


def compute_cgMLST(matrix, sorted_genomes, threshold, step):
    """Prune allele call results based on presence/absence threshold.

    Removes columns from an allele calling matrix based on
    threshold for loci presence/absence.

    Parameters
    ----------
    matrix : pandas.core.frame.DataFrame
        Pandas dataframe with allelic profiles.
        Each row has the allelic profile of a genome
        and each column has the allele identifiers
        determined for a single gene.
    sorted_genomes : list
        List of genome identifiers sorted in order of
        decresing number of loci.
    threshold : float
        Core genome determination threshold.
    step : int
        Number of genomes added to the cgMLST computation at
        each step.

    Returns
    -------
    pruned_df : pandas.core.frame.DataFrame
        Dataframe with the cgMLST profiles for the last step value
        (the last step includes all genomes in the input matrix).
    cgMLST_size : dict
        Dictionary with the number of genomes used to compute the
        cgMLST as keys and the size of the core-genome as values.
    """
    # determine genes at or above threshold
    cgMLST_size = {}
    for i in im.inclusive_range(1, len(sorted_genomes), step):
        # get subdataframe for current genomes
        current_df = matrix.loc[sorted_genomes[:i]]
        pa_rows, _ = current_df.shape
        is_above_threshold = current_df.apply(above_threshold,
                                              args=(pa_rows, threshold,))
        above = current_df.columns[is_above_threshold]
        cgMLST_size[pa_rows] = len(above)
        print('\r', 'Computed for...{0} genomes.'.format(i), end='')

    # return list of genes in cgMLST and cgMLST count per genome threshold
    return [above, cgMLST_size]


def presAbs(matrix, output_directory):
    """Create a presence/absence matrix.

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
    pa_path : str
        Path to the output TSV file that contains the
        presence-absence matrix.
    """
    presence_absence = matrix.apply(binarize_matrix)

    pa_path = fo.join_paths(output_directory, [ct.PRESENCE_ABSENCE_BASENAME])
    presence_absence.to_csv(pa_path, sep='\t')

    return [presence_absence, pa_path]


def missing_data_table(presence_absence):
    """Determine missing data per genome.

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


def main(input_file, output_directory, threshold, step,
         genes2remove, genomes2remove):
    """Determine the cgMLST based on allele calling results.

    Parameters
    ----------
    input_file : str
        Path a TSV file with allelic profiles for a set
        of genomes.
    output_directory : str
        Path to the directory where the process will
        store output files.
    threshold : list
        Core genome determination thresholds.
    step : int
        Number of genomes added to the cgMLST computation at
        each step.
    genes2remove : str
        Path to TXT file with the list of genomes to remove.
    genomes2remove : str
        Path to TXT file with the list of loci to remove.

    Returns
    -------
    List with the paths to three files:

    - Path a TSV file with the cgMLST matrix.
    - Path to a TXT file with the list of genes that
      constitute the core genome.
    - Path to a TSV file with the information about
      missing data per genome.
    """
    fo.create_directory(output_directory)

    # import matrix with allelic profiles
    matrix = pd.read_csv(input_file, header=0, index_col=0,
                         sep='\t', low_memory=False)

    total_genomes, total_loci = matrix.shape
    print('Input matrix has {0} profiles for {1} '
          'loci.'.format(total_genomes, total_loci))

    # read list of genes and list of genomes to exclude
    genesToRemove = []
    if genes2remove:
        genesToRemove = fo.read_lines(genes2remove)
    print('{0} loci to exclude.'.format(len(genesToRemove)))
    genomesToRemove = []
    if genomes2remove:
        genomesToRemove = fo.read_lines(genomes2remove)
    print('{0} genomes to exclude.'.format(len(genomesToRemove)))

    cgMLST_thresholds = threshold

    # remove genomes
    if len(genomesToRemove) > 0:
        matrix = remove_genomes(matrix, genomesToRemove)
        total_genomes -= len(genomesToRemove)

    # remove genes
    if len(genesToRemove) > 0:
        matrix = remove_genes(matrix, genesToRemove)
        total_loci -= len(genesToRemove)

    # mask missing data
    print('Masking matrix with {0} profiles for {1} '
          'loci...'.format(total_genomes, total_loci))
    masked_matrix = matrix.apply(im.replace_chars)
    print('Masked {0} profiles.'.format(total_genomes))

    # build presence/absence matrix
    print('Building presence-absence matrix...')
    pa_matrix, pa_outfile = presAbs(masked_matrix, output_directory)
    print('Presence-absence matrix saved to {0}'.format(pa_outfile))

    # count number of missing data per genome
    print('Determining missing data per genome...')
    missing_data_df = missing_data_table(pa_matrix)
    missing_stats = missing_data_df['missing'].describe()

    # sort genomes based on order of decreasing missing data
    missing_data_df = missing_data_df.sort_values('missing', ascending=True)
    sorted_genomes = missing_data_df.index.tolist()

    # write table with missing data stats
    mdata_path = os.path.join(output_directory, ct.MISSING_LOCI_BASENAME)
    missing_data_df.to_csv(mdata_path, sep='\t', index=False)
    print('Missing data table saved to {0}'.format(mdata_path))

    # compute cgMLST
    line_traces = []
    for i, t in enumerate(cgMLST_thresholds):
        print('Determining cgMLST for loci presence threshold of {0}...'.format(t))
        cgMLST_genes, cgMLST_counts = compute_cgMLST(pa_matrix,
                                                     sorted_genomes,
                                                     t,
                                                     step)

        retained = len(cgMLST_genes)
        print('\ncgMLST for loci presence threshold of {0} composed '
              'of {1}/{2} genes.'.format(t, retained, total_loci))

        # write cgMLST matrix
        # get subset from masked matrix
        cgMLST_matrix = masked_matrix[cgMLST_genes]
        cgmlst_path = os.path.join(output_directory, 'cgMLST{0}.tsv'.format(int(t*100)))
        cgMLST_matrix.to_csv(cgmlst_path, sep='\t')
        print('cgMLST profiles saved to {0}'.format(cgmlst_path))

        # write genes in cgMLST to file
        loci_path = os.path.join(output_directory, 'cgMLSTschema{0}.txt'.format(int(t*100)))
        fo.write_lines(list(cgMLST_genes), loci_path)
        print('List of loci in the cgMLST saved to {0}'.format(loci_path))

        # create line trace
        cgMLST_trace = go.Scattergl(x=list(cgMLST_counts.keys()),
                                    y=list(cgMLST_counts.values()),
                                    mode='lines',
                                    name='cgMLST{0}'.format(int(t*100)),
                                    hovertemplate=('%{y}'))
        line_traces.append(cgMLST_trace)

    # create trace with present loci
    present = [total_loci-i for i in missing_data_df['missing']]
    genomes_index = list(range(1, len(present)+1))
    miss_trace = go.Scattergl(x=genomes_index,
                              y=present,
                              mode='lines',
                              name='Present in genome',
                              line=dict(color='#000000'),
                              hovertemplate=('%{y}<br>'
                                             'Genome: %{text}<br>'),
                              text=sorted_genomes)
    line_traces.append(miss_trace)

    fig = go.Figure(data=line_traces)
    fig.update_layout(title={'text': 'Number of loci in cgMLST',
                             'font_size': 30},
                      xaxis_title='Number of genomes',
                      yaxis_title='Number of loci',
                      template='simple_white',
                      hovermode='x')
    fig.update_xaxes(range=[0, len(sorted_genomes)],
                     tickfont=dict(size=18),
                     titlefont=dict(size=20),
                     showgrid=True)
    fig.update_yaxes(tickfont=dict(size=18),
                     titlefont=dict(size=20),
                     showgrid=True)
    output_html = os.path.join(output_directory, 'cgMLST.html')
    plot(fig, filename=output_html, auto_open=False)
    print('HTML file with cgMLST per loci presence threshold '
          'and per step saved to {0}'.format(output_html))


if __name__ == '__main__':

    main()
