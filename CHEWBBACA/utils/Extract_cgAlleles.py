#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd


def presAbs(calls, listgenomesRemove, outputfile, cgPercent):

    def binarize(column):

        fix_inf = column.replace(to_replace='INF-', value='', regex=True)
        col = fix_inf.replace(to_replace='\D+.*', value='0', regex=True)
        coln = pd.to_numeric(col)

        return np.int64(coln > 0)

    def above_threshold(column, column_length):

        return (np.sum(column) / column_length) >= cgPercent

    # remove genomes
    for genome in calls.index:

        if genome in listgenomesRemove:

            print("removed genome : {}".format(genome))

    to_remove_bool = calls.index.isin(listgenomesRemove)

    calls_removed = calls.loc[~ to_remove_bool]

    print ("building the presence and absence matrix...")

    presence_absence = calls_removed.apply(binarize)


    pa_path = os.path.join(outputfile, 'Presence_Absence.tsv')
    presence_absence.to_csv(pa_path, sep='\t')

    print ("presence and abscence matrix built")

    pa_rows, _ = presence_absence.shape
    is_above_threshold = presence_absence.apply(above_threshold,
                                                args=(pa_rows,))

    columns_to_remove = calls_removed.columns[~ is_above_threshold]

    return presence_absence, calls_removed, columns_to_remove

def missing_data_table(cleaned_calls):

    n_genes = len(cleaned_calls.columns)
    genes_present = cleaned_calls.apply(np.count_nonzero, axis=1)

    missing_data = {'FILE': cleaned_calls.index,
                    'number of missing data': n_genes - genes_present,
                    'percentage': 1 - (genes_present / n_genes)}

    missing_data_df = pd.DataFrame(missing_data,
                                   columns=['FILE',
                                            'number of missing data',
                                            'percentage'])

    return missing_data_df


def clean(inputfile, outputfile, toremovegenes, toremovegenomes, cgPercent):

    # open the raw file to be clean
    calls = pd.read_csv(inputfile,
                        header=0, index_col=0, sep='\t', low_memory=False)

    # get presence abscence matrix
    presence_absence, calls_removed, cols_to_remove = presAbs(calls,
                                                              toremovegenomes,
                                                              outputfile,
                                                              cgPercent)


    total_genes_to_delete = set(toremovegenes).union(set(cols_to_remove))

    # clean the original matrix,
    # using the information on the presence/abscence matrix
    print ("processing the matrix")

    to_keep = ~calls_removed.columns.isin(total_genes_to_delete)

    cleaned = calls_removed.loc[:, to_keep]

    # remove INF and other missing data tags from the profile
    words = ['LNF', 'PLOT3', 'PLOT5', 'ASM', 'ALM', 'NIPHEM', 'NIPH', 'LOTSC']
    for word in words:
        cleaned = cleaned.replace(to_replace='{}.*'.format(word),
                                  value='0', regex=True)

    cleaned = cleaned.replace(to_replace='INF-', value='', regex=True)
    cleaned = cleaned.apply(pd.to_numeric)
    # cleaned = cleaned.reindex(cleaned.index.rename(['FILE']))

    # count number of missing data per genome
    missing_data_df = missing_data_table(cleaned)

    # write the output file
    cgmlst_path = os.path.join(outputfile, 'cgMLST.tsv')
    cleaned.to_csv(cgmlst_path, sep='\t')

    schema_path = os.path.join(outputfile, 'cgMLSTschema.txt')
    pd.Series(list(cleaned.index)).to_csv(schema_path, index=False)

    mdata_path = os.path.join(outputfile, 'mdata_stats.tsv')
    missing_data_df.to_csv(mdata_path, sep='\t', index=False)

    cgmlst_genes = len(cleaned.columns)
    totaldeletedgenes = len(calls.columns) - cgmlst_genes

    print ("deleted : {} loci".format(totaldeletedgenes))
    print ("total loci remaining : {}".format(cgmlst_genes))


def main(pathOutputfile,newfile,percent,genesToRemoveFile,genomesToRemoveFile):

    if not os.path.exists(newfile):
        os.makedirs(newfile)

    genesToRemove = []
    genomesToRemove = []

    if genesToRemoveFile:

        genesToRemove = pd.read_csv(genesToRemoveFile,
                                    header=None, index_col=None)[0]

    if genomesToRemoveFile:

        genomesToRemove = pd.read_csv(genomesToRemoveFile,
                                      header=None, index_col=None)[0]

    clean(pathOutputfile, newfile, genesToRemove, genomesToRemove, percent)


if __name__ == "__main__":
    main()
