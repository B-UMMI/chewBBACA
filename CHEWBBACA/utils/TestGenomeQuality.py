#!/usr/bin/env python3


import os
from collections import Counter

import plotly
import numpy as np
import plotly.graph_objs as go

try:
    from utils import constants as ct
except:
    from CHEWBBACA.utils import constants as ct


def presAbs(d2c):

    row = 1
    while row < d2c.shape[0]:
        column = 1
        while column < d2c.shape[1]:
            try:
                aux = int(d2c[row, column])
                if aux > 0:
                    d2c[row, column] = 1
                else:
                    d2c[row, column] = 0
            except:
                try:
                    # remove 'INF-' prefix
                    aux = str((d2c[row, column])).replace("INF-", "")
                    # remove '*' added to chewie-NS schemas
                    aux = str((d2c[row, column])).replace("*", "")
                    aux = int(aux)
                    d2c[row, column] = 1
                except Exception as e:
                    d2c[row, column] = 0

            column += 1
        row += 1

    return d2c


# function to report the most problematic locus/genomes and report the number of good loci by % of genomes
# a list of genomes with a sum of problems higher than the given threshold will be returned
def presence3(d2, ythreshold, vector, abscenceMatrix):

    classifications = ct.ALLELECALL_CLASSIFICATIONS + [ct.PROBABLE_LNF]

    d2d = d2
    genomeslist = d2d[1:, :1]
    geneslist = d2d[0]

    column = 1
    totals = []
    allbadgenomes = []
    reallybadgenomes = []
    allverygood = []
    allverybad = []
    plus95 = 0
    plus99 = 0
    plus995 = 0
    listgenes2show = []

    # if a locus has one problematic call, check how many genomes have problem in that call
    while column < d2d.shape[1]:
        row = 1
        notfound = 0
        badgenomes = []
        while row < d2d.shape[0]:
            if not abscenceMatrix:
                if any([c in d2d[row][column] for c in classifications]):
                    d2d[row, column] = 0
                    notfound += 1
                    badgenomes.append(row - 1)
                else:
                    d2d[row, column] = 1
                    pass
            else:
                if int(d2d[row, column]) == 0:
                    notfound += 1
                    badgenomes.append(row - 1)

                else:
                    d2d[row, column] = 1
                    pass

            row += 1

        value = float(d2d.shape[0] - notfound) / float(d2d.shape[0])

        if len(genomeslist) > 500:
            xthreshold = 0.95
        elif len(genomeslist) < 20:
            xthreshold = 0.90
        else:
            xthreshold = 0.95

        if value > xthreshold:
            listgenes2show.append(geneslist[column])

        if value > xthreshold and value < 1:
            for badgenome in badgenomes:
                allbadgenomes.append((genomeslist[int(badgenome)])[0])
        elif value == 1:
            allverygood.append(geneslist[column])
        elif value == 0:
            allverybad.append(geneslist[column])
        if value >= 0.95:
            plus95 += 1
        if value >= 0.99:
            plus99 += 1
        if value >= 0.995:
            plus995 += 1

        totals.append(float(float(d2d.shape[0] - notfound) / float(d2d.shape[0])))
        column += 1

    counter = Counter(allbadgenomes)

    for elem in counter.most_common():
        if int(elem[1]) > ythreshold:
            reallybadgenomes.append((elem)[0])

    d2d = d2d.T

    # number of used genomes
    vector[0].append(len(genomeslist))

    # number of loci at 95%
    vector[1].append(plus95)

    # number of loci at 99%
    vector[2].append(plus99)

    # number of loci at 100%
    vector[6].append(len(allverygood))

    # number of loci at 0%
    vector[4].append(len(allverybad))

    # number of to be removed genomes%
    try:
        vector[5].append(len(reallybadgenomes) + ((vector[5])[-1]))
    except:
        vector[5].append(len(reallybadgenomes))

    # number of loci at 99.5%
    vector[3].append(plus995)

    return d2d, reallybadgenomes, vector, True, listgenes2show


def removegenomes(d2a, bagenomeslist):

    deleted = 0
    genomesList = (d2a[1:, :1]).tolist()

    badgenomeslist2 = []
    for elem in bagenomeslist:
        badgenomeslist2.append(elem)

    if len(bagenomeslist) > 0:
        i = 0
        for genome in genomesList:
            i += 1
            if genome[0] in badgenomeslist2:
                d2a = np.delete(d2a, i, axis=0)
                deleted += 1
                i -= 1

    return d2a


def clean(d2, iterations, ythreshold, out_folder):

    abscencematrix = True

    i = 0
    removedlistgenomes = []
    toremovegenomes = []

    statsvector = [[] for x in range(7)]

    # run a function that gives information about the usable locus and the possible bad genomes
    lastremovedgenomesCount = 0
    iterationStabilizedat = None
    isStable = False
    listgenes2show = []
    listgenes2showtotal = []

    while i <= iterations:

        for elem in toremovegenomes:
            removedlistgenomes.append(elem)

        if len(removedlistgenomes) > lastremovedgenomesCount and i > 0:
            lastremovedgenomesCount = len(removedlistgenomes)
        elif iterationStabilizedat is None and i > 0:
            iterationStabilizedat = i
            isStable = True
            listgenes2showtotal.append(listgenes2show)
        if not isStable:
            d2 = removegenomes(d2, toremovegenomes)
            matrix3, toremovegenomes, statsvector, abscencematrix, listgenes2show = presence3(d2, ythreshold,
                                                                                              statsvector,
                                                                                              abscencematrix)
        else:
            for vector in statsvector:
                vector.append(vector[-1])

        i += 1

    with open(os.path.join(out_folder, "removedGenomes.txt"), "a") as f:
        f.write((str(ythreshold) + "\t" + (' '.join(map(str, removedlistgenomes))) + "\n"))

    with open(os.path.join(out_folder, "Genes_95%.txt"), "a") as f:
        f.write(str(ythreshold) + "\t")

        for x in listgenes2showtotal:
            f.write((' '.join(map(str, x))) + "\n")

    return statsvector, iterationStabilizedat


def main(input_file, output_directory, max_iteration, max_threshold, step):

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    allresults = []
    threshold = 0
    thresholdlist = []
    listStableIter = []

    # for each threshold run a clean function on the dataset, using the previous run output (to be removed genomes) as input for the new one
    with open(os.path.join(output_directory, "removedGenomes.txt"), "w") as f:
        f.write("Threshold\tRemoved_genomes\n")
    with open(os.path.join(output_directory, "Genes_95%.txt"), "w") as f:
        f.write("Threshold\tPresent_genes\n")

    d2original = np.genfromtxt(input_file, delimiter='\t', dtype=None, encoding=None)
    d2copy = np.copy(d2original)
    # create abscence/presence matrix
    print('Creating presence/absence matrix.')
    d2copy = presAbs(d2copy)

    while threshold <= max_threshold:
        print('\r', 'Calculating for threshold {0}...'.format(threshold), end='')
        thresholdlist.append(threshold)
        result, stabilizedIter = clean(d2copy, max_iteration, threshold, output_directory)
        listStableIter.append(stabilizedIter)
        allresults.append(result)
        threshold += step

    labels = ["number of genomes",
              "Number of Loci present in 95% genomes",
              "Number of Loci present in 99% genomes",
              "Number of Loci present in 99.5% genomes",
              "Number of Loci present in 100% genomes"]

    i = max_iteration
    while i <= max_iteration:
        aux2 = []
        for resultPerThresh in allresults:
            aux = []

            for result in resultPerThresh:
                aux.append(result[i])
            aux2.append(aux)

        d2 = np.asarray(aux2)

        d2 = d2.T
        linenmbr = 0

        d2 = np.delete(d2, (4), axis=0)
        d2 = np.delete(d2, (4), axis=0)

        listtraces = []
        for line in d2:
            if (linenmbr == 0):
                trace = go.Scatter(
                    x=thresholdlist,
                    y=line,
                    name="Number of genomes used",
                    mode='lines+markers',
                    yaxis='y2',
                    marker=dict(symbol='diamond-dot', size=10)
                )
            else:
                trace = go.Scatter(
                    x=thresholdlist,
                    name=labels[linenmbr],
                    mode='lines+markers',
                    y=line,
                    marker=dict(symbol='star-dot', size=10),
                    line=dict(dash='dash')
                )

            listtraces.append(trace)
            linenmbr += 1

        layout = go.Layout(
            title='Test genomes quality',
            xaxis=dict(
                title='Threshold'
            ),
            yaxis=dict(
                title='Number of loci'
            ),
            yaxis2=dict(
                title='Number of genomes',
                overlaying='y',
                side='right'
            )
        )

        fig = go.Figure({"data": listtraces, "layout": layout})
        plotly.offline.plot(fig, filename=os.path.join(output_directory, 'GenomeQualityPlot.html'))

        i += 1


if __name__ == "__main__":
    main()
