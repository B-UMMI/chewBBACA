#!/usr/bin/env python3
import os
from Bio import SeqIO
import argparse
import plotly
import plotly.graph_objs as go

import math

def calculateN_50(input_file,ratio):

    lengths=[]
    sums=0

    for seq_record in SeqIO.parse(input_file, "fasta"):
        lengths.append(len(seq_record.seq))

    lengths=sorted(lengths, reverse=True)

    contigsMore1000=0
    totalLengthMore1000=0

    for elem in lengths:
        if elem >=10000:
            contigsMore1000+=1
            totalLengthMore1000+=int(elem)
    N_half=sum(lengths)*ratio

    for i in range(0, len(lengths)):
        sums=sums+lengths[i]
        if sums >= N_half:
            return lengths[i],len(lengths),N_half/ratio,contigsMore1000,totalLengthMore1000


def addAnnotations(aux1,aux2,annotationsList,listFiles2Check):
    auxDict={}
    for k in aux1:
        if k in listFiles2Check:
            x=aux1.index(k)
            y=aux2[x]
            aux=dict(
                y=math.log(y,10),
                x=x,
                xref='x',
                yref='y',
                #~ text=k,
                showarrow=True,
                ax=20,
                ay=-40,

            )
            annotationsList.append(aux)
            auxDict[k]=[x,y]
    return auxDict

def main():

    parser = argparse.ArgumentParser(description="do a plot of number of pairs of isolates vs number of allelic differences")
    parser.add_argument('-i', nargs='?', type=str, help='path where genomes are', required=True)
    parser.add_argument('-g', nargs='?', type=str, help='List of genomes to annotate', required=False)


    args = parser.parse_args()
    contigsfile = args.i
    try:
        listGenomes = args.g
    except:
        listGenomes=[]

    results=[]
    dictContigs={}
    dictContigs1000={}
    dictBases={}
    dictBases1000={}
    dictn50={}
    dictn90={}

    try:
        f=open( contigsfile, 'r')
        f.close()
    except IOError:
        listbasename=os.path.basename(os.path.normpath(contigsfile))
        for assembly in os.listdir(contigsfile):
            print("processing: "+assembly)
            try:
                assemblypath=os.path.join(contigsfile,assembly)
                #~ n90,numberContigs,numberBases,contigsMore1000=calculateN_50(assemblypath,0.9)
                n50,numberContigs,numberBases,contigsMore1000,numberBases1000=calculateN_50(assemblypath,0.5)
                results.append([assembly,numberContigs,numberBases,n50])
                dictContigs[assembly]=numberContigs
                dictContigs1000[assembly]=contigsMore1000
                dictBases[assembly]=numberBases
                dictBases1000[assembly]=numberBases1000
                dictn50[assembly]=n50
                #~ dictn90[assembly]=n90
            except:
                print(assembly+ " is not a fasta file")

    listFiles = []

    with open(listGenomes, 'r') as genomes_fp:
        for genome in genomes_fp:
            genome = genome.rstrip('\n')
            genome = genome.rstrip('\r')
            listFiles.append( genome )



    with open("AssemblyStats.tsv", "w") as f:
        f.write("File\tNo.Contigs\tNo.Bases\tN50\n")
        for elem in results:
            f.write('\t'.join(map(str,elem))+"\n")

    #~ sortedKeys90=sorted(dictn90.keys(),key=dictn90.__getitem__,reverse=True)

    if len(listFiles)==0:
        typeScatter=go.Scattergl
    else:
        typeScatter=go.Scatter

        sortedKeys50=sorted(dictn50.keys(),key=dictn50.__getitem__,reverse=True)
        sortedKeysContigs=sorted(dictContigs.keys(),key=dictContigs.__getitem__,reverse=True)
        sortedKeysContigs1000=sorted(dictContigs1000.keys(),key=dictContigs1000.__getitem__,reverse=True)
        sortedKeysBases=sorted(dictBases.keys(),key=dictBases.__getitem__,reverse=True)
        sortedKeysBases1000=sorted(dictBases1000.keys(),key=dictBases1000.__getitem__,reverse=True)
        #~ sortedvalues90=sorted(dictn90.values(),reverse=True)
        sortedvalues50=sorted(dictn50.values(),reverse=True)
        sortedvaluesContigs=sorted(dictContigs.values(),reverse=True)
        sortedvaluesContigs1000=sorted(dictContigs1000.values(),reverse=True)
        sortedvaluesBases=sorted(dictBases.values(),reverse=True)
        sortedvaluesBases1000=sorted(dictBases1000.values(),reverse=True)


    annotations=[]

    annotDict=addAnnotations(sortedKeys50,sortedvalues50,annotations,listFiles)
    annotDict2=addAnnotations(sortedKeysContigs,sortedvaluesContigs,annotations,listFiles)
    #~ addAnnotations(sortedKeysBases,sortedvaluesContigs,annotations,listFiles)
    annotDict3=addAnnotations(sortedKeysBases1000,sortedvaluesBases1000,annotations,listFiles)
    annotDict4=addAnnotations(sortedKeysContigs1000,sortedvaluesContigs1000,annotations,listFiles)

    listTraces=[]

    for k,v in annotDict2.items():
        x1=v[0]
        y1=v[1]
        v=annotDict4[k]
        x2=v[0]
        y2=v[1]

        trace = typeScatter(
            x= [x1,x2],
            y = [y1,y2],
            name= k,
            mode = 'lines',
            hoverinfo= "none",
            showlegend=False,
            line = dict(
                color = 'grey',
                width = 1,
                dash = 'dash')
        )
        listTraces.append(trace)

    for k,v in annotDict.items():
        x1=v[0]
        y1=v[1]
        v=annotDict3[k]
        x2=v[0]
        y2=v[1]

        trace = typeScatter(
            x= [x1,x2],
            y = [y1,y2],
            name= k,
            mode = 'lines',
            hoverinfo= "none",
            showlegend=False,
            line = dict(
                color = 'grey',
                width = 1,
                dash = 'dash')
        )
        listTraces.append(trace)

        #~ for k,v in annotDict.iteritems():
        #~ x1=v[0]
        #~ y1=v[1]
        #~ v=annotDict2[k]
        #~ x2=v[0]
        #~ y2=v[1]
        #~
        #~ trace = go.Scatter(
        #~ x= [x1,x2],
        #~ y = [y1,y2],
        #~ name= k,
        #~ mode = 'lines',
        #~ hoverinfo= "none",
        #~ showlegend=False,
        #~ line = dict(
        #~ color = 'grey',
        #~ width = 1,
        #~ dash = 'dash')
        #~ )
        #~ listTraces.append(trace)


    trace1 = typeScatter(
        y = sorted(dictContigs.values(),reverse=True),
        name= "No. contigs",
        mode = 'markers',
        text = sorted(dictContigs.keys(),key=dictContigs.__getitem__,reverse=True),
        hoverinfo= "y+text",
        marker = dict(
            color = 'pink',
        ),
    )
    listTraces.append(trace1)
    #~ trace2 = go.Scatter(
    #~ y = sorted(dictBases.values(),reverse=True),
    #~ name= "No. bases",
    #~ mode = 'markers',
    #~ text = sorted(dictBases.keys(),key=dictBases.__getitem__,reverse=True),
    #~ )
    #~ trace1 = go.Scatter(
    #~ y = sorted(dictn90.values(),reverse=True),
    #~ name= "N90",
    #~ text = sorted(dictn90.keys(),key=dictn90.__getitem__,reverse=True),
    #~ )
    trace3 = typeScatter(
        y = sorted(dictn50.values(),reverse=True),
        name= "N50",
        mode = 'markers',
        marker = dict(
            color = 'red',
        ),
        text = sorted(dictn50.keys(),key=dictn50.__getitem__,reverse=True),
        hoverinfo= "y+text",
    )
    listTraces.append(trace3)
    trace4 = typeScatter(
        y = sorted(dictContigs1000.values(),reverse=True),
        name= "No. contigs > 10k bp",
        mode = 'markers',
        marker = dict(
            color = 'green',
        ),
        text = sorted(dictContigs1000.keys(),key=dictContigs1000.__getitem__,reverse=True),
        hoverinfo= "y+text",
    )
    listTraces.append(trace4)
    trace5 = typeScatter(
        y = sorted(dictBases1000.values(),reverse=True),
        name= "No. bases >10k bp",
        mode = 'markers',
        hoverinfo= "y+text",
        marker = dict(
            color = 'blue',
        ),
        text = sorted(dictBases1000.keys(),key=dictBases1000.__getitem__,reverse=True),
    )
    listTraces.append(trace5)




    layout = go.Layout(
        hovermode= 'closest',

        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showline=False,
            showticklabels=False
        ),
        yaxis=dict(
            type='log',
            autorange=True
        ),
        annotations=annotations
    )


    fig = go.Figure(data=listTraces, layout=layout)
    plot_url = plotly.offline.plot(fig, filename='AssemblyStatsStack.html')


if __name__ == "__main__":
    main()
