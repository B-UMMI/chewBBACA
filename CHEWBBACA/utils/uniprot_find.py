#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO
import os
import argparse
from SPARQLWrapper import SPARQLWrapper, JSON
import csv
from collections import defaultdict
import multiprocessing
virtuoso_server=SPARQLWrapper('http://sparql.uniprot.org/sparql')

class Result():
    def __init__(self):
        self.val = []

    def update_result(self, val):
        lala=val
        self.val.append(lala)

    def get_result(self):
        return self.val

def translateSeq(DNASeq, verbose):
    if verbose:
        def verboseprint(*args):
            for arg in args:
                print (arg),
            print
    else:
        verboseprint = lambda *a: None  # do-nothing function

    seq = DNASeq
    tableid = 11
    inverted = False

    myseq = Seq(seq)
    protseq = Seq.translate(myseq, table=tableid, cds=True)


    return str(protseq)

def check_if_list_or_folder(folder_or_list):
    list_files = []
    # check if given a list of genomes paths or a folder to create schema
    try:
        f = open(folder_or_list, 'r')
        f.close()
        list_files = folder_or_list
    except IOError:

        for gene in os.listdir(folder_or_list):
            try:
                genepath = os.path.join(folder_or_list, gene)
                for allele in SeqIO.parse(genepath, "fasta"):
                    break
                list_files.append(os.path.abspath(genepath))
            except Exception as e:
                print (e)
                pass

    return list_files


def get_data(sparql_query):
    virtuoso_server.setQuery(sparql_query)
    virtuoso_server.setReturnFormat(JSON)
    virtuoso_server.setTimeout(10)
    try:
        result = virtuoso_server.query().convert()
    except:
        print ("A request to uniprot timed out, trying new request")
        time.sleep(5)
        result = virtuoso_server.query().convert()
    return result

def get_protein_info(proteinSequence):
    proteinSequence=proteinSequence.replace("*","")

    name=''
    url=''
    prevName=''
    prevUrl=''

    query='PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>  PREFIX up: <http://purl.uniprot.org/core/> select ?seq ?fname ?fname2 ?fname3  where {{?b a up:Simple_Sequence; rdf:value "'+proteinSequence+'". ?seq up:sequence ?b. OPTIONAL{?seq up:submittedName ?sname. ?sname up:fullName ?fname2} OPTIONAL{?seq up:recommendedName ?rname.?rname up:fullName ?fname} }UNION{?seq a up:Sequence; rdf:value "'+proteinSequence+'"; rdfs:label ?fname3. }}'
    result = get_data(query)
    #~ print (query)

    try:
        result["results"]["bindings"][0]
        aux=result["results"]["bindings"]
        for elem in aux:
            if 'fname' in elem.keys():
                name=str(elem['fname']['value'])
                url=str(elem['seq']['value'])
            elif 'fname2' in elem.keys():
                name=str(elem['fname2']['value'])
                url=str(elem['seq']['value'])
            elif 'fname3' in elem.keys():
                name=str(elem['fname3']['value'])
                url=str(elem['seq']['value'])


            if not "Uncharacterized protein" in name:
                break

            if prevName=='' and (not "Uncharacterized protein" in name or not "hypothetical" in name or not "DUF" in name):
                prevName=name
                prevUrl=url
            else:
                name=prevName
                url=prevUrl

            #~ print (name)
            #~ print (url)

    except Exception as e:
        return False




    return str(name),str(url)

def proc_gene(gene,auxBar):

    #~ print gene
    name=''
    url=''
    prevName=''
    prevUrl=''
    for allele in SeqIO.parse(gene, "fasta"):
        params = {}
        sequence=str(allele.seq)
        try:
            proteinSequence=translateSeq(sequence,False)
        except:
            continue
        try:
            name,url=get_protein_info(proteinSequence)
            if "Uncharacterized protein" in name or "hypothetical" in name or "DUF" in name :
                if not prevName=="":
                    name=prevName
                    url=prevUrl
                #~ print("trying next allele")
                continue
            else:
                prevName=name
                prevUrl=url
                #~ print (name)
                #~ print (url)
                break
        except Exception as e:
            #~ print (e)
            #~ print("trying next allele")
            continue

    if gene in auxBar:
        auxlen=len(auxBar)
        index=auxBar.index(gene)
        print ( "["+"="*index+">"+" "*(auxlen-index)+"] Querying "+str(int((float(index)/auxlen)*100))+"%")

    return [gene, name, url]

def main(geneFiles,proteinid2genome,cpu2use):

    geneFiles = check_if_list_or_folder(geneFiles)
    if isinstance(geneFiles, list):
        with open("listGenes.txt", "w") as f:
            for genome in geneFiles:
                f.write(genome + "\n")
        geneFiles = "listGenes.txt"

    listGenes=[]
    gene_fp = open( geneFiles, 'r')
    for gene in gene_fp:
        gene = gene.rstrip('\n')
        listGenes.append(gene)
    gene_fp.close()

    try:
        os.remove("listGenes.txt")
    except:
        pass

    dictaux= defaultdict(list)

    with open(proteinid2genome) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        headers = next(reader)
        for row in reader:
            for elem in row:
                dictaux[row[0]+row[-1]].append(str(elem))



    print ("Processing the fastas")
    got=0
    notgot=[]
    uncharacterized=[]
    selected_prots=[]
    counter=0


    #test down bar
    auxBar=[]
    #~ orderedkeys=sorted(newDict.keys())
    step=(int((len(listGenes))/10))+1
    counter=0
    while counter < len(listGenes):
        auxBar.append(listGenes[counter])
        counter+=step

    pool = multiprocessing.Pool(cpu2use)
    result = Result()
    for gene in listGenes:
        p=pool.apply_async(proc_gene,args=[gene,auxBar],callback=result.update_result)

    pool.close()
    pool.join()
    listResults=result.get_result()


    for result in listResults:
        #~ lala=process_locus(gene,path2schema,newDict,auxBar)
        gene=result[0]
        name=result[1]
        url=result[2]

        print ("final name: "+name)
        if "Uncharacterized protein" in name or "hypothetical" in name:
            uncharacterized.append(gene)
            got+=1

        elif name=="":
            notgot.append(gene)
        else:
            got+=1

        aux=gene.split("-protein")
        protid=aux[-1].replace(".fasta","")
        aux2=aux[0].split("/")[-1].replace("-","_")
        dictaux[aux2+protid].append(str(name))
        dictaux[aux2+protid].append(str(url))
        dictaux[aux2+protid].insert(0,os.path.basename(gene))
        selected_prots.append(aux2+protid)


    newProfileStr="Fasta\t"+('\t'.join(map(str, headers)))+"\tname\turl\n"
    for key in set(sorted(selected_prots,key=str)):
        newProfileStr += ('\t'.join(map(str, dictaux[key])))+"\n"

    print ("Found : " +str(got)+ ", "+str(len(uncharacterized))+" of them uncharacterized/hypothetical. Nothing found on :"+str(len(notgot)))
    with open("new_protids.tsv", "w") as f:
        f.write(newProfileStr)

    print ("Done")


if __name__ == "__main__":
    main()
