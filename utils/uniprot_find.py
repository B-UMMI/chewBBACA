#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import os
import argparse
import urllib
import datetime
from SPARQLWrapper import SPARQLWrapper, JSON
import requests
import csv
from collections import defaultdict
virtuoso_server=SPARQLWrapper('http://sparql.uniprot.org/sparql')


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
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        try:
            seq = reverseComplement(seq)
            myseq = Seq(seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
            inverted = True
        except:
            try:
                seq = seq[::-1]
                myseq = Seq(seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
                inverted = True
            except:
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                    inverted = False
                except Exception as e:
                    verboseprint("translation error")
                    verboseprint(e)
                    raise

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
                for allele in SeqIO.parse(genepath, "fasta", generic_dna):
                    break
                list_files.append(os.path.abspath(genepath))
            except Exception as e:
                print (e)
                pass

    return list_files


def get_data(sparql_query):
    virtuoso_server.setQuery(sparql_query)
    virtuoso_server.setReturnFormat(JSON)
    result = virtuoso_server.query().convert()
    return result

def get_protein_info(proteinSequence,notgot,uncharacterized,got,selected_prots):
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


def main():
    parser = argparse.ArgumentParser(
        description="This program get names for a schema")
    parser.add_argument('-i', nargs='?', type=str, help='path to folder containg the schema fasta files ( alternative a list of fasta files)', required=True)
    parser.add_argument('-t', nargs='?', type=str, help='path to tsv file)', required=True)

    args = parser.parse_args()

    geneFiles = args.i
    proteinid2genome = args.t
    
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
    for gene in listGenes:
        counter+=1
        print (gene)
        print (str(counter)+"/"+str(len(listGenes)))
        name=""
        url=""
        prevName=""
        prevUrl=""
        for allele in SeqIO.parse(gene, "fasta", generic_dna):
            params = {}
            sequence=str(allele.seq)
            try:
                proteinSequence=translateSeq(sequence,False)
            except:
                continue
            try:
                name,url=get_protein_info(proteinSequence,notgot,uncharacterized,got,selected_prots)
                if "Uncharacterized protein" in name or "hypothetical" in name or "DUF" in name :
                    if not prevName=="":
                        name=prevName
                        url=prevUrl
                    print("trying next allele")
                    continue
                else:
                    prevName=name
                    prevUrl=url
                    #~ print (name)
                    #~ print (url)
                    break
            except:
                #~ print("trying next allele")
                continue
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
