#!/usr/bin/env python3


import os
import csv

from Bio import SeqIO

try:
    from utils import (
        file_operations as fo,
        uniprot_requests as ur,
        parameters_validation as pv,
        sequence_manipulation as sm,
        multiprocessing_operations as mo
    )
except:
    from CHEWBBACA.utils import (
        file_operations as fo,
        uniprot_requests as ur,
        parameters_validation as pv,
        sequence_manipulation as sm,
        multiprocessing_operations as mo
    )


def get_protein_info(proteinSequence):

    proteinSequence = proteinSequence.replace("*", "")

    name = ''
    url = ''
    prevName = ''
    prevUrl = ''

    query = ('PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>  '
             'PREFIX up: <http://purl.uniprot.org/core/> '
             'select ?seq ?fname ?fname2 ?fname3  where {'
             '{?b a up:Simple_Sequence; rdf:value '
             '"'+proteinSequence+'". ?seq up:sequence ?b. '
             'OPTIONAL{?seq up:submittedName ?sname. ?sname up:fullName ?fname2} '
             'OPTIONAL{?seq up:recommendedName ?rname.?rname up:fullName ?fname} }'
             'UNION{?seq a up:Sequence; rdf:value "'+proteinSequence+'"; '
             'rdfs:label ?fname3. }}')

    result = ur.get_data(query)
    try:
        result["results"]["bindings"][0]
        aux = result["results"]["bindings"]
        for elem in aux:
            if 'fname' in elem.keys():
                name = str(elem['fname']['value'])
                url = str(elem['seq']['value'])
            elif 'fname2' in elem.keys():
                name = str(elem['fname2']['value'])
                url = str(elem['seq']['value'])
            elif 'fname3' in elem.keys():
                name = str(elem['fname3']['value'])
                url = str(elem['seq']['value'])

            if "Uncharacterized protein" not in name:
                break

            if prevName == '' and (not "Uncharacterized protein" in name or not "hypothetical" in name or not "DUF" in name):
                prevName = name
                prevUrl = url
            else:
                name = prevName
                url = prevUrl

    except Exception as e:
        return False

    return [name, url]


def proc_gene(gene):

    name = ''
    url = ''
    prevName = ''
    prevUrl = ''
    for allele in SeqIO.parse(gene, "fasta"):
        sequence = str(allele.seq)
        proteinSequence = str(sm.translate_sequence(sequence, table_id=11))

        info = get_protein_info(proteinSequence)
        if info is not False:
            name, url = info
            if "Uncharacterized protein" in name or "hypothetical" in name or "DUF" in name:
                if not prevName == "":
                    name = prevName
                    url = prevUrl
                continue
            else:
                prevName = name
                prevUrl = url
                break
        else:
            continue

    return [gene, name, url]


def main(input_files, protein_table, output_directory, cpu_cores):

    # define output paths
    outpath = os.path.abspath(output_directory)

    # create output directories
    # check if they exist first
    fo.create_directory(outpath)

    genes_list = os.path.join(output_directory, 'listGenes.txt')
    genes_list = pv.check_input_type(input_files, genes_list)

    with open(genes_list, 'r') as infile:
        listGenes = [l.rstrip('\n') for l in infile.readlines()]

    os.remove(genes_list)

    with open(protein_table, 'r') as infile:
        table_lines = list(csv.reader(infile, delimiter='\t'))
        header = table_lines[0]
        loci_info = {l[0]+l[-2]: l for l in table_lines[1:]}

    print('Searching for annotations...\n')

    uniprot_args = [[gene, proc_gene] for gene in listGenes]

    listResults = mo.map_async_parallelizer(uniprot_args,
                                            mo.function_helper,
                                            cpu_cores,
                                            show_progress=True)

    # create table lines
    got = 0
    notgot = []
    uncharacterized = []
    selected_prots = []
    for result in listResults:
        gene, name, url = result

        if "Uncharacterized protein" in name or "hypothetical" in name:
            uncharacterized.append(gene)
            got += 1

        elif name == "":
            notgot.append(gene)
        else:
            got += 1

        aux = gene.split("-protein")
        protid = aux[-1].replace(".fasta", "")
        aux2 = aux[0].split("/")[-1].replace("-", "_")
        loci_info[aux2+protid].append(str(name))
        loci_info[aux2+protid].append(str(url))
        loci_info[aux2+protid].insert(0, os.path.basename(gene))
        selected_prots.append(aux2+protid)

    new_header = ['locus_id'] + header + ['name', 'url']
    new_lines = [new_header]
    for key in set(sorted(selected_prots, key=str)):
        new_lines.append(loci_info[key])

    new_lines = ['\t'.join(l) for l in new_lines]
    table_text = '\n'.join(new_lines)

    print('\n\nFound: {0}, {1} of them uncharacterized/hypothetical. '
          'Nothing found on: {2} '.format(got,
                                          len(uncharacterized),
                                          len(notgot)))

    output_table = os.path.join(outpath, 'new_protids.tsv')
    with open(output_table, 'w') as outfile:
        outfile.write(table_text+'\n')

    print('The table with new information can be '
          'found at:\n{0}'.format(output_table))


if __name__ == "__main__":

    main()
