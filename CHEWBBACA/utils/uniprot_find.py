#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module

Expected input
--------------


Code documentation
------------------
"""


import os
import gzip
import argparse
import urllib.request

from Bio import SeqIO

try:
    from utils import (
        constants as ct,
        blast_wrapper as bw,
        file_operations as fo,
        uniprot_requests as ur,
        fasta_operations as fao,
        parameters_validation as pv,
        sequence_manipulation as sm,
        iterables_manipulation as im,
        multiprocessing_operations as mo
    )
except:
    from CHEWBBACA.utils import (
        constants as ct,
        blast_wrapper as bw,
        file_operations as fo,
        uniprot_requests as ur,
        fasta_operations as fao,
        parameters_validation as pv,
        sequence_manipulation as sm,
        iterables_manipulation as im,
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


def download_file(file_url, outfile):
    """
    """

    urllib.request.urlretrieve(file_url, outfile)


def species_lines(uniprot_readme, species_name):
    """
    """

    proteomes = []
    for line in uniprot_readme:
        if species_name in line:
            proteomes.append(line)

    proteomes = [line.strip('\n') for line in proteomes]
    proteomes = [line.split('\t') for line in proteomes]

    return proteomes


def get_proteomes(proteomes, output_dir):
    """
    """

    # construct FTP URLs for each proteome
    proteomes_files = []
    for l in proteomes:
        domain = '{0}{1}'.format(l[3][0].upper(), l[3][1:])
        proteome_id = '{0}_{1}'.format(l[0], l[1])
        proteome_file = '{0}.fasta.gz'.format(proteome_id)
        local_proteome_file = os.path.join(output_dir, proteome_file)
        proteome_url = os.path.join(ct.UNIPROT_PROTEOMES_FTP, domain, l[0], proteome_file)
        proteomes_files.append(local_proteome_file)
        download_file(proteome_url, local_proteome_file)

    return proteomes_files


def unzip_file(compressed_file):
    """
    """

    lines = []
    with gzip.open(compressed_file, 'rb') as f:
        for line in f:
            lines.append(line.decode())

    return lines


def main(input_files, protein_table, output_directory, species,
         blast_score_ratio, cpu_cores):

    # create output directories
    # check if they exist first
    fo.create_directory(output_directory)

    # create temp directory
    temp_directory = os.path.join(output_directory, 'temp')
    fo.create_directory(temp_directory)

    # validate input files
    genes_list = os.path.join(temp_directory, 'listGenes.txt')
    genes_list = pv.check_input_type(input_files, genes_list)
    listGenes = fo.read_lines(genes_list)

    # get paths to short directory
    reps_paths = []
    for path in listGenes:
        basename = os.path.basename(path)
        rep_basename = basename.replace('.fasta', '_short.fasta')
        parent_dir = os.path.dirname(path)
        short_path = os.path.join(parent_dir, 'short', rep_basename)
        reps_paths.append(short_path)

    # translate representatives for all loci
    translated_reps = fo.join_paths(temp_directory, ['translated_reps'])
    fo.create_directory(translated_reps)

    reps_protein_files = []
    for path in reps_paths:
        prot_basename = os.path.basename(path).replace('.fasta', '_protein.fasta')
        records = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(path, 'fasta')]
        translated_records = [(l[0], sm.translate_dna(l[1], 11, 0)) for l in records]
        prot_records = ['>{0}\n{1}'.format(r[0], r[1][0][0]) for r in translated_records]
        new_file = fo.join_paths(translated_reps, [prot_basename])
        recs_lines = im.join_list(prot_records, '\n')
        fo.write_to_file(recs_lines, new_file, 'a', '\n')
        reps_protein_files.append(new_file)

    # read "cds_info.tsv" file created by CreateSchema
    table_lines = fo.read_tabular(protein_table)
    header = table_lines[0]
    loci_info = {}
    for l in table_lines[1:]:
        locus_id = l[0].replace('_', '-')
        locus_id = locus_id + '-protein{0}'.format(l[-2])
        loci_info[locus_id] = l

    print('\nLoci to annotate: {0}'.format(len(listGenes)))

    # find annotations based on reference proteomes for species
    found = False
    if species is not None:

        proteomes_readme = os.path.join(ct.UNIPROT_PROTEOMES_FTP, 'README')
        local_readme = os.path.join(temp_directory, 'reference_proteomes_readme.txt')

        # get README file with list of reference proteomes
        download_file(proteomes_readme, local_readme)

        # get lines with proteomes info for species of interest
        readme_lines = fo.read_lines(local_readme, strip=False)

        proteomes = species_lines(readme_lines, species)
        if len(proteomes) > 0:
            # create directory to store proteomes
            proteomes_directory = os.path.join(temp_directory, 'proteomes')
            fo.create_directory(proteomes_directory)

            proteomes_files = get_proteomes(proteomes, proteomes_directory)

        # uncompress files and join into single FASTA
        proteomes_records = []
        for file in proteomes_files:
            records = unzip_file(file)
            proteomes_records.extend(records)

        full_reference = os.path.join(proteomes_directory, 'full_proteome.fasta')
        fo.write_to_file(''.join(proteomes_records), full_reference, 'w', '')

        # get self-score
        # concatenate protein files
        self_file = fo.concatenate_files(reps_protein_files, os.path.join(temp_directory, 'self.fasta'))

        integer_seqids = os.path.join(temp_directory, 'self_int.fasta')
        ids_dict = fao.integer_headers(self_file, integer_seqids)

        self_blastdb = os.path.join(temp_directory, 'self_db')
        stderr = bw.make_blast_db('makeblastdb', integer_seqids,
                                  self_blastdb, 'prot')

        self_results = bw.run_blast('blastp', self_blastdb, integer_seqids,
                                    os.path.join(temp_directory, 'self_blastout.tsv'),
                                    threads=cpu_cores, max_targets=1)

        self_lines = fo.read_tabular(os.path.join(temp_directory, 'self_blastout.tsv'))
        self_lines_ids = {ids_dict[l[0]]: l[-1] for l in self_lines}

        # create BLASTdb with proteome sequences
        proteome_blastdb = os.path.join(proteomes_directory, 'proteomes_db')
        stderr = bw.make_blast_db('makeblastdb', full_reference,
                                  proteome_blastdb, 'prot')

        # BLASTp to determine annotations
        blast_inputs = [['blastp', proteome_blastdb, file, file+'_blastout.tsv',
                         1, 1, None, None, 1, None, bw.run_blast] for file in reps_protein_files]

        blast_results = mo.map_async_parallelizer(blast_inputs,
                                                  mo.function_helper,
                                                  cpu_cores,
                                                  show_progress=True)

        # process results for each BLASTp
        # index proteome file
        indexed_proteome = SeqIO.index(full_reference, 'fasta')
        proteome_results = {}
        blastout_files = [fo.join_paths(translated_reps, [file]) for file in os.listdir(translated_reps) if 'blastout' in file]
        for file in blastout_files:
            locus_id = os.path.basename(file).split('_short')[0]
            res = fo.read_tabular(file)
            if len(res) > 0:
                # add BSR values
                for r in res:
                    r.append(float(r[2])/float(self_lines_ids[r[0]]))
                # sort based on BSR
                sorted_res = sorted(res, key=lambda x: x[-1])
                if sorted_res[0][-1] > blast_score_ratio:
                    proteome_results[locus_id] = sorted_res[0]
                    # get record and extract relevant info
                    hit = indexed_proteome[sorted_res[0][1]]
                    hit_id = hit.id
                    hit_description = hit.description
                    hit_product = (hit_description.split(hit_id+' ')[1]).split(' OS=')[0]
                    hit_gene_name = (hit_description.split('GN=')[1]).split(' PE=')[0]
                    proteome_results[locus_id].extend([hit_id, hit_product, hit_gene_name])
                else:
                    proteome_results[locus_id] = []
            else:
                proteome_results[locus_id] = res

        found = True
############################################################################
# find annotations in SPARQL endpoint

    # create inputs to multiprocessing
    uniprot_args = [[gene, proc_gene] for gene in listGenes]

    # this should work with all alleles in the loci to maximize chance of finding info
    workers = cpu_cores if cpu_cores <= 4 else 4
    listResults = mo.map_async_parallelizer(uniprot_args,
                                            mo.function_helper,
                                            workers,
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

        gene_basename = os.path.basename(gene)
        gene_basename = gene_basename.split('.fasta')[0]

        loci_info[gene_basename].append(str(name))
        loci_info[gene_basename].append(str(url))
        loci_info[gene_basename].insert(0, os.path.basename(gene))
        selected_prots.append(gene_basename)

    new_header = ['Locus_ID'] + header + ['Uniprot_Name', 'UniProt_URL']
    if found is True:
        new_header.extend(['Proteome_ID', 'Proteome_Product', 'Proteome_Gene_Name', 'Proteome_BSR'])

    new_lines = [new_header]
    for key in set(sorted(selected_prots, key=str)):
        if found is False:
            new_lines.append(loci_info[key])
        else:
            if len(proteome_results[key]) > 0:
                proteome_info = proteome_results[key][4:]+[str(proteome_results[key][3])]
            else:
                proteome_info = ['']*4

            new_lines.append(loci_info[key]+proteome_info)

    new_lines = ['\t'.join(l) for l in new_lines]
    table_text = '\n'.join(new_lines)

    print('\n\nFound: {0}, {1} of them uncharacterized/hypothetical. '
          'Nothing found on: {2} '.format(got,
                                          len(uncharacterized),
                                          len(notgot)))

    output_table = os.path.join(output_directory, 'new_protids.tsv')
    with open(output_table, 'w') as outfile:
        outfile.write(table_text+'\n')

    print('The table with new information can be '
          'found at:\n{0}'.format(output_table))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', type=str,
                        required=True, dest='input_files',
                        help='Path to the schema\'s directory or to a file with '
                             'a list of paths to loci FASTA files, one per line.')

    parser.add_argument('-t', '--protein-table', type=str,
                        required=True, dest='protein_table',
                        help='Path to the "proteinID_Genome.tsv" file created by '
                             'the CreateSchema process.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='The directory where the output files will be '
                             'saved (will create the directory if it does not '
                             'exist).')

    parser.add_argument('-s', type=str, required=False, dest='species',
                        help='')

    parser.add_argument('--bsr', type=float, required=False, dest='blast_score_ratio',
                        default=0.6, help='')

    parser.add_argument('--cpu', '--cpu-cores', type=int,
                        required=False, dest='cpu_cores',
                        default=1, help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
