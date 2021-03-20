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
import shutil
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


def species_lines(uniprot_readme, taxa):
    """
    """

    proteomes = []
    for line in uniprot_readme:
        if any([taxon in line for taxon in taxa]) is True:
            proteomes.append(line)

    proteomes = [line.strip('\n') for line in proteomes]
    proteomes = [line.split('\t') for line in proteomes]

    return proteomes


def get_proteomes(proteomes, output_dir):
    """
    """

    print('Downloading reference proteomes...')
    # construct FTP URLs for each proteome
    downloaded = 0
    proteomes_files = []
    for l in proteomes:
        domain = '{0}{1}'.format(l[3][0].upper(), l[3][1:])
        proteome_id = '{0}_{1}'.format(l[0], l[1])
        proteome_file = '{0}.fasta.gz'.format(proteome_id)
        local_proteome_file = os.path.join(output_dir, proteome_file)
        proteome_url = os.path.join(ct.UNIPROT_PROTEOMES_FTP, domain, l[0], proteome_file)
        proteomes_files.append(local_proteome_file)
        download_file(proteome_url, local_proteome_file)
        downloaded += 1
        print('\r', 'Downloaded {0}/{1}'.format(downloaded, len(proteomes)), end='')

    return proteomes_files


def unzip_file(compressed_file, archive_type='.gz'):
    """
    """

    lines = []
    with gzip.open(compressed_file, 'rb') as f:
        for line in f:
            lines.append(line.decode())

    # save uncompressed contents
    uncompressed_file = compressed_file.rstrip('.gz')
    fo.write_lines(lines, uncompressed_file, joiner='')

    return uncompressed_file


def translate_fastas(fasta_paths, output_directory):
    """
    """

    protein_files = []
    for path in fasta_paths:
        records = fao.import_sequences(path)
        translated_records = {seqid: str(sm.translate_dna(seq, 11, 0)[0][0])
                              for seqid, seq in records.items()}
        translated_lines = fao.fasta_lines(list(translated_records.keys()),
                                           translated_records)

        basename = fo.file_basename(path).replace('.fasta', '_protein.fasta')
        prot_file = fo.join_paths(output_directory, [basename])

        fo.write_lines(translated_lines, prot_file)
        protein_files.append(prot_file)

    return protein_files


def get_self_scores(fasta_file, temp_directory, threads):
    """
    """

    basename = fo.file_basename(fasta_file, suffix=False)

    integer_seqids = fo.join_paths(temp_directory,
                                   ['{0}_int.fasta'.format(basename)])
    ids_dict = fao.integer_headers(fasta_file, integer_seqids)

    blastdb = fo.join_paths(temp_directory, ['{0}_db'.format(basename)])
    stderr = bw.make_blast_db('makeblastdb', integer_seqids,
                              blastdb, 'prot')

    blastout = fo.join_paths(temp_directory, ['self_blastout.tsv'])
    self_results = bw.run_blast('blastp', blastdb, integer_seqids,
                                blastout, threads=threads, max_targets=1)

    self_lines = fo.read_tabular(blastout)
    self_lines_ids = {ids_dict[l[0]]: l[-1] for l in self_lines}

    return self_lines_ids


def proteome_annotations(input_files, temp_directory, taxa, blast_score_ratio,
                         cpu_cores, proteome_matches):
    """
    """

    # get paths to files with representative sequences
    short_directory = fo.join_paths(input_files, ['short'])
    reps_paths = [fo.join_paths(short_directory, [file])
                  for file in os.listdir(short_directory)
                  if file.endswith('.fasta') is True]

    print('Translating representative sequences...')
    # translate representatives for all loci
    translated_reps = fo.join_paths(temp_directory, ['translated_reps'])
    fo.create_directory(translated_reps)

    reps_protein_files = translate_fastas(reps_paths, translated_reps)

    print('Downloading list of reference proteomes...')
    remote_readme = os.path.join(ct.UNIPROT_PROTEOMES_FTP, 'README')
    local_readme = os.path.join(temp_directory, 'reference_proteomes_readme.txt')

    # get README file with list of reference proteomes
    download_file(remote_readme, local_readme)

    # get lines with proteomes info for species of interest
    readme_lines = fo.read_lines(local_readme, strip=False)

    selected_proteomes = species_lines(readme_lines, taxa)
    print('Found {0} reference proteomes for {1}.'.format(len(selected_proteomes), taxa))
    if len(selected_proteomes) > 0:
        # create directory to store proteomes
        proteomes_directory = os.path.join(temp_directory, 'proteomes')
        fo.create_directory(proteomes_directory)

        proteomes_files = get_proteomes(selected_proteomes, proteomes_directory)

    # uncompress files and concatenate into single FASTA
    uncompressed_proteomes = [unzip_file(file) for file in proteomes_files]
    proteomes_concat = fo.join_paths(proteomes_directory,
                                     ['full_proteome.fasta'])
    proteomes_concat = fo.concatenate_files(uncompressed_proteomes,
                                            proteomes_concat)

    # get self-scores
    # concatenate protein files
    reps_concat = fo.concatenate_files(reps_protein_files,
                                       fo.join_paths(temp_directory, ['reps_concat.fasta']))

    print('\nDetermining self-score of representatives...')
    self_scores = get_self_scores(reps_concat, temp_directory, cpu_cores)

    # create BLASTdb with proteome sequences
    proteome_blastdb = fo.join_paths(proteomes_directory,
                                     ['proteomes_db'])
    stderr = bw.make_blast_db('makeblastdb', proteomes_concat,
                              proteome_blastdb, 'prot')

    # BLASTp to determine annotations
    blast_inputs = [['blastp', proteome_blastdb, file, file+'_blastout.tsv',
                     1, 1, None, None, proteome_matches, None, bw.run_blast] for file in reps_protein_files]

    print('\nBLASTing representatives against proteomes...')
    blast_results = mo.map_async_parallelizer(blast_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    blastout_files = [fo.join_paths(translated_reps, [file])
                      for file in os.listdir(translated_reps)
                      if 'blastout' in file]

    # index proteome file
    indexed_proteome = SeqIO.index(proteomes_concat, 'fasta')

    # process results for each BLASTp
    proteome_results = {}
    hit_info = 0
    for file in blastout_files:
        locus_id = fo.file_basename(file).split('_short')[0]
        results = fo.read_tabular(file)
        proteome_results[locus_id] = []
        if len(results) > 0:
            # compute BSR values
            for r in results:
                r.append(float(r[2])/float(self_scores[r[0]]))
            # sort based on BSR
            sorted_results = sorted(results, key=lambda x: x[-1])
            # get results equal or above BSR
            high_bsr_results = [r for r in sorted_results if r[-1] >= blast_score_ratio]
            for res in high_bsr_results[0:proteome_matches]:
                proteome_results[locus_id].append(res)
                # get record and extract relevant info
                hit = indexed_proteome[res[1]]
                hit_dict = vars(hit)
                # some tags might be missing
                hit_id = hit_dict.get('id', '')
                hit_description = hit_dict.get('description', '')
                if hit_description != '':
                    if 'OS=' in hit_description:
                        # get organism name
                        hit_species = (hit_description.split('OS=')[1]).split(' OX=')[0]
                        hit_product = (hit_description.split(hit_id+' ')[1]).split(' OS=')[0]
                    else:
                        hit_product = ''
                    if 'GN=' in hit_description:
                        hit_gene_name = (hit_description.split('GN=')[1]).split(' PE=')[0]
                    else:
                        hit_gene_name = ''
                proteome_results[locus_id][-1].extend([hit_id, hit_product, hit_gene_name, hit_species])
        # sort from highest to lowest BSR
        proteome_results[locus_id] = sorted(proteome_results[locus_id], key=lambda x: x[3], reverse=True)

    return proteome_results


def sparql_annotations(listGenes, cpu_cores):
    """
    """

    # create inputs to multiprocessing
    uniprot_args = [[gene, proc_gene] for gene in listGenes]

    # this should work with all alleles in the loci to maximize chance of finding info
    workers = cpu_cores if cpu_cores <= 4 else 4
    listResults = mo.map_async_parallelizer(uniprot_args,
                                            mo.function_helper,
                                            workers,
                                            show_progress=True)

    return listResults


def join_annotations(listResults, proteome_results, loci_info):
    """
    """

    selected = {}
    for result in listResults:
        gene, name, url = result

        gene_basename = fo.file_basename(gene, suffix=False)

        loci_info[gene_basename].append(str(name))
        loci_info[gene_basename].append(str(url))
        loci_info[gene_basename].insert(0, os.path.basename(gene))

        proteome_info = []
        if len(proteome_results) > 0:
            if len(proteome_results[gene_basename]) > 0:
                proteome_info = proteome_results[gene_basename]
            else:
                proteome_info = []

            loci_info[gene_basename].append(proteome_info)

        selected[gene_basename] = loci_info[gene_basename]

    return selected


def create_annotations_table(annotations, output_directory, header,
                             schema_name):
    """
    """

    new_lines = [header]
    for locus, data in annotations.items():
        new_line = [locus] + data[1:9]
        if len(data[-1]) > 0:
            relevant_data = [d[4:]+[str(round(d[3],2))] for d in data[-1]]
            proteome_data = list(zip(*relevant_data))
            proteome_data = [';'.join(list(map(str, d))) for d in proteome_data]
            new_line.extend(proteome_data)
        new_lines.append(new_line)

    new_lines = ['\t'.join(l) for l in new_lines]
    table_text = '\n'.join(new_lines)

    output_table = os.path.join(output_directory, '{0}_annotations.tsv'.format(schema_name))
    with open(output_table, 'w') as outfile:
        outfile.write(table_text+'\n')

    return output_table


input_files = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/spyogenes_schema/schema_seed'
protein_table = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/spyogenes_schema/cds_info.tsv'
output_directory = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/spyogenes_schema/annotations'
blast_score_ratio = 0.6
cpu_cores = 4
no_cleanup = True
taxa = ['Streptococcus agalactiae', 'Streptococcus pyogenes']
proteome_matches = 2
def main(input_files, protein_table, output_directory, blast_score_ratio,
         cpu_cores, no_cleanup, taxa):

    # create output directory
    fo.create_directory(output_directory)

    # create temp directory
    temp_directory = os.path.join(output_directory, 'temp')
    fo.create_directory(temp_directory)

    # validate input files
    genes_list = os.path.join(temp_directory, 'listGenes.txt')
    genes_list = pv.check_input_type(input_files, genes_list)
    listGenes = fo.read_lines(genes_list)

    print('Schema: {0}'.format(os.path.dirname(listGenes[0])))
    print('Number of loci: {0}'.format(len(listGenes)))
    schema_basename = fo.file_basename(os.path.dirname(listGenes[0]))
    # check if loci identifiers are in cds_info table

    # find annotations based on reference proteomes for species
    proteome_results = {}
    if taxa is not None:
        proteome_results = proteome_annotations(input_files, temp_directory,
                                                taxa, blast_score_ratio,
                                                cpu_cores, proteome_matches)

    # find annotations in SPARQL endpoint
    listResults = sparql_annotations(listGenes, cpu_cores)

    # read "cds_info.tsv" file created by CreateSchema
    table_lines = fo.read_tabular(protein_table)
    loci_info = {}
    for l in table_lines[1:]:
        locus_id = l[0].replace('_', '-')
        locus_id = locus_id + '-protein{0}'.format(l[-2])
        loci_info[locus_id] = l

    annotations = join_annotations(listResults, proteome_results, loci_info)

    # table header
    header = ['Locus_ID'] + table_lines[0] + ['Uniprot_Name', 'UniProt_URL']
    if len(proteome_results) > 0:
        header.extend(['Proteome_ID', 'Proteome_Product',
                       'Proteome_Gene_Name', 'Proteome_Species',
                       'Proteome_BSR'])

    output_table = create_annotations_table(annotations, output_directory,
                                            header, schema_basename)

    print('\n\nFound: {0}, {1} of them uncharacterized/hypothetical. '
          'Nothing found on: {2} '.format(got,
                                          len(uncharacterized),
                                          len(notgot)))

    if no_cleanup is False:
        shutil.rmtree(temp_directory)

    print('The table with new information can be '
          'found at:\n{0}'.format(output_table))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', type=str,
                        required=True, dest='input_files',
                        help='Path to the schema\'s directory or to a file '
                             'with a list of paths to loci FASTA files, one '
                             'per line.')

    parser.add_argument('-t', '--protein-table', type=str,
                        required=True, dest='protein_table',
                        help='Path to the "cds_info.tsv" file created by '
                             'the CreateSchema process.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='The directory where the output files will be '
                             'saved (will create the directory if it does not '
                             'exist).')

    parser.add_argument('--bsr', type=float, required=False, dest='blast_score_ratio',
                        default=0.6, help='')

    parser.add_argument('--cpu', '--cpu-cores', type=int,
                        required=False, dest='cpu_cores',
                        default=1, help='')

    parser.add_argument('--no-cleanup', action='store_true',
                        required=False, dest='no_cleanup',
                        help='')

    parser.add_argument('--taxa', nargs='+', type=str,
                        required=False, dest='taxa',
                        help='')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
