#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the creation of a TSV file with annotation
terms for the loci in a schema.

The process queries UniProt's SPARQL endpoint to find exact
matches and retrieve the product name and page URL for those
matches. If users provide a taxon/taxa name/s, the process
will also search for reference proteomes for the specified
taxon/taxa and use BLASTp to align local sequences against
reference sequences to assign annotation terms based on the
BSR value.

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_files`` : Path to the schema's directory or to a file
  with a list of paths to loci FASTA files, one per line.

    - e.g.: ``/home/user/schemas/my_schema``

- ``-t``, ``protein_table`` : Path to the "cds_info.tsv" file created by
  the CreateSchema process.

    - e.g.: ``/home/user/schemas/my_schema/cds_info.tsv``

- ``-o``, ``output_directory`` : Output directory where the process will
  store intermediate files and save the final TSV file with the annotations.

    - e.g.: ``/home/user/schemas/my_schema/annotations``

- ``--bsr``, ``blast_score_ratio`` : BLAST Score Ratio value. This
  value is only used when a taxon/taxa is provided and local sequences
  are aligned against reference proteomes.

    - e.g.: ``0.6``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores used to run the process.

    - e.g.: ``4``

- ``--taxa``, ``taxa`` : List of scientific names for a set of taxa. The
  process will search for and download reference proteomes with terms that
  match any of the provided taxa.

    - e.g.: ``"Streptococcus pyogenes"``

- ``--pm``, ``proteome_matches`` : Maximum number of proteome matches to
  report.

    - e.g.: ``2``

- ``--no-cleanup``, ``no_cleanup`` : If provided, intermediate files
  generated during process execution are not removed at the end.

Code documentation
------------------
"""


import os
import shutil
import argparse

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


def extract_annotations(blastout_files, indexed_proteome, self_scores,
                        blast_score_ratio, proteome_matches):
    """
    """

    proteome_results = {}
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
                hit_id = hit_dict.get('id', 'not_found')
                hit_description = hit_dict.get('description', 'not_found')
                if hit_description != '':
                    if 'OS=' in hit_description:
                        # get organism name
                        hit_species = (hit_description.split('OS=')[1]).split(' OX=')[0]
                        hit_product = (hit_description.split(hit_id+' ')[1]).split(' OS=')[0]
                    else:
                        hit_species = 'not_found'
                        hit_product = 'not_found'

                    if 'GN=' in hit_description:
                        hit_gene_name = (hit_description.split('GN=')[1]).split(' PE=')[0]
                    else:
                        hit_gene_name = 'not_found'

                proteome_results[locus_id][-1].extend([hit_id, hit_product, hit_gene_name, hit_species])
        # sort from highest to lowest BSR
        proteome_results[locus_id] = sorted(proteome_results[locus_id], key=lambda x: x[3], reverse=True)

    return proteome_results


def proteome_annotations(schema_directory, temp_directory, taxa, blast_score_ratio,
                         cpu_cores, proteome_matches, blast_path):
    """
    """

    # get paths to files with representative sequences
    short_directory = fo.join_paths(schema_directory, ['short'])
    reps_paths = [fo.join_paths(short_directory, [file])
                  for file in os.listdir(short_directory)
                  if file.endswith('.fasta') is True]

    print('Translating representative sequences...', end='')
    # translate representatives for all loci
    translated_reps = fo.join_paths(temp_directory, ['translated_reps'])
    fo.create_directory(translated_reps)

    reps_protein_files = fao.translate_fastas(reps_paths, translated_reps, 11)
    print('done.')

    print('Downloading list of reference proteomes...', end='')
    remote_readme = fo.join_paths(ct.UNIPROT_PROTEOMES_FTP, ['README'])
    local_readme = fo.join_paths(temp_directory,
                                 ['reference_proteomes_readme.txt'])

    # get README file with list of reference proteomes
    res = fo.download_file(remote_readme, local_readme)
    print('done.')

    # get lines with proteomes info for species of interest
    readme_lines = fo.read_lines(local_readme, strip=False)

    selected_proteomes = im.contained_terms(readme_lines, taxa)
    selected_proteomes = [line.strip('\n') for line in selected_proteomes]
    selected_proteomes = [line.split('\t') for line in selected_proteomes]
    print('Found {0} reference proteomes for {1}.'.format(len(selected_proteomes), taxa))
    if len(selected_proteomes) > 0:
        # create directory to store proteomes
        proteomes_directory = fo.join_paths(temp_directory, ['proteomes'])
        fo.create_directory(proteomes_directory)

        proteomes_files = ur.get_proteomes(selected_proteomes, proteomes_directory)

    # uncompress files and concatenate into single FASTA
    uncompressed_proteomes = [fo.unzip_file(file) for file in proteomes_files]
    proteomes_concat = fo.join_paths(proteomes_directory,
                                     ['full_proteome.fasta'])
    proteomes_concat = fo.concatenate_files(uncompressed_proteomes,
                                            proteomes_concat)

    # get self-scores
    # concatenate protein files
    reps_concat = fo.concatenate_files(reps_protein_files,
                                       fo.join_paths(temp_directory, ['reps_concat.fasta']))

    print('\nDetermining self-score of representatives...', end='')
    blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
    makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)
    self_scores = fao.get_self_scores(reps_concat, temp_directory, cpu_cores,
                                      blastp_path, makeblastdb_path)
    print('done.')

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
    print('\n')

    blastout_files = [fo.join_paths(translated_reps, [file])
                      for file in os.listdir(translated_reps)
                      if 'blastout' in file]

    # index proteome file
    indexed_proteome = SeqIO.index(proteomes_concat, 'fasta')

    # process results for each BLASTp
    proteome_results = extract_annotations(blastout_files,
                                           indexed_proteome,
                                           self_scores,
                                           blast_score_ratio,
                                           proteome_matches)

    return proteome_results


def sparql_annotations(loci_files, cpu_cores):
    """ Retrieves annotations from UniProt's SPARQL endpoint.

        Parameters
        ----------
        loci_files : list
            List with the paths to the loci FASTA files.
        cpu_cores : int
            Number of files to process in parallel.

        Returns
        -------
        annotations : list
            List with sublists. Each sublist contains
            the path to the FASTA file of a locus, the
            product name found for that locus and the
            URL to the page of the record that matched the
            locus.
    """

    # create inputs to multiprocessing
    uniprot_args = [[gene, proc_gene] for gene in loci_files]

    # this works with all alleles in the loci to maximize
    # chance of finding info
    workers = cpu_cores if cpu_cores <= 4 else 4
    annotations = mo.map_async_parallelizer(uniprot_args,
                                            mo.function_helper,
                                            workers,
                                            show_progress=True)

    return annotations


def join_annotations(sparql_results, proteome_results, loci_info):
    """
    """

    selected = {}
    for result in sparql_results:
        gene_basename = fo.file_basename(result[0], suffix=False)
        locus_info = [gene_basename]
        locus_info += loci_info.get(gene_basename, ['']*6)
        locus_info += result[1:3]
        locus_info += [proteome_results.get(gene_basename, [])]
        selected[gene_basename] = locus_info

    return selected


def create_annotations_table(annotations, output_directory, header,
                             schema_name, loci_info):
    """
    """

    new_lines = [header]
    for locus, data in annotations.items():
        new_line = [locus]
        if len(loci_info) > 0:
            new_line += data[1:9]
        else:
            new_line += data[7:9]

        if len(data[-1]) > 0:
            relevant_data = [d[4:]+[str(round(d[3],2))] for d in data[-1]]
            proteome_data = list(zip(*relevant_data))
            proteome_data = [';'.join(list(map(str, d))) for d in proteome_data]
            proteome_data = ['' if set(d) == {';'} else d for d in proteome_data]
            new_line.extend(proteome_data)
        new_lines.append(new_line)

    new_lines = ['\t'.join(l) for l in new_lines]
    table_text = '\n'.join(new_lines)

    table_basename = '{0}_annotations.tsv'.format(schema_name)
    output_table = fo.join_paths(output_directory, [table_basename])
    with open(output_table, 'w') as outfile:
        outfile.write(table_text+'\n')

    return output_table


def main(input_files, output_directory, protein_table, blast_score_ratio,
         cpu_cores, taxa, proteome_matches, no_cleanup, blast_path):

    # create output directory
    fo.create_directory(output_directory)

    # create temp directory
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # validate input files
    genes_list = fo.join_paths(temp_directory, ['listGenes.txt'])
    genes_list = pv.check_input_type(input_files, genes_list)
    loci_paths = fo.read_lines(genes_list)

    schema_directory = os.path.dirname(loci_paths[0])
    schema_basename = fo.file_basename(schema_directory)
    print('Schema: {0}'.format(schema_directory))
    print('Number of loci: {0}'.format(len(loci_paths)))

    # find annotations based on reference proteomes for species
    proteome_results = {}
    if taxa is not None:
        proteome_results = proteome_annotations(schema_directory, temp_directory,
                                                taxa, blast_score_ratio,
                                                cpu_cores, proteome_matches,
                                                blast_path)

    # find annotations in SPARQL endpoint
    print('Querying UniProt\'s SPARQL endpoint...')
    sparql_results = sparql_annotations(loci_paths, cpu_cores)

    loci_info = {}
    if protein_table is not None:
        # read cds_info table
        # read "cds_info.tsv" file created by CreateSchema
        table_lines = fo.read_tabular(protein_table)
        for l in table_lines[1:]:
            # create locus identifier based on genome identifier and
            # cds identifier in file
            locus_id = l[0].replace('_', '-')
            locus_id = locus_id + '-protein{0}'.format(l[-2])
            loci_info[locus_id] = l

    annotations = join_annotations(sparql_results, proteome_results, loci_info)

    # table header
    header = ['Locus_ID']
    if len(loci_info) > 0:
        header += table_lines[0]

    header += ['Uniprot_Name', 'UniProt_URL']

    if len(proteome_results) > 0:
        header.extend(['Proteome_ID', 'Proteome_Product',
                       'Proteome_Gene_Name', 'Proteome_Species',
                       'Proteome_BSR'])

    output_table = create_annotations_table(annotations, output_directory,
                                            header, schema_basename,
                                            loci_info)

    if no_cleanup is False:
        shutil.rmtree(temp_directory)

    print('\n\nThe table with new information can be found at:'
          '\n{0}'.format(output_table))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', type=str,
                        required=True, dest='input_files',
                        help='Path to the schema\'s directory or to a file '
                             'with a list of paths to loci FASTA files, one '
                             'per line.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Output directory where the process will '
                             'store intermediate files and save the final '
                             'TSV file with the annotations.')

    parser.add_argument('-t', '--protein-table', type=str,
                        required=False, dest='protein_table',
                        help='Path to the "cds_info.tsv" file created by '
                             'the CreateSchema process.')

    parser.add_argument('--bsr', type=float, required=False,
                        dest='blast_score_ratio',
                        default=0.6,
                        help='BLAST Score Ratio value. This value is only '
                             'used when a taxon/taxa is provided and local '
                             'sequences are aligned against reference '
                             'proteomes.')

    parser.add_argument('--cpu', '--cpu-cores', type=int,
                        required=False, dest='cpu_cores',
                        default=1,
                        help='Number of CPU cores used to run the process.')

    parser.add_argument('--taxa', nargs='+', type=str,
                        required=False, dest='taxa',
                        help='List of scientific names for a set of taxa. The '
                             'process will search for and download reference '
                             'proteomes with terms that match any of the '
                             'provided taxa.')

    parser.add_argument('--pm', type=int, required=False,
                        default=1, dest='proteome_matches',
                        help='Maximum number of proteome matches to report.')

    parser.add_argument('--no-cleanup', action='store_true',
                        required=False, dest='no_cleanup',
                        help='If provided, intermediate files generated '
                             'during process execution are not removed '
                             'at the end.')

    parser.add_argument('--b', '--blast-path', type=pv.check_blast,
                        required=False, default='', dest='blast_path',
                        help='Path to the BLAST executables.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
