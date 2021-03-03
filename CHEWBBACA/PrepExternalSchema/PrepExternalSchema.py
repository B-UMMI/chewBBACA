#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module enables the adaptation of external schemas so that the loci and
alleles present in those schemas can be used with chewBBACA. During the
process, alleles that do not correspond to a complete CDS or that cannot be
translated are discarded from the final schema. One or more alleles of each
gene/locus will be chosen as representatives and included in the 'short'
directory.

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_files`` : Path to the folder containing the fasta files,
  one fasta file per gene/locus (alternatively, a file with a list of paths
  can be given).

    - e.g.: ``/home/user/chewie_schemas/schema_dir``

- ``-o``, ``output_directory`` : The directory where the output files will
  be saved (will create the directory if it does not exist).

    - e.g.: ``/home/user/adapted_schema``

- ``ptf``, ``ptf_path`` : Path to the Prodigal training file that will
  be associated with the adapted schema.

    - e.g.: ``/home/user/training_files/training_file``

- ``--cpu``, ``cpu_cores`` : The number of CPU cores to use (default=1).

    - e.g.: ``4``

- ``--bsr``, ``blast_score_ratio`` : The BLAST Score Ratio value that
  will be used to adapt the external schema (default=0.6).

    - e.g.: ``0.6``

- ``--l``, ``minimum_length`` : Minimum sequence length accepted.
  Sequences with a length value smaller than the value passed to this
  argument will be discarded (default=0).

    - e.g.: ``201``

- ``--t``, ``translation_table`` : Genetic code to use for CDS
  translation (default=11, for Bacteria and Archaea).

    - e.g.: ``11``

- ``--st``, ``size_threshold`` : CDS size variation threshold. At the
  default value of 0.2, alleles with size variation +-20 percent when
  compared to the representative will not be included in the final schema.

    - e.g.: ``0.2``

Code documentation
------------------
"""

import os
import shutil
import argparse
import itertools
import multiprocessing

try:
    from utils import (constants as ct,
                       blast_wrapper as bw,
                       file_operations as fo,
                       fasta_operations as fao,
                       sequence_manipulation as sm,
                       parameters_validation as pv,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 blast_wrapper as bw,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_manipulation as sm,
                                 parameters_validation as pv,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


def bsr_categorizer(blast_results, representatives,
                    representatives_scores, min_bsr, max_bsr):
    """ Determines the BLAST hits that have a BSR below a minimum threshold
        and the BLAST hits that have a BSR above a maximum threshold.

        Parameters
        ----------
        blast_results : list of list
            A list with sublists, each sublist contains information
            about a BLAST hit.
        representatives : list
            List with sequence identifiers of representative
            sequences.
        representatives_scores : dict
            Dictionary with self BLAST raw score for every
            representative.
        min_bsr : float
            Minimum BSR value accepted to consider a sequence as
            a possible new representative.
        max_bsr : float
            Maximum BSR value accepted to consider a sequence as
            a possible new representative.

        Returns
        -------
        List with the following elements:
            high_bsr : list
                list with all sequence identifiers of subject
                sequences that had hits with a BSR higher than the
                maximum defined threshold.
            low_bsr : list
                list with all sequence identifiers of subject
                sequences that had hits with a BSR lower than the
                minimum defined threshold.
    """

    high_bsr = []
    hotspot_bsr = []
    low_bsr = []

    high_reps = {}
    hot_reps = {}
    low_reps = {}

    filtered_results = [res for res in blast_results
                        if res[0] != res[1] and res[1] not in representatives]
    bsr_values = [float(res[2])/float(representatives_scores[res[0]])
                  for res in filtered_results]

    high_bsr = [res[1] for ind, res in enumerate(filtered_results)
                if bsr_values[ind] >= max_bsr]
    low_bsr = [res[1] for ind, res in enumerate(filtered_results)
               if bsr_values[ind] < min_bsr]
    hotspot_bsr = [res[1] for ind, res in enumerate(filtered_results)
                   if bsr_values[ind] >= min_bsr and bsr_values[ind] < max_bsr]

    for ind, res in enumerate(filtered_results):
        if bsr_values[ind] >= min_bsr:
            high_reps.setdefault(res[0], []).append(res[1])
        if bsr_values[ind] < min_bsr:
            low_reps.setdefault(res[0], []).append(res[1])
        if bsr_values[ind] >= min_bsr and bsr_values[ind] < max_bsr:
            hot_reps.setdefault(res[0], []).append(res[1])

    # determine representatives that only led to low BSR
    low_reps = list(set(low_reps) - set(high_reps))

    return [high_bsr, low_bsr, hotspot_bsr, high_reps, low_reps, hot_reps]


def select_candidate(candidates, proteins, seqids,
                     representatives, final_representatives):
    """ Chooses a new representative sequence.

        Parameters
        ----------
        candidates : list
            List with the sequence identifiers of all candidates.
        proteins : dict
            A dictionary with protein identifiers as keys and
            protein sequences as values.
        seqids : list
            A list with the sequence identifiers that still have
            no representative (representatives identifiers are
            included because they have to be BLASTed in order to
            determine their self score).
        representatives : list
            The sequence identifiers of all representatives.

        Returns
        -------
        representatives : list
            The set of all representatives, including the new
            representative that was chosen by the function.
    """

    # with more than one sequence as candidate, select longest
    if len(candidates) > 1:

        # determine length of all candidates
        candidates_len = [(seqid, len(proteins[seqid]))
                          for seqid in candidates]

        # order representative candidates by length descending order
        candidates_len = sorted(candidates_len, key=lambda x: x[1],
                                reverse=True)

        # longest allele is the new representative
        representatives.append(candidates_len[0][0])
        final_representatives.append(candidates_len[0][0])

    # if there is only one candidate, keep that
    elif len(candidates) == 1:

        representatives.append(candidates[0])
        final_representatives.append(candidates[0])

    # if no hit qualifies and there are still sequences
    # without representative
    elif len(candidates) == 0 and \
            len(seqids) > len(representatives):

        # determine length of remaining sequences
        # (representatives not included)
        candidates_len = [(seqid, len(proteins[seqid]))
                          for seqid in seqids
                          if seqid not in representatives]

        # sort by descending length
        candidates_len = sorted(candidates_len, key=lambda x: x[1],
                                reverse=True)

        # longest of remaining sequences is new representative
        representatives.append(candidates_len[0][0])
        final_representatives.append(candidates_len[0][0])

    return [representatives, final_representatives]


def adapt_loci(genes, schema_path, schema_short_path, bsr, min_len,
               table_id, size_threshold, blastp_path, makeblastdb_path):
    """ Adapts a set of genes/loci from an external schema so that
        that schema  can be used with chewBBACA. Removes invalid alleles
        and selects representative alleles to include in the "short" directory.

        Parameters
        ----------
        genes_list : list
            A list with the following elements:

            - List with paths to the files to be processed.
            - Path to the schema directory.
            - Path to the "short" directory.
            - BLAST Score Ratio value.
            - Minimum sequence length value.
            - Genetic code.
            - Sequence size variation threshold.

        Returns
        -------
        invalid_alleles : list
            List with the identifiers of the alleles that were
            determined to be invalid.
        invalid_genes : list
            List with the identifiers of the genes that had no
            valid alleles.
        summary_stats : list of list
            List with one sublist per processed locus. Each
            sublist has four elements:

            - The identifier of the locus.
            - The number of alleles in the external file.
            - The number of alleles that were a valid CDS.
            - The number of representatives determined determined
              by the process.

        The function writes the schema files .
    """

    # divide input list into variables
    summary_stats = []
    invalid_genes = []
    invalid_alleles = []
    for gene in genes:

        representatives = []
        final_representatives = []

        # get gene basename and identifier
        gene_basename = os.path.basename(gene)
        gene_id = gene_basename.split('.f')[0]

        # create paths to gene files in new schema
        gene_file = fo.join_paths(schema_path,
                                  ['{0}{1}'.format(gene_id, '.fasta')])

        gene_short_file = fo.join_paths(schema_short_path,
                                        ['{0}{1}'.format(gene_id, '_short.fasta')])

        # create path to temp working directory for current gene
        gene_temp_dir = fo.join_paths(schema_path,
                                      ['{0}{1}'.format(gene_id, '_temp')])

        # create temp directory for the current gene
        fo.create_directory(gene_temp_dir)

        # dictionaries mapping gene identifiers to DNA sequences
        # and Protein sequences
        gene_seqs, prot_seqs, gene_invalid, seqids_map, total_sequences = \
            sm.get_seqs_dicts(gene, gene_id, table_id, min_len, size_threshold)
        invalid_alleles.extend(gene_invalid)

        # if locus has no valid CDS sequences,
        # continue to next locus
        if len(prot_seqs) == 0:
            shutil.rmtree(gene_temp_dir)
            invalid_genes.append(gene_id)
            summary_stats.append([gene_id, str(total_sequences), '0', '0'])
            continue

        if len(gene_seqs) > 1:
            # identify DNA sequences that code for same protein
            equal_prots = sm.determine_duplicated_seqs(prot_seqs)

            # get only one identifier per protein
            ids_to_blast = [protids[0] for protein, protids in equal_prots.items()]

            # get longest sequence as first representative
            longest = sm.determine_longest(ids_to_blast, prot_seqs)
            representatives.append(longest)
            final_representatives.append(longest)

            # create FASTA file with distinct protein sequences
            protein_file = fo.join_paths(gene_temp_dir,
                                         ['{0}_protein.fasta'.format(gene_id)])
            protein_lines = fao.fasta_lines(ids_to_blast, prot_seqs)
            fo.write_list(protein_lines, protein_file)

            # create blastdb with all distinct proteins
            blastp_db = os.path.join(gene_temp_dir, gene_id)
            bw.make_blast_db(makeblastdb_path, protein_file, blastp_db, 'prot')

            # determine appropriate blastp task (proteins < 30aa need blastp-short)
            blastp_task = bw.determine_blast_task(equal_prots)

            # cycles to BLAST representatives against non-representatives until
            # all non-representatives have a representative
            while len(set(ids_to_blast) - set(representatives)) != 0:

                # create FASTA file with representative sequences
                rep_file = fo.join_paths(gene_temp_dir,
                                         ['{0}_rep_protein.fasta'.format(gene_id)])
                rep_protein_lines = fao.fasta_lines(representatives, prot_seqs)
                fo.write_list(rep_protein_lines, rep_file)

                # create file with seqids to BLAST against
                ids_str = im.concatenate_list([str(i) for i in ids_to_blast], '\n')
                ids_file = fo.join_paths(gene_temp_dir,
                                         ['{0}_ids.txt'.format(gene_id)])
                fo.write_to_file(ids_str, ids_file, 'w', '')

                # BLAST representatives against non-represented
                blast_output = fo.join_paths(gene_temp_dir,
                                             ['{0}_blast_out.tsv'.format(gene_id)])
                # set max_target_seqs to huge number because BLAST only
                # returns 500 hits by default

                blast_stderr = bw.run_blast(blastp_path, blastp_db, rep_file,
                                            blast_output, 1, 1, ids_file,
                                            blastp_task, 100000)
                if len(blast_stderr) > 0:
                    raise ValueError(blast_stderr)

                # import BLAST results
                blast_results = fo.read_tabular(blast_output)

                # get self-score for representatives
                rep_self_scores = {res[1]: res[2] for res in blast_results
                                   if res[0] == res[1]}

                # divide results into high, low and hot BSR values
                hitting_high, hitting_low, hotspots, high_reps, low_reps, hot_reps = \
                    bsr_categorizer(blast_results, representatives,
                                    rep_self_scores, bsr, bsr+0.1)

                excluded_reps = []

                # remove high BSR hits that have representative
                hitting_high = set(hitting_high)
                ids_to_blast = [i for i in ids_to_blast if i not in hitting_high]

                # remove representatives that led to high BSR with subjects that were removed
                prunned_high_reps = {k: [r for r in v if r in ids_to_blast] for k, v in high_reps.items()}
                reps_to_remove = [k for k, v in prunned_high_reps.items() if len(v) == 0]

                excluded_reps.extend(reps_to_remove)

                # determine smallest set of representatives that allow to get all cycle candidates
                excluded = []
                hotspot_reps = set(im.flatten_list(list(hot_reps.values())))
                for rep, hits in hot_reps.items():
                    common = hotspot_reps.intersection(set(hits))
                    if len(common) > 0:
                        hotspot_reps = hotspot_reps - common
                    else:
                        excluded.append(rep)

                excluded_reps.extend(excluded)

                # remove representatives that only led to low BSR
                excluded_reps.extend(low_reps)

                representatives = [rep for rep in representatives if rep not in excluded_reps]
                ids_to_blast = [i for i in ids_to_blast if i not in excluded_reps]

                # determine next representative from candidates
                rep_candidates = list(set(hotspots) - hitting_high)
                # sort to guarantee reproducible results with same datasets
                rep_candidates = sorted(rep_candidates, key=lambda x: int(x))
                representatives, final_representatives = select_candidate(rep_candidates,
                                                                          prot_seqs,
                                                                          ids_to_blast,
                                                                          representatives,
                                                                          final_representatives)

                # remove files created for current gene iteration
                os.remove(rep_file)
                os.remove(blast_output)
                os.remove(ids_file)

        else:
            final_representatives = list(prot_seqs.keys())

        # write schema file with all alleles
        gene_lines = fao.fasta_lines(list(gene_seqs.keys()), gene_seqs)
        fo.write_list(gene_lines, gene_file)

        # get total number of valid sequences
        valid_sequences = len(gene_lines)

        # write schema file with representatives
        final_representatives = [seqids_map[rep] for rep in final_representatives]
        gene_rep_lines = fao.fasta_lines(final_representatives, gene_seqs)
        fo.write_list(gene_rep_lines, gene_short_file)

        # get number of representatives
        representatives_number = len(gene_rep_lines)

        summary_stats.append([gene_id,
                              str(total_sequences),
                              str(valid_sequences),
                              str(representatives_number)])

        shutil.rmtree(gene_temp_dir)

    return [invalid_alleles, invalid_genes, summary_stats]


def main(input_files, output_directory, cpu_cores, blast_score_ratio,
         minimum_length, translation_table, ptf_path, size_threshold,
         blast_path):

    print('Adapting schema in the following '
          'directory:\n{0}'.format(os.path.abspath(input_files)))
    print('Prodigal training file:\n{0}'.format(ptf_path))
    print('Number of cores: {0}'.format(cpu_cores))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum accepted sequence length: {0}'.format(minimum_length))
    print('Size threshold: {0}'.format(size_threshold))

    # define output paths
    schema_path = os.path.abspath(output_directory)
    schema_short_path = fo.join_paths(schema_path, ['short'])

    # create output directories
    # check if they exist first
    fo.create_directory(schema_path)
    fo.create_directory(schema_short_path)

    # list schema gene files
    genes_file = pv.check_input_type(input_files,
                                     os.path.join(output_directory, 'schema_genes.txt'))

    # import list of schema files
    with open(genes_file, 'r') as gf:
        genes_list = [line.rstrip('\n') for line in gf]
    os.remove(genes_file)

    print('Number of genes to adapt: {0}\n'.format(len(genes_list)))

    print('Determining the total number of alleles and '
          'allele mean length per gene...\n'.format())

    # count number of sequences and mean length per gene
    genes_info = []
    genes_pools = multiprocessing.Pool(processes=cpu_cores)
    gp = genes_pools.map_async(fao.gene_seqs_info, genes_list,
                               callback=genes_info.extend)
    gp.wait()

    # split files according to number of sequences and sequence mean length
    # in each file to pass even groups of sequences to all cores
    even_genes_groups = mo.split_genes_by_core(genes_info, cpu_cores*4,
                                               'seqcount')
    # with few inputs, some sublists might be empty
    even_genes_groups = [i for i in even_genes_groups if len(i) > 0]

    # add common arguments
    blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
    makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)
    even_genes_groups = [[i, schema_path, schema_short_path,
                          blast_score_ratio, minimum_length,
                          translation_table, size_threshold,
                          blastp_path, makeblastdb_path,
                          adapt_loci] for i in even_genes_groups]

    print('Adapting {0} genes...\n'.format(len(genes_list)))

    invalid_data = mo.map_async_parallelizer(even_genes_groups,
                                             mo.function_helper,
                                             cpu_cores,
                                             show_progress=True)

    # define paths and write files with list of invalid
    # alleles and invalid genes
    output_schema_basename = os.path.basename(output_directory.rstrip('/'))
    schema_parent_directory = os.path.dirname(schema_path)

    # write file with alleles that were determined to be invalid
    invalid_alleles = [sub[0] for sub in invalid_data]
    invalid_alleles = list(itertools.chain.from_iterable(invalid_alleles))
    invalid_alleles_file = os.path.join(schema_parent_directory,
                                        '{0}_{1}'.format(output_schema_basename, 'invalid_alleles.txt'))

    with open(invalid_alleles_file, 'w') as inv:
        lines = ['{0}: {1}\n'.format(allele[0], allele[1]) for allele in invalid_alleles]
        inv.writelines(lines)

    # write file with identifiers of genes that had no valid alleles
    invalid_genes = [sub[1] for sub in invalid_data]
    invalid_genes = list(itertools.chain.from_iterable(invalid_genes))
    invalid_genes_file = os.path.join(schema_parent_directory,
                                      '{0}_{1}'.format(output_schema_basename, 'invalid_genes.txt'))

    with open(invalid_genes_file, 'w') as inv:
        invalid_geqids = '\n'.join(invalid_genes)
        inv.write(invalid_geqids)

    stats_lines = [sub[2] for sub in invalid_data]
    stats_lines = list(itertools.chain.from_iterable(stats_lines))
    stats_lines = ['\t'.join(line) for line in stats_lines]
    stats_genes_file = '{0}/{1}_{2}'.format(schema_parent_directory,
                                            output_schema_basename,
                                            'summary_stats.txt')

    with open(stats_genes_file, 'w') as stats:
        summary_stats_text = '\n'.join(stats_lines)
        stats.write('Gene\tTotal_alleles\tValid_alleles\tNumber_representatives\n')
        stats.write(summary_stats_text)

    print('\n\nNumber of invalid genes: {0}'.format(len(invalid_genes)))
    print('Number of invalid alleles: {0}'.format(len(invalid_alleles)))

    print('\nSuccessfully adapted {0}/{1} genes present in the '
          'input schema.'.format(len(genes_list)-len(invalid_genes),
                                 len(genes_list)))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True, dest='input_files',
                        help='Path to the folder containing the fasta files, '
                             'one fasta file per gene/locus (alternatively, '
                             'a file with a list of paths can be given).')

    parser.add_argument('-o', type=str, required=True, dest='output_directory',
                        help='The directory where the output files will be '
                             'saved (will create the directory if it does not '
                             'exist).')

    parser.add_argument('--ptf', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file that '
                             'will be associated with the adapted schema.')

    parser.add_argument('--cpu', type=int, required=False, default=1,
                        dest='cpu_cores',
                        help='The number of CPU cores to use (default=1).')

    parser.add_argument('--bsr', type=float, required=False, default=0.6,
                        dest='blast_score_ratio',
                        help='The BLAST Score Ratio value that will be '
                             'used to adapt the external schema (default=0.6).')

    parser.add_argument('--l', type=int, required=False, default=0,
                        dest='minimum_length',
                        help='Minimum sequence length accepted. Sequences with'
                        ' a length value smaller than the value passed to this'
                        ' argument will be discarded (default=0).')

    parser.add_argument('--t', type=int, required=False, default=11,
                        dest='translation_table',
                        help='Genetic code to use for CDS translation.'
                        ' (default=11, for Bacteria and Archaea)')

    parser.add_argument('--st', type=float, required=True,
                        default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent when compared to the representative '
                             'will not be included in the final schema.')

    parser.add_argument('--b', type=pv.check_blast, required=False,
                        default='', dest='blast_path',
                        help='Path to the BLAST executables.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    main(**vars(args))
