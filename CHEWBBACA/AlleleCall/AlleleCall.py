#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This module enables

Expected input
--------------

The process expects the following variables whether through command line
execution or invocation of the :py:func:`main` function:

- ``-i``, ``input_files`` : Path to the directory that contains the input
  FASTA files. Alternatively, a single file with a list of paths to FASTA
  files, one per line.

    - e.g.: ``/home/user/genomes``

- ``-o``, ``output_directory`` : Output directory where the process will
  store intermediate files and create the schema's directory.

    - e.g.: ``/home/user/schemas/new_schema``

- ``--ptf``, ``ptf_path`` : Path to the Prodigal training file.

    - e.g.: ``/home/user/training_files/species.trn``

- ``--bsr``, ``blast_score_ratio`` : BLAST Score Ratio value.

    - e.g.: ``0.6``

- ``--l``, ``minimum_length`` : Minimum sequence length. Coding sequences
  shorter than this value are excluded.

    - e.g.: ``201``

- ``--t``, ``translation_table`` : Genetic code used to predict genes and
  to translate coding sequences.

    - e.g.: ``11``

- ``--st``, ``size_threshold`` : CDS size variation threshold. Added to the
  schema's config file and used to identify alleles with a length value that
  deviates from the locus length mode during the allele calling process.

    - e.g.: ``0.2``

- ``--w``, ``word_size`` : word size used to generate k-mers during the
  clustering step.

    - e.g.: ``5``

- ``--ws``, ``window_size`` : window size value. Number of consecutive
  k-mers included in each window to determine a minimizer.

    - e.g.: ``5``

- ``--cs``, ``clustering_sim`` : clustering similarity threshold. Minimum
  decimal proportion of shared distinct minimizers for a sequence to be
  added to a cluster.

    - e.g.: ``0.2``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores used to run the process.

    - e.g.: ``4``

- ``--b``, ``blast_path`` : Path to the BLAST executables.

    - e.g.: ``/home/software/blast``

- ``--pm``, ``prodigal_mode`` : Prodigal running mode.

    - e.g.: ``single``

- ``--CDS``, ``cds_input`` : If provided, input is a single or several FASTA
  files with coding sequences (skips gene prediction and CDS extraction).

    - e.g.: ``/home/user/coding_sequences_files``

- ``--no-cleanup``, ``no_cleanup`` : If provided, intermediate files
  generated during process execution are not removed at the end.

Code documentation
------------------
"""


import os
import re
import csv
import sys
import time
import argparse
from collections import Counter

from Bio import SeqIO

try:
    from utils import (constants as ct,
                       blast_wrapper as bw,
                       core_functions as cf,
                       file_operations as fo,
                       fasta_operations as fao,
                       sequence_manipulation as sm,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 blast_wrapper as bw,
                                 core_functions as cf,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


# import module to determine variable size
import get_varSize_deep as gs


def match_regex(text, pattern):
    """
    """

    match_span = re.search(pattern, text).span()
    match_str = text[match_span[0]:match_span[1]]

    return match_str


###############################
# Translating everything, maybe we should only create Fasta files with distinct records
def translate_fastas(fasta_files, output_directory,
                     translation_table, cpu_cores):
    """
    """

    translation_inputs = im.divide_list_into_n_chunks(fasta_files,
                                                      len(fasta_files))

    common_args = [output_directory, translation_table]

    # add common arguments to all sublists
    translation_inputs = im.multiprocessing_inputs(translation_inputs,
                                                   common_args,
                                                   fao.translate_fasta)

    # run Prodigal to predict genes
    protein_files = mo.map_async_parallelizer(translation_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    return protein_files


def select_high_score(input_file, blast_score_ratio,
                      ids_map, output_directory):
    """
    """

    # import BLASTp results
    current_results = process_blast_results(input_file,
                                            blast_score_ratio,
                                            ids_map)

    # only save and evaluate results for loci that
    # had at least one high BSR match
    if len(current_results[0]) > 0:
        locus = fo.get_locus_id(input_file)
        current_results.append(locus)
        # save data to pickle
        pickle_out = fo.join_paths(output_directory, [locus+'_bsr'])
        fo.pickle_dumper(current_results, pickle_out)

        return pickle_out


def group_by_genome(file, sequences, cds_hashtable, prot_hashtable, cds_index):
    """
    """

    high_bsr, representatives_info, locus = fo.pickle_loader(file)
    genomes_hits = [{}, representatives_info]
    for k, v in high_bsr.items():
        target_prot_hash = im.hash_sequence(sequences[k])
        target_seqids = prot_hashtable[target_prot_hash][1:]
        for g in target_seqids:
            target_dna = str(cds_index.get(g).seq)
            target_dna_len = len(target_dna)
            target_dna_hash = im.hash_sequence(target_dna)

            # get ids for all genomes with same CDS as representative
            all_genomes_with_cds = cds_hashtable[target_dna_hash][1:]
            for a in all_genomes_with_cds:
                genomes_hits[0].setdefault(a, []).append((k, target_prot_hash,
                                                          target_dna_hash, target_dna_len,
                                                          high_bsr[k]))

    return [locus, genomes_hits]


def update_classification(genome_id, locus_results, classification,
                          representative_id, blast_score_ratio,
                          multiple_classification='NIPH'):
    """
    """

    if genome_id in locus_results:
        locus_results[genome_id][0] = multiple_classification
        locus_results[genome_id].append((representative_id,
                                         classification,
                                         blast_score_ratio))
    else:
        locus_results[genome_id] = [classification,
                                    (representative_id,
                                     classification,
                                     blast_score_ratio)]

    return locus_results


def schema_exact_matches(loci_files, output_directory,
                         distinct_htable, loci_basenames):
    """
    """

    dna_exact_hits = 0
    dna_matches_ids = []
    classification_files = []
    for file in loci_files:
        pickle_out = fo.join_paths(output_directory,
                                   [loci_basenames[file]+'_results'])
        em_results = locus_exact_matches(file, distinct_htable)
        # save classifications
        fo.pickle_dumper(em_results[0], pickle_out)
        classification_files.append(pickle_out)
        # extend list of matched seqids
        dna_matches_ids.extend(em_results[1])
        dna_exact_hits += em_results[2]

    return [classification_files, dna_matches_ids, dna_exact_hits]


def locus_exact_matches(fasta_file, distinct_table):
    """
    """

    seq_generator = SeqIO.parse(fasta_file, 'fasta')

    exact_hits = {}
    matches_ids = []
    total_matches = 0
    for rec in seq_generator:
        seqid = (rec.id).split('_')[0]
        sequence = str(rec.seq.upper())
        sequence_hash = im.hash_sequence(sequence)
        # file cannot have duplicated alleles or it will
        # classify based on same DNA sequence twice
        if sequence_hash in distinct_table:
            # exclude element in index 0, it is the seqid chosen as
            # representative and is repeated in index 1
            current_matches = distinct_table[sequence_hash]
            total_matches += len(current_matches) - 1
            for m in current_matches[1:]:
                exact_hits = update_classification(m, exact_hits, 'EXC',
                                                   seqid, 1.0, 'NIPHEM')

            matches_ids.append(current_matches[0])

    return [exact_hits, matches_ids, total_matches]



# fasta_file = protein_files[0]
# distinct_prot_table = distinct_pseqids
# distinct_dna_table = dna_distinct_htable
# locus = fo.get_locus_id(fasta_file)
# pickle_out = fo.join_paths(preprocess_dir, [locus+'_results'])
# output_file = pickle_out
# protein_index = protein_index
# dna_index = dna_index
def protein_exact_matches(fasta_file, distinct_prot_table,
                          distinct_dna_table, output_file,
                          protein_index, dna_index):
    """
    """

    seq_generator = SeqIO.parse(fasta_file, 'fasta')

    seen_dna = []
    seen_prot = []
    total_cds = 0
    total_prots = 0
    exact_hits = {}
    exact_prot_hashes = []
    for rec in seq_generator:
        protid = (rec.id).split('_')[0]
        protein = str(rec.seq.upper())
        prot_hash = im.hash_sequence(protein)
        # do not proceed if distinct protein sequence has already been seen
        # different alleles might code for the same protein
        if prot_hash in distinct_prot_table and prot_hash not in seen_prot:
            # get protids for distinct CDSs
            protids = distinct_prot_table[prot_hash][1:]
            total_prots += len(protids)
            exact_prot_hashes.extend(protids)
            for p in protids:
                cds = str(dna_index.get(p).seq)
                cds_hash = im.hash_sequence(cds)
                # get all CDS ids associated associated to hash
                seqids = distinct_dna_table[cds_hash]
                total_cds += len(seqids) - 1
                for s in seqids[1:]:
                    if s not in exact_hits:
                        if cds_hash not in seen_dna:
                            exact_hits = update_classification(s, exact_hits,
                                                               'INF', p, 1.0)
                            seen_dna.append(cds_hash)
                        else:
                            exact_hits = update_classification(s, exact_hits,
                                                               'EXC', p, 1.0)
                    else:
                        if cds_hash not in seen_dna:
                            exact_hits = update_classification(s, exact_hits,
                                                               'INF', p, 1.0)
                            exact_hits[s].append((seqids[0], 'INF'))
                        else:
                            exact_hits = update_classification(s, exact_hits,
                                                               'EXC', p, 1.0)

            seen_prot.append(prot_hash)

    # merge results with DNA exact matches
    locus_results = fo.pickle_loader(output_file)
    for k, v in exact_hits.items():
        locus_results.setdefault(k, [v[0]]).extend(v[1:])
        if len(v[1:]) > 1:
            locus_results[k][0] = 'NIPH'

    fo.pickle_dumper(locus_results, output_file)

    return [exact_prot_hashes, total_prots, total_cds]


def process_blast_results(blast_output, blast_score_ratio, ids_mapping=None):
    """
    """

    current_results = fo.read_tabular(blast_output)

    # get allele identifier for the representative
    if ids_mapping is not None:
        representatives_info = {ids_mapping[r[0]]: (((int(r[3])*3)+3), float(r[5]))
                                for r in current_results
                                if r[0] == r[4]}
    else:
        representatives_info = {r[0]: (((int(r[3])*3)+3), float(r[5]))
                                for r in current_results
                                if r[0] == r[4]}

    # exclude representative self-alignment
    current_results = [r
                       for r in current_results
                       if r[0] != r[4]]

    # only keep the match with the greatest score for each target
    # several representatives from the same locus might match the target
    highest_scores = {}
    for r in current_results:
        if r[4] not in highest_scores:
            highest_scores[r[4]] = r
        elif r[4] in highest_scores:
            current_score = float(r[5])
            previous_score = float(highest_scores[r[4]][5])
            if current_score > previous_score:
                highest_scores[r[4]] = r

    # determine BSR values
    if ids_mapping is not None:
        bsr_values = {ids_mapping[k]: (float(v[5])/representatives_info[ids_mapping[v[0]]][1],
                                       (int(v[1])-1)*3,  # subtract 1 to exclude start position
                                       (int(v[2])+1)*3,  # add 1 to include stop codon
                                       ids_mapping[v[0]],
                                       ids_mapping[v[0]].split('_')[-1])
                      for k, v in highest_scores.items()}
    else:
        bsr_values = {k: (float(v[5])/representatives_info[v[0]][1],
                          (int(v[1])-1)*3,  # subtract 1 to exclude start position
                          (int(v[2])+1)*3,  # add 1 to include stop codon
                          v[0],
                          v[0].split('_')[-1])
                      for k, v in highest_scores.items()}

    # only keep matches above BSR threshold
    high_bsr = {k: v
                for k, v in bsr_values.items()
                if v[0] >= blast_score_ratio}

    return [high_bsr, representatives_info]


def classify_alleles(locus, genomes_matches, representatives_info,
                     inv_map, contigs_lengths, locus_results_file,
                     locus_mode, temp_directory, size_threshold):
    """
    """

    locus_results = fo.pickle_loader(locus_results_file)

    seen_dna = []
    seen_prot = []
    excluded = []
    inf_matches = 0
    genome_coordinates = {}
    representative_candidates = []
    for genome, matches in genomes_matches.items():
        current_g = inv_map[genome]
        clen = contigs_lengths[current_g]

        # open pickle for genome and get coordinates
        genome_cds_file = fo.join_paths(temp_directory, ['2_cds_extraction', current_g+'_cds_hash'])
        genome_cds_coordinates = fo.pickle_loader(genome_cds_file)
        for m in matches:
            target_seqid = m[0]
            excluded.append(target_seqid)
            target_dna_len = m[3]
            bsr = m[4][0]
            rep_length = representatives_info[m[4][3]][0]
            rep_alleleid = m[4][4]
            target_prot_hash = m[1]
            target_dna_hash = m[2]

            # careful about offset!!!
            # determine right match position on representative allele
            rep_right_pos = rep_length - m[4][2]
            # determine left match position on representative allele
            rep_left_pos = m[4][1]

            if target_dna_hash in seen_dna:
                locus_results = update_classification(genome, locus_results, 'EXC',
                                                      rep_alleleid, 1.0)
                inf_matches += 1
                continue

            # check if CDS DNA sequence is not in schema but matches a translated allele
            if target_prot_hash in seen_prot:
                locus_results = update_classification(genome, locus_results, 'INF',
                                                      rep_alleleid, 1.0)
                inf_matches += 1
                # add DNA hash to classify the next match as exact
                seen_dna.append(target_dna_hash)
                continue

            # there is no exact match in the schema, perform full evaluation
            # get contig lengths
            genome_coordinates[target_seqid] = genome_cds_coordinates[target_dna_hash][0]
            contig_left_pos = int(genome_coordinates[target_seqid][1])
            contig_right_pos = int(genome_coordinates[target_seqid][2])
            matched_contig_len = clen[genome_coordinates[target_seqid][0]]
            ###########################################################
            # need to invert when necessary
            # check if this is necessary
            # get total number of bases left outside match in contig
            if contig_right_pos > contig_left_pos:
                contig_right_rest = matched_contig_len - contig_right_pos
                contig_left_rest = contig_left_pos
                rep_right_rest = rep_right_pos
                rep_left_rest = rep_left_pos
            # reverse values because CDS was identified in reverse strand
            elif contig_right_pos < contig_left_pos:
                print('lol')
                contig_left_rest = matched_contig_len - contig_right_pos
                contig_right_rest = contig_left_pos
                # also need to reverse values for representative
                rep_right_rest = rep_left_pos
                rep_left_rest = rep_right_pos
                sys.exit()

            # check LOTSC
            if contig_left_rest < rep_left_rest and contig_right_rest < rep_right_rest:
                locus_results = update_classification(genome, locus_results, 'LOTSC',
                                                      rep_alleleid, bsr)
                continue
            # check if PLOT
            elif contig_left_rest < rep_left_rest:
                locus_results = update_classification(genome, locus_results, 'PLOT3',
                                                      rep_alleleid, bsr)
                continue
            elif contig_right_rest < rep_right_rest:
                locus_results = update_classification(genome, locus_results, 'PLOT5',
                                                      rep_alleleid, bsr)
                continue

            # check if ASM or ALM
            # we only need to evaluate one of the genomes, if they are ASM/ALM we can classify all of them as the same!
            if target_dna_len < (locus_mode[0]-(locus_mode[0])*size_threshold):
                locus_results = update_classification(genome, locus_results, 'ASM',
                                                      rep_alleleid, bsr)
                continue
            elif target_dna_len > (locus_mode[0]+(locus_mode[0])*size_threshold):
                locus_results = update_classification(genome, locus_results, 'ALM',
                                                      rep_alleleid, bsr)
                continue

            # add INF
            locus_results = update_classification(genome, locus_results, 'INF',
                                                  rep_alleleid, bsr)

            # also add sequence as representative candidate
            if bsr >= blast_score_ratio and bsr < blast_score_ratio+0.1:
                representative_candidates.append((genome, target_seqid))
                excluded.remove(target_seqid)

            # add hash of newly inferred allele
            seen_dna.append(target_dna_hash)
            seen_prot.append(target_prot_hash)

        # update locus mode value if a new allele was inferred
        if locus_results[genome][0] == 'INF':
            # append length of inferred allele to list with allele sizes
            locus_mode[1].append(target_dna_len)
            #loci_modes[locus][1].append(target_dna_len)
            # compute mode
            #loci_modes[locus][0] = sm.determine_mode(loci_modes[locus][1])[0]
            locus_mode[0] = sm.determine_mode(locus_mode[1])[0]

    # save updated results
    fo.pickle_dumper(locus_results, locus_results_file)

    return {locus: [locus_results_file, locus_mode, inf_matches,
                    excluded, representative_candidates]}


#input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids.txt'
#input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids32.txt'
input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids320.txt'
#input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids_2058spyogenes.txt'
output_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/test_allelecall'
ptf_path = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed/Streptococcus_agalactiae.trn'
#ptf_path = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/spyogenes_schema_processed/Streptococcus_pyogenes.trn'
blast_score_ratio = 0.6
minimum_length = 201
translation_table = 11
size_threshold = 0.2
word_size = 5
window_size = 5
clustering_sim = 0.2
representative_filter = 0.9
intra_filter = 0.9
cpu_cores = 6
blast_path = '/home/rfm/Software/anaconda3/envs/spyder/bin'
prodigal_mode = 'single'
cds_input = False
only_exact = False
#schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed'
schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema_called_320/schema_seed'
#schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/spyogenes_schema_processed'
def allele_calling(input_files, schema_directory, output_directory, ptf_path,
                   blast_score_ratio, minimum_length, translation_table,
                   size_threshold, word_size, window_size, clustering_sim,
                   representative_filter, intra_filter, cpu_cores, blast_path,
                   prodigal_mode, cds_input, only_exact):
    """
    """

    start = time.time()
    # define directory for temporary files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # read file with paths to input files
    fasta_files = fo.read_lines(input_files, strip=True)

    # sort paths to FASTA files
    fasta_files = im.sort_data(fasta_files, sort_key=lambda x: x.lower())
    # map input basenames to
    inputs_basenames = fo.mapping_function(fasta_files,
                                           fo.file_basename, [False])

    # create dictionary with integer mapping for inputs
    # this reduces memory usage during sequence deduplication
    # because the integer identifier is added, not the string
    map_ids = im.integer_mapping(list(inputs_basenames.values()))
    inv_map = im.invert_dictionary(map_ids)

    if cds_input is False:

        print('Number of inputs: {0}'.format(len(fasta_files)))

        # create directory to store files with Prodigal results
        prodigal_path = fo.join_paths(temp_directory, ['1_gene_prediction'])
        fo.create_directory(prodigal_path)

        # run Prodigal to determine CDSs for all input genomes
        print('\nPredicting gene sequences...\n')

        # gene prediction step
        gp_results = cf.predict_genes(fasta_files, ptf_path,
                                      translation_table, prodigal_mode,
                                      cpu_cores, prodigal_path,
                                      output_directory)

        failed, failed_file = gp_results

        if len(failed) > 0:
            print('\nFailed to predict genes for {0} genomes'
                  '.'.format(len(failed)))
            print('Make sure that Prodigal runs in meta mode (--pm meta) '
                  'if any input file has less than 100kbp.')
            print('Info for failed cases stored in: {0}'.format(failed_file))

        # remove failed genomes from paths
        fasta_files = im.filter_list(fasta_files, failed)

        if len(fasta_files) == 0:
            sys.exit('\nCould not predict gene sequences from any '
                     'of the input files.\nPlease provide input files '
                     'in the accepted FASTA format.')

        # CDS extraction step
        print('\n\nExtracting coding sequences...\n')
        # create output directory
        cds_extraction_path = fo.join_paths(temp_directory,
                                            ['2_cds_extraction'])
        fo.create_directory(cds_extraction_path)
        eg_results = cf.extract_genes(fasta_files, prodigal_path,
                                      cpu_cores, cds_extraction_path,
                                      output_directory)
        cds_files, total_extracted = eg_results

        print('\n\nExtracted a total of {0} coding sequences from {1} '
              'genomes.'.format(total_extracted, len(fasta_files)))
    else:
        cds_files = fasta_files
        print('Number of inputs: {0}'.format(len(cds_files)))

    # create directory to store files from pre-process steps
    preprocess_dir = fo.join_paths(temp_directory, ['3_cds_preprocess'])
    fo.create_directory(preprocess_dir)

    # DNA sequences deduplication step
    # keep hash of unique sequences and a list with the integer
    # identifiers of genomes that have those sequences
    print('\nRemoving duplicated DNA sequences...', end='')
    distinct_dna_template = 'distinct_cds_{0}.fasta'
    dna_dedup_results = cf.exclude_duplicates(cds_files, preprocess_dir,
                                              cpu_cores, distinct_dna_template,
                                              map_ids, False)

    dna_distinct_htable, distinct_file, repeated = dna_dedup_results
    print('removed {0} sequences.'.format(repeated))

    # get size of hash table
    size1 = gs.convert_bytes(dna_distinct_htable, set())
    print(size1.split('\n')[2])

    # remove small sequences!?

    ## Check if reverse complement matches?
    # get list of loci files
    print('Getting list of loci...', end='')
    loci_files = fo.listdir_fullpath(schema_directory,
                                     substring_filter='.fasta')

    loci_basenames = fo.mapping_function(loci_files, fo.get_locus_id, [])
    print('schema has {0} loci.'.format(len(loci_files)))

    print('Finding DNA exact matches...', end='')
    # find exact DNA matches
    # Add multiprocessing?
    # use psutil library to get free memory and decide maximum number of processes to spawn?
    dna_em_results = schema_exact_matches(loci_files, preprocess_dir,
                                          dna_distinct_htable, loci_basenames)

    classification_files, dna_matches_ids, dna_exact_hits = dna_em_results

    print('found {0} exact matches (matching {1} alleles).'
          ''.format(dna_exact_hits, len(dna_matches_ids)))

    if only_exact is True:
        # return classification files for creation of output files
        end = time.time()
        delta = end - start
        print('\n', delta/60)
        sys.exit('Finished allele calling.')

    # create new Fasta file without the DNA sequences that were exact matches
    dna_index = SeqIO.index(distinct_file, 'fasta')
    unique_fasta = fo.join_paths(preprocess_dir, ['dna_non_exact.fasta'])
    total_selected, selected_ids = fao.exclude_sequences_by_id(distinct_file,
                                                               dna_matches_ids,
                                                               dna_index,
                                                               unique_fasta)

    # translate DNA sequences and identify duplicates
    # sequence translation step
    print('\nTranslating {0} DNA sequences...'.format(total_selected))

    ts_results = cf.translate_sequences(selected_ids, unique_fasta,
                                        preprocess_dir, translation_table,
                                        minimum_length, cpu_cores)

    dna_file, protein_file, ut_seqids, ut_lines = ts_results

    print('Removed {0} DNA sequences that could not be '
          'translated.'.format(len(ut_seqids)))

    # write info about invalid alleles to file
    invalid_alleles_file = fo.join_paths(output_directory,
                                         ['invalid_cds.txt'])
    invalid_alleles = im.join_list(ut_lines, '\n')
    fo.write_to_file(invalid_alleles, invalid_alleles_file, 'w', '\n')
    print('Info about untranslatable and small sequences '
          'stored in {0}'.format(invalid_alleles_file))

    # protein sequences deduplication step
    print('\nRemoving duplicated protein sequences...', end='')
    distinct_prot_template = 'distinct_prots_{0}.fasta'
    ds_results = cf.exclude_duplicates([protein_file], preprocess_dir, 1,
                                       distinct_prot_template, map_ids, True)
    print('removed {0} sequences.'.format(ds_results[2]))
    distinct_pseqids = ds_results[0]

    # translate loci files
    print('Translating schema alleles...')
    protein_dir = os.path.join(temp_directory, '4_protein_dir')
    fo.create_directory(protein_dir)
    protein_files = translate_fastas(loci_files, protein_dir,
                                     translation_table, cpu_cores)

    # create protein file index
    protein_index = SeqIO.index(ds_results[1], 'fasta')
    # identify exact matches at protein level
    # exact matches are new alleles that can be added to the schema
    print('\nFinding protein exact matches...', end='')
    exc_cds = 0
    exc_prot = 0
    exact_phashes = []
    for file in protein_files:
        locus = fo.get_locus_id(file)
        pickle_out = fo.join_paths(preprocess_dir, [locus+'_results'])
        em_results = protein_exact_matches(file, distinct_pseqids,
                                           dna_distinct_htable, pickle_out,
                                           protein_index, dna_index)
        exact_phashes.extend(em_results[0])
        exc_prot += em_results[1]
        exc_cds += em_results[2]

    # need to determine number of distinct proteins to print correct info!!!
    print('found {0} protein exact matches (matching {1} distinct proteins).'
          ''.format(exc_prot, len(exact_phashes)))

    print('Protein matches correspond to {0} DNA matches.'.format(exc_cds))

    # create new Fasta file without the Protein sequences that were exact matches
    unique_pfasta = fo.join_paths(preprocess_dir, ['protein_non_exact.fasta'])
    protein_index = SeqIO.index(ds_results[1], 'fasta')
    selected_ids = fao.exclude_sequences_by_id(ds_results[1], exact_phashes,
                                               protein_index, unique_pfasta)

    # translate schema representatives
    print('Translating schema representatives...')
    rep_dir = os.path.join(schema_directory, 'short')
    rep_list = [os.path.join(rep_dir, f) for f in os.listdir(rep_dir) if f.endswith('.fasta')]
    protein_repfiles = translate_fastas(rep_list, protein_dir,
                                        translation_table, cpu_cores)

    # cluster protein sequences
    # protein clustering step
    # read protein sequences
    proteins = fao.import_sequences(unique_pfasta)

    # create directory to store clustering data
    clustering_dir = fo.join_paths(temp_directory, ['5_clustering'])
    fo.create_directory(clustering_dir)

    # BLASTp clusters step
    blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
    makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)

    # concatenate all representative
    concat_reps = os.path.join(protein_dir, 'concat_reps.fasta')
    fo.concatenate_files(protein_repfiles, concat_reps)

    prot_index = SeqIO.index(fo.join_paths(preprocess_dir, ['protein.fasta']), 'fasta')
    dna_cds_index = SeqIO.index(unique_fasta, 'fasta')

    # create index for representative sequences
    rep_proteins = fao.import_sequences(concat_reps)

    representatives = im.kmer_index(rep_proteins, 5)[0]
    cs_results = cf.cluster_sequences(proteins, word_size, window_size,
                                      clustering_sim, representatives, False,
                                      1, 30, clustering_dir, cpu_cores,
                                      'clusters', True, False)

    # remove singletons
    clusters = {k: v for k, v in cs_results.items() if len(v) > 0}

    # create file with all proteins, including loci representatives
    if len(clusters) > 0:
        blasting_dir = fo.join_paths(clustering_dir, ['cluster_blaster'])
        fo.create_directory(blasting_dir)

        # create Fasta file with remaining proteins and representatives
        fo.concatenate_files([unique_pfasta, concat_reps], os.path.join(blasting_dir, 'all_prots.fasta'))

        all_prots = os.path.join(blasting_dir, 'all_prots.fasta')
        all_proteins = fao.import_sequences(all_prots)

        blast_results, ids_dict = cf.blast_clusters(clusters, all_proteins,
                                                    blasting_dir, blastp_path,
                                                    makeblastdb_path, cpu_cores,
                                                    'blast', True)

        blast_files = im.flatten_list(blast_results)

    # for ASM and ALM calssification I only need to determine it for the one of the ids associated to that protein
    # the other proteins have the same length and the classification will be the same
    # for PLOT cases I need to check all cases

    # add id of representative CDS, get hash for that sequence and then open file with genome CDS to get CDS coordinates
    # hash file with DNA sequences to get DNA sequence to hash?

    # add hash of inferred alleles to dict with results per locus
    # we can check if hash is already in the dict to avoid calculating everyhting, effectivelly detecting new exact matches!
    # we can store the DNA hash and the protein hash. DNA hash is an exact match and protein hash is a new allele.

    # group files for representatives of the same locus
    loci_results = {}
    for f in blast_files:
        # create function to accept regex pattern and extract string?
        locus_rep = ids_dict[match_regex(f, r'seq_[0-9]+')]
        locus_id = match_regex(locus_rep, r'^.+-protein[0-9]+')
        loci_results.setdefault(locus_id, []).append(f)

    # concatenate results for representatives of the same locus
    concatenated_files = []
    concatenate_file_template = '{0}_concatenated_blastout.tsv'
    for locus, files in loci_results.items():
        outfile = fo.join_paths(clustering_dir,
                                [concatenate_file_template.format(locus)])
        fo.concatenate_files(files, outfile)
        concatenated_files.append(outfile)

    # import results and determine BSR
    high_scoring = []
    for f in concatenated_files:
        # exclude results in the BSR+0.1 threshold
        # that might be representative candidates
        outfile = select_high_score(f, blast_score_ratio+0.1,
                                    ids_dict, clustering_dir)
        if outfile is not None:
            high_scoring.append(outfile)

    # group results for each locus per genome
    loci_results = {}
    for f in high_scoring:
        results = group_by_genome(f, all_proteins, dna_distinct_htable,
                                  distinct_pseqids, dna_cds_index)
        loci_results[results[0]] = results[1]

    # get contig length for all genomes
    contigs_lengths = {inputs_basenames[f]: fao.sequences_lengths(f)
                       for f in fasta_files}

    # determine size mode for each locus
    print('\nDetermining sequence length mode for all loci...', end='')
    loci_modes = {}
    for file in loci_files:
        alleles_sizes = list(fao.sequences_lengths(file).values())
        loci_modes[loci_basenames[file]] = [sm.determine_mode(alleles_sizes)[0], alleles_sizes]

    print('done.')

    # process results per genome and per locus
    classification_inputs = []
    for locus, hits in loci_results.items():
        genomes_matches, representatives_info = hits

        # get locus length mode
        locus_mode = loci_modes[locus]

        # import file with locus classifications
        locus_results_file = fo.join_paths(preprocess_dir, [locus+'_results'])
        
        classification_inputs.append([locus, genomes_matches,
                                      representatives_info,
                                      inv_map, contigs_lengths,
                                      locus_results_file, locus_mode,
                                      temp_directory, size_threshold,
                                      classify_alleles])

    class_results = mo.map_async_parallelizer(classification_inputs,
                                              mo.function_helper,
                                              cpu_cores,
                                              show_progress=True)

    inf_matches = 0
    excluded = []
    for r in class_results:
        locus = list(r.keys())[0]
        results = r[locus]
        inf_matches += results[2]
        excluded.extend(results[3])
        loci_modes[locus] = results[1]

    #######################################
    # Determine representatives iteratively
    # convert to set, some ids might be duplicated and operations are faster with set
    excluded = set(excluded)
    print('\nExcluded: {0}'.format(len(excluded)))
    # BLASTp to align representatives against remaining sequences
    # create file with remaining sequences
    remaining_ids = [rec.id
                     for rec in SeqIO.parse(unique_pfasta, 'fasta')
                     if rec.id not in excluded]
    reps_ids = [rec.id for rec in SeqIO.parse(concat_reps, 'fasta')]
    print('Remaining: {0}'.format(len(remaining_ids)))

    # create directory to store temp results
    iterative_rep_dir = fo.join_paths(temp_directory, ['6_iterative_reps'])
    fo.create_directory(iterative_rep_dir)

    # create Fasta file with protein sequences of reminaing alleles and
    # loci representatives
    new_prot_file = fo.join_paths(iterative_rep_dir, ['iterative_proteins.fasta'])
    fo.concatenate_files([fo.join_paths(preprocess_dir, ['protein.fasta']), concat_reps],
                          new_prot_file)
    prot_index = SeqIO.index(new_prot_file, 'fasta')

    iteration = 1
    exausted = False
    # add schema representatives to get self-score
    #remaining_ids += reps_ids
    while exausted is False:
        ######################################################
        # Need to remove schema representatives from ids!!!! #
        ######################################################
        remaining_seqs_file = fo.join_paths(iterative_rep_dir,
                                            ['remaining_prots_iter{0}.fasta'.format(iteration)])
        # create Fasta with sequences that were not classified
        fao.get_sequences_by_id(prot_index, remaining_ids+reps_ids,
                                remaining_seqs_file, limit=20000)

        # change identifiers to shorten and avoid BLASTp error?
        blast_db = fo.join_paths(iterative_rep_dir, ['blastdb_iter{0}'.format(iteration)])
        db_stderr = bw.make_blast_db(makeblastdb_path, remaining_seqs_file,
                                     blast_db, 'prot')

        # BLAST representatives against remaining sequences
        # iterative process until the process does not detect new representatives
        # pre-compute self-score for representatives?
        print('Representative sets to BLAST against remaining '
              'sequences: {0}\n'.format(len(protein_repfiles)))

        # create BLASTp inputs
        output_files = []
        blast_inputs = []
        for file in protein_repfiles:
            locus_base = fo.file_basename(file, False)
            outfile = fo.join_paths(iterative_rep_dir,
                                    [locus_base+'_blast_results_iter{0}.tsv'.format(iteration)])
            output_files.append(outfile)
            # create file with ids to BLAST against
            # get rep ids
            current_repids = [rec.id
                              for rec in SeqIO.parse(file, 'fasta')]
            # add remaining ids
            current_repids.extend(remaining_ids)
            # save file with ids
            current_file = fo.join_paths(iterative_rep_dir,
                                         [locus_base+'_ids_iter{0}.txt'.format(iteration)])
            with open(current_file, 'w') as outkadgk:
                joined_ids = '\n'.join(current_repids)
                outkadgk.write(joined_ids)

            blast_inputs.append([blastp_path, blast_db, file, outfile,
                                 1, 1, current_file, bw.run_blast])

        blastp_results = mo.map_async_parallelizer(blast_inputs,
                                                   mo.function_helper,
                                                   cpu_cores,
                                                   show_progress=True)

        # prepare results for allele calling
        high_scoring = []
        for f in output_files:
            outfile = select_high_score(f, blast_score_ratio,
                                        None, iterative_rep_dir)
            if outfile is not None:
                high_scoring.append(outfile)

        loci_results = {}
        for f in high_scoring:
            results = group_by_genome(f, all_proteins, dna_distinct_htable,
                                      distinct_pseqids, dna_cds_index)
            loci_results[results[0]] = results[1]

        print('\nLoci with new representative candidates: '
              '{0}'.format(len(loci_results)))

        if len(loci_results) == 0:
            exausted = True
            continue

        # process results per genome and per locus
        classification_inputs = []
        for locus, hits in loci_results.items():
            genomes_matches, representatives_info = hits

            # get locus length mode
            locus_mode = loci_modes[locus]

            # import file with locus classifications
            locus_results_file = fo.join_paths(preprocess_dir, [locus+'_results'])

            classification_inputs.append([locus, genomes_matches,
                                          representatives_info,
                                          inv_map, contigs_lengths,
                                          locus_results_file, locus_mode,
                                          temp_directory, size_threshold,
                                          classify_alleles])

        class_results = mo.map_async_parallelizer(classification_inputs,
                                                  mo.function_helper,
                                                  cpu_cores,
                                                  show_progress=True)

        excluded = []
        representative_candidates = {}
        for r in class_results:
            locus = list(r.keys())[0]
            results = r[locus]
            inf_matches += results[2]
            excluded.extend(results[3])
            loci_modes[locus] = results[1]
            representative_candidates[locus] = results[4]

        print('Classified {0} alleles.'.format(len(excluded)))
        # exclude sequences that were classified
        remaining_ids = list(set(remaining_ids) - set(excluded))

        print('Selecting representatives for next iteration.')
        # determine representatives and then BLAST new representatives
        # against remaining again. If only one representative candidate
        # there is no need to determine which are the best representatives!
        total_selected = 0
        representatives = {}
        # select representatives for loci that have representative candidates
        for k, v in representative_candidates.items():
            if len(v) == 1:
                representatives.setdefault(k, []).append(v[0][1])
            # select for loci that have multiple candidates
            elif len(v) > 1:
                # create Fasta file with all candidates
                reps_candidates_ids = [s[1] for s in v]
                tmp_file = fo.join_paths(iterative_rep_dir, ['{0}_reps_selection_iter{1}.fasta'.format(k, iteration)])
                fao.get_sequences_by_id(prot_index, reps_candidates_ids, tmp_file)
                # BLAST all against all
                # change identifiers to shorten and avoid BLASTp error?
                blast_db = fo.join_paths(iterative_rep_dir,
                                         ['{0}_reps_iter{1}_blastdb'.format(k, iteration)])
                db_stderr = bw.make_blast_db(makeblastdb_path, tmp_file,
                                             blast_db, 'prot')

                outfile = fo.join_paths(iterative_rep_dir,
                                        ['{0}_reps_iter{1}_blast_results'.format(k, iteration)])
                bw.run_blast(blastp_path, blast_db, tmp_file, outfile)

                blast_results = fo.read_tabular(outfile)
                self_scores = {r[0]: float(r[5]) for r in blast_results if r[0] == r[4]}

                # filter self-alignment lines
                filtered_results = [r for r in blast_results if r[0] != r[4]]

                # sort based on input order
                sorted_results = sorted(filtered_results,
                                        key= lambda x: x[4].split('-protein')[0])

                # compute BSR
                for r in sorted_results:
                    r.append(float(r[-1])/self_scores[r[0]])

                to_remove = []
                for r in sorted_results:
                    if r[0] not in to_remove:
                        if r[-1] >= blast_score_ratio+0.1:
                            to_remove.append(r[4])

                selected = [r[0] for r in sorted_results if r[0] not in to_remove]
                selected = list(set(selected))
                total_selected += len(selected)
                if len(selected) > 0:
                    representatives.setdefault(k, []).extend(selected)

                # remove candidates that were not selected
                remaining_ids = list(set(remaining_ids) - set(to_remove))

        print('Remaining unclassified alleles: {0}'.format(len(remaining_ids)))

        # stop iterating if it is not possible to identify representatives
        if len(representatives) == 0:
            exausted = True
        else:
            iteration += 1
            # create files with representative sequences
            reps_ids = []
            protein_repfiles = []
            for k, v in representatives.items():
                reps_ids += v
                rep_file = fo.join_paths(iterative_rep_dir,
                                         ['{0}_reps_iter{1}.fasta'.format(k, iteration)])
                fao.get_sequences_by_id(prot_index, v, rep_file)
                protein_repfiles.append(rep_file)

    end = time.time()
    delta = end - start
    print('\n', delta/60)
    # create inputs to feed into allele calling function!

    # BLAST schema representatives against remaining sequences
    # then select high-BSR results
    # allele calling function will classify matches and return list of representative candidates
    # Remove sequences excluded by allele calling function
    # BLAST representative candidates and check if that allows to detect more alleles for some loci
    # stop when it is not possible to detect more representative candidates
    # in the end BLAST representative candidates for a locus and determine which should be added? (some representatives might have high BSR)
    # or BLAST detected representatives before/after (better after...) allele calling just to determine if both should be used for BLAST?

    # count number of cases or each locus
    exc = 0
    # we will get more NIPH/EM classifications because the current implementation stops when it detects an exact match
    # it does not perform BLASTp to detect high scoring hits that can change the classification to NIPH
    niphem = 0
    niph = 0
    asm = 0
    alm = 0
    inf = 0
    lotsc = 0
    plot3 = 0
    plot5 = 0
    total_niph = 0
    total_niphem = 0
    total_plot3 = 0
    total_plot5 = 0
    total_asm = 0
    total_alm = 0
    # get total number of classified CDSs! NIPHEM and NIPH encompass multiple CDSs
    ################################################################
    # there are a lot of CDSs that are not being classified!
    # probably losing some sequences along the way, need to check for that
    # the new algorithm classifies more or less CDSs per genome?
    total_cds = 0
    for file in classification_files:
        locus_results = fo.pickle_loader(file)
        total_cds += sum([len(c)-1 for g, c in locus_results.items()])
        all_classifications = [c[0] for g, c in locus_results.items()]
        exc += sum([1 for c in all_classifications if 'EXC' in c])
        niphem += all_classifications.count('NIPHEM')
        total_niphem += sum([len(c)-1 for g, c in locus_results.items() if c[0] == 'NIPHEM'])
        niph += all_classifications.count('NIPH')
        total_niph += sum([len(c)-1 for g, c in locus_results.items() if c[0] == 'NIPH'])
        inf += all_classifications.count('INF')
        asm += all_classifications.count('ASM')
        total_asm += sum([len([r for r in c[1:] if r[1] == 'ASM']) for g, c in locus_results.items()])
        alm += all_classifications.count('ALM')
        lotsc += all_classifications.count('LOTSC')
        plot3 += all_classifications.count('PLOT3')
        plot5 += all_classifications.count('PLOT5')

    print('')
    print('exc: {}\n'
          'niphem: {}\n'
          'niph: {}\n'
          'inf: {}\n'
          'asm: {}\n'
          'alm: {}\n'
          'lotsc: {}\n'
          'plot3: {}\n'
          'plot5: {}'.format(exc, niphem, niph,
                             inf, asm, alm,
                             lotsc, plot3, plot5))

    print('Classified a total of {0} CDSs.'.format(total_cds))
    print(inf_matches)

    return classification_files


def create_outputs(classification_files, inv_map, output_directory):
    """
    """

    # create output folder
    results_dir = fo.join_paths(output_directory,
                                ['results_{0}'.format(str(time.strftime("%Y%m%dT%H%M%S")))])
    fo.create_directory(results_dir)

    # create results_alleles.tsv
    total_inputs = len(inv_map)
    file_cols = []
    inputs_col = ['FILE'] + [inv_map[i] for i in range(1, total_inputs+1)]
    file_cols.append(inputs_col)
    for file in classification_files:
        locus = fo.get_locus_id(file)
        locus_results = fo.pickle_loader(file)
        col = [locus]
        col += [locus_results[i][0] if i in locus_results else 'LNF'
                for i in range(1, total_inputs+1)]

        file_cols.append(col)

    lines_generator = zip(*file_cols)

    results_alleles_outfile = fo.join_paths(results_dir, ['results_alleles.tsv'])
    with open(results_alleles_outfile, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for row in lines_generator:
            writer.writerow(row)

    # create results_statsitics.tsv
    class_counts = {i: [] for i in range(1, total_inputs+1)}
    for file in classification_files:
        locus = fo.get_locus_id(file)
        locus_results = fo.pickle_loader(file)

        for i in range(1, total_inputs+1):
            if i in locus_results:
                class_counts[i].append(locus_results[i][0])
            else:
                class_counts[i].append('LNF')

    total_counts = {i: Counter(v).most_common() for i, v in class_counts.items()}
    classes = ['EXC', 'INF', 'LNF', 'PLOT3', 'PLOT5', 'LOTSC', 'NIPH', 'NIPHEM', 'ALM', 'ASM']
    template_counts = {i: {c: 0 for c in classes} for i in total_counts}

    for k, v in total_counts.items():
        for r in v:
            template_counts[k][r[0]] = r[1]

    final_counts = {inv_map[i]: v for i, v in template_counts.items()}

    header = ['FILE'] + classes
    lines = [header]
    for k, v in final_counts.items():
        new_line = [k] + [str(v[c]) for c in classes]
        lines.append(new_line)

    outlines = ['\t'.join(l) for l in lines]

    results_statistics_outfile = fo.join_paths(results_dir, ['results_statistics.tsv'])
    with open(results_statistics_outfile, 'w') as outfile:
        text = '\n'.join(outlines)
        outfile.write(text+'\n')

# implement mode that only detects exact matches (ultrafast)
# implement mode that infers new alleles and that only adds inferred alleles to schema if requested
#######################################################################################
# The protein exact matches, that are not detected at DNA level, are inferred alleles!?
# They have BSR=1, but the DNA sequences are different than the ones in the schema.
# We can classify them as inferred, but we do not need to check if they are representatives
# because the BSR=1. We can classify as INF and it might be altered later
# The protein exact matches are not being properly counted? I think I am not addin/counting the
# number of DNA alleles that match those protein exact matches
# need to make sure that protein exact matches that are inferred do not count all as INF, only one can count as INF
# other genomes with same protein/CDS need to have EXC
##########################
# Check all steps that remove sequences/seqids to verify that we are not removing sequences that
# should not be removed or keeping some that we should not!!!
def main(input_files, schema_directory, output_directory, ptf_path,
         blast_score_ratio, minimum_length, translation_table,
         size_threshold, word_size, window_size, clustering_sim,
         representative_filter, intra_filter, cpu_cores, blast_path,
         cds_input, prodigal_mode, only_exact, no_cleanup):

    print('Prodigal training file: {0}'.format(ptf_path))
    print('CPU cores: {0}'.format(cpu_cores))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum sequence length: {0}'.format(minimum_length))
    print('Size threshold: {0}'.format(size_threshold))
    print('Word size: {0}'.format(word_size))
    print('Window size: {0}'.format(window_size))
    print('Clustering similarity: {0}'.format(clustering_sim))

    results = allele_calling(input_files, schema_directory, output_directory,
                             ptf_path, blast_score_ratio, minimum_length,
                             translation_table, size_threshold, word_size,
                             window_size, clustering_sim, representative_filter,
                             intra_filter, cpu_cores, blast_path,
                             prodigal_mode, cds_input, only_exact)

    # remove temporary files
#    if no_cleanup is False:
#        fo.delete_directory(results[1])


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', nargs='?', type=str,
                        required=True, dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line.')

    parser.add_argument('-g', '--schema-directory', type=str,
                        required=True, dest='schema_directory',
                        help='')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Output directory where the process will store '
                             'intermediate files and create the schema\'s '
                             'directory.')

    parser.add_argument('--ptf', '--training-file', type=str,
                        required=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('--bsr', '--blast-score-ratio', type=float,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', '--minimum-length', type=int,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length value. Coding sequences '
                             'shorter than this value are excluded.')

    parser.add_argument('--t', '--translation-table', type=int,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences.')

    parser.add_argument('--st', '--size-threshold', type=float,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. Added to the '
                             'schema\'s config file and used to identify '
                             'alleles with a length value that deviates from '
                             'the locus length mode during the allele calling '
                             'process.')

    parser.add_argument('--cpu', '--cpu-cores', type=int,
                        required=False, default=1, dest='cpu_cores',
                        help='Number of CPU cores that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores).')

    parser.add_argument('--b', '--blast-path', type=str,
                        required=False, default='', dest='blast_path',
                        help='Path to the BLAST executables.')

    parser.add_argument('--pm', '--prodigal-mode', required=False,
                        choices=['single', 'meta'],
                        default='single', dest='prodigal_mode',
                        help='Prodigal running mode.')

    parser.add_argument('--CDS', required=False, action='store_true',
                        dest='cds_input',
                        help='If provided, input is a single or several FASTA '
                             'files with coding sequences.')

    parser.add_argument('--only-exact', required=False, action='store_true',
                        dest='only_exact',
                        help='If provided, the process will only determine '
                             'exact matches.')

    parser.add_argument('--no-cleanup', required=False, action='store_true',
                        dest='no_cleanup',
                        help='If provided, intermediate files generated '
                             'during process execution are not removed at '
                             'the end.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main()
