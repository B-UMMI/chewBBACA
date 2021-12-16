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
import csv
import sys
import argparse
from collections import Counter

from Bio import SeqIO

try:
    from utils import (ParalogPrunning,
                       constants as ct,
                       blast_wrapper as bw,
                       core_functions as cf,
                       file_operations as fo,
                       fasta_operations as fao,
                       process_datetime as pdt,
                       sequence_manipulation as sm,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (ParalogPrunning,
                                 constants as ct,
                                 blast_wrapper as bw,
                                 core_functions as cf,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 process_datetime as pdt,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


# import module to determine variable size
import get_varSize_deep as gs


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

    results = mo.map_async_parallelizer(translation_inputs,
                                        mo.function_helper,
                                        cpu_cores,
                                        show_progress=True)

    protein_files = {r[0]: r[1] for r in results}

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


def update_classification(genome_id, locus_results, match_info,
                          multiple_classification='NIPH'):
    """
    """

    # add data about match
    locus_results.setdefault(genome_id, [match_info[3]]).append(match_info)
    # change classification if genome has multiple matches for same locus
    if len(locus_results[genome_id]) > 2:
        locus_results[genome_id][0] = multiple_classification

    return locus_results


# locus = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed/GCA-000007265-protein1.fasta'
# results_file = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/test_allelecall/temp/3_cds_preprocess/GCA-000007265-protein1_results'
# locus_classifications = fo.pickle_loader(results_file)
# fasta_file = locus
# distinct_table = dna_distinct_htable
def dna_exact_matches(fasta_file, distinct_table, locus_classifications):
    """
    """

    # create generator to get Fasta sequences
    seq_generator = SeqIO.parse(fasta_file, 'fasta')

    matches_ids = []
    total_matches = 0
    for rec in seq_generator:
        seqid = (rec.id).split('_')[-1]
        sequence = str(rec.seq.upper())
        sequence_hash = im.hash_sequence(sequence)
        # file cannot have duplicated alleles or it will
        # classify based on same DNA sequence twice
        if sequence_hash in distinct_table:
            current_matches = distinct_table[sequence_hash]
            # exclude element in index 0, it is the seqid chosen as
            # representative
            total_matches += len(current_matches) - 1
            for m in current_matches[1:]:
                # for DNA exact matches, the allele ID,
                # the seqid chosen as repsentative during sequence deduplication,
                # the hash of the schema allele,
                # the classification (EXC) and
                # the BSR value (1.0) are stored
                match_info = (seqid, current_matches[0],
                              sequence_hash, 'EXC', 1.0)
                locus_classifications = update_classification(m, locus_classifications,
                                                              match_info, 'NIPHEM')

            matches_ids.append(current_matches[0])

    return [locus_classifications, matches_ids, total_matches]


def protein_exact_matches(fasta_file, distinct_prot_table,
                          distinct_dna_table, locus_classifications,
                          dna_index):
    """
    """

    seq_generator = SeqIO.parse(fasta_file, 'fasta')

    seen_dna = []
    seen_prot = []
    total_cds = 0
    total_prots = 0
    exact_prot_hashes = []
    for rec in seq_generator:
        protid = rec.id
        protein = str(rec.seq.upper())
        prot_hash = im.hash_sequence(protein)
        # do not proceed if distinct protein sequence has already been seen
        # different alleles might code for the same protein
        if prot_hash in distinct_prot_table and prot_hash not in seen_prot:
            # get protids for distinct DNA CDSs
            protids = distinct_prot_table[prot_hash][1:]
            total_prots += len(protids)
            exact_prot_hashes.extend(protids)
            # for each distinct CDS that codes for the protein
            for p in protids:
                cds = str(dna_index.get(p).seq)
                cds_hash = im.hash_sequence(cds)
                # get IDs of genomes that contain the CDS
                genome_ids = distinct_dna_table[cds_hash]
                total_cds += len(genome_ids) - 1
                # for each genome ID that contains the CDS
                for gid in genome_ids[1:]:
                    # first time seeing CDS
                    if cds_hash not in seen_dna:
                        current_class = 'INF'
                        seen_dna.append(cds_hash)
                    else:
                        current_class = 'EXCP'
                    # for protein exact matches, the seqid of the translated allele,
                    # the seqid of the protein chosen as representative during sequence deduplication,
                    # the hash of the schema allele,
                    # the classification (INF or EXCP) and
                    # the BSR value (1.0) are stored
                    match_info = (protid, p, cds_hash, current_class, 1.0)
                    locus_classifications = update_classification(gid, locus_classifications,
                                                                  match_info)

            seen_prot.append(prot_hash)

    return [locus_classifications, exact_prot_hashes, total_prots, total_cds]


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
                     locus_mode, temp_directory, size_threshold,
                     blast_score_ratio):
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

            ##### This might lead to NIPH if a genome has several exact matches after the CDS
            # is detected in a previous genome!!!
            if target_dna_hash in seen_dna:
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'EXC', 1.0))
                inf_matches += 1
                continue

            # check if CDS DNA sequence is not in schema but matches a translated allele
            if target_prot_hash in seen_prot:
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'INF', 1.0))
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
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'LOTSC', bsr))
                continue
            # check if PLOT
            elif contig_left_rest < rep_left_rest:
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'PLOT3', bsr))
                continue
            elif contig_right_rest < rep_right_rest:
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'PLOT5', bsr))
                continue

            # check if ASM or ALM
            # we only need to evaluate one of the genomes, if they are ASM/ALM we can classify all of them as the same!
            if target_dna_len < (locus_mode[0]-(locus_mode[0])*size_threshold):
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'ASM', bsr))
                continue
            elif target_dna_len > (locus_mode[0]+(locus_mode[0])*size_threshold):
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                      target_dna_hash, 'ALM', bsr))
                continue

            # add INF
            locus_results = update_classification(genome, locus_results,
                                                  (rep_alleleid, target_seqid,
                                                  target_dna_hash, 'INF', bsr))

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


# file = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/test_allelecall/temp/3_cds_preprocess/GCA-000007265-protein1_results'
def assign_ids(classification_files, schema_directory):
    """
    """

    # assign allele identifiers
    new_alleles = {}
    for locus, results in classification_files.items():
        locus_id = fo.get_locus_id(locus)
        # get loci records
        records = fao.import_sequences(locus)
        # get id of max record
        max_alleleid = max([int(rec.split('_')[-1]) for rec in records])
        locus_results = fo.pickle_loader(results)
        # iterate and assign allele identifiers
        seen = {}
        # sort to get INF classifications first
        sorted_r = sorted(locus_results.items(),
                          key=lambda x: x[1][0] == 'INF',
                          reverse=True)
        for k in sorted_r:
            genome_id = k[0]
            current_results = k[1]
            cds_hash = current_results[1][2]
            if current_results[0] == 'EXC':
                locus_results[genome_id].append(current_results[1][0])
            elif current_results[0] == 'EXCP':
                # this should be True every single case, but INF cases might be converted
                # to NIPH due to a second match after being used to detect EXCP
                if cds_hash in seen:
                    locus_results[genome_id].append(seen[cds_hash])
                # the exact match might not be found because the INF that
                # added the match was converted to NIPH
                # we need to add the new allele based on the EXCP case...
                else:
                    max_alleleid += 1
                    locus_results[genome_id].append('INF-{0}'.format(max_alleleid))
                    seen[cds_hash] = str(max_alleleid)
                    new_alleles.setdefault(locus, []).append(current_results[1][1])
            elif current_results[0] == 'INF':
                if cds_hash not in seen:
                    max_alleleid += 1
                    locus_results[genome_id].append('INF-{0}'.format(max_alleleid))
                    seen[cds_hash] = str(max_alleleid)
                    new_alleles.setdefault(locus, []).append(current_results[1][1])
                else:
                    print('INF case was already seen!')
                    locus_results[genome_id].append(seen[cds_hash])

        # save updated info
        fo.pickle_dumper(locus_results, results)

    return new_alleles


def write_logfile(start_time, end_time, total_inputs,
                  total_loci, cpu_cores, blast_score_ratio,
                  output_directory):
    """
    """

    start_time_str = pdt.datetime_str(start_time,
                                      date_format='%Y-%m-%d - %H:%M:%S')
    
    end_time_str = pdt.datetime_str(end_time,
                                    date_format='%Y-%m-%d - %H:%M:%S')


    log_outfile = fo.join_paths(output_directory, ['logging_info.txt'])
    with open(log_outfile, 'w') as outfile:
        outfile.write('Started Script at: {0}'.format(start_time_str))
        outfile.write('\nFinished Script at: {0}'.format(end_time_str))
        outfile.write('\nNumber of genomes: {0}'.format(total_inputs))
        outfile.write('\nNumber of loci: {0}'.format(total_loci))
        outfile.write('\nUsed this number of CPU cores: {0}'.format(cpu_cores))
        outfile.write('\nUsed a bsr of: {0}\n'.format(blast_score_ratio))

    return log_outfile


######
# Locus file with results for testing
#results_test = 'GCA-000007265-protein1_results'


#input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids.txt'
input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids32.txt'
#input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids320.txt'
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
schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed'
#schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema_called_320/schema_seed'
#schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/spyogenes_schema_processed'
def allele_calling(input_files, schema_directory, output_directory, ptf_path,
                   blast_score_ratio, minimum_length, translation_table,
                   size_threshold, word_size, window_size, clustering_sim,
                   cpu_cores, blast_path, prodigal_mode, cds_input,
                   only_exact):
    """
    """

    # define directory for temporary files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # read file with paths to input files
    fasta_files = fo.read_lines(input_files, strip=True)

################
# Genomes are not being classified in the order they were sorted
# leading to incorrect order of INF and EXC in the matrix!

    # sort paths to FASTA files
    fasta_files = im.sort_data(fasta_files, sort_key=str.lower)

    # map full paths to basename
    inputs_basenames = fo.mapping_function(fasta_files,
                                           fo.file_basename, [False])

    # map input identifiers to integers
    # use the mapped integers to refer to each input
    # this reduces memory usage compared to using the string identifiers
    map_ids = im.integer_mapping(list(inputs_basenames.values()))
    inv_map = im.invert_dictionary(map_ids)

    # inputs are genome assemblies
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
    # inputs are Fasta files with the predicted CDSs
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

    # get list of loci files
    print('Getting list of loci...', end='')
    loci_files = fo.listdir_fullpath(schema_directory,
                                     substring_filter='.fasta')

    # get mapping between locus file path and locus identifier
    loci_basenames = fo.mapping_function(loci_files, fo.get_locus_id, [])
    print('schema has {0} loci.'.format(len(loci_files)))

    # create files with empty results data structure
    classification_files = {}
    empty_results = {}
    for file in loci_files:
        pickle_out = fo.join_paths(preprocess_dir,
                                   [loci_basenames[file]+'_results'])

        # create file with empty results structure
        fo.pickle_dumper(empty_results, pickle_out)
        classification_files[file] = pickle_out

    print('Finding DNA exact matches...', end='')
    # find exact DNA matches
    # Add multiprocessing?
    # use psutil library to get free memory and decide maximum number of processes to spawn?
    dna_exact_hits = 0
    dna_matches_ids = []
    for locus, results_file in classification_files.items():
        locus_classifications = fo.pickle_loader(results_file)
        em_results = dna_exact_matches(locus, dna_distinct_htable, locus_classifications)
        # save classifications
        fo.pickle_dumper(em_results[0], results_file)
        # extend list of matched seqids
        dna_matches_ids.extend(em_results[1])
        dna_exact_hits += em_results[2]

    print('found {0} exact matches (matching {1} alleles).'
          ''.format(dna_exact_hits, len(dna_matches_ids)))

    if only_exact is True:
        # return classification files for creation of output files
        return [classification_files, inv_map, []]

    # create new Fasta file without the DNA sequences that were exact matches
    dna_index = SeqIO.index(distinct_file, 'fasta')
    unique_fasta = fo.join_paths(preprocess_dir, ['dna_non_exact.fasta'])
    # this step is taking too long? Is it also using too much memory?
    # Optimize to run faster with a great number of input assemblies!
    # multiprocessing?
    total_selected = fao.exclude_sequences_by_id(dna_index,
                                                 dna_matches_ids,
                                                 unique_fasta)

    # verify if having these identifiers in memory is problematic for huge datasets
    selected_ids = (rec for rec in dna_index if rec not in dna_matches_ids)
    selected_ids = list(selected_ids)

    # translate DNA sequences and identify duplicates
    print('\nTranslating {0} DNA sequences...'.format(total_selected))

    # this step excludes small sequences
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

    # identify exact matches at protein level
    # exact matches are new alleles that can be added to the schema
    print('\nFinding protein exact matches...', end='')
    exc_cds = 0
    exc_prot = 0
    exact_phashes = []
    for locus, pfile in protein_files.items():
        results_file = classification_files[locus]
        locus_classifications = fo.pickle_loader(results_file)
        em_results = protein_exact_matches(pfile, distinct_pseqids,
                                           dna_distinct_htable, locus_classifications,
                                           dna_index)

        fo.pickle_dumper(em_results[0], results_file)
        exact_phashes.extend(em_results[1])
        exc_prot += em_results[2]
        exc_cds += em_results[3]

    # need to determine number of distinct proteins to print correct info!!!
    print('found {0} protein exact matches (matching {1} distinct proteins).'
          ''.format(exc_prot, len(exact_phashes)))

    print('Protein matches correspond to {0} DNA matches.'.format(exc_cds))

    # create new Fasta file without the Protein sequences that were exact matches
    unique_pfasta = fo.join_paths(preprocess_dir, ['protein_non_exact.fasta'])
    # create protein file index
    protein_index = SeqIO.index(ds_results[1], 'fasta')
    total_selected = fao.exclude_sequences_by_id(protein_index, exact_phashes, unique_pfasta)

    # translate schema representatives
    print('Translating schema representatives...')
    rep_dir = os.path.join(schema_directory, 'short')
    rep_list = fo.listdir_fullpath(rep_dir, '.fasta')
    protein_repfiles = translate_fastas(rep_list, protein_dir,
                                        translation_table, cpu_cores)

    # cluster protein sequences
    proteins = fao.import_sequences(unique_pfasta)

    # create directory to store clustering data
    clustering_dir = fo.join_paths(temp_directory, ['5_clustering'])
    fo.create_directory(clustering_dir)

    # BLASTp clusters step
    blastp_path = os.path.join(blast_path, ct.BLASTP_ALIAS)
    makeblastdb_path = os.path.join(blast_path, ct.MAKEBLASTDB_ALIAS)

    # concatenate all representative
    concat_reps = os.path.join(protein_dir, 'concat_reps.fasta')
    fo.concatenate_files(protein_repfiles.values(), concat_reps)

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
        all_prots = os.path.join(blasting_dir, 'all_prots.fasta')
        
        # create Fasta file with remaining proteins and representatives
        fo.concatenate_files([unique_pfasta, concat_reps], all_prots)

        all_proteins = fao.import_sequences(all_prots)

        blast_results, ids_dict = cf.blast_clusters(clusters, all_proteins,
                                                    blasting_dir, blastp_path,
                                                    makeblastdb_path, cpu_cores,
                                                    'blast', True)

        blast_files = im.flatten_list(blast_results)

    # for ASM and ALM classifications we only need to determine it for the one of the ids associated to that protein
    # the other proteins have the same length and the classification will be the same
    # for PLOT cases we need to check all cases

    # add id of representative CDS, get hash for that sequence and then open file with genome CDS to get CDS coordinates
    # hash file with DNA sequences to get DNA sequence to hash?

    # add hash of inferred alleles to dict with results per locus
    # we can check if hash is already in the dict to avoid calculating everyhting, effectivelly detecting new exact matches!
    # we can store the DNA hash and the protein hash. DNA hash is an exact match and protein hash is a new allele.

    # group files for representatives of the same locus
    loci_results = {}
    for f in blast_files:
        # create function to accept regex pattern and extract string?
        locus_rep = ids_dict[im.match_regex(f, r'seq_[0-9]+')]
        locus_id = im.match_regex(locus_rep, r'^.+-protein[0-9]+')
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
                                      blast_score_ratio,
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

    # Determine representatives iteratively
    # convert to set, some ids might be duplicated
    excluded = set(excluded)
    print('\nExcluded: {0}'.format(len(excluded)))

    # get seqids of remaining unclassified sequences
    unclassified_ids = [rec.id
                        for rec in SeqIO.parse(unique_pfasta, 'fasta')
                        if rec.id not in excluded]
    print('Remaining: {0}'.format(len(unclassified_ids)))

    # get seqids of schema representatives
    reps_ids = [rec.id for rec in SeqIO.parse(concat_reps, 'fasta')]
    print('Schema has a total of {0} representative alleles.'
          ''.format(len(reps_ids)))

    # create directory to store data for each iteration
    iterative_rep_dir = fo.join_paths(temp_directory, ['6_iterative_reps'])
    fo.create_directory(iterative_rep_dir)

    # create Fasta file with remaining unclassified sequences
    # loci representatives
    new_prot_file = fo.join_paths(iterative_rep_dir, ['iterative_proteins.fasta'])
    ### Using a Fasta file with more than unclassified alleles!
    fo.concatenate_files([fo.join_paths(preprocess_dir, ['protein.fasta']), concat_reps],
                          new_prot_file)
    prot_index = SeqIO.index(new_prot_file, 'fasta')

    protein_repfiles = list(protein_repfiles.values())

    iteration = 1
    exausted = False
    while exausted is False:
        ######################################################
        # Need to remove schema representatives from ids!!!! #
        ######################################################
        remaining_seqs_file = fo.join_paths(iterative_rep_dir,
                                            ['remaining_prots_iter{0}.fasta'.format(iteration)])
        # create Fasta with sequences that were not classified
        fao.get_sequences_by_id(prot_index, unclassified_ids+reps_ids,
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
            current_repids.extend(unclassified_ids)
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
                                          blast_score_ratio,
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
        unclassified_ids = list(set(unclassified_ids) - set(excluded))

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
                unclassified_ids = list(set(unclassified_ids) - set(to_remove))

        print('Remaining unclassified alleles: {0}'.format(len(unclassified_ids)))

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

    exc = 0
    # we will get more NIPH/EM classifications because the current
    # implementation stops when it detects an exact match it does
    # not perform BLASTp to detect high scoring hits that can change
    # the classification to NIPH
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
    # get total number of classified CDSs
    total_cds = 0
    for file in classification_files.values():
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

    return [classification_files, inv_map]


start_time = pdt.get_datetime()
end_time = pdt.get_datetime()
def create_outputs(classification_files, inv_map, output_directory,
                   start_time, end_time, cpu_cores, blast_score_ratio):
    """
    """

    # create output folder
    results_dir = fo.join_paths(output_directory,
                                ['results_{0}'.format(end_time)])
    fo.create_directory(results_dir)

    print('Writing logging_info.txt...', end='')
    write_logfile(start_time, end_time, len(inv_map),
                  len(classification_files), cpu_cores,
                  blast_score_ratio, results_dir)
    print('done.')

    # create results_alleles.tsv
    print('Writing results_alleles.tsv...', end='')
    total_inputs = len(inv_map)
    file_cols = []
    inputs_col = ['FILE'] + [inv_map[i] for i in range(1, total_inputs+1)]
    file_cols.append(inputs_col)
    for file in classification_files.values():
        locus = fo.get_locus_id(file)
        locus_results = fo.pickle_loader(file)
        col = [locus]
        for i in range(1, total_inputs+1):
            if i in locus_results:
                current_result = locus_results[i]
                if current_result[0] in ['EXC', 'EXCP', 'INF']:
                    col.append(current_result[-1])
                else:
                    col.append(current_result[0])
            else:
                col.append('LNF')

        file_cols.append(col)

    lines_generator = zip(*file_cols)

    results_alleles_outfile = fo.join_paths(results_dir, ['results_alleles.tsv'])
    with open(results_alleles_outfile, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        for row in lines_generator:
            writer.writerow(row)

    print('done.')

    # create results_statsitics.tsv
    print('Writing results_statsitics.tsv...', end='')
    class_counts = {i: [] for i in range(1, total_inputs+1)}
    for locus, results in classification_files.items():
        locus_id = fo.get_locus_id(locus)
        locus_results = fo.pickle_loader(results)

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
    
    print('done.')

    # create results_contigsInfo.tsv
    # using too much memory???
    print('Writing results_contigsInfo.tsv...', end='')
    invalid_classes = ['LNF', 'PLOT3', 'PLOT5',
                       'LOTSC', 'NIPH', 'NIPHEM',
                       'ALM', 'ASM']

    file_cols = []
    inputs_col = ['FILE'] + [inv_map[i] for i in range(1, total_inputs+1)]
    file_cols.append(inputs_col)
    for file in classification_files:
        locus = fo.get_locus_id(file)
        locus_results = fo.pickle_loader(file)
        col = [locus]
        col += [locus_results[i][1][2]
                if i in locus_results and locus_results[i][0] not in invalid_classes
                else locus_results.get(i, ['LNF'])[0]
                for i in range(1, total_inputs+1)]

        file_cols.append(col)

    lines_generator = zip(*file_cols)

    # fetch coordinates for each genome
    coordinates_dir = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/test_allelecall/temp/2_cds_extraction'
    coordinates_files = fo.listdir_fullpath(coordinates_dir, 'cds_hash')
    coordinates_files = {fo.file_basename(f, True).split('_cds_hash')[0]: f
                         for f in coordinates_files}

    final_lines = []
    header = list(lines_generator.__next__())
    final_lines.append(header)
    for l in lines_generator:
        genome_id = l[0]
        # open file with coordinates
        coordinates = fo.pickle_loader(coordinates_files[genome_id])
        new_line = [coordinates[c][0] if c not in invalid_classes else c for c in l[1:]]
        new_line = ['{0}&{1}-{2}&{3}'.format(c[0], c[1], c[2], c[4])
                    if c not in invalid_classes else c
                    for c in new_line]
        final_lines.append([genome_id]+new_line)

    #########################
    # Start positions have offset of 1 and do not match start positions reported by previous implementation?
    #########################

    # write file
    results_contigs_outfile = fo.join_paths(results_dir, ['results_contigsInfo.tsv'])
    with open(results_contigs_outfile, 'w') as outfile:
        outlines = ['\t'.join(l) for l in final_lines]
        text = '\n'.join(outlines)
        outfile.write(text+'\n')
        
    print('done.')

    # determine paralogous loci and write RepeatedLoci.txt file
    print('Writing RepeatedLoci.txt...', end='')
    ParalogPrunning.main(results_contigs_outfile, results_dir)
    print('\nResults available in {0}'.format(results_dir))


# need to make sure that protein exact matches that are inferred do not count all as INF, only one can count as INF
# other genomes with same protein/CDS need to have EXC
##########################
# Check all steps that remove sequences/seqids to verify that we are not removing sequences that
# should not be removed or keeping some that we should not!!!
results = [classification_files, inv_map]
def main(input_files, schema_directory, output_directory, ptf_path,
         blast_score_ratio, minimum_length, translation_table,
         size_threshold, word_size, window_size, clustering_sim,
         cpu_cores, blast_path, cds_input, prodigal_mode, only_exact,
         add_inferred, no_cleanup):

    print('Prodigal training file: {0}'.format(ptf_path))
    print('CPU cores: {0}'.format(cpu_cores))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum sequence length: {0}'.format(minimum_length))
    print('Size threshold: {0}'.format(size_threshold))
    print('Word size: {0}'.format(word_size))
    print('Window size: {0}'.format(window_size))
    print('Clustering similarity: {0}'.format(clustering_sim))

    start_time = pdt.get_datetime()

    results = allele_calling(input_files, schema_directory, output_directory,
                             ptf_path, blast_score_ratio, minimum_length,
                             translation_table, size_threshold, word_size,
                             window_size, clustering_sim, cpu_cores, blast_path,
                             prodigal_mode, cds_input, only_exact)

    new_alleles = assign_ids(results[0], schema_directory)

    if add_inferred is True:
        if len(new_alleles) > 0:
            # add inferred alleles to schema
            pass
        else:
            print('No new alleles to add to schema.')

    end_time = pdt.get_datetime()

    # create output files
    create_outputs(results[0], results[1], output_directory,
                   start_time, end_time, cpu_cores, blast_score_ratio)

    # remove temporary files
    if no_cleanup is False:
        fo.delete_directory(fo.join_paths(output_directory, ['temp']))


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
    
    parser.add_argument('--gl', '--genes-list', type=str,
                        required=False, default=False, dest='genes_list',
                        help='Path to a file with the list of genes '
                             'in the schema that the process should '
                             'identify alleles for.')

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

    parser.add_argument('--add-inferred', required=False, action='store_true',
                        dest='add_inferred',
                        help='If provided, the process will add the sequences '
                             'of inferred alleles to the schema.')

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
