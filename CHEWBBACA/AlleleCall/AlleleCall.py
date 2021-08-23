#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables ...

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
import sys
import time
import argparse

from Bio import SeqIO

try:
    from utils import (constants as ct,
                       blast_wrapper as bw,
                       core_functions as cf,
                       gene_prediction as gp,
                       file_operations as fo,
                       fasta_operations as fao,
                       sequence_clustering as sc,
                       sequence_manipulation as sm,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except:
    from CHEWBBACA.utils import (constants as ct,
                                 blast_wrapper as bw,
                                 core_functions as cf,
                                 gene_prediction as gp,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 sequence_clustering as sc,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


# import module to determine variable size
import get_varSize_deep as gs


def get_mode(locus_fasta):
    """
    """

    seq_generator = SeqIO.parse(locus_fasta, 'fasta')
    alleles_sizes = [len(rec.seq)
                     for rec in seq_generator]
    modes = sm.determine_mode(alleles_sizes)

    # only return first mode
    return modes[0]


def exact_matches(locus_file, locus_basename, distinct_table, output_dir, map_ids=False, match_type='EXCDNA'):
    """
    """

    seq_generator = SeqIO.parse(locus_file, 'fasta')
    pickle_out = fo.join_paths(output_dir, [locus_basename+'_results'])

    total = 0
    exact_hashes = []
    exact_matches = {}
    for rec in seq_generator:
        seqid = (rec.id).split('_')[-1]
        sequence = str(rec.seq.upper())
        seq_hash = im.hash_sequence(sequence)
        if seq_hash in distinct_table:
            current_matches = distinct_table[seq_hash]
            if map_ids is not False:
                current_matches = [current_matches[0]] + [map_ids[c.split('-protein')[0]] for c in current_matches[1:]]
            total += len(current_matches) - 1
            for m in current_matches[1:]:
                if m not in exact_matches:
                    exact_matches[m] = ['EXC', (seqid, match_type)]
                else:
                    if all([i[1] == 'EXCDNA' for i in exact_matches[m][1:]]) is True:
                        exact_matches[m][0] = 'NIPHEM'
                        exact_matches[m].append((seqid, match_type))
                    else:
                        exact_matches[m][0] = 'NIPH'
                        exact_matches[m].append((seqid, match_type))
            exact_hashes.append(current_matches[0])

    # update classifications
    if os.path.isfile(pickle_out) is False:
        fo.pickle_dumper(exact_matches, pickle_out)
    else:
        locus_results = fo.pickle_loader(pickle_out)
        for k, v in exact_matches.items():
            locus_results.setdefault(k, [v[0]]).extend(v[1:])
            if len(v[1:]) > 1:
                if all([i[1] == 'EXCDNA' for i in locus_results[k][1:]]) is True:
                    locus_results[k][0] = 'NIPHEM'
                else:
                    locus_results[k][0] = 'NIPH'

        fo.pickle_dumper(locus_results, pickle_out)

    return [exact_hashes, pickle_out, total]


input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids.txt'
#input_files = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/ids32.txt'
output_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/test_allelecall'
ptf_path = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed/Streptococcus_agalactiae.trn'
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
blast_path = '/home/rfm/Software/anaconda3/envs/ns/bin'
prodigal_mode = 'single'
cds_input = False
schema_directory = '/home/rfm/Desktop/rfm/Lab_Software/AlleleCall_tests/sagalactiae32_schema/schema_seed'
def allele_calling(input_files, schema_directory, output_directory, ptf_path,
                   blast_score_ratio, minimum_length, translation_table,
                   size_threshold, word_size, window_size, clustering_sim,
                   representative_filter, intra_filter, cpu_cores, blast_path,
                   prodigal_mode, cds_input):
    """
    """

    # define directory for temporary files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    # read file with paths to input files
    fasta_files = fo.read_lines(input_files, strip=True)

    # sort paths to FASTA files
    fasta_files = im.sort_data(fasta_files, sort_key=lambda x: x.lower())
    inputs_basenames = fo.mapping_function(fasta_files,
                                           fo.file_basename, [False])
    inv_basenames = im.invert_dictionary(inputs_basenames)

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
        cds_extraction_path = os.path.join(temp_directory,
                                           '2_cds_extraction')
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
    start = time.time()
    print('\nRemoving duplicated DNA sequences...', end='')
    distinct_dna_template = 'distinct_cds_{0}.fasta'
    dna_dedup_results = cf.exclude_duplicates(cds_files, preprocess_dir,
                                              cpu_cores, distinct_dna_template,
                                              map_ids, False)

    distinct_htable, distinct_file, repeated = dna_dedup_results
    print('removed {0} sequences.'.format(repeated))
    end = time.time()
    delta = end - start
    print(delta)

    # get size of hash table
    size1 = gs.convert_bytes(distinct_htable, set())
    print(size1)

    # remove small sequences!?
    #####################
    ## Add classification information as genome_id --> classification, BSR e.g.: {4: ['EXC', ('1', 1)]}
    print('Finding DNA exact matches...', end='')
    #####################
    ## Check if reverse complement matches?
    # find DNA exact matches

    # get list of loci files
    loci_files = fo.listdir_fullpath(schema_directory,
                                     substring_filter='.fasta')

    loci_basenames = fo.mapping_function(loci_files, fo.file_basename, [False])

    # determine size mode for each locus
    loci_modes = {}
    for file in loci_files:
        loci_modes[loci_basenames[file]] = get_mode(file)

    # find exact DNA matches
    exc_dna = 0
    exact_hashes = []
    classification_files = []
    for file in loci_files:
        em_results = exact_matches(file, loci_basenames[file],
                                   distinct_htable, preprocess_dir)
        exact_hashes.extend(em_results[0])
        classification_files.append(em_results[1])
        exc_dna += em_results[2]

    print('found {0} exact matches (matching {1} alleles).'.format(exc_dna, len(exact_hashes)))

    # remove DNA sequences that were exact matches from the file with all distinct CDSs
    seq_generator = SeqIO.parse(distinct_file, 'fasta')
    selected_ids = [rec.id
                    for rec in seq_generator
                    if rec.id not in exact_hashes]
    cds_index = SeqIO.index(distinct_file, 'fasta')
    unique_fasta = fo.join_paths(preprocess_dir, ['dna_non_exact.fasta'])
    fao.get_sequences_by_id(cds_index, selected_ids, unique_fasta, limit=20000)

    # translate DNA sequences and identify duplicates
    # sequence translation step
    print('\nTranslating {0} DNA sequences...'.format(len(selected_ids)))

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

    # translate loci files and identify exact matches at protein level
    # exact matches are new alleles that can be added to the schema
    # find exact matches
    # key can be reduced to basename and allele match to allele identifier

    # translate loci files
    protein_dir = os.path.join(temp_directory, '4_protein_dir')
    fo.create_directory(protein_dir)
    protein_files = fao.translate_fastas(loci_files, protein_dir, 11)

    # get the integer identifier of the last allele for all loci
    loci_topid = {}
    for f in loci_files:
        records = SeqIO.parse(f, 'fasta')
        allele_ids = [int((rec.id).split('_')[-1]) for rec in records]
        topid = max(allele_ids)
        loci_topid[loci_basenames[f]] = topid

    print('Finding protein exact matches...', end='')
    exc_prot = 0
    exact_phashes = []
    for file in protein_files:
        em_results = exact_matches(file, loci_basenames[os.path.join(schema_directory, os.path.basename(file)).replace('_protein', '')],
                                   distinct_pseqids, preprocess_dir, map_ids, 'EXCPROT')
        exact_phashes.extend(em_results[0])
        exc_prot += em_results[2]

    print('found {0} exact matches (matching {1} proteins).'.format(exc_prot, len(exact_phashes)))

    # remove Protein sequences that were exact matches from the file with all distinct proteins
    seq_generator = SeqIO.parse(ds_results[1], 'fasta')
    selected_pids = [rec.id
                    for rec in seq_generator
                    if rec.id not in exact_phashes]
    prot_index = SeqIO.index(ds_results[1], 'fasta')
    unique_pfasta = fo.join_paths(preprocess_dir, ['protein_non_exact.fasta'])
    fao.get_sequences_by_id(prot_index, selected_pids, unique_pfasta, limit=20000)

    # translate schema representatives
    rep_dir = os.path.join(schema_directory, 'short')
    rep_list = [os.path.join(rep_dir, f) for f in os.listdir(rep_dir) if f.endswith('.fasta')]
    protein_repfiles = fao.translate_fastas(rep_list, protein_dir, 11)

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

    ##########################################
    # iterate and perform allele calling until the process cannot identify new representatives
    dna_cds_index = SeqIO.index(unique_fasta, 'fasta')
    new_reps = True
    current_iteration = 1
    while new_reps is True:
        # create directory for current iteration
        iteration_directory = fo.join_paths(protein_dir, ['iteration_{0}'.format(current_iteration)])
        fo.create_directory(iteration_directory)

        # create index for representative sequences
        # concatenate all representative
        concat_reps = os.path.join(protein_dir,
                                   'concat_reps_{0}.fasta'.format(current_iteration))
        fo.concatenate_files(protein_repfiles, concat_reps)
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

        print('')
        # import results and determine classifications
        processed = 0
        inf_matches = 0
        representative_candidates = {}
        for f in blast_files:
            # import BLASTp results for single representative cluster
            current_results = fo.read_tabular(f)
            # get allele identifier for the representative
            rep_id = ids_dict[current_results[0][0]]
            rep_alleleid = rep_id.split('_')[-1]
            # get length of representative DNA sequence
            rep_len = (int(current_results[0][3])*3) + 3
            # get locus identifier
            locus = rep_id.split('_')[0]
            # get locus length mode
            locus_mode = loci_modes[locus]
    
            # get BLASTp self-score for representative
            self_score = float([r[5] for r in current_results if r[0] == r[4]][0])
    
            # exclude representative self-alignment
            current_results = [r for r in current_results if r[0] != r[4]]
    
            # determine BSR values
            bsr_values = {s[4]: (float(s[5])/self_score, int(s[1]), int(s[2]))
                          for s in current_results}
    
            # only keep matches above BSR threshold
            high_bsr = {k: v
                        for k, v in bsr_values.items()
                        if v[0] >= blast_score_ratio}
    
            # start classifying high BSR hits
            if len(high_bsr) > 0:
                # instantiate list to store hashes of protein sequences that are classified as inferred
                seen_dna = []
                seen_prot = []
                # import allele calling results for locus
                locus_results_file = fo.join_paths(preprocess_dir, [locus+'_results'])
                locus_results = fo.pickle_loader(locus_results_file)
                for k, v in high_bsr.items():
                    # get CDS string original identifier
                    target_id = ids_dict[k]
                    bsr = v[0]
                    # get Protein hash
                    prot_hash = im.hash_sequence(all_proteins[target_id])
    
                    # careful about offset!!!
                    # determine right match position on representative allele
                    rep_right_pos = rep_len - ((v[2]+1)*3)
                    # determine left match position on representative allele
                    rep_left_pos = ((v[1]-1)*3)
    
                    # use hash to get all CDS identifiers that coded for the same protein
                    prot_seqids = distinct_pseqids[prot_hash][1:]
    
                    # need to get the identifiers for all genomes that had the CDSs!
                    # We are only getting the CDSs for the original DNA CDSs representatives...
    
                    # open file with CDS coordinates for each genome and get coordinates
                    genome_coordinates = {}
                    for g in prot_seqids:
                        # get genome identifier for the representative
                        #genome_strid = g.split('-protein')[0]
                        #genome_id = map_ids[genome_strid]
    
                        # get CDS DNA sequence
                        cds_dna = str(dna_cds_index.get(g).seq)
                        cds_len = len(cds_dna)
                        cds_hash = im.hash_sequence(cds_dna)
    
                        # get ids for all genomes with same CDS as representative
                        all_genomes_with_cds = distinct_htable[cds_hash][1:]
    
                        # check if CDS DNA sequence was previously identified as new allele
                        for og in all_genomes_with_cds:
                            if cds_hash in seen_dna:
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'EXCBLASTp', 1.0))
                                else:
                                    locus_results[og] = ['EXC', (rep_alleleid, 'EXCBLASTp', 1.0)]
                                inf_matches += 1
                                continue
    
                            # check if CDS DNA sequence is not in schema but matches a translated allele
                            if prot_hash in seen_prot:
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'EXCBLASTp', 1.0))
                                else:
                                    locus_results[og] = ['INF', (rep_alleleid, 'EXCBLASTp', 1.0)]
                                    seen_dna.append(cds_hash)
                                inf_matches += 1
                                continue
    
                            current_g = inv_map[og]
                            # there is no exact match in the schema, perform full evaluation
                            # get contig lengths
                            #contigs = SeqIO.parse(inv_basenames[genome_strid], 'fasta')
                            contigs = SeqIO.parse(inv_basenames[current_g], 'fasta')
                            contigs_lengths = {rec.id: len(rec.seq) for rec in contigs}
    
                            # open pickle for genome and get coordinates
                            #genome_cds_file = fo.join_paths(temp_directory, ['2_cds_extraction', genome_strid+'_cds_hash'])
                            genome_cds_file = fo.join_paths(temp_directory, ['2_cds_extraction', current_g+'_cds_hash'])
                            genome_cds_coordinates = fo.pickle_loader(genome_cds_file)
    
                            genome_coordinates[g] = genome_cds_coordinates[cds_hash][0]
                            ###########################################################
                            # need to invert when necessary
                            if genome_cds_coordinates[cds_hash][0][1] > genome_cds_coordinates[cds_hash][0][2]:
                                print('lol')
    
                            matched_contig_len = contigs_lengths[genome_coordinates[g][0]]
                            contig_left_pos = int(genome_coordinates[g][1])
                            contig_right_pos = matched_contig_len - int(genome_coordinates[g][2])
    
                            # check LOTSC
                            if contig_left_pos < rep_left_pos and contig_right_pos < rep_right_pos:
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'LOTSC', bsr))
                                else:
                                    locus_results[og] = ['LOTSC', (rep_alleleid, 'LOTSC', bsr)]
                                continue
                            # check if PLOT
                            elif contig_left_pos < rep_left_pos:
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'PLOT3', bsr))
                                else:
                                    locus_results[og] = ['PLOT3', (rep_alleleid, 'PLOT3', bsr)]
                                continue
                            elif contig_right_pos < rep_right_pos:
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'PLOT5', bsr))
                                else:
                                    locus_results[og] = ['PLOT5', (rep_alleleid, 'PLOT5', bsr)]
                                continue
    
                            # check if ASM or ALM
                            if cds_len < (locus_mode-(locus_mode)*size_threshold):
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'ASM', bsr))
                                else:
                                    locus_results[og] = ['ASM', (rep_alleleid, 'ASM', bsr)]
                                continue
                            elif cds_len > (locus_mode+(locus_mode)*size_threshold):
                                if og in locus_results:
                                    locus_results[og][0] = 'NIPH'
                                    locus_results[og].append((rep_alleleid, 'ALM', bsr))
                                else:
                                    locus_results[og] = ['ALM', (rep_alleleid, 'ALM', bsr)]
                                continue
    
                            # add INF
                            if og in locus_results:
                                locus_results[og][0] = 'NIPH'
                                locus_results[og].append((rep_alleleid, 'BLASTp', bsr))
                            else:
                                locus_results[og] = ['INF', (rep_alleleid, 'BLASTp', bsr)]
                                # add hash of newly inferred allele
                                seen_dna.append(cds_hash)
                                seen_prot.append(prot_hash)
                                # add as representative candidate based on BSR value
                                if bsr >= 0.6 and bsr <= 0.7:
                                    representative_candidates.setdefault(locus, []).append(g)
    
            # save updated results
            fo.pickle_dumper(locus_results, locus_results_file)
    
            processed += 1
            print('\r', 'Processed: {}'.format(processed), end='')

        # select representatives

    # after the first round of allele calling, it is necessary to identify new representatives based on hits with 0.6>=BSR<=0.7
    # we need to select representatives, cluster remaining sequences against new representatives and BLAST clusters
    # we need to do this until the process does not find any new representatives
    # it should be possible to keep updating the files with the classifications forthe loci

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
    # get total number of classified CDSs! NIPHEM and NIPH encompass multiple CDSs
    total_cds = 0
    for file in classification_files:
        locus_results = fo.pickle_loader(file)
        total_cds += sum([len(c)-1 for g, c in locus_results.items()])
        all_classifications = [c[0] for g, c in locus_results.items()]
        exc += sum([1 for c in all_classifications if 'EXC' in c])
        niphem += all_classifications.count('NIPHEM')
        niph += all_classifications.count('NIPH')
        inf += all_classifications.count('INF')
        asm += all_classifications.count('ASM')
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
          'plot5: {}'.format(exc, niphem, niph, inf, asm, alm, lotsc, plot3, plot5))

    print(inf_matches)

    # create results_statistics.tsv file
    


def main(input_files, schema_directory, output_directory, ptf_path,
         blast_score_ratio, minimum_length, translation_table,
         size_threshold, word_size, window_size, clustering_sim,
         representative_filter, intra_filter, cpu_cores, blast_path,
         cds_input, prodigal_mode):#, no_cleanup):

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
                             prodigal_mode, cds_input)

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
