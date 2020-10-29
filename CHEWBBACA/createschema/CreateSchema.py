#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the creation of a whole genome multi locus sequence
typing schema seed.

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

Code documentation
------------------
"""


import os
import sys
import time
import argparse
import itertools
from itertools import repeat

from Bio import SeqIO

try:
    from utils import runProdigal
    from utils import auxiliary_functions as aux
except:
    from CHEWBBACA.utils import runProdigal
    from CHEWBBACA.utils import auxiliary_functions as aux


def gene_prediction_component(fasta_files, ptf_path, translation_table,
                              prodigal_mode, cpu_cores, temp_directory,
                              output_directory):
    """ Runs Prodigal to predict coding sequences in FASTA
        files with genomic sequence.

        Parameters
        ----------
        fasta_files : list
            List with paths to FASTA files with genomic
            sequences.
        ptf_path : str
            Path to the Prodigal training file. Should
            be False if a training file is not provided.
        translation_table : int
            Genetic code used to predict and translate
            coding sequences.
        prodigal_mode : str
            Prodigal execution mode.
        cpu_cores : int
            Number of process that will run Prodigal in
            parallel.
        temp_directory : str
            Path to the directory where temporary files
            should be stored in.
        output_directory : str
            Path to the main directory of the process.

        Returns
        -------
        A list with the following elements:
            fasta_files : list
                Input list without the paths to the files
                that Prodigal could not predict genes for.
            prodigal_path : str
                Path to the directory with the files with
                Prodigal results.

        Raises
        ------
        SystemExit
            If Prodigal could not predict genes for any of
            the input files.
    """

    # create directory to store files with Prodigal results
    prodigal_path = aux.join_paths(temp_directory, ['prodigal_cds_prediction'])
    aux.create_directory(prodigal_path)

    # run Prodigal to determine CDSs for all input genomes
    print('\nPredicting gene sequences...\n')
    # divide input genomes into equal number of sublists for
    # maximum process progress resolution
    prodigal_inputs = aux.divide_list_into_n_chunks(fasta_files,
                                                    len(fasta_files))
    # add common arguments to all sublists
    common_args = [prodigal_path, ptf_path, translation_table,
                   prodigal_mode, runProdigal.main]
    prodigal_inputs = [i+common_args for i in prodigal_inputs]

    # run Prodigal to predict genes
    prodigal_results = aux.map_async_parallelizer(prodigal_inputs,
                                                  aux.function_helper,
                                                  cpu_cores,
                                                  show_progress=True)

    # determine if Prodigal predicted genes for all genomes
    failed, failed_file = aux.check_prodigal_results(prodigal_results,
                                                     output_directory)

    if len(failed) > 0:
        print('Failed to predict genes for {0} genomes.'.format(len(failed)))
        print('Info for failed cases stored in: {0}'.format(failed_file))

    # remove failed genomes from paths
    for f in failed:
        fasta_files.remove(f[0])

    if len(fasta_files) == 0:
        sys.exit('\nCould not predict gene sequences from any '
                 'of the input files.\nPlease provide input files '
                 'in the accepted FASTA format.')

    return [fasta_files, prodigal_path]


def cds_extraction_component(fasta_files, prodigal_path, cpu_cores,
                             temp_directory, output_directory):
    """ Extracts coding sequences from FASTA files with genomic
        sequences and saves coding sequences and info about coding
        sequences to files.

        Parameters
        ----------
        fasta_files : list
            List with paths to FASTA files with genomic sequences.
        prodigal_path : str
            Path to the directory with the files with Prodigal
            results.
        cpu_cores : int
            Number of process that will extract coding sequences
            in parallel.
        temp_directory : str
            Path to the directory where temporary files should
            be stored in.
        output_directory : str
            Path to the main directory of the process.

        Returns
        -------
        cds_files : list
            List with paths to the FASTA files that contain
            extracted coding sequences.
    """

    # divide inputs into at least 20 sublists for 5% process
    # progress resolution
    num_chunks = 20 if cpu_cores < 20 else cpu_cores
    extractor_inputs = aux.divide_list_into_n_chunks(fasta_files, num_chunks)
    # add common arguments and unique index/identifier
    extractor_inputs = [[extractor_inputs[i-1], prodigal_path,
                         temp_directory, i, aux.cds_batch_extractor]
                        for i in range(1, len(extractor_inputs)+1)]

    # extract coding sequences
    print('\n\nExtracting coding sequences...\n')
    extracted_cdss = aux.map_async_parallelizer(extractor_inputs,
                                                aux.function_helper,
                                                cpu_cores,
                                                show_progress=True)

    total_extracted = sum([f[2] for f in extracted_cdss])
    print('\n\nExtracted a total of {0} coding sequences from {1} '
          'genomes.'.format(total_extracted, len(fasta_files)))

    # create full table file
    table_files = [f[0] for f in extracted_cdss]
    table_file = aux.join_paths(output_directory, ['protein_info.tsv'])
    table_header = 'Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n'
    aux.concatenate_files(table_files, table_file, table_header)
    aux.remove_files(table_files)

    cds_files = [f[1] for f in extracted_cdss]

    return cds_files


def deduplication_component(fasta_files, temp_directory, cpu_cores,
                            outfile_template):
    """ Identifies duplicated sequences in FASTA files and
        selects a distinct set of sequences.

        Parameters
        ----------
        fasta_files : list
            List with paths to FASTA files.
        temp_directory : str
            Path to the directory where new files will be
            created.
        cpu_cores : int
            Number of deduplication processes to run in
            parallel.
        outfile_template : str
            Template for the name of output files.

        Returns
        -------
        A list with the following elements:
            distinct_seqids : list
                List with the sequence identifiers of distinct
                sequences.
            distinct_seqs : str
                Path to the FASTA file with distinct sequences.
    """

    dedup_inputs = [[file,
                     aux.join_paths(temp_directory,
                                    [outfile_template.format(i+1)]),
                     aux.determine_distinct]
                    for i, file in enumerate(fasta_files)]

    # determine distinct sequences (keeps 1 seqid per sequence)
    print('\nRemoving duplicated sequences...')
    dedup_results = aux.map_async_parallelizer(dedup_inputs,
                                               aux.function_helper,
                                               cpu_cores,
                                               show_progress=False)

    # determine number of duplicated sequences
    repeated = sum([d[0] for d in dedup_results])

    # one last round after concatenating files
    if len(dedup_inputs) > 1:
        dedup_files = [f[1] for f in dedup_inputs]
        cds_file = aux.join_paths(temp_directory, ['sequences.fasta'])
        cds_file = aux.concatenate_files(dedup_files, cds_file)
        distinct_seqs = aux.join_paths(temp_directory, ['distinct_seqs.fasta'])
        dedup_results = aux.determine_distinct(cds_file, distinct_seqs)
        repeated_seqs, distinct_seqids = dedup_results

        repeated += repeated_seqs

        print('Removed {0} repeated sequences.'.format(repeated))

        return [distinct_seqids, distinct_seqs]
    else:
        print('Removed {0} repeated sequences.'.format(repeated))

        return [dedup_results[0][1], dedup_inputs[0][1]]


def small_sequences_component(fasta_file, minimum_length):
    """ Identifies sequences smaller that a specified length
        value.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file.
        minimum_length : int
            Sequences with a length value below this value are
            considered small.

        Returns
        -------
        small_seqids : list
            List with the sequence identifiers of small
            sequences.
    """

    # determine small sequences and keep their seqids
    print('\nRemoving sequences smaller than {0}...'.format(minimum_length))
    small_seqids = aux.determine_small(fasta_file, minimum_length)

    ss_lines = ['{0}: smaller than {1} chars'.format(seqid, minimum_length)
                for seqid in small_seqids]

    print('Removed {0} sequences shorter than '
          '{1}.'.format(len(small_seqids), minimum_length))

    return [small_seqids, ss_lines]


def translation_component(sequence_ids, sequences_file, temp_directory,
                          translation_table, minimum_length, cpu_cores):
    """ Translates DNA sequences, returns information about
        sequences that are untranslatable and saves
        translatable DNA sequences and proteins to FASTA
        files.

        Parameters
        ----------
        sequence_ids : list
            List with the identifiers of the sequences that
            should be translated.
        sequences_file : str
            Path to the FASTA files with the DNA sequences.
        temp_directory : str
            Path to the directory where new files will be
            created.
        translation_table : int
            Genetic code used to predict and translate
            coding sequences.
        minimum_length : int
            Sequences with a length value below this value are
            considered small and are excluded.
        cpu_cores : int
            Number of translation processes to run in
            parallel.

        Returns
        -------
        A list with the following elements:
            dna_file : str
                Path to the FASTA with all DNA sequences
                that could be translated.
            protein_file : str
                Path to the FASTA file with the translated
                sequences.
            untrans_seqids : list
                List with the identifiers of the sequences
                that could not be translated.
            untrans_lines : list
                List with one string per untranslatable
                sequence. Each string has a sequence
                identifier and a small description that
                indicates why the sequence could not be
                translated.
    """

    # translate distinct DNA sequences
    print('\nTranslating {0} DNA sequences...'.format(len(sequence_ids)))

    # divide inputs into sublists
    translation_inputs = aux.divide_list_into_n_chunks(sequence_ids, cpu_cores)

    # create paths to files with translatable DNA sequences
    dna_template = aux.join_paths(temp_directory, ['dna_{0}.fasta'])
    dna_files = [dna_template.format(i+1)
                 for i in range(len(translation_inputs))]

    # create paths to files with protein sequences
    protein_template = aux.join_paths(temp_directory, ['protein_{0}.fasta'])
    protein_files = [protein_template.format(i+1)
                     for i in range(len(translation_inputs))]

    # add common args to sublists
    common_args = [sequences_file, translation_table, minimum_length]
    translation_inputs = [[translation_inputs[i], *common_args, dna_files[i],
                           protein_files[i], aux.translate_coding_sequences]
                          for i in range(0, len(dna_files))]

    # translate sequences
    translation_results = aux.map_async_parallelizer(translation_inputs,
                                                     aux.function_helper,
                                                     cpu_cores,
                                                     show_progress=False)

    # concatenate files
    dna_file = aux.join_paths(temp_directory, ['dna.fasta'])
    dna_file = aux.concatenate_files(dna_files, dna_file)
    protein_file = aux.join_paths(temp_directory, ['protein.fasta'])
    protein_file = aux.concatenate_files(protein_files, protein_file)

    # determine sequences that could not be translated
    untrans_lines = []
    untrans_seqids = []
    for res in translation_results:
        if len(res[0]) > 0:
            untrans_lines.extend(['{0}: {1}'.format(r[0], r[1])
                                  for r in res[0]])
            untrans_seqids.extend([r[0] for r in res[0]])

    print('Removed {0} DNA sequences that could not be '
          'translated.'.format(len(untrans_seqids)))

    return [dna_file, protein_file, untrans_seqids, untrans_lines]


def clustering_helper(protein_file, cpu_cores, word_size, clustering_mode, temp_directory, intra_filter, representative_filter):
    """
    """

    # Use clustering to reduce number of BLAST comparisons
    prots = {rec.id: str(rec.seq) for rec in SeqIO.parse(protein_file, 'fasta')}

    # sort proteins by length and alphabetically
    sorted_prots = {k: v for k, v in sorted(prots.items(),
                                            key=lambda item: len(item[1]),
                                            reverse=True)}
    divided_prots = aux.split_iterable(sorted_prots, int(len(sorted_prots)/cpu_cores+10))
    divided_prots = [[d, word_size] for d in divided_prots]
    print(int(len(sorted_prots)/cpu_cores+10), len(divided_prots))
    # cluster proteins
    print('\nClustering protein sequences...')
    start_cluster = time.time()

    # need to pass word size!!!
    # clustering_results = aux.map_async_parallelizer(divided_prots, aux.cluster_sequences,
    #                                                 cpu_cores)
    clustering_results = aux.map_async_parallelizer(divided_prots, aux.cluster_helper,
                                                    cpu_cores)

    #print('Clustered {0} proteins into {1} clusters'.format(len(distinct_protein_seqs),
    #                                                        len(clusters)))
    #print(len(clustering_results))
    #for r in clustering_results:
        #print(len(r))

    # join representatives of each cluster
    rep_sequences = {}
    for r in clustering_results:
        rep_sequences = {**rep_sequences, **r[1]}

    # cluster final round
    clusters = aux.cluster_sequences(rep_sequences, word_size, 0.8,
                                     clustering_mode, minimizer=True)
    #print(len(clusters[0]))
    #p = 0
    #for k, v in clusters[0].items():
        #if len(v) > 1:
            #p += 1

    end_cluster = time.time()
    print(end_cluster-start_cluster)

    #print(p)

    # merge clusters
    all_clusters = {}
    for c in clustering_results:
        all_clusters = {**all_clusters, **c[0]}

    merged_clusters = {}
    for k, v in clusters[0].items():
        new_cluster = {k: []}
        #if len(v) > 1:
            #print(v)
        for n in v:
            new_cluster[k] += all_clusters[n[0]]

        merged_clusters[k] = new_cluster[k]

    #print(len(merged_clusters))
    #seen = []
    #total = 0
    #for j, b in merged_clusters.items():
        #total += len(b)
        #repeated = [d for d in b if d in seen]
        #if len(repeated) > 0:
            #print(j, b)
            #print('repe', repeated)
        #seen.extend(b)
    #print(total, len(set([s[0] for s in seen])))

    # write file with clustering results
    clusters_out = os.path.join(temp_directory, 'clustered_results.txt')
    aux.write_clusters(merged_clusters, clusters_out)

    # remove sequences that are very similar to representatives
    prunned_clusters, excluded_alleles = aux.representative_prunner(merged_clusters,
                                                                    representative_filter)

    # write file with prunning results
    prunned_out = os.path.join(temp_directory, 'clustered_prunned.txt')
    aux.write_clusters(prunned_clusters, prunned_out)

    excluded_alleles = [e[0] for e in excluded_alleles]

    print('Removed {0} sequences based on high similarity with '
          'cluster representative.'.format(len(excluded_alleles)))
    # determine clusters that only have the representative
    singletons = aux.determine_singletons(prunned_clusters)
    print('Found {0} singletons.'.format(len(singletons)))
    # remove singletons and keep clusters that need to be BLASTed
    final_clusters = aux.remove_entries(prunned_clusters, singletons)

    # determine number of sequences that still need to be evaluated
    # +1 to include representative
    clustered_sequences = sum([len(v)+1 for k, v in final_clusters.items()])
    print('Remaining sequences after representative and singleton prunning: {0}'.format(clustered_sequences))

    # identify clusters with more than 1 sequence besides the representative
    intra_clusters = {k: v for k, v in final_clusters.items() if len(v) > 1}

    # create kmer profile for each cluster and determine similarity
    indexed_prots = SeqIO.index(protein_file, 'fasta')

    excluded_dict = aux.intra_cluster_sim(intra_clusters, indexed_prots, word_size, intra_filter)

    # write results to file
    intrasim_out = os.path.join(temp_directory, 'clustered_intrasim.txt')
    intrasim_lines = []
    for k, v in excluded_dict.items():
        if len(v[0]) > 0:
            header_line = '>{0}'.format(k)
            intrasim_lines.append(header_line)
            sims_cases = v[1]
            sims_lines = ['\t{0}, {1}, {2}'.format(s1, s2[0], s2[1]) for s1, s2 in sims_cases.items()]
            intrasim_lines.extend(sims_lines)

    aux.write_to_file('\n'.join(intrasim_lines), intrasim_out, 'w', '\n')

    intra_excluded = [v[0] for k, v in excluded_dict.items()]
    intra_excluded = list(itertools.chain.from_iterable(intra_excluded))

    excluded_alleles = excluded_alleles + intra_excluded

    for k, v in excluded_dict.items():
        if len(v[0]) > 0:
            final_clusters[k] = [e for e in final_clusters[k] if e[0] not in v[0]]

    # add key because it is representative identifier
    clustered_sequences2 = [[k]+[e[0] for e in v] for k, v in final_clusters.items()]
    clustered_sequences2 = list(itertools.chain.from_iterable(clustered_sequences2))
    print('Remaining sequences after intra cluster prunning: {0}'.format(len(clustered_sequences2)))

    print('Clusters to BLAST: {0}'.format(len(final_clusters)))

    return [final_clusters, clustered_sequences2, indexed_prots, excluded_alleles]


def blaster_helper(temp_directory, clustered_sequences2, final_clusters, blastp_path, cpu_cores, indexed_prots):
    """
    """

    # create BLASTdb with all protein sequences from the protogenome, and with files to
    # possibilitate Blasting only against certain database sequences
    clustered_seqs_file = os.path.join(temp_directory, 'clustered_proteins.fasta')
    aux.get_sequences_by_id(indexed_prots, clustered_sequences2, clustered_seqs_file)

    integer_clusters = os.path.join(temp_directory, 'clustered_proteins_int.fasta')
    ids_dict = aux.integer_headers(clustered_seqs_file, integer_clusters)

    # create BLAST DB
    blast_db = aux.join_paths(temp_directory, ['clustered_proteins_int'])
    aux.make_blast_db(integer_clusters, blast_db, 'prot')

    blast_results_dir = os.path.join(temp_directory, 'blast_results')
    os.mkdir(blast_results_dir)

    seqids_to_blast = aux.blast_inputs(final_clusters, blast_results_dir, ids_dict)
    print(len(seqids_to_blast))
    # distribute clusters per available cores
    splitted_seqids = aux.split_blast_inputs_by_core(seqids_to_blast,
                                                     20,
                                                     blast_results_dir)

    splitted_seqids = list(map(aux.extend_list, splitted_seqids, repeat(blastp_path),
                               repeat(blast_db), repeat(blast_results_dir),
                               repeat(integer_clusters)))

    # create the FASTA files with the protein sequences before BLAST?
    print('BLASTing protein sequences in each cluster...\n')

    # BLAST each sequences in a cluster against every sequence in that cluster
    blast_results = aux.map_async_parallelizer(splitted_seqids, aux.cluster_blaster,
                                               cpu_cores, show_progress=True)

    print('\n\nFinished BLASTp. Determining schema representatives...')

    return [blast_results, ids_dict]


def main(input_files, output_directory, schema_name, ptf_path,
         blast_score_ratio, minimum_length, translation_table,
         size_threshold, clustering_mode, word_size, clustering_sim,
         representative_filter, intra_filter, cpu_cores, blastp_path,
         cds_input, prodigal_mode, verbose, cleanup):

    start_date = aux.get_datetime()
    print('Started at: {0}\n'.format(aux.datetime_str(start_date)))

    # read file with paths to input files
    fasta_files = aux.read_lines(input_files, strip=True)

    # sort paths to FASTA files
    fasta_files = aux.sort_data(fasta_files, sort_key=lambda x: x.lower())

    print('Number of genomes/assemblies: {0}'.format(len(fasta_files)))
    print('Training file: {0}'.format(ptf_path))
    print('Number of cores: {0}'.format(cpu_cores))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum sequence length: {0}'.format(minimum_length))
    print('Clustering mode: {0}'.format(clustering_mode))
    print('Word size: {0}'.format(word_size))
    print('Clustering similarity: {0}'.format(clustering_sim))
    print('Representative filter: {0}'.format(representative_filter))
    print('Intra filter: {0}'.format(intra_filter))

    # define directory for temporary files
    temp_directory = aux.join_paths(output_directory, ['temp'])

    # gene prediction step
    gp_results = gene_prediction_component(fasta_files, ptf_path,
                                           translation_table, prodigal_mode,
                                           cpu_cores, temp_directory,
                                           output_directory)

    fasta_files, prodigal_path = gp_results

    # CDS extraction step
    cds_files = cds_extraction_component(fasta_files, prodigal_path,
                                         cpu_cores, temp_directory,
                                         output_directory)

    # DNA sequences deduplication step
    distinct_dna_template = 'distinct_seqs_{0}.fasta'
    ds_results = deduplication_component(cds_files, temp_directory, cpu_cores,
                                         distinct_dna_template)

    schema_seqids, distinct_seqs_file = ds_results

    # determine small sequences step
    ss_results = small_sequences_component(distinct_seqs_file, minimum_length)

    small_seqids, ss_lines = ss_results

    # exclude seqids of small sequences
    schema_seqids = list(set(schema_seqids) - set(small_seqids))
    schema_seqids = aux.sort_data(schema_seqids, sort_key=lambda x: x.lower())

    # sequence translation step
    ts_results = translation_component(schema_seqids, distinct_seqs_file,
                                       temp_directory, translation_table,
                                       minimum_length, cpu_cores)

    dna_file, protein_file, ut_seqids, ut_lines = ts_results

    # write info about invalid alleles to fle
    invalid_alleles_file = aux.join_paths(output_directory,
                                          ['invalid_alleles.txt'])
    invalid_alleles = aux.join_list(ut_lines+ss_lines, '\n')
    aux.write_to_file(invalid_alleles, invalid_alleles_file, 'w', '\n')
    print('Info about untranslatable and small sequences '
          'stored in {0}'.format(invalid_alleles_file))

    # protein sequences deduplication step
    print('\nRemoving repeated Protein sequences...')
    distinct_prot_template = 'distinct_seqs_{0}.fasta'
    ds_results = deduplication_component([protein_file], temp_directory, 1,
                                         distinct_prot_template)

    distinct_protein_seqs, distinct_prots_file = ds_results

    distinct_protein_seqs.sort(key=lambda y: y.lower())

    # write protein FASTA file
    qc_protein_file = os.path.join(temp_directory, 'filtered_proteins.fasta')
    indexed_protein_file = SeqIO.index(protein_file, 'fasta')
    aux.get_sequences_by_id(indexed_protein_file, distinct_protein_seqs,
                            qc_protein_file)

    # write DNA FASTA file
    qc_dna_file = os.path.join(temp_directory, 'filtered_dna.fasta')
    indexed_dna_file = SeqIO.index(dna_file, 'fasta')
    aux.get_sequences_by_id(indexed_dna_file, distinct_protein_seqs,
                            qc_dna_file)

    print('\nKept {0} sequences after filtering the initial '
          'sequences.'.format(len(distinct_protein_seqs)))

    final_clusters, clustered_sequences2, indexed_prots, excluded_alleles = clustering_helper(qc_protein_file, cpu_cores,
                                                                                              word_size, clustering_mode,
                                                                                              temp_directory, intra_filter,
                                                                                              representative_filter)

    blast_results, ids_dict = blaster_helper(temp_directory, clustered_sequences2, final_clusters, blastp_path, cpu_cores, indexed_prots)

    blast_files = aux.flatten_list(blast_results)

    # index fasta file
    indexed_fasta = SeqIO.index(qc_dna_file, 'fasta')

    splitted_results = [[file, indexed_fasta, blast_score_ratio, ids_dict] for file in blast_files]

    blast_excluded_alleles = [aux.apply_bsr(i) for i in splitted_results]

    # merge bsr results
    blast_excluded_alleles = aux.flatten_list(blast_excluded_alleles)
    blast_excluded_alleles = [ids_dict[seqid] for seqid in blast_excluded_alleles]
    excluded_alleles.extend(blast_excluded_alleles)

    # perform final BLAST to avoid creating a schema with paralogs
    print('Performing a final BLAST to check for paralogs...')
    schema_seqids = list(set(distinct_protein_seqs) - set(excluded_alleles))
    beta_file = os.path.join(temp_directory, 'beta_schema.fasta')
    print('Total of {0} sequences to compare in final BLAST.'.format(len(schema_seqids)))
    aux.get_sequences_by_id(indexed_protein_file, schema_seqids, beta_file)

    integer_seqids = os.path.join(temp_directory, 'int_proteins_int.fasta')
    ids_dict2 = aux.integer_headers(beta_file, integer_seqids)

    blast_db = aux.join_paths(temp_directory, ['int_proteins_int'])
    aux.make_blast_db(integer_seqids, blast_db, 'prot')

    blast_output = '{0}/{1}_blast_out.tsv'.format(temp_directory, 'beta_schema')
    stderr = aux.run_blast(blastp_path, blast_db, integer_seqids, blast_output, 1,
                           cpu_cores)

    final_excluded = aux.apply_bsr([blast_output, indexed_fasta, blast_score_ratio, ids_dict2])
    final_excluded = [ids_dict2[seqid] for seqid in final_excluded]
    excluded_alleles = excluded_alleles + final_excluded
    schema_seqids_final = list(set(schema_seqids) - set(excluded_alleles))
    print('Removed {0} loci that were too similar with other loci in the schema.'.format(len(final_excluded)))

    output_schema = os.path.join(temp_directory, 'schema_seed.fasta')

    # create file with the schema representative sequences
    aux.get_sequences_by_id(indexed_fasta, schema_seqids_final, output_schema)

    schema_dir = aux.join_paths(output_directory, [schema_name])
    aux.create_directory(schema_dir)

    # create directory and schema files
    filenames = (record.name for record in SeqIO.parse(output_schema, 'fasta'))
    schema_files = aux.split_fasta(output_schema, schema_dir, 1, filenames)
    aux.create_short(schema_files, schema_dir)

    # remove temporary files
    if cleanup is True:
        aux.delete_directory(temp_directory)

    end_date = aux.get_datetime()
    end_date_str = aux.datetime_str(end_date)

    minutes, seconds = aux.datetime_diff(start_date, end_date)

    print('Created schema with {0} genes based on {1} genomes in'
          '{2: .0f}m{3: .0f}s.'.format(len(schema_files), len(fasta_files),
                                       minutes, seconds))
    print('\nFinished at: {0}'.format(end_date_str))


def parse_arguments():


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_directory',
                        help='Output directory where the process will store '
                             'intermediate files and create the schema\'s directory.')

    parser.add_argument('--n', type=str, required=False,
                        default='schema_seed', dest='schema_name',
                        help='Name to give to folder that will store the schema files.')

    parser.add_argument('--ptf', type=str, required=False,
                        default=False, dest='ptf_path',
                        help='Path to the Prodigal training file.')

    parser.add_argument('--bsr', type=float,
                        required=False, default=0.6, dest='blast_score_ratio',
                        help='BLAST Score Ratio value. Sequences with '
                             'alignments with a BSR value equal to or '
                             'greater than this value will be considered '
                             'as sequences from the same gene.')

    parser.add_argument('--l', type=int,
                        required=False, default=201, dest='minimum_length',
                        help='Minimum sequence length accepted for a '
                             'coding sequence to be included in the schema.')

    parser.add_argument('--t', type=int,
                        required=False, default=11, dest='translation_table',
                        help='Genetic code used to predict genes and'
                             ' to translate coding sequences.')

    parser.add_argument('--st', type=float,
                        required=False, default=0.2, dest='size_threshold',
                        help='CDS size variation threshold. At the default '
                             'value of 0.2, alleles with size variation '
                             '+-20 percent will be classified as ASM/ALM.')

    parser.add_argument('--cpu', type=int, required=False,
                        default=1, dest='cpu_cores',
                        help='Number of CPU cores that will be '
                             'used to run the CreateSchema process '
                             '(will be redefined to a lower value '
                             'if it is equal to or exceeds the total'
                             'number of available CPU cores).')

    parser.add_argument('--cm', type=str, required=False,
                        default='greedy', dest='clustering_mode',
                        help='The clustering mode. There are two modes: '
                             'greedy and full. Greedy will add sequences '
                             'to a single cluster. Full will add sequences '
                             'to all clusters they share high similarity with.')
    
    parser.add_argument('--ws', type=int, required=False,
                        default=4, dest='word_size',
                        help='Value of k used to decompose protein sequences '
                             'into k-mers.')

    parser.add_argument('--cs', type=float, required=False,
                        default=0.20, dest='clustering_sim',
                        help='Similarity threshold value necessary to '
                             'consider adding a sequence to a cluster. This '
                             'value corresponds to the percentage of shared k-mers.')

    parser.add_argument('--rf', type=float, required=False,
                        default=0.80, dest='representative_filter',
                        help='Similarity threshold value that is considered '
                             'to determine if a sequence belongs to the same '
                             'gene as the cluster representative purely based '
                             'on the percentage of shared k-mers.')
    
    parser.add_argument('--if', type=float, required=False,
                        default=0.80, dest='intra_filter',
                        help='Similarity threshold value that is considered '
                             'to determine if sequences in the same custer '
                             'belong to the same gene. Only one of those '
                             'sequences is kept.')

    parser.add_argument('--b', type=str, required=False,
                        default='blastp', dest='blastp_path',
                        help='Path to the BLASTp executables.')

    parser.add_argument('--CDS', required=False, action='store_true',
                        dest='cds_input',
                        help='Input is a FASTA file with one representative '
                             'sequence per gene in the schema.')

    parser.add_argument('--pm', required=False, choices=['single', 'meta'],
                        default='single', dest='prodigal_mode',
                        help='Prodigal running mode.')

    parser.add_argument('--v', required=False, action='store_true',
                        dest='verbose',
                        help='Increased output verbosity during execution.')

    parser.add_argument('--c', '--cleanup', required=False, action='store_true',
                        dest='cleanup',
                        help='Delete intermediate files at the end.')

    args = parser.parse_args()

    return [args.input_files, args.output_directory, args.schema_name,
            args.ptf_path, args.blast_score_ratio, args.minimum_length,
            args.translation_table, args.size_threshold,
            args.clustering_mode, args.word_size, args.clustering_sim,
            args.representative_filter, args.intra_filter, args.cpu_cores,
            args.blastp_path, args.cds_input, args.prodigal_mode, args.verbose,
            args.cleanup]


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2], args[3], args[4], args[5],
         args[6], args[7], args[8], args[9], args[10], args[11],
         args[12], args[13], args[14], args[15], arg[16], args[17],
         args[18])
