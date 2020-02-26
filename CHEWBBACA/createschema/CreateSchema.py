#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import time
import pickle
import shutil
import argparse
import itertools
from collections import Counter
from multiprocessing import Pool

from Bio import SeqIO

import new_utils
import runProdigal
import CreateSchema_aux as rfm


#input_files = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/ref32_genomes'
#output_directory = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/dup_ref32_schema_seed'
#prodigal_training_file = '/home/rfm/Desktop/rfm/Lab_Software/chewBBACA/CHEWBBACA/prodigal_training_files/Streptococcus_agalactiae.trn'
#schema_name = 'dup_ref32_schema_seed'
#cpu_count = 6
#blastp_path = shutil.which('blastp')
#blast_score_ratio = 0.6
#minimum_cds_length = 201
#cleanup = 'yes'
## cluster parameters
#clustering_mode = 'greedy'
#filtering_sim = 0.15
#word_size = 4
#clustering_sim = 0.20
#cluster_filter = 0.80


def main(input_files, output_directory, prodigal_training_file, schema_name, cpu_count,
         blastp_path, blast_score_ratio, minimum_cds_length, translation_table, clustering_mode,
         word_size, clustering_sim, rep_filter, intra_filter, cleanup):

    global_start = time.time()

    cpu_to_apply = new_utils.verify_cpu_usage(cpu_count)

    chosen_taxon = os.path.abspath(prodigal_training_file)

    print('\nCreating schema based on the genomes '
          'in the following directory:\n{0}'.format(os.path.abspath(input_files)))
    print('Training file: {0}'.format(chosen_taxon))
    print('Number of cores: {0}'.format(cpu_count))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum accepted sequence length: {0}'.format(minimum_cds_length))
    print('Clustering mode: {0}'.format(clustering_mode))
    print('Word size: {0}'.format(word_size))
    print('Clustering similarity: {0}'.format(clustering_sim))
    print('Representative filter: {0}'.format(rep_filter))
    print('Intra filter: {0}'.format(intra_filter))

    # path to directory with genomes FASTA files
    genomes_path = input_files

    # list files in genomes FASTA files directory and determine absolute path for each
    fasta_files = os.listdir(genomes_path)
    for f in range(len(fasta_files)):
        fasta_files[f] = os.path.join(genomes_path, fasta_files[f])

    # maintain genome order to assign identifiers correctly
    fasta_files.sort(key=lambda y: y.lower())
    print('Number of genomes/assemblies: {0}'.format(len(fasta_files)))

    # determine and store genome identifiers
    genomes_identifiers = []
    for genome in fasta_files:
        genomes_identifiers.append(rfm.genome_id(genome))

    # define parent directory
    parent_directory = os.path.abspath(output_directory)
    if not os.path.exists(parent_directory):
        os.makedirs(parent_directory)

    # define directory for temporary files
    temp_directory = os.path.join(parent_directory, 'temp')

    # define output directory where Prodigal files will be stored
    prodigal_path = os.path.join(temp_directory, 'prodigal_cds_prediction')

    # create directory for Prodigal files, if it does not exist
    if not os.path.exists(prodigal_path):
        os.makedirs(prodigal_path)

    # run Prodigal to determine CDSs for all input genomes
    prodigal_start = time.time()
    print('\nPredicting CDS sequences...')
    print ('Started Prodigal at: {0}'.format(time.strftime('%H:%M:%S - %d/%m/%Y')))

    # run Prodigal with multiprocessing
    pool = Pool(cpu_to_apply)
    for genome in fasta_files:
        pool.apply_async(runProdigal.main, (genome, prodigal_path, chosen_taxon, translation_table))
    pool.close()
    pool.join()

    print('\nChecking if Prodigal created all the necessary files...')
    new_utils.check_prodigal_output_files(prodigal_path, fasta_files)
    print('Finished Prodigal at: {0}'.format(time.strftime("%H:%M:%S - %d/%m/%Y")))
    prodigal_end = time.time()
    prodigal_delta = prodigal_end - prodigal_start
    print('Prodigal delta: {0}'.format(prodigal_delta))

    # get CDSs for each genome
    protein_table = os.path.join(parent_directory, 'protein_info.tsv')
    with open(protein_table, 'w') as file:
        file.write('Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n')

    # add multiprocessing!!!
    print('Extracting coding sequences from all genomes...')
    cds_start = time.time()
    protid = 1
    cds_file = os.path.join(temp_directory, 'coding_sequences.fasta')
    # divide list of input genomes per core and write into separate files to be capable of
    # multiprocessing?
    for g in range(len(fasta_files)):
        # determine Prodigal ORF file path for current genome
        orf_file_path = os.path.join(prodigal_path, '{0}_ORF.txt'.format(genomes_identifiers[g]))
        # import contigs for current genome/assembly (importing for all genomes into memory might occupy too much memory)
        contigs = rfm.import_sequences(fasta_files[g])
        # extract coding sequences from contigs
        genome_info = rfm.extract_coding_sequences(orf_file_path, contigs, protid)
        # save coding sequences to file
        genome_id = genomes_identifiers[g]

        # create records and write them to file
        cds_lines = rfm.create_fasta_lines(genome_info[0], genome_id)
        rfm.write_fasta(cds_lines, cds_file)

        rfm.write_protein_table(protein_table, genome_id, genome_info[1])

        # keep track of CDSs identifiers to assign them sequentially
        protid = 1

    cds_end = time.time()
    cds_delta = cds_end - cds_start
    print('CDS extraction delta: {0}:'.format(cds_delta))
    # determine repeated sequences and keep only one representative
    r1_start = time.time()
    repeated_seqs_file = os.path.join(temp_directory, 'repeated_dna_seqids.txt')
    unique_dna_seqids = os.path.join(temp_directory, 'unique_dna_seqids.fasta')
    repeated_dna_cds = rfm.determine_repeated(cds_file, repeated_seqs_file, unique_dna_seqids)
    r1_end = time.time()
    r1_delta = r1_end - r1_start

    print('Repeated DNA sequences removal delta: {0}:'.format(r1_delta))
    print('\nRemoved {0} repeated DNA sequences.'.format(repeated_dna_cds[0]))

    # determine small DNA sequences and remove those seqids
    small_dna_cds = rfm.determine_small(unique_dna_seqids, minimum_cds_length)

    valid_dna_seqs = small_dna_cds[1]
    print('Removed {0} DNA sequences shorter than {1} nucleotides.'.format(len(small_dna_cds[0]), minimum_cds_length))
    valid_dna_seqs.sort(key=lambda y: y.lower())

    # convert to protein
    trans_start = time.time()
    dna_valid_file = os.path.join(temp_directory, 'valid_dna.fasta')
    protein_valid_file = os.path.join(temp_directory, 'valid_protein.fasta')
    print('\nTranslating {0} DNA sequences...'.format(len(valid_dna_seqs)))
    untranslatable_cds, total_seqs = rfm.translate_coding_sequences(unique_dna_seqids, valid_dna_seqs,
                                                                    dna_valid_file, protein_valid_file, translation_table)
    untranslatable_seqids = [t[0] for t in untranslatable_cds]

    trans_end = time.time()
    trans_delta = trans_end - trans_start
    print('Translation delta: {0}:'.format(trans_delta))
    # write file with invalid alleles info
    invalid_alleles_file = os.path.join(parent_directory, 'invalid_alleles.txt')
    with open(invalid_alleles_file, 'w') as inv:
        lines = []
        for allele in untranslatable_cds:
            seqid = allele[0]
            error = allele[1]
            line = '{0}: {1}\n'.format(seqid, error)
            lines.append(line)

        inv.writelines(lines)
        
    # add multiprocessing!!!
    # remove DNA sequences that could not be translated
    translatable_cds = list(set(valid_dna_seqs) - set(untranslatable_seqids))
    print('Removed {0} DNA sequences that could not be translated.'.format(len(untranslatable_seqids)))
    print('Info about untranslatable alleles stored in {0}'.format(invalid_alleles_file))
    translatable_cds.sort(key=lambda y: y.lower())

    # next round of finding repeated sequences, but for proteins
    repeated_prots_file = os.path.join(temp_directory, 'repeated_prots_seqids.txt')
    unique_prots_file = os.path.join(temp_directory, 'unique_prots_seqs.fasta')
    repeated_protein_cds = rfm.determine_repeated(protein_valid_file, repeated_prots_file, unique_prots_file)

    # remove seqids that are from repeated protein sequences
    unique_protein_seqs = repeated_protein_cds[1]
    #print('Removed {0} repeated Protein sequences.'.format(len(repeated_protein_cds[0])))
    unique_protein_seqs.sort(key=lambda y: y.lower())

    final_seqids = unique_protein_seqs
    final_seqids.sort(key=lambda y: y.lower())

    # write protein FASTA file
    protein_file = os.path.join(temp_directory, 'filtered_proteins.fasta')
    indexed_protein_valid_file = SeqIO.index(protein_valid_file, 'fasta')
    rfm.get_sequences_by_id(indexed_protein_valid_file, final_seqids, protein_file)

    # write DNA FASTA file
    dna_file = os.path.join(temp_directory, 'filtered_dna.fasta')
    indexed_dna_valid_file = SeqIO.index(dna_valid_file, 'fasta')
    rfm.get_sequences_by_id(indexed_dna_valid_file, final_seqids, dna_file)

    print('Kept {0} sequences after filtering the initial sequences.'.format(len(final_seqids)))

    # Use clustering to reduce number of BLAST comparisons
    # should not import all proteins into memory you noob!!!
    # proteins should be imported one by one and clustered that way?
    # maybe sort them in new file and then use SeqIO.parse and cluster one by one
    prots = {}
    for record in SeqIO.parse(protein_file, 'fasta'):
        seqid = record.id
        prot_seq = str(record.seq)
        prots[seqid] = prot_seq

    # sort proteins by length and alphabetically
    sorted_prots = sorted(list(prots.items()),
                          key=lambda x: (-len(x[1]), x[0]))

    # cluster proteins
    cluster_start = time.time()
    clusters = rfm.cluster_sequences(sorted_prots, word_size,
                                     clustering_sim, clustering_mode)
    print('Clustered proteins into {0} clusters'.format(len(clusters)))
    cluster_end = time.time()
    cluster_delta = cluster_end - cluster_start
    print('Cluster delta: {0}'.format(cluster_delta))

    # remove sequences that are very similar to representatives
    prunned_clusters, excluded_alleles = rfm.cluster_prunner(clusters,
                                                             rep_filter)
    excluded_alleles = [e[0] for e in excluded_alleles]
    # determine clusters that only have the representative
    singletons = rfm.determine_singletons(prunned_clusters)
    print('Singletons: {0}'.format(len(singletons)))
    # remove singletons and keep clusters that need to be BLASTed
    final_clusters = rfm.remove_clusters(prunned_clusters, singletons)
    print('Clusters to BLAST: {0}'.format(len(final_clusters)))
    
    # determine number of sequences that still need to be evaluated
    # +1 to include representative
    clustered_sequences = sum([len(v)+1 for k, v in final_clusters.items()])
    print('Remaining sequences after clustering and prunning: {0}'.format(clustered_sequences))

    # identify clusters with more than 1 sequence besides the representative
    print('Determining intra cluster similarity...')

    intra_start = time.time()
    intra_clusters = {k: v for k, v in final_clusters.items() if len(v) > 1}

    # create kmer profile for each cluster and determine similarity
    indexed_prots = SeqIO.index(protein_file, 'fasta')

    excluded_dict = rfm.intra_cluster_sim(intra_clusters, indexed_prots, word_size, intra_filter)

    intra_excluded = [v[0] for k, v in excluded_dict.items()]
    intra_excluded = list(itertools.chain.from_iterable(intra_excluded))

    intra_end = time.time()
    intra_delta = intra_end - intra_start
    print('Intra cluster similarity: {0}'.format(intra_delta))

    excluded_alleles = excluded_alleles + intra_excluded

    for k, v in excluded_dict.items():
        if len(v[0]) > 0:
            final_clusters[k] = [e for e in final_clusters[k] if e[0] not in v[0]]

    # add key because it is representative identifier
    clustered_sequences2 = [[k]+[e[0] for e in v] for k, v in final_clusters.items()]
    clustered_sequences2 = list(itertools.chain.from_iterable(clustered_sequences2))
    print('Remaining sequences after intra cluster prunning: {0}'.format(len(clustered_sequences2)))

    # create BLASTdb with all protein sequences from the protogenome, and with files to
    # possibilitate Blasting only against certain database sequences
    clustered_seqs_file = os.path.join(temp_directory, 'clustered_proteins.fasta')
    rfm.get_sequences_by_id(indexed_prots, clustered_sequences2, clustered_seqs_file)
    
    makedb_cmd = 'makeblastdb -in {0} -parse_seqids -dbtype prot >/dev/null 2>&1'.format(clustered_seqs_file)
    os.system(makedb_cmd)

    # BLAST necessary sequences against only the sequences from the same cluster
    blast_db = os.path.join(temp_directory, 'clustered_proteins.fasta')

    blast_results_dir = os.path.join(temp_directory, 'blast_results')
    os.mkdir(blast_results_dir)

    seqids_to_blast = rfm.blast_inputs(final_clusters, blast_results_dir)

    # distribute clusters per available cores, try to group inputs into
    # even groups in terms of number of clusters and sum of number of sequences
    # per inputs group
    splitted_seqids = rfm.split_blast_inputs_by_core(seqids_to_blast,
                                                     cpu_to_apply,
                                                     blast_results_dir)

    for s in range(len(splitted_seqids)):
        splitted_seqids[s].append(blastp_path)
        splitted_seqids[s].append(blast_db)
        splitted_seqids[s].append(blast_results_dir)
        splitted_seqids[s].append(clustered_seqs_file)
    
    # create the FASTA files with the protein sequences before BLAST?
    
    
    blast_start = time.time()
    print('BLASTing protein sequences in each cluster...')
    # ADD progress bar!!!
    # BLAST each sequences in a cluster against every sequence in that cluster
    p = Pool(processes = cpu_to_apply)
    r = p.map_async(rfm.cluster_blaster, splitted_seqids)
    r.wait()

    print('Finished BLASTp. Determining schema representatives...')
    blast_end = time.time()
    blast_delta = blast_end - blast_start
    print('BLAST delta: {0}'.format(blast_delta))

    # Import BLAST results
    blast_files = os.listdir(blast_results_dir)
    blast_files = [os.path.join(blast_results_dir, file) for file in blast_files if 'blast_out.tsv' in file]

    # index fasta file
    indexed_fasta = SeqIO.index(dna_file, 'fasta')

    splitted_results = []
    for file in blast_files:
        splitted_results.append([file, indexed_fasta, blast_score_ratio])

    blast_excluded_alleles = []
    for i in splitted_results:
        try:
            e = rfm.apply_bsr(i)
            blast_excluded_alleles.append(e)
        # proteins with low complexity regions might
        # fail to align 
        except Exception:
            print('Could not apply BSR to {0} results. Probably low complexity proteins. FAGGOTS!'.format(i))

    # merge bsr results
    blast_excluded_alleles = list(itertools.chain.from_iterable(blast_excluded_alleles))
    excluded_alleles.extend(blast_excluded_alleles)

    # perform final BLAST to avoid creating a schema with paralogs
    print('Performing a final BLAST to check for paralogs...')
    schema_seqids = list(set(final_seqids) - set(excluded_alleles))
    beta_file = os.path.join(temp_directory, 'beta_schema.fasta')
    rfm.get_sequences_by_id(indexed_protein_valid_file, schema_seqids, beta_file)

    makedb_cmd = 'makeblastdb -in {0} -dbtype prot >/dev/null 2>&1'.format(beta_file)
    os.system(makedb_cmd)

    blast_db = os.path.join(temp_directory, beta_file)
    
    second_blast_start = time.time()
    blast_output = '{0}/{1}_blast_out.tsv'.format(temp_directory, 'beta_schema')
    blast_command = ('{0} -db {1} -query {2} -out {3} -outfmt "6 qseqid sseqid score" '
                     '-max_hsps 1 -num_threads {4} -evalue 0.001'.format(blastp_path, beta_file, beta_file, blast_output, cpu_to_apply))
    os.system(blast_command)
    second_blast_end = time.time()
    second_blast_delta = second_blast_end - second_blast_start
    print('Total BLAST time: {0}'.format(second_blast_delta))

    final_excluded = rfm.apply_bsr([blast_output, indexed_fasta, blast_score_ratio])
    excluded_alleles = excluded_alleles + final_excluded
    schema_seqids_final = list(set(schema_seqids) - set(excluded_alleles))
    print('Removed {0} loci that were too similar with other loci in the schema.'.format(len(final_excluded)))

    # decide which identifiers to keep and redo BLAST to check for residual paralogs
    output_schema = os.path.join(temp_directory, 'cd_hit_chewie.fasta')

    # create file with the schema representative sequences
    rfm.get_schema_reps(indexed_fasta, schema_seqids_final,
                        protein_file, output_schema)

    schema_dir = os.path.join(parent_directory, schema_name)
    if not os.path.exists(schema_dir):
        os.makedirs(schema_dir)

    # create directory and schema files
    rfm.build_schema(output_schema, schema_dir)

    # write hidden config files and genes list
    # write hidden config file with parameters
    config_file = os.path.join(schema_dir, '.schema_config')
    shutil.copy(prodigal_training_file, schema_dir)
    ptf_name = os.path.basename(prodigal_training_file)

    # create dictionary with parameters values
    parameters = {'bsr': blast_score_ratio,
                  'prodigal_training_file': ptf_name,
                  'translation_table': 11,
                  'minimum_locus_length': minimum_cds_length,
                  'chewBBACA_version': '2.1.0'}

    with open(config_file, 'wb') as cf:
        pickle.dump(parameters, cf)

    # create hidden file with genes/loci list
    schema_files = [file for file in os.listdir(schema_dir) if '.fasta' in file]
    schema_list_file = os.path.join(schema_dir, '.genes_list')
    with open(schema_list_file, 'wb') as sl:
        pickle.dump(schema_files, sl)
    
    # remove temporary files
    if cleanup == 'yes':
        shutil.rmtree(temp_directory)

    global_end = time.time()
    global_delta = global_end - global_start
    print('Created schema based on {0} assemblies/genomes of {1} in {2}m{3}s.'.format(len(fasta_files),
                                                                  os.path.basename(chosen_taxon).rstrip('.trn'),
                                                                  int(global_delta/60), int(((global_delta/60)%1)*60)))


def parse_arguments():


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input_files', type=str, required=True, dest='input_files',
                        help='List of FASTA files, with complete genomes or assemblies'
                             ' given as input to create a new schema.')

    parser.add_argument('-o', '--output_directory', type=str, required=True, dest='output_directory',
                        help='Name or path to the output directory.')

    parser.add_argument('-ptf', '-prodigal_training_file', type=str, required=True, dest='prodigal_training_file',
                        help='Full path to user-provided Prodigal training file.')

    parser.add_argument('-n', '--schema_name', type=str, required=False, dest='schema_name',
                        help='Name to give to folder that will store the final schema.')

    parser.add_argument('--cpu', '--cpu_count', type=int, required=False, dest='cpu_count',
                        default=1, help='Number of CPU cores to use (default=1). If the provided value '
                                        'exceeds the maximum number of available cores, will'
                                        'maximum -2.')

    parser.add_argument('--bp', '--blastp_path', type=str, required=False, dest='blastp_path',
                        default='blastp', help='Path to BLASTp executable (default=blastp, assumes '
                                          'that BLAST was added to the Path).')

    parser.add_argument('--bsr', '--blast_score_ratio', type=float, required=False, dest='blast_score_ratio',
                        default=0.6, help='The BLAST Score Ratio value that will be used to compare CDSs and '
                                          'create the new schema (default=0.6).')

    parser.add_argument('--l', '--minimum_cds_length', type=int, required=False, dest='minimum_cds_length',
                        default=201, help='Minimum acceptable CDS length (default=201).')

    parser.add_argument('--t', '--translation_table', type=int, required=False, dest='translation_table',
                        default=11, help='Genetic code to use for gene prediction with Prodigal and CDS '
                        'translation (default=11, Bacteria and Archaea).')

    parser.add_argument('--cm', '--clustering_mode', type=str, required=False, dest='clustering_mode',
                        default='greedy', help='')

    parser.add_argument('--wf', '--word_filter', type=int, required=False, dest='word_filter',
                        default=4, help='')
    
    parser.add_argument('--fs', '--filtering_sim', type=float, required=False, dest='filtering_sim',
                        default=0.15, help='')

    parser.add_argument('--cs', '--clustering_sim', type=float, required=False, dest='clustering_sim',
                        default=0.20, help='')

    parser.add_argument('--ws', '--word_size', type=int, required=False, dest='word_size',
                        default=4, help='')

    parser.add_argument('--cf', '--cluster_filter', type=float, required=False, dest='cluster_filter',
                        default=0.80, help='')

    parser.add_argument('--c', '--cleanup', type=str, required=False, dest='cleanup',
                        default='yes', help='If the temporary directory should be deleted at the end '
                        '(default=yes).')

    args = parser.parse_args()

    return [args.input_files, args.output_directory, args.prodigal_training_file,
            args.schema_name, args.cpu_count, args.blastp_path, args.blast_score_ratio,
            args.minimum_cds_length, args.translation_table, args.clustering_mode, args.word_filter,
            args.filtering_sim, args.clustering_sim, args.word_size,
            args.cluster_filter, args.cleanup]


if __name__ == '__main__':

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4],
         args[5], args[6], args[7], args[8], args[9],
         args[10], args[11], args[12], args[13], args[14],
         args[15])
