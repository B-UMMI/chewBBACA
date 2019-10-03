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
import shutil
import argparse
import itertools
from multiprocessing import Pool

from Bio import SeqIO

import new_utils
import runProdigal
import CreateSchema_aux as rfm


input_files = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/GBS_680_genomes'
output_directory = '/home/rfm/Desktop/rfm/Lab_Software/CreateSchema_tests/new_create_schema_scripts/sagalactiae_680ref_schema'
prodigal_training_file = '/home/rfm/Desktop/rfm/Lab_Software/chewBBACA/CHEWBBACA/prodigal_training_files/Streptococcus_agalactiae.trn'
schema_name = 'sagalactiae_schema_seed'
cpu_count = 6
blastp_path = shutil.which('blastp')
blast_score_ratio = 0.6
minimum_cds_length = 201
cleanup = 'yes'
# cluster parameters
clustering_mode = 'greedy'
word_filter = 10
word_size = 4
clustering_sim = 0.0
cluster_filter_sim = 0.80


def main(input_files, output_directory, prodigal_training_file, schema_name, cpu_count, 
         blastp_path, blast_score_ratio, minimum_cds_length, cd_hit_sim, cd_hit_word,
         cutoff_sim, cleanup):
    """
    """

    start = time.time()

    cpu_to_apply = new_utils.verify_cpu_usage(cpu_count)
    print('\nNumber of cores: {0}'.format(cpu_to_apply))

    chosen_taxon = os.path.abspath(prodigal_training_file)
    print('Training file: {0}'.format(chosen_taxon))

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
    print('\nPredicting CDS sequences...')
    print ('Started Prodigal at: {0}'.format(time.strftime('%H:%M:%S - %d/%m/%Y')))

    # run Prodigal with multiprocessing
    pool = Pool(cpu_to_apply)
    for genome in fasta_files:
        pool.apply_async(runProdigal.main, (genome, prodigal_path, chosen_taxon))
    pool.close()
    pool.join()

    print('\nChecking if Prodigal created all the necessary files...')
    new_utils.check_prodigal_output_files(prodigal_path, fasta_files)
    print('Finished Prodigal at: {0}'.format(time.strftime("%H:%M:%S - %d/%m/%Y")))

    # get CDSs for each genome
    protein_table = os.path.join(parent_directory, 'protein_info.tsv')
    with open(protein_table, 'w') as file:
        file.write('Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n')

    # add multiprocessing!!!
    protid = 1
    cds_ids = []
    cds_file = os.path.join(temp_directory, 'coding_sequences.fasta')
    for g in range(len(fasta_files)):
        # determine Prodigal ORF file path for current genome
        orf_file_path = os.path.join(prodigal_path, '{0}_ORF.txt'.format(genomes_identifiers[g]))
        # import contigs for current genome/assembly (importing for all genomes into memory might occupy too much memory)
        contigs = rfm.import_contigs(fasta_files[g])
        # extract coding sequences from contigs
        genome_info = rfm.extract_coding_sequences(orf_file_path, contigs, protid)
        # save coding sequences to file
        genome_id = genomes_identifiers[g]
        # create sequence ids
        genome_ids = [genome_id+'-protein'+str(seqid) for seqid in genome_info[0].keys()]

        cds_lines = rfm.create_fasta_lines(genome_info[0], genome_id)
        rfm.write_fasta(cds_lines, cds_file)
        cds_ids.extend(genome_ids)
        rfm.write_protein_table(protein_table, genome_id, genome_info[1])

        # keep track of CDSs identifiers to assign them sequentially
        protid = 1

    # determine seqids of repeated sequences
    repeated_dna_cds = rfm.determine_repeated(cds_file)

    # remove seqids that are from repeated sequences
    unique_dna_seqs = list(set(cds_ids) - set(repeated_dna_cds))
    print('\nRemoved {0} repeated DNA sequences.'.format(len((repeated_dna_cds))))
    unique_dna_seqs.sort(key=lambda y: y.lower())

    # determine small DNA sequences and remove those seqids
    small_dna_cds = rfm.determine_small(cds_file, minimum_cds_length)

    valid_dna_seqs = list(set(unique_dna_seqs) - set(small_dna_cds))
    print('Removed {0} DNA sequences shorter than {1} nucleotides.'.format(len(small_dna_cds), minimum_cds_length))
    valid_dna_seqs.sort(key=lambda y: y.lower())

    # convert to protein
    dna_valid_file = os.path.join(temp_directory, 'valid_dna.fasta')
    protein_valid_file = os.path.join(temp_directory, 'valid_protein.fasta')
    print('\nTranslating {0} DNA sequences...'.format(len(valid_dna_seqs)))
    untranslatable_cds, total_seqs = rfm.translate_coding_sequences(cds_file, valid_dna_seqs,
                                                                    dna_valid_file, protein_valid_file, 11)
    untranslatable_seqids = [t[0] for t in untranslatable_cds]

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
    repeated_protein_cds = rfm.determine_repeated(protein_valid_file)

    # remove seqids that are from repeated protein sequences
    unique_protein_seqs = list(set(translatable_cds) - set(repeated_protein_cds))
    print('Removed {0} repeated Protein sequences.'.format(len(repeated_protein_cds)))
    unique_protein_seqs.sort(key=lambda y: y.lower())

    final_seqids = unique_protein_seqs
    final_seqids.sort(key=lambda y: y.lower())

    # write protein FASTA file
    protein_file = os.path.join(temp_directory, 'filtered_proteins.fasta')
    rfm.get_sequences_by_id(protein_valid_file, final_seqids, protein_file)

    # write DNA FASTA file
    dna_file = os.path.join(temp_directory, 'filtered_dna.fasta')
    rfm.get_sequences_by_id(dna_valid_file, final_seqids, dna_file)

    print('Kept {0} sequences after filtering the initial {1} sequences.'.format(len(final_seqids), len(cds_ids)))

    # Use clustering to reduce number of BLAST comparisons
    # import proteins to cluster
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
    clusters = rfm.cluster_sequences(sorted_prots, word_filter,
                                     word_size, clustering_sim,
                                     clustering_mode)
    cluster_end = time.time()
    cluster_delta = cluster_end - cluster_start

    # remove sequences that are very similar to representatives
    prunned_clusters, excluded_alleles = rfm.cluster_prunner(clusters,
                                                             cluster_filter_sim)
    excluded_alleles = [e[0] for e in excluded_alleles]
    # determine clusters that only have the representative
    singletons = rfm.determine_singletons(prunned_clusters)
    # remove singletons and keep clusters that need to be BLASTed
    final_clusters = rfm.remove_clusters(prunned_clusters, singletons)

    # create dictionary mapping seqids to sequences
    protein_dict_file = os.path.join(temp_directory, 'protein_dictionary')
    rfm.create_protein_dict(protein_file, protein_dict_file)

    # create BLASTdb with all protein sequences from the protogenome, and with files to
    # possibilitate Blasting only against certain database sequences
    makedb_cmd = 'makeblastdb -in {0} -parse_seqids -dbtype prot >/dev/null 2>&1'.format(protein_file)
    os.system(makedb_cmd)

    # BLAST necessary sequences against only the sequences from the same cluster
    blast_db = os.path.join(temp_directory, 'filtered_proteins.fasta')

    blast_results_dir = os.path.join(temp_directory, 'blast_results')
    os.mkdir(blast_results_dir)

    seqids_to_blast = rfm.blast_inputs(final_clusters, blast_db, 
                                       protein_dict_file, blast_results_dir)

    # distribute clusters per available cores, try to group inputs into
    # even groups in terms of number of clusters and sum of number of sequences
    # per inputs group
    splitted_seqids = rfm.split_blast_inputs_by_core(seqids_to_blast,
                                                     cpu_to_apply, blast_results_dir)
    for s in range(len(splitted_seqids)):
        splitted_seqids[s].append(blastp_path)
        splitted_seqids[s].append(blast_results_dir)

    print('BLASTing protein sequences in each cluster determined by CD-HIT...')
    # ADD progress bar!!!
    # BLAST each sequences in a cluster against every sequence in that cluster
    p = Pool(processes = cpu_to_apply)
    r = p.map_async(rfm.cluster_blaster, splitted_seqids)
    r.wait()

    print('Finished BLASTp. Determining schema representatives...')

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
        e = rfm.apply_bsr(i)
        blast_excluded_alleles.append(e)

    # merge bsr results
    blast_excluded_alleles = list(itertools.chain.from_iterable(blast_excluded_alleles))
    excluded_alleles.extend(blast_excluded_alleles)

    # perform final BLAST to avoid creating a schema with paralogs
    print('Performing a final BLAST to check for paralogs...')
    schema_seqids = list(set(final_seqids) - set(excluded_alleles))
    beta_file = os.path.join(temp_directory, 'beta_schema.fasta')
    rfm.get_sequences_by_id(protein_valid_file, schema_seqids, beta_file)

    makedb_cmd = 'makeblastdb -in {0} -dbtype prot >/dev/null 2>&1'.format(beta_file)
    os.system(makedb_cmd)

    blast_db = os.path.join(temp_directory, beta_file)

    blast_output = '{0}/{1}_blast_out.tsv'.format(temp_directory, 'beta_schema')
    blast_command = ('{0} -db {1} -query {2} -out {3} -outfmt "6 qseqid sseqid score" '
                     '-max_hsps 1 -num_threads {4} -evalue 0.001'.format(blastp_path, beta_file, beta_file, blast_output, cpu_to_apply))
    os.system(blast_command)

    final_excluded = rfm.apply_bsr([blast_output, indexed_fasta, blast_score_ratio])
    excluded_alleles = excluded_alleles + final_excluded
    print('Removed {0} loci that were too similar with other loci in the schema.'.format(len(final_excluded)))

    # decide which identifiers to keep and redo BLAST to check for residual paralogs
    output_schema = os.path.join(temp_directory, 'cd_hit_chewie.fasta')

    # filter protogenome to keep only representative sequences
    # optimize? seems to take too long!
    rfm.remove_same_locus_alleles(dna_file, excluded_alleles,
                                  protein_file, output_schema, minimum_cds_length)

    schema_dir = os.path.join(parent_directory, schema_name)

    # create directory and schema files
    rfm.build_schema(output_schema, schema_dir)

    if cleanup == 'yes':
        shutil.rmtree(temp_directory)

    end = time.time()
    delta = end - start
    print('Created schema based on {0} assemblies/genomes of {1} in {2}m{3}s.'.format(len(fasta_files),
                                                                  os.path.basename(chosen_taxon).rstrip('.trn'),
                                                                  int(delta/60), int(((delta/60)%1)*60)))

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

    parser.add_argument('--cd_hit_sim', '--cd_hit_similarity', type=float, required=False, dest='cd_hit_similarity',
                        default=0.6, help='Clustering similarity used to cluster proteins (default=0.6).')

    parser.add_argument('--cd_hit_word', '--cd_hit_word_size', type=int, required=False, dest='cd_hit_word_size',
                        default=4, help='Word size used in CD-HIT to cluster proteins (default=4).')

    parser.add_argument('--cutoff_sim', '--cutoff_similarity', type=float, required=False, dest='cutoff_similarity',
                        default=70.00, help='Similarity cutoff that will be used to consider that sequences in the same cluster'
                                            'correspond to alleles of the same gene (default=70.00).')

    parser.add_argument('--c', '--cleanup', type=str, required=False, dest='cleanup',
                        default='yes', help='If the temporary directory should be deleted at the end '
                                          '(default=yes.')

    args = parser.parse_args()

    return [args.input_files, args.output_directory, args.prodigal_training_file, 
            args.schema_name, args.cpu_count, args.blastp_path, args.blast_score_ratio,
            args.minimum_cds_length, args.cd_hit_similarity, args.cd_hit_word_size,
            args.cutoff_similarity, args.cleanup]


if __name__ == '__main__':

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4],
         args[5], args[6], args[7], args[8], args[9],
         args[10], args[11])
