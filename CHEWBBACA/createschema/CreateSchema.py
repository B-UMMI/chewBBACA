#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""


import os
import sys
import time
import pickle
import shutil
import argparse
import itertools
import datetime as dt
from collections import Counter
from multiprocessing import Pool

from Bio import SeqIO

try:
    from createschema import new_utils
    from utils import runProdigal
    from utils import auxiliary_functions as aux
except:
    from CHEWBBACA.createschema import new_utils
    from CHEWBBACA.utils import runProdigal
    from CHEWBBACA.utils import auxiliary_functions as aux


#input_files = '/home/rfm/Desktop/NS_tutorial_data/tutorial_data/sagalactiae_tutorial/sagalactiae_genomes/subset1'
#output_directory = '/home/rfm/Desktop/NS_tutorial_data/tutorial_data/sagalactiae_tutorial/sagalactiae_genomes/beep'
#ptf_path = '/home/rfm/Desktop/NS_tutorial_data/tutorial_data/sagalactiae_tutorial/sagalactiae_schema/Streptococcus_agalactiae.trn'
#schema_name = 'sagalactiae_test'
#cpu_count = 6
#blastp_path = '/home/rfm/Software/anaconda3/envs/ns/bin/blastp'
#blast_score_ratio = 0.6
#minimum_length = 201
#cleanup = False
#translation_table = 11
#size_threshold = 0.2
## cluster parameters
#clustering_mode = 'greedy'
#representative_filter = 0.80
#intra_filter = 0.80
#word_size = 4
#clustering_sim = 0.20
#cds_input = False
#verbose = False

def main(input_files, output_directory, schema_name, ptf_path, blast_score_ratio,
         minimum_length, translation_table, size_threshold, clustering_mode,
         word_size, clustering_sim, representative_filter, intra_filter, cpu_count,
         blastp_path, cds_input, verbose, cleanup):

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
    print('Started at: {0}\n'.format(start_date_str))

    cpu_to_apply = aux.verify_cpu_usage(cpu_count)

    print('\nCreating schema based on the genomes '
          'in the following directory:\n{0}'.format(os.path.abspath(input_files)))
    print('Training file: {0}'.format(ptf_path))
    print('Number of cores: {0}'.format(cpu_count))
    print('BLAST Score Ratio: {0}'.format(blast_score_ratio))
    print('Translation table: {0}'.format(translation_table))
    print('Minimum accepted sequence length: {0}'.format(minimum_length))
    print('Clustering mode: {0}'.format(clustering_mode))
    print('Word size: {0}'.format(word_size))
    print('Clustering similarity: {0}'.format(clustering_sim))
    print('Representative filter: {0}'.format(representative_filter))
    print('Intra filter: {0}'.format(intra_filter))

    # read file with paths to input files
    with open(input_files, 'r') as infile:
        fasta_files = [file.strip() for file in infile.readlines()]

    # maintain genome order to assign identifiers correctly
    fasta_files.sort(key=lambda y: y.lower())
    print('Number of genomes/assemblies: {0}'.format(len(fasta_files)))

    # determine and store genome identifiers
    genomes_identifiers = [aux.file_basename(g, False) for g in fasta_files]

    # define directory for temporary files
    temp_directory = os.path.join(output_directory, 'temp')

    # define output directory where Prodigal files will be stored
    prodigal_path = os.path.join(temp_directory, 'prodigal_cds_prediction')
    if not os.path.exists(prodigal_path):
        os.makedirs(prodigal_path)

    # run Prodigal to determine CDSs for all input genomes
    prodigal_start = time.time()
    print('\nPredicting CDS sequences...')
    print ('Started Prodigal at: {0}'.format(time.strftime('%H:%M:%S - %d/%m/%Y')))

    # run Prodigal with multiprocessing
    inputs = [[file, prodigal_path, ptf_path, translation_table, 'single'] for file in fasta_files]
    results = []
    completed = False
    tickval = 5
    ticknum = 20
    pool = Pool(cpu_to_apply)
    rawr = pool.map_async(runProdigal.main, inputs, callback=results.extend, chunksize=1)

    while completed is False:
        completed = aux.progress_bar(rawr, len(inputs), tickval, ticknum, completed)

    rawr.wait()

    print('\nChecking if Prodigal created all the necessary files...')
    fasta_files, genomes_identifiers = new_utils.check_prodigal_output_files(prodigal_path, fasta_files, input_files, results, genomes_identifiers, output_directory)
    print('Finished Prodigal at: {0}'.format(time.strftime("%H:%M:%S - %d/%m/%Y")))
    prodigal_end = time.time()
    prodigal_delta = prodigal_end - prodigal_start
    print('Prodigal delta: {0}'.format(prodigal_delta))

    # divide inputs
    inputs2 = aux.divide_list_into_n_chunks(fasta_files, cpu_to_apply)
    inputs2 = [i for i in inputs2 if len(i) > 0]
    for i in range(len(inputs2)):
        inputs2[i].append(prodigal_path)
        inputs2[i].append(output_directory)
        inputs2[i].append(i+1)

    # extract coding sequences
    print('Extracting coding sequences from all genomes...')
    cds_start = time.time()
    results2 = []
    pool = Pool(cpu_to_apply)
    rawr = pool.map_async(aux.process_cdss, inputs2, callback=results2.extend)
    rawr.wait()
    cds_end = time.time()
    cds_delta = cds_end - cds_start
    print('CDS extraction delta: {0}:'.format(cds_delta))

    table_files = [f[0] for f in results2]
    table_file = os.path.join(output_directory, 'protein_info.tsv')
    table_header = 'Genome\tContig\tStart\tStop\tProtein_ID\tCoding_Strand\n'
    aux.concatenate_files(table_files, table_file, table_header)
    for f in table_files:
        os.remove(f)

    cds_files = [f[1] for f in results2]
    inputs3 = []
    unique_files = []
    for i in range(len(cds_files)):
        inputs3.append([cds_files[i], os.path.join(temp_directory, 'unique_dna_seqids_{0}.fasta'.format(i+1))])
        unique_files.append(os.path.join(temp_directory, 'unique_dna_seqids_{0}.fasta'.format(i+1)))

    # determine repeated sequences and keep only one representative
    r1_start = time.time()
    results3 = []
    pool = Pool(cpu_to_apply)
    rawr = pool.map_async(aux.determine_repeated, inputs3, callback=results3.extend)
    rawr.wait()
    # one last round after concatenating files
    cds_file = os.path.join(temp_directory, 'coding_sequences_all.fasta')
    cds_file = aux.concatenate_files(unique_files, cds_file)
    unique_dna_seqids = os.path.join(temp_directory, 'unique_dna_seqids.fasta')
    repeated_dna_cds = aux.determine_repeated([cds_file, unique_dna_seqids])
    r1_end = time.time()
    r1_delta = r1_end - r1_start

    print('Repeated DNA sequences removal delta: {0}:'.format(r1_delta))
    print('\nRemoved {0} repeated DNA sequences.'.format(repeated_dna_cds[0]))

    # we should remove small DNA sequences from multiple files, then join and remove duplicated
    # sequences again

    # determine small DNA sequences and remove those seqids
    small_dna_cds = aux.determine_small(unique_dna_seqids, minimum_length)

    valid_dna_seqs = small_dna_cds[1]
    print('Removed {0} DNA sequences shorter than {1} nucleotides.'.format(len(small_dna_cds[0]), minimum_length))
    valid_dna_seqs.sort(key=lambda y: y.lower())

    # convert to protein
    print('\nTranslating {0} DNA sequences...'.format(len(valid_dna_seqs)))
    trans_start = time.time()
    results4 = []
    inputs4 = aux.divide_list_into_n_chunks(valid_dna_seqs, cpu_to_apply)
    dna_files = []
    protein_files = []
    for i in range(len(inputs4)):
        inputs4[i].append(unique_dna_seqids)
        inputs4[i].append(translation_table)
        inputs4[i].append(minimum_length)
        inputs4[i].append(os.path.join(temp_directory, 'valid_dna_{0}.fasta'.format(i+1)))
        inputs4[i].append(os.path.join(temp_directory, 'valid_protein_{0}.fasta'.format(i+1)))
    pool = Pool(cpu_to_apply)
    rawr = pool.map_async(aux.translate_coding_sequences, inputs4, callback=results4.extend)
    rawr.wait()

    # concatenate files
    dna_files = [i[-2] for i in inputs4]
    protein_files = [i[-1] for i in inputs4]
    dna_valid_file = os.path.join(temp_directory, 'valid_dna.fasta')
    dna_valid_file = aux.concatenate_files(dna_files, dna_valid_file)
    protein_valid_file = os.path.join(temp_directory, 'valid_protein.fasta')
    protein_valid_file = aux.concatenate_files(protein_files, protein_valid_file)

    untranslatable_cds = []
    for r in results4:
        untranslatable_cds.extend(r[0])

    untranslatable_seqids = [t[0] for t in untranslatable_cds]

    trans_end = time.time()
    trans_delta = trans_end - trans_start
    print('Translation delta: {0}:'.format(trans_delta))
    # write file with invalid alleles info
    invalid_alleles_file = os.path.join(output_directory, 'invalid_alleles.txt')
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
    unique_prots_file = os.path.join(temp_directory, 'unique_prots_seqs.fasta')
    repeated_protein_cds = aux.determine_repeated([protein_valid_file, unique_prots_file])

    # remove seqids that are from repeated protein sequences
    unique_protein_seqs = repeated_protein_cds[1]
    print('Removed {0} repeated Protein sequences.'.format(repeated_protein_cds[0]))
    unique_protein_seqs.sort(key=lambda y: y.lower())

    final_seqids = unique_protein_seqs
    final_seqids.sort(key=lambda y: y.lower())

    # write protein FASTA file
    protein_file = os.path.join(temp_directory, 'filtered_proteins.fasta')
    indexed_protein_valid_file = SeqIO.index(protein_valid_file, 'fasta')
    aux.get_sequences_by_id(indexed_protein_valid_file, final_seqids, protein_file)

    # write DNA FASTA file
    dna_file = os.path.join(temp_directory, 'filtered_dna.fasta')
    indexed_dna_valid_file = SeqIO.index(dna_valid_file, 'fasta')
    aux.get_sequences_by_id(indexed_dna_valid_file, final_seqids, dna_file)

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
    clusters = aux.cluster_sequences(sorted_prots, word_size,
                                     clustering_sim, clustering_mode, 1)
    print('Clustered proteins into {0} clusters'.format(len(clusters)))
    cluster_end = time.time()
    cluster_delta = cluster_end - cluster_start
    print('Cluster delta: {0}'.format(cluster_delta))

    # remove sequences that are very similar to representatives
    prunned_clusters, excluded_alleles = aux.cluster_prunner(clusters,
                                                             representative_filter)
    excluded_alleles = [e[0] for e in excluded_alleles]
    # determine clusters that only have the representative
    singletons = aux.determine_singletons(prunned_clusters)
    print('Singletons: {0}'.format(len(singletons)))
    # remove singletons and keep clusters that need to be BLASTed
    final_clusters = aux.remove_clusters(prunned_clusters, singletons)
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

    excluded_dict = aux.intra_cluster_sim(intra_clusters, indexed_prots, word_size, intra_filter)

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
    aux.get_sequences_by_id(indexed_prots, clustered_sequences2, clustered_seqs_file)

    integer_clusters = os.path.join(temp_directory, 'clustered_proteins_int.fasta')
    ids_dict = aux.integer_headers(clustered_seqs_file, integer_clusters)

    makedb_cmd = 'makeblastdb -in {0} -parse_seqids -dbtype prot >/dev/null 2>&1'.format(integer_clusters)
    #makedb_cmd = 'makeblastdb -in {0} -parse_seqids -dbtype prot'.format(integer_clusters)
    os.system(makedb_cmd)

    # BLAST necessary sequences against only the sequences from the same cluster
    blast_db = os.path.join(temp_directory, 'clustered_proteins_int.fasta')

    blast_results_dir = os.path.join(temp_directory, 'blast_results')
    os.mkdir(blast_results_dir)

    seqids_to_blast = aux.blast_inputs(final_clusters, blast_results_dir, ids_dict)

    # distribute clusters per available cores, try to group inputs into
    # even groups in terms of number of clusters and sum of number of sequences
    # per inputs group
    splitted_seqids = aux.split_blast_inputs_by_core(seqids_to_blast,
                                                     cpu_to_apply,
                                                     blast_results_dir)

    for s in range(len(splitted_seqids)):
        splitted_seqids[s].append(blastp_path)
        splitted_seqids[s].append(blast_db)
        splitted_seqids[s].append(blast_results_dir)
        splitted_seqids[s].append(integer_clusters)
    
    # create the FASTA files with the protein sequences before BLAST?
    blast_start = time.time()
    print('BLASTing protein sequences in each cluster...')
    # ADD progress bar!!!
    # BLAST each sequences in a cluster against every sequence in that cluster
    p = Pool(processes = cpu_to_apply)
    r = p.map_async(aux.cluster_blaster, splitted_seqids)
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
        splitted_results.append([file, indexed_fasta, blast_score_ratio, ids_dict])

    blast_excluded_alleles = []
    for i in splitted_results:
        try:
            e = aux.apply_bsr(i)
            blast_excluded_alleles.append(e)
        # proteins with low complexity regions might
        # fail to align 
        except Exception:
            print('Could not apply BSR to {0} results. Probably low complexity proteins.'.format(i))

    # merge bsr results
    blast_excluded_alleles = list(itertools.chain.from_iterable(blast_excluded_alleles))
    blast_excluded_alleles = [ids_dict[seqid] for seqid in blast_excluded_alleles]
    excluded_alleles.extend(blast_excluded_alleles)

    # perform final BLAST to avoid creating a schema with paralogs
    print('Performing a final BLAST to check for paralogs...')
    schema_seqids = list(set(final_seqids) - set(excluded_alleles))
    beta_file = os.path.join(temp_directory, 'beta_schema.fasta')
    aux.get_sequences_by_id(indexed_protein_valid_file, schema_seqids, beta_file)

    integer_seqids = os.path.join(temp_directory, 'int_proteins_int.fasta')
    ids_dict2 = aux.integer_headers(beta_file, integer_seqids)

    makedb_cmd = 'makeblastdb -in {0} -dbtype prot >/dev/null 2>&1'.format(integer_seqids)
    os.system(makedb_cmd)

    blast_db = os.path.join(temp_directory, integer_seqids)
    
    second_blast_start = time.time()
    blast_output = '{0}/{1}_blast_out.tsv'.format(temp_directory, 'beta_schema')
    blast_command = ('{0} -db {1} -query {2} -out {3} -outfmt "6 qseqid sseqid score" '
                     '-max_hsps 1 -num_threads {4} -evalue 0.001'.format(blastp_path, integer_seqids, integer_seqids, blast_output, cpu_to_apply))
    os.system(blast_command)
    second_blast_end = time.time()
    second_blast_delta = second_blast_end - second_blast_start
    print('Total BLAST time: {0}'.format(second_blast_delta))

    final_excluded = aux.apply_bsr([blast_output, indexed_fasta, blast_score_ratio, ids_dict2])
    final_excluded = [ids_dict2[seqid] for seqid in final_excluded]
    excluded_alleles = excluded_alleles + final_excluded
    schema_seqids_final = list(set(schema_seqids) - set(excluded_alleles))
    print('Removed {0} loci that were too similar with other loci in the schema.'.format(len(final_excluded)))

    # decide which identifiers to keep and redo BLAST to check for residual paralogs
    output_schema = os.path.join(temp_directory, 'schema_seed.fasta')

    # create file with the schema representative sequences
    aux.get_schema_reps(indexed_fasta, schema_seqids_final,
                        protein_file, output_schema)

    schema_dir = os.path.join(output_directory, schema_name)
    if not os.path.exists(schema_dir):
        os.makedirs(schema_dir)

    # create directory and schema files
    aux.build_schema(output_schema, schema_dir)

    # remove temporary files
    if cleanup is True:
        shutil.rmtree(temp_directory)

    end_date = dt.datetime.now()
    end_date_str = dt.datetime.strftime(end_date, '%Y-%m-%dT%H:%M:%S')

    delta = end_date - start_date
    minutes, seconds = divmod(delta.total_seconds(), 60)

    print('\nFinished at: {0}'.format(end_date_str))
    print('Created schema based on {0} genomes in'
          '{1: .0f}m{2: .0f}s.'.format(len(fasta_files),
                                       minutes, seconds))


def parse_arguments():


    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', nargs='?', type=str, required=True,
                        dest='input_files',
                        help='Path to the directory that contains the input '
                             'FASTA files. Alternatively, a single file with '
                             'a list of paths to FASTA files, one per line')

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
            args.representative_filter, args.intra_filter, args.cpu_count,
            args.blastp_path, args.cds_input, args.verbose, args.cleanup]


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2], args[3], args[4], args[5],
         args[6], args[7], args[8], args[9], args[10], args[11],
         args[12], args[13], args[14], args[15], arg[16], args[17])
