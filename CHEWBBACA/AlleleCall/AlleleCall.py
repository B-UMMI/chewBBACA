#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module performs allele calling to determine the allelic profiles
for a set of bacterial strains.

Code documentation
------------------
"""


import os
import csv
import sys
import math
from collections import Counter

from Bio import SeqIO

try:
    from utils import (constants as ct,
                       blast_wrapper as bw,
                       profile_hasher as ph,
                       core_functions as cf,
                       file_operations as fo,
                       fasta_operations as fao,
                       process_datetime as pdt,
                       sequence_manipulation as sm,
                       iterables_manipulation as im,
                       multiprocessing_operations as mo)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 blast_wrapper as bw,
                                 profile_hasher as ph,
                                 core_functions as cf,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 process_datetime as pdt,
                                 sequence_manipulation as sm,
                                 iterables_manipulation as im,
                                 multiprocessing_operations as mo)


def loci_mode_file(loci_files, output_file):
    """Determine the allele size mode for a set of loci.

    Parameters
    ----------
    loci_files : list
        List with the full paths to loci FASTA files.
    output_file : str
        Path to the output file created to store the sequence size mode
        values (created with the Pickle module).

    Returns
    -------
    loci_modes : str
        Path to the output file with the allele size mode values (a
        dictionary with loci identifiers as keys and the allele size
        mode and the list of allele sizes as values is saved with the
        Pickle module).
    """
    loci_modes = {}
    for file in loci_files:
        locus_id = fo.file_basename(file, False)
        allele_sizes = list(fao.sequence_lengths(file).values())
        # select first if there are several values with same frequency
        loci_modes[locus_id] = [sm.determine_mode(allele_sizes)[0], allele_sizes]

    fo.pickle_dumper(loci_modes, output_file)

    return loci_modes


def allele_hash_table(fasta_files, table_id, translation_table,
                      output_directory, file_prefix):
    """Create a hash table with allele hashes mapped to allele identifiers.

    Parameters
    ----------
    fasta_files : list
        List with paths to loci FASTA files.
    table_id : int
        Unique integer identifier added to the name of the output file.
    translation_table : int
        Genetic code used to translate the alleles if the hashes should be
        determined based on the aminoacid sequences. Will use the DNA
        sequences if this argument is not provided.
    output_directory : str
        Path to the output directory.
    file_prefix : str
        Prefix added to the name of the file before the unique integer
        identifier.

    Returns
    -------
    table_file : str
        Path to the output file (a dictionary with allele hashes as keys
        and the locus integer index and allele identifier as values is saved
        with the Pickle module).
    """
    hashtable = {}
    for file in fasta_files:
        locus_index = file[1]
        alleles = [(rec.id, str(rec.seq))
                   for rec in fao.sequence_generator(file[0])]
        if translation_table is not None:
            alleles = [(rec[0], str(sm.translate_sequence(rec[1], translation_table)))
                       for rec in alleles]

        for record in alleles:
            # needs to be string because of '*' added to schemas from Chewie-NS
            allele_id = record[0].split('_')[-1]
            allele_hash = im.hash_sequence(record[1])
            hashtable.setdefault(allele_hash, []).append((locus_index, allele_id))

    table_file = fo.join_paths(output_directory, ['{0}{1}'.format(file_prefix, table_id)])
    fo.pickle_dumper(hashtable, table_file)

    return table_file


def pre_compute_hash_tables(output_directory, loci_files, translation_table,
                            cpu_cores):
    """Create allele and translated allele hash tables.

    Parameters
    ----------
    output_directory : str
        Path to the output directory.
    loci_files : list
        List with paths to loci FASTA files.
    translation_table : int
        Genetic code used to translate the alleles if the hashes should be
        determined based on the aminoacid sequences. Will use the DNA
        sequences if this argument is not provided.
    cpu_cores : int
        Number of CPU cores used to compute the hash tables.

    Returns
    -------
    hash_tables : list
        List with the paths for the files that contain the hash tables (
        files are created in the output diretory with the Pickle module).
    """
    # define maximum number of sequence hashes per hash table
    max_sequences = ct.HASH_TABLE_MAXIMUM_ALLELES

    input_groups = []
    current_group = []
    total_alleles = 0
    for file in loci_files:
        # count number of sequences per file
        with open(file, 'r') as infile:
            num_alleles = sum(1 for _ in infile) / 2
        # increment total sequences and reset if it reaches limit
        total_alleles += num_alleles
        current_group.append((file, loci_files.index(file)))
        if total_alleles >= max_sequences or file == loci_files[-1]:
            input_groups.append(current_group)
            current_group = []
            total_alleles = 0

    # create inputs to parallelize hash table creation
    # hash tables for DNA sequences
    inputs = [[g, i+1, None, output_directory, 'DNAtable', allele_hash_table]
              for i, g in enumerate(input_groups)]

    hash_tables = mo.map_async_parallelizer(inputs,
                                            mo.function_helper,
                                            cpu_cores,
                                            show_progress=False)

    # hash tables for translated alleles
    inputs = [[g, i+1, translation_table, output_directory, 'PROTEINtable', allele_hash_table]
              for i, g in enumerate(input_groups)]

    hash_tables = mo.map_async_parallelizer(inputs,
                                            mo.function_helper,
                                            cpu_cores,
                                            show_progress=False)

    return hash_tables


def update_hash_tables(loci_files, loci_to_call, translation_table,
                       pre_computed_dir):
    """Update pre-computed hash tables.

    Parameters
    ----------
    loci_files : dict
        Dictionary with paths to schema loci FASTA files as keys and
        a list with the paths to the temporary FASTA files that contain
        the newly inferred alleles for the loci as values.
    loci_to_call : dict
        Dictionary with paths to schema loci FASTA files as keys and
        loci integer identifiers as values.
    translation_table : int
        Genetic code used to translate the alleles.
    pre_computed_dir : str
        Path to the directory that contains the files with the pre-computed
        hash tables.

    Returns
    -------
    Total number of new allele hashes added to the pre-computed hash tables.
    """
    # create hash tables with data for new alleles
    novel_dnatable = {}
    novel_proteintable = {}
    for key, value in loci_files.items():
        # get locus integer identifier
        locus_index = loci_to_call[key]
        # import new alleles
        records = fao.sequence_generator(value[0])
        for rec in records:
            # compute allele SHA256 hash
            allele_id = (rec.id).split('_')[-1]
            sequence = str(rec.seq)
            seq_hash = im.hash_sequence(sequence)
            prot_hash = im.hash_sequence(str(sm.translate_sequence(sequence, translation_table)))
            # add to hash tables that will be used to update pre-computed data
            novel_dnatable.setdefault(seq_hash, []).append((locus_index, allele_id))
            novel_proteintable.setdefault(seq_hash, []).append(prot_hash)

    # update pre-computed hash tables
    # list files with pre-computed hash tables and select last created
    dna_tables = fo.listdir_fullpath(pre_computed_dir, 'DNAtable')
    prot_tables = fo.listdir_fullpath(pre_computed_dir, 'PROTEINtable')
    # load hash tables to update
    latest_dna_table = sorted(dna_tables,
                              key=lambda x: int(x.split('table')[-1]))[-1]
    latest_prot_table = sorted(prot_tables,
                               key=lambda x: int(x.split('table')[-1]))[-1]
    current_dna_table = fo.pickle_loader(latest_dna_table)
    current_prot_table = fo.pickle_loader(latest_prot_table)
    for key, value in novel_dnatable.items():
        # add entry for each new allele
        current_dna_table.setdefault(key, []).extend(value)
        current_prot_table.setdefault(novel_proteintable[key][0], []).extend(value)
        # check if it reached the maximum number of alleles per table file
        if len(current_dna_table) >= ct.HASH_TABLE_MAXIMUM_ALLELES:
            # save current hash tables and create new files
            fo.pickle_dumper(current_dna_table, latest_dna_table)
            fo.pickle_dumper(current_prot_table, latest_prot_table)
            new_index = int(latest_dna_table.split('DNAtable')[-1]) + 1
            latest_dna_table = fo.join_paths(pre_computed_dir,
                                             ['DNAtable{0}'.format(new_index)])
            current_dna_table = {}
            latest_prot_table = fo.join_paths(pre_computed_dir,
                                              ['PROTEINtable{0}'.format(new_index)])
            current_prot_table = {}

    if len(current_dna_table) > 0:
        fo.pickle_dumper(current_dna_table, latest_dna_table)
        fo.pickle_dumper(current_prot_table, latest_prot_table)

    return len(novel_dnatable)


def create_classification_file(locus_id, output_directory, locus_results):
    """Create file to store classifications for a locus.

    Parameters
    ----------
    locus_id : str
        The identifier of the locus.
    output_directory : str
        Path to the output directory where the file will
        be created.
    locus_results : dict
        Results to save to file.

    Return
    ------
    pickle_out : str
        Path to the file created to store classification
        results.
    """
    pickle_out = fo.join_paths(output_directory,
                               [locus_id+'_results'])

    # create file with empty results structure
    fo.pickle_dumper(locus_results, pickle_out)

    return pickle_out


def update_classification(genome_id, locus_results, match_info):
    """Update locus classification for an input.

    Parameters
    ----------
    genome_id : int
        Integer identifier attributed to the input.
    locus_results : dict
        Dictionary with the matches found for the locus
        in the inputs.
    match_info : list
        List with information about the match found for
        the locus.

    Returns
    -------
    locus_results : dict
        Updated results.
    """
    # add data about match
    locus_results.setdefault(genome_id, [match_info[3]]).append(match_info)

    # get all classifications
    classes_list = [c[3] for c in locus_results[genome_id][1:]]
    # evaluate classification for genomes with multiple matches
    if len(classes_list) > 1:
        classes_counts = Counter(classes_list)
        # multiple matches, single class
        if len(classes_counts) == 1:
            if 'EXC' in classes_counts:
                locus_results[genome_id][0] = 'NIPHEM'
            # multiple INF, ASM, ALM, etc classes are classified as NIPH
            else:
                locus_results[genome_id][0] = 'NIPH'
        # multiple matches and classes
        elif len(classes_counts) > 1:
            # inputs that include both EXC and INF are classified as NIPH
            if 'EXC' and 'INF' in classes_counts:
                locus_results[genome_id][0] = 'NIPH'
            # any class with PLOT3, PLOT5 or LOTSC are classified as NIPH
            elif any([c in ['PLOT3', 'PLOT5', 'LOTSC'] for c in classes_counts]) is True:
                locus_results[genome_id][0] = 'NIPH'
            # EXC or INF with ASM/ALM
            elif 'EXC' in classes_counts or 'INF' in classes_counts:
                match_count = classes_counts.get('EXC', classes_counts['INF'])
                # Single EXC or INF classified as EXC or INF even if there are ASM/ALM
                if match_count == 1:
                    locus_results[genome_id][0] = 'EXC' if 'EXC' in classes_counts else 'INF'
                # multiple EXC or INF classified as NIPH
                else:
                    locus_results[genome_id][0] = 'NIPH'
            # multiple ASM and ALM are classified as NIPH
            else:
                locus_results[genome_id][0] = 'NIPH'

    return locus_results


def count_classifications(classification_files, classification_labels):
    """Determine counts for each classification type except LNF.

    Parameters
    ----------
    classification files : list
        List of paths to pickled files that contain the
        classifications for a set of loci.

    Returns
    -------
    global_counts : dict
        Dicitonary with classification types as keys
        and the total number of inputs classified per
        type as values.
    total_cds : int
        The total number of coding sequences that
        have been classified.
    """
    classification_counts = Counter()
    # get total number of classified CDSs
    total_cds = 0
    for file in classification_files:
        locus_results = fo.pickle_loader(file)
        total_cds += sum([len([r for r in c if isinstance(r, tuple)])
                          for g, c in locus_results.items()])
        locus_classifications = [c[0] for g, c in locus_results.items()]
        locus_counts = Counter(locus_classifications)
        classification_counts += locus_counts

    # add classification that might be missing
    classification_counts.update(Counter({k: 0 for k in classification_labels[:-1]
                                          if k not in classification_counts}))

    return [classification_counts, total_cds]


def dna_exact_matches(table_file, presence_dnahashtable, loci_files,
                      classification_files, input_ids):
    """Find exact matches between input CDSs and alleles from a locus.

    Parameters
    ----------
    locus_file : str
        Path to the locus FASTA file that contains the locus
        alleles.
    presence_dnahashtable : dict
        Dictionary with SHA-256 hashes for distinct DNA
        sequences extracted from the inputs and lists of
        genome integer identifiers enconded with the
        polyline algorithm as values.
    locus_classifications : str
        Dictionary with the matches found for the locus
        in the inputs.
    input_ids : dict
        Dictionary with input integer identifiers as keys
        and input sequence identifiers as values.

    Returns
    -------
    locus_classifications : dict
        Updated results with the exact matches found for the
        locus.
    matched_seqids : list
        Sequence identifiers of the distinct CDSs that were
        matched.
    total_matches : int
        Number of exact matches.
    """
    variant_table = fo.pickle_loader(table_file)
    matches = set(variant_table.keys()).intersection(set(presence_dnahashtable.keys()))
    matches_data = {m: variant_table[m] for m in matches}
    # sort based on locus integer identifier
    sorted_matches = sorted(list(matches_data.items()), key=lambda x: x[1][0][0])
    groups = {}
    for match in sorted_matches:
        # only get information about the first locus that contains alleles
        # does not get additional info if allele is present in several loci
        groups.setdefault(match[1][0][0], []).append((match[0], match[1][0][1]))

    groups = {loci_files[k]: v
              for k, v in groups.items()
              if k in loci_files}

    matched_seqids = []
    dna_exact_hits = 0
    dna_matches_ids = 0
    for g, results in groups.items():
        class_file = classification_files[g]
        locus_classifications = fo.pickle_loader(class_file)

        locus_matched_seqids = []
        locus_total_matches = 0
        for r in results:
            # decode list of inputs that contain allele
            matched_inputs = im.polyline_decoding(presence_dnahashtable[r[0]])
            # seqid chosen as representative during sequence deduplication
            representative_seqid = '{0}-protein{1}'.format(input_ids[matched_inputs[1]], matched_inputs[0])
            match_data = (r[1], representative_seqid, r[0], 'EXC', 1.0)
            # classify as exact matches
            # skip first value, it is the protein id
            for gid in matched_inputs[1:]:
                locus_classifications = update_classification(gid,
                                                              locus_classifications,
                                                              match_data)

            locus_total_matches += len(matched_inputs[1:])
            # store representative id for the sequences
            locus_matched_seqids.append(representative_seqid)

        # save updated classifications
        fo.pickle_dumper(locus_classifications, class_file)
        # extend list of matched seqids
        matched_seqids.extend(locus_matched_seqids)
        dna_exact_hits += locus_total_matches
        dna_matches_ids += len(locus_matched_seqids)

    return [matched_seqids, dna_exact_hits, dna_matches_ids]


def protein_exact_matches(table_file, presence_PROThashtable, loci_files,
                          classification_files, input_ids, dna_index,
                          presence_DNAhashtable):
    """Find exact matches at protein level.

    Parameters
    ----------
    locus_file : str
        Path to the locus FASTA file that contains the locus
        translated alleles.
    presence_PROThashtable : dict
        Dictionary with SHA-256 hashes for distinct protein
        sequences extracted from the inputs and lists of
        sequence identifiers enconded with the polyline
        algorithm as values.
    presence_DNAhashtable : dict
        Dictionary with SHA-256 hashes for distinct DNA
        sequences extracted from the inputs and lists of
        genome integer identifiers enconded with the
        polyline algorithm as values.
    locus_classifications : dict
        Dictionary with the matches found for the locus
        in the inputs.
    dna_index : Bio.File._IndexedSeqFileDict
        Fasta file index created with BioPython.
    input_ids : dict
        Dictionary with input integer identifiers as keys
        and input sequence identifiers as values.

    Returns
    -------
    locus_classifications : dict
        Updated results with the exact matches found for the
        locus.
    exact_prot_hashes : list
        Sequence identifiers of the distinct CDSs that
        were matched.
    total_prots : int
        Number of matched distinct CDSs.
    total_cds : int
        Number of matched CDSs.
    total_distinct_prots : int
        Number of matched distinct prots.
    """
    variant_table = fo.pickle_loader(table_file)
    common = set(variant_table.keys()).intersection(set(presence_PROThashtable.keys()))
    common_data = {c: variant_table[c] for c in common}
    sorted_common = sorted(list(common_data.items()), key=lambda x: x[1][0][0])
    groups = {}
    for r in sorted_common:
        # only gets information about the first locus that contains alleles
        # does not get addiitonal info if allele is present in several loci
        groups.setdefault(r[1][0][0], []).append((r[0], r[1][0][1]))

    groups = {loci_files[k]: v for k, v in groups.items() if k in loci_files}

    exact_phashes = []
    exc_prot = 0
    exc_cds = 0
    exc_distinct_prot = 0
    loci_modes = []
    for g, results in groups.items():
        class_file = classification_files[g]
        locus_classifications = fo.pickle_loader(class_file)

        locus_total_cds = 0
        locus_total_prots = 0
        locus_total_distinct_prots = 0
        locus_matched_dna = {}
        locus_matched_proteins = set()
        locus_exact_prot_hashes = []
        locus_lengths = []
        for r in results:
            if r[0] in presence_PROThashtable:
                # get protids for distinct DNA CDSs
                matched_protids = im.polyline_decoding(presence_PROThashtable[r[0]])
                matched_protids = ['{0}-protein{1}'.format(input_ids[matched_protids[i+1]], matched_protids[i])
                                   for i in range(0, len(matched_protids), 2)]
                locus_total_prots += len(matched_protids)
                locus_exact_prot_hashes.extend(matched_protids)
                locus_total_distinct_prots += 1
                # for each distinct CDS that codes for the protein
                for m in matched_protids:
                    cds = str(dna_index.get(m).seq)
                    cds_hash = im.hash_sequence(cds)
                    # get IDs of genomes that contain the CDS
                    matched_inputs = im.polyline_decoding(presence_DNAhashtable[cds_hash])
                    locus_total_cds += len(matched_inputs)-1
                    # for each genome ID that contains the CDS
                    for gid in matched_inputs[1:]:
                        # first time seeing CDS
                        if cds_hash not in locus_matched_dna:
                            current_class = 'INF'
                            # representative_seqid = r[1]
                            locus_matched_dna[cds_hash] = m
                            locus_lengths.append(len(cds))
                        else:
                            current_class = 'EXC'
                            # if it matches a INF, change protid to seqid of INF
                            # that will be assigned a new allele id later
                            # representative_seqid = locus_matched_dna[cds_hash]
                        # for protein exact matches, the seqid of the translated allele,
                        # the seqid of the protein chosen as representative during sequence deduplication,
                        match_data = (r[1], m, cds_hash, current_class, 1.0)
                        locus_classifications = update_classification(gid, locus_classifications,
                                                                      match_data)
                locus_matched_proteins.add(r[0])

        fo.pickle_dumper(locus_classifications, class_file)
        exact_phashes.extend(locus_exact_prot_hashes)
        exc_prot += locus_total_prots
        exc_cds += locus_total_cds
        exc_distinct_prot += locus_total_distinct_prots
        # update locus mode
        if len(locus_lengths) > 0:
            locus_id = fo.file_basename(g, False)
            loci_modes.append([locus_id, locus_lengths])

    return [exact_phashes, exc_prot, exc_cds, exc_distinct_prot, loci_modes]


def contig_position_classification(representative_length, representative_leftmost_pos,
                                   representative_rightmost_pos, contig_length,
                                   contig_leftmost_pos, contig_rightmost_pos):
    """Determine classification based on the alignment position on the contig.

    Parameters
    ----------
    representative_length : int
        Length of the representative allele that matched a
        coding sequence identified in the input contig.
    representative_leftmost_pos : int
        Representative sequence leftmost aligned position.
    representative_rightmost_pos : int
        Representative sequence rightmost aligned position.
    contig_length : int
        Length of the contig that contains the coding sequence
        that matched with the representative allele.
    contig_leftmost_pos : int
        Contig leftmost aligned position.
    contig_rightmost_pos : int
        Contig rightmost aligned position.

    Returns
    -------
    'LOTSC' if the contig is smaller than the matched representative
    allele, 'PLOT5' or 'PLOT3' if the matched allele unaligned part
    exceeds one of the contig ends, None otherwise.
    """

    # check if it is LOTSC because the contig is smaller than matched allele
    if contig_length < representative_length:
        return 'LOTSC'

    # check if it is PLOT
    # match in sense strand
    if contig_rightmost_pos > contig_leftmost_pos:
        # determine rightmost aligned position in contig
        contig_rightmost_rest = contig_length - contig_rightmost_pos
        # determine leftmost aligned position in contig
        contig_leftmost_rest = contig_leftmost_pos
        # determine number of rightmost bases in the target that did not align
        representative_rightmost_rest = representative_length - representative_rightmost_pos
        # determine number of leftmost bases in the target that did not align
        representative_leftmost_rest = representative_leftmost_pos
    # reverse values because CDS was identified in reverse strand
    elif contig_rightmost_pos < contig_leftmost_pos:
        contig_leftmost_rest = contig_rightmost_pos
        contig_rightmost_rest = contig_length - contig_leftmost_pos
        # also need to reverse values for representative
        representative_leftmost_rest = representative_rightmost_pos
        representative_rightmost_rest = representative_length - representative_leftmost_pos

    # check if the unaligned region of the matched allele exceeds
    # one of the contig ends
    if contig_leftmost_rest < representative_leftmost_rest:
        return 'PLOT5'
    elif contig_rightmost_rest < representative_rightmost_rest:
        return 'PLOT3'


def allele_size_classification(sequence_length, locus_mode, size_threshold):
    """Determine classification based on sequence size.

    Parameters
    ----------
    sequence_length : int
        Length of the DNA sequence.
    locus_mode : int
        Locus allele size mode.
    size_threshold : float
        Sequence size variation threshold.

    Returns
    -------
    'ASM' if sequence size value is below computed sequence
    size interval, 'ALM' if it is above and None if it is
    contained in the interval.
    """
    if size_threshold is not None:
        if sequence_length < (locus_mode[0]-(locus_mode[0])*size_threshold):
            return 'ASM'
        elif sequence_length > (locus_mode[0]+(locus_mode[0])*size_threshold):
            return 'ALM'


def write_loci_summary(classification_files, output_directory, total_inputs,
                       classification_labels):
    """Write a TSV file with classification counts per locus.

    Parameters
    ----------
    classification_files : dict
        Dictionary with the paths to loci FASTA files as keys
        and paths to loci classification files as values.
    output_directory : str
        Path to the output directory where the TSV file will
        be created.
    total_inputs : int
        Total number of inputs.
    """
    loci_stats = [ct.LOCI_STATS_HEADER]
    for k, v in classification_files.items():
        locus_id = fo.get_locus_id(k)
        if locus_id is None:
            locus_id = fo.file_basename(v).split('_results')[0]
        locus_results = fo.pickle_loader(v)

        # count locus classifications
        current_counts = count_classifications([v], classification_labels)
        counts_list = [locus_id]
        for c in classification_labels[:-1]:
            counts_list.append(str(current_counts[0][c]))
        # add LNF count
        counts_list.append(str(total_inputs-len(locus_results)))
        counts_list.append(str(current_counts[1]))
        locus_line = im.join_list(counts_list, '\t')
        loci_stats.append(locus_line)

    output_file = fo.join_paths(output_directory, [ct.LOCI_STATS_BASENAME])
    fo.write_lines(loci_stats, output_file)


def write_logfile(start_time, end_time, total_inputs,
                  total_loci, cpu_cores, blast_score_ratio,
                  output_directory):
    """Write the log file.

    Parameters
    ----------
    start_time : datetime.datetime
        Datetime object with the date and hour
        determined when the process started running.
    end_time : datetime.datetime
        Datetime object with the date and hour
        determined when the process concluded.
    total_inputs : int
        Number of inputs passed to the process.
    cpu_cores : int
        Number of CPU cores/threads used by the
        process.
    blast_score_ratio : float
        BLAST Score Ratio value used by the
        process.
    output_directory : str
        Path to the output directory where the
        log file will be created.

    Returns
    -------
    log_outfile : str
        Path to the log file.
    """
    start_time_str = pdt.datetime_str(start_time,
                                      date_format='%H:%M:%S-%d/%m/%Y')

    end_time_str = pdt.datetime_str(end_time,
                                    date_format='%H:%M:%S-%d/%m/%Y')

    log_outfile = fo.join_paths(output_directory, [ct.LOGFILE_BASENAME])
    logfile_text = ct.LOGFILE_TEMPLATE.format(start_time_str, end_time_str,
                                              total_inputs, total_loci,
                                              cpu_cores, blast_score_ratio)

    fo.write_to_file(logfile_text, log_outfile, 'w', '')

    return log_outfile


def write_results_alleles(classification_files, input_identifiers,
                          output_directory, missing_class):
    """Write a TSV file with the allelic profiles for the input samples.

    Parameters
    ----------
    classification_files : list
        List with the paths to loci classification files.
    input_identifiers : list
        Sorted list that contains input string identifiers.
    output_directory : str
        Path to the output directory.
    """
    # add first column with input identifiers
    columns = [['FILE'] + input_identifiers]
    for file in classification_files:
        # get locus identifier to add as column header
        locus_id = fo.get_locus_id(file)
        if locus_id is None:
            locus_id = fo.file_basename(file).split('_results')[0]
        locus_results = fo.pickle_loader(file)
        locus_column = [locus_id]
        for i in range(1, len(input_identifiers)+1):
            # determine if locus was found in each input
            if i in locus_results:
                current_result = locus_results[i]
                # exact or inferred, append assigned allele id
                if current_result[0] in ['EXC', 'INF']:
                    locus_column.append(current_result[-1])
                # missing data (PLOT, ASM, ALM, ...)
                else:
                    locus_column.append(current_result[0])
            # locus was not identified in the input
            else:
                locus_column.append(missing_class)

        columns.append(locus_column)

    # group elements with same list index
    lines = im.aggregate_iterables(columns)
    lines = ['\t'.join(line) for line in lines]

    output_file = fo.join_paths(output_directory, [ct.RESULTS_ALLELES_BASENAME])
    fo.write_lines(lines, output_file)

    return output_file


def write_results_statistics(classification_files, input_identifiers,
                             output_directory, classification_labels):
    """Write a TSV file with classification counts per input.

    Parameters
    ----------
    classification_files : dict
        Dictionary with the paths to loci FASTA files as keys
        and paths to loci classification files as values.
    input_identifiers : dict
        Dictionary with input integer identifiers as keys
        and input string identifiers as values.
    output_directory : str
        Path to the output directory where the TSV file will
        be created.
    """
    # initialize classification counts per input
    class_counts = {i: {c: 0 for c in classification_labels}
                    for i in input_identifiers}
    for file in classification_files.values():
        locus_id = fo.get_locus_id(file)
        if locus_id is None:
            locus_id = fo.file_basename(file).split('_results')[0]
        locus_results = fo.pickle_loader(file)

        for i in class_counts:
            if i in locus_results:
                class_counts[i][locus_results[i][0]] += 1
            else:
                class_counts[i][classification_labels[-1]] += 1

    # substitute integer identifiers by string identifiers
    class_counts = {input_identifiers[i]: v for i, v in class_counts.items()}

    # initialize with header line
    lines = [['FILE'] + classification_labels]
    for k, v in class_counts.items():
        input_line = [k] + [str(v[c]) for c in classification_labels]
        lines.append(input_line)

    outlines = ['\t'.join(line) for line in lines]

    output_file = fo.join_paths(output_directory, ['results_statistics.tsv'])
    fo.write_lines(outlines, output_file)


def write_results_contigs(classification_files, input_identifiers,
                          output_directory, cds_coordinates_files,
                          classification_labels):
    """Write a TSV file with the CDS coordinates for each input.

    Writes a TSV file with coding sequence coordinates (contig
    identifier, start and stop positions and coding strand) for
    EXC and INF classifications or with the classification type
    if it is not EXC or INF.

    Parameters
    ----------
    classification_files : list
        List with the paths to loci classification files.
    input_identifiers : dict
        Dictionary with input integer identifiers as keys
        and input string identifiers as values.
    output_directory : str
        Path to the output directory where the TSV file will
        be created.
    cds_coordinates_files : dict
        Dictionary with input string identifiers as keys
        and paths to pickled files with coding sequence
        coordinates as values.
    classification_labels : list
        List with the class labels attributed by chewBBACA.

    Returns
    -------
    output_file : str
        Path to the output file that contains the sequence
        coordinates.
    repeated : dict
        Dictionary with hashes for the CDSs that matched multiple loci
        as keys and a list with information about the genome of origin
        and matched loci as values.
    """
    # do not include LNF class
    invalid_classes = classification_labels[2:-1]
    intermediate_file = fo.join_paths(output_directory,
                                      ['inter_results_contigsInfo.tsv'])
    columns = [['FILE'] + list(input_identifiers.values())]
    # limit the number of lines to store in memory
    line_limit = 500
    # get hash if coordinates are available, seqid otherwise
    id_index = 2 if cds_coordinates_files is not None else 1
    for i, file in enumerate(classification_files):
        locus_id = fo.get_locus_id(file)
        if locus_id is None:
            locus_id = fo.file_basename(file).split('_results')[0]
        locus_results = fo.pickle_loader(file)
        column = [locus_id]
        # get sequence hash for exact and inferred
        # get classification for other cases or LNF for no classification
        column += [locus_results[i][1][id_index]
                   if i in locus_results and locus_results[i][0] not in invalid_classes
                   else locus_results.get(i, [classification_labels[-1]])[0]
                   for i in input_identifiers]

        columns.append(column)

        if len(columns) >= line_limit or (i+1) == len(classification_files):
            inter_lines = [im.join_list(c, '\t') for c in columns]
            fo.write_lines(inter_lines, intermediate_file, write_mode='a')
            columns = []

    # transpose intermediate file
    transposed_file = fo.transpose_matrix(intermediate_file, output_directory)

    # include LNF class in list to exclude from counting
    # and coordinate retrieval
    invalid_classes.append(classification_labels[-1])

    # get all hashes that match more than one locus per input
    # alleles that match more than one locus are only identified if they are
    # classified as EXC or INF for both loci
    # if they are classified as EXC/INF for one locus and NIPH for the other,
    # they will not be identified based on results for that genome
    # we need to get all alleles that match multiple loci before deciding if
    # they should be added to the schema.
    repeated = {}
    with open(transposed_file, 'r') as infile:
        csv_reader = csv.reader(infile, delimiter='\t')
        header = csv_reader.__next__()
        for i, l in enumerate(csv_reader):
            genome_id = l[0]
            # count number of ocurrences for each hash/seqid
            hashes = [h for h in l[1:] if h not in invalid_classes]
            hash_counts = Counter(hashes)
            repeated_hashes = [count
                               for count in hash_counts.most_common()
                               if count[1] > 1]
            for h in repeated_hashes:
                repeated.setdefault(h[0], [])

    # create final output file
    output_file = fo.join_paths(output_directory, ['results_contigsInfo.tsv'])
    with open(transposed_file, 'r') as infile:
        csv_reader = csv.reader(infile, delimiter='\t')
        header = csv_reader.__next__()
        output_lines = [header]
        # convert hashes to CDS coordinates
        if cds_coordinates_files is not None:
            for i, l in enumerate(csv_reader):
                genome_id = l[0]
                # open file with loci coordinates
                coordinates = fo.pickle_loader(cds_coordinates_files[genome_id])[0]
                # start position is 0-based, stop position is upper-bound exclusive
                # convert to PAMA CDSs that matched multiple loci
                cds_coordinates = []
                for j, c in enumerate(l[1:]):
                    if c not in repeated:
                        cds_coordinates.append(coordinates.get(c, [c])[0])
                    else:
                        cds_coordinates.append(classification_labels[-2])
                        # add CDS coordinates to list of PAMA classifications
                        repeated_coordinates = coordinates[c][0]
                        repeated[c].append([genome_id, header[j+1],
                                            '{0}&{1}-{2}&{3}'.format(*repeated_coordinates[:3],
                                                                     repeated_coordinates[4])])

                # contig identifier, start and stop positions and strand
                # 1 for sense, 0 for antisense
                cds_coordinates = ['{0}&{1}-{2}&{3}'.format(*c[:3], c[4])
                                   if c not in invalid_classes else c
                                   for c in cds_coordinates]

                output_lines.append([genome_id]+cds_coordinates)

                if len(output_lines) >= line_limit or (i+1) == len(input_identifiers):
                    output_lines = ['\t'.join(l) for l in output_lines]
                    fo.write_lines(output_lines, output_file, write_mode='a')
                    output_lines = []
        # use distinct sequence identifier if coordinates are not available
        else:
            for i, l in enumerate(csv_reader):
                genome_id = l[0]
                cds_coordinates = []
                for j, c in enumerate(l[1:]):
                    if c not in repeated:
                        cds_coordinates.append(c)
                    else:
                        cds_coordinates.append(classification_labels[-2])
                        # add CDS coordinates to list of PAMA classifications
                        repeated[c].append([genome_id, header[j+1], c])

                output_lines.append([genome_id]+cds_coordinates)

                if len(output_lines) >= line_limit or (i+1) == len(input_identifiers):
                    output_lines = ['\t'.join(l) for l in output_lines]
                    fo.write_lines(output_lines, output_file, write_mode='a')
                    output_lines = []

    # delete intermediate files
    fo.remove_files([intermediate_file, transposed_file])

    return repeated


def create_unclassified_fasta(fasta_file, prot_file, unclassified_protids,
                              protein_hashtable, output_directory, inv_map):
    """Write the coding sequences that were not classified to a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to FASTA file that contains the distinct coding
        sequences identified in the inputs.
    prot_file : str
        Path to FASTA file that contains the distinct translated
        coding sequences identified in the inputs.
    unclassified_protids : list
        List with the sequence identifiers of the representative
        sequences that were not classified.
    protein_hashtable : dict
        Dictionary with SHA-256 hashes for distinct DNA
        sequences extracted from the inputs and lists of
        genome integer identifiers enconded with the
        polyline algorithm as values.
    output_directory : str
        Path to the output directory where the file will be
        created.
    inv_map : dict
        Dictionary with input integer identifiers as keys
        and input string identifiers as values.
    """
    if fasta_file != prot_file:
        unclassified_seqids = []
        prot_distinct_index = fao.index_fasta(prot_file)
        for protid in unclassified_protids:
            prot_seq = str(prot_distinct_index[protid].seq)
            # determine hash
            prot_hash = im.hash_sequence(prot_seq)
            # get all seqids for DNA sequences that code for protein
            seqids = im.polyline_decoding(protein_hashtable[prot_hash])
            # pairs of protein_id, input_id
            seqids = ['{0}-protein{1}'.format(inv_map[seqids[i+1]], seqids[i])
                      for i in range(0, len(seqids), 2)]
            unclassified_seqids.extend(seqids)
    else:
        unclassified_seqids = unclassified_protids

    output_file = fo.join_paths(output_directory, [ct.UNCLASSIFIED_BASENAME])
    dna_index = fao.index_fasta(fasta_file)
    # create FASTA file with unclassified CDSs
    fao.get_sequences_by_id(dna_index, unclassified_seqids, output_file)


def assign_allele_ids(locus_files, ns, repeated):
    """Assign allele identifiers to coding sequences classified as EXC or INF.

    Parameters
    ----------
    classification_files : dict
        Dictionary with the paths to loci FASTA files as keys
        and paths to loci classification files as values.
    ns : bool
        If the schema was downloaded from Chewie-NS. If True,
        adds '*' before allele integer identifiers.
    repeated : list
        List of allele SHA256 hashes for the alleles that matched
        multiple loci. Those alleles are not assigned identifers
        and the classification is changes to PAMA (PAralogous MAtch).

    Returns
    -------
    novel_alleles : dict
        Dictionary with paths to loci FASTA files as keys and
        lists with SHA-256 hashes and allele integer identifiers
        for each novel allele.
    """
    # assign allele identifiers
    novel_alleles = {}
    # import allele calling results and sort to get INF first
    locus_results = fo.pickle_loader(locus_files[1])
    # sort by input order
    sorted_results = sorted(locus_results.items(), key=lambda x: x[0])

    # only keep INF and EXC classifications
    sorted_results = [r for r in sorted_results
                      if r[1][0] in ['EXC', 'INF']]

    # sort to get INF classifications first
    sorted_results = sorted(sorted_results, key=lambda x: x[1][0] == 'INF',
                            reverse=True)

    if len(sorted_results) > 0:
        # import locus records
        records = fao.import_sequences(locus_files[0])
        # determine hash for all locus alleles
        matched_alleles = {im.hash_sequence(v): k.split('_')[-1]
                           for k, v in records.items()}
        # get greatest allele integer identifier
        max_alleleid = max([int(rec.replace('*', '').split('_')[-1])
                            for rec in records])

        for k in sorted_results:
            genome_id = k[0]
            current_results = k[1]
            # get match that was EXC or INF
            current_match = [c for c in current_results[1:]
                             if c[3] in ['EXC', 'INF']][0]
            cds_hash = current_match[2]
            # do not add to schema CDSs that matched several loci
            if cds_hash not in repeated:
                if cds_hash in matched_alleles:
                    locus_results[genome_id].append(matched_alleles[cds_hash])
                else:
                    max_alleleid += 1
                    if ns is True:
                        locus_results[genome_id].append('INF-*{0}'.format(max_alleleid))
                        matched_alleles[cds_hash] = '*' + str(max_alleleid)
                        # add the unique SHA256 value
                        novel_alleles.setdefault(locus_files[0], []).append([cds_hash, '*' + str(max_alleleid)])
                    else:
                        locus_results[genome_id].append('INF-{0}'.format(max_alleleid))
                        matched_alleles[cds_hash] = str(max_alleleid)
                        # add the unique SHA256 value
                        novel_alleles.setdefault(locus_files[0], []).append([cds_hash, str(max_alleleid)])

                    # EXC to INF to enable accurate count of INF classifications
                    # some INF classifications might be converted to NIPH based on similar
                    # matches on the same genome. Matches to the new INF/NIPH will be
                    # classified as EXC and need to be added as new alleles and converted to INF
                    if current_results[0] == 'EXC':
                        locus_results[genome_id][0] = 'INF'
            # classify as PAMA when a CDS matches multiple loci
            else:
                locus_results[genome_id][0] = ct.ALLELECALL_CLASSIFICATIONS[9]

        # save updated results
        fo.pickle_dumper(locus_results, locus_files[1])

    return novel_alleles


def create_novel_fastas(inferred_alleles, inferred_representatives,
                        sequences_file, output_directory):
    """Create FASTA files with the novel alleles for each locus.

    Parameters
    ----------
    inferred_alleles : dict
        Dictionary with paths to loci FASTA files as keys and
        lists with SHA-256 hashes, allele integer identifiers and
        sequence identifiers for each novel allele.
    inferred_representatives : dict
        Dictionary with loci identifiers as keys and lists with
        sequence identifiers, SHA-256 hashes and allele integer
        identifiers for each novel representative allele.
    sequences_file : str
        Path to FASTA file that contains the distinct coding
        sequences identified in the inputs.
    output_directory : str
        Path to the directory where the FASTA files that contain
        the novel alleles will be saved to.

    Returns
    -------
    total_inferred : int
        Total number of inferred alleles added to the schema.
    total_representatives : int
        Total number of representative alleles added to the
        schema.
    updated_novel : dict
        Dictionary with loci identifiers as keys and paths to
        the FASTA files that contain the novel alleles as values.
    """
    # create index for Fasta file with distinct CDSs
    sequence_index = fao.index_fasta(sequences_file)

    # count number of novel and representative alleles added to schema
    total_inferred = 0
    total_representative = 0
    updated_novel = {}
    for locus, alleles in inferred_alleles.items():
        locus_id = fo.get_locus_id(locus)
        if locus_id is None:
            locus_id = fo.file_basename(locus, False)

        updated_novel[locus] = []
        # get novel alleles through indexed Fasta file
        novel_alleles = ['>{0}_{1}\n{2}'.format(locus_id, a[1],
                                                str(sequence_index.get(a[2]).seq))
                         for a in alleles]
        # create Fasta file with novel alleles
        novel_file = fo.join_paths(output_directory, ['{0}.fasta'.format(locus_id)])
        fo.write_lines(novel_alleles, novel_file)
        updated_novel[locus].append(novel_file)
        total_inferred += len(novel_alleles)

        # add representatives
        novel_representatives = inferred_representatives.get(locus_id, None)
        if novel_representatives is not None:
            reps_sequences = ['>{0}_{1}\n{2}'.format(locus_id, a[2], str(sequence_index.get(a[0]).seq))
                              for a in novel_representatives]
            # create Fasta file with novel representative alleles
            novel_rep_file = fo.join_paths(output_directory, ['short', '{0}_short.fasta'.format(locus_id)])
            fo.write_lines(reps_sequences, novel_rep_file)
            updated_novel[locus].append(novel_rep_file)
            total_representative += len(reps_sequences)

    return [total_inferred, total_representative, updated_novel]


def add_inferred_alleles(inferred_alleles):
    """Add inferred alleles to a schema.

    Prameters
    ---------
    inferred_alleles . dict
        Dictionary with loci identifiers as keys and paths
        to the FASTA files that contain the novel alleles as values.

    Returns
    -------
    added : list
        List that contains the number of novel alleles added for each
        locus (same order as `inferred_alleles` keys)
    """
    added = []
    for locus, files in inferred_alleles.items():
        locus_id = fo.get_locus_id(locus)
        if locus_id is None:
            locus_id = fo.file_basename(locus, False)

        # get novel alleles
        novel_alleles = fo.read_lines(files[0])
        # append novel alleles to locus FASTA file
        fo.write_lines(novel_alleles, locus, write_mode='a')
        added.append(len(novel_alleles)//2)

        if len(files) > 1:
            # add representatives
            novel_representatives = fo.read_lines(files[1])
            # append novel alleles to file in 'short' directory
            locus_short_path = fo.join_paths(os.path.dirname(locus),
                                             ['short', locus_id+'_short.fasta'])
            fo.write_lines(novel_representatives, locus_short_path, write_mode='a')

    return added


def select_highest_scores(blast_outfile):
    """Select the highest-scoring match for each BLAST target.

    Parameters
    ----------
    blast_outfile : str
        Path to the TSV file created by BLAST.

    Returns
    -------
    best_matches : list
        List with the highest-scoring match/line for each
        distinct target.
    """
    blast_results = fo.read_tabular(blast_outfile)
    # sort results based on decreasing raw score
    blast_results = im.sort_iterable(blast_results,
                                     lambda x: int(x[5]), reverse=True)

    # select matches with highest score for each target
    best_matches = []
    for r in blast_results:
        # only get the best raw score for each target
        if r[4] not in best_matches:
            best_matches.append(r)

    return best_matches


def process_blast_results(blast_results, bsr_threshold, query_scores,
                          inputids_mapping):
    """Process BLAST results to get data relevant for classification.

    Parameters
    ----------
    blast_results : list
        List with one sublist per BLAST match (must have one
        sublist per target with the highest-scoring match).
    bsr_threshold : float
        BLAST Score Ratio (BSR) value to select matches. Matches
        will be kept if the computed BSR is equal or greater than
        this value.
    query_scores :  dict
        Dictionary with loci representative sequence identifiers
        as values and a tuple with sequence length and the raw
        score for the self-alignment as value.
    inputids_mapping : dict
        Maping between short sequence identifiers and original
        sequence identifiers.

    Returns
    -------
    match_info : dict
        Dictionary with the distinct target sequence identifiers
        as keys and a tuple with the BSR value, target sequence
        length, query sequence length and query sequence identifier
        for the highest-scoring match for each target as values.
    """
    # replace query and target identifiers if they were simplified
    # to avoid BLAST warnings/errors
    if inputids_mapping is not None:
        blast_results = [im.replace_list_values(r, inputids_mapping)
                         for r in blast_results]

    # determine BSR values
    match_info = {}
    for r in blast_results:
        query_id = r[0]
        target_id = r[4]
        raw_score = float(r[6])
        bsr = cf.compute_bsr(raw_score, query_scores[query_id][1])
        # only keep matches above BSR threshold
        if bsr >= bsr_threshold:
            # BLAST has 1-based positions
            qstart = (int(r[1])-1)*3  # subtract 1 to exclude start position
            qend = (int(r[2])*3)+3  # add 3 to count stop codon
            target_length = (int(r[5])*3)+3
            query_length = query_scores[query_id][0]

            match_info[target_id] = (bsr, qstart, qend,
                                     target_length, query_length, query_id)

    return match_info


def expand_matches(match_info, pfasta_index, dfasta_index, dhashtable,
                   phashtable, inv_map):
    """Expand distinct matches to create matches for all matched inputs.

    Parameters
    ----------
    match_info : dict
        Dictionary with the distinct target sequence identifiers
        as keys and a tuple with the BSR value, target sequence
        length, query sequence length and query sequence identifier
        for the highest-scoring match for each target as values.
    pfasta_index : Bio.File._IndexedSeqFileDict
        Fasta file index created with BioPython. Index for the
        distinct protein sequences.
    dfasta_index : Bio.File._IndexedSeqFileDict
        Fasta file index created with BioPython. Index for the
        distinct DNA coding sequences.
    dhashtable : dict
        Dictionary with SHA-256 hashes for distinct DNA
        sequences extracted from the inputs and lists of
        genome integer identifiers enconded with the
        polyline algorithm as values.
    phashtable : dict
        Dictionary with SHA-256 hashes for distinct protein
        sequences extracted from the inputs and lists of
        sequence identifiers enconded with the polyline
        algorithm as values.
    inv_map : dict
        Dictionary with input integer identifiers as keys
        and input string identifiers as values.

    Returns
    -------
    input_matches : dict
        Dictionary with input integer identifiers as keys
        and tuples with information about matches identified
        in the inputs as values.
    """
    input_matches = {}
    for target_id in match_info:
        target_protein = str(pfasta_index.get(target_id).seq)
        target_phash = im.hash_sequence(target_protein)
        target_integers = im.polyline_decoding(phashtable[target_phash])
        target_seqids = ['{0}-protein{1}'.format(inv_map[target_integers[i+1]], target_integers[i])
                         for i in range(0, len(target_integers), 2)]
        for seqid in target_seqids:
            target_cds = str(dfasta_index.get(seqid).seq)
            target_dhash = im.hash_sequence(target_cds)
            # get ids for all genomes with same CDS as representative
            target_inputs = im.polyline_decoding(dhashtable[target_dhash])[1:]
            for i in target_inputs:
                input_matches.setdefault(i, []).append((target_id, target_phash,
                                                        target_dhash, *match_info[target_id]))

    return input_matches


def identify_paralogous(repeated, output_directory):
    """Identifiy groups of paralogous loci in the schema.

    Parameters
    ----------
    repeated : dict
        Dictionary with hashes for the CDSs that matched multiple loci
        as keys and a list with information about the genome of origin
        and matched loci as values.
    output_directory : str
        Path to the output directory where the file with
        the list of paralogus loci will be created.

    Returns
    -------
    The total number of paralogous loci detected.
    """
    paralogous_data = {}
    paralogous_counts = {}
    for value in repeated.values():
        for p in value:
            paralogous_data.setdefault(p[2], {})
            paralogous_data[p[2]].setdefault(p[0], []).append(p[1])
            paralogous_counts.setdefault(p[1], []).append(1)

    # write paralogous loci counts
    paralogous_counts_lines = [ct.PARALOGOUS_COUNTS_HEADER]
    paralogous_counts_lines.extend(['{0}\t{1}'.format(k, sum(v))
                                    for k, v in paralogous_counts.items()])
    paralogous_counts_outfile = fo.join_paths(output_directory,
                                              [ct.PARALOGOUS_COUNTS_BASENAME])
    fo.write_lines(paralogous_counts_lines, paralogous_counts_outfile)

    # write groups of paralogous loci per input
    paralogous_lines = [ct.PARALOGOUS_LIST_HEADER]
    for cds, g in paralogous_data.items():
        paralogous_lines.extend(['{0}\t{1}\t{2}'.format(k, '|'.join(v), cds)
                                 for k, v in g.items()])

    paralogous_loci_outfile = fo.join_paths(output_directory,
                                            [ct.PARALOGOUS_LOCI_BASENAME])
    fo.write_lines(paralogous_lines, paralogous_loci_outfile)

    return len(paralogous_counts)


def classify_inexact_matches(locus, genomes_matches, inv_map,
                             locus_results_file, locus_mode, temp_directory,
                             size_threshold, blast_score_ratio, output_directory,
                             cds_input):
    """Classify inexact matches found for a locus.

    Parameters
    ----------
    locus : str
        Locus identifier.
    genomes_matches : str
        Path to file with data about matches found in the inputs.
    inv_map : dict
        Dictionary with input integer identifeirs as keys and
        input string identifiers as values.
    locus_results_file : str
        Path to file with classification results for the locus.
    locus_mode : list
        List where wthe first element is the locus allele size mode
        and the second element is a list with the length values for
        all alleles.
    temp_directory : str
        Path to the directory where temporary files will be stored.
    size_threshold : float or None
        Sequence size variation threshold.
    blast_score_ratio : float
        BLAST Score Ratio value.
    output_directory : str
        Path to the directory where the files with information for
        each locus will be saved to.
    cds_input : bool
        True if the input files contained CDSs, False if they contained
        contigs.

    Returns
    -------
    locus_info_file : str
        Path to pickle file with the data about the classification
        of inexact matches (contains dictionary with locus identifier
        as key and a list with the path to the pickle file with the
        locus classifications, locus allele size mode, sequence
        identifiers of the distinct sequences that were classified
        and a list with data about representative candidates as value).
    """
    # import classifications
    locus_results = fo.pickle_loader(locus_results_file)

    # import matches
    genomes_matches = fo.pickle_loader(genomes_matches)

    # initialize lists to store hashes of CDSs that have been classified
    seen_dna = {}
    seen_prot = []
    # initialize list to store sequence identifiers that have been classified
    excluded = []
    representative_candidates = []
    for genome, matches in genomes_matches.items():
        current_g = inv_map[genome]

        for match in matches:
            # get sequence identifier for the representative CDS
            target_seqid = match[0]
            # get allele identifier for the schema representative
            rep_alleleid = match[8]
            # determine if schema representative has allele id
            try:
                # split to get allele identifier
                # need to replace '*' for novel alleles added to schemas from Chewie-NS
                int(rep_alleleid.replace('*', '').split('_')[-1])
                rep_alleleid = rep_alleleid.split('_')[-1]
            except Exception as e:
                pass

            # get hash of the CDS DNA sequence
            target_dna_hash = match[2]
            # get hash of the translated CDS sequence
            target_prot_hash = match[1]

            # get the BSR value
            bsr = match[3]

            # CDS DNA sequence was identified in one of the previous inputs
            # This will change classification to NIPH if the input
            # already had a classification for the current locus
            if target_dna_hash in seen_dna:
                locus_results = update_classification(genome, locus_results,
                                                      (seen_dna[target_dna_hash], target_seqid,
                                                       target_dna_hash, 'EXC', 1.0))
                continue

            # translated CDS matches other translated CDS that was classified
            if target_prot_hash in seen_prot:
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                       target_dna_hash, 'INF', 1.0))
                # add DNA hash to classify the next match as EXC
                seen_dna[target_dna_hash] = target_seqid
                continue

            if cds_input is False:
                # there is no DNA or Protein exact match, perform full evaluation
                # open pickle for genome and get coordinates
                genome_cds_file = fo.join_paths(temp_directory, ['2_cds_extraction', current_g+'.cds_hash'])
                genome_cds_coordinates = fo.pickle_loader(genome_cds_file)
                # classifications based on position on contig (PLOT3, PLOT5 and LOTSC)
                # get CDS start and stop positions
                genome_coordinates = genome_cds_coordinates[0][target_dna_hash][0]
                contig_leftmost_pos = int(genome_coordinates[1])
                contig_rightmost_pos = int(genome_coordinates[2])
                # get contig length
                contig_length = genome_cds_coordinates[1][genome_coordinates[0]]
                # get representative length
                representative_length = match[7]
                # get target left and right positions that aligned
                representative_leftmost_pos = match[4]
                representative_rightmost_pos = match[5]
                # determine if it is PLOT3, PLOT5 or LOTSC
                relative_pos = contig_position_classification(representative_length,
                                                              representative_leftmost_pos,
                                                              representative_rightmost_pos,
                                                              contig_length,
                                                              contig_leftmost_pos,
                                                              contig_rightmost_pos)

                if relative_pos is not None:
                    locus_results = update_classification(genome, locus_results,
                                                          (rep_alleleid, target_seqid,
                                                           target_dna_hash, relative_pos, bsr))
                    # need to exclude so that it does not duplicate ASM/ALM classifications later
                    excluded.append(target_seqid)
                    continue

            target_dna_len = match[6]
            # check if ASM or ALM
            relative_size = allele_size_classification(target_dna_len, locus_mode, size_threshold)
            if relative_size is not None:
                locus_results = update_classification(genome, locus_results,
                                                      (rep_alleleid, target_seqid,
                                                       target_dna_hash, relative_size, bsr))
                # need to exclude so that it does not duplicate PLOT3/5 classifications later
                excluded.append(target_seqid)
                continue

            # add INF
            # this will turn into NIPH if there are multiple hits for the same input
            locus_results = update_classification(genome, locus_results,
                                                  (rep_alleleid, target_seqid,
                                                   target_dna_hash, 'INF', bsr))

            seen_dna[target_dna_hash] = target_seqid
            excluded.append(target_seqid)
            seen_prot.append(target_prot_hash)

        # update locus mode value if classification for genome is INF
        if genome in locus_results and locus_results[genome][0] == 'INF':
            # append length of inferred allele to list with allele sizes
            locus_mode[1].append(target_dna_len)
            # compute mode
            locus_mode[0] = sm.determine_mode(locus_mode[1])[0]
            # only add as representative candidate if classification is not NIPH
            inf_bsr = locus_results[genome][1][4]
            if inf_bsr >= blast_score_ratio and inf_bsr < blast_score_ratio+0.1:
                representative_candidates.append((genome, target_seqid,
                                                  match[8], target_dna_hash))

    # save updated results
    fo.pickle_dumper(locus_results, locus_results_file)

    # save info about updated mode, excluded ids and representative candidates
    locus_info = {locus: [locus_results_file, locus_mode,
                          excluded, representative_candidates]}
    locus_info_file = fo.join_paths(output_directory, ['{0}_classification_info'.format(locus)])
    fo.pickle_dumper(locus_info, locus_info_file)

    return locus_info_file


def create_missing_fasta(class_files, fasta_file, input_map, dna_hashtable,
                         output_directory, coordinates_files, classification_labels):
    """Create Fasta file with sequences for missing data classes.

    Parameters
    ----------
    class_files : dict
        Dictionary with paths to loci files as keys and paths to
        pickled files with classification results as values.
    fasta_file : str
        Path to Fasta file with the distinct CDS extracted from
        the input genomes.
    input_map : dict
        Dictionary with the mapping between the input integer
        identifiers and input string identifiers.
    dna_hashtable : dict
        Dictionary with hashes of the distinct CDS extracted from
        input genomes as keys and lists containing the integer
        identifiers for te inputs that contained the CDS encoded
        with the polyline algorithm.
    output_directory : str
        Path to the output directory where the Fasta file will
        be saved to.
    coordinates_files : dict
        Dictionary with the mapping between input string identifiers
        and paths to pickled files that contain a dictionary with the
        coordinates of the CDS identified in each input.
    classification_labels : list
        List with the class labels attributed by chewBBACA.
    """
    invalid_cases = classification_labels[2:-1]

    # get information about missing cases for each input genome
    missing_cases = {}
    # get hash if coordinates are available, seqid otherwise
    id_index = 2 if coordinates_files is not None else 1
    for locus, file in class_files.items():
        locus_id = fo.get_locus_id(locus)
        if locus_id is None:
            locus_id = fo.file_basename(locus, False)
        locus_classifications = fo.pickle_loader(file)
        # get data for genomes that do not have EXC or INF classifications
        # it will not get invalid classes if a genome is classified as EXC|INF
        for gid, v in locus_classifications.items():
            if v[0] in invalid_cases:
                genome_info = [locus_id, v[0], [[e[id_index], e[3]] for e in v[1:]]]
                missing_cases.setdefault(input_map[gid], []).append(genome_info)

    # add integer index to FASTA header and first TSV column
    index = 1
    tsv_lines = ['Index\tGenome\tLocus\tLocus_classification\tCDS\tCDS_classification']
    # fetch genome coordinates and create FASTA headers
    # input files contained genome assemblies
    if coordinates_files is not None:
        for k, v in missing_cases.items():
            genome_coordinates = fo.pickle_loader(coordinates_files[k])[0]
            # genomes may have duplicated CDSs
            # store hash and increment i to get correct positions
            hashes = {}
            for c in v:
                locus_id = c[0]
                classification = c[1]
                for h in c[2]:
                    current_hash = h[0]
                    coordinates = genome_coordinates[current_hash]

                    if current_hash not in hashes:
                        hashes[current_hash] = {locus_id: 0}
                    else:
                        # multiple matches to the same locus
                        if locus_id in hashes[current_hash]:
                            hashes[current_hash][locus_id] += 1
                        # multiple matches to multiple loci
                        else:
                            hashes[current_hash][locus_id] = 0

                    current_index = hashes[current_hash][locus_id]
                    protid = coordinates[current_index][3]
                    h.append('{0}|{1}|{2}&{3}|{1}-protein{4}&{5}'.format(index,
                                                                         k,
                                                                         locus_id,
                                                                         classification,
                                                                         protid,
                                                                         h[1]))
                    tsv_lines.append('{0}\t{1}\t{2}\t{3}\t{1}-protein{4}\t{5}'.format(index,
                                                                         k,
                                                                         locus_id,
                                                                         classification,
                                                                         protid,
                                                                         h[1]))
                    index += 1
    # input files contained CDSs and there are no genome coordinates
    else:
        for k, v in missing_cases.items():
            for c in v:
                locus_id = c[0]
                classification = c[1]
                for h in c[2]:
                    seqid = h[0]
                    seqid_classification = h[1]
                    h.append('{0}|{1}|{2}&{3}|{4}&{5}'.format(index, k,
                                                              locus_id,
                                                              classification,
                                                              seqid,
                                                              seqid_classification))
                    tsv_lines.append('{0}\t{1}\t{2}\t{3}\t{1}-protein{4}\t{5}'.format(index, k,
                                                              locus_id,
                                                              classification,
                                                              seqid,
                                                              seqid_classification))
                    index += 1

    # get FASTA sequences
    missing_records = []
    dna_index = fao.index_fasta(fasta_file)
    if coordinates_files is not None:
        for value in missing_cases.values():
            current_records = []
            for c in value:
                for h in c[2]:
                    hash_entry = im.polyline_decoding(dna_hashtable[h[0]])
                    seqid = '{0}-protein{1}'.format(input_map[hash_entry[1]], hash_entry[0])
                    new_rec = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
                                                   [h[2],
                                                    str(dna_index[seqid].seq)])
                    current_records.append(new_rec)

            missing_records.extend(current_records)
    else:
        for value in missing_cases.values():
            current_records = []
            for c in value:
                for h in c[2]:
                    new_rec = fao.fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
                                                   [h[2],
                                                    str(dna_index[h[0]].seq)])
                    current_records.append(new_rec)

            missing_records.extend(current_records)

    # write FASTA file
    output_file = fo.join_paths(output_directory, ['missing_classes.fasta'])
    fo.write_lines(missing_records, output_file)

    # write TSV file
    output_file = fo.join_paths(output_directory, ['missing_classes.tsv'])
    fo.write_lines(tsv_lines, output_file)


def select_representatives(representative_candidates, locus, fasta_file,
                           iteration, output_directory, blastp_path,
                           blast_db, blast_score_ratio, threads):
    """Select new representative alleles for a locus.

    Parameters
    ----------
    representative_candidates : dict
        Dictionary with sequence identifiers as keys and sequence
        hashes as values.
    locus : str
        Locus identifier.
    fasta_file : path
        Path to Fasta file that contains the translated sequences
        of the representative candidates.
    iteration : int
        Iteration number to add to generated files.
    output_directory : str
        Path to the output directory.
    blastp_path : str
        Path to the BLASTp executable.
    blast_db : str
        Path to the BLAST database.
    blast_score_ratio : float
        BLAST Score Ratio value.
    threads : int
        Number of threads passed to BLAST.

    Returns
    -------
    locus : str
        Locus identifier.
    selected : list
        List that contains one tuple per selected candidate (tuples
        contain the sequence identifier and the sequence hash for
        each new representative).
    """
    # create file with candidate ids
    ids_file = fo.join_paths(output_directory,
                             ['{0}_candidates_ids_{1}.fasta'.format(locus, iteration)])
    fo.write_lines(list(representative_candidates.keys()), ids_file)

    # BLASTp to compare all candidates
    blast_output = fo.join_paths(output_directory,
                                 ['{0}_candidates_{1}_blastout.tsv'.format(locus, iteration)])
    # pass number of max targets per query to reduce execution time
    blastp_stderr = bw.run_blast(blastp_path, blast_db, fasta_file,
                                 blast_output, threads=threads,
                                 ids_file=ids_file, max_targets=30)

    blast_results = fo.read_tabular(blast_output)
    # get self scores
    candidates_self_scores = {line[0]: ((int(line[3])*3)+3, float(line[6]))
                              for line in blast_results if line[0] == line[4]}
    # select results between different candidates
    blast_results = [line for line in blast_results if line[0] != line[4]]

    # compute bsr
    for line in blast_results:
        line.append(cf.compute_bsr(candidates_self_scores[line[4]][1], candidates_self_scores[line[0]][1]))
    # sort by sequence length to process longest candidates first
    blast_results = sorted(blast_results, key=lambda x: int(x[3]))

    excluded_candidates = []
    for r in blast_results:
        if r[7] >= blast_score_ratio+0.1:
            if r[4] not in excluded_candidates:
                excluded_candidates.append(r[0])

    selected_candidates = list(set([line[0]
                                    for line in blast_results
                                    if line[0] not in excluded_candidates]))

    selected = [(line, representative_candidates[line])
                for line in selected_candidates
                if line not in excluded_candidates]

    return [locus, selected]


def allele_calling(fasta_files, schema_directory, temp_directory,
                   loci_modes, loci_files, config, pre_computed_dir):
    """Perform allele calling for a set of inputs.

    Parameters
    ----------
    fasta_files : list
        List with the full paths to the input FASTA files, one per
        strain.
    schema_directory : str
        Path to the schema directory.
    temp_directory : str
        Path to the temporary directory in which the intermediate files
        are created.
    loci_modes : dict
        Dictionary with the allele size mode for each locus.
    loci_files : dict
        Dictionary with the paths to the loci FASTA files as keys and
        the loci integer identifiers as values.
    config : dict
        Dictionary with the config values that will be used to perform
        allele calling.
    pre_computed_dir : str
        Path to the the directory that contains the files with the
        pre-computed hash tables.

    Returns
    -------
    template_dict : dict
        Dictionary used to store data and paths to files with results
        (paths to files with the classification results, dictionary with
        the mapping between input identifiers and unique integer identifiers,
        paths to files with the CDSs coordinates for each input genome,
        path to a FASTA file with all distinct DNA sequences, path to a
        FASTA file with all distinct protein sequences, a hash table with
        the mapping between the SHA256 hash for each distinct CDS extracted
        from the input genomes and the list of input genomes that contain
        each CDS, a hash table with the mapping between the SHA256 hash
        for each translated CDS and the distinct DNA CDSs that code for each
        protein, list with the paths to input files for which it was not
        possible to predict CDSs with prodigal, path to file with invalid
        CDSs, list with the sequence identifiers for the sequences that
        were not classified, dictionary with the mapping between the
        sequence identifiers of the representative alleles and their
        self-alignment BLASTp raw score, dictionary with information
        about new representatives for each locus).
    """
    # get dictionary template to store variables to return
    template_dict = ct.ALLELECALL_DICT

    # map full paths to basename
    inputs_basenames = im.mapping_function(fasta_files,
                                           fo.file_basename, [False])

    # map input identifiers to integers
    # use the mapped integers to refer to each input
    # this reduces memory usage compared to using string identifiers
    basename_map = im.integer_mapping(inputs_basenames.values())
    basename_inverse_map = im.invert_dictionary(basename_map)
    template_dict['basename_map'] = basename_inverse_map

    # detect if some inputs share the same unique prefix
    if len(basename_inverse_map) < len(fasta_files):
        basename_counts = [[basename, list(inputs_basenames.values()).count(basename)]
                           for basename in basename_map]
        repeated_basenames = ['{0}: {1}'.format(*l)
                              for l in basename_counts if l[1] > 1]
        # only delete temp directory created for each run
        # do not delete output directory because it might include other files
        fo.delete_directory(temp_directory)
        sys.exit('\nSome input files share the same filename prefix '
                 '(substring before the first "." in the filename). '
                 'Please make sure that every input file has a unique '
                 'filename prefix.\n{0}'.format('\n'.join(repeated_basenames)))

    # inputs are genome assemblies
    if config['CDS input'] is False:
        # create directory to store files with Prodigal results
        prodigal_path = fo.join_paths(temp_directory, ['1_cds_prediction'])
        fo.create_directory(prodigal_path)

        # run Prodigal to determine CDSs for all input genomes
        print('\n== CDS prediction ==\n')

        # gene prediction step
        print('Predicting CDS for {0} inputs...'.format(len(fasta_files)))
        failed = cf.predict_genes(fasta_files,
                                  config['Prodigal training file'],
                                  config['Translation table'],
                                  config['Prodigal mode'],
                                  config['CPU cores'],
                                  prodigal_path)

        if len(failed) > 0:
            print('\nFailed to predict CDS for {0} inputs'
                  '.'.format(len(failed)))
            print('Make sure that Prodigal runs in meta mode (--pm meta) '
                  'if any input file has less than 100kbp.')

            # remove failed genomes from paths
            fasta_files = im.filter_list(fasta_files, failed)

        if len(fasta_files) == 0:
            sys.exit('\nCould not predict CDS for any '
                     'of the input files.\nPlease provide input files '
                     'in the accepted FASTA format.')

        # CDS extraction step
        print('\n\n== CDS extraction ==\n')
        # create output directory
        cds_extraction_path = fo.join_paths(temp_directory,
                                            ['2_cds_extraction'])
        fo.create_directory(cds_extraction_path)
        print('Extracting predicted CDS for {0} inputs...'.format(len(fasta_files)))
        eg_results = cf.extract_genes(fasta_files,
                                      prodigal_path,
                                      config['CPU cores'],
                                      cds_extraction_path)
        cds_files, total_extracted, cds_coordinates = eg_results

        print('\nExtracted a total of {0} CDS from {1} '
              'inputs.'.format(total_extracted, len(fasta_files)))
    # inputs are Fasta files with the predicted CDSs
    else:
        failed = []
        cds_files = fasta_files
        cds_coordinates = {}

    if len(failed) > 0:
        template_dict['invalid_inputs'] = failed
    if len(cds_coordinates) > 0:
        template_dict['cds_coordinates'] = cds_coordinates

    # create directory to store files from pre-process steps
    preprocess_dir = fo.join_paths(temp_directory, ['3_cds_preprocess'])
    fo.create_directory(preprocess_dir)

    # DNA sequences deduplication step
    # keep hash of unique sequences and a list with the integer
    # identifiers of genomes that have those sequences
    # lists of integers are encoded with polyline algorithm
    print('\n== CDS deduplication ==')
    # create directory to store files from DNA deduplication
    dna_dedup_dir = fo.join_paths(preprocess_dir, ['cds_deduplication'])
    fo.create_directory(dna_dedup_dir)
    distinct_dna_template = 'distinct_cds_{0}'
    print('\nIdentifying distinct CDS...', end='')
    dna_dedup_results = cf.exclude_duplicates(cds_files, dna_dedup_dir,
                                              config['CPU cores'],
                                              distinct_dna_template,
                                              [basename_map, basename_inverse_map],
                                              False, False)

    dna_distinct_htable, distinct_file, repeated = dna_dedup_results
    print('identified {0} distinct CDS.'.format(len(dna_distinct_htable)))
    template_dict['dna_fasta'] = distinct_file
    template_dict['dna_hashtable'] = dna_distinct_htable

    # get mapping between locus file path and locus identifier
    loci_basenames = im.mapping_function(loci_files, fo.file_basename, [False])

    # create files with empty results data structure
    classification_dir = fo.join_paths(temp_directory, ['classification_files'])
    fo.create_directory(classification_dir)
    empty_results = {}
    inputs = [[loci_basenames[file], classification_dir, empty_results]
              for file in loci_files]
    classification_files = {file: create_classification_file(*inputs[i])
                            for i, file in enumerate(loci_files)}

    print('\n== CDS exact matches ==')
    matched_seqids = []
    dna_exact_hits = 0
    dna_matches_ids = 0
    print('\nSearching for DNA exact matches...', end='')
    # list files in pre-computed dir
    dna_tables = fo.listdir_fullpath(pre_computed_dir, 'DNAtable')
    for file in dna_tables:
        dna_matches = dna_exact_matches(file,
                                        dna_distinct_htable,
                                        im.invert_dictionary(loci_files),
                                        classification_files,
                                        basename_inverse_map)
        # might have duplicate IDs
        matched_seqids.extend(dna_matches[0])
        dna_exact_hits += dna_matches[1]
        dna_matches_ids += dna_matches[2]

    print('found {0} exact matches (matching {1} distinct '
          'alleles).'.format(dna_exact_hits, dna_matches_ids))

    # save seqids that matched
    dna_exact_outfile = fo.join_paths(preprocess_dir, ['cds_exact_matches.txt'])
    fo.write_lines(matched_seqids, dna_exact_outfile)

    # get sequence identifiers for unclassified sequences
    # reading to get lines with '>' is faster that reading with BioPython
    # and filtering based on sequence identifiers
    matched_lines = fo.matching_lines(distinct_file, '>')
    matched_lines = [line.strip()[1:] for line in matched_lines]
    selected_ids = im.filter_list(matched_lines, matched_seqids)

    print('Unclassified CDS: {0}'.format(len(selected_ids)))

    # user only wants to determine exact matches
    if config['Mode'] == 1:
        template_dict['classification_files'] = classification_files
        template_dict['protein_fasta'] = distinct_file
        template_dict['unclassified_ids'] = selected_ids

        return template_dict

    # create Fasta file without distinct sequences that were exact matches
    dna_index = fao.index_fasta(distinct_file)

    # translate DNA sequences and identify duplicates
    print('\n== CDS translation ==\n')

    # this step excludes small sequences
    # create directory to store translation results
    cds_translation_dir = fo.join_paths(preprocess_dir, ['cds_translation'])
    fo.create_directory(cds_translation_dir)
    print('Translating {0} CDS...'.format(len(selected_ids)))
    ts_results = cf.translate_sequences(selected_ids, distinct_file,
                                        cds_translation_dir,
                                        config['Translation table'],
                                        config['Minimum sequence length'],
                                        config['CPU cores'])

    protein_file, ut_seqids, ut_lines = ts_results

    print('\nIdentified {0} CDS that could not be translated.'.format(len(ut_seqids)))

    # write info about invalid alleles to file
    invalid_alleles_file = fo.join_paths(temp_directory,
                                         ['invalid_cds.txt'])
    invalid_alleles = im.join_list(im.sort_iterable(ut_lines), '\n')
    fo.write_to_file(invalid_alleles, invalid_alleles_file, 'w', '\n')
    print('Information about untranslatable and small sequences '
          'stored in {0}'.format(invalid_alleles_file))
    template_dict['invalid_alleles'] = invalid_alleles_file
    print('Unclassified CDS: {0}'.format(len(selected_ids)-len(ut_seqids)))

    # protein sequences deduplication step
    print('\n== Protein deduplication ==')
    # create directory to store files from protein deduplication
    protein_dedup_dir = fo.join_paths(preprocess_dir, ['protein_deduplication'])
    fo.create_directory(protein_dedup_dir)
    distinct_prot_template = 'distinct_proteins_{0}'
    print('\nIdentifying distinct proteins...', end='')
    ds_results = cf.exclude_duplicates([protein_file], protein_dedup_dir, 1,
                                       distinct_prot_template, [basename_map, basename_inverse_map], True, False)
    distinct_pseqids = ds_results[0]
    print('identified {0} distinct proteins.'.format(len(distinct_pseqids)))
    template_dict['protein_hashtable'] = distinct_pseqids

    # identify exact matches at protein level
    # exact matches are novel alleles that can be added to the schema
    print('\n== Protein exact matches ==')
    exc_cds = 0
    exc_prot = 0
    exc_distinct_prot = 0
    exact_phashes = []
    print('\nSearching for Protein exact matches...', end='')
    # list files in pre-computed dir
    protein_tables = fo.listdir_fullpath(pre_computed_dir, 'PROTEINtable')
    for file in protein_tables:
        protein_matches = protein_exact_matches(file,
                                                distinct_pseqids,
                                                im.invert_dictionary(loci_files),
                                                classification_files,
                                                basename_inverse_map,
                                                dna_index,
                                                dna_distinct_htable)

        exact_phashes.extend(protein_matches[0])
        exc_prot += protein_matches[1]
        exc_cds += protein_matches[2]
        exc_distinct_prot += protein_matches[3]
        # update locus mode
        if len(protein_matches[4]) > 0:
            for r in protein_matches[4]:
                loci_modes[r[0]][1].extend(r[1])
                loci_modes[r[0]][0] = sm.determine_mode(loci_modes[r[0]][1])[0]

    print('found {0} exact matches ({1} distinct CDS, {2} total CDS).'
          ''.format(exc_distinct_prot, exc_prot, exc_cds))

    # save seqids that matched
    protein_exact_outfile = fo.join_paths(preprocess_dir, ['protein_exact_matches.txt'])
    fo.write_lines(exact_phashes, protein_exact_outfile)

    # create new Fasta file without the Protein sequences that were exact matches
    unique_pfasta = fo.join_paths(preprocess_dir, ['protein_distinct.fasta'])
    # create protein file index
    protein_index = fao.index_fasta(ds_results[1])
    # the list of "exact_phases" corresponds to the seqids for the DNA sequences
    # this means that it can have more elements that the number of protein exact matches
    # because different alleles might code for same protein
    matched_lines = fo.matching_lines(ds_results[1], '>')
    matched_lines = [line.strip()[1:] for line in matched_lines]
    selected_ids = im.filter_list(matched_lines, exact_phashes)
    total_selected = fao.get_sequences_by_id(protein_index, selected_ids, unique_pfasta)

    print('Unclassified proteins: {0}'.format(total_selected))

    if config['Mode'] == 2:
        template_dict['classification_files'] = classification_files
        template_dict['protein_fasta'] = unique_pfasta
        template_dict['unclassified_ids'] = selected_ids

        return template_dict

    # translate schema representatives
    print('\n== Clustering ==')
    print('\nTranslating schema\'s representative alleles...', end='')
    rep_dir = fo.join_paths(schema_directory, ['short'])
    rep_list = fo.listdir_fullpath(rep_dir, '.fasta')
    # filter to get only files in list of loci
    rep_basenames = {file: fo.file_basename(file, False).replace('_short', '') for file in rep_list}
    rep_list = [k for k, v in rep_basenames.items() if v in loci_basenames.values()]

    reps_protein_dir = fo.join_paths(temp_directory, ['4_translated_representatives'])
    fo.create_directory(reps_protein_dir)
    protein_files = mo.parallelize_function(fao.translate_fasta, rep_list,
                                            [reps_protein_dir,
                                             config['Translation table']],
                                            config['CPU cores'],
                                            False)
    protein_repfiles = [r[1] for r in protein_files]
    print('done.')

    # cluster protein sequences
    proteins = fao.import_sequences(unique_pfasta)

    # create directory to store clustering data
    clustering_dir = fo.join_paths(temp_directory, ['5_clustering'])
    fo.create_directory(clustering_dir)

    # define BLASTp and makeblastdb paths
    blastp_path = fo.join_paths(config['BLAST path'], [ct.BLASTP_ALIAS])
    makeblastdb_path = fo.join_paths(config['BLAST path'], [ct.MAKEBLASTDB_ALIAS])

    # concatenate all schema representative
    concat_reps = fo.join_paths(reps_protein_dir, ['concat_reps.fasta'])
    fo.concatenate_files(protein_repfiles, concat_reps)

    # create dict with mapping between allele header and locus id
    rep_recs = fao.sequence_generator(concat_reps)
    rep_map = {}
    for rec in rep_recs:
        rep_map[rec.id] = '_'.join((rec.id).split('_')[:-1])

    # determine self-score for representatives if file is missing
    self_score_file = fo.join_paths(schema_directory, ['short', 'self_scores'])
    if os.path.isfile(self_score_file) is False:
        print('Determining BLASTp raw score for each representative...', end='')
        self_score_dir = fo.join_paths(reps_protein_dir, ['self_scores'])
        fo.create_directory(self_score_dir)
        self_scores = cf.determine_self_scores(concat_reps, self_score_dir,
                                               makeblastdb_path, blastp_path,
                                               'prot',
                                               config['CPU cores'])
        fo.pickle_dumper(self_scores, self_score_file)
        print('done.')
    else:
        self_scores = fo.pickle_loader(self_score_file)

    # create Kmer index for representatives
    print('Creating minimizer index for representative alleles...', end='')
    representatives = im.kmer_index(concat_reps, 5)
    print('done.')
    print('Created index with {0} distinct minimizers for {1} loci.'.format(len(representatives), len(loci_files)))

    # cluster CDSs into representative clusters
    print('Clustering proteins...')
    # define input group size based on number of available CPU cores
    group_size = math.ceil(len(proteins)/config['CPU cores'])
    cs_results = cf.cluster_sequences(proteins,
                                      config['Word size'],
                                      config['Window size'],
                                      config['Clustering similarity'],
                                      representatives,
                                      False,
                                      1,
                                      30,
                                      clustering_dir,
                                      config['CPU cores'],
                                      group_size,
                                      False)

    # exclude singletons
    clusters = {k: v for k, v in cs_results.items() if len(v) > 0}
    print('\nClustered {0} proteins into {1} clusters.'
          ''.format(len(proteins), len(clusters)))

    # create Fasta file with remaining proteins and representatives
    all_prots = fo.join_paths(clustering_dir, ['distinct_proteins.fasta'])
    fo.concatenate_files([unique_pfasta, concat_reps], all_prots)

    # create index for distinct protein sequences
    prot_index = fao.index_fasta(all_prots)

    # BLASTp if there are clusters with n>1
    excluded = []
    if len(clusters) > 0:
        blasting_dir = fo.join_paths(clustering_dir, ['cluster_BLASTer'])
        fo.create_directory(blasting_dir)

        all_proteins = fao.import_sequences(all_prots)
        # BLAST clustered sequences against cluster representatives
        print('Clusters to BLAST: {0}'.format(len(clusters)))
        blast_results, ids_dict = cf.blast_clusters(clusters, all_proteins,
                                                    blasting_dir, blastp_path,
                                                    makeblastdb_path,
                                                    config['CPU cores'],
                                                    True)

        blast_files = im.flatten_list(blast_results)

        # concatenate results for representatives of the same locus
        loci_results = {}
        for f in blast_files:
            locus_rep = ids_dict[im.match_regex(f, r'seq_[0-9]+')]
            locus_id = im.match_regex(locus_rep, r'^.+-protein[0-9]+')
            # for schemas that do not have 'protein' in the identifier
            # would fail for schemas from Chewie-NS
            if locus_id is None:
                locus_id = '_'.join(locus_rep.split('_')[0:-1])
            loci_results.setdefault(locus_id, []).append(f)

        concatenated_files = []
        concatenate_file_template = '{0}.concatenated_blastout.tsv'
        for locus, files in loci_results.items():
            outfile = fo.join_paths(blasting_dir,
                                    ['BLAST_results', concatenate_file_template.format(locus)])
            fo.concatenate_files(files, outfile)
            concatenated_files.append(outfile)

        loci_results = {}
        blast_matches_dir = fo.join_paths(clustering_dir, ['matches'])
        fo.create_directory(blast_matches_dir)
        for f in concatenated_files:
            locus_id = fo.get_locus_id(f)
            if locus_id is None:
                locus_id = fo.file_basename(f).split('.concatenated')[0]
            # exclude results in the BSR+0.1 threshold
            # process representative candidates in later stage
            best_matches = select_highest_scores(f)
            match_info = process_blast_results(best_matches, config['BLAST Score Ratio']+0.1,
                                               self_scores, ids_dict)
            locus_results = expand_matches(match_info, prot_index, dna_index,
                                           dna_distinct_htable, distinct_pseqids, basename_inverse_map)

            if len(locus_results) > 0:
                # save results to file
                locus_file = fo.join_paths(blast_matches_dir, ['{0}_locus_matches'.format(locus_id)])
                fo.pickle_dumper(locus_results, locus_file)
                loci_results[locus_id] = locus_file

        if len(loci_results) > 0:
            # process results per genome and per locus
            print('\nClassifying clustered proteins...')
            classification_inputs = []
            blast_clusters_results_dir = fo.join_paths(clustering_dir, ['results'])
            fo.create_directory(blast_clusters_results_dir)
            for locus, file in loci_results.items():
                # get locus length mode
                locus_mode = loci_modes[locus]

                # import file with locus classifications
                locus_results_file = fo.join_paths(classification_dir, [locus+'_results'])

                classification_inputs.append([locus, file,
                                              basename_inverse_map,
                                              locus_results_file, locus_mode,
                                              temp_directory,
                                              config['Size threshold'],
                                              config['BLAST Score Ratio'],
                                              blast_clusters_results_dir,
                                              config['CDS input'],
                                              classify_inexact_matches])

            class_results = mo.map_async_parallelizer(classification_inputs,
                                                      mo.function_helper,
                                                      config['CPU cores'],
                                                      show_progress=True)

            for r in class_results:
                current_results = fo.pickle_loader(r)
                for locus, v in current_results.items():
                    loci_modes[locus] = v[1]
                    excluded.extend(v[2])

            # may have repeated elements due to same CDS matching different loci
            excluded = set(excluded)

        print('\nClassified {0} distinct proteins.'.format(len(excluded)))

    # get seqids of remaining unclassified sequences
    unclassified_ids = [rec.id
                        for rec in SeqIO.parse(unique_pfasta, 'fasta')
                        if rec.id not in excluded]
    print('Unclassified proteins: {0}'.format(len(unclassified_ids)))

    if config['Mode'] == 3:
        template_dict['classification_files'] = classification_files
        template_dict['protein_fasta'] = all_prots
        template_dict['unclassified_ids'] = unclassified_ids

        return template_dict

    print('\n== Representative determination ==\n')
    # create directory to store data for each iteration
    iterative_rep_dir = fo.join_paths(temp_directory, ['6_representative_determination'])
    fo.create_directory(iterative_rep_dir)

    remaining_seqs_file = fo.join_paths(iterative_rep_dir, ['unclassified_proteins.fasta'])
    # create Fasta with unclassified sequences
    fao.get_sequences_by_id(prot_index, unclassified_ids,
                            remaining_seqs_file, limit=50000)

    # shorten ids to avoid BLASTp error?
    blast_db = fo.join_paths(iterative_rep_dir, ['blastdb'])
    # will not work if file contains duplicated seqids
    db_stderr = bw.make_blast_db(makeblastdb_path, remaining_seqs_file,
                                 blast_db, 'prot')

    # get seqids of schema representatives
    reps_ids = [rec.id for rec in SeqIO.parse(concat_reps, 'fasta')]

    # BLAST schema representatives against remaining unclassified CDSs
    new_reps = {}
    iteration = 1
    exausted = False
    # keep iterating while new representatives are discovered
    # this step can run faster if we align several representatives per core, instead of distributing 1 per core.
    while exausted is False:
        iter_string = 'Iteration {0}'.format(iteration)
        print(iter_string)
        print('='*len(iter_string))
        # create directory for current iteration
        iteration_directory = fo.join_paths(iterative_rep_dir, ['iteration_{0}'.format(iteration)])
        fo.create_directory(iteration_directory)
        # create text file with unclassified seqids
        remaining_seqids_file = fo.join_paths(iteration_directory, ['unclassified_seqids_{0}.txt'.format(iteration)])
        fo.write_lines(unclassified_ids, remaining_seqids_file)
        # BLAST representatives against remaining sequences
        # iterative process until the process does not detect new representatives
        print('Loci: {0}'.format(len(protein_repfiles)))

        # concatenate to create groups of 100 loci
        concat_repfiles = []
        for i in range(0, len(protein_repfiles), 100):
            concat_file = fo.join_paths(iteration_directory, ['concat_reps{0}-{1}.fasta'.format(i+1, i+len(protein_repfiles[i:i+100]))])
            fo.concatenate_files(protein_repfiles[i:i+100], concat_file)
            concat_repfiles.append(concat_file)

        # create BLASTp inputs
        output_files = []
        blast_inputs = []
        # create directory to store BLASTp results
        iteration_blast_dir = fo.join_paths(iteration_directory, ['BLAST_results'])
        fo.create_directory(iteration_blast_dir)
        for file in concat_repfiles:
            concat_rep_basename = fo.file_basename(file, False)
            outfile = fo.join_paths(iteration_blast_dir,
                                    [concat_rep_basename+'_blast_results_iter{0}.tsv'.format(iteration)])
            output_files.append(outfile)

            blast_inputs.append([blastp_path, blast_db, file, outfile,
                                 1, 1, remaining_seqids_file, 'blastp', 10,
                                 ct.IGNORE_RAISED, bw.run_blast])

        print('BLASTing loci representatives against unclassified proteins...', end='')
        # BLAST representatives against unclassified sequences
        blastp_results = mo.map_async_parallelizer(blast_inputs,
                                                   mo.function_helper,
                                                   config['CPU cores'],
                                                   show_progress=False)
        print('done.')

        loci_output_files = []
        for f in output_files:
            concat_results = fo.read_tabular(f)
            loci_separate_results = {}
            for line in concat_results:
                loci_separate_results.setdefault(rep_map[line[0]], []).append(line)
            for k, v in loci_separate_results.items():
                loci_separate_outfile = fo.join_paths(iteration_blast_dir, ['{0}_blast_results_iter{1}.tsv'.format(k, iteration)])
                lines = ['\t'.join(l) for l in v]
                fo.write_lines(lines, loci_separate_outfile)
                loci_output_files.append(loci_separate_outfile)

        loci_results = {}
        # create directory to store files with matches
        iteration_matches_dir = fo.join_paths(iteration_directory, ['matches'])
        fo.create_directory(iteration_matches_dir)
        for f in loci_output_files:
            locus_id = fo.get_locus_id(f)
            if locus_id is None:
                locus_id = fo.file_basename(f).split('_blast')[0]
            best_matches = select_highest_scores(f)
            match_info = process_blast_results(best_matches, config['BLAST Score Ratio'],
                                               self_scores, None)
            locus_results = expand_matches(match_info, prot_index, dna_index,
                                           dna_distinct_htable, distinct_pseqids, basename_inverse_map)

            if len(locus_results) > 0:
                locus_file = fo.join_paths(iteration_matches_dir, ['{0}_locus_matches'.format(locus_id)])
                fo.pickle_dumper(locus_results, locus_file)
                loci_results[locus_id] = locus_file

        print('Loci with high-scoring matches: {0}'.format(len(loci_results)))

        if len(loci_results) == 0:
            exausted = True
            continue

        # process results per genome and per locus
        classification_inputs = []
        blast_iteration_results_dir = fo.join_paths(iteration_directory, ['results'])
        fo.create_directory(blast_iteration_results_dir)
        for locus, file in loci_results.items():
            # get locus length mode
            locus_mode = loci_modes[locus]

            # import file with locus classifications
            locus_results_file = fo.join_paths(classification_dir, [locus+'_results'])

            classification_inputs.append([locus, file,
                                          basename_inverse_map,
                                          locus_results_file, locus_mode,
                                          temp_directory,
                                          config['Size threshold'],
                                          config['BLAST Score Ratio'],
                                          blast_iteration_results_dir,
                                          config['CDS input'],
                                          classify_inexact_matches])

        print('Classifying proteins...', end='')
        class_results = mo.map_async_parallelizer(classification_inputs,
                                                  mo.function_helper,
                                                  config['CPU cores'],
                                                  show_progress=False)

        # may have repeated elements due to same CDS matching different loci
        excluded = []
        representative_candidates = {}
        for r in class_results:
            current_results = fo.pickle_loader(r)
            for locus, v in current_results.items():
                loci_modes[locus] = v[1]
                excluded.extend(v[2])
                if len(v[3]) > 0:
                    representative_candidates[locus] = v[3]

        # remove representative candidates ids from excluded
        excluded = set(excluded)

        # include new representatives
        print('classified {0} proteins.'.format(len(excluded)))

        # exclude sequences that were excluded
        unclassified_ids = set(unclassified_ids) - excluded

        # create directory to store new representatives
        new_reps_directory = fo.join_paths(iteration_directory, ['representative_candidates'])
        fo.create_directory(new_reps_directory)

        representatives = {}
        representative_inputs = []
        if len(representative_candidates) > 0:
            candidates_dir = fo.join_paths(new_reps_directory, ['candidates'])
            fo.create_directory(candidates_dir)
            selection_dir = fo.join_paths(new_reps_directory, ['selection'])
            fo.create_directory(selection_dir)
            blast_selection_dir = fo.join_paths(selection_dir, ['BLAST_results'])
            fo.create_directory(blast_selection_dir)
            print('Selecting representatives for next iteration...', end='')
            for k, v in representative_candidates.items():
                if len(v) > 1:
                    current_candidates = {e[1]: e[3] for e in v}
                    fasta_file = fo.join_paths(candidates_dir,
                                               ['{0}_candidates.fasta'.format(k)])
                    # create file with sequences
                    fao.get_sequences_by_id(prot_index, list(current_candidates.keys()), fasta_file)
                    representative_inputs.append([current_candidates, k, fasta_file,
                                                  iteration, blast_selection_dir, blastp_path,
                                                  blast_db, config['BLAST Score Ratio'], 1,
                                                  select_representatives])
                else:
                    representatives[k] = [(v[0][1], v[0][3])]

            selected_candidates = mo.map_async_parallelizer(representative_inputs,
                                                            mo.function_helper,
                                                            config['CPU cores'],
                                                            show_progress=False)

            for c in selected_candidates:
                representatives[c[0]] = c[1]

            for k, v in representatives.items():
                new_reps.setdefault(k, []).extend(v)

            total_representatives = sum([len(v) for k, v in representatives.items()])
            print('selected {0} representatives.'.format(total_representatives))

        # new representatives and alleles that amtch in other genomes should have been all classified
        print('Unclassified proteins: {0}\n'.format(len(unclassified_ids)))

        # stop iterating if there are no new representatives
        if len(representatives) == 0:
            exausted = True
        else:
            selected_dir = fo.join_paths(new_reps_directory, ['selected'])
            fo.create_directory(selected_dir)
            # create files with representative sequences
            rep_map = {}
            reps_ids = []
            protein_repfiles = []
            for k, v in representatives.items():
                # get new representative for locus
                current_new_reps = [e[0] for e in v]
                reps_ids.extend(current_new_reps)
                for c in current_new_reps:
                    rep_map[c] = k

                # need to add 'short' or locus id will not be split
                rep_file = fo.join_paths(selected_dir,
                                         ['{0}_short_reps_iter.fasta'.format(k)])
                fao.get_sequences_by_id(prot_index, current_new_reps, rep_file)
                protein_repfiles.append(rep_file)

            # concatenate reps
            concat_repy = fo.join_paths(new_reps_directory, ['concat_reps.fasta'])
            fao.get_sequences_by_id(prot_index, set(reps_ids), concat_repy, limit=50000)
            # determine self-score for new reps
            candidates_blast_dir = fo.join_paths(new_reps_directory, ['representatives_self_score'])
            fo.create_directory(candidates_blast_dir)
            new_self_scores = cf.determine_self_scores(concat_repy, candidates_blast_dir,
                                                       makeblastdb_path, blastp_path,
                                                       'prot', config['CPU cores'])

            self_scores = {**self_scores, **new_self_scores}

        iteration += 1

    template_dict['classification_files'] = classification_files
    template_dict['protein_fasta'] = all_prots
    template_dict['unclassified_ids'] = unclassified_ids
    template_dict['self_scores'] = self_scores
    template_dict['representatives'] = new_reps

    return template_dict


def main(input_file, loci_list, schema_directory, output_directory,
         no_inferred, output_unclassified, output_missing,
         no_cleanup, hash_profiles, ns, config):

    for k, v in config.items():
        print('{0}: {1}'.format(k, v))

    # create directory to store intermediate files
    temp_directory = fo.join_paths(output_directory, ['temp'])
    fo.create_directory(temp_directory)

    start_time = pdt.get_datetime()

    # read file with paths to input files
    input_files = fo.read_lines(input_file, strip=True)
    # sort paths to FASTA files
    input_files = im.sort_iterable(input_files, sort_key=str.lower)
    print('Number of inputs: {0}'.format(len(input_files)))

    # get list of loci to call
    loci_to_call = fo.read_lines(loci_list)
    print('Number of loci: {0}'.format(len(loci_to_call)))

    # get list of schema loci
    schema_loci = fo.pickle_loader(fo.join_paths(schema_directory,
                                                 [ct.LOCI_LIST_FILE]))
    schema_loci_fullpath = [fo.join_paths(schema_directory, [file])
                            for file in schema_loci]

    # get size mode for all loci
    loci_modes_file = fo.join_paths(schema_directory, ['loci_modes'])
    if os.path.isfile(loci_modes_file) is True:
        loci_modes = fo.pickle_loader(loci_modes_file)
    else:
        print('\nDetermining sequence length mode for all loci...', end='')
        loci_modes = loci_mode_file(schema_loci_fullpath, loci_modes_file)
        print('done.')

    # check if schema contains folder with pre-computed hash tables
    pre_computed_dir = fo.join_paths(schema_directory, ['pre_computed'])
    if os.path.isdir(pre_computed_dir) is False:
        print('\nCreating pre-computed hash tables...', end='')
        # create hash tables to speedup exact matching and avoid schema
        # translation in each run
        fo.create_directory(pre_computed_dir)
        pre_compute_hash_tables(pre_computed_dir,
                                schema_loci_fullpath,
                                config['Translation table'],
                                config['CPU cores'])
        print('done.')

    # get index for loci to call
    loci_to_call = {file: schema_loci_fullpath.index(file)
                    for file in loci_to_call}

    results = allele_calling(input_files, schema_directory, temp_directory,
                             loci_modes.copy(), loci_to_call, config,
                             pre_computed_dir)

    print('\n== Wrapping up ==\n')

    # adjust missing locus classification based on mode
    classification_labels = ct.ALLELECALL_CLASSIFICATIONS
    if config['Mode'] != 4:
        classification_labels[-1] = ct.PROBABLE_LNF

    # sort for order similar to v2.0
    results['classification_files'] = dict(sorted(results['classification_files'].items()))

    # list files with CDSs coordinates
    if config['CDS input'] is False:
        coordinates_dir = fo.join_paths(temp_directory, ['2_cds_extraction'])
        coordinates_files = fo.listdir_fullpath(coordinates_dir, '.cds_hash')
        coordinates_files = {fo.file_basename(f, True).split('.cds_hash')[0]: f
                             for f in coordinates_files}
    else:
        coordinates_files = None

    print('Writing results_contigsInfo.tsv...', end='')
    repeated = write_results_contigs(list(results['classification_files'].values()),
                                     results['basename_map'],
                                     output_directory,
                                     coordinates_files,
                                     classification_labels)
    print('done.')

    # determine paralogous loci and write RepeatedLoci.txt file
    print('Writing paralogous_loci.tsv and paralogous_counts.tsv...', end='')
    total_paralogous = identify_paralogous(repeated, output_directory)
    print('done.')
    print('Detected number of paralogous loci: '
          '{0}'.format(total_paralogous))

    # assign allele identifiers to novel alleles
    assignment_inputs = list(results['classification_files'].items())
    assignment_inputs = [[g, ns, set(list(repeated.keys())), assign_allele_ids]
                         for g in assignment_inputs]

    novel_alleles = mo.map_async_parallelizer(assignment_inputs,
                                              mo.function_helper,
                                              config['CPU cores'],
                                              show_progress=False)

    novel_alleles = im.merge_dictionaries(novel_alleles)

    updated_files = {}
    if config['Mode'] != 1:
        # get seqids that match hashes
        for k, v in novel_alleles.items():
            for r in v:
                rep_seqid = im.polyline_decoding(results['dna_hashtable'][r[0]])[0:2]
                rep_seqid = '{0}-protein{1}'.format(results['basename_map'][rep_seqid[1]], rep_seqid[0])
                r.append(rep_seqid)

        reps_info = {}
        if config['Mode'] == 4:
            # get info for new representative alleles that must be added to files in the short directory
            for k, v in novel_alleles.items():
                locus_id = fo.get_locus_id(k)
                if locus_id is None:
                    locus_id = fo.file_basename(k, False)
                current_results = results['representatives'].get(locus_id, None)
                if current_results is not None:
                    for e in current_results:
                        allele_id = [line[1] for line in v if line[0] == e[1]]
                        # we might have representatives that were converted to NIPH but still appear in the list
                        if len(allele_id) > 0:
                            reps_info.setdefault(locus_id, []).append(list(e)+allele_id)

            if no_inferred is False:
                # update self_scores
                reps_to_del = set()
                for k, v in reps_info.items():
                    for r in v:
                        new_id = k+'_'+r[-1]
                        results['self_scores'][new_id] = results['self_scores'][r[0]]
                        # delete old entries
                        if r[0] not in reps_to_del:
                            reps_to_del.add(r[0])

                for r in reps_to_del:
                    del results['self_scores'][r]

                # save updated self-scores
                self_score_file = fo.join_paths(schema_directory, ['short', 'self_scores'])
                fo.pickle_dumper(results['self_scores'], self_score_file)

        if len(novel_alleles) > 0:
            # create Fasta files with novel alleles
            novel_directory = fo.join_paths(temp_directory, ['novel_alleles'])
            novel_rep_directory = fo.join_paths(novel_directory, ['short'])
            fo.create_directory(novel_rep_directory)
            added = create_novel_fastas(novel_alleles, reps_info, results['dna_fasta'], novel_directory)
            updated_files = added[2]
            if no_inferred is False:
                # add inferred alleles to schema
                added2 = add_inferred_alleles(added[2])
                # recompute mode for loci with novel alleles
                for file in novel_alleles:
                    alleles_sizes = list(fao.sequence_lengths(file).values())
                    # select first value in list if there are several values with same frequency
                    loci_modes[fo.file_basename(file, False)] = [sm.determine_mode(alleles_sizes)[0], alleles_sizes]
                fo.pickle_dumper(loci_modes, loci_modes_file)

                # add novel alleles hashes to pre-computed hash tables
                total_hashes = update_hash_tables(added[2], loci_to_call,
                                   config['Translation table'], pre_computed_dir)

    end_time = pdt.get_datetime()

    # create output files
    print('Writing logging_info.txt...', end='')
    write_logfile(start_time,
                  end_time,
                  len(results['basename_map']),
                  len(results['classification_files']),
                  config['CPU cores'],
                  config['BLAST Score Ratio'],
                  output_directory)
    print('done.')

    print('Writing results_alleles.tsv...', end='')
    profiles_table = write_results_alleles(list(results['classification_files'].values()),
                                           list(results['basename_map'].values()),
                                           output_directory,
                                           classification_labels[-1])
    print('done.')

    print('Writing results_statistics.tsv...', end='')
    write_results_statistics(results['classification_files'],
                             results['basename_map'],
                             output_directory,
                             classification_labels)
    print('done.')

    print('Writing loci_summary_stats.tsv...', end='')
    write_loci_summary(results['classification_files'],
                       output_directory,
                       len(input_files),
                       classification_labels)
    print('done.')

    if output_unclassified is True:
        print('Writing FASTA file with unclassified CDSs...', end='')
        # create Fasta file with the distinct CDSs that were not classified
        create_unclassified_fasta(results['dna_fasta'],
                                  results['protein_fasta'],
                                  results['unclassified_ids'],
                                  results['protein_hashtable'],
                                  output_directory,
                                  results['basename_map'])
        print('done.')

    if output_missing is True:
        # Create Fasta file with CDS that were classified as ASM, ALM, ...
        print('Writing FASTA file with CDSs classified as ASM, ALM, NIPH, '
              'NIPHEM, PAMA, PLOT3, PLOT5 and LOTSC...', end='')
        create_missing_fasta(results['classification_files'],
                             results['dna_fasta'],
                             results['basename_map'],
                             results['dna_hashtable'],
                             output_directory,
                             coordinates_files,
                             classification_labels)
        print('done.')

    if hash_profiles is not None:
        # create TSV file with hashed profiles
        print('Writing file with hashed profiles...', end='')
        ph.main(profiles_table, schema_directory, output_directory,
                hash_profiles, 4, 1000, updated_files,
                no_inferred)
        print('done.')

    # move file with CDSs coordinates
    # will not be created if input files contain predicted CDS
    if config['CDS input'] is False:
        fo.move_file(results['cds_coordinates'], output_directory)

    # move file with list of excluded CDS
    # file is not created if we only search for exact matches
    if config['Mode'] != 1:
        fo.move_file(results['invalid_alleles'], output_directory)

    # write Prodigal stderr for inputs that failed gene prediction
    if results['invalid_inputs'] is not None:
        failed_file = fo.join_paths(output_directory, ['prodigal_stderr.tsv'])
        lines = ['{0}\t{1}'.format(line[0], line[1])
                 for line in results['invalid_inputs']]
        fo.write_lines(lines, failed_file)

    # count total for each classification type
    global_counts, total_cds = count_classifications(results['classification_files'].values(),
                                                     classification_labels)

    print('Classified a total of {0} CDSs.'.format(total_cds))
    print('\n'.join(['{0}: {1}'.format(k, v)
                     for k, v in global_counts.items()]))
    if no_inferred is False and config['Mode'] != 1 and len(novel_alleles) > 0:
        print('Added {0} novel alleles to schema.'.format(sum(added2)))
        print('Added {0} representative alleles to schema.'.format(added[1]))
    elif no_inferred is True or len(novel_alleles) == 0:
        print('No new alleles to add to schema.')

    # remove temporary files
    if no_cleanup is False:
        fo.delete_directory(temp_directory)

    print('\nResults available in {0}'.format(output_directory))
