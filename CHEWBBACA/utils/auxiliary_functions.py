#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


import os
import sys
import time
import json
import shutil
import traceback
import multiprocessing
from collections import Counter
from multiprocessing import Pool
from SPARQLWrapper import SPARQLWrapper, JSON

from Bio import SeqIO

try:
    from utils import runProdigal
    from utils import io_utils as io
    from utils import str_utils as su
    from utils import list_utils as lu
    from utils import constants as cnst
    from utils import files_utils as fu
    from utils import fasta_utils as fau
    from utils import translation_utils as tu
except:
    from CHEWBBACA.utils import runProdigal
    from CHEWBBACA.utils import io_utils as io
    from CHEWBBACA.utils import str_utils as su
    from CHEWBBACA.utils import list_utils as lu
    from CHEWBBACA.utils import constants as cnst
    from CHEWBBACA.utils import files_utils as fu
    from CHEWBBACA.utils import fasta_utils as fau
    from CHEWBBACA.utils import translation_utils as tu


UNIPROT_SERVER = SPARQLWrapper("http://sparql.uniprot.org/sparql")


def extract_genome_cds(reading_frames, contigs, starting_id):
    """ Extracts CDSs from contigs based on the start
        and stop codon positions determined by Prodigal.

        Parameters
        ----------
        reading_frames : str
            Path to the ORF file created by Prodigal.
        contigs : dict
            Dictionary with contig ids as keys and contig
            sequences as values.
        starting_id : int
            Integer identifier to give to the first CDS extracted
            and that will be incremented to serve as identifier
            for subsequent CDSs.

        Returns
        -------
        coding_sequences : dict
            Dictionary with coding sequences ids as keys and
            coding sequences as values.
        coding_sequences_info : list
            List with a sublist for each extracted CDS. Sublists
            have information about the extracted CDS (identifier
            of the contig where the CDS was found, start position
            in the contig, stop position in the contig, sequence
            identifier attributed to that CDS and the strand that
            coded for that CDS).
    """

    seqid = starting_id
    coding_sequences = {}
    coding_sequences_info = []
    for contig_id, frames in reading_frames.items():
        sequence = contigs[contig_id]
        # for each start and stop codon in the contig
        for cds in frames:
            start_pos = cds[0]
            stop_pos = cds[1]
            strand = cds[2]
            # extract CDS sequence
            cds_sequence = su.extract_single_cds(sequence, *cds).upper()

            # store CDS with unique id
            coding_sequences[seqid] = cds_sequence

            # store CDS information
            coding_sequences_info.append([contig_id, str(start_pos),
                                          str(stop_pos), str(seqid),
                                          str(strand)])

            # increment seqid
            seqid += 1

    return [coding_sequences, coding_sequences_info]


def write_protein_table(output_file, genome_id, cds_info):
    """ Writes information about coding sequences in a
        genome to a file.

        Parameters
        ----------
        output_file : str
            Path to the output file to which info will
            be saved.
        genome_id : str
            Identifier of the genome to add to first field
            of every new line.
        cds_info : list
            List with information about each coding sequence
            identified in the genome (contig identifier,
            cds start position, cds stop position, cds
            identifier and cds coding strand).
    """

    table_lines = [[genome_id] + protein_info
                   for protein_info in cds_info]
    table_lines = [lu.join_list(line, '\t') for line in table_lines]
    table_text = lu.join_list(table_lines, '\n')
    io.write_to_file(table_text, output_file, 'a', '\n')


def save_extracted_cds(genome, identifier, orf_file, protein_table, cds_file):
    """ Extracts coding sequences from genome based on Prodigal
        gene predictions. Writes coding sequences to FASTA file
        and information about coding sequences to TSV file.

        Parameters
        ----------
        genome : str
            Path to the FASTA file with the FASTA sequences for
            a genome.
        identifier : str
            Genome identifier to add to FASTA records headers
            and to the first field in the TSV file.
        orf_file : str
            Path to the file with Prodigal results.
        protein_table : str
            Path to the TSV file to which coding sequences
            information will be written.
        cds_file : str
            Path to the FASTA file to which coding sequences
            will be written.

        Returns
        -------
        total_cds : int
            Total number of coding sequences extracted from
            the genome.
    """

    # import contigs for current genome/assembly
    contigs = fau.import_sequences(genome)
    # extract coding sequences from contigs
    reading_frames = io.pickle_loader(orf_file)
    genome_info = extract_genome_cds(reading_frames,
                                           contigs, 1)
    # save coding sequences to file
    # create records and write them to file
    cds_lines = fau.create_fasta_lines(genome_info[0], identifier)
    io.write_lines(cds_lines, cds_file)

    write_protein_table(protein_table, identifier, genome_info[1])

    total_cds = len(genome_info[0])

    return total_cds


def function_helper(input_args):
    """
    """

    try:
        results = input_args[-1](*input_args[0:-1])
    except Exception as e:
        func_name = (input_args[-1]).__name__
        traceback_lines = traceback.format_exception(etype=type(e), value=e,
                                                     tb=e.__traceback__)
        traceback_text = ''.join(traceback_lines)
        print('Error on {0}:\n{1}\n'.format(func_name, traceback_text))

    return results


def cds_batch_extractor(genomes, prodigal_path, temp_directory, index):
    """ Extracts coding sequences from a set of genomes.

        Parameters
        ----------
        input_data : list
            List with a set of paths for FASTA files with
            genomic sequences, followed by the path to the
            directory with files with Prodigal resutls, the
            path to the temporary directory for all files and
            directories that will be read and written and
            an index/identifier to add to the output files
            with coding sequences and coding sequences info.

        Returns
        -------
        A list with the following elements:
            protein_table : str
                Path to the TSV file to which coding sequences
                info was written.
            cds_file : str
                Path to the FASTA file to which coding sequences
                were written.
            batch_total : int
                Total number of coding sequences extracted from
                the set of input genomes.
    """

    protein_table = fu.join_paths(temp_directory,
                               ['protein_info_{0}.tsv'.format(index)])

    cds_file = fu.join_paths(temp_directory,
                          ['coding_sequences_{0}.fasta'.format(index)])

    batch_total = 0
    for g in genomes:
        # determine Prodigal ORF file path for current genome
        identifier = fu.file_basename(g, False)
        orf_file_path = fu.join_paths(prodigal_path,
                                   ['{0}_ORF.txt'.format(identifier)])
        total = save_extracted_cds(g, identifier, orf_file_path,
                                   protein_table, cds_file)
        batch_total += total

    return [protein_table, cds_file, batch_total]


def write_gene_list(schema_dir):
    """ Creates list with gene files in a schema and
        uses the pickle module to save the list to a file.

        Parameters
        ----------
        schema_dir : str
            Path to the directory with schema files.

        Returns
        -------
        A list with two elements. A boolean value that
        is True if the file with the list of genes was
        created and False otherwise. The second element
        is the path to the created file.
    """

    schema_files = [file for file in os.listdir(schema_dir) if '.fasta' in file]
    schema_list_file = fu.join_paths(schema_dir, ['.genes_list'])
    io.pickle_dumper(schema_files, schema_list_file)

    return [os.path.isfile(schema_list_file), schema_list_file]


def write_schema_config(blast_score_ratio, ptf_hash, translation_table,
                        minimum_sequence_length, chewie_version, size_threshold,
                        word_size, clustering_sim, representative_filter,
                        intra_filter, output_directory):
    """ Writes chewBBACA's parameters values used to create
        a schema to a file.

        Parameters
        ----------
        blast_score_ratio : float
            BLAST Score Ratio value used to create the
            schema.
        ptf_hash : str
            BLAKE2 hash of the Prodigal training file
            content.
        translation_table : int
            Genetic code used to predict and translate
            coding sequences.
        minimum_sequence_length : int
            Minimum sequence length, sequences with a
            length value lower than this value are not
            included in the schema.
        chewie_version : str
            Version of the chewBBACA suite used to create
            the schema.
        size_threshold : float
            Sequence size variation percentage threshold,
            new alleles cannot have a length value that
            deviates +/- than this value in relation to the
            locus's representative sequence.
        word_size : int
            Word/k value used to cluster protein sequences
            during schema creation and allele calling.
        clustering_sim : float
            Proportion of k-mers/minimizers that two proteins
            need to have in common to be clustered together.
        representative_filter : float
            Proportion of k-mers/minimizers that a clustered
            protein has to have in common with the representative
            protein of the cluster to be considered the same gene.
        intra_filter : float
            Proportion of k-mers/minimizers that clustered
            proteins have to have in common to be considered
            of the same gene.
        output_directory : str
            Path to the output directory where the file with
            schema parameters values will be created.

        Returns
        -------
        A list with two elements. A boolean value that
        is True if the file with the parameters values was
        created and False otherwise. The second element
        is the path to the created file.
    """

    size_threshold = None if size_threshold in [None, 'None'] else float(size_threshold)

    params = {}
    params['bsr'] = [float(blast_score_ratio)]
    params['prodigal_training_file'] = [ptf_hash]
    params['translation_table'] = [int(translation_table)]
    params['minimum_locus_length'] = [int(minimum_sequence_length)]
    params['chewBBACA_version'] = [chewie_version]
    params['size_threshold'] = [size_threshold]
    params['word_size'] = [word_size]
    params['cluster_sim'] = [clustering_sim]
    params['representative_filter'] = [representative_filter]
    params['intraCluster_filter'] = [intra_filter]

    config_file = os.path.join(output_directory, '.schema_config')
    io.pickle_dumper(params, config_file)

    return [os.path.isfile(config_file), config_file]


def select_name(result):
    """ Extracts the annotation description from the result
        of a query to the UniProt SPARQL endpoint.

        Parameters
        ----------
        result : dict
            A dictionary with the results from querying
            the UniProt SPARQL endpoint.

        Returns
        -------
        A list with the following elements:
            name : str
                The annotation descrition.
            url : str
                The URI to the UniProt page for the protein.
            label : str
                A label that has descriptive value.
    """

    url = ''
    name = ''
    label = ''

    i = 1
    found = False
    # get the entries with results
    aux = result['results']['bindings']
    total_res = len(aux)
    # only check results that are not empty
    if total_res > 0:
        # iterate over all results to find suitable
        while found is False:
            current_res = aux[i]
            res_keys = aux[i].keys()

            # annotation name can be associated
            # to different keys
            if 'fname' in res_keys:
                name = str(current_res['fname']['value'])
                found = True
            elif 'sname2' in res_keys:
                name = str(current_res['sname2']['value'])
                found = True
            elif 'label' in res_keys:
                name = str(current_res['label']['value'])
                found = True

            if 'label' in res_keys:
                label = str(current_res['label']['value'])
            else:
                label = name

            # get UniProt URL
            if 'uri' in res_keys:
                url = str(current_res['seq']['value'])
            elif 'seq' in res_keys:
                url = str(current_res['seq']['value'])

            if i == total_res:
                found = True

    return [name, url, label]


def uniprot_query(sequence):
    """ Constructs a SPARQL query to search for exact matches in
        UniProt's SPARQL endpoint.

        Parameters
        ----------
        sequence : str
            The Protein sequence that will be added to
            the query/searched for.

        Returns
        -------
        query : str
            The SPARQL query that will allow to search for
            exact matches in the UniProt database.
    """

    query = ('PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>  '
             'PREFIX up: <http://purl.uniprot.org/core/> '
             'select ?seq ?fname ?sname2 ?label  where {'
             '{?b a up:Simple_Sequence; rdf:value '
             '"'+sequence+'". ?seq up:sequence ?b. '
             'OPTIONAL{?seq up:submittedName ?sname. ?sname up:fullName ?sname2} '
             'OPTIONAL{?seq up:recommendedName ?rname.?rname up:fullName ?fname} }'
             'UNION{?seq a up:Sequence; rdf:value "'+sequence+'"; '
             'rdfs:label ?label. }}')

    return query


def verify_cpu_usage(cpu_to_use):
    """ Verify if the value provided for the number of CPU
        cores/threads does not exceed system limit or affect
        system performance.

        Parameters
        ----------
        cpu_to_use : int
            Value provided for the number of CPU cores/threads.

        Returns
        -------
        cpu_to_use : int
            Value of CPU cores/threads that will be used after
            determining if the provided value was safe.

    """
    total_cpu = multiprocessing.cpu_count()

    # do not allow a value greater than the number of cores
    if cpu_to_use >= total_cpu:
        print('Warning! You have provided a CPU core count value '
              'that is equal to or exceeds the number of CPU '
              'cores in your system!')
        # define a value that is safe according to the number of
        # available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1
        print('Resetting to: {0}'.format(cpu_to_use))
    elif cpu_to_use == (total_cpu - 1):
        print('Warning! You have provided a CPU core count value '
              'that is close to the maximum core count of your '
              'machine ({0}/{1}). This may affect your system '
              'responsiveness.'.format(cpu_to_use, total_cpu))

    return cpu_to_use


def check_ptf(ptf_path):
    """ Determines if path to Prodigal training file exists.

        Parameters
        ----------
        ptf_path : str
            Path to the Prodigal training file.

        Returns
        -------
        A list with a bool value that is True if the Prodigal
        training file exists or False otherwise and the path
        to the file if it exists or a message if it does not
        exist.
    """

    if os.path.isfile(ptf_path) is False:
        message = ('Cannot find specified Prodigal training file.'
                   '\nPlease provide a valid training file.\n\nYou '
                   'can create a training file for a species of '
                   'interest with the following command:\n  prodigal '
                   '-i <reference_genome> -t <training_file.trn> -p '
                   'single\n\nIt is strongly advised to provide a '
                   'high-quality and closed genome for the training '
                   'process.')
        return [False, message]
    else:
        return [True, ptf_path]


def check_prodigal_results(prodigal_results, output_directory):
    """ Determine if Prodigal could not predict genes for any input
        assembly.

        Parameters
        ----------
        prodigal_results : list
            List with gene prediction results from Prodigal.
        output_directory : str
            Path to the output directory where the file with information
            about failed cases will be written to.

        Returns
        -------
        A list with the following elements:
            failed : list
                List with the stderr for the cases that Prodigal
                failed to predict genes for.
            failed_file : str
                Path to the file with information about the failed
                cases.
    """

    no_cds = [l for l in prodigal_results if l[1] == 0]
    errors = [l for l in prodigal_results if isinstance(l[1], str) is True]
    failed = no_cds + errors

    failed_file = os.path.join(output_directory, 'prodigal_fails.tsv')
    if len(failed) > 0:
        lines = ['{0}\t{1}'.format(l[0], l[1]) for l in failed]
        io.write_lines(lines, failed_file)

    return [failed, failed_file]


def map_async_parallelizer(inputs, function, cpu, callback='extend',
                           chunksize=1, show_progress=False):
    """ Parallelizes function calls by creating several processes
        and distributing inputs.

        Parameters
        ----------
        inputs : list
            List with inputs to process.
        function
            Function to be parallelized.
        cpu : int
            Number of processes to create (based on the
            number of cores).
        callback : str
            Results can be appended, 'append', to the
            list that stores results or the list of results
            can be extended, 'extend'.
        chunksize : int
            Size of input chunks that will be passed to
            each process. The function will create groups
            of inputs with this number of elements.
        show_progress: bool
            True to show a progress bar with the percentage
            of results that have been processed, False
            otherwise.

        Returns
        -------
        results : list
            List with the results returned for each function
            call.
    """

    results = []
    pool = Pool(cpu)
    if callback == 'extend':
        rawr = pool.map_async(function, inputs,
                              callback=results.extend, chunksize=chunksize)
    elif callback == 'append':
        rawr = pool.map_async(function, inputs,
                              callback=results.append, chunksize=chunksize)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = progress_bar(rawr, len(inputs))

    rawr.wait()

    return results


def check_input_type(input_path, output_file):
    """ Checks if the input path is for a file or for a
        directory. If the path is for a directory, the
        function creates a file with the list of paths
        to FASTA files in the directory.

        Parameters
        ----------
        input_path : str
            Path to file or directory.
        output_file : str
            Path to the output file with the list of FASTA
            files.

        Returns
        -------
        list_files : str
            Path to a file with a list of paths for FASTA
            files.

        Raises
        ------
        SystemExit
            - If there were no FASTA files in the directory.
            - If input path is not a valid path for a file or
              for a directory.
    """

    # check if input argument is a file or a directory
    if os.path.isfile(input_path):
        list_files = input_path
    elif os.path.isdir(input_path):
        # we need to get only files with FASTA extension
        files = os.listdir(input_path)
        files = fu.filter_files(files, cnst.FASTA_SUFFIXES)
        # get absolute paths
        files = [os.path.join(input_path, file) for file in files]
        # filter any directories that migh end with FASTA extension
        files = [file for file in files if os.path.isdir(file) is False]

        # only keep files whose content is typical of a FASTA file
        fasta_files = fau.filter_non_fasta(files)

        # if there are FASTA files
        if len(fasta_files) > 0:
            # store full paths to FASTA files
            with open(output_file, 'w') as f:
                for file in fasta_files:
                    f.write(file + '\n')
        else:
            sys.exit('\nCould not get input files. Please '
                     'provide a directory with FASTA files '
                     'or a file with the list of full paths '
                     'to the FASTA files and ensure that '
                     'filenames end with one of the '
                     'following suffixes: {0}.'
                     ''.format(cnst.FASTA_SUFFIXES))

        list_files = output_file
    else:
        sys.exit('\nInput argument is not a valid directory or '
                 'file with a list of paths. Please provide a '
                 'valid input, either a folder with FASTA files '
                 'or a file with the list of full paths to FASTA '
                 'files (one per line).')

    return list_files


def determine_duplicated_seqs(sequences):
    """ Creates a dictionary with sequences as keys and all sequence
        identifiers associated with a sequence as values.

        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.

        Returns
        -------
        equal_seqs : dict
            Dictionary with sequences as keys and sequence
            identifiers that are associated with each
            sequence as values.
    """

    equal_seqs = {}
    for seqid, seq in sequences.items():
        # if protein sequence was already added as key
        if seq in equal_seqs:
            # append new protid
            equal_seqs[seq].append(seqid)
        # else add new protein sequence as key and protid
        # as value
        else:
            equal_seqs[seq] = [seqid]

    return equal_seqs


def determine_longest(seqids, sequences):
    """ Determines which sequence is the longest among
        sequences with the specified identifiers.

        Parameters
        ----------
        seqids : list
            List with sequence identifiers.
        sequences : dict
            Dictionary with sequence identifiers as keys
            and sequences as values.

        Returns
        -------
        chosen : str
            Sequence identifier of the longest sequence.
    """

    seqids_tups = [(seqid, sequences[seqid]) for seqid in seqids]
    sorted_tups = sorted(seqids_tups, key=lambda x: len(x[1]), reverse=True)
    chosen = sorted_tups[0][0]

    return chosen


def determine_mode(int_values):
    """ Determines the mode value from a list of integer values.
        Returns list with multiple values if distribution is
        multimodal.

        Parameters
        ----------
        int_values : list
            List with integer values.

        Returns
        -------
        modes : list
            The most frequent integer values.
    """

    # determine frequency of each length value
    counts = Counter(int_values)

    # order by most common first
    most_common = counts.most_common()

    # get most common
    modes = [most_common[0][0]]

    # determine if there are more length values that are as common
    modes += [m[0] for m in most_common[1:] if m[1] == most_common[0][1]]

    return modes


def mode_filter(sequences, size_threshold):
    """ Determines the mode from a set of input sequences
        and identifies sequences that have a length value
        smaller or greater than the mode based on a threshold.

        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        size_threshold : float
            Sequences with +/- this value * mode will be
            reported as above or below the mode.

        Returns
        -------
        A list with the following variables:
            modes : list
                List with mode values determined based on the
                length of input sequences.
            alm : list
                List with the sequence identifiers of the
                sequences that are above the mode value by
                mode*size_threshold.
            asm : list
                List with the sequence identifiers of the
                sequences that are below the mode value by
                mode*size_threshold.
            seqs_lengths : dict
                Dictionary with sequence identifiers as keys
                and sequence lengths as values.
    """

    # determine length value of all sequences
    seqs_lengths = {seqid: len(seq) for seqid, seq in sequences.items()}

    # determine mode/s
    modes = determine_mode(list(seqs_lengths.values()))

    # determine top and bot length value limits
    max_mode = max(modes)
    top_limit = max_mode + (max_mode*size_threshold)
    min_mode = min(modes)
    bot_limit = min_mode - (min_mode*size_threshold)

    # determine sequences that are below or above limits
    alm = [seqid for seqid, length in seqs_lengths.items()
           if length > top_limit]
    asm = [seqid for seqid, length in seqs_lengths.items()
           if length < bot_limit]

    return [modes, alm, asm, seqs_lengths]


def add_prefix(ids, prefix):
    """
    """

    ids_map = {}
    for i in ids:
        new_id = '{0}_{1}'.format(prefix, i.split('_')[-1])
        ids_map[i] = new_id

    return ids_map


def get_seqs_dicts(fasta_path, gene_id, table_id, min_len, size_threshold):
    """ Creates a dictionary mapping seqids to DNA sequences and
        another dictionary mapping protids to protein sequences.

        Parameters
        ----------
            gene_file (str): path/name of the FASTA file with
            DNA sequences.
            table_id (int): translation table identifier.

        Returns:
            List with following elements:
                dna_seqs (dict): dictionary with sequence identifiers as keys
                and DNA sequences as values.
                prot_seqs (dict): dictionary with protein identifiers as keys
                and Protein sequences as values.
                invalid_alleles (list): list with sequence identifiers of
                alleles that are not valid because they could not be
                translated.
    """

    sequences = fau.import_sequences(fasta_path)

    # translate sequences
    translated_seqs = {k: tu.translate_dna(v, table_id, min_len)
                       for k, v in sequences.items()}

    # add locus identifier to headers
    # some headers might only have the allele identifier
    seqids = list(translated_seqs.keys())
    new_seqids = add_prefix(seqids, gene_id)
    # switch ids
    sequences = {new_seqids[k]: v for k, v in translated_seqs.items()}

    valid = {k: v[0][1] for k, v in sequences.items() if isinstance(v, list) is True}
    invalid = [[k, v] for k, v in sequences.items() if isinstance(v, list) is False]

    seqid = 1
    seqids_map = {}
    dna_seqs = {}
    prot_seqs = {}
    for k, v in valid.items():
        seqids_map[str(seqid)] = k
        dna_seqs[k] = v[0][1]
        prot_seqs[str(seqid)] = str(v[0][0])
        seqid += 1

    if size_threshold is not None and len(prot_seqs) > 0:
        # remove alleles based on length mode and size threshold
        modes, alm, asm, alleles_lengths = mode_filter(dna_seqs, size_threshold)
        excluded = set(asm + alm)

        dna_seqs = {seqid: seq for seqid, seq in dna_seqs.items() if seqid not in excluded}
        prot_seqs = {seqid: seq for seqid, seq in prot_seqs.items() if seqids_map[seqid] not in excluded}

        modes_concat = ':'.join(map(str, modes))
        st_percentage = int(size_threshold*100)
        invalid += [[s, 'allele greater than {0}% locus length mode '
                                '({1}>{2})'.format(st_percentage, alleles_lengths[s], modes_concat)] for s in alm]
        invalid += [[s, 'allele smaller than {0}% locus length mode '
                                '({1}<{2})'.format(st_percentage, alleles_lengths[s], modes_concat)] for s in asm]

    total_seqs = len(translated_seqs)

    return [dna_seqs, prot_seqs,
            invalid, seqids_map, total_seqs]


def split_genes_by_core(inputs, cores, method):
    """ Creates balanced lists of loci to distribute per number
        of available cores. Loci lists can be created based
        on the number of sequence per locus (seqcount), the mean
        length of the sequences in each locus or the product of
        both values.

        Parameters
        ----------
        inputs : list
            List with one sublist per locus. Each sublist has
            a locus identifier, the total number of sequences
            and sequence mean legth for that locus.
        cores : int
            The number of loci groups that should be created.
            Based on the number of CPU cores that will be
            used to process the inputs.
        method : str
            "seqcount" to create loci lists based on the total
            number of sequences, "length" to split based
            on mean length of sequences and "seqcount+length" to
            split based on both criteria.

        Returns
        -------
        splitted_ids : list
            List with sublists that contain loci identifiers
            Sublists are balanced based on the chosen method.
    """

    # initialize list with sublists to store inputs
    splitted_ids = [[] for cpu in range(cores)]
    # initialize list with chosen criterion values
    # for each sublist of inputs
    splitted_values = [0 for cpu in range(cores)]
    i = 0
    for locus in inputs:
        if method == 'seqcount':
            splitted_values[i] += locus[1]
        elif method == 'length':
            splitted_values[i] += locus[2]
        elif method == 'seqcount+length':
            splitted_values[i] += locus[1] * locus[2]
        splitted_ids[i].append(locus[0])
        # at the end of each iteration, choose the sublist
        # with lowest criterion value
        i = splitted_values.index(min(splitted_values))

    return splitted_ids


def read_configs(schema_path, filename):
    """ Reads file with schema config values.

        Parameters
        ----------
        schema_path : str
            Path to the schema's directory.
        filename : str
            Name of the file that contains the config values.

        Returns
        -------
        configs : dict
            Dictionary with config names as keys and config
            values as values.
    """

    config_file = os.path.join(schema_path, filename)
    if os.path.isfile(config_file):
        # Load configs dictionary
        configs = io.pickle_loader(config_file)
    else:
        sys.exit('Could not find a valid config file.')

    return configs


def get_data(sparql_query):
    """ Sends request to query UniProts's SPARQL endpoint.

        Parameters
        ----------
        sparql_query : str
            SPARQL query.

        Returns
        -------
        result : dict
            Dictionary with data retrieved from UniProt.
    """

    tries = 0
    max_tries = 5
    success = False
    while success is False and tries < max_tries:
        try:
            UNIPROT_SERVER.setQuery(sparql_query)
            UNIPROT_SERVER.setReturnFormat(JSON)
            UNIPROT_SERVER.setTimeout(60)
            result = UNIPROT_SERVER.query().convert()
            success = True
        except Exception as e:
            tries += 1
            result = e
            time.sleep(1)

    return result


def progress_bar(process, total, tickval=5, ticknum=20, completed=False):
    """ Creates and prints progress bar to stdout.

        Parameters
        ----------
        process : multiprocessing.pool.MapResult
            Multiprocessing object.
        total : int
            Total number of inputs that have to be processed.
        tickval : int
            Progress completion percentage value for each
            tick.
        ticknum : int
            Total number of ticks in progress bar.
        completed : bool
            Boolean indicating if process has completed.

        Returns
        -------
        completed : bool
            Boolean indicating if process has completed.
    """

    # check if process has finished
    if (process.ready()):
        # print full progress bar and satisfy stopping condition
        progress_bar = '[{0}] 100%'.format('='*ticknum)
        completed = True

    # check how many inputs have been processed
    remaining = process._number_left
    if remaining == total:
        # print empty progress bar
        progress_bar = '[{0}] 0%'.format(' '*ticknum)
    else:
        # print progress bar, incremented by 5%
        progress = int(100-(remaining/total)*100)
        progress_tick = progress//tickval
        progress_bar = '[{0}{1}] {2}%'.format('='*progress_tick,
                                              ' '*(ticknum-progress_tick),
                                              progress)

    print('\r', progress_bar, end='')
    time.sleep(0.5)

    return completed


def translate_coding_sequences(seqids, sequences_file, translation_table,
                               minimum_length, dna_file, protein_file):
    """ Translates CDSs into protein sequences.

        Parameters
        ----------
        input_data : list
            A list with the sequence identifiers of the
            sequences that should be translated, the path
            to the FASTA file that contains the DNA sequences,
            the translation table identifier, the minimum
            sequence length value, the path to a file to
            save DNA sequences and the path to a file to
            save protein sequences.

        Returns
        -------
        A list with following elements:
            invalid_alleles : list
                List with one sublist per invalid allele.
                Each sublist contains a sequence identifer
                and the exception message returned after
                attempting translation.
            total_seqs : int
                Total number of DNA sequences that were
                translated.
    """

    # define limit of records to keep in memory
    dna_lines = []
    total_seqs = 0
    prot_lines = []
    line_limit = 5000
    invalid_alleles = []

    cds_index = SeqIO.index(sequences_file, 'fasta')

    for i, seqid in enumerate(seqids):
        try:
            sequence = str(cds_index.get(seqid).seq)
        except Exception as e:
            print(e)

        translation = tu.translate_dna(sequence, translation_table, minimum_length)
        if isinstance(translation, list):
            dna_lines.append('>{0}'.format(seqid))
            dna_lines.append(translation[0][1])
            prot_lines.append('>{0}'.format(seqid))
            prot_lines.append(str(translation[0][0]))
            total_seqs += 1
        # if returned value is a string, translation failed and
        # string contains exceptions
        elif isinstance(translation, str):
            invalid_alleles.append([seqid, translation])

        if len(dna_lines)//2 == line_limit or i+1 == len(seqids):

            dna_lines = lu.join_list(dna_lines, '\n')
            io.write_to_file(dna_lines, dna_file, 'a', '\n')
            dna_lines = []

            prot_lines = lu.join_list(prot_lines, '\n')
            io.write_to_file(prot_lines, protein_file, 'a', '\n')
            prot_lines = []

    return [invalid_alleles, total_seqs]


def determine_distinct(sequences_file, unique_fasta):
    """ Identifies duplicated sequences in a FASTA file.
        Returns a single sequence identifier per distinct
        sequence and saves distinct sequences to a FASTA
        file.

        Parameters
        ----------
        sequences_file : str
            Path to a FASTA file.
        unique_fasta : str
            Path to a FASTA file that will be created to
            store distinct sequences.

        Returns
        -------
        List with following elements:
            total : int
                Total number of times sequences were repeated.
            unique_seqids : list
                List with one sequence identifier per distinct
                sequence. The first identifier observed for a
                distinct sequence is the one stored in the list.
    """

    total = 0
    seqs_dict = {}
    out_limit = 10000
    out_seqs = []
    exausted = False
    seq_generator = SeqIO.parse(sequences_file, 'fasta')
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            # seq object has to be converted to string
            sequence = str(record.seq.upper())
            seqid = record.id
            seq_hash = su.hash_sequence(sequence)

            # store only the hash for distinct sequences
            if seq_hash not in seqs_dict:
                seqs_dict[seq_hash] = seqid
                recout = fau.fasta_str_record(seqid, sequence)
                out_seqs.append(recout)
            elif seq_hash in seqs_dict:
                total += 1
        else:
            exausted = True

        if len(out_seqs) == out_limit or exausted is True:
            if len(out_seqs) > 0:
                out_seqs = lu.join_list(out_seqs, '\n')
                io.write_to_file(out_seqs, unique_fasta, 'a', '\n')
                out_seqs = []

    unique_seqids = list(seqs_dict.values())

    return [total, unique_seqids]


def determine_small(sequences_file, minimum_length):
    """ Find protein sequences that are shorter than
        desired length.

        Parameters
        ----------
        sequences_file : str
            Path to a FASTA file.
        minimum_length : int
            Sequences with a length value below this value
            are considered small.

        Returns
        -------
        small_seqs : list
            List with the identifiers of small sequences.
    """

    small_seqids = []
    for record in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        sequence = str(record.seq)
        seqid = record.id

        if len(sequence) < minimum_length:
            small_seqids.append(seqid)

    return small_seqids


def create_short(schema_files, schema_dir):
    """ Creates the 'short' directory for a schema.
        Creates the directory and copies schema files
        to the directory (should be used when the schema
        only has 1 sequence per gene/locus).

        Parameters
        ----------
        schema_files : list
            List with paths to all FASTA files in the schema.
        schema_dir : str
            Path to the schema's directory.

        Returns
        -------
        True on completion.
    """

    short_path = fu.join_paths(schema_dir, ['short'])
    fu.create_directory(short_path)

    for file in schema_files:
        short_file = fu.join_paths(short_path, [fu.file_basename(file)])
        short_file = short_file.replace('.fasta', '_short.fasta')
        shutil.copy(file, short_file)

    return True


def apply_bsr(blast_results, fasta_file, bsr, ids_dict):
    """ Computes the BLAST Score Ratio value for BLAST
        alignments and returns the identifiers of the
        sequences that are similar to sequences with
        the same size or that are larger.

        Parameters
        ----------
        inputs : list
            List with the path to a file with BLAST
            results, the path to a FASTA file that
            contains the sequences that were aligned,
            the BSR value to use as threshold and a
            dictionary with the mapping between
            sequence identifiers used for BLAST and
            the original sequence identifiers.

        Returns
        -------
        excluded_alleles : list
            List with the identifiers of the sequences
            that were highly similar to larger sequences
            or sequences of the same size.
    """

    self_scores = {r[0]: r[2] for r in blast_results if r[0] == r[1]}
    # do not include self-scores lines, no need to evaluate those hits
    blast_results = [r for r in blast_results if r[0] != r[1]]

    lengths = {}
    for k in self_scores:
        record = fasta_file.get(ids_dict[k])
        sequence = str(record.seq)
        lengths[k] = len(sequence)

    excluded_alleles = []
    for res in blast_results:

        query = res[0]
        hit = res[1]
        score = res[2]

        if query not in excluded_alleles:
            # try to apply BSR strategy
            try:
                self_blast_score = self_scores[query]

                query_length = lengths[query]
                hit_length = lengths[hit]
                blast_score_ratio = float(score) / float(self_blast_score)

                # BSR has to be greater than threshold, just as in the original function
                if blast_score_ratio >= bsr and hit not in excluded_alleles:

                    if hit_length > query_length and query not in excluded_alleles:
                        excluded_alleles.append(query)

                    elif hit_length <= query_length:
                        excluded_alleles.append(hit)
            # it might not work because there is no self score for
            # some sequences due to low complexity regions...
            except Exception:
                excluded_alleles.append(query)

    return excluded_alleles


def sort_data(data, sort_key=None, reverse=False):
    """ Sorts an iterable.

        Parameters
        ----------
        data : iter
            Iterable to sort.
        sort_key
            If provided, data will be sorted based
            on this function.
        reverse : bool
            If sorting order should be inverted.

        Returns
        -------
        sorted_data
            List with sorted elements.
    """

    if sort_key is None:
        sorted_data = sorted(data, reverse=reverse)
    elif sort_key is not None:
        sorted_data = sorted(data, key=sort_key, reverse=reverse)

    return sorted_data


def kmer_index(sequences, word_size):
    """
    """

    kmers_mapping = {}
    for seqid, seq in sequences.items():
        minimizers = determine_minimizers(seq, word_size, word_size, position=False)
        kmers = set(minimizers)

        # create dict with kmers as keys and list
        # of sequences with given kmers as values
        for kmer in kmers:
            kmers_mapping.setdefault(kmer, []).append(seqid)

    return kmers_mapping
