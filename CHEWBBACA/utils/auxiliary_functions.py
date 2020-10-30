#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


import os
import re
import sys
import csv
import time
import json
import shutil
import pickle
import hashlib
import zipfile
import requests
import itertools
import subprocess
import datetime as dt
import multiprocessing
from getpass import getpass
from itertools import islice
from collections import Counter
from multiprocessing import Pool
from multiprocessing import TimeoutError
from multiprocessing.pool import ThreadPool
from SPARQLWrapper import SPARQLWrapper, JSON
from urllib.parse import urlparse, urlencode, urlsplit, parse_qs

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

try:
    from utils import runProdigal
    from utils import constants as cnst
except:
    from CHEWBBACA.utils import runProdigal
    from CHEWBBACA.utils import constants as cnst


UNIPROT_SERVER = SPARQLWrapper("http://sparql.uniprot.org/sparql")


def read_lines(input_file, strip=True):
    """ Reads lines in an input file and stores those lines
        in a list.

        Parameters
        ----------
        input_file : str
            Path to the input file.
        strip : bool
            Specify if lines should be stripped of leading
            and trailing white spaces and new line characters.

        Returns
        -------
        lines : list
            List with the lines read from the input file.
    """

    with open(input_file, 'r') as infile:
        if strip is True:
            lines = [file.strip() for file in infile.readlines()]
        else:
            lines = [file for file in infile.readlines()]

    return lines


def join_list(lst, link):
    """ Joins all elements in a list into a single string.

        Parameters
        ----------
        lst : list
            List with elements to be joined.
        link : str
            Character used to join list elements.

        Returns
        -------
        joined_list : str
            A single string with all elements in the input
            list joined by the character chosen as link.
    """

    joined_list = link.join(lst)

    return joined_list


def file_basename(file_path, suffix=True):
    """ Extract file basename from path.

        Parameters
        ----------
        file_path : str
            Path to the file.
        suffix : bool
            Specify if the basename should include the file
            extension.

        Returns
        -------
        basename : str
            File basename extracted from input path.
    """

    basename = os.path.basename(file_path)

    if suffix is False:
        basename = basename.split('.')[0]

    return basename


def check_connection(ns_url, headers=cnst.HEADERS_GET_JSON):
    """ Verifies connection to a chewie-NS instance.

        Parameters
        ----------
        ns_url : str
           The base URL for a Chewie Nomenclature Server
           instance.
        headers : dict
            HTTP headers for GET requests.

        Returns
        -------
        conn : bool
            True if it was possible to return data from the
            server, False otherwise.
    """

    url = make_url(ns_url, *['stats', 'summary'])

    try:
        res = requests.get(url, headers=headers, timeout=30, verify=False)
        server_status = res.status_code
        if server_status in [200, 201]:
            conn = True
        else:
            conn = False
    except Exception:
        conn = False

    return conn


def pickle_dumper(content, output_file):
    """ Use the Pickle module to serialize an object.

        Parameters
        ----------
        content : type
            Variable that refers to the object that will
            be serialized and written to the output file.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'wb') as po:
        pickle.dump(content, po)


def pickle_loader(input_file):
    """ Use the Pickle module to de-serialize an object.

        Parameters
        ----------
        input_file : str
            Path to file with byte stream to be de-serialized.

        Returns
        -------
        content : type
            Variable that refers to the de-serialized
            object.
    """

    with open(input_file, 'rb') as pi:
        content = pickle.load(pi)

    return content


def file_zipper(input_file, zip_file):
    """ Zips (compresses) a file.

        Parameters
        ----------
        input_file : str
            Path to the file that will be compressed.
        zip_file : str
            Path to the ZIP file that will be created.

        Returns
        -------
        zip_file : str
            Path to the ZIP file that was created by
            compressing the input file.
    """

    with zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(input_file, os.path.basename(input_file))

    return zip_file


def remove_files(files):
    """ Deletes a list of files.

        Parameters
        ----------
        files : list
            List with paths to the files to be deleted.
    """

    for f in files:
        os.remove(f)


def count_sequences(fasta_file):
    """ Counts the number of sequences in a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file.

        Returns
        -------
        total_seqs : int
            Number of sequences in the input FASTA file.
    """

    records = SeqIO.parse(fasta_file, 'fasta')
    total_seqs = len(list(records))

    return total_seqs


def simple_get_request(base_url, headers, endpoint_list):
    """ Constructs an URI for an endpoint and requests
        data through a GET method sent to that endpoint.

        Parameters
        ----------
        base_url : str
            The base URL for a Chewie Nomenclature Server
            instance.
        headers : dict
            HTTP headers for GET requests.
        endpoint_list : list
            List with elements that will be concatenated
            to the base URI to create the URI for the API endpoint.

        Returns
        -------
        res : requests.models.Response
            Response object from the GET method.
    """

    # unpack list of sequential endpoints and pass to create URI
    url = make_url(base_url, *endpoint_list)
    res = requests.get(url, headers=headers, timeout=30, verify=False)

    return res


def simple_post_request(base_url, headers, endpoint_list, data):
    """ Constructs an URI for an endpoint and sends
        data through a POST method sent to that endpoint.

        Parameters
        ----------
        base_url : str
            The base URL for a Chewie Nomenclature Server
            instance.
        headers : dict
            HTTP headers for GET requests.
        endpoint_list : list
            List with elements that will be concatenated
            to the base URI to create the URI for the API endpoint.
        data
            Data to send. Must be an object that can be JSON
            serialized.

        Returns
        -------
        res : requests.models.Response
            Response object from the POST method.
    """

    # unpack list of sequential endpoints and pass to create URI
    url = make_url(base_url, *endpoint_list)
    res = requests.post(url, data=json.dumps(data),
                        headers=headers, verify=False)

    return res


def concatenate_files(files, output_file, header=None):
    """ Concatenates the contents of a set of files.

        Parameters
        ----------
        files : list
            List with the paths to the files to concatenate.
        output_file : str
            Path to the output file that will store the
            concatenation of input files.
        header : str or NoneType
            Specify a header that should be written as the
            first line in the output file.

        Returns
        -------
        output_file : str
            Path to the output file that was created with
            the concatenation of input files.
    """

    with open(output_file, 'w') as of:
        if header is not None:
            of.write(header)
        for f in files:
            with open(f, 'r') as fd:
                shutil.copyfileobj(fd, of)

    return output_file


def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single string to a file.

        Parameters
        ----------
        text : str
            A single string to write to the output file.
        output_file : str
            Path to the output file.
        write_mode : str
            Write mode can be 'w', writes text and overwrites
            any text in file, or 'a', appends text to text
            already in file.
        end_char : str
            Character added to the end of the file.
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def write_lines(lines, output_file):
    """ Writes a list of strings to a file. The strings
        are joined with newlines before being written to
        file.

        Parameters
        ----------
        lines : list
            List with the lines/strings to write to the
            output file.
        output_file : str
            Path to the output file.
    """

    joined_lines = join_list(lines, '\n')

    write_to_file(joined_lines, output_file, 'a', '\n')


def write_records(records, output_file):
    """ Writes FASTA records (BioPython SeqRecord) to a file.

        Parameters
        ----------
        records : list
            List with BioPython SeqRecord objects.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'w') as output_handle:
        fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(records)


def integer_headers(input_fasta, output_fasta, start=1, limit=5000):
    """ Switches FASTA records headers in a file by integer
        values.

        Parameters
        ----------
        input_fasta : str
            Path to the a FASTA file.
        output_fasta : str
            Path to the output file with modified headers.

        Returns
        -------
        ids_map : dict
            Dictionary with mapping between integer and original
            headers.
    """

    seqs = []
    ids_map = {}
    exausted = False
    seq_generator = SeqIO.parse(input_fasta, 'fasta')
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            new_id = 'seq_{0}'.format(start)
            ids_map[new_id] = record.id
            sequence = str(record.seq)
            new_rec = '>{0}\n{1}'.format(new_id, sequence)
            seqs.append(new_rec)
            start += 1
        elif record is None:
            exausted = True

        if len(seqs) == limit or exausted is True:
            write_lines(seqs, output_fasta)
            seqs = []

    return ids_map


def create_fasta_lines(sequences, prefix):
    """ Creates FASTA records in string format.

        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys
            and sequences as values.
        prefix : str
            Prefix to include in sequences headers.

        Returns
        -------
        lines : list
            List with Fasta records in string format.
    """

    template = '>{0}-protein{1}\n{2}'

    lines = [template.format(prefix, seqid, sequence)
             for seqid, sequence in sequences.items()]

    return lines


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.

        Returns
        -------
        seqs_dict : dict
            Dictionary that has sequences ids as keys and
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta')
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def extract_subsequence(sequence, start, stop):
    """ Extract substring from string.

        Parameters
        ----------
        sequence : str
            Input string.
        start : int
            Substring start position in input string.
        stop : int
            Substring stop position in input string.

        Returns
        -------
        subsequence : str
            Substring extracted from input string.
    """

    subsequence = sequence[start:stop]

    return subsequence


def extract_single_cds(sequence, start, stop, strand):
    """ Extract coding sequence from contig.

        Parameters
        ----------
        sequence : str
            Contig that contains the coding sequence.
        start : int
            Coding sequence start position in contig.
        stop : int
            Coding sequence stop position in contig.
        strand : int
            Conding sequence orientation.

        Returns
        -------
        coding_sequence : str
            Coding sequence extracted from input contig
            in sense orientation.
    """

    coding_sequence = extract_subsequence(sequence, start, stop)

    if strand == 0:
        coding_sequence = reverse_complement(coding_sequence)

    return coding_sequence


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
            cds_sequence = extract_single_cds(sequence, *cds).upper()

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
    table_lines = [join_list(line, '\t') for line in table_lines]
    table_text = join_list(table_lines, '\n')
    write_to_file(table_text, output_file, 'a', '\n')


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
    contigs = import_sequences(genome)
    # extract coding sequences from contigs
    reading_frames = pickle_loader(orf_file)
    genome_info = extract_genome_cds(reading_frames,
                                           contigs, 1)
    # save coding sequences to file
    # create records and write them to file
    cds_lines = create_fasta_lines(genome_info[0], identifier)
    write_lines(cds_lines, cds_file)

    write_protein_table(protein_table, identifier, genome_info[1])

    total_cds = len(genome_info[0])

    return total_cds


def function_helper(input_args):
    """
    """

    results = input_args[-1](*input_args[0:-1])

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

    protein_table = join_paths(temp_directory,
                               ['protein_info_{0}.tsv'.format(index)])

    cds_file = join_paths(temp_directory,
                          ['coding_sequences_{0}.fasta'.format(index)])

    batch_total = 0
    for g in genomes:
        # determine Prodigal ORF file path for current genome
        identifier = file_basename(g, False)
        orf_file_path = join_paths(prodigal_path,
                                   ['{0}_ORF.txt'.format(identifier)])
        total = save_extracted_cds(g, identifier, orf_file_path,
                                   protein_table, cds_file)
        batch_total += total

    return [protein_table, cds_file, batch_total]


def hash_file(file, read_mode):
    """ Creates a hash based on the contents of a file.

        Parameters
        ----------
        file : str
            Path to a file.
        read_mode : str
            File read mode.

        Returns
        -------
        hash_str : str
            Hash computed from file contents.
    """

    with open(file, read_mode) as f:
        hash_obj = hashlib.blake2b()
        file_content = f.read()
        hash_obj.update(file_content)
        hash_str = hash_obj.hexdigest()

    return hash_str


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
    schema_list_file = join_paths(schema_dir, ['.genes_list'])
    pickle_dumper(schema_files, schema_list_file)

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
    pickle_dumper(params, config_file)

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


def user_info(base_url, headers_get):
    """ Get the identifier and role for a user registered
        in a Chewie-NS instance and determine if user has
        permissions based on role.

        Parameters
        ----------
        base_url : str
            The base URL for a Chewie Nomenclature Server
            instance.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        A list with the following elements:
            user_id : str
                User identifier in Chewie-NS.
            user_role : str
                User role in Chewie-NS.
            permission : bool
                True if the user is the Admin or a
                Contributor, False otherwise.
    """

    # verify user role to check permission
    user_info = simple_get_request(base_url, headers_get,
                                   ['user', 'current_user'])
    user_info = user_info.json()

    user_id = str(user_info['id'])
    user_role = user_info['roles'].split(' ')[-1][:-1]
    permission = any(role in user_role for role in ['Admin', 'Contributor'])

    return [user_id, user_role, permission]


def species_ids(species_id, base_url, headers_get):
    """ Gets species identifier or name in a Chewie-NS
        instance.

        Parameters
        ----------
        species_id : str
            The identifier of the species in Chewie-NS or
            the name of the species.
        base_url : str
            The base URL for a Chewie Nomenclature Server
            instance.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        A list with the species identifier and name if the
        species exists in Chewie-NS or a 404 HTTP code if
        the species was not found.
    """

    try:
        int(species_id)
        species_info = simple_get_request(base_url, headers_get,
                                          ['species', species_id])
        if species_info.status_code in [200, 201]:
            species_name = species_info.json()[0]['name']['value']
            return [species_id, species_name]
        else:
            return 404
    except ValueError:
        species_name = species_id
        ns_species = species_list(base_url, headers_get, ['species', 'list'])
        species_id = ns_species.get(species_name, 'not_found')
        if species_id != 'not_found':
            return [species_id, species_name]
        else:
            return 404


def species_list(base_url, headers_get, endpoint_list):
    """ Gets the list of species in a Chewie-NS instance.

        Parameters
        ----------
        base_url : str
            The base URL for a Chewie Nomenclature Server
            instance.
        headers_get : dict
            HTTP headers for GET requests.
        endpoint_list : list
            List with strings to add to tht base_url to
            create the URI for the endpoint.

        Returns
        -------
        species_lst : dict
            Dictionary with species names as keys and
            species identifiers as values.
    """

    res = simple_get_request(base_url, headers_get, endpoint_list)
    res = res.json()
    species_lst = {}
    for sp in res:
        species = sp['name']['value']
        species_url = sp['species']['value']
        species_id = species_url.split('/')[-1]

        species_lst[species] = species_id

    return species_lst


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
        write_lines(lines, failed_file)

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


def extend_list(input_list, *elements):
    """
    """

    input_list.extend(elements)

    return input_list


def is_fasta(file_path):
    """ Checks if a file is a FASTA file.

        Parameters
        ----------
        file_path : str
            Path to the file.

        Returns
        -------
        True if file has valid FASTA format,
        False otherwise.
    """

    with open(file_path, 'r') as handle:
        try:
            fasta = SeqIO.parse(handle, 'fasta')
        except:
            fasta = [False]

        # returns True if FASTA file, False otherwise
        return any(fasta)


def filter_files(files, suffixes, reverse=False):
    """ Filters files names based on a list of suffixes.

        Parameters
        ----------
        files : list
            A list with filenames of file paths.
        suffixes : list
            List with suffixes.
        reverse : bool
            True if files should be filtered out from
            input list or False otherwise.

        Returns
        -------
        filtered : list
            List with files that passed filtering.
    """

    if reverse is False:
        filtered = [file for file in files
                    if any([True for suffix in suffixes if suffix in file])]
    elif reverse is True:
        filtered = [file for file in files
                    if not any([True for suffix in suffixes if suffix in file])]

    return filtered


def filter_non_fasta(files):
    """ Creates a new list of files names/paths that only contains
        FASTA files.

        Parameters
        ----------
        files : list
            A list with files names/paths.

        Returns
        -------
        fasta_files : list
            List with files names/paths that have a FASTA
            format.
    """

    fasta_files = [file for file in files if is_fasta(file) is True]

    return fasta_files


def gene_seqs_info(fasta_file):
    """ Determines the total number of sequences and the mean
        length of sequences in a fasta file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file.

        Returns
        -------
        stats : list
            A list with the path to the FASTA file, the
            total number of records in that file and
            the sequence mean length for the sequences
            in the file.
    """

    seq_generator = SeqIO.parse(fasta_file, 'fasta')
    seqs_lengths = [len(seq) for seq in seq_generator]
    mean_length = sum(seqs_lengths)/len(seqs_lengths)
    total_seqs = len(seqs_lengths)
    stats = [fasta_file, total_seqs, mean_length]

    return stats


def make_blast_db(input_fasta, output_path, db_type):
    """ Creates a BLAST database.

        Parameters
        ----------
        input_fasta : str
            Path to the FASTA file that contains the sequences
            that should be added to the BLAST database.
        output_path : str
            Path to the directory where the database files
            will be created. Database files will have names
            with the path's basemane.
        db_type : str
            Type of the database, nucleotide (nuc) or
            protein (prot).

        Returns
        -------
        Creates a BLAST database with the input sequences.
    """

    blastdb_cmd = ['makeblastdb', '-in', input_fasta, '-out', output_path,
                   '-parse_seqids', '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    makedb_cmd.wait()


def determine_blast_task(sequences):
    """ Determine the type of task that should be used to
        run BLAST (alignments with short sequences require
        definition of different task).

        Parameters
        ----------
        sequences : str
            Path to a file with sequences.

        Returns
        -------
        blast_task : str
            A string that indicates the type of BLAST
            task to run.
    """

    blast_task = 'blastp'
    proteins_lengths = [len(p) for p in sequences]
    minimum_length = min(proteins_lengths)
    if minimum_length < 30:
        blast_task = 'blastp-short'

    return blast_task


def create_directory(directory_path):
    """ Creates a diretory if it does not exist."""

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


def join_paths(parent_path, child_paths):
    """ Creates a new path by joining a parent directory
        and a list with child paths."""

    joined_paths = os.path.join(parent_path, *child_paths)

    return joined_paths


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
        files = filter_files(files, cnst.FASTA_SUFFIXES)
        # get absolute paths
        files = [os.path.join(input_path, file) for file in files]
        # filter any directories that migh end with FASTA extension
        files = [file for file in files if os.path.isdir(file) is False]

        # only keep files whose content is typical of a FASTA file
        fasta_files = filter_non_fasta(files)

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


def escape_special_characters(input_string):
    """ Escapes strings to use in regex.

        Parameters
        ----------
        input_string : str
            String containing characters to escape.

        Returns
        -------
        escaped_string : str
            Escaped string.
    """

    escaped_string = re.escape(input_string)

    return escaped_string


def replace_multiple_characters(input_string):
    """ Replaces multiple characters in a string.

        Parameters
        ----------
        input_string : str
            String with characters to be replaced.

        Returns
        -------
        replaced : str
            Input string without replaced characters.
    """

    replaced = input_string.replace("|", "_")\
                           .replace("_", "-")\
                           .replace("(", "")\
                           .replace(")", "")\
                           .replace("'", "")\
                           .replace("\"", "")\
                           .replace(":", "")

    return replaced


def listdir_fullpath(directory_path):
    """ Gets the full path of the files from a directory.

        Parameters
        ----------
        directory_path : str
            Path to a directory.

        Returns
        -------
        List containing the full path for every file
        in the input directory.
    """

    return [os.path.join(directory_path, f)
            for f in os.listdir(directory_path)]


def flatten_list(list_to_flatten):
    """ Flattens one level of a nested list.

        Parameters
        ----------
        list_to_flatten : list
            List with nested lists.

        Returns
        -------
        flattened_list : str
            Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list


def invert_dictionary(dictionary):
    """ Inverts a dictionary (Keys become values and vice-versa).

        Parameters
        ----------
        dictionary : dict
            Dictionary to be inverted.

        Returns
        -------
        inverted_dict : dict
            Inverted dictionary.
    """

    inverted_dict = {value: key for key, value in dictionary.items()}

    return inverted_dict


def isListEmpty(input_list):
    """ Checks if a nested list is empty. """

    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False


def read_tabular(input_file, delimiter='\t'):
    """ Read tabular file.

        Parameters
        ----------
        input_file : str
            Path to a tabular file.
        delimiter : str
            Delimiter used to separate file fields.

        Returns
        -------
        lines : list
            A list with a sublist per line in the input file.
            Each sublist has the fields that were separated by
            the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def fasta_lines(identifiers, sequences_dictionary):
    """ Creates list with line elements for a FASTA file based
        on the sequence identifiers passed.

        Parameters
        ----------
        identifiers : list
            A list with the identifiers of sequences that will be
            included in the list.
        sequences_dictionary : dict
            A dictionary with sequence identifiers as keys and
            sequences as values.

        Returns
        -------
        seqs_lines : list
            A list with strings representing the header of
            the sequence and the sequence for each of the specified
            sequence identifiers.
    """

    seqs_lines = ['>{0}\n{1}\n'.format(seqid, sequences_dictionary[seqid])
                  for seqid in identifiers]

    return seqs_lines


def write_list(lines, output_file):
    """ Writes list elements to file.

        Parameters
        ----------
        lines : list
            List with the lines that will be written to
            the output file.
        output_file : str
            Path to the output file.

        Returns
        -------
        Writes contents of input list to the output file
        (function does not add any character between lines).
    """

    with open(output_file, 'w') as file:
        file.writelines(lines)


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


def sequences_lengths(fasta_file):
    """ Determines the length of all sequences in a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file with sequences.

        Returns
        -------
        lengths : dict
            Dictionary with the `fasta_file` basename as key and
            a nested dictionary with sequences hashes as keys and
            sequences lengths as values.
    """

    basename = os.path.basename(fasta_file)
    lengths = {basename: {hashlib.sha256(str(rec.seq).encode('utf-8')).hexdigest(): len(rec.seq)
                          for rec in SeqIO.parse(fasta_file, 'fasta')}}

    return lengths


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

    sequences = import_sequences(fasta_path)

    # translate sequences
    translated_seqs = {k: translate_dna(v, table_id, min_len)
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
            the path to the FASTA file with the locus sequences,
            the total number of sequences and sequence mean legth
            for a locus.
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
            List with sublists that contain paths to loci FASTA
            files. Sublists are balanced based on the chosen
            method.
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


# this function should receive an iterable with values per cluster
# possibly convert it so that it can be used for the same as previous function!
def split_blast_inputs_by_core(blast_inputs, cores, blast_files_dir):
    """
    """

    splitted_ids = [[] for cpu in range(cores)]
    cluster_sums = [0 for cpu in range(cores)]
    i = 0
    for cluster in blast_inputs:
        cluster_file = os.path.join(blast_files_dir,
                                    '{0}_ids.txt'.format(cluster))
        with open(cluster_file, 'r') as infile:
            cluster_seqs = [line.strip() for line in infile.readlines()]
        splitted_ids[i].append(cluster)
        cluster_sums[i] += len(cluster_seqs)
        i = cluster_sums.index(min(cluster_sums))

    return splitted_ids


def split_iterable(iterable, size):
    """ Splits a dictionary.

        Parameters
        ----------
        iterable : dict
            Dictionary to split.
        size : int
            Size of dictionaries created from the input
            dictionary.

        Returns
        -------
        chunks : list
            List with dictionaries of defined size
            resulting from splitting the input dictionary.
    """

    chunks = []
    it = iter(iterable)
    for i in range(0, len(iterable), size):
        chunks.append({k: iterable[k] for k in islice(it, size)})

    return chunks


def concatenate_list(str_list, join_char):
    """ Concatenates list elements with specified
        character between each original list element.

        Parameters
        ----------
        str_list : list
            List with strings that will be concatenated.
        join_char : str
            Character that will be used to join list
            elements.

        Returns
        -------
        ids_str : str
            String resulting from the concatenation of
            all strings in the input list.
    """

    concat = join_char.join(str_list)

    return concat


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.

        Returns
        -------
        reverse_complement_strand : str
            The reverse complement of the DNA sequence (lowercase
            is converted to uppercase).
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A',
                       'n': 'N', 'N': 'N'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    # default to 'N' if nucleotide is not in base_complement dictionary
    bases = [base_complement.get(base, 'N') for base in bases]

    complement_strand = ''.join(bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def reverse_str(string):
    """ Reverse character order in input string.

        Parameters
        ----------
        string : str
         String to be reversed.

        Returns
        -------
        revstr : str
            Reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Parameters
        ----------
        dna_str : str
            String representing a DNA sequence.
        table_id : int
            Translation table identifier.

        Returns
        -------
        protseq : Bio.Seq.Seq
            Protein sequence created by translating the
            input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Attempts to translate an input DNA sequence in specified
        orientation and stores exceptions when the input sequence
        cannot be translated.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.
        method : str
            Sequence orientation to attempt translation.
        table_id : int
            Translation table identifier.

        Returns
        -------
        If the sequence can be translated:
            protseq : Bio.Seq.Seq
                Translated DNA sequence.
            myseq : str
                String representing the DNA sequence in the
                orientation used to translate it.
        Otherwise, returns string with the description of the
        exception that was raised.
    """

    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def check_str_alphabet(string, alphabet):
    """ Determine if a string only contains characters from
        specified alphabet.

        Parameters
        ----------
        string : str
            Input string.
        alphabet : str
            String that has all characters from desired
            alphabet.

        Returns
        -------
        "True" if sequence only has characters from specified
        alphabet and string "ambiguous or invalid characters" if
        it any of its characters is not in the alphabet.
    """

    valid_chars = alphabet
    if all(n in valid_chars for n in string) is True:
        return True
    else:
        return 'ambiguous or invalid characters'


def check_str_multiple(string, number):
    """ Determine if length of input string is multiple of
        a specified number.

        Parameters
        ----------
        string : str
            Input string.
        number : int
            Length value should be a multiple of this number.

        Returns
        -------
        "True" if the length of the sequence is a multiple of the
        specified number and "sequence length is not a multiple of number"
        if condition is not satisfied.
    """

    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def translate_dna(dna_sequence, table_id, min_len):
    """ Checks if sequence is valid and attempts to translate
        it, calling several functions to ensure that the sequence
        only has 'ACTG', is multiple of 3 and that it can be
        translated in any of 4 different orientations. Stores
        exceptions so that it is possible to understand why the
        sequence could not be translated.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.
        table_id : int
            Translation table identifier.

        Returns
        -------
        If the sequence can be translated:
            sequence : list
                List with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
            coding_strand : str
                The sequence orientation that codes for the
                protein.
        Otherwise:
            exception_str : str
                A string containing the exceptions that
                explain why the the sequence could not be
                translated.
    """

    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the string is DNA, without ambiguous bases
    valid_dna = check_str_alphabet(original_seq, 'ACTG')
    if valid_dna is not True:
        return valid_dna

    # check if sequence is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return valid_length

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        sequence, exception_collector = retranslate(original_seq,
                                                    translating_methods[i],
                                                    table_id, strands[i],
                                                    exception_collector)

        i += 1
        if i == len(strands) or isinstance(sequence, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(sequence, list):
        return [sequence, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str


def retranslate(sequence, method, table_id, strands, exception_collector):
    """ Sends sequence for translation and collects exceptions when
        the sequence cannot be translated.

        Parameters
        ----------
        sequence : str
            String representing a DNA sequence.
        method : str
            Sequence orientation to attempt translation.
        table_id : int
            Translation table identifier.
        strands : list
            List with maximum of 4 different orientations
            to attempt translation, 'original', 'revcomp',
            'rev' and 'revrevcomp'.
        exception_collector : list
            List used to store all exceptions arising from
            translation attempts.

        Returns
        -------
        If the sequence can be translated:
            translated_seq : list
                List with the protein sequence and with the
                DNA sequence in the orientation used for translation.
            exception_collector : list
                List with the exceptions that were captured when the
                sequence could not be translated.
        Otherwise:
            translated_seq : str
                String with the description of the last exception
                captured because sequence could not be translated.
            exception_collector : list
                List with all exception that have been captured
                for translation attempts.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def is_url(url):
    """ Checks if a URL is valid.

        Parameters
        ----------
        url : str
            The url to be checked.

        Returns
        -------
        True if url is valid, False otherwise.
    """

    try:
        result = urlparse(url)
        return all([result.scheme, result.netloc, result.path])
    except:
        return False


def make_url(base_url, *res, **params):
    """ Creates a URL.

        Parameters
        ----------
        base_url : str
            The base URL.
        res : str
            Endpoint(s) to add to the base URL.
        params : str
            Addtional parameters (WIP).

        Returns
        -------
        url : str
            URL constructed by adding endpoints and
            additional parameters. If the provided
            base URL is invalid returns 'An invalid
            URL was provided.'
    """

    url = base_url
    # Check if the url is valid
    if is_url(url):

        if url[-1] == '/':
            url = url[:-1]

        # Add the endpoints
        for r in res:
            url = '{0}/{1}'.format(url, r)

        # Add params if they are provided
        if params:
            url = '{0}?{1}'.format(url, urlencode(params))

        return url
    else:
        return 'An invalid URL was provided.'


def get_sequence_from_url(url):
    """ Extracts sequence from URL.

        Parameters
        ----------
        url : str
            URL that contains the sequence.

        Returns
        -------
        seq : str
            Sequence extracted from the URL.
    """

    seq = parse_qs(urlsplit(url).query)['sequence'][0]

    return seq


def login_user_to_NS(server_url, email, password):
    """ Login user to a Chewie-NS instance.

        Parameters
        ----------
        server_url : str
            Base URL of the Chewie-NS instance.
        email : str
            Email linked to user's account in Chewie-NS.
        password : str
            Password for the user's account in Chewie-NS.

        Returns
        -------
        If Login is successful:
            token : str
                Authorization token to perform requests to
                Chewie-NS.
        Otherwise, returns False.
    """

    auth_params = {}
    auth_params['email'] = email
    auth_params['password'] = password

    auth_headers = {}
    auth_headers['Content-Type'] = 'application/json'
    auth_headers['accepts'] = 'application/json'

    auth_url = make_url(server_url, 'auth', 'login')

    auth_r = requests.post(auth_url, data=json.dumps(auth_params),
                           headers=auth_headers, verify=False)

    auth_result = auth_r.json()
    if auth_result['status'] == 'success':
        token = auth_result['access_token']
    else:
        token = False

    return token


def capture_login_credentials(base_url):
    """ Captures login credentials and logs user into
        Chewie-NS to get authorization token.

        Parameters
        ----------
        base_url : str
            Base URL of the Chewie-NS instance.

        Returns
        -------
        token : str
            Authorization token to perform requests to
            Chewie-NS.
    """

    print('\nPlease provide login credentials:')
    user = input('USERNAME: ')
    password = getpass('PASSWORD: ')
    print()
    # get token
    token = login_user_to_NS(base_url, user, password)
    # if login was not successful, stop the program
    if token is False:
        sys.exit('Invalid credentials.')

    return token


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
        configs = pickle_loader(config_file)
    else:
        sys.exit('Could not find a valid config file.')

    return configs


def get_species_schemas(schema_id, species_id, base_url, headers_get):
    """ Determines if a species in Chewie-NS has a schema
        with specified identifier.

        Parameters
        ----------
        schema_id : str
            The identifier of the schema in the Chewie-NS.
        species_id : str
            The identifier of the schema's species in the
            Chewie-NS.
        base_url : str
            Base URL of the Chewie Nomenclature server.
        headers_get : dict
            HTTP headers for GET requests.

        Returns
        -------
        A list with the schema identifier, the schema URI,
        the schema name and a dictionary with all properties
        for the schema in Chewie-NS.

        Raises
        ------
        SystemExit
            - If the schema with the specified identifier does
              not exist.
            - If the process cannot retrieve the list of schemas
              for the species.
    """

    # get the list of schemas for the species
    schema_get = simple_get_request(base_url, headers_get,
                                    ['species', species_id, 'schemas'])
    schema_get_status = schema_get.status_code
    if schema_get_status in [200, 201]:
        species_schemas = schema_get.json()

        # extract schemas identifiers, URIs and names from response
        schemas_info = []
        for s in species_schemas:
            schema_uri = s['schemas']['value']
            scid = schema_uri.split('/')[-1]
            schema_name = s['name']['value']
            schemas_info.append([scid, schema_uri, schema_name])

        # select schema with specified identifier
        schema = [s for s in schemas_info if schema_id in s]
        if len(schema) > 0:
            # get schema parameters
            schema = schema[0]
            schema_params = requests.get(schema[1], headers=headers_get,
                                         verify=False)
            schema_params = schema_params.json()[0]
            schema.append(schema_params)
            return schema
        else:
            sys.exit('Could not find a schema with provided identifier.')
    else:
        sys.exit('Could not retrieve schemas for current species.')


def upload_file(file, filename, url, headers, verify_ssl):
    """ Uploads a file to the NS.

        Parameters
        ----------
        file : str
            Path to the file to upload.
        filename : str
            Name used to save the file in the NS.
        url : str
            Endpoint URL that receives the POST request.
        headers : dict
            HTTP POST request headers.
        verify_sll : bool
            If the SSL certificates should be verified in
            HTTPS requests (False for no verification, True
            otherwise).

        Returns
        -------
        response : requests.models.Response
            Response object from the 'requests' module.
    """

    file_handle = open(file, 'rb')
    files = {'file': (filename, file_handle)}
    response = requests.post(url,
                             headers=headers,
                             files=files,
                             verify=verify_ssl)
    file_handle.close()

    return response


def upload_data(data, url, headers, verify_ssl):
    """ Uploads data to Chewie-NS.

        Parameters
        ----------
        data
            The data that will be sent to Chewie-NS (any data
            type accepted by requests 'data' argument).
        url : str
            Endpoint URL that receives the POST request.
        headers : dict
            HTTP POST request headers.
        verify_sll : bool
            If the SSL certificates should be verified in
            HTTPS requests (False for no verification, True otherwise).

        Returns
        -------
        response : requests.models.Response
            Response object from the 'requests' module.
    """

    response = requests.post(url,
                             headers=headers,
                             data=data,
                             verify=verify_ssl)

    return response


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


def input_timeout(prompt, timeout):
    """ Adds timeout feature when requesting user input.

        Parameters
        ----------
        prompt : str
            Message to print to stdout to request for user
            input.
        timeout : int
            Maximum number of seconds that the process will
            wait for input.

        Returns
        -------
        String with user input.

        Raises
        ------
        SystemExit
            - If there is no user input before timeout.
    """

    pool = ThreadPool(processes=1)
    answer = pool.apply_async(input, args=[prompt])

    try:
        return answer.get(timeout=timeout)
    except TimeoutError as e:
        sys.exit('Timed out.')


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

        translation = translate_dna(sequence, translation_table, minimum_length)
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

            dna_lines = join_list(dna_lines, '\n')
            write_to_file(dna_lines, dna_file, 'a', '\n')
            dna_lines = []

            prot_lines = join_list(prot_lines, '\n')
            write_to_file(prot_lines, protein_file, 'a', '\n')
            prot_lines = []

    return [invalid_alleles, total_seqs]


def hash_sequence(string):
    """ Compute SHA256 for an input string.

        Parameters
        ----------
        string : str
            Input string to hash.

        Returns
        -------
        sha256 : str
            String representation of the sha256 HASH object.
    """

    sha256 = hashlib.sha256(string.encode('utf-8')).hexdigest()

    return sha256


def fasta_str_record(seqid, sequence):
    """ Creates the string representation of a FASTA record.

        Parameters
        ----------
        seqid : str
            Sequence identifier to include in the header.
        sequence : str
            String representing DNA or Protein sequence.

        Returns
        -------
        record : str
            String representation of the FASTA record.
    """

    record = '>{0}\n{1}'.format(seqid, sequence)

    return record


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
            seq_hash = hash_sequence(sequence)

            # store only the hash for distinct sequences
            if seq_hash not in seqs_dict:
                seqs_dict[seq_hash] = seqid
                recout = fasta_str_record(seqid, sequence)
                out_seqs.append(recout)
            elif seq_hash in seqs_dict:
                total += 1
        else:
            exausted = True

        if len(out_seqs) == out_limit or exausted is True:
            if len(out_seqs) > 0:
                out_seqs = join_list(out_seqs, '\n')
                write_to_file(out_seqs, unique_fasta, 'a', '\n')
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


def get_sequences_by_id(sequences_index, seqids, out_file, limit=5000):
    """ Retrieves sequences from an indexed FASTA file.

        Parameters
        ----------
        sequences_index : Bio.File._IndexedSeqFileDict
            Fasta file index.
        seqids : list
            List with the identifiers of the sequences
            that should be retrieved.
        out_file : str
            Path to the FASTA file to which selected
            sequences will be saved.
        limit : int
            Maximum number of sequences that will be
            kept in memory at a time (to avoid keeping
            huge datasets in memory).

        Returns
        -------
        Creates a file with the sequences that have the
        identifiers in the input list.
    """

    records = []
    for seqid in seqids:
        sequence = str(sequences_index[seqid].seq)
        record = fasta_str_record(seqid, sequence)
        records.append(record)

        if len(records) == limit or seqid == seqids[-1]:
            lines = join_list(records, '\n')
            write_to_file(lines, out_file, 'a', '\n')
            records = []


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

    short_path = join_paths(schema_dir, ['short'])
    create_directory(short_path)

    for file in schema_files:
        short_file = join_paths(short_path, [file_basename(file)])
        short_file = short_file.replace('.fasta', '_short.fasta')
        shutil.copy(file, short_file)

    return True


def split_fasta(fasta_path, output_path, num_seqs, filenames):
    """ Splits a FASTA file.

        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.
        output_path : str
            Path to the output directory where new FASTA
            files will be created.
        num_seqs : int
            Split FASTA file into files with this number
            of sequences.
        filenames : list
            List with names to attribute to new files.

        Returns
        -------
        splitted_files : list
            List with paths to the new files that were
            created by splitting the input FASTA file.
    """

    splitted_files = []
    current_recs = []
    records = [rec for rec in SeqIO.parse(fasta_path, 'fasta')]
    for record in records:
        current_recs.append(record)
        if len(current_recs) == num_seqs or record.id == records[-1].id:
            file_name = filenames.__next__()
            file_name = replace_multiple_characters(file_name)

            new_file = join_paths(output_path,
                                  ['{0}{1}'.format(file_name, '.fasta')])

            splitted_files.append(new_file)

            write_records(current_recs, new_file)

            current_recs = []

    return splitted_files


def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None):
    """ Execute BLAST to align sequences in a FASTA file
        against a BLAST database.

        Parameters
        ----------
        blast_path : str
            Path to BLAST executables.
        blast_db : str
            Path to the BLAST database.
        fasta_file : str
            Path to the FASTA file with sequences to
            align against the BLAST database.
        blast_output : str
            Path to the file that will be created to
            store BLAST results.
        max_hsps : int
            Maximum number of High Scoring Pairs per
            pair of aligned sequences.
        threads : int
            Number of threads/cores used to run BLAST.
        ids_file : str
            Path to a file with sequence identifiers,
            one per line. Sequences will only be aligned
            to the sequences in the BLAST database that
            have any of the identifiers in this file.

        Returns
        -------
        stderr : str
            String with errors raised during BLAST execution.
    """

    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid score',
                  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
                  '-evalue', '0.001']

    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])

    blast_proc = subprocess.Popen(blast_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()

    return stderr


def cluster_blaster(inputs):
    """ Aligns sequences in the same cluster with BLAST.

        Parameters
        ----------
        inputs : list
            List with clusters identifiers, the path to
            BLAST executables, the path to a BLAST database,
            the path to the directory that contains files with
            the identifiers of the sequences in each cluster
            and where FASTA files with the sequences in each
            cluster and files with BLAST results will be
            written to and the path to the FASTA file with
            the protein sequences in all clusters.

        Returns
        -------
        out_files : list
            List with the paths to the files with BLAST
            results for each cluster.
    """

    blast_inputs = inputs[0:-4]
    blastp_path = inputs[-4]
    output_directory = inputs[-2]
    proteins_file = inputs[-1]
    blast_db = inputs[-3]

    indexed_fasta = SeqIO.index(proteins_file, 'fasta')

    out_files = []
    try:
        for cluster in blast_inputs:

            cluster_id = cluster
            ids_file = os.path.join(output_directory,
                                    '{0}_ids.txt'.format(cluster_id))

            with open(ids_file, 'r') as clstr:
                cluster_ids = [l.strip() for l in clstr.readlines()]

            fasta_file = os.path.join(output_directory,
                                      '{0}_protein.fasta'.format(cluster_id))
            # create file with protein sequences
            get_sequences_by_id(indexed_fasta, cluster_ids, fasta_file)

            blast_output = os.path.join(output_directory,
                                        '{0}_blast_out.tsv'.format(cluster_id))

            # Use subprocess to capture errors and warnings
            stderr = run_blast(blastp_path, blast_db, fasta_file, blast_output,
                               1, 1, ids_file)

            out_files.append(blast_output)
    except Exception as e:
        print(e)

    return out_files


def blast_inputs(clusters, output_directory, ids_dict):
    """ Creates files with the identifiers of the sequences
        in each cluster.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of cluster
            representatives as keys and a list with tuples
            as values (each tuple has the identifier of a
            sequence that is in the cluster, the percentage
            of shared minimizers and the length os that
            sequence).
        output_directory : str
            Path to the directory where files with identifiers
            will be created.
        ids_dict : dict
            Dictionary that maps sequence identifiers to
            shorter and unique identifiers that will be
            saved in the files and used as sequence
            identifiers during BLAST to avoid errors
            related with sequence headers/identifiers
            that exceed length limit allowed by BLAST.

        Returns
        -------
        ids_to_blast : list
            List with the identifiers of all clusters.
    """

    rev_ids = {v: k for k, v in ids_dict.items()}

    ids_to_blast = []
    for i in clusters:

        cluster_file = os.path.join(output_directory,
                                    '{0}_ids.txt'.format(rev_ids[i]))
        cluster_ids = [rev_ids[i]] + [rev_ids[seq[0]] for seq in clusters[i]]
        cluster_lines = join_list(cluster_ids, '\n')
        write_to_file(cluster_lines, cluster_file, 'w', '')
        ids_to_blast.append(rev_ids[i])

    return ids_to_blast


def divide_list_into_n_chunks(list_to_divide, n):
    """ Divides a list into a defined number of sublists.

        Parameters
        ----------
        list_to_divide : list
            List to divide into sublists.
        n : int
            Number of sublists to create.

        Returns
        -------
        sublists : list
            List with the sublists created by dividing
            the input list.
    """

    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    # exclude lists that are empty due to small number of elements
    sublists = [i for i in sublists if len(i) > 0]

    return sublists


def apply_bsr(inputs):
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

    ids_dict = inputs[3]
    blast_file = inputs[0]
    blast_results = read_tabular(blast_file)
    self_scores = {r[0]: r[2] for r in blast_results if r[0] == r[1]}
    # do not include self-scores lines, no need to evaluate those hits
    blast_results = [r for r in blast_results if r[0] != r[1]]

    fasta_file = inputs[1]
    lengths = {}
    for k in self_scores:
        record = fasta_file.get(ids_dict[k])
        sequence = str(record.seq)
        lengths[k] = len(sequence)

    bsr = inputs[2]

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


def sequence_kmerizer(sequence, k_value, offset=1, position=False):
    """ Decomposes a sequence into kmers.

        Parameters
        ----------
        sequence : str
            Sequence to divide into kmers.
        k_value : int
            Value for the size k of kmers.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        kmers : list
            List with the kmers determined for the input
            sequence. The list will contain strings if
            it is not specified that positions should be
            stored and tuples of kmer and start position
            if the position is stored.
    """

    if position is False:
        kmers = [sequence[i:i+k_value]
                 for i in range(0, len(sequence)-k_value+1, offset)]
    elif position is True:
        kmers = [(sequence[i:i+k_value], i)
                 for i in range(0, len(sequence)-k_value+1, offset)]

    return kmers


def determine_minimizers(sequence, adjacent_kmers, k_value, position=False):
    """ Determine the minimizers for a sequence based on
        lexicographical order. Skips windows that
        cannot have a minimizer based on the minimizer
        computed in the previous iteration.

        Parameters
        ----------
        sequence : str
            String representing the sequence.
        adjacent_kmers : int
            Window size value. Number of adjacent kmers per group.
        k_value : int
            Value of k for kmer size.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        minimizers : list
            A list with the set of minimizers determined
            for the input sequence.
    """

    # break sequence into kmers
    kmers = sequence_kmerizer(sequence, k_value, position=position)

    i = 0
    previous = None
    sell = False
    minimizers = []
    # determine total number of windows
    last_window = (len(kmers)-adjacent_kmers)
    while i <= last_window:
        # get kmers in current window
        window = kmers[i:i+adjacent_kmers]
        # pick smallest kmer as minimizer
        minimizer = [min(window)]
        # get position in window of smallest minimizer
        minimizer_idx = window.index(minimizer[0])
        # sliding window that does not included last minimizer
        if previous is None:
            # simply store smallest minimizer
            minimizers.extend(minimizer)
        # sliding window includes last minimizer because we
        # skipped some sliding windows
        else:
            # check if minimizer is the same as the one picked
            # in the last window
            # Do not store minimizer if it is the same
            if minimizer[0] != previous:
                # get kmers smaller than last minimizer
                skipped = window[1:minimizer_idx]
                # determine if any of the smaller kmers is
                # the minimizer of a skipped window
                minimal = previous
                for m in skipped:
                    if m < minimal:
                        minimizer.append(m)
                        minimal = m
                minimizers.extend(minimizer)

        # slide by 1 if minimizer has index 0 in window
        if minimizer_idx == 0:
            i += 1
            previous = None
        # skip sliding windows based on minimizer position
        else:
            i += minimizer_idx
            # if adding minimizer index surpasses last window value we
            # might miss one last minimizer because it will fail the condition
            # find a better way to control this condition!
            if i > last_window and sell is False:
                i = last_window
                sell = True
            previous = minimizer[0]

    return minimizers


# for AlleleCall this function needs to receive the variable 'reps_groups'
# and to accept an argument that controls if it is allowed to create new
# clusters based on new genes that it finds.
def cluster_sequences(sorted_sequences, word_size=4, clustering_sim=0.2,
                      mode='greedy', representatives=None, grow=True,
                      offset=1, minimizer=True):
    """ Cluster sequences based on shared percentage of kmers/minimizers.

        Parameters
        ----------
        sorted_sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values. Sorted by decreasing sequence
            length.
        word_size : int
            Value k for the kmer size.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.
        mode : str
            Clustering mode.
        representatives : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values.
        grow : bool
            If it is allowed to create new clusters.
        offset : int
            Value to indicate offset of consecutive kmers.
        minimizer : bool
            If clustering should be based on shared minimizers.

        Returns
        -------
        A list with the following elements:
            clusters : dict
                Dictionary with the identifiers of sequences
                that are clusters representatives as keys and
                a list with tuples as values. Each tuple has
                the identifier of a sequence that was added to
                the cluster, the percentage of shared
                kmers/minimizers and the length of the clustered
                sequence.
            reps_sequences : dict
                Dictionary with the identifiers of sequences
                that are clusters representatives as keys and
                their sequences as values.
    """

    print(word_size)
    clusters = {}
    reps_sequences = {}
    if representatives is None:
        reps_groups = {}
    else:
        reps_groups = representatives

    # repetitive = 0
    for protid, protein in sorted_sequences.items():
        if minimizer is True:
            minimizers = determine_minimizers(protein, word_size,
                                              word_size, position=False)
            kmers = set(minimizers)
            # if len(kmers) < (0.98*len(minimizers)):
            #     repetitive += 1
        # check if set of distinct kmers is much smaller than the
        # set of minimizers to understand if sequence has too much redundancy
        elif minimizer is False:
            kmers = sequence_kmerizer(protein, word_size, offset, False)

        current_reps = [reps_groups[k] for k in kmers if k in reps_groups]
        current_reps = flatten_list(current_reps)

        # count number of kmer hits per representative
        counts = Counter(current_reps)
        selected_reps = [(k, v/len(kmers))
                         for k, v in counts.items()
                         if v/len(kmers) >= clustering_sim]

        # sort to get most similar at index 0
        if len(selected_reps) > 0:
            for s in selected_reps:
                clusters[s[0]].append((protid, s[1], len(protein)))
        else:
            for k in kmers:
                reps_groups.setdefault(k, []).append(protid)

            clusters[protid] = [(protid, 1.0, len(protein))]
            reps_sequences[protid] = protein

    # print(repetitive)
    return [clusters, reps_sequences]


def write_clusters(clusters, outfile):
    """ Writes clusters to file.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        outfile : str
            Path to the file that will be
            created to save clusters.
    """

    cluster_num = 0
    cluster_lines = []
    for rep, seqids in clusters.items():
        cluster_lines.append('>Cluster_{0}'.format(cluster_num))
        clustered = ['\t{0}, {1}, {2}'.format(s[0], s[1], s[2])
                     for s in seqids]
        cluster_lines.extend(clustered)
        cluster_num += 1
    cluster_text = join_list(cluster_lines, '\n')

    write_to_file(cluster_text, outfile, 'w', '\n')


def representative_prunner(clusters, sim_cutoff):
    """ Removes sequences from clusters based on a similarity
        threshold.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        sim_cutoff : float
            Similarity threshold value. Sequences with
            equal or greater similarity value are excluded
            from clusters.

        Returns
        -------
        A list with the following elements:
            prunned_clusters : dict
                Input dictionary without values for the
                sequences that had similarity values
                equal or greater that defined threshold.
            excluded : list
                List with a list per cluster. Each list
                has the sequences that were excluded from
                a cluster.
    """

    excluded = []
    prunned_clusters = {}
    for rep, seqids in clusters.items():
        # get high scoring elements
        # determine if distribution is multimodal
        # determine if rep length is considerably larger than first high scoring element

        # this removes the representative from the cluster
        # rep_info = clusters[rep][0]
        # length_cutoff = rep_info[2] - (rep_info[2]*0.2)
        # keep = []
        # remove = []
        # for s in seqids:
        #     query_sim = s[1]
        #     query_len = s[2]
        #     if query_sim < sim_cutoff:
        #         keep.append(s)
        #     elif query_sim >= sim_cutoff and query_len < length_cutoff:
        #         keep.append(s)
        #     elif query_sim >= sim_cutoff and query_len > length_cutoff:
        #         if s[0] != rep:
        #             remove.append(s)

        # prunned_clusters[rep] = keep
        # excluded.extend(remove)

        prunned_clusters[rep] = [seqid
                                 for seqid in seqids
                                 if seqid[1] < sim_cutoff]
        excluded.extend([seqid
                         for seqid in seqids
                         if seqid[1] >= sim_cutoff and seqid[0] != rep])

    return [prunned_clusters, excluded]


def determine_singletons(clusters):
    """ Determines clusters that only have one sequence.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.

        Returns
        -------
        singletons : dict
            Dictionary entries for clusters that only have
            one sequence.
    """

    singletons = {k: v for k, v in clusters.items() if len(v) == 0}

    return singletons


def remove_entries(dictionary, keys):
    """ Creates new dictionary without entries with
        specified keys.

        Parameters
        ----------
        dictionary : dict
            Input dictionary.
        keys : list
            List of keys for the entries that should
            not be included in the new dictionary.

        Returns
        -------
        new_dict : dict
            Dictionary without entries with keys in
            the input list.
    """

    new_dict = {k: v for k, v in dictionary.items() if k not in keys}

    return new_dict


def intra_cluster_sim(clusters, protein_file, word_size, intra_filter):
    """ Determines the percentage of shared kmers/minimizers
        between sequences in the same cluster and excludes
        sequences from a clusters sequences that are similar
        to other sequences in the cluster.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        protein_file : str
            Path to a FASTA file with protein sequences.
        word_size : int
            Value k for the kmer size.
        intra_filter : float
            Similarity threshold value. If two sequences in
            the same cluster have a similarity value equal
            or greater to this value, one of the sequences
            will be excluded from the cluster.

        Returns
        -------
        excluded_dict : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
            lists as values. Each list has two elements:
            a list with the identifiers of the sequences
            that were excluded from the clusters and a
            dictionary with sequences identifiers as keys
            and tuples with sequence identifiers and the
            similarity value for each match that led to
            an eclusion.
    """

    excluded_dict = {}
    for k, v in clusters.items():
        # get identifiers of sequences in the cluster
        cluster_ids = v

        # get protein sequences
        cluster_sequences = {}
        for seqid in cluster_ids:
            # get sequences through indexed FASTA file
            cluster_sequences[seqid[0]] = str(protein_file[seqid[0]].seq)

        # get all kmers per sequence
        kmers_mapping = {}
        cluster_kmers = {}
        for seqid, prot in cluster_sequences.items():
            minimizers = determine_minimizers(prot, word_size, word_size)
            kmers = set([m[0] for m in minimizers])
            # dict with sequence indentifiers and kmers
            cluster_kmers[seqid] = kmers

            # create dict with kmers as keys and list
            # of sequences with given kmers as values
            for kmer in kmers:
                kmers_mapping.setdefault(kmer, []).append(seqid)

        sims_cases = {}
        excluded = []
        # for each sequence in the cluster
        for seqid, kmers in cluster_kmers.items():
            if seqid not in excluded:
                query_kmers = kmers
                # determine sequences that also have the same kmers
                current_reps = [kmers_mapping[kmer] for kmer in query_kmers if kmer in kmers_mapping]
                current_reps = flatten_list(current_reps)

                # count number of common kmers with other sequences
                counts = Counter(current_reps)
                # determine sequences that are equal or above a similarty threshold
                current_reps = [(s, v/len(kmers)) for s, v in counts.items() if v/len(kmers) >= intra_filter]
                # sort to get most similar first
                sims = sorted(current_reps, key=lambda x: x[1], reverse=True)

                if len(sims) > 1:
                    # exclude current sequence
                    candidates = [s for s in sims if s[0] != seqid]
                    for c in candidates:
                        # if query sequence is longer, keep it and exclude candidate
                        if len(query_kmers) >= len(cluster_kmers[c[0]]):
                            sims_cases[c[0]] = (seqid, c[1])
                            excluded.append(c[0])
                        # otherwise, exclude query sequence
                        elif len(cluster_kmers[c[0]]) > len(query_kmers):
                            sims_cases[seqid] = (c[0], c[1])
                            excluded.append(seqid)

        # convert into set first to remove possible duplicates
        excluded_dict[k] = [list(set(excluded)), sims_cases]

    return excluded_dict


def get_datetime():
    """ Returns datetime module object with
        information about current date and hour.
    """

    current_datetime = dt.datetime.now()

    return current_datetime


def datetime_str(datetime_obj, date_format='%Y-%m-%dT%H:%M:%S'):
    """ Converts datetime module object to formatted string.

        Parameters
        ----------
        datetime_obj : datetime.datetime
            Datetime module object.
        date_format : str
            Format for string representation of the date
            object.

        Returns
        -------
        dt_str : str
            String representation of the date according
            to specified format.
    """

    dt_str = dt.datetime.strftime(datetime_obj, date_format)

    return dt_str


def datetime_diff(sdate, edate):
    """ Returns the difference in minutes and the
        remaining in seconds between two dates.

        Parameters
        ----------
        sdate : datetime.datetime
            Datetime module object corresponding to
            the oldest date.
        edate : datetime.datetime
            Datetime module object corresponding to
            the most recent date.

        Returns
        -------
        A list with the following elements:
            minutes : float
                Time difference between the two dates
                in minutes.
            seconds : float
                The remaining time difference in seconds.
    """

    delta = edate - sdate
    minutes, seconds = divmod(delta.total_seconds(), 60)

    return [minutes, seconds]


def delete_directory(directory_path):
    """ Deletes a directory.
    """

    shutil.rmtree(directory_path)


def validate_date(date):
    """ 
    """

    valid = False
    try:
        date = dt.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
        valid = date
    except ValueError:
        date = dt.datetime.strptime(date+'.0', '%Y-%m-%dT%H:%M:%S.%f')
        valid = date

    return valid
