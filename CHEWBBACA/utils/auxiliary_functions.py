#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHORS

    Mickael Silva
    github: @

    Pedro Cerqueira
    github: @pedrorvc

    Rafael Mamede
    github: @rfm-targa

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
import threading
import itertools
import datetime as dt
import multiprocessing
import concurrent.futures
from getpass import getpass
from collections import Counter
from multiprocessing import TimeoutError
from multiprocessing.pool import ThreadPool
from SPARQLWrapper import SPARQLWrapper, JSON
from urllib.parse import urlparse, urlencode, urlsplit, parse_qs

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

try:
    from utils import constants as cnst
    from createschema import new_utils
    from createschema import init_schema_4_bbaca
except:
    from CHEWBBACA.utils import constants as cnst
    from CHEWBBACA.createschema import new_utils
    from CHEWBBACA.createschema import init_schema_4_bbaca

UNIPROT_SERVER = SPARQLWrapper("http://sparql.uniprot.org/sparql")


def join_list(lst, link):
    """
    """

    return link.join(lst)


def file_basename(file_path, suffix=True):
    """ Get file basename without extension suffix.

        Args:
            file_path (str): the path to the file.
            Full path or relative path.

        Returns:
            file (str): the file basename without the
            file extension suffix.
    """

    file = os.path.basename(file_path)

    if suffix is False:
        file = file.split('.')[0]

    return file


def check_connection(ns_url, headers=cnst.HEADERS_GET_JSON):
    """
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


def pickle_dumper(pickle_out, content):
    """
    """

    with open(pickle_out, 'wb') as po:
        pickle.dump(content, po)


def pickle_loader(pickle_in):
    """
    """

    with open(pickle_in, 'rb') as pi:
        data = pickle.load(pi)

    return data


def file_zipper(ori_file, zip_file):
    """
    """

    with zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(ori_file, os.path.basename(ori_file))

    return zip_file


def remove_files(files):
    """
    """
        
    for f in files:
        os.remove(f)


def count_sequences(fasta_file):
    """
    """

    records = SeqIO.parse(fasta_file, 'fasta')
    total_seqs = len(list(records))

    return total_seqs


def simple_get_request(base_url, headers, endpoint_list):
    """ Constructs an endpoint URI and uses a GET method to retrive
        information from the endpoint.

        Args:
            base_url (str): the base URI for the NS, used to concatenate
            with a list of elements and obtain endpoints URL.
            headers (dict): headers for the GET method used to
            get data from the API endpoints.
            endpoint_list (list): list with elements that will be
            concatenated to the base URL to obtain the URL for
            the API endpoint.
        Returns:
            res (requests.models.Response): response object from
            the GET method.
    """

    # unpack list of sequential endpoints and pass to create URI
    url = make_url(base_url, *endpoint_list)

    res = requests.get(url, headers=headers, timeout=30, verify=False)

    return res


def simple_post_request(base_url, headers, endpoint_list, data):
    """ Constructs an endpoint URI and uses a POST method to insert
        information into the NS structure.

        Args:
            base_url (str): the base URL for the NS, used to concatenate
            with a list of elements and obtain endpoints URL.
            headers (dict): headers for the POST method used to
            insert data into the NS.
            endpoint_list (list): list with elements that will be
            concatenated to the base URL to obtain the URL for
            the API endpoint.
        Returns:
            res (requests.models.Response): response object from
            the POST method.
    """

    # unpack list of sequential endpoints and pass to create URI
    url = make_url(base_url, *endpoint_list)
    res = requests.post(url, data=json.dumps(data), headers=headers, verify=False)

    return res


def concatenate_files(files, output_file, header=None):
    """
    """

    with open(output_file, 'w') as of:
        if header is not None:
            of.write(header)
        for f in files:
            with open(f, 'r') as fd:
                shutil.copyfileobj(fd, of)

    return output_file


def write_to_file(text, output_file, write_mode, end_char):
    """
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def write_fasta(fasta_lines, output_file):
    """
    """

    joined_lines = join_list(fasta_lines, '\n')

    write_to_file(joined_lines, output_file, 'a', '\n')


def integer_headers(input_fasta, output_fasta):
    """
    """

    start = 1
    ids_dict = {}
    seqs = []
    limit = 5000
    for rec in SeqIO.parse(input_fasta, 'fasta'):
        new_id = 'seq_{0}'.format(start)
        ids_dict[new_id] = rec.id
        sequence = str(rec.seq)
        new_rec = '>{0}\n{1}'.format(new_id, sequence)
        seqs.append(new_rec)
        if len(seqs) == limit:
            write_fasta(seqs, output_fasta)
            seqs = []
        start += 1

    if len(seqs) > 0:
        write_fasta(seqs, output_fasta)

    return ids_dict


def create_fasta_lines(sequences, genome_id):
    """

        Args:

        Returns:

        Example:
    """

    lines = []
    for seqid in sequences:
        sequence = sequences[seqid]
        record = '>{0}-protein{1}\n{2}'.format(genome_id, seqid, sequence)

        lines.append(record)

    return lines


def import_sequences(fasta_path):
    """ Imports sequences from a FASTA file.

        Args:
            fasta_path (str): full path to the FASTA file.

        Returns:
            dictionary that has sequences ids as keys and DNA
            sequences as values.
    """

    records = SeqIO.parse(fasta_path, 'fasta', generic_dna)
    seqs_dict = {rec.id: str(rec.seq.upper()) for rec in records}

    return seqs_dict


def extract_coding_sequences(reading_frames, contigs, starting_id):
    """ Extracts CDSs from contigs based on the start
        and stop codon positions determined by Prodigal.

        Args:
            reading_frames (str): full path to the ORF file derived
            from Prodigal.
            contigs (dict): a dictionary with contigs ids as keys
            and contigs sequences as values.
            starting_id (int): integer identifier to give to the
            first CDS extracted and that will be incremented to
            serve as identifier for subsequent CDSs.

        Returns:
            coding_sequences (dict): dictionary with CDSs ids as keys and
            CDSs DNA sequences as values.
            coding_sequences_info (list): list with a sublist for each
            extracted CDS. Sublists have information about the extracted CDS
            (identifier of the contig where the CDS was found, start position
            in the contig, stop position in the contig, sequence identifier
            attributed to that CDS and the strand that coded for that CDS.)
    """

    # load binary file with a list of lists for each contig
    # each sublist has a start codon and stop codon positions
    # in the contig
    with open(reading_frames, 'rb') as orf_file:
        rfs = pickle.load(orf_file)

    seqid = starting_id
    coding_sequences = {}
    coding_sequences_info = []
    for contig_id, frames in rfs.items():
        # for each start and stop codon in the contig
        for cds in frames:
            start_codon = cds[0]
            stop_codon = cds[1]
            strand = cds[2]
            # extract CDS sequence
            cds_sequence = contigs[contig_id][start_codon:stop_codon].upper()
            # check coding strand to change sequence orientation
            # if needed
            if strand == 0:
                cds_sequence = reverse_complement(cds_sequence)

            # store CDS with unique id
            coding_sequences[seqid] = cds_sequence

            # store CDS information
            coding_sequences_info.append([contig_id, str(start_codon),
                                          str(stop_codon), str(seqid),
                                          str(strand)])

            # increment seqid
            seqid += 1

    return [coding_sequences, coding_sequences_info]


def write_protein_table(file_name, genome_id, cds_info):
    """
    """

    table_lines = [[genome_id] + protein_info
                   for protein_info in cds_info]
    table_lines = [join_list(line, '\t') for line in table_lines]
    table_text = join_list(table_lines, '\n')
    write_to_file(table_text, file_name, 'a', '\n')


def process_cdss(input_data):
    """
    """

    genomes = input_data[0:-3]
    prodigal_path = input_data[-3]
    parent_directory = input_data[-2]
    index = input_data[-1]
    temp_directory = os.path.join(parent_directory, 'temp')
    protein_table = os.path.join(parent_directory,
                                 'protein_info_{0}.tsv'.format(index))

    cds_file = os.path.join(temp_directory,
                            'coding_sequences_{0}.fasta'.format(index))
    protid = 1
    for i, g in enumerate(genomes):
        try:
            # determine Prodigal ORF file path for current genome
            identifier = os.path.basename(g).split('.')[0]
            orf_file_path = os.path.join(prodigal_path,
                                         '{0}_ORF.txt'.format(identifier))
            # import contigs for current genome/assembly
            contigs = import_sequences(g)
            # extract coding sequences from contigs
            genome_info = extract_coding_sequences(orf_file_path,
                                                       contigs, protid)
            # save coding sequences to file
            # create records and write them to file
            cds_lines = create_fasta_lines(genome_info[0], identifier)
            write_fasta(cds_lines, cds_file)

            write_protein_table(protein_table, identifier, genome_info[1])

            # keep track of CDSs identifiers to assign them sequentially
            protid = 1
        except Exception as e:
            print(e, identifier, orf_file_path, g)

    return [protein_table, cds_file]


def hash_file(file, read_mode):
    """
    """

    with open(file, read_mode) as f:
        hash_obj = hashlib.blake2b()
        file_content = f.read()
        hash_obj.update(file_content)
        hash_str = hash_obj.hexdigest()

    return hash_str


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


def write_gene_list(schema_dir):
    """
    """

    schema_files = [file for file in os.listdir(schema_dir) if '.fasta' in file]
    schema_list_file = os.path.join(schema_dir, '.genes_list')
    pickle_dumper(schema_list_file, schema_files)

    return [os.path.isfile(schema_list_file), schema_list_file]


def write_schema_config(blast_score_ratio, ptf_hash, translation_table,
                        minimum_sequence_length, chewie_version, size_threshold,
                        word_size, clustering_sim, representative_filter,
                        intra_filter, output_directory):
    """
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
    pickle_dumper(config_file, params)

    return [os.path.isfile(config_file), config_file]


def select_name(result):
    """ Extracts the annotation description from the result
        of a query to the UniProt SPARQL endpoint.

        Args:
            result (dict): a dictionary with the results
            from querying the UniProt SPARQL endpoint.
        Returns:
            A list with the following elements:
                - the annotation descrition;
                - the URI to the UniProt page for the protein;
                - a label that has descriptive value.
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
    """ Constructs a SPARQL query to search for exact matches in the
        UniProt endpoint.

        Args:
            sequence (str): the Protein sequence that will be added
            to the query.
        Returns:
            query (str): the SPARQL query that will allow to seaarch for
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
    """
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
    """
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
    """
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
    """ Verify the cpu usage for chewBBACA.

        Args:
            cpu_to_use (int): the number of cpu provided to chewBBACA

        Returns:
            cpu_to_use (int): the number of cpu to use after verification

        Example:

            >>> verify_cpu_usage(6)
            6
    """
    total_cpu = multiprocessing.cpu_count()

    # do not allow a value of cpuToUse greater than the number of cores/threads
    if cpu_to_use >= total_cpu:
        print('Warning! You have provided a CPU core count value '
              'that is equal to or exceeds the number of CPU '
              'cores/threads in your machine!')
        print('Setting a different value...')
        # define a value that is safe according to the number of
        # available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1
        print('CPU cores/threads value set to: {0}'.format(cpu_to_use))

    elif cpu_to_use == (total_cpu - 1):
        print('Warning! You have provided a CPU core count value '
              'that is close to the maximum core count of your '
              'machine ({0}/{1}). This may affect your system '
              'responsiveness.'.format(cpu_to_use, total_cpu))

    return cpu_to_use


def check_ptf(ptf_path):
    """
    """

    if os.path.isfile(ptf_path) is False:
        message = ('Cannot find specified Prodigal training file.\nPlease provide a '
                   'valid training file.\n\nYou can create a training '
                   'file for a species of interest with the following command:\n  '
                   'prodigal -i <reference_genome> -t <training_file.trn> -p single\n\n'
                   'It is strongly advised to provide a high-quality and closed genome '
                   'for the training process.')

        return [False, message]
    else:
        return [True, ptf_path]


def check_prodigal_output_files(path_to_temp, list_of_genomes):
    """ Checks if Prodigal created ORF files
        equal to the number of genome files provided.

        Args:
            path_to_temp (str): the full path to the 'temp'
            directory created by chewBBACA.
            list_of_genomes (list): list containing the full
            path to the input genomes.

        Returns:
            prints a message if Prodigal created all the
            necessary files, otherwise raises a ValueError
    """

    # list ORF files created in the 'temp' directory
    listOfORFCreated = []
    for orffile in os.listdir(path_to_temp):
        if orffile.endswith('_ORF.txt'):
            listOfORFCreated.append(orffile)

    # raise exception if the number of ORF files is not equal
    # to the number of genome files provided
    if len(list_of_genomes) > len(listOfORFCreated):
        message = ('Missing some ORF files from Prodigal run. '
                   'Missing {0} ORF files out of {1} genome files '
                   'provided.'.format(len(list_of_genomes) - len(listOfORFCreated),
                                      len(list_of_genomes)))
        # remove 'temp' directory
        shutil.rmtree(path_to_temp)
        raise ValueError(message)
    else:
        print('Prodigal created all the necessary files.')


def is_fasta(filename):
    """ Checks if a file is a FASTA file.

        Args:
            filename (str): the full path to the FASTA file

        Returns:
            True if FASTA file,
            False otherwise
    """

    with open(filename, 'r') as handle:
        try:
            fasta = SeqIO.parse(handle, 'fasta')
        except:
            fasta = [False]

        # returns True if FASTA file, False otherwise
        return any(fasta)



def filter_files(files_list, suffixes):
    """ Checks if files names contain any suffix from a list of suffixes.

        Args:
            files_list (list): a list with all files names.
        Returns:
            suffixes (list): a list with all suffixes to search for in
            the files names.
    """

    accepted = [file for file in files_list
                if any([True for suffix in suffixes if suffix in file])]

    return accepted


def filter_non_fasta(files_list):
    """ Creates a new list of files names/paths that only contains FASTA files.

        Args:
            files_list (list): a list with files names/paths.
        Returns:
            real_fasta (list): a list with files names/paths that correspond
            to FASTA files.
    """

    real_fasta = [file for file in files_list if is_fasta(file) is True]

    return real_fasta


def gene_seqs_info(gene):
    """ Determines the total number of alleles and the mean length
        of allele sequences per gene.

        Args:
            genes_list (list): a list with names/paths for FASTA
            files.
        Returns:
            genes_info (list): a list with a sublist for each input
            gene file. Each sublist contains a gene identifier, the
            total number of alleles for that gene and the mean length
            of allele sequences for that gene.
    """

    seq_generator = SeqIO.parse(gene, 'fasta', generic_dna)
    alleles_lengths = [len(allele) for allele in seq_generator]
    mean_length = sum(alleles_lengths)/len(alleles_lengths)
    total_seqs = len(alleles_lengths)
    genes_info = [gene, total_seqs, mean_length]

    return genes_info


def make_blast_db(input_fasta, output_path, db_type):
    """ Creates a BLAST database.

        Args:
            input_fasta (str): path to the input file with sequences.
            output_path (str): path to the output database.
            db_type (str): type of the database, nucleotide (nuc) or
            protein (prot).

        Returns:
            Creates a BLAST database with the input sequences.
    """

    makedb_cmd = ('makeblastdb -in {0} -out {1} -parse_seqids '
                  '-dbtype {2} > /dev/null'.format(input_fasta,
                                                   output_path,
                                                   db_type))
    os.system(makedb_cmd)


def determine_blast_task(proteins):
    """ Determine the type of task that should be used to run BLAST.

        Args:
            proteins (str): path to a file with sequences.

        Returns:
            blast_task (str): a string that indicates the type of BLAST
            task to run.
    """

    blast_task = 'blastp'
    proteins_lengths = [len(p) for p in proteins]
    minimum_length = min(proteins_lengths)
    if minimum_length < 30:
        blast_task = 'blastp-short'

    return blast_task


def create_directory(directory_path):
    """ Creates a diretory if it does not exist."""

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


def join_paths(parent_path, child_path):
    """ Joins a parent directory and a subdirectory."""

    joined_paths = os.path.join(parent_path, child_path)

    return joined_paths


def check_input_type(input_path, output_file):
    """ Checks if the input is a file or a directory.

        Args:
            folder_or_list (str): the full path to the file or directory

        Returns:
            list_files (str) if folder_or_list is a path to a file,
            list_files (list) if folder_or_list is a path to a directory,
            Raises Exception otherwise
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


def escape_special_characters(a_string):
    """ Escapes strings to use in regex

        Args:
            a_string (str): string containing characters to escape

        Returns:
            escaped (str): escaped string
    """

    escaped = re.escape(a_string)

    return escaped


def replace_multiple_characters(namefile):
    """ Replaces multiple characters in a string

        Args:
            namefile (str): string containing the name of the contig
            with characters to replace.

        Returns:
            replaced (str): string containing the name of the contig
            without characters to replace.
    """

    replaced = namefile.replace("|", "_")\
                       .replace("_", "-")\
                       .replace("(", "")\
                       .replace(")", "")\
                       .replace("'", "")\
                       .replace("\"", "")\
                       .replace(":", "")

    return replaced


def listdir_fullpath(path):
    """ Gets the full path of the files from a directory

        Args:
            path (str): full path to a directory

        Returns:
            list containing the full path of every file
            contained in the input directory.
    """

    return [os.path.join(path, f) for f in os.listdir(path)]


def flatten_list(list_to_flatten):
    """Flattens one level of a nested list

        Args:
            list_to_flatten (list)

        Returns:
            flattened list

        Example:

            >>> flatten_list([[[1,2],[3,4]]])
            [[1, 2], [3, 4]]

    """

    return list(itertools.chain(*list_to_flatten))


def invert_dictionary(dictionary):
    """ Inverts a dictionary. Keys become values and vice-versa

        Args:
            dictionary (dict)

        Returns:
            inverted (dict): inverted dictionary

        Example:

            >>> inverted_dictionary({key:value})
            {value:key}
    """

    inverted = {value: key for key, value in dictionary.items()}

    return inverted


def threads_for_blast(files_to_blast, cpu_to_apply):
    """ Define the number of threads for BLAST

        Args:
            files_to_blast (list): list containing the full
            path to the files to BLAST
            cpu_to_apply (int): number of cpu to use

        Returns:
            blast_threads (list): list contaning the number
            of threads to use for each file.
            proc (int): Number of processes to use in multiprocessing
    """

    # define number of processes and available cores for each BLAST
    if len(files_to_blast) >= cpu_to_apply:
        blast_threads = [1 for protogenome in files_to_blast]
        proc = cpu_to_apply

    elif cpu_to_apply % len(files_to_blast) == 0:
        blast_threads = [int(cpu_to_apply / len(files_to_blast))
                         for protogenome in files_to_blast]
        proc = len(blast_threads)

    elif cpu_to_apply % len(files_to_blast) == 1:
        blast_threads = [2] + [1 for protogenome in range(0,len(files_to_blast)-1)]
        proc = len(blast_threads)

    elif cpu_to_apply % len(files_to_blast) > 1:
        base_cpu = int(cpu_to_apply / len(files_to_blast))
        blast_threads = [base_cpu
                         for protogenome in range(0, len(files_to_blast))]
        extra_cpu = cpu_to_apply - sum(blast_threads)
        i = 0
        while extra_cpu > 0:
            blast_threads[i] += 1
            extra_cpu -= 1
            i += 1
        proc = len(blast_threads)

    return blast_threads, proc


def isListEmpty(inList):
    """ Checks if a nested list is empty
    """
    if isinstance(inList, list): # Is a list
        return all(map(isListEmpty, inList)) if isinstance(inList, list) else False


def read_blast_tabular(blast_tabular_file):
    """ Read a file with BLAST results in tabular format

        Args:
            blast_tabular_file (str): path to output file of BLAST.

        Returns:
            blasting_results (list): a list with a sublist per line
            in the input file.
    """

    with open(blast_tabular_file, 'r') as blastout:
        reader = csv.reader(blastout, delimiter='\t')
        blasting_results = [row for row in reader]

    return blasting_results


def fasta_lines(identifiers, sequences_dictionary):
    """ Creates list with line elements for a FASTA file based on the sequence
        identifiers passed.

        Args:
            identifiers (list): a list with the identifiers of sequences that
            will be included in the list.
            sequences_dictionary (dict): a dictionary with sequence identifeirs
            as keys and sequences as values.

        Returns:
            seqs_lines (list): a list with strings representing the header of
            the sequence and the sequence for each of the specified sequence
            identifiers.
    """

    seqs_lines = ['>{0}\n{1}\n'.format(seqid, sequences_dictionary[seqid])
                  for seqid in identifiers]

    return seqs_lines


def write_list(lines, output_file):
    """ Writes list elements to file.

        Args:
            lines (list): list with the ordered lines that will be written
            to the output file.
            output_file (str): name/path of the output file.

        Returns:
            Writes contents of 'lines' argument into 'output_file'.
    """

    with open(output_file, 'w') as file:
        file.writelines(lines)


def determine_duplicated_prots(proteins):
    """ Creates a dictionary with protein sequences as keys and all sequence
        identifiers associated with that protein as values.

        Args:
            proteins (dict): dictionary with protein sequence identifiers as
            keys and protein sequences as values.

        Returns:
            equal_prots (dict): dictionary with protein sequence as keys and
            sequence identifiers that are associated with each protein sequence
            as values.
    """
    # use sequences hashes as keys and protids as values
    # read file and process generator to save memory???
    equal_prots = {}
    for protid, protein in proteins.items():
        # if protein sequence was already added as key
        if protein in equal_prots:
            # append new protid
            equal_prots[protein].append(protid)
        # else add new protein sequence as key and protid
        # as value
        else:
            equal_prots[protein] = [protid]

    return equal_prots


def sequences_lengths(fasta_file):
    """ Determines the length of all DNA sequences in a FASTA file.

        Parameters
        ----------
        fasta_file : str
            Path to a FASTA file with DNA sequences.

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


def determine_longest(seqids, proteins):
    """ Determines which sequence is the longest among
        sequences with the specified identifiers.
    """

    seqids_prots = [(seqid, proteins[seqid]) for seqid in seqids]
    sorted_prots = sorted(seqids_prots, key= lambda x: len(x[1]), reverse=True)
    chosen = sorted_prots[0][0]

    return chosen


def locus_mode(alleles):
    """ Determines the mode value from a set of sequence length values.

        Args:
            alleles (dict): dictionary with alleles identifiers as keys
            and the allele length as value.
        Returns:
            modes (list): The most frequent length values. The distribution
            of length values for a locus might have more than one mode.
    """

    # determine frequency of each length value
    counts = Counter(alleles.values())
    # order by most common first
    most_common = counts.most_common()

    # get most common
    modes = [most_common[0][0]]
    # determine if there are more length values that are as common
    modes += [m[0] for m in most_common[1:] if m[1] == most_common[0][1]]

    return modes


def mode_filter(alleles, size_threshold):
    """ Determines the mode from a set of input sequences
        and identifies sequences that have a length value
        smaller or greater than the mode based on a threshold.

        Args:
            alleles (dict):
            size_threshold (float):
        Returns:
            A list with the following variables:
                - modes (list):
                - alm (list):
                - asm (list):
                - alleles_lengths (dict):
    """

    alm = []
    asm = []

    # determine length value of all sequences
    alleles_lengths = {seqid: len(seq) for seqid, seq in alleles.items()}

    # determine mode/s
    modes = locus_mode(alleles_lengths)
    # determine top and bot length value limits
    max_mode = max(modes)
    top_limit = max_mode + (max_mode*size_threshold)
    min_mode = min(modes)
    bot_limit = min_mode - (min_mode*size_threshold)

    # determine sequences that are below or above limits
    alm = [seqid for seqid, length in alleles_lengths.items() if length > top_limit]
    asm = [seqid for seqid, length in alleles_lengths.items() if length < bot_limit]

    return [modes, alm, asm, alleles_lengths]


def get_seqs_dicts(gene_file, gene_id, table_id, min_len, size_threshold, max_proteins=None):
    """ Creates a dictionary mapping seqids to DNA sequences and
        another dictionary mapping protids to protein sequences.

        Args:
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

    seqid = 1
    dna_seqs = {}
    prot_seqs = {}
    seqids_map = {}
    invalid_alleles = []
    seq_generator = SeqIO.parse(gene_file, 'fasta', generic_dna)
    if max_proteins is None:
        translated_seqs = [(rec.id, translate_dna(str(rec.seq), table_id, min_len)) for rec in seq_generator]
    else:
        translated_seqs = []
        exausted = False
        invalid = 0
        seen = []
        while (len(translated_seqs)-invalid) < max_proteins and exausted is False:
            current_rec = next(seq_generator, None)
            if current_rec is not None:
                recid = current_rec.id
                sequence = str(current_rec.seq)
                prot = (recid, translate_dna(sequence, table_id, min_len))
                if isinstance(prot[1], str) is True:
                    invalid += 1
                    translated_seqs.append(prot)
                else:
                    if prot[1] not in seen:
                        translated_seqs.append(prot)
                        seen.append(prot[1])
            else:
                exausted = True

    total_seqs = len(translated_seqs)
    for rec in translated_seqs:
        # if the allele identifier is just an integer
        # add gene identifier as prefix
        try:
            int_seqid = int(rec[0])
            # Python converts '2_1' to 21
            if '_' in rec[0]:
                int_seqid = int(rec[0].split('_')[-1])
            new_seqid = '{0}_{1}'.format(gene_id, int_seqid)
        except Exception:
            new_seqid = rec[0]

        # if returned value is a list, translation was successful
        if isinstance(rec[1], list):
            # we need to assign simple integers as sequence identifiers
            # because BLAST will not work if sequence identifiers are
            # too long
            seqids_map[str(seqid)] = new_seqid
            dna_seqs[new_seqid] = rec[1][0][1]
            prot_seqs[str(seqid)] = str(rec[1][0][0])
            seqid += 1
        # if returned value is a string, translation failed and
        # string contains exceptions
        elif isinstance(rec[1], str):
            invalid_alleles.append([new_seqid, rec[1]])

    if size_threshold is not None and len(prot_seqs) > 0:
        # remove alleles based on length mode and size threshold
        modes, alm, asm, alleles_lengths = mode_filter(dna_seqs, size_threshold)
        excluded = set(asm + alm)

        dna_seqs = {seqid: seq for seqid, seq in dna_seqs.items() if seqid not in excluded}
        prot_seqs = {seqid: seq for seqid, seq in prot_seqs.items() if seqids_map[seqid] not in excluded}

        modes_concat = ':'.join(map(str, modes))
        st_percentage = int(size_threshold*100)
        invalid_alleles += [[s, 'allele greater than {0}% locus length mode '
                                '({1}>{2})'.format(st_percentage, alleles_lengths[s], modes_concat)] for s in alm]
        invalid_alleles += [[s, 'allele smaller than {0}% locus length mode '
                                '({1}<{2})'.format(st_percentage, alleles_lengths[s], modes_concat)] for s in asm]

    return [dna_seqs, prot_seqs,
            invalid_alleles, seqids_map, total_seqs]


def split_genes_by_core(inputs, threads, method):
    """ Splits list with loci inputs into several sublists based
        on the number of sequence per locus (seqcount), the mean
        length of the sequences in each locus or the product of
        both variables.

        Args:
            inputs (list): list with information about the data of
            each locus that needs to be processed.
            threads (int): the number of sublists with inputs that
            should be created, based on the number of CPU threads
            that will be used to process the inputs.
            method (str): "seqcount" to split inputs into sublists
            with even number of sequences, "length" to split based
            on mean length of sequences and "seqcount+length" to
            split based on both criteria.

        Returns:
            splitted_ids (list): subslists with paths to loci, each
            sublist containing paths for a set of loci that should
            not differ much from other sublists based on the criterion
            used to separate the inputs.
    """

    # initialize list with sublists to store inputs
    splitted_ids = [[] for cpu in range(threads)]
    # initialize list with chosen criterion values
    # for each sublist of inputs
    splitted_values = [0 for cpu in range(threads)]
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


def concatenate_list(str_list, join_char):
    """ Concatenates list elements with specified
        character between each original list element.

        Args:
            sequence_ids (list): list with strings that will be concatenated.
            join_char (str): character that will be used to join all list
            elements.

        Returns:
            ids_str (str): concatenation of all strings in the input list.
    """

    concat = join_char.join(str_list)

    return concat


def write_text_chunk(output_file, text):
    """ Write single string to file.

        Args:
            output_file (str): path/name of the file that will store
            the input text.
            text (str): single string to write to file.

        Returns:
            Writes input text to output file.
    """

    with open(output_file, 'w') as out:
        out.write(text)


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Args:
            strDNA (str): string representing a DNA sequence.

        Returns:
            revC_dna (str): the reverse complement of the DNA sequence, without
            lowercase letters.

        Example:
            >>> reverse_complement('ATCGgcaNn')
            'NNTGCCGAT'
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

        Args:
            string (str): string to be reversed.

        Returns:
            revstr (str): reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.

        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Tries to translate an input DNA sequence in specified orientation
        and stores exceptions when the input sequence cannot be translated.

        Args:
            dna_sequence (str): string representing a DNA sequence.
            method (str): a string specifying the way the sequence will
            be oriented to attempt translation.
            table_id (int): translation table identifier.

        Returns:
            List with following elements if translation is successful:
                protseq (str): string representing the translated DNA sequence.
                myseq (str): string representing the DNA sequence in the
                orientation used to translate it.
            Otherwise, returns string derived from captured exception.
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
    """ Determine if a string only has characters from specified
        alphabet.

        Args:
            string (str): input string.
            alphabet (str): string that has all characters from desired
            alphabet.

        Returns:
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

        Args:
            string (str): input string.
            number (int): integer that will be used to check if sequence
            length is multiple of.

        Returns:
            "True" if the length of the sequence is a multiple of the
            specified number and "sequence length is not a multiple of number"
            if condition is not satisfied.
    """

    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def translate_dna(dna_sequence, table_id, min_len):
    """ Checks if sequence is valid and attempts to translate it,
        calling several functions to ensure that the sequence only has
        'ACTG', is multiple of 3 and that it can be translated in any of 4
        different orientations. Stores exceptions so that it is possible to
        understand the sequence could not be translated.

        Args:
            dna_sequence (str):
            table_id (int):

        Returns:
            If the sequence can be translated,
            a list with following elements:
                sequence (list): a list with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
                coding_strand (str): the strand orientation that had could be
                translated.
            Otherwise:
                exception_str (str): a string containing the exceptions that
                determined that the sequence could not be translated.
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

        Args:
            sequence (str): string representing a DNA sequence.
            method (str): a string specifying the sequence orientation
            that should be used to attempt translation.
            table_id (int): translation table identifier.
            strands (list): list with 4 different orientations that can
            be checked.
            exception_collector (list): list used to store all exceptions
            arising from translation attempts.

        Returns:
            A list with following elements, if the sequence can be translated:
                translated_seq (list): a list with the protein sequence and
                with the DNA sequence in the orientation used for translation.
                exception_collector (list): a list with the exceptions that are
                captured when the sequence could not be translated.
            Otherwise:
                translated_seq (str): a string with the exception/reason why
                the sequence could not be translated.
                exception_collector (list): list with all exception that have
                been captured during translation attempts of the current
                sequence.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def is_url(url):
    """ Checks if a url is valid
    
        Args: 
        url (str): the url to be checked

        Returns:
        True if url is valid, False otherwise.
    
    """
    
    try:
        
        result = urlparse(url)
        return all([result.scheme, result.netloc, result.path])
    
    except:
        return False


def make_url(base_url, *res, **params):
    """ Creates a url. 
    
        Args: 
            base_url (str): the base url
            res (str): endpoint(s) to add to the base url
            params (str): addtional parameters (WIP)

        Returns:
            url (str) with the provided parameters.
            Otherwise, returns base_url.

    """
    
    url = base_url
    
    # Check if the url is valid
    if is_url(url):
        
        if url[-1] == "/":
            url = url[:-1]
    
        # Add the endpoints
        for r in res:
            url = '{0}/{1}'.format(url, r)
        
        # Add params if they are provided
        if params:
            url = '{0}?{1}'.format(url, urlencode(params))
        
        return url
    
    else:
        return "An invalid URL was provided."


def get_sequence_from_url(url):
    """
    """
    
    seq = parse_qs(urlsplit(url).query)["sequence"][0]
    
    return seq


def login_user_to_NS(server_url, email, password):
    """ Logs a user in Nomenclature Server
    
        Args:
            server_url (str): url of Nomeclature Server API
            email (str): email of the user in NS
            password (str): password of the user in NS
            
        Returns:
            token (str): authorization token to perform requests to NS
    """
    
    auth_params = {}
    auth_params["email"] = email 
    auth_params["password"] = password
    
    auth_headers = {}
    auth_headers["Content-Type"] = "application/json"
    auth_headers["accepts"] = "application/json"
    
    auth_url = make_url(server_url, "auth", "login")
    
    auth_r = requests.post(auth_url, data=json.dumps(auth_params), headers=auth_headers, verify=False)
    
    auth_result = auth_r.json()
    if auth_result['status'] == 'success':
        token = auth_result["access_token"]
    else:
        token = False
    
    return token


def capture_login_credentials(base_url):
    """
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
    """
    """

    config_file = os.path.join(schema_path, filename)
    if os.path.isfile(config_file):
        # Load configs dictionary
        configs = pickle_loader(config_file)
    else:
        sys.exit('Could not find a valid config file.')

    return configs


def get_species_schemas(schema_id, species_id, base_url, headers_get):
    """ Determines if a species in the Chewie-NS has a schema
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
        list
            A list with the following elements:

            - The schema id (str).
            - The schema URI (str).
            - The schema name (str).

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
            schema_params = requests.get(schema[1], headers=headers_get, verify=False)
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
            HTTPS requests (False for no verification, True otherwise).

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
    """ Uploads data to the NS.

        Parameters
        ----------
        data
            The data that will be sent to the NS (any data
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
    """ Gets data from Virtuoso """

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


def progress_bar(process, total, tickval, ticknum, completed):
    """
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
    """
    """

    pool = ThreadPool(processes=1)
    answer = pool.apply_async(input, args=[prompt])

    try:
        return answer.get(timeout=timeout)
    except TimeoutError as e:
        sys.exit('Timed out.')


#def translate_coding_sequences(sequences_file, valid_seqs, dna_valid_file,
#                               protein_valid_file, table_id):
def translate_coding_sequences(input_data):
    """ Translates CDSs into protein sequences.

        Args:
            cdss (dict): a dictionary with CDSs ids as keys and CDSs DNA
            sequences as values.

        Returns:
            prots (dict): a dictionary with CDSs/proteins ids as keys and
            protein sequences as values.
            trans_state (dict): a dictionary with the CDSs/proteins ids as
            keys and the DNA strand that coded for those proteins.
            ambiguous (dict): a dictionary with CDSs/proteins ids as keys
            and CDSs DNA sequences that had ambiguous bases as values.
            untranslatable (dict): a dictionary with CDSs/proteins ids as
            keys and CDSs DNA sequences that could not be translated.
    """

    sequences_file = input_data[-4]
    valid_seqs = input_data[0:-4]
    dna_valid_file = input_data[-2]
    protein_valid_file = input_data[-1]
    table_id = input_data[-3]

    # define limit of records to keep in memory
    dna_lines = []
    total_seqs = 0
    prot_lines = []
    line_limit = 5000
    invalid_alleles = []
    dna_seqs_file = dna_valid_file
    protein_seqs_file = protein_valid_file

    cds_index = SeqIO.index(sequences_file, 'fasta')

    for i, seqid in enumerate(valid_seqs):
        sequence = str(cds_index.get(seqid).seq)

        translation = translate_dna(sequence, table_id)
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

        if len(dna_lines)//2 == line_limit or i+1 == len(valid_seqs):

            dna_lines = join_list(dna_lines, '\n')
            write_to_file(dna_lines, dna_seqs_file, 'a', '\n')
            dna_lines = []

            prot_lines = join_list(prot_lines, '\n')
            write_to_file(prot_lines, protein_seqs_file, 'a', '\n')
            prot_lines = []

    return [invalid_alleles, total_seqs]


def determine_repeated(input_data):
    """
    """

    sequences_file, unique_fasta = input_data

    total = 0
    seqs_dict = {}
    out_limit = 5000
    out_seqs = []
    for record in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        sequence = str(record.seq.upper())
        seqid = record.id
        seq_hash = hashlib.sha256(sequence.encode('utf-8')).hexdigest()

        # store each sequence and first seqid found
        # with that sequence
        # write sequence and seqid to file and use hash as key
        if seq_hash not in seqs_dict:
            seqs_dict[seq_hash] = seqid
            recout = '>{0}\n{1}'.format(seqid, sequence)
            out_seqs.append(recout)
            if len(out_seqs) == out_limit:
                out_seqs = join_list(out_seqs, '\n')
                write_to_file(out_seqs, unique_fasta, 'a', '\n')
                out_seqs = []
        elif seq_hash in seqs_dict:
            total += 1

    if len(out_seqs) > 0:
        out_seqs = join_list(out_seqs, '\n')
        write_to_file(out_seqs, unique_fasta, 'a', '\n')

    unique_seqids = list(seqs_dict.values())

    return [total, unique_seqids]


def determine_small(sequences_file, minimum_length):
    """ Find protein sequences that are shorter than desired length.

        Args:
            prots (dict): a dictionary with protein ids as keys and protein
            sequences as values.
            min_len (int): Proteins with a number of amino acids lower than
            this value are considered small.

        Returns:
            small_proteins (dict): a dictionary with the ids of small proteins
            as keys and their amino acid sequence as values.
    """

    small_seqs = []
    valid_seqs = []
    for record in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        #sequence = sequences_index.get(seqid).seq
        sequence = str(record.seq)
        seqid = record.id

        if len(sequence) < minimum_length:
            small_seqs.append(seqid)
        else:
            valid_seqs.append(seqid)

    return [small_seqs, valid_seqs]


def get_sequences_by_id(sequences_index, seqids, out_file):
    """
    """

    selected = []
    seq_limit = 5000
    for i, seqid in enumerate(seqids):
        identifier = sequences_index[seqid].id
        header = '>{0}'.format(identifier)
        sequence = str(sequences_index[seqid].seq)
        selected.append(header)
        selected.append(sequence)

        if len(selected) // 2 == seq_limit or i+1 == len(seqids):

            lines = join_list(selected, '\n')
            write_to_file(lines, out_file, 'a', '\n')

            selected = []


def get_schema_reps(dna_index, schema_seqids, protein_file, output_file):
    """
    """

    schema_lines = []
    for seqid in schema_seqids:
        header = '>{0}'.format(seqid)
        sequence = str(dna_index[seqid].seq)

        schema_lines.append(header)
        schema_lines.append(sequence)

    schema_text = join_list(schema_lines, '\n')
    write_to_file(schema_text, output_file, 'w', '')


def build_schema(schema_file, output_path):
    """
    """

    total_genes = 0
    schema_files = []
    for record in SeqIO.parse(schema_file, 'fasta', IUPAC.unambiguous_dna):
        file_name = record.name
        file_name = new_utils.replace_multiple_characters(file_name)

        new_file = '{0}{1}'.format(file_name, '.fasta')
        new_file_path = os.path.join(output_path, new_file)
        schema_files.append(new_file_path)

        header = '>{0}_1'.format(file_name)
        sequence = str(record.seq).upper()
        file_text = join_list([header, sequence], '\n')
        write_to_file(file_text, new_file_path, 'w', '\n')

        total_genes += 1

    init_schema_4_bbaca.get_Short(schema_files)
    print('\nTotal of {0} loci that constitute the schema.'.format(total_genes))


def cluster_blaster(inputs):
    """
    """

    blast_inputs = inputs[0:-4]
    blastp_path = inputs[-4]
    output_directory = inputs[-2]
    proteins_file = inputs[-1]
    blast_db = inputs[-3]

    indexed_fasta = SeqIO.index(proteins_file, 'fasta')

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
        blast_command = ('{0} -db {1} -query {2} -out {3} '
                         '-outfmt "6 qseqid sseqid score" '
                         '-max_hsps 1 -num_threads {4} -evalue '
                         '0.001 -seqidlist {5}'.format(blastp_path,
                                                       blast_db,
                                                       fasta_file,
                                                       blast_output,
                                                       1, ids_file))

        os.system(blast_command)


def blast_inputs(clusters, output_directory, ids_dict):
    """
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


def split_blast_inputs_by_core(blast_inputs, threads, blast_files_dir):
    """
    """

    splitted_ids = [[] for cpu in range(threads)]
    splitted_values = [[] for cpu in range(threads)]
    cluster_sums = [0] * threads
    i = 0
    for cluster in blast_inputs:
        cluster_file = os.path.join(blast_files_dir,
                                    '{0}_ids.txt'.format(cluster))
        with open(cluster_file, 'r') as infile:
            cluster_seqs = [line.strip() for line in infile.readlines()]
        splitted_values[i].append(len(cluster_seqs))
        splitted_ids[i].append(cluster)
        cluster_sums[i] += len(cluster_seqs)
        i = cluster_sums.index(min(cluster_sums))

    return splitted_ids


def divide_list_into_n_chunks(list_to_divide, n):
    """
    """

    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    return sublists


def apply_bsr(inputs):
    """
    """

    ids_dict = inputs[3]
    blast_file = inputs[0]
    blast_results = read_blast_tabular(blast_file)
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
            #print(res)
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
                        #print(query)

                    elif hit_length <= query_length:
                        excluded_alleles.append(hit)
                        #print(hit)
            # it might not work because there is no self score for
            # some sequences due to low complexity regions...
            except Exception:
                excluded_alleles.append(query)

    return excluded_alleles


def sequence_kmerizer(sequence, k_value, offset=1):
    """
    """

    kmers = [sequence[i:i+k_value] for i in range(0, len(sequence)-k_value+1, offset)]

    return kmers


# for AlleleCall this function need to receive the variable 'reps_groups'
# and to accept a argument that controls if it is allowed to create new
# clusters based on new genes that it finds.
def cluster_sequences(sorted_sequences, word_size, clustering_sim, mode, offset=1):#, interval=0.25):
    """
    """

    clusters = {}
    reps_groups = {}
    for prot in sorted_sequences:
        protid = prot[0]
        protein = prot[1]
        kmers = sequence_kmerizer(protein, word_size, offset)

        #first_value = len(kmers) * interval
        #intervals = [first_value, first_value*2, first_value*3, len(kmers)]
        current_reps = []
        for k in kmers:
        #for i, k in enumerate(kmers):
            if k in reps_groups:
                current_reps.extend(list(reps_groups[k]))
            #if i in intervals:
            #    counts = Counter(current_reps)
            #    current_reps = [(k, v/len(kmers)) for k, v in counts.items() if v/len(kmers) >= clustering_sim]

        # count number of kmer hits per representative
        counts = Counter(current_reps)
        current_reps = [(k, v/len(kmers)) for k, v in counts.items() if v/len(kmers) >= clustering_sim]

        # sort to get most similar at index 0
        sims = sorted(current_reps, key=lambda x: x[1], reverse=True)

        if len(sims) > 0:

            if mode == 'greedy':
                clusters[sims[0][0]].append((protid, sims[0][1]))
            elif mode == 'full':
                for s in sims:
                    clusters[s[0]].append((protid, s[1]))
        else:
            for k in kmers:
                if k in reps_groups:
                    reps_groups[k].add(protid)
                else:
                    reps_groups[k] = set([protid])
            clusters[protid] = [(protid, 1.0)]

    return clusters


def cluster_prunner(clusters, sim_cutoff):
    """
    """

    excluded = []
    prunned_clusters = {}
    for rep, seqids in clusters.items():
        # this removes the representative from the cluster
        prunned_clusters[rep] = [seqid for seqid in seqids if seqid[1] < sim_cutoff]
        excluded.extend([seqid for seqid in seqids if seqid[1] >= sim_cutoff and seqid[0] != rep])

    return [prunned_clusters, excluded]


def determine_singletons(clusters):
    """
    """

    singletons = {k: v for k, v in clusters.items() if len(v) == 0}

    return singletons


def remove_clusters(clusters, cluster_ids):
    """
    """

    new_clusters = {k: v for k, v in clusters.items() if k not in cluster_ids}

    return new_clusters


def intra_cluster_sim(clusters, protein_file, word_size, intra_filter):
    """
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
            prot_kmers = sequence_kmerizer(prot, word_size)
            # dict with sequence indentifiers and kmers
            cluster_kmers[seqid] = prot_kmers
            
            # create dict with kmers as keys and list
            # of sequences with given kmers as values
            for kmer in prot_kmers:
                if kmer in kmers_mapping:
                    kmers_mapping[kmer].add(seqid)
                else:
                    kmers_mapping[kmer] = set([seqid])

        sims_cases = {}
        excluded = []
        # for each sequence in the cluster
        for seqid, kmers in cluster_kmers.items():
            if seqid not in excluded:
                query_kmers = kmers
                current_reps = []
                # determine sequences that also have the same kmers
                for kmer in query_kmers:
                    if kmer in kmers_mapping:
                        current_reps.extend(list(kmers_mapping[kmer]))
                
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
