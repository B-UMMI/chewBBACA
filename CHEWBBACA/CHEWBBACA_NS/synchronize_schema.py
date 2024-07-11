#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module allows users to synchronize a local schema, previously
downloaded from the Chewie-NS, with the remote schema in the Chewie-NS.

The synchronization process will retrieve new alleles added to the remote
schema after the last synchronization date and update the local schema,
ensuring that local and remote alleles have the same integer identifiers.
The process also allows users to submit novel alleles and update the remote
schema. Novel alleles that are not submitted keep a '*' in the identifier.

Code documentation
------------------
"""


import os
import sys
import json
import time
import shutil
import requests
from collections import defaultdict
from urllib3.exceptions import InsecureRequestWarning

try:
    from PrepExternalSchema import adapt_schema
    from utils import (constants as ct,
                       file_operations as fo,
                       fasta_operations as fao,
                       chewiens_requests as cr,
                       # profiles_sqlitedb as ps,
                       parameters_validation as pv,
                       iterables_manipulation as im)
except ModuleNotFoundError:
    from CHEWBBACA.PrepExternalSchema import adapt_schema
    from CHEWBBACA.utils import (constants as ct,
                                 file_operations as fo,
                                 fasta_operations as fao,
                                 chewiens_requests as cr,
                                 # profiles_sqlitedb as ps,
                                 parameters_validation as pv,
                                 iterables_manipulation as im)


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


def create_lengths_files(upload, temp_dir):
    """Create files with sequence length values to send to Cheiwe-NS.

    Parameters
    ----------
    upload : dict
        Dictionary with loci identifiers as keys and lists
        as values. Each list contains the following elements:

        - A dictionary with sequences hashes as keys and a list
          with the sequence identifier, DNA sequence and sequence
          length value.
        - The locus URI.
    temp_dir : str
        Path to the directory where the output files with
        the length values of the sequences will be created.

    Returns
    -------
    length_files : list of str
        List with paths to the pickled files that contain
        the sequence length values. Each pickled file
        contains a dictionary with hashes (sequence hashes
        computed with sha256) as keys and sequence lengths
        as values.
    """
    length_files = []
    for locus, v in upload.items():
        records = v[0]
        locus_id = locus.rstrip('.fasta')
        lengths = {locus: {h: info[2] for h, info in records.items()}}

        lengths_file = fo.join_paths(temp_dir,
                                     ['{0}_lengths'.format(locus_id)])
        fo.pickle_dumper(lengths, lengths_file)
        length_files.append(lengths_file)

    return length_files


def create_alleles_files(upload, base_url, user_id, species_name,
                         species_id, schema_id, temp_dir):
    """Create files with the essential data to insert alleles in Chewie-NS.

    Parameters
    ----------
    upload : dict
        Dictionary with loci identifiers as keys and lists
        as values. Each list contains the following elements:

        - A dictionary with sequences hashes as keys and a list
          with the sequence identifier, DNA sequence and sequence
          length value.
        - The locus URI.
    base_url : str
        Base URL of the Chewie-NS.
    user_id : int
        Current user identifier.
    species_name : str
        Name of the species.
    species_id : int
        The identifier of the schema's species in the NS.
    schema_id : int
        The identifier of the schema in the NS.
    temp_dir : str
        Path to the directory where the output files  will
        be created.

    Returns
    -------
    list of list
        List with the following elements:

        - List with paths to files that contain the data that
          is necessary to upload and insert new alleles.
        - List with the loci identifiers in the NS.
        - List with loci basenames.
    """
    loci_ids = []
    loci_names = []
    alleles_files = []
    user_uri = '{0}users/{1}'.format(base_url, user_id)
    for locus, recs in upload.items():
        locus_uri = recs[1]
        alleles_sequences = [v[1] for k, v in recs[0].items()]

        post_inputs = [locus_uri, species_name,
                       user_uri, tuple(alleles_sequences)]

        locus_id = locus_uri.split('/')[-1]
        loci_ids.append(locus_id)

        locus_name = locus.rstrip('.fasta')
        loci_names.append(locus_name)

        alleles_file = fo.join_paths(temp_dir,
                                     ['{0}_{1}_{2}'.format(species_id,
                                                           schema_id,
                                                           locus_id)])
        alleles_files.append(alleles_file)

        fo.pickle_dumper(post_inputs, alleles_file)

    return [alleles_files, loci_ids, loci_names]


def upload_alleles_data(alleles_data, length_files, base_url, headers_post,
                        headers_post_bytes, species_id, schema_id):
    """Upload files with the data to insert alleles in Chewie-NS.

    Parameters
    ----------
    alleles_data : list
        List with tuples, one per locus, that contain the path
        to the ZIP archive with the data to insert alleles,
        the identifier of the locus, the locus file hash and
        the basename of the locus file.
    length_files : list
        List with paths to the pickled files that contain a
        dictionary with sequences hashes as keys and sequence
        length as values.
    base_url : str
        Base URL of the Nomenclature server.
    headers_post : dict
        HTTP headers for POST requests that accept JSON
        formatted data.
    headers_post_bytes : dict
        HTTP headers for POST requests that support file
        upload.
    species_id : int
        The identifier of the schema's species in the NS.
    schema_id : int
        The identifier of the schema in the NS.

    Returns
    -------
    failed : list of str
        List with the identifiers of the loci whose alleles
        data could not be fully uploaded.
    zip_res : dict
        A dictionary with the response returned by the last
        POST method. It has loci identifiers as keys and
        lists with two dictionaries as values (the dictionaries
        have sequences hashes as keys and sequence identifiers in
        the Chewie-NS as values. The first dictionary has the hashes
        of the sequences that were sent to the Chewie-NS but that were
        already present in the loci and the identifiers of those repeated
        alleles that were sent to the Chewie-NS. The second dictionary
        has the same structure but for the sequences that were accepted and
        inserted into each locus).
    """
    uploaded = 0
    failed = []
    for i, a in enumerate(alleles_data):
        locus_id = a[1]
        # get length of alleles from current locus
        current_len = length_files[i]
        data = fo.pickle_loader(current_len)
        data = {locus_id: data[next(iter(data))]}
        data = {'content': data}

        # send data to the NS
        send_url = cr.make_url(base_url, 'species', species_id,
                               'schemas', schema_id, 'loci',
                               locus_id, 'lengths')

        lengths_res = cr.simple_post_request(send_url, headers_post,
                                             data=json.dumps(data))[1]
        length_status = lengths_res.status_code

        # get path to ZIP archive with data to insert alleles
        current_zip = a[0]

        # send data to insert alleles in the NS
        zip_url = cr.make_url(base_url, 'species', species_id,
                              'schemas', schema_id, 'loci',
                              locus_id, 'update')

        if alleles_data[i] == alleles_data[-1]:
            headers_post_bytes['complete'] = 'True'

        zip_res = cr.upload_file(current_zip, os.path.basename(current_zip),
                                 zip_url, headers_post_bytes,
                                 False)

        # determine if upload was successful
        zip_status = zip_res.status_code

        # determine if upload was successful
        if length_status not in [200, 201] or zip_status not in [200, 201]:
            failed.append(locus_id)
        elif length_status in [200, 201] and zip_status in [200, 201]:
            uploaded += 1
            print('\r', '    Sent data for alleles of '
                  '{0}/{1} loci.'.format(uploaded, len(alleles_data)), end='')

    return [failed, zip_res]


def pickle_to_fasta(locus, pickled_file, temp_dir, identifiers, reassigned):
    """Create FASTA files from pickled files.

    Parameters
    ----------
    locus : str
        The identifier of the locus with '.fasta' suffix.
    pickled_file : str
        Path to the pickled file with a dictionary that
        has integer identifiers as keys and a tuple with
        two elements: the identifier that should be assigned
        to the allele (might differ from the key if the allele
        is new, in which case it starts with '*') and the DNA
        sequence of the allele.
    temp_dir : str
        Path to the directory where the output FASTA file will
        be created.
    identifiers : dict
        The `zip_res` variable returned by the
        :py:func:`upload_alleles_data` function. It will be used
        to change allele identifiers that were successfully
        inserted into the Chewie-NS.

    Returns
    -------
    fasta_path : str
        Path to the FASTA file created by this function
    """
    locus_id = locus.rstrip('.fasta')
    locus_int = locus_id.split('-')[-1].lstrip('0')
    if locus_int in identifiers:
        repeated = identifiers[locus_int][0]
        attributed = identifiers[locus_int][1]
    else:
        repeated = {}
        attributed = {}

    inv_reassigned = {}
    if locus in reassigned:
        inv_reassigned = {v: k for k, v in reassigned[locus].items()}

    locus_sequences = fo.pickle_loader(pickled_file)

    natsorted_locus = sorted(locus_sequences)

    fasta_path = fo.join_paths(temp_dir, [locus])
    records = []
    for seqid in natsorted_locus:
        recid = locus_sequences[seqid][0]
        seq = locus_sequences[seqid][1]
        seq_hash = im.hash_sequence(seq)
        # switch by the identifier attributed by the Chewie-NS
        if seq_hash in attributed:
            new_recid = attributed[seq_hash]
            if recid in inv_reassigned:
                old_id = inv_reassigned[recid]
                reassigned[locus][old_id] = new_recid
            else:
                if locus not in reassigned:
                    reassigned[locus] = {recid: new_recid}
                else:
                    reassigned[locus][recid] = new_recid
            recid = new_recid
        elif seq_hash in repeated:
            new_recid = repeated[seq_hash]
            if recid in inv_reassigned:
                old_id = inv_reassigned[recid]
                reassigned[locus][old_id] = new_recid
            else:
                if locus not in reassigned:
                    reassigned[locus] = {recid: new_recid}
                else:
                    reassigned[locus][recid] = new_recid
            recid = new_recid

        record = '>{0}_{1}\n{2}'.format(locus_id, recid, seq)
        records.append(record)

    fasta_text = '\n'.join(records)

    with open(fasta_path, 'w') as fp:
        fp.write(fasta_text)

    os.remove(pickled_file)

    return reassigned


def retrieve_alleles(loci_new_alleles, server_time, schema_uri,
                     count, headers_get, ns_date):
    """Retrieve alleles added to Chewie-NS since last sync date.

    Retrieves alleles added to a schema in Chewie-NS
    during a time interval, up to the maximum number of
    alleles that the server returns at a time (50000).

    Parameters
    ----------
    loci_new_alleles : dict
        A dictionary with loci identifiers as keys and
        dictionaries with alleles identifiers and DNA
        sequences as values.
    server_time : str
        The function will return alleles added to the
        schema after this date (format Y-%m-%dT%H:%M:%S).
    schema_uri : str
        The URI of the schema in the Chewie-NS.
    count : int
        The cumulative number of sequences that has been
        returned.
    headers_get : dict
        HTTP headers for GET requests.
    ns_date : str
        The function will return alleles added to the
        schema up to this date (format Y-%m-%dT%H:%M:%S).

    Returns
    -------
    A list with the following variables:

    loci_new_alleles :  dict
        Input `loci_new_alleles` dictionary with alleles
        returned in the current and previous iterations.
    server_time : str
        The date of insertion of the last allele returned by
        the Chewie-NS.
    count : int
        The cumulative number of sequences that has been
        returned.
    """
    # request the new alleles starting on the date given
    url = cr.make_url(schema_uri, 'loci')
    payload = {'local_date': server_time, 'ns_date': ns_date}
    # get the new alleles
    response = cr.simple_get_request(url, headers_get, parameters=payload)[1]
    response_content = response.json()

    # get info about sequences that were added since last date
    new_alleles = response_content['newAlleles']

    # get headers info
    response_headers = response.headers
    if len(new_alleles) > 0:
        # get date of last added allele
        server_time = response_headers['Last-Allele']
        # group retrieved alleles by locus
        for allele in new_alleles:
            locus = '{0}{1}'.format(allele['locus_name']['value'], '.fasta')
            allele_id = allele['allele_id']['value']
            sequence = allele['nucSeq']['value']

            loci_new_alleles[locus][allele_id] = sequence

        # keep count of the number of retrieved alleles
        count += len(new_alleles)
    else:
        # get current server date if no alleles were added
        # since last server time
        server_time = response_headers['Server-Date']

    return (loci_new_alleles, server_time, count)


def altered_loci(loci, schema_dir, pickled_loci, not_in_ns,
                 temp_dir, rearranged):
    """Reassign local allele identifiers based on Chewie-NS nomenclature.

    Reassigns allele identifiers of loci that were altered
    in the Chewie-NS since last sync process. Identifier
    reassignment will ensure that local and remote alleles
    have the same identifiers or that alleles that are only
    local are assigned an identifier that does not match any
    identifier in the Chewie-NS.

    Parameters
    ----------
    loci : dict
        A dictionary with loci identifiers as keys and nested
        dictionaries as values (one nested dictionary per entry,
        with alleles identifiers as keys and DNA sequences as values).
    schema_dir : str
        Path to the schema's diretory.
    pickled_loci : dict
        A dictionary that will be used to store paths to
        pickled files with the records for each locus that is
        processed.
    not_in_ns : dict
        A dictionary that will be used to store information
        about local alleles that are not in the Chewie-NS.
    temp_dir : str
        Path to the directory where the output FASTA file will
        be created.

    Returns
    -------
    A list with three elements:

    - pickled_loci : dict
        Input `pickled_loci` variable updated with the paths
        to the pickled files with the records of the loci that
        were altered in the Chewie-NS since last sync process.
    - not_in_ns : dict
        Input `not_in_ns` variable that was updated with the
        information about all alleles that are not in the
        Chewie-NS.
    - rearranged : dict
        A dictionary with loci identifiers as keys and lists
        as values. Each list contains one or more two-element
        tuples that represent identifiers reassignments
        (e.g.: (*2, 3), *2 --> 3 ).
    """
    for locus, alleles in loci.items():
        locus_id = locus.rstrip('.fasta')
        # paths for current and temp locus file
        locus_file = fo.join_paths(schema_dir, [locus])

        # get the latest alleles retrieved from Chewie-NS
        ns_seqs = {seq: seqid for seqid, seq in alleles.items()}

        local_seqs = {str(rec.seq): [rec.id, (rec.id).split('_')[-1]]
                      for rec in fao.sequence_generator(locus_file)}

        # create dictionary with updated records
        # add local records without '*'
        updated_records = {v[0]: [v[1], k] for k, v in local_seqs.items() if '*' not in v[1]}

        # add records that were added to the remote schema since last sync
        switched = {}
        for seq, seqid in ns_seqs.items():
            updated_records[seqid] = [seqid.split('_')[-1], seq]
            if seq in local_seqs:
                switched[local_seqs[seq][1]] = seqid
                # remove sequence from local sequences to avoid sending
                # duplicated sequences to Chewie-NS
                del local_seqs[seq]

        # identify largest allele identifier in Chewie-NS
        max_id = max([int(v[0]) for k, v in updated_records.items()])

        local_ids = [(v[1], k) for k, v in local_seqs.items() if '*' in v[0]]
        sorted_novel = sorted(local_ids, key=lambda x: int(x[0][1:]))

        for si in sorted_novel:
            max_id += 1
            new_id = '*{0}'.format(max_id)
            new_seqid = '{0}_{1}'.format(locus_id, new_id)
            updated_records[new_seqid] = [new_id, si[1]]
            switched[si[0]] = new_id

        if len(updated_records) > 0:
            novel_alleles = {im.hash_sequence(v[1]):
                             [v[0], v[1], len(v[1])]
                             for k, v in updated_records.items() if '*' in v[0]}
            not_in_ns[locus] = [novel_alleles]

        if len(switched) > 0:
            rearranged[locus] = switched

        updated_records2 = {}
        for seqid, v in updated_records.items():
            int_seqid = int(v[0]) if '*' not in v[0] else int(v[0][1:])
            updated_records2[int_seqid] = (seqid, v[1])

        temp_file = fo.join_paths(temp_dir, ['{0}_pickled'.format(locus_id)])
        fo.pickle_dumper(updated_records2, temp_file)

        pickled_loci[locus] = temp_file

    return [pickled_loci, not_in_ns, rearranged]


def unaltered_loci(loci, schema_dir, pickled_loci, not_in_ns, temp_dir):
    """Determine set of loci that were altered locally only.

    Determines if local schema has new alleles for loci that
    were not updated in the Chewie-NS. Creates a pickled
    file for each locus that has new alleles.

    Parameters
    ----------
    loci : dict
        A dictionary with loci identifiers as keys and nested
        dictionaries as values (one nested dictionary per entry,
        with alleles identifiers as keys and DNA sequences as values).
    schema_dir : str
        Path to the schema's diretory.
    pickled_loci : dict
        A dictionary that will be used to store paths to
        pickled files with the records for each locus that is
        processed.
    not_in_ns : dict
        A dictionary that will be used to store information
        about local alleles that are not in the Chewie-NS.
    temp_dir : str
        Path to the directory where the output FASTA file will
        be created.

    Returns
    -------
    A list with two elements:

    - pickled_loci : dict
        Input `pickled_loci` variable updated with the paths
        to the pickled files with the records of the loci that
        were not altered in the Chewie-NS since last sync process
        but were altered locally.
    - not_in_ns : dict
        Input `not_in_ns` variable that was updated with the
        information about all alleles that are not in the
        Chewie-NS.
    """
    for gene in loci:
        locus_id = gene.rstrip('.fasta')
        locus_file = fo.join_paths(schema_dir, [gene])
        # get local locus sequences
        records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                   for rec in fao.sequence_generator(locus_file)}

        # determine if there are sequences with '*'
        novel_local = {im.hash_sequence(v[1]):
                       [v[0], v[1], len(v[1])]
                       for k, v in records.items() if '*' in v[0]}
        if len(novel_local) > 0:
            not_in_ns[gene] = [novel_local]
            updated_records = {}
            for seqid, seq in records.items():
                int_seqid = int(seq[0].replace('*', ''))
                updated_records[int_seqid] = (seq[0], seq[1])

            temp_file = fo.join_paths(temp_dir, ['{0}_pickled'.format(locus_id)])
            fo.pickle_dumper(updated_records, temp_file)

            pickled_loci[gene] = temp_file

    return [pickled_loci, not_in_ns]


def update_loci_files(loci, local_loci, schema_dir, temp_dir):
    """Determine set loci/alleles to update in the local and remote schemas.

    Determines which loci were or were not updated in the
    remote schema, alleles that are not in the Chewie-NS,
    reassigns identifiers to ensure that alleles that are
    common to local and remote schema have the same
    identifiers and saves updated records to pickled files
    to allow the update of the local schema and of the
    remote schema by sending novel alleles to the Chewie-NS.

    Parameters
    ----------
    loci : dict
        Dictionary with the alleles added to the remote
        schema since last sync date.
    local_loci : list
        List with the names of the schema's loci files.
    schema_dir : str
        Path to the schema's diretory.
    temp_dir : str
        Path to the directory where the output pickled files
        will be created.

    Returns
    -------
    A list with five elements:

    - not_in_ns : dict
        Dictionary with information about local alleles that are not
        in the Chewie-NS. It has loci identifiers as keys and lists
        as values. Each list contains the following elements:

        - A dictionary with sequences hashes as keys and a list
          with the sequence identifier, DNA sequence and sequence
          length value.
        - The locus URI.
    - pickled_loci : dict
        A dictionary with paths to pickled files that contain the
        records for each locus that has been altered in the remote
        schema and/or locally with new alleles added during the
        AlleleCall process.
    - updated : dict
        Dictionary with the alleles added to the remote
        schema since last sync date.
    - not_updated : list
        List of loci that were not altered in the remote schema
        since last sync process.
    - rearranged : dict
        A dictionary with loci identifiers as keys and lists
        as values. Each list contains one or more two-element
        tuples that represent identifiers reassignments
        (e.g.: (*2, 3), *2 --> 3 ).
    """
    rearranged = {}
    not_in_ns = {}
    pickled_loci = {}
    # loci in the remote schema that were altered
    updated = {gene: loci[gene] for gene in local_loci if gene in loci}
    # loci that were not altered in the remote schema
    not_updated = [gene for gene in local_loci if gene not in loci]

    # update local identifiers based on identifiers in the Chewie-NS
    pickled_loci, not_in_ns, \
        rearranged = altered_loci(updated, schema_dir, pickled_loci,
                                  not_in_ns, temp_dir, rearranged)

    # determine loci that were not altered in the Chewie-NS but
    # have been altered locally
    pickled_loci, not_in_ns = unaltered_loci(not_updated, schema_dir,
                                             pickled_loci, not_in_ns,
                                             temp_dir)

    return [not_in_ns, pickled_loci, rearranged]


def retrieve_latest(local_date, schema_uri, headers_get, ns_date):
    """Retrieve alleles added to a schema in Chewie-NS after a specified date.

    Parameters
    ----------
    local_date : str
        Last sync date. The function will return alleles added
        to the remote schema after this date (format Y-%m-%dT%H:%M:%S).
    schema_uri : str
        The URI of the schema in the Chewie-NS.
    headers_get : dict
        HTTP headers for GET requests.
    ns_date : str
        Last modification date of the remote schema in the Chewie-NS.
        The function will return alleles added to the schema up to
        this date (format Y-%m-%dT%H:%M:%S).

    Returns
    -------
    A list with the following variables:

    new_alleles : dict
        A dictionary with loci identifiers as keys and
        dictionaries with alleles identifiers and DNA
        sequences as values.
    server_time : str
        The last modification date of the remote schema.
    count : int
        The total number of sequences that were retrieved from
        the Chewie-NS.
    """
    # get list of sequences that are new considering the last date
    count = 0
    server_time = local_date
    new_alleles = defaultdict(dict)

    # get new alleles from Chewie-NS
    while server_time != ns_date:
        new_alleles, server_time, count = retrieve_alleles(new_alleles,
                                                           server_time,
                                                           schema_uri,
                                                           count,
                                                           headers_get,
                                                           ns_date)

    return [new_alleles, server_time, count]


def main(schema_directory, cpu_cores, nomenclature_server,
         submit, blast_path):  # update_profiles):

    # get ns configs
    local_date, schema_uri = pv.read_configs(schema_directory, '.ns_config')
    # get schema and species identifiers
    schema_id = schema_uri.split('/')[-1]
    species_id = schema_uri.split('/')[-3]
    if nomenclature_server is None:
        nomenclature_server = schema_uri.split('species/')[0]

    if submit is True and 'tutorial' not in nomenclature_server:
        print('\nOnly authorized registered users may submit new alleles.')
        token = cr.capture_login_credentials(nomenclature_server)
    else:
        token = ''

    # GET request headers
    headers_get = ct.HEADERS_GET_JSON
    headers_get['Authorization'] = token

    # determine current user ID and Role
    if submit is True and 'tutorial' not in nomenclature_server:
        user_id, user_role, user_auth = cr.user_info(nomenclature_server, headers_get)
        # verify if user has authorization to submit
        url = cr.make_url(nomenclature_server, 'auth', 'check')
        response = cr.simple_get_request(url, headers_get)[1]
        if response.status_code == 200:
            user_auth = True
        else:
            sys.exit('Current user has no authorization to submit novel '
                     'alleles.\nYou can request authorization to submit '
                     'novel alleles by sending an e-mail to: '
                     'imm-bioinfo@medicina.ulisboa.pt')
        print('User id: {0}'.format(user_id))
        print('User role: {0}\n'.format(user_role))
    else:
        user_id = ''
        user_role = ''
        user_auth = True if 'tutorial' in nomenclature_server else False

    # POST requests headers
    headers_post = ct.HEADERS_POST_JSON
    headers_post['Authorization'] = token
    headers_post['user_id'] = user_id
    # POST headers to send binary data
    headers_post_bytes = ct.HEADERS_POST
    headers_post_bytes['Authorization'] = token
    headers_post_bytes['user_id'] = user_id

    schema_params = pv.read_configs(schema_directory, ct.SCHEMA_CONFIG_BASENAME)

    # verify that local configs have a single value per parameter
    if all([len(schema_params[k]) == 1
            for k in schema_params
            if k != 'chewBBACA_version']) is not True:
        sys.exit('Cannot sync schema with multiple values per parameter.')

    # check if schema exists in Chewie-NS
    schema_name, ns_params = cr.get_species_schemas(schema_id,
                                                    species_id,
                                                    nomenclature_server,
                                                    headers_get)[2:]

    # verify that local configs match NS configs
    if all([str(schema_params[k][0]) == ns_params[k]['value']
            for k in schema_params
            if k not in ['chewBBACA_version', 'window_size']]) is not True:
        sys.exit('Local configs do not match Chewie-NS configs.')

    # Get the name of the species from the provided id
    # or vice-versa
    species_id, species_name = cr.species_ids(species_id, nomenclature_server,
                                              headers_get)

    print('Schema id: {0}'.format(schema_id))
    print('Schema name: {0}'.format(schema_name))
    print("Schema's species: {0} (id={1})".format(species_name, species_id))
    print('Last synced: {0}'.format(local_date))

    # get last modification date
    # setting syncing date to last modification date will allow
    # all users to sync even when the schema is locked and being
    # updated by another user
    ns_date = ns_params['last_modified']['value']
    print('\nRemote schema was last modified on: {0}'.format(ns_date))

    # exit if remote schema has not been updated since last
    # sync date and current user does not wish to submit new alleles
    if local_date == ns_date and submit is False:
        sys.exit('\nRemote schema has not been updated since last sync '
                 'process. Local schema is up-to-date.')

    # Create a temporary dir for the new alleles
    temp_dir = fo.join_paths(os.path.dirname(schema_directory), ['temp'])
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # retrieve alleles added to schema after last sync date
    print('\nRetrieving alleles added to remote schema '
          'after {0}...'.format(local_date))
    loci_alleles, server_time, count = retrieve_latest(local_date, schema_uri,
                                                       headers_get, ns_date)

    print('Retrieved {0} alleles for {1} loci.'
          ''.format(count, len(loci_alleles)))

    # Get schema files from genes list file
    genes_list = fo.join_paths(schema_directory, [ct.GENE_LIST_BASENAME])
    genes = fo.pickle_loader(genes_list)

    # update loci structure
    not_in_ns, pickled_loci, rearranged = update_loci_files(loci_alleles, genes,
                                                            schema_directory, temp_dir)

    total_local = sum([len(v[0]) for k, v in not_in_ns.items()])
    print('Local schema has {0} novel alleles for {1} '
          'loci.'.format(total_local, len(not_in_ns)))

    # check if there are any changes to make
    if len(pickled_loci) == 0:
        shutil.rmtree(temp_dir)
        sys.exit('Remote schema has not been altered and local schema '
                 'does not have novel alleles.')

    results = {}
    attributed = 0
    if submit is True and user_auth is True and len(not_in_ns) > 0:

        # attempt to lock schema
        lock_res = cr.simple_post_request(nomenclature_server, headers_post,
                                          ['species', species_id,
                                           'schemas', schema_id,
                                           'lock'],
                                          data=json.dumps({'action': 'lock'}))[1]
        # if schema is already locked user cannot send alleles
        lock_status = lock_res.status_code
        if lock_status == 403:
            print('Schema is already locked. Another user might be updating '
                  'the schema. Please repeat the syncing process after a '
                  'while to add your new alleles to the Chewie-NS.\n The '
                  'process will now update your local schema with the alleles '
                  'retrieved from the Chewie-NS.')
        else:

            # after locking, check if date matches ns_date
            date_res = cr.simple_get_request(nomenclature_server, headers_get,
                                             ['species', species_id, 'schemas',
                                              schema_id, 'modified'])[1]

            date_value = (date_res.json()).split(' ')[-1]

            if date_value != ns_date:
                print('Data retrieved from the Chewie-NS has an older '
                      'timestamp than current schema timestamp. Schema '
                      'might have been updated before this syncing process. '
                      'Please repeat the syncing process in order to add '
                      'your new alleles to the schema. The process will now '
                      'update your local schema with the alleles retrieved '
                      'from the Chewie-NS.')

                # unlock schema
                lock_res = cr.simple_post_request(nomenclature_server, headers_post,
                                                  ['species', species_id, 'schemas',
                                                   schema_id, 'lock'],
                                                  data=json.dumps({'action': 'unlock'}))[1]
            else:
                print('Collecting data and creating files to submit local alleles...')
                # get list of loci for schema in the NS
                loci_res = cr.simple_get_request(nomenclature_server, headers_get,
                                                 ['species', species_id, 'schemas',
                                                  schema_id, 'loci'])[1]
                # get loci files names from response
                for l in loci_res.json()['Loci']:
                    locus_name = l['name']['value'] + '.fasta'
                    locus_uri = l['locus']['value']
                    if locus_name in not_in_ns:
                        not_in_ns[locus_name].append(locus_uri)

                # create files with length values to update
                length_files = create_lengths_files(not_in_ns, temp_dir)

                # create new alleles data
                alleles_files, \
                    loci_ids, \
                    loci_names = create_alleles_files(not_in_ns, nomenclature_server,
                                                      user_id, species_name,
                                                      species_id, schema_id,
                                                      temp_dir)

                # compress files with new alleles
                zipped_files = ['{0}.zip'.format(file) for file in alleles_files]
                list(map(fo.file_zipper, alleles_files, zipped_files))
                alleles_data = list(zip(zipped_files, loci_ids, loci_names))

                print('Sending and inserting new alleles...')
                failed, \
                    start_count = upload_alleles_data(alleles_data, length_files,
                                                      nomenclature_server, headers_post,
                                                      headers_post_bytes, species_id,
                                                      schema_id)

                # track progress through endpoint
                # set time limit for task completion (seconds)
                print()
                time_limit = 2100
                current_time = 0
                status = 'Updating'
                start_count = int(start_count.json()['nr_alleles'])
                while status != 'Complete' and (current_time < time_limit):
                    insertion_status = cr.simple_get_request(nomenclature_server,
                                                             headers_get,
                                                             ['species', species_id,
                                                              'schemas', schema_id,
                                                              'loci', 'locus', 'update'])[1]
                    insertion_status = insertion_status.json()
                    if 'message' in insertion_status:
                        status = 'Complete'
                        results = insertion_status['identifiers']

                    current_count = int(insertion_status['nr_alleles'])

                    inserted = current_count - start_count
                    print('\r', '    Inserted {0} alleles.'.format(inserted), end='')
                    time.sleep(2)
                    current_time += 2

                if current_time != time_limit:
                    # determine alleles that were attributed an identifier
                    repeated = sum([len(r[0]) for l, r in results.items()])
                    attributed = sum([len(r[1]) for l, r in results.items()])

                    print('\nThe Chewie-NS inserted {0} new alleles and detected '
                          '{1} repeated alleles.'.format(attributed, repeated))
                else:
                    print('\nCould not retrieve allele identifiers assigned by '
                          'Chewie-NS. Will adapt schema with retrieved alleles. '
                          'Please repeat the syncing process in order to assign '
                          'the new identifiers for the submitted alleles.')

                # remove files in temp folder
                fo.remove_files(length_files)
                fo.remove_files(alleles_files)
                fo.remove_files(zipped_files)

    # Change pickled files to FASTA files
    for locus, pick in pickled_loci.items():
        rearranged = pickle_to_fasta(locus, pick, temp_dir, results,
                                     rearranged)

    # update pre-computed tables
    # pre_computed_dir = fo.join_paths(schema_directory, ['pre_computed'])
    # if os.path.isdir(pre_computed_dir):
    #     schema_loci_fullpath = [fo.join_paths(schema_directory, [file])
    #                             for file in genes]
    #     loci_to_call = {file: schema_loci_fullpath.index(file)
    #                     for file in schema_loci_fullpath}
    #     update_hash_tables(rearranged, loci_files, loci_to_call,
    #                        schema_params['translation_table'][0],
    #                        pre_computed_dir)
    # Delete pre-computed directory
    # Future versions will change the pre-computed hash tables
    pre_computed_dir = fo.join_paths(schema_directory, ['pre_computed'])
    if os.path.isdir(pre_computed_dir):
        fo.delete_directory(pre_computed_dir)

    # change identifiers in SQLite DB
    # if len(rearranged) > 0 and update_profiles is True:
    #     print('\nUpdating local allele identifiers...')
    #     altered = ps.update_profiles(schema_directory, rearranged)
    #     if altered is not None:
    #         print('Updated {0} profiles.\n'.format(altered))
    #     else:
    #         print('Could not find local SQLite database to upload profiles.\n')

    # Re-determine the representative sequences
    schema_short_directory = fo.join_paths(schema_directory, ['short'])
    if attributed > 0 or count > 0:
        loci_list = fo.join_paths(schema_directory, [ct.LOCI_LIST])
        loci_list, total_loci = pv.check_input_type(temp_dir, loci_list)
        adapt_schema.main(loci_list, [schema_directory, schema_short_directory],
                          cpu_cores, float(schema_params['bsr'][0]),
                          0, 11, None, blast_path)

        # Delete invalid alleles and genes files
        parent_dir = os.path.dirname(schema_directory)
        files = [fo.join_paths(parent_dir, [file])
                 for file in os.listdir(parent_dir)
                 if 'invalid' in file]

        fo.remove_files(files)

        # get last modification date
        last_modified = cr.simple_get_request(nomenclature_server, headers_get,
                                              ['species', species_id,
                                               'schemas', schema_id,
                                               'modified'])[1]
        last_modified = (last_modified.json()).split(' ')[-1]
        server_time = last_modified

        # update NS config file with latest server time
        ns_configs = fo.join_paths(schema_directory, ['.ns_config'])
        fo.pickle_dumper([server_time, schema_uri], ns_configs)

        # Delete file with list of loci to adapt
        os.remove(loci_list)

    print('Received {0} new alleles for {1} loci and sent '
          '{2} for {3} loci. '.format(count, len(pickled_loci),
                                      attributed, len(not_in_ns)))

    # delete temp directory
    shutil.rmtree(temp_dir)

    # delete pre-computed loci mode values
    mode_file = fo.join_paths(schema_directory, ['loci_modes'])
    fo.remove_files([mode_file])

    # delete pre-computed BSR values from 'short' directory
    # representatives might have changed and BSR values are outdated
    bsr_file = fo.join_paths(schema_directory, ['short', 'self_scores'])
    fo.remove_files([bsr_file])
