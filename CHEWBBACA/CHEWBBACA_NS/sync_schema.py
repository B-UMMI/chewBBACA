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

Expected input
--------------

The process expects the following variables wheter through command line
execution or invocation of the :py:func:`main` function:

- ``-sc``, ``schema_directory`` : Path to the directory with the schema to be
  synced.

    - e.g.: ``/home/user/chewie_schemas/ecoli_schema``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores that will be used to
  determine new representatives if the process downloads new alleles from
  the Chewie-NS.

    - e.g.: ``4``

- ``--ns_url``, ``nomenclature_server_url`` : The base URL for the
  Nomenclature Server.

    - e.g.: ``http://127.0.0.1:5000/NS/api/`` (localhost)

- ``--submit`` : If the process should identify new alleles in the local
  schema and send them to the Chewie-NS (only registered users can submit
  new alleles).

Code documentation
------------------
"""


import os
import sys
import json
import shutil
import pickle
import hashlib
import argparse
import requests
import datetime as dt
from copy import deepcopy
from collections import defaultdict
from SPARQLWrapper import SPARQLWrapper
from urllib3.exceptions import InsecureRequestWarning

from Bio import SeqIO

try:
    from utils import constants as cnst
    from utils import auxiliary_functions as aux
    from PrepExternalSchema import PrepExternalSchema
except:
    from CHEWBBACA.utils import constants as cnst
    from CHEWBBACA.utils import auxiliary_functions as aux
    from CHEWBBACA.PrepExternalSchema import PrepExternalSchema


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


uniprot_sparql = SPARQLWrapper(cnst.UNIPROT_SPARQL)


def create_lengths_files(upload, temp_dir):
    """ Creates files with the length values of the sequences
        that will be sent to the Chewie-NS.

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

        lengths_file = os.path.join(temp_dir, '{0}_lengths'.format(locus_id))
        aux.pickle_dumper(lengths_file, lengths)
        length_files.append(lengths_file)

    return length_files


def create_alleles_files(upload, base_url, user_id, species_name,
                         species_id, schema_id, temp_dir):
    """ Creates files with the essential data to insert alleles in the NS.

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

        alleles_file = os.path.join(temp_dir,
                                    '{0}_{1}_{2}'.format(species_id,
                                                         schema_id,
                                                         locus_id))
        alleles_files.append(alleles_file)

        aux.pickle_dumper(alleles_file, post_inputs)

    return [alleles_files, loci_ids, loci_names]


def upload_alleles_data(alleles_data, length_files, base_url, headers_post,
                        headers_post_bytes, species_id, schema_id):
    """ Uploads files with the data to insert alleles and the
        length values for the sequences of each locus.

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

    failed = []
    for i, a in enumerate(alleles_data):

        locus_id = a[1]

        # get length of alleles from current locus
        current_len = length_files[i]
        data = aux.pickle_loader(current_len)
        data = {locus_id: data[next(iter(data))]}
        data = json.dumps({'content': data})

        # send data to the NS
        send_url = aux.make_url(base_url, 'species', species_id,
                                'schemas', schema_id, 'loci',
                                locus_id, 'lengths')

        lengths_res = aux.upload_data(data, send_url, headers_post, False)
        length_status = lengths_res.status_code

        # get path to ZIP archive with data to insert alleles
        current_zip = a[0]

        # send data to insert alleles in the NS
        zip_url = aux.make_url(base_url, 'species', species_id,
                               'schemas', schema_id, 'loci',
                               locus_id, 'update')

        if alleles_data[i] == alleles_data[-1]:
            headers_post_bytes['complete'] = 'True'

        zip_res = aux.upload_file(current_zip, os.path.basename(current_zip),
                                  zip_url, headers_post_bytes,
                                  False)
        zip_status = zip_res.status_code

        # determine if upload was successful
        if length_status not in [200, 201] or zip_status not in [200, 201]:
            failed.append(locus_id)

    return [failed, zip_res.json()]


def pickle_to_fasta(locus, pickled_file, temp_dir, identifiers):
    """ Creates FASTA files with the information contained in
        a pickled file.

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
            The `zip_res` variable returned by the :py:func:`upload_alleles_data`
            function. It will be used to change allele identifiers that were
            successfully inserted into the Cewie-NS.

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

    locus_sequences = aux.pickle_loader(pickled_file)

    natsorted_locus = sorted(locus_sequences)

    fasta_path = os.path.join(temp_dir, locus)
    records = []
    for seqid in natsorted_locus:
        recid = locus_sequences[seqid][0]
        seq = locus_sequences[seqid][1]
        seq_hash = hashlib.sha256(seq.encode('utf-8')).hexdigest()
        # switch by the identifier attributed by the Chewie-NS
        if seq_hash in attributed:
            recid = attributed[seq_hash]
        elif seq_hash in repeated:
            recid = repeated[seq_hash]

        record = '>{0}_{1}\n{2}'.format(locus_id, recid, seq)
        records.append(record)

    fasta_text = '\n'.join(records)

    with open(fasta_path, 'w') as fp:
        fp.write(fasta_text)

    os.remove(pickled_file)

    return fasta_path


def retrieve_alleles(loci_new_alleles, server_time, schema_uri,
                     count, headers_get, ns_date):
    """ Retrieves alleles added to a schema in the Chewie-NS
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
    url = aux.make_url(schema_uri, 'loci')
    payload = {'local_date': server_time, 'ns_date': ns_date}
    # get the new alleles
    response = requests.get(url, headers=headers_get,
                            timeout=30, params=payload)

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


def process_common(ns_id, local_id, local_seqid, seq, local_seqs,
                   ns_seqs, id_map, locus_id, switched):
    """ Determines the integer identifier that must be attributed
        to a local allele that has been added to the schema through
        allele calling and is in the Chewie-NS.

        Parameters
        ----------
        ns_id : int
            The integer identifier of the allele in the remote schema.
        local_id : int
            The integer identifier of the allele in the local schema.
        local_seqid : str
            The full identifier (locus prefix plus integer identifier)
            of the allele in the local schema.
        seq : str
            The DNA sequence of the allele.
        local_seqs : dict
            A dictionary with alleles full identifiers as keys and lists,
            that contain integer identifiers and DNA sequences, as values.
        ns_seqs : dict
            A dictionary with DNA sequences as keys and alleles full
            identifiers as values.
        id_map : dict
            A dictionary with alleles integer identifiers as keys and
            full alleles identifiers as values.
        locus_id : str
            The locus prefix.
        switched : list
            A list with all the allele identifiers reassignments
            performed for the sequences of the current locus.

        Returns
        -------
        A list with three elements:

        - local_seqs : dict
            Input `local_seqs` with the allele identifier for
            the current allele updated (as well as the identifier
            of any other allele that had to be switched to be able
            to assign the correct identifier to the current allele).
        - ns_seqs : dict
            Input `ns_seqs` without the entry for the sequence
            that had its allele identifier evaluated.
        - switched : dict
            Input `switched` variable with a new entry for the
            reassignment performed.
    """

    # invert reassignments dict
    inv_switched = {v: k for k, v in switched.items()}
    # local and remote integer identifiers are equal
    if ns_id == local_id:
        # simply remove '*' from identifier
        new_id = local_seqid.replace('*', '')
        local_seqs[new_id] = [str(local_id), seq]

        if local_seqid in inv_switched:
            switched[inv_switched[local_seqid]] = new_id
        else:
            switched[local_seqid] = new_id

        # delete record with '*' and record retrieved from the NS
        del local_seqs[local_seqid]
        del ns_seqs[seq]
    # local and remote integer identifiers differ
    elif ns_id != local_id:
        # identifier in the NS was assigned to other sequence
        if ns_id in id_map:
            # NS identifier will be the new identifier
            new_id = '{0}_{1}'.format(locus_id, ns_id)

            # current sequence identifier
            old_id = local_seqid
            # identifier of sequence with same identifier as NS
            beep_id = '{0}_*{1}'.format(locus_id, ns_id)
            if beep_id in local_seqs:
                old_seq = local_seqs[beep_id][1]

            local_seqs[new_id] = [str(ns_id), seq]
            if old_id in inv_switched:
                switched[inv_switched[old_id]] = new_id
            else:
                switched[old_id] = new_id

            if beep_id in local_seqs:
                local_seqs[old_id] = [str(local_id), old_seq]
                if beep_id in inv_switched:
                    switched[inv_switched[beep_id]] = old_id
                else:
                    switched[beep_id] = old_id
                del local_seqs[beep_id]
            else:
                del local_seqs[old_id]
            del ns_seqs[seq]
        # ns identifier is greater than maximum local identifier
        elif ns_id not in id_map:
            new_id = '{0}_{1}'.format(locus_id, ns_id)
            local_seqs[new_id] = [str(ns_id), seq]

            if local_seqid in inv_switched:
                switched[inv_switched[local_seqid]] = new_id
            else:
                switched[local_seqid] = new_id

            # simply delete the old record with the same allele
            del local_seqs[local_seqid]
            del ns_seqs[seq]

    return [local_seqs, ns_seqs, switched]


def process_novel(local_id, local_seqid, seq, local_seqs, ns_seqs,
                  max_id, locus_id, not_in_ns, locus, switched):
    """ Determines the integer identifier that must be attributed
        to a local allele that is not present in the Chewie-NS.

        Parameters
        ----------
        local_id : int
            The integer identifier of the allele in the local schema.
        local_seqid : str
            The full identifier (locus prefix plus integer identifier)
            of the allele in the local schema.
        seq : str
            The DNA sequence of the allele.
        local_seqs : dict
            A dictionary with alleles full identifiers as keys and lists,
            that contain integer identifiers and DNA sequences, as values.
        ns_seqs : dict
            A dictionary with DNA sequences as keys and alleles full
            identifiers as values.
        max_id : int
            The highest integer identifier attributed to local and remote
            alleles.
        locus_id : str
            The locus prefix.
        not_in_ns : dict
            Dictionary with information about local alleles that are not
            in the Chewie-NS. It has loci identifiers as keys and lists
            as values. Each list contains the following elements:

            - A dictionary with sequences hashes as keys and a list
              with the sequence identifier, DNA sequence and sequence
              length value.
            - The locus URI.
        locus : str
            Full locus identifier with '.fasta' suffix.

        Returns
        -------
        A list with four elements:

        - local_seqs : dict
            Input `local_seqs` variable with the allele identifier for
            the current allele updated (as well as the identifier
            of any other allele that had to be switched to be able
            to assign the correct identifier to the current allele).
        - ns_seqs : dict
            Input `ns_seqs` variable without the entry for the sequence
            that had its allele identifier evaluated.
        - not_in_ns : dict
            Input `not_in_ns` variable that was updated with the info
            of the allele that is not in the Chewie-NS.
        - max_id : int
            Updated value for the maximum integer identifier
            attributed to an allele.
    """

    seq_hash = hashlib.sha256(seq.encode('utf-8')).hexdigest()
    # local sequence may have an identifier
    # that was attributed to another allele in the NS
    if str(local_id) in list(ns_seqs.values()):
        # invert local_seqs
        inv_local = {v[1]: [k, v[0]] for k, v in local_seqs.items()}

        # attribute identifier greater than total number
        # of alleles
        max_id += 1
        new_id = '{0}_*{1}'.format(locus_id, max_id)
        local_seqs[new_id] = ['*{0}'.format(max_id), seq]

        # get sequence in the NS that has that allele identifier
        ns_seq = [k for k, v in ns_seqs.items() if int(v) == local_id]
        updated_seqid = local_seqid.replace('*', '')
        # attribute the NS sequence to the allele identifier
        local_seqs[updated_seqid] = [str(local_id), ns_seq[0]]
        switched[local_seqid] = new_id
        # check if sequence from NS is in local sequences
        if ns_seq[0] in inv_local:
            old_id = inv_local[ns_seq[0]][0]
            inv_switched = {v: k for k, v in switched.items()}
            if old_id in inv_switched:
                switched[inv_switched[old_id]] = updated_seqid
            else:
                switched[old_id] = updated_seqid
            del local_seqs[old_id]

        # delete identifier with '*'
        del local_seqs[local_seqid]
        del ns_seqs[ns_seq[0]]

        # store info about sequences that are not in the NS
        not_in_ns[locus][0][seq_hash] = ['{0}_*{1}'.format(locus_id, max_id),
                                         seq, len(seq)]
    # local sequence has an identifier that is not attributed
    # to any sequence in the NS
    else:
        not_in_ns[locus][0][seq_hash] = [local_seqid, seq, len(seq)]

    return [local_seqs, ns_seqs, not_in_ns, max_id, switched]


def alternative_alter(loci, schema_dir, pickled_loci, not_in_ns, temp_dir, rearranged):
    """
    """

    for locus, alleles in loci.items():

        locus_id = locus.rstrip('.fasta')

        # get latest alleles retrieved from the Chewie-NS
        ns_seqs = {seq: seqid for seqid, seq in alleles.items()}

        # paths for current and temp locus file
        locus_file = os.path.join(schema_dir, locus)

        # get local locus sequences
        records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                   for rec in SeqIO.parse(locus_file, 'fasta')}

        # check if the NS and local schema have the same set of
        # sequence identifiers
        ns_ids = set(list(ns_seqs.values()))
        local_ids = set([rec[0] for seqid, rec in records.items()])

        # proceed to next locus if set of identifiers is equa
        if ns_ids == local_ids:
            continue

        # invert local dict
        inv_local = {v[1]: [k, v[0]] for k, v in records.items()}

        # alter identifiers of local alleles that were added to the NS
        switched = {}
        for seq, seqid in ns_seqs.items():
            records[seqid] = [seqid.split('_')[-1], seq]
            if seq in inv_local:
                switched[inv_local[seq][0]] = seqid
                del records[inv_local[seq][0]]

        # identify alleles with '*' and move them to top
        max_id = max([int(v[0]) for k, v in records.items() if '*' not in k])

        novel_ids = [k for k in records if '*' in k]
        sorted_novel = sorted(novel_ids, key=lambda x: int(x.split('*')[-1]))

        for si in sorted_novel:
            max_id += 1
            new_id = '{0}*{1}'.format(si.split('*')[0], max_id)
            records[new_id] = [new_id.split('_')[-1], records[si][1]]
            if si != new_id:
                del records[si]
                switched[si] = new_id

        # determine records that are not in the NS
        not_in_ns[locus] = {hashlib.sha256(v[1].encode('utf-8')).hexdigest(): [v[0], v[1], len(v[1])]
                            for k, v in records.items() if '*' in v[0]}

        rearranged[locus] = switched

        updated_records = {}
        for seqid, seq in records.items():
            int_seqid = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
            updated_records[int_seqid] = (seq[0], seq[1])

        temp_file = os.path.join(temp_dir, '{0}_pickled'.format(locus_id))
        aux.pickle_dumper(temp_file, updated_records)

        pickled_loci[locus] = temp_file

    return [pickled_loci, not_in_ns, rearranged]


def altered_loci(loci, schema_dir, pickled_loci, not_in_ns, temp_dir, rearranged):
    """ Reassigns alleles identifiers of loci that were altered
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
            processed..
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

        # get latest alleles retrieved from the Chewie-NS
        ns_seqs = {seq: seqid for seqid, seq in alleles.items()}

        # paths for current and temp locus file
        locus_file = os.path.join(schema_dir, locus)
        temp_file = os.path.join(temp_dir, '{0}_pickled'.format(locus_id))

        # get local locus sequences
        records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                   for rec in SeqIO.parse(locus_file, 'fasta')}

        # check if the NS and local schema have the same set of
        # sequence identifiers
        ns_ids = set(list(ns_seqs.values()))
        local_ids = set([rec[0] for seqid, rec in records.items()])

        # proceed to next locus if set of identifiers is equa
        if ns_ids == local_ids:
            continue

        # deepcopy local records to alter and sync with the NS locus
        local_seqs = deepcopy(records)
        id_map = {}
        for seqid, seq in local_seqs.items():
            allele_id = seq[0]
            id_map[int(allele_id.replace('*', ''))] = seqid

        # determine max integer identifier to later increment
        identifiers = list(ns_ids) + list(id_map.keys())
        max_id = max([int(i) for i in identifiers])

        # initalize list
        pre = 0
        switched = {}
        not_in_ns[locus] = [{}]
        # for each local record
        for seqid, seq in records.items():
            inv_local = {v[1]: [k, v[0]] for k, v in local_seqs.items()}
            # get full and integer identifiers
            current_seqid = inv_local[seq[1]][0]
            #current_seqid = seqid
            allele_id = inv_local[seq[1]][1]
            allele_id = int(allele_id.replace('*', ''))
            #allele_id = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
            current_seq = seq[1]
            # local allele is in the Chewie-NS
            if current_seq in ns_seqs and '*' in current_seqid:
                # get allele identifier in the Chewie-NS
                ns_id = int(ns_seqs[current_seq])

                local_seqs, \
                    ns_seqs, \
                    switched = process_common(ns_id, allele_id, current_seqid,
                                              current_seq, local_seqs, ns_seqs,
                                              id_map, locus_id, switched)

            # local allele is not in the Chewie-NS
            elif current_seq not in ns_seqs and '*' in current_seqid:

                local_seqs, \
                    ns_seqs, \
                    not_in_ns, \
                    max_id, \
                    switched = process_novel(allele_id, current_seqid,
                                             current_seq, local_seqs,
                                             ns_seqs, max_id, locus_id,
                                             not_in_ns, locus, switched)

            # local allele was previously synced
            elif current_seq in ns_seqs and '*' not in current_seqid:
                del ns_seqs[current_seq]
                pre += 1

        if len(switched) > 0:
            rearranged[locus] = switched

        # get local records that were synced in terms of identifiers
        updated_records = {}
        for seqid, seq in local_seqs.items():
            int_seqid = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
            updated_records[int_seqid] = (seq[0], seq[1])

        # check if there are any completely new sequences that are only present in the NS
        if len(ns_seqs) > 0:
            ns_rest = {int(seqid): (seqid, seq) for seq, seqid in ns_seqs.items()}
            updated_records = {**updated_records, **ns_rest}

        aux.pickle_dumper(temp_file, updated_records)

        pickled_loci[locus] = temp_file

        if len(not_in_ns[locus][0]) == 0:
            del not_in_ns[locus]

    return [pickled_loci, not_in_ns, rearranged]


def unaltered_loci(loci, schema_dir, pickled_loci, not_in_ns, temp_dir):
    """ Determines if local schema has new alleles for loci that
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
        locus_file = os.path.join(schema_dir, gene)
        # get local locus sequences
        records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                   for rec in SeqIO.parse(locus_file, 'fasta')}

        # determine if there are sequences with '*'
        novel_local = {hashlib.sha256(v[1].encode('utf-8')).hexdigest(): [v[0], v[1], len(v[1])]
                       for k, v in records.items() if '*' in v[0]}
        if len(novel_local) > 0:
            not_in_ns[gene] = [novel_local]
            updated_records = {}
            for seqid, seq in records.items():
                int_seqid = int(seq[0].replace('*', ''))
                updated_records[int_seqid] = (seq[0], seq[1])

            temp_file = os.path.join(temp_dir, '{0}_pickled'.format(locus_id))
            aux.pickle_dumper(temp_file, updated_records)

            pickled_loci[gene] = temp_file

    return [pickled_loci, not_in_ns]

loci = loci_alleles
local_loci = genes


def update_loci_files(loci, local_loci, schema_dir, temp_dir):
    """ Determines which loci were or were not updated in the
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
    pickled_loci, not_in_ns, rearranged = altered_loci(updated, schema_dir,
                                                       pickled_loci, not_in_ns,
                                                       temp_dir, rearranged)

    # determine loci that were not altered in the Chewie-NS but
    # have been altered locally
    pickled_loci, not_in_ns = unaltered_loci(not_updated, schema_dir,
                                             pickled_loci, not_in_ns,
                                             temp_dir)

    return [not_in_ns, pickled_loci, updated, not_updated, rearranged]


def retrieve_latest(local_date, schema_uri, headers_get, ns_date):
    """ Retrieves alleles added to a schema in the Chewie-NS after
        a specified date.

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


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-sc', type=str, dest='schema_directory', required=True,
                        help='Path to the directory with the schema to be'
                             'synced.')

    parser.add_argument('--cpu', type=int, required=False,
                        dest='cpu_cores', default=1,
                        help='Number of CPU cores that will '
                             'be used to determine new representatives '
                             'if the process downloads new alleles from '
                             'the Chewie-NS.')

    parser.add_argument('--ns_url', type=str, required=False,
                        dest='nomenclature_server_url',
                        default=cnst.HOST_NS,
                        help='The base URL for the Nomenclature Server.')

    parser.add_argument('--submit', required=False,
                        action='store_true', dest='submit',
                        help='If the process should identify new alleles '
                             'in the local schema and send them to the '
                             'NS. (only users with permissions level of '
                             'Contributor can submit new alleles).')

    args = parser.parse_args()

    return [args.schema_directory, args.cpu_cores,
            args.nomenclature_server_url, args.submit]


def main(schema_dir, core_num, base_url, submit):

    if submit is True:
        print('\nOnly registered users may submit new alleles.')
        token = aux.capture_login_credentials(base_url)
    else:
        token = ''

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
    print('Started at: {0}\n'.format(start_date_str))

    # GET request headers
    headers_get = cnst.HEADERS_GET_JSON
    headers_get['Authorization'] = token

    # determine current user ID and Role
    if submit is True:
        user_id, user_role, user_auth = aux.user_info(base_url, headers_get)
        user_auth = True
        print('User id: {0}'.format(user_id))
        print('User role: {0}\n'.format(user_role))
    else:
        user_id = ''
        user_role = ''
        user_auth = False

    # POST requests headers
    headers_post = cnst.HEADERS_POST_JSON
    headers_post['Authorization'] = token
    headers_post['user_id'] = user_id
    # POST headers to send binary data
    headers_post_bytes = cnst.HEADERS_POST
    headers_post_bytes['Authorization'] = token
    headers_post_bytes['user_id'] = user_id

    # get schema configs
    local_date, schema_uri = aux.read_configs(schema_dir, '.ns_config')
    # get schema and species identifiers
    schema_id = schema_uri.split('/')[-1]
    species_id = schema_uri.split('/')[-3]

    schema_params = aux.read_configs(schema_dir, '.schema_config')

    # check if schema exists in the NS
    schema_name, ns_params = aux.get_species_schemas(schema_id,
                                                     species_id,
                                                     base_url,
                                                     headers_get)[2:]

    # Get the name of the species from the provided id
    # or vice-versa
    species_id, species_name = aux.species_ids(species_id, base_url, headers_get)

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
    temp_dir = os.path.join(os.path.dirname(schema_dir), 'temp')
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
    genes_list = os.path.join(schema_dir, '.genes_list')
    with open(genes_list, 'rb') as gl:
        genes = pickle.load(gl)

    # update loci structure
    not_in_ns, pickled_loci, \
        updated, not_update, \
        rearranged = update_loci_files(loci_alleles, genes,
                                       schema_dir, temp_dir)
    print(rearranged)
    total_local = sum([len(v) for k, v in not_in_ns.items()])
    print('Local schema has {0} novel alleles.'.format(total_local))

    # check if there are any changes to make
    if len(pickled_loci) == 0:
        shutil.rmtree(temp_dir)
        sys.exit('Remote schema has not been altered and local schema '
                 'does not have novel alleles.')

    results = {}
    attributed = 0
    if submit is True and user_auth is True and len(not_in_ns) > 0:

        # attempt to lock schema
        lock_res = aux.simple_post_request(base_url, headers_post,
                                           ['species', species_id,
                                            'schemas', schema_id,
                                            'lock'], {'action': 'lock'})
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
            date_res = aux.simple_get_request(base_url, headers_get,
                                              ['species', species_id,
                                               'schemas', schema_id,
                                               'modified'])

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
                lock_res = aux.simple_post_request(base_url, headers_post,
                                                   ['species', species_id,
                                                    'schemas', schema_id,
                                                    'lock'],
                                                   {'action': 'unlock'})
            else:
                print('Determining local alleles to submit...')
                # get list of loci for schema in the NS
                loci_res = aux.simple_get_request(base_url, headers_get,
                                                  ['species', species_id,
                                                   'schemas', schema_id,
                                                   'loci'])
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
                    loci_names = create_alleles_files(not_in_ns, base_url, user_id,
                                                      species_name, species_id,
                                                      schema_id, temp_dir)

                # compress files with new alleles
                zipped_files = ['{0}.zip'.format(file) for file in alleles_files]
                list(map(aux.file_zipper, alleles_files, zipped_files))
                alleles_data = list(zip(zipped_files, loci_ids, loci_names))

                print('Sending and inserting new alleles...')
                failed, results = upload_alleles_data(alleles_data, length_files, base_url,
                                                      headers_post, headers_post_bytes,
                                                      species_id, schema_id)

                # determine alleles that were attributed an identifier
                repeated = 0
                attributed = 0
                for locus, r in results.items():
                    repeated += len(r[0])
                    attributed += len(r[1])

                print('The Chewie-NS inserted {0} new alleles and detected '
                      '{1} repeated alleles.'.format(attributed, repeated))

                # get last modification date
                last_modified = aux.simple_get_request(base_url, headers_get,
                                                       ['species', species_id,
                                                        'schemas', schema_id,
                                                        'modified'])
                last_modified = (last_modified.json()).split(' ')[-1]
                server_time = last_modified

                # remove files in temp folder
                aux.remove_files(length_files)
                aux.remove_files(alleles_files)
                aux.remove_files(zipped_files)

    # change pickled files to FASTA files
    for locus, pick in pickled_loci.items():
        pickle_to_fasta(locus, pick, temp_dir, results)

    # Re-determine the representative sequences
    if attributed > 0 or count > 0:
        PrepExternalSchema.main(temp_dir, schema_dir,
                                core_num, float(schema_params['bsr'][0]),
                                int(schema_params['minimum_locus_length'][0]),
                                11, '', None)

        # delete invalid alleles and genes files
        parent_dir = os.path.dirname(schema_dir)
        files = [os.path.join(parent_dir, file)
                 for file in os.listdir(parent_dir)
                 if 'invalid' in file]

        for f in files:
            os.remove(f)

        # update NS config file with latest server time
        ns_configs = os.path.join(schema_dir, '.ns_config')
        with open(ns_configs, 'wb') as nc:
            pickle.dump([server_time, schema_uri], nc)

    print('Received {0} new alleles for {1} loci and sent '
          '{2} for {3} loci. '.format(count, len(pickled_loci),
                                      attributed, len(not_in_ns)))

    # delete temp directory
    shutil.rmtree(temp_dir)

    end_date = dt.datetime.now()
    end_date_str = dt.datetime.strftime(end_date, '%Y-%m-%dT%H:%M:%S')

    delta = end_date - start_date
    minutes, seconds = divmod(delta.total_seconds(), 60)

    print('\nFinished at: {0}'.format(end_date_str))
    print('Elapsed time: {0:.0f}m{1:.0f}s'.format(minutes, seconds))


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3])
