#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the download of chewBBACA's schemas from the
Chewie-NS.

The process enables the download of ZIP archives that contain ready-to-use
versions of any schema in the Chewie-NS. It also allows users to download
any schema with the structure it had at a specific time point.



Expected input
--------------

The process expects the following variables wheter through command line
execution or invocation of the :py:func:`main` function:

- ``-sc``, ``schema_id`` : The schema identifier in the Chewie-NS.

    - e.g.: ``1``

- ``-sp``, ``species_id`` : The integer identifier or name of the species
  that the schema will be associated to in the Chewie-NS.

    - e.g.: ``1`` or ``'Yersinia pestis'``

- ``-o``, ``download_folder`` : Path to the parent directory of the folder
  that will store the downloaded schema. The process will create a folder
  with the schema's name inside the directory specified through this argument.

    - e.g.: ``/home/user/chewie_schemas``

- ``--cpu``, ``cpu_cores`` : Number of CPU cores that will be used to
  construct the schema if the process downloads FASTA files instead of
  the compressed version.

    - e.g.: ``4``

- ``--ns_url``, ``nomenclature_server_url`` : The base URL for the
  Nomenclature Server.

    - e.g.: ``http://127.0.0.1:5000/NS/api/`` (local host)

- ``--d``, ``date`` : Download schema with structure it had at
  specified date. Must be in the format "Y-m-dTH:M:S" or "Y-m-dTH:M:S.f".

    - e.g.: ``2020-03-27T11:38:00`` or ``2020-03-27T11:38:01.100``

- ``--latest`` : If the compressed version that is available is not the
  latest, downloads all loci and constructs schema locally.

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


def determine_upload(ns_loci, schema_path, temp_path,
                     base_url, headers_get):
    """
    """

    complete = []
    local_uniq = {}
    total_seqs = 0
    for locus, url in ns_loci.items():

        locus_id = locus.rstrip('.fasta')
        pickled_file = os.path.join(temp_path, '{0}_pickled'.format(locus_id))
        with open(pickled_file, 'rb') as pf:
            locus_sequences = pickle.load(pf)

        local_uniq[locus] = []
        local_uniq[locus].append(url)

        # after sync, the allele identifiers of the sequences
        # exclusive to the local schema should be ok to use
        # as identifiers for new alleles in the NS
        local_sequences = {k: v for k, v in locus_sequences.items()
                           if '*' in v[0]}
        for seqid, rec in local_sequences.items():
            local_uniq[locus].append([seqid, rec[1], len(rec[1])])
            # delete entry with '*'
            locus_sequences[seqid] = (str(seqid), rec[1])
            total_seqs += 1

        # save updated pickle
        os.remove(pickled_file)
        with open(pickled_file, 'wb') as pf:
            pickle.dump(locus_sequences, pf)

        complete.append(locus)

    print('Found a total of {0} new sequences in'
          ' local schema.'.format(total_seqs))

    return [local_uniq, complete]


def retrieve_alleles(loci_new_alleles, server_time, schema_uri,
                     count, headers_get, ns_date):
    """
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
    print('Got {0}'.format(len(new_alleles)))
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


def update_loci_files(loci, local_loci, schema_dir, temp_dir):
    """
    """

    # Check if new alleles are already on schema
    not_in_ns = {}
    pickled_loci = {}
    for gene in local_loci:
        if gene in loci:
            locus = gene
            alleles = loci[locus]

            not_in_ns[locus] = []

            # get locus file identifier without '.fasta'
            locus_id = locus.rstrip('.fasta')
            # get latest locus alleles retrieved from the NS
            locus_ns_seqs = {seq: seqid
                             for seqid, seq in alleles.items()}

            # full path for local locus file
            locus_file = os.path.join(schema_dir, locus)
            # full path for temporary file that will
            # substitute the current file
            temp_file = os.path.join(temp_dir, '{0}_pickled'.format(locus_id))

            # get local locus sequences
            with open(locus_file, 'r') as sf:
                # read sequences in local file
                records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                           for rec in SeqIO.parse(sf, 'fasta')}

            # check if the NS and local schema have the same set of
            # sequence identifiers
            ns_ids = set(list(locus_ns_seqs.values()))
            local_ids = set([rec[0] for seqid, rec in records.items()])
            # if the set of identifiers is equal, do not create a pickled file

            if ns_ids == local_ids:
                del not_in_ns[locus]
                continue

            # deepcopy local records to alter and sync with the NS locus
            altered_locus = deepcopy(records)
            id_map = {}
            for seqid, seq in altered_locus.items():
                if '*' in seq[0]:
                    id_map[int(seq[0][1:])] = seqid
                else:
                    id_map[int(seq[0])] = seqid

            # determine max integer identifier to later increment
            #max_id = max(id_map.keys())
            max_id = max([int(i) for i in ns_ids])

            # initalize list
            switched = []
            pre = 0
            # for each local record
            for seqid, seq in records.items():
                # get full sequence idetifier
                current_seqid = seqid
                # get integer identifier
                allele_id = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
                # get DNA sequence
                current_seq = seq[1]
                # if the DNA sequence is in the NS but was added to the
                # local schema through allele call
                if current_seq in locus_ns_seqs and '*' in current_seqid:
                    # get allele identifier in the NS
                    ns_id = int(locus_ns_seqs[current_seq])
                    # allele has the same integer identifier
                    # in local and NS schemas
                    if ns_id == allele_id:
                        # simply remove '*'
                        new_id = current_seqid.replace('*', '')
                        # add record without '*'
                        altered_locus[new_id] = [str(allele_id), seq[1]]
                        # delete record with '*'
                        del altered_locus[current_seqid]
                        # delete from latest NS sequences
                        del locus_ns_seqs[current_seq]
                    # allele may have a different identifier between the local
                    # and NS schemas
                    elif ns_id != allele_id:
                        # the same allele identifier might have been attributed
                        # to different allele sequences
                        if ns_id in id_map:
                            # NS identifier will be the new local identifier
                            new_id = '{0}_{1}'.format(locus.rstrip('.fasta'),
                                                      ns_id)
                            # the local allele identifier will replace the
                            # identifier of the local sequence that currently has
                            # the same identifier as the NS allele
                            old_id = current_seqid
                            # where is the '*'
                            old_seq = altered_locus['{0}_*{1}'.format(locus_id, ns_id)][1]
                            altered_locus[new_id] = [str(ns_id), current_seq]
                            switched.append(current_seq)
                            if old_seq not in switched:
                                altered_locus[old_id] = [allele_id, old_seq]

                            del altered_locus['{0}_*{1}'.format(locus_id, ns_id)]
                            del locus_ns_seqs[current_seq]
                        # the identifier attributed to the allele in the NS
                        # might be greater than the number of alleles in the
                        # local schema
                        elif ns_id not in id_map:
                            new_id = '{0}_{1}'.format(locus_id, ns_id)
                            altered_locus[new_id] = [str(ns_id), current_seq]
                            # simply delete the old record with the same allele
                            del altered_locus[current_seqid]
                            del locus_ns_seqs[current_seq]

                # local allele may not be in the NS
                elif current_seq not in locus_ns_seqs and '*' in current_seqid:
                    # local sequence may have an identifier
                    # that was attributed to another allele in the NS
                    if str(allele_id) in list(locus_ns_seqs.values()):
                        # attribute identifier greater than total number
                        # of alleles
                        max_id += 1
                        new_id = '{0}_*{1}'.format(locus_id, max_id)
                        altered_locus[new_id] = ['*{0}'.format(max_id), current_seq]

                        # get sequence in the NS that has that allele
                        ns_seq = [k for k, v in locus_ns_seqs.items() if int(v) == allele_id]
                        updated_seqid = current_seqid.replace('*', '')
                        # attribute the NS sequence to the allele identifier
                        altered_locus[updated_seqid] = [str(allele_id), ns_seq[0]]

                        # delete identifier with '*'
                        del altered_locus[current_seqid]
                        del locus_ns_seqs[ns_seq[0]]

                        # store info about sequences that are not in the NS
                        not_in_ns[locus].append(('{0}_*{1}'.format(locus_id,
                                                                   max_id),
                                                                   current_seq))
                    # local sequence has an identifier that is not attributed
                    # to any sequence in the NS
                    else:
                        not_in_ns[locus].append((current_seqid, current_seq))
                # local sequence was previously synced
                elif current_seq in locus_ns_seqs and '*' not in current_seqid:
                    del locus_ns_seqs[current_seq]
                    pre += 1

            # get local records that were synced in terms of identifiers
            updated_records = {}
            for seqid, seq in altered_locus.items():
                int_seqid = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
                updated_records[int_seqid] = (seq[0], seq[1])

            # check if there are any completely new sequences that are only present in the NS
            if len(locus_ns_seqs) > 0:
                ns_rest = {int(seqid): (seqid, seq) for seq, seqid in locus_ns_seqs.items()}
                updated_records = {**updated_records, **ns_rest}

            with open(temp_file, 'wb') as pl:
                pickle.dump(updated_records, pl)

            pickled_loci[locus] = temp_file
            if len(not_in_ns[locus]) == 0:
                del not_in_ns[locus]
        else:
            # there are no new alleles in the NS for this locus
            # check if local schema has new alleles
            locus_file = os.path.join(schema_dir, gene)
            # get local locus sequences
            with open(locus_file, 'r') as sf:
                # read sequences in local file
                records = {rec.id: [(rec.id).split('_')[-1], str(rec.seq)]
                           for rec in SeqIO.parse(sf, 'fasta')}

            # determine if there are sequences with '*'
            not_in_ns[gene] = [(v[0], v[1]) for k, v in records.items() if '*' in v[0]]

            if len(not_in_ns[gene]) > 0:
                updated_records = {}
                for seqid, seq in records.items():
                    int_seqid = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
                    updated_records[int_seqid] = (seq[0], seq[1])

                locus_id = gene.rstrip('.fasta')
                temp_file = os.path.join(temp_dir, '{0}_pickled'.format(locus_id))
                with open(temp_file, 'wb') as pl:
                    pickle.dump(updated_records, pl)

                pickled_loci[gene] = temp_file
            else:
                del not_in_ns[gene]

    return [not_in_ns, pickled_loci]


def retrieve_latest(local_date, schema_uri, schema_dir, temp_dir,
                    ns_url, headers_get, ns_date):
    """
    """

    # get list of sequences that are new considering the last date
    count = 0
    server_time = local_date
    new_alleles = defaultdict(dict)

    # get new alleles from Chewie-NS, maximum of 10k at a time
    print('Determining if there are new sequences in the NS...')
    print('Getting new alleles added to the NS since '
          '{0}'.format(str(server_time)))
    while server_time != ns_date:
        new_alleles, server_time, count = retrieve_alleles(new_alleles,
                                                           server_time,
                                                           schema_uri,
                                                           count,
                                                           headers_get,
                                                           ns_date)

    print('Retrieved {0} new alleles added since '
          '{1}'.format(count, local_date))

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

    token = aux.capture_login_credentials(base_url)

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
    print('Started at: {0}\n'.format(start_date_str))

    # GET request headers
    headers_get = cnst.HEADERS_GET_JSON
    headers_get['Authorization'] = token

    # determine current user ID and Role
    user_id, user_role, user_auth = aux.user_info(base_url, headers_get)
    print('User id: {0}'.format(user_id))
    print('User role: {0}'.format(user_role))

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
    species_info = aux.species_ids(species_id, base_url, headers_get)
    species_id, species_name = species_info
    print('\nNS species with identifier {0} '
          'is {1}.'.format(species_id, species_name))

    # get last modification date
    # setting syncing date to last modification date will allow
    # all users to sync even when the schema is locked and being
    # updated by another user
    ns_date = ns_params['last_modified']['value']

    # Create a temporary dir for the new alleles
    temp_dir = os.path.join(os.path.dirname(schema_dir), 'temp')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # retrieve alleles added to schema after last sync date
    loci_alleles, server_time, count = retrieve_latest(local_date, schema_uri,
                                                       schema_dir, temp_dir,
                                                       base_url, headers_get,
                                                       ns_date)

    # Get schema files from genes list file
    genes_list = os.path.join(schema_dir, '.genes_list')
    with open(genes_list, 'rb') as gl:
        genes = pickle.load(gl)

    # update loci structure
    not_in_ns, pickled_loci = update_loci_files(loci_alleles,
                                                genes,
                                                schema_dir,
                                                temp_dir)

    # check if there are any changes to make
    if len(pickled_loci) == 0:
        shutil.rmtree(temp_dir)
        sys.exit('Retrieved alleles are common to local and NS schema. '
                 'Local schema is up to date.')

    num_sent = 0
    if submit is True and user_auth is True and len(not_in_ns) > 0:

        # attempt to lock schema
        lock_res = aux.simple_post_request(base_url, headers_post,
                                           ['species', species_id,
                                            'schemas', schema_id,
                                            'lock'], {'action': 'lock'})
        # if schema is already locked user cannot send alleles
        lock_status = lock_res.status_code
        if lock_status == 403:
            print('Schema is already locked. Another user '
                  'might be updating the schema. Please repeat '
                  'the syncing process after a while to add your '
                  'new alleles to the Chewie-NS.\n The process '
                  'will now update your local schema with the '
                  'alleles retrieved from the Chewie-NS.')
        else:

            # after locking, check if date matches ns_date
            date_res = aux.simple_get_request(base_url, headers_get,
                                              ['species', species_id,
                                               'schemas', schema_id,
                                               'modified'])

            date_value = (date_res.json()).split(' ')[-1]

            if date_value != ns_date:
                print('Data retrieved from the Chewie-NS has an '
                      'older timestamp than current schema timestamp. '
                      'Schema might have been updated before this '
                      'syncing process. Please repeat the syncing '
                      'process in order to add your new alleles to the '
                      'schema. The process will now update your '
                      'local schema with the alleles retrieved '
                      'from the Chewie-NS.')

                # unlock schema
                lock_res = aux.simple_post_request(base_url, headers_post,
                                                   ['species', species_id,
                                                    'schemas', schema_id,
                                                    'lock'],
                                                   {'action': 'unlock'})
            else:

                print('Determining local alleles to submit...')

                # compare list of genes, if they do not intersect, halt process
                # get list of loci for schema in the NS
                loci_res = aux.simple_get_request(base_url, headers_get,
                                                  ['species', species_id,
                                                   'schemas', schema_id,
                                                   'loci'])
                # get loci files names from response
                ns_loci = {}
                for l in loci_res.json()['Loci']:
                    locus_name = l['name']['value'] + '.fasta'
                    locus_uri = l['locus']['value']
                    if locus_name in not_in_ns:
                        ns_loci[locus_name] = locus_uri

                # determine sequences that are not in the NS
                # we synced before so we just need to go to each file and
                # identify the alleles with '*' as the alleles to be added!
                # we have to submit them and alter locally too...
                upload, completed = determine_upload(ns_loci, schema_dir,
                                                     temp_dir, base_url,
                                                     headers_get)

                # create files with length values to send!
                # create an endpooint to send new alleles, it receives files
                # and calls a script that inserts the new alleles and updates
                # the pre-computed files!
                length_files = []
                for locus, recs in upload.items():
                    lengths = {locus: {hashlib.sha256(rec[1].encode('utf-8')).hexdigest(): rec[2]
                               for rec in recs[1:]}}
                    lengths_file = os.path.join(temp_dir,
                                                '{0}_lengths'.format(locus.rstrip('.fasta')))

                    aux.pickle_dumper(lengths_file, lengths)
                    length_files.append(lengths_file)

                # create new alleles data
                loci_ids = []
                loci_names = []
                alleles_files = []
                user_uri = '{0}users/{1}'.format(base_url, user_id)
                for locus, recs in upload.items():
                    locus_uri = recs[0]
                    alleles_sequences = [r[1] for r in recs[1:]]

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

                # compress files with new alleles
                zipped_files = ['{0}.zip'.format(file) for file in alleles_files]
                list(map(aux.file_zipper, alleles_files, zipped_files))
                alleles_data = list(zip(zipped_files, loci_ids, loci_names))

                # send data
                failed = []
                uploaded = 0
                print('Sending and inserting new alleles...')
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
                    elif length_status in [200, 201] and zip_status in [200, 201]:
                        uploaded += 1

                # get last modification date
                last_modified = aux.simple_get_request(base_url, headers_get,
                                                       ['species', species_id,
                                                        'schemas', schema_id,
                                                        'modified'])
                last_modified = (last_modified.json()).split(' ')[-1]
                server_time = last_modified

                num_sent = sum([len(v) for k, v in not_in_ns.items()])

    # change profiles in the master file so that local
    # profiles have updated allele identifiers

                # remove files in temp folder
                aux.remove_files(length_files)
                aux.remove_files(alleles_files)
                aux.remove_files(zipped_files)

    # change pickled files to FASTA files
    for locus, pick in pickled_loci.items():
        locus_id = locus.rstrip('.fasta')
        pickled_file = pick
        with open(pickled_file, 'rb') as pf:
            locus_sequences = pickle.load(pf)

        natsorted_locus = sorted(locus_sequences)

        fasta_path = os.path.join(temp_dir, locus)
        records = ['>{0}_{1}\n{2}'.format(locus_id,
                                          locus_sequences[seqid][0],
                                          locus_sequences[seqid][1])
                   for seqid in natsorted_locus]

        fasta_text = '\n'.join(records)

        with open(fasta_path, 'w') as fp:
            fp.write(fasta_text)

        os.remove(pickled_file)

    # Re-determine the representative sequences
    if num_sent > 0 or count > 0:
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
                                      num_sent, len(not_in_ns)))

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
