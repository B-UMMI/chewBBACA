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
import shutil
import pickle
import argparse
import requests
import datetime as dt
from copy import deepcopy
import concurrent.futures
from itertools import repeat
from collections import defaultdict
from SPARQLWrapper import SPARQLWrapper
from urllib3.exceptions import InsecureRequestWarning

from Bio import SeqIO

try:
    from utils import constants as cnst
    from utils import auxiliary_functions as aux
    from PrepExternalSchema import PrepExternalSchema
except ImportError:
    from CHEWBBACA.utils import constants as cnst
    from CHEWBBACA.utils import auxiliary_functions as aux
    from CHEWBBACA.PrepExternalSchema import PrepExternalSchema


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


uniprot_sparql = SPARQLWrapper(cnst.UNIPROT_SPARQL)


def determine_upload(local_schema_loci, ns_schema_loci,
                     ns_schema_locid_map, schema_path,
                     temp_path, base_url, headers_get):
    """
    """

    local_uniq = {}
    comp = []
    for locus in local_schema_loci:

        local_locus = locus
        locus_id = local_locus.rstrip('.fasta')
        pickled_file = os.path.join(temp_path, '{0}_pickled'.format(locus_id))
        with open(pickled_file, 'rb') as pf:
            locus_sequences = pickle.load(pf)
        local_sequences = local_schema_loci[locus]

        if local_locus in ns_schema_loci:
            local_uniq[locus] = []

            # get uniprot info for locus
            ns_uri = ns_schema_locid_map[local_locus][0]

            local_uniq[local_locus].append(ns_uri)

            # after sync, the allele identifiers of the sequences
            # exclusive to the local schema should be ok to use
            # as identifiers for new alleles in the NS
            for rec in local_sequences:
                seqid = rec[0]
                seqid = seqid.replace('*', '')
                seq = rec[1]

                local_uniq[local_locus].append([seqid,
                                                seq, len(seq)])

                # delete entry with '*'
                allele_id = int(seqid.split('_')[-1])
                locus_sequences[allele_id] = (str(allele_id), locus_sequences[allele_id][1])
            
            # save updated pickle
            os.remove(pickled_file)
            with open(pickled_file, 'wb') as pf:
                pickle.dump(locus_sequences, pf)
            
            comp.append(locus)

    print('Found a total of {0} new sequences in'
          ' local schema.'.format(len(local_uniq)))

    return [local_uniq, comp]


def build_fasta_files(new_allele_seq_dict, path2schema, path_new_alleles):
    """
    """
    
    for name, allele_dict in new_allele_seq_dict.items():
        
        for allele_id, seq in allele_dict.items():
            
            auxDict = {}
            auxname = name.split(".")
            listIndexAux = []

		# create a dictionary with the info from fasta to build the fasta file; {id:sequence}
        listIndexAux.append(int(allele_id))
        auxDict[int(allele_id)] = f"{seq}"

		# sort by allele id, create the fasta string and write the file
        listIndexAux.sort()
        auxString = ''
        for index in listIndexAux:
            auxString += f">{auxname[0]}_{index}\n{auxDict[index]}\n"
        
        # if the fasta already exists, add new allele
        if os.path.exists(os.path.join(path2schema, name)):
            
            with open(os.path.join(path2schema, name), 'a') as f:
                f.write(auxString)
            
            # Copy the updated file to new_alleles dir 
            new_alleles_path = shutil.copy(os.path.join(path2schema, name), path_new_alleles)
        
        # write new file to schema dir and new_alleles dir
        else:
            
            # Schema dir
            with open(os.path.join(path2schema, name), 'w') as f1:
                f1.write(auxString)

            # New alleles dir
            with open(os.path.join(path_new_alleles, name), 'w') as f2:
                f2.write(auxString)

    print(f"{len(new_allele_seq_dict)} files have been written")

    return True


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
    for locus, alleles in loci.items():
        not_in_ns[locus] = []
        if locus in local_loci:
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

        # get local records that were synced in terms of identifiers
        updated_records = {}
        for seqid, seq in altered_locus.items():
            int_seqid = int(seq[0]) if '*' not in seq[0] else int(seq[0][1:])
            updated_records[int_seqid] = (seq[0], seq[1])

        # check if there are any completely new sequences that are only present in the NS
        if len(locus_ns_seqs) > 0:
            ns_rest = {int(seqid): (seqid, seq) for seq, seqid in locus_ns_seqs.items()}
            updated_records = {**updated_records, **ns_rest}

        # sort by integer identifier
        #updated_records = sorted(updated_records, key=lambda x: x[0])

        with open(temp_file, 'wb') as pl:
            pickle.dump(updated_records, pl)

        pickled_loci[locus] = temp_file
        if len(not_in_ns[locus]) == 0:
            del not_in_ns[locus]

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

    parser.add_argument('-sc', type=str, dest='schema_dir', required=True,
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
                             'NS. (only users with permissons level of '
                             'Contributor can submit new alleles).')

    # add waiting time for the step that adds alleles to the NS? If schema is locked wait this time.

    args = parser.parse_args()

    return [args.schema_dir, args.cpu_cores,
            args.nomenclature_server_url, args.submit]


schema_dir = '/home/rfm/Desktop/ns_test/test_download/sagalactiae_sagalactiae3'
core_num = 6
base_url = 'http://127.0.0.1:5000/NS/api/'
submit = 'yes'


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

    # get last modification date
    # setting syncing date to last modification date will allow
    # all users to sync even when the schema is locked and being
    # updated by another user
    ns_date = ns_params['last_modified']['value']

    # Create a temporary dir for the new alleles
    temp_dir = os.path.join(os.path.dirname(schema_dir), 'temp')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    # change for testing!!!
    local_date = '2020-05-20T21:44:31.08'
    #local_date = '2020-05-24T19:06:17.408700'

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

    if submit == 'yes' and user_auth is True:

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

                # Get the name of the species from the provided id
                # or vice-versa
                species_info = aux.species_ids(species_id, base_url, headers_get)
                species_id, species_name = species_info
                print('\nNS species with identifier {0} '
                      'is {1}.'.format(species_id, species_name))

                # compare list of genes, if they do not intersect, halt process
                # get list of loci for schema in the NS
                ns_loci_get = aux.simple_get_request(base_url, headers_get,
                                                     ['species', species_id,
                                                      'schemas', schema_id,
                                                      'loci'])
                # get loci files names from response
                ns_schema_loci = []
                ns_schema_locid_map = {}
                for l in ns_loci_get.json()['Loci']:
                    locus_file = l['original_name']['value']
                    locus_name = l['name']['value'] + '.fasta'
                    ns_schema_loci.append(locus_name)
                    locus_uri = l['locus']['value']
                    ns_schema_locid_map[locus_name] = (locus_uri, locus_file)
        
                completed = list(pickled_loci.keys())
                incomplete = []
                for locus in not_in_ns:
                    if locus in ns_schema_loci:
                        if len(not_in_ns[locus]) > 0:
                            incomplete.append(locus)
                        else:
                            if locus not in completed:
                                completed.append(locus)
        
                if len(incomplete) > 0:
        
                    uniq_local = {k: v for k, v in not_in_ns.items()
                                  if k in incomplete}
                    loci_uris = {k: v for k, v in ns_schema_locid_map.items()
                                 if k in incomplete}
        
                    # determine sequences that are not in the NS
                    # we synced before so we just need to go to each file and
                    # identify the alleles with '*' as the alleles to be added!
                    # we have to submit them and alter locally too...
                    upload, comp = determine_upload(uniq_local, incomplete,
                                                    loci_uris, schema_dir,
                                                    temp_dir, base_url, headers_get)
        
                    completed.extend(comp)
                    completed = list(set(completed))

                    # create files with length values to send!
                    # create an endpooint to send new alleles, it receives files
                    # and calls a script that inserts the new alleles and updates
                    # the pre-computed files!

                    # add sequences to the NS
                    # use multiprocessing
                    alleles_data = []
                    for locus, info in upload.items():
                        new_loci_url = info[0]
                        name = info[1]
                        label = info[2]
                        url = info[3]
                        allele_seq_list = [rec[1] for rec in info[4:]]
                        start_id = min([int(rec[0].split('_')[-1])
                                        for rec in info[4:]])
                        post_data = aux.create_allele_data(allele_seq_list,
                                                           new_loci_url,
                                                           name,
                                                           label,
                                                           url,
                                                           species_name,
                                                           True,
                                                           headers_post,
                                                           user_id,
                                                           start_id)
                        alleles_data.append(post_data)
        
                    for inlist in alleles_data:
                        post_results = []
                        total_inserted = 0
                        total_alleles = len(inlist)
        
                        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                            # Start the load operations and mark each future
                            # with its URL
                            for res in executor.map(aux.post_allele, inlist):
                                post_results.append(res)
                                total_inserted += 1
                                print('\r',
                                      'Processed {0}/{1} '
                                      'alleles.'.format(total_inserted,
                                                        total_alleles),
                                      end='')
#        else:
#            print('There are no new alleles in the local schema that '
#                  'need to be sent to the NS.')

    elif submit == 'no':
        completed = list(not_in_ns.keys())

    # change profiles in the master file so that local
    # profiles have updated allele identifiers

    # change pickled files to FASTA files
    for locus in completed:
        locus_id = locus.rstrip('.fasta')
        pickled_file = pickled_loci[locus]
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
    PrepExternalSchema.main(temp_dir, schema_dir,
                            core_num, bsr, 0, 11)

    # delete invalid alleles and genes files
    parent_dir = os.path.dirname(schema_dir)
    files = [os.path.join(parent_dir, file)
             for file in os.listdir(parent_dir)
             if 'invalid' in file]

    for f in files:
        os.remove(f)

    # delete temp directory
    shutil.rmtree(temp_dir)

    # update NS config file with latest server time
    ns_configs = os.path.join(schema_dir, '.ns_config')
    with open(ns_configs, 'wb') as nc:
        pickle.dump([server_time, schema_uri], nc)

    if user_role is True and submit == 'yes':
        print('Updated local schema with {0} alleles '
              'and sent {1} novel alleles to the '
              'NS.'.format(new_count, len(alleles_data)))
    else:
        print('Updated local schema with {0} alleles.'.format(new_count))

    # unlock schema if it was locked
    if locked is True:
        schema_unlock = aux.simple_post_request(base_url, headers_post,
                                                ['species', species_id,
                                                 'schemas', schema_id, 'lock'],
                                                {'action': 'unlock'})

    print('Done!')


if __name__ == "__main__":

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3])
