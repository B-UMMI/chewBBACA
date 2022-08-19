#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import py
import os
import sys
import shutil
import pickle
import pytest
import filecmp
from unittest.mock import patch
# from contextlib import nullcontext as does_not_raise

from CHEWBBACA import chewBBACA


def pickle_loader(pickle_in):
    """
    """

    with open(pickle_in, 'rb') as pi:
        data = pickle.load(pi)

    return data


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_genome_dir',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      'data/createschema_data/expected_results'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_genome_list/mock_genomes.txt',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      'data/createschema_data/expected_results'),
     ])
def test_createschema_valid(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        capture = py.io.StdCapture()
        chewBBACA.main()
        stdout, stderr = capture.reset()

    # check output files
    schema_seed = os.path.join(test_args[5], 'schema_seed')
    output_files = [os.path.join(schema_seed, file)
                    for file in os.listdir(schema_seed)
                    if 'short' != file]
    output_files.sort()

    expected_seed = os.path.join(expected, 'schema_seed')
    expected_files = [os.path.join(expected_seed, file)
                      for file in os.listdir(expected_seed)
                      if 'short' != file]
    expected_files.sort()

    # get config files
    genes_lists = [output_files.pop(0), expected_files.pop(0)]
    schemas_configs = [output_files.pop(0), expected_files.pop(0)]

    # compare configs
    assert pickle_loader(genes_lists[0]).sort() == pickle_loader(genes_lists[1]).sort()
    assert pickle_loader(schemas_configs[0]) == pickle_loader(schemas_configs[1])

    files = output_files + expected_files
    basename_dict = {}
    for f in files:
        basename = os.path.basename(f)
        basename_dict.setdefault(basename, []).append(f)

    # assert that files in each pair are equal
    file_cmps = []
    for k, v in basename_dict.items():
        file_cmps.append(filecmp.cmp(v[0], v[1], shallow=False))

    assert all(file_cmps) is True

    # delete results
    try:
        shutil.rmtree(test_args[5])
    except Exception as e2:
        pass


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/genome_dir_with_empty_genomes',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      '\nCould not get input files. Please provide a directory'
      ' with FASTA files or a file with the list of full '
      'paths to the FASTA files and ensure that filenames end '
      'with one of the following suffixes: '
      '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\'].'),
     (['chewBBACA.py', 'CreateSchema',
         '-i', 'data/createschema_data/zero_bytes_pair',
         '-o', 'createschema_results',
         '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      '\nCould not get input files. Please provide a directory'
      ' with FASTA files or a file with the list of full '
      'paths to the FASTA files and ensure that filenames end '
      'with one of the following suffixes: '
      '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\'].'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_schema_dir',
       '-o', 'createschema_results',
       '--ptf', 'path/does/not/exist'],
      'Invalid path for Prodigal training file.'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/genome_dir_with_header_fasta_only',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      '\nCould not predict gene sequences from any '
      'of the input files.\nPlease provide input files '
      'in the accepted FASTA format.'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/invalid_genome_dir',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      '\nCould not predict gene sequences from any '
      'of the input files.\nPlease provide input files '
      'in the accepted FASTA format.')
     ])
def test_createschema_empty_pairs(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert expected in e.value.code
