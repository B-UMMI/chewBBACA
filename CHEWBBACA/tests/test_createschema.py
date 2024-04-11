#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the CreateSchema module works as expected.
"""


import os
import sys
import shutil
import pickle
import pytest
import filecmp
from unittest.mock import patch
# from contextlib import nullcontext as does_not_raise

from CHEWBBACA import chewBBACA


def pickle_loader(input_file):
    """Use the Pickle module to de-serialize an object.

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
    with open(input_file, 'rb') as pinfile:
        content = pickle.load(pinfile)

    return content


# Test successful process
# 1) With a path to a folder that contains FASTA files
# 2) With a path to a file that contains a list of paths to FASTA files
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
	 (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_cds_dir',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn',
	   '--cds-input'],
      'data/createschema_data/expected_results'),
     ])
def test_createschema_valid(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        chewBBACA.main()

    # Check schema files
	# Schema created by test
    schema_seed = os.path.join(test_args[5], 'schema_seed')
    output_files = [os.path.join(schema_seed, file)
                    for file in os.listdir(schema_seed)
                    if 'short' not in file]
    output_files.sort()
	# Expected schema data
    expected_seed = os.path.join(expected, 'schema_seed')
    expected_files = [os.path.join(expected_seed, file)
                      for file in os.listdir(expected_seed)
                      if 'short' not in file]
    expected_files.sort()

    # Get config values
    genes_lists = [output_files.pop(0), expected_files.pop(0)]
    schemas_configs = [output_files.pop(0), expected_files.pop(0)]

    # Compare lists of loci
    assert pickle_loader(genes_lists[0]).sort() == pickle_loader(genes_lists[1]).sort()
    # Read config values
    # Ignore chewBBACA version value
    configs1 = pickle_loader(schemas_configs[0])
    del configs1['chewBBACA_version']
    configs2 = pickle_loader(schemas_configs[1])
    del configs2['chewBBACA_version']
    assert configs1 == configs2

	# Group test schema and expected schema files based on basename
    files = output_files + expected_files
    basename_dict = {}
    for f in files:
        basename = os.path.basename(f)
        basename_dict.setdefault(basename, []).append(f)

    # Assert that files in each pair are equal
    file_cmps = []
    for k, v in basename_dict.items():
        file_cmps.append(filecmp.cmp(v[0], v[1], shallow=False))

    assert all(file_cmps) is True

    # Delete output folder or next test might fail
    try:
        shutil.rmtree(test_args[5])
    except Exception as exc_msg:
        print(exc_msg)


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/genome_dir_with_empty_genomes',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      'Could not get input files. Please provide a directory'
      ' with FASTA files or a file with the list of full '
      'paths to the FASTA files and ensure that filenames end '
      'with one of the following extensions: '
      '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
     (['chewBBACA.py', 'CreateSchema',
         '-i', 'data/createschema_data/zero_bytes_pair',
         '-o', 'createschema_results',
         '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      'Could not get input files. Please provide a directory'
      ' with FASTA files or a file with the list of full '
      'paths to the FASTA files and ensure that filenames end '
      'with one of the following extensions: '
      '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_schema_dir',
       '-o', 'createschema_results',
       '--ptf', 'path/does/not/exist'],
      'Invalid path for Prodigal training file.'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/genome_dir_with_header_fasta_only',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      'Could not predict CDSs from any of the input files.'
	  '\nPlease provide input files in the accepted FASTA format.'),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/invalid_genome_dir',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'],
      'Could not predict CDSs from any of the input files.'
	  '\nPlease provide input files in the accepted FASTA format.')
     ])
def test_createschema_empty_pairs(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
	# Check that the exit message includes expected message
    assert expected in e.value.code

    # Delete output folder or next test might fail
    try:
        shutil.rmtree(test_args[5])
    except Exception as exc_msg:
        print(exc_msg)
