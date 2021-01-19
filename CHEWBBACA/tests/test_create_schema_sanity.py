#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import py
import os
import sys
import pickle
import pytest
import filecmp
from unittest.mock import patch
#from contextlib import nullcontext as does_not_raise

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
      'data/createschema_data/expected_results')
     ])
def test_createschema_valid(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        capture = py.io.StdCapture()
        chewBBACA.main()
        stdout, stderr = capture.reset()

    # check text printed to stdout

    # check if Prodigal created files
    assert 'All files were created.' in stdout

    # check output files
    output_files = [os.path.join(test_args[5], file)
                    for file in os.listdir(test_args[5])
                    if 'short' != file]
    output_files.sort()

    expected_files = [os.path.join(expected, file)
                      for file in os.listdir(expected)
                      if 'short' != file]
    expected_files.sort()

    # get config files
    # genes_lists = [output_files.pop(0), expected_files.pop(0)]
    # schemas_configs = [output_files.pop(0), expected_files.pop(0)]

    # compare configs
    # assert pickle_loader(genes_lists[0]).sort() == pickle_loader(genes_lists[1]).sort()
    # assert pickle_loader(schemas_configs[0]) == pickle_loader(schemas_configs[1])

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


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/invalid_genome_dir',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'], 1)
     ])
def test_createschema_invalid_pairs(test_args, expected):
    with pytest.raises(ValueError) as e:
        with patch.object(sys, 'argv', test_args):
            capture = py.io.StdCapture()
            chewBBACA.main()
            stdout, stderr = capture.reset()

        assert 'BLAST Database error: No alias or index file found for protein database' in stdout

        assert e.type == ValueError


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/genome_dir_with_empty_genomes',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'], "Could not get input files."),
     (['chewBBACA.py', 'CreateSchema',
         '-i', 'data/createschema_data/zero_bytes_pair',
         '-o', 'createschema_results',
         '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'], "Could not get input files."),
     (['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_schema_dir',
       '-o', 'createschema_results',
       '--ptf', 'path/does/not/exist'], "Cannot find specified Prodigal training file.")
     ])
def test_createschema_empty_pairs(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert expected in e.value.code


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/genome_dir_with_header_fasta_only',
       '-o', 'createschema_results',
       '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn'], 1)
     ])
def test_createschema_header_only_pairs(test_args, expected):
    with pytest.raises(ValueError) as e:
        with patch.object(sys, 'argv', test_args):
            capture = py.io.StdCapture()
            chewBBACA.main()
            stdout, stderr = capture.reset()

        assert 'Error:  no input sequences to analyze.' in stdout

        assert 'BLAST Database error: No alias or index file found for protein database' in stdout

        assert e.type == ValueError
