#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the AlleleCall module works as expected.
"""


import py
import os
import sys
import shutil
import pytest
import filecmp
from unittest.mock import patch
# from contextlib import nullcontext as does_not_raise

from CHEWBBACA import chewBBACA


args_template = ['chewBBACA.py', 'AlleleCall',
                 '-i', 'data/allelecall_data/test_genome',
                 '-g', 'data/allelecall_data/sagalactiae_schema',
                 '-o', 'allelecall_results']


@pytest.mark.parametrize(
        'test_args, expected',
        [(args_template,
         'data/allelecall_data/test_results/mode4'),
		(args_template+['--mode', '1'],
          'data/allelecall_data/test_results/mode1'),
		(args_template+['--mode', '2'],
          'data/allelecall_data/test_results/mode2'),
		(args_template+['--mode', '3'],
          'data/allelecall_data/test_results/mode3'),
		(args_template[0:3]+['data/allelecall_data/test_cds_input']+args_template[4:]+['--cds-input'],
          'data/allelecall_data/test_results/mode4'),
        (args_template[0:3]+['data/allelecall_data/test_genomes_list/test_genomes.txt']+args_template[4:],
          'data/allelecall_data/test_results/mode4'),
        (args_template+['--gl', 'data/allelecall_data/test_genes_list/test_genes_extension.txt'],
         'data/allelecall_data/test_genes_list/test_genes_results'),
        (args_template+['--gl', 'data/allelecall_data/test_genes_list/test_genes_no_extension.txt'],
         'data/allelecall_data/test_genes_list/test_genes_results'),
        (args_template+['--gl', 'data/allelecall_data/test_genes_list/test_genes_path.txt'],
         'data/allelecall_data/test_genes_list/test_genes_results'),
        ])
def test_allelecall_valid(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        capture = py.io.StdCapture()
        chewBBACA.main()
        stdout, stderr = capture.reset()

    # Get paths to output files
    for root, dirs, files in os.walk(test_args[7]):
        output_files = [os.path.join(root, file)
                        for file in files
						if 'logging_info.txt' != file]

	# Get paths to files with expected results
    expected_files = [os.path.join(expected, file)
                      for file in os.listdir(expected)]

	# Group test results and expected results based on basename
    files = output_files + expected_files
    basename_dict = {}
    for f in files:
        basename = os.path.basename(f)
        basename_dict.setdefault(basename, []).append(f)

    # Assert that files in each pair are equal
    file_cmps = []
    for k, v in basename_dict.items():
		assert filecmp.cmp(v[0], v[1], shallow=False) is True
        #file_cmps.append(filecmp.cmp(v[0], v[1], shallow=False))

    #assert all(file_cmps) is True

	# Delete output folder or next test might fail
    try:
        shutil.rmtree(test_args[7])
    except Exception as exc_msg:
        print(exc_msg)


@pytest.mark.parametrize(
        'test_args, expected',
        [(['chewBBACA.py', 'AlleleCall',
           '-i', 'data/prep_data/empty_dir',
           '-g', 'data/allelecall_data/sagalactiae_schema',
           '-o', 'allelecall_results'],
         'Could not get input files. Please provide a directory'
         ' with FASTA files or a file with the list of full '
         'paths to the FASTA files and ensure that filenames end '
         'with one of the following extensions: '
         '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
         (['chewBBACA.py', 'AlleleCall',
           '-i', 'data/createschema_data/genome_dir_with_empty_genomes',
           '-g', 'data/allelecall_data/sagalactiae_schema',
           '-o', 'allelecall_results'],
         'Could not get input files. Please provide a directory'
         ' with FASTA files or a file with the list of full '
         'paths to the FASTA files and ensure that filenames end '
         'with one of the following extensions: '
         '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
         (['chewBBACA.py', 'AlleleCall',
           '-i', 'data/createschema_data/zero_bytes_pair',
           '-g', 'data/allelecall_data/sagalactiae_schema',
           '-o', 'allelecall_results'],
         'Could not get input files. Please provide a directory'
         ' with FASTA files or a file with the list of full '
         'paths to the FASTA files and ensure that filenames end '
         'with one of the following extensions: '
         '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
         (['chewBBACA.py', 'AlleleCall',
           '-i', 'this/path/aint/real',
           '-g', 'data/allelecall_data/sagalactiae_schema',
           '-o', 'allelecall_results'],
          'Input argument is not a valid directory or '
          'file with a list of paths to FASTA files. Please provide a '
          'valid input, either a folder with FASTA files '
          'or a file with the list of full paths to FASTA '
          'files (one per line and ending with one of the '
          'following file extensions: [\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\']).')
        ])
def test_invalid_input(test_args, expected):

    # create empty dir for empty dir test
    if 'empty_dir' in test_args[3] and os.path.isdir(test_args[3]) is False:
        os.mkdir(test_args[3])

    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert expected in e.value.code


@pytest.mark.parametrize(
        'test_args, expected',
        [(args_template+['--bsr', '-1'],
         '\nBSR value is not contained in the [0.0, 1.0] interval.'),
         (args_template+['--bsr', '1.1'],
         '\nBSR value is not contained in the [0.0, 1.0] interval.'),
         (args_template+['--bsr', 'sus'],
          '\nInvalid BSR value of sus. BSR value must be '
          'contained in the [0.0, 1.0] interval.'),
         (args_template+['--l', '-1'],
          '\nInvalid minimum sequence length value. '
          'Must be equal or greater than 0.'),
         (args_template+['--l', 'sus'],
          '\nInvalid minimum sequence length value. '
          'Value must be a positive integer.'),
         (args_template+['--st', '-1'],
          '\nInvalid size threshold value. '
          'Must be contained in the [0.0, 1.0] interval.'),
         (args_template+['--st', '1.1'],
          '\nInvalid size threshold value. '
          'Must be contained in the [0.0, 1.0] interval.'),
         (args_template+['--st', 'sus'],
          '\nInvalid size threshold value used to '
          'create schema. Value must be None or a '
          'positive float in the [0.0, 1.0] interval.'),
         ])
def test_invalid_args(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert e.value.code == expected
