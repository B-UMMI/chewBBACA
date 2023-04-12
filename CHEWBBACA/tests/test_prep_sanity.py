#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import py
import sys
import pickle
import pytest
import shutil
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
        [(['chewBBACA.py', 'PrepExternalSchema',
           '-i', 'data/prep_data/empty_dir',
           '-o', 'adapted_schema'],
         'Could not get input files. Please provide a directory'
         ' with FASTA files or a file with the list of full '
         'paths to the FASTA files and ensure that filenames end '
         'with one of the following extensions: '
         '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
         (['chewBBACA.py', 'PrepExternalSchema',
           '-i', 'data/prep_data/empty_files',
           '-o', 'adapted_schema'],
         'Could not get input files. Please provide a directory'
         ' with FASTA files or a file with the list of full '
         'paths to the FASTA files and ensure that filenames end '
         'with one of the following extensions: '
         '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
         (['chewBBACA.py', 'PrepExternalSchema',
           '-i', 'data/prep_data/zero_bytes_pair',
           '-o', 'adapted_schema'],
         'Could not get input files. Please provide a directory'
         ' with FASTA files or a file with the list of full '
         'paths to the FASTA files and ensure that filenames end '
         'with one of the following extensions: '
         '[\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\'].'),
         (['chewBBACA.py', 'PrepExternalSchema',
           '-i', 'this/path/aint/real',
           '-o', 'adapted_schema'],
          'Input argument is not a valid directory or '
          'file with a list of paths to FASTA files. Please provide a '
          'valid input, either a folder with FASTA files '
          'or a file with the list of full paths to FASTA '
          'files (one per line and ending with one of the '
          'following file extensions: [\'.fasta\', \'.fna\', \'.ffn\', \'.fa\', \'.fas\']).')
         ])
def test_prep_invalid_input(test_args, expected):

    # create empty dir for empty dir test
    if 'empty_dir' in test_args[3] and os.path.isdir(test_args[3]) is False:
        os.mkdir(test_args[3])

    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    try:
        shutil.rmtree(test_args[5])
    except Exception as e2:
        pass

    assert e.type == SystemExit
    assert expected in e.value.code


@pytest.mark.parametrize(
        'test_args, expected',
        [
            (['chewBBACA.py', 'PrepExternalSchema',
              '-i', 'data/prep_data/valid_input',
              '-o', 'preped_schema'],
             'data/prep_data/expected_results_valid'),
            (['chewBBACA.py', 'PrepExternalSchema',
              '-i', 'data/prep_data/file_extensions',
              '-o', 'preped_schema'],
             'data/prep_data/expected_results_extension'),
        ]
)
def test_prep_valid_input(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        capture = py.io.StdCapture()
        chewBBACA.main()
        stdout, stderr = capture.reset()

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
    genes_lists = [output_files.pop(0), expected_files.pop(0)]
    schemas_configs = [output_files.pop(0), expected_files.pop(0)]

    # compare configs
    assert pickle_loader(genes_lists[0]).sort() == pickle_loader(genes_lists[1]).sort()
    assert pickle_loader(schemas_configs[0]) == pickle_loader(schemas_configs[1])

    # compare FASTA files
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

    # Delete output directory before next test
    try:
        shutil.rmtree(test_args[5])
    except Exception as e2:
        pass
