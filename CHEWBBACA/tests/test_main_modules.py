#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import py
import os
import sys
import pytest
import filecmp
from unittest.mock import patch
#from contextlib import nullcontext as does_not_raise

from CHEWBBACA import chewBBACA


@pytest.mark.parametrize(
        'test_args, expected',
        [(['chewBBACA.py', '-h'], 0),
         (['chewBBACA.py', '--help'], 0),
         (['chewBBACA.py', '-v'], 0),
         (['chewBBACA.py', '--version'], 0),
         (['chewBBACA.py', 'CreateSchema', '-h'], 0),
         (['chewBBACA.py', 'CreateSchema', '--help'], 0),
         (['chewBBACA.py', 'AlleleCall', '-h'], 0),
         (['chewBBACA.py', 'AlleleCall', '--help'], 0),
         (['chewBBACA.py', 'TestGenomeQuality', '-h'], 0),
         (['chewBBACA.py', 'TestGenomeQuality', '--help'], 0),
         (['chewBBACA.py', 'ExtractCgMLST', '-h'], 0),
         (['chewBBACA.py', 'ExtractCgMLST', '--help'], 0),
         (['chewBBACA.py', 'SchemaEvaluator', '-h'], 0),
         (['chewBBACA.py', 'SchemaEvaluator', '--help'], 0),
         (['chewBBACA.py', 'DownloadSchema', '-h'], 0),
         (['chewBBACA.py', 'DownloadSchema', '--help'], 0),
         (['chewBBACA.py', 'LoadSchema', '-h'], 0),
         (['chewBBACA.py', 'LoadSchema', '--help'], 0),
         (['chewBBACA.py', 'SyncSchema', '-h'], 0),
         (['chewBBACA.py', 'SyncSchema', '--help'], 0),
         (['chewBBACA.py', 'NSStats', '-h'], 0),
         (['chewBBACA.py', 'NSStats', '--help'], 0),
         ])
def test_modules_exit_codes(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert e.value.code == expected


@pytest.mark.parametrize(
        'test_args, expected',
        [(['chewBBACA.py', 'AlleleCall',
           '-i', 'data/allelecall_data/test_genome',
           '-g', 'data/allelecall_data/sagalactiae_schema',
           '-o', 'allelecall_results'],
         'data/allelecall_data/test_results')
         ])
def test_allelecall_valid(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        capture = py.io.StdCapture()
        chewBBACA.main()
        stdout, stderr = capture.reset()

    # check text printed to stdout
    assert 'Writing output files' in stdout

    # check output files
    for root, dirs, files in os.walk(test_args[7]):
        output_files = [os.path.join(root, file)
                        for file in files
                        if 'logging_info.txt' != file]

    expected_files = [os.path.join(expected, file)
                      for file in os.listdir(expected)
                      if 'logging_info.txt' != file]

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
