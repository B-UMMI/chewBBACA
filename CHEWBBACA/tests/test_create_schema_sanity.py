#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    [(['chewBBACA.py', 'CreateSchema',
       '-i', 'data/createschema_data/mock_schema_dir',
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

    expected_files = [os.path.join(expected, file)
                      for file in os.listdir(expected)
                      if 'short' != file]

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
