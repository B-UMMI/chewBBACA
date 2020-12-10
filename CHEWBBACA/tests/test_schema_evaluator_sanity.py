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


@pytest.mark.parametrize(
    'test_args, expected',
    [(['chewBBACA.py', 'SchemaEvaluator',
       '-i', 'data/schemaEvaluatorData/test_yersinia',
       '-l', './schemaEvaluatorData_results', ],
      'data/schemaEvaluatorData/expected_results')
     ])
def test_schemaEvaluator_valid(test_args, expected):
    with patch.object(sys, 'argv', test_args):
        capture = py.io.StdCapture()
        chewBBACA.main()
        stdout, stderr = capture.reset()

    # check text printed to stdout

    # check if Prodigal created files
    assert 'The report has been created.' in stdout

    # check output html files
    output_html_files = [os.path.join(test_args[5], "html_files", file)
                         for file in os.listdir(os.path.join(test_args[5], "html_files"))
                         ]
    output_html_files.sort()

    expected_html_files = [os.path.join(expected, "html_files", file)
                           for file in os.listdir(os.path.join(expected, "html_files"))
                           ]
    expected_html_files.sort()

    html_files = output_html_files + expected_html_files
    basename_html_dict = {}
    for f in html_files:
        basename = os.path.basename(f)
        basename_html_dict.setdefault(basename, []).append(f)

    # assert that files in each pair are equal
    file_cmps = []
    for k, v in basename_html_dict.items():
        file_cmps.append(filecmp.cmp(v[0], v[1], shallow=False))

    assert all(file_cmps) is True

    # check output main files
    output_main_files = [os.path.join(test_args[5], "SchemaEvaluator_pre_computed_data", file)
                         for file in os.listdir(os.path.join(test_args[5], "SchemaEvaluator_pre_computed_data"))
                         if "prot_files" != file]
    output_main_files.sort()

    expected_main_files = [os.path.join(expected, "SchemaEvaluator_pre_computed_data", file)
                           for file in os.listdir(os.path.join(expected, "SchemaEvaluator_pre_computed_data"))
                           if "prot_files" != file]
    expected_main_files.sort()

    main_files = output_main_files + expected_main_files
    basename_main_dict = {}
    for f in main_files:
        basename = os.path.basename(f)
        basename_main_dict.setdefault(basename, []).append(f)

    # assert that files in each pair are equal
    file_cmps_main = []
    for k, v in basename_main_dict.items():
        file_cmps_main.append(filecmp.cmp(v[0], v[1], shallow=False))

    assert all(file_cmps_main) is True

    # check output protein files
    output_prot_files = [os.path.join(test_args[5], "SchemaEvaluator_pre_computed_data", "prot_files", file)
                         for file in os.listdir(os.path.join(test_args[5], "SchemaEvaluator_pre_computed_data", "prot_files"))
                         if "exceptions" != file]
    output_prot_files.sort()

    expected_prot_files = [os.path.join(expected, "SchemaEvaluator_pre_computed_data", "prot_files", file)
                           for file in os.listdir(os.path.join(expected, "SchemaEvaluator_pre_computed_data", "prot_files"))
                           if "exceptions" != file]
    expected_prot_files.sort()

    prot_files = output_prot_files + expected_prot_files
    basename_prot_dict = {}
    for f in prot_files:
        basename = os.path.basename(f)
        basename_prot_dict.setdefault(basename, []).append(f)

    # assert that files in each pair are equal
    file_cmps_prot = []
    for k, v in basename_prot_dict.items():
        file_cmps_prot.append(filecmp.cmp(v[0], v[1], shallow=False))

    assert all(file_cmps_prot) is True

    # check output exception files
    output_exc_files = [os.path.join(test_args[5], "SchemaEvaluator_pre_computed_data", "prot_files", "exceptions", file)
                         for file in os.listdir(os.path.join(test_args[5], "SchemaEvaluator_pre_computed_data", "prot_files", "exceptions"))
                        ]
    output_exc_files.sort()

    expected_exc_files = [os.path.join(expected, "SchemaEvaluator_pre_computed_data", "prot_files", "exceptions", file)
                           for file in os.listdir(os.path.join(expected, "SchemaEvaluator_pre_computed_data", "prot_files", "exceptions"))
                          ]
    expected_exc_files.sort()

    exc_files = output_exc_files + expected_exc_files
    basename_exc_dict = {}
    for f in exc_files:
        basename = os.path.basename(f)
        basename_exc_dict.setdefault(basename, []).append(f)

    # assert that files in each pair are equal
    file_cmps_exc = []
    for k, v in basename_exc_dict.items():
        file_cmps_exc.append(filecmp.cmp(v[0], v[1], shallow=False))

    assert all(file_cmps_exc) is True
