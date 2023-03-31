#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import py
import os
import sys
import pickle
import pytest
import shutil
import filecmp
from unittest.mock import patch

# from contextlib import nullcontext as does_not_raise

from CHEWBBACA import chewBBACA


@pytest.mark.parametrize(
    "test_args, expected",
    [
        (
            [
                "chewBBACA.py",
                "SchemaEvaluator",
                "-g",
                "data/schemaevaluator_data/empty_files",
                "-o",
                "schema_report",
            ],
            "Could not get input files.",
        ),
        (
            [
                "chewBBACA.py",
                "SchemaEvaluator",
                "-g",
                "data/schemaevaluator_data/zero_bytes_pair",
                "-o",
                "schema_report",
            ],
            "Could not get input files.",
        ),
        (
            [
                "chewBBACA.py",
                "SchemaEvaluator",
                "-g",
                "this/path/aint/real",
                "-o",
                "schema_report",
            ],
            "Input argument is not a valid directory. Exiting...",
        ),
    ],
)
def test_schemaEvaluator_invalid_input(test_args, expected):

    with pytest.raises(SystemExit) as e:
        with patch.object(sys, "argv", test_args):
            chewBBACA.main()

    try:
        shutil.rmtree(test_args[5])
    except Exception as e2:
        pass

    assert e.type == SystemExit
    assert expected in e.value.code


# @pytest.mark.parametrize(
#     "test_args, expected",
#     [
#         (
#             [
#                 "chewBBACA.py",
#                 "SchemaEvaluator",
#                 "-g",
#                 "data/schemaevaluator_data/test_schema",
#                 "-o",
#                 "schema_report",
#                 "--loci-reports",
#                 "--add-sequences",
#             ],
#             "data/schemaevaluator_data/expected_results",
#         )
#     ],
# )
# def test_schemaEvaluator_valid(test_args, expected):
#     with patch.object(sys, "argv", test_args):
#         capture = py.io.StdCapture()
#         chewBBACA.main()
#         stdout, stderr = capture.reset()

#     # check Schema Report HTML file
#     schema_report_html = os.path.join(test_args[5], "schema_report.html")

#     expected_schema_file = os.path.join(expected, "schema_report.html")

#     schema_report_cmp = filecmp.cmp(schema_report_html, expected_schema_file, shallow=True)
#     assert(schema_report_cmp) is True

#     # check Locus Report HTML files and JS bundle
#     locus_report_files = [
#         os.path.join(test_args[5], "loci_reports", file)
#         for file in os.listdir(os.path.join(test_args[5], "loci_reports"))
#     ]
#     locus_report_files.sort()

#     expected_report_files = [
#         os.path.join(expected, "loci_reports", file)
#         for file in os.listdir(os.path.join(expected, "loci_reports"))
#     ]
#     expected_report_files.sort()

#     # assert that files in each pair are equal
#     file_cmps = []
#     for i, file in enumerate(expected_report_files):
#         file_cmps.append(filecmp.cmp(file, locus_report_files[i], shallow=False))

#     assert all(file_cmps) is True
