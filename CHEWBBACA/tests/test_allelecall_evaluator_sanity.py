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
                "AlleleCallEvaluator",
                "-i",
                "data/allelecall_data/fake_results",
                "-g",
                "data/allelecall_data/sagalactiae_schema",
                "-o",
                "results_report",
            ],
            "Input argument is not a valid directory. Exiting...",
        ),
        # (
        #     [
        #         "chewBBACA.py",
        #         "AlleleCallEvaluator",
        #         "-g",
        #         "this/path/aint/real",
        #         "-o",
        #         "schema_report",
        #     ],
        #     "Input argument is not a valid directory. Exiting...",
        # ),
    ],
)
def test_allelecallevaluator_invalid_input(test_args, expected):

    with pytest.raises(SystemExit) as e:
        with patch.object(sys, "argv", test_args):
            chewBBACA.main()

    # Delete output directory
    try:
        shutil.rmtree(test_args[7])
    except Exception as e2:
        pass

    assert e.type == SystemExit
    assert expected in e.value.code


@pytest.mark.parametrize(
    "test_args, expected",
    [
        (
            [
                "chewBBACA.py",
                "AlleleCallEvaluator",
                "-i",
                "data/allelecall_data/test_results",
                "-g",
                "data/allelecall_data/sagalactiae_schema",
                "-o",
                "results_report",
            ],
            0,
        ),
        # (
        #     [
        #         "chewBBACA.py",
        #         "AlleleCallEvaluator",
        #         "-i",
        #         "data/",
        #         "-g",
        #         "data/",
        #         "-o",
        #         "results_report",
        #     ],
        #     0,
        # ),
    ],
)
def test_allelecallevaluator_valid(test_args, expected):
    with pytest.raises(SystemExit) as e:
        chewBBACA.main()

    try:
        shutil.rmtree(test_args[7])
    except Exception as e2:
        pass

    # Check exit code
    assert e.type == SystemExit
    assert expected == e.value.code
