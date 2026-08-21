#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the ComputeDistances module works as expected.
"""


import os
import sys
import shutil
import pytest
import filecmp

from CHEWBBACA import chewBBACA
from CHEWBBACA.utils import constants as ct
from CHEWBBACA.tests import test_arguments as ta
from CHEWBBACA.utils import file_operations as fo


# Use the tmp_path fixture to create a tmp directory for each test
@pytest.fixture
def args_fixture(request, tmp_path):
	# Setup step
	args = request.param
	args[0][5] = os.path.join(tmp_path, args[0][5])
	yield args
	# Teardown step only runs after yield
	shutil.rmtree(tmp_path)


# Test with valid inputs
@pytest.mark.parametrize(
    "args_fixture",
    [
	 (ta.COMPUTEDISTANCES_TEST_DEFAULT, 'data/computedistances_data/expected_results/default'),
     (ta.COMPUTEDISTANCES_TEST_DEFAULT_SYMMETRIC, 'data/computedistances_data/expected_results/default_symmetric'),
	 (ta.COMPUTEDISTANCES_TEST_DEFAULT_SIMILARITY, 'data/computedistances_data/expected_results/default_similarity'),
	 (ta.COMPUTEDISTANCES_TEST_DEFAULT_SYMMETRIC_SIMILARITY, 'data/computedistances_data/expected_results/default_symmetric_similarity'),
	 (ta.COMPUTEDISTANCES_TEST_DEFAULT_TABLE, 'data/computedistances_data/expected_results/default_table'),
	 (ta.COMPUTEDISTANCES_TEST_JACCARD_SYMMETRIC, 'data/computedistances_data/expected_results/jaccard_symmetric'),
	 (ta.COMPUTEDISTANCES_TEST_JACCARD_SYMMETRIC_SIMILARITY, 'data/computedistances_data/expected_results/jaccard_symmetric_similarity'),
	 (ta.COMPUTEDISTANCES_TEST_LOCI_SYMMETRIC, 'data/computedistances_data/expected_results/loci_symmetric'),
	 (ta.COMPUTEDISTANCES_TEST_LOCI_SYMMETRIC_SIMILARITY, 'data/computedistances_data/expected_results/loci_symmetric_similarity'),
	 (ta.COMPUTEDISTANCES_TEST_CORE_SYMMETRIC, 'data/computedistances_data/expected_results/core_symmetric')
	 (ta.COMPUTEDISTANCES_TEST_CORE_SYMMETRIC_SIMILARITY, 'data/computedistances_data/expected_results/core_symmetric_similarity')
    ],
	indirect=True
)
def test_computedistances_valid_input(monkeypatch, args_fixture):
    # Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		# List output files
		output_dir = args_fixture[0][5]
		output_files = [os.path.join(output_dir, file)
						for file in os.listdir(output_dir)]
		output_files.sort()
		# List expected results files
		expected_dir = args_fixture[1]
		expected_files = [os.path.join(expected_dir, file)
						  for file in os.listdir(expected_dir)]
		expected_files.sort()

		# Group test results and expected results files based on basename
		files = output_files + expected_files
		basename_dict = {}
		for f in files:
			basename = os.path.basename(f)
			basename_dict.setdefault(basename, []).append(f)

		# Assert that files in each pair are equal
		for k, v in basename_dict.items():
			assert filecmp.cmp(v[0], v[1], shallow=False) is True
