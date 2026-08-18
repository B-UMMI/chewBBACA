#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the MergeResults module works as expected.
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
	 (ta.MERGERESULTS_TEST_DEFAULT, 'data/mergeresults_data/expected_results_default'),
     (ta.MERGERESULTS_TEST_COMMON, 'data/mergeresults_data/expected_results_common'),
    ],
	indirect=True
)
def test_mergeresults_valid_input(monkeypatch, args_fixture):
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


@pytest.mark.parametrize(
    "args_fixture",
    [
	 (ta.MERGERESULTS_TEST_INPUT_FILE, ct.MERGERESULTS_INPUTFILE_EXCEPTION[:33]),
	 (ta.MERGERESULTS_TEST_NOCOMMONLOCI, ct.MERGERESULTS_NOCOMMONLOCI_EXCEPTION),
	 (ta.MERGERESULTS_TEST_LOCIDIFFER, ct.MERGERESULTS_LOCIDIFFER_EXCEPTION),
    ],
	indirect=True
)
def test_mergeresults_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			# Check that the exit message includes expected message
			assert args_fixture[1] in e.value.code
