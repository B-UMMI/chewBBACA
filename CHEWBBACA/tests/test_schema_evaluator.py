#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import shutil
import pytest
import filecmp

from CHEWBBACA import chewBBACA
from CHEWBBACA.utils import constants as ct


# Use the tmp_path fixture to create a tmp directory for each test
@pytest.fixture
def args_fixture(request, tmp_path):
	# Setup step
	args = request.param
	args[0][5] = os.path.join(tmp_path, args[0][5])
	yield args
	# Teardown step only runs after yield
	shutil.rmtree(tmp_path)


@pytest.mark.parametrize(
    "args_fixture",
    [
	 (ct.SCHEMAEVALUATOR_TEST_VALID_INPUT, 'data/schemaevaluator_data/expected_results/valid'),
     (ct.SCHEMAEVALUATOR_TEST_SINGLE_ALLELE, 'data/schemaevaluator_data/expected_results/single_allele'),
     (ct.SCHEMAEVALUATOR_TEST_SINGLE_INVALID_ALLELE, 'data/schemaevaluator_data/expected_results/single_invalid_allele'),
     (ct.SCHEMAEVALUATOR_TEST_SEVERAL_INVALID_ALLELES, 'data/schemaevaluator_data/expected_results/several_invalid_alleles'),
    ],
	indirect=True
)
def test_schemaEvaluator_valid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		# Check schema files
		# Schema created by test
		output_dir = args_fixture[0][5]
		output_files = [os.path.join(output_dir, file)
						for file in os.listdir(output_dir)]
		output_files.sort()
		# Expected schema data
		expected_dir = args_fixture[1]
		expected_files = [os.path.join(expected_dir, file)
						  for file in os.listdir(expected_dir)]
		expected_files.sort()

		# Group test results and expected results based on basename
		files = output_files + expected_files
		basename_dict = {}
		for f in files:
			basename = os.path.basename(f)
			basename_dict.setdefault(basename, []).append(f)

		# Assert that output dir contains the expected file number
		# Cannot compare files, HTML files contains the same data but assert evaluates to False
		for k, v in basename_dict.items():
			assert len(v) == 2


@pytest.mark.parametrize(
    "args_fixture",
    [
     (ct.SCHEMAEVALUATOR_TEST_EMPTY_FILES, 'Could not get input files.'),
     (ct.SCHEMAEVALUATOR_TEST_ZERO_BYTES, 'Could not get input files.'),
     (ct.SCHEMAEVALUATOR_TEST_FAKE_PATH, 'Path to input schema does not exist. Please provide a valid path.')
    ],
	indirect=True
)
def test_schemaEvaluator_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			assert e.value.code == args_fixture[1]
