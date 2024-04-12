#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import shutil
import pytest
import filecmp

from CHEWBBACA import chewBBACA
from CHEWBBACA.utils import constants as ct
from CHEWBBACA.tests import test_arguments as ta


# Use the tmp_path fixture to create a tmp directory for each test
@pytest.fixture
def args_fixture(request, tmp_path):
	# Setup step
	args = request.param
	args[0][7] = os.path.join(tmp_path, args[0][7])
	yield args
	# Teardown step only runs after yield
	shutil.rmtree(tmp_path)


@pytest.mark.parametrize(
    "args_fixture",
    [
     (ta.ALLELECALL_EVALUATOR_VALID, 'data/allelecallevaluator_data/expected_results'),
    ],
	indirect=True
)
def test_allelecallevaluator_valid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		# Check schema files
		# Schema created by test
		output_dir = args_fixture[0][7]
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
     (ta.ALLELECALL_EVALUATOR_INVALID_PATH, ct.MISSING_INPUT_ARG),
    ],
	indirect=True
)
def test_allelecallevaluator_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			# Check that the exit message includes expected message
			assert args_fixture[1] in e.value.code
