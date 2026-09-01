#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the PredictGenes module works as expected.
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
	args[0][5] = os.path.join(tmp_path, args[0][7])
	yield args
	# Teardown step only runs after yield
	shutil.rmtree(tmp_path)


@pytest.mark.parametrize(
	"args_fixture",
	[
	 (ta.PREDICTGENES_TEST_DEFAULT, 'data/predictgenes_data/expected_results/default'),
	 (ta.PREDICTGENES_TEST_TRAINING_FILE, 'data/predictgenes_data/expected_results/training_file'),
	 (ta.PREDICTGENES_TEST_JUST_TRAINING, 'data/predictgenes_data/expected_results/just_training'),
	 (ta.PREDICTGENES_TEST_MULTIPLE_OUTFORMATS, 'data/predictgenes_data/expected_results/multiple_outformats'),
	],
	indirect=True # Pass parameters through args_fixture fixture
)
def test_predictgenes_valid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		output_dir = args_fixture[0][5]
		output_files = fo.list_files_recursively(output_dir)

		# Get paths to files with expected results
		expected_dir = args_fixture[1]
		expected_files = fo.list_files_recursively(expected_dir)

		# Group test results and expected results based on basename
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
	 (ta.PREDICTGENES_TEST_NO_CDS, ct.CANNOT_PREDICT)
	],
	indirect=True
)
def test_predictgenes_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			assert args_fixture[1] in e.value.code
