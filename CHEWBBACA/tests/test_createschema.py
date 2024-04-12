#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the CreateSchema module works as expected.
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


# Test successful process
# 1) With a path to a folder that contains FASTA files
# 2) With a path to a file that contains a list of paths to FASTA files
@pytest.mark.parametrize(
    "args_fixture",
    [
	 (ta.CREATESCHEMA_TEST_GENOME_TEMPLATE, 'data/createschema_data/expected_results'),
     (ta.CREATESCHEMA_TEST_GENOME_LIST, 'data/createschema_data/expected_results'),
	 (ta.CREATESCHEMA_TEST_CDS, 'data/createschema_data/expected_results'),
    ],
	indirect=True
)
def test_createschema_valid_input(monkeypatch, args_fixture):
    # Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		# Check schema files
		# Schema created by test
		output_dir = args_fixture[0][5]
		schema_dir = os.path.join(output_dir, 'schema_seed')
		output_files = [os.path.join(schema_dir, file)
						for file in os.listdir(schema_dir)
						if 'short' not in file]
		output_files.sort()
		# Expected schema data
		expected_dir = os.path.join(args_fixture[1], 'schema_seed')
		expected_files = [os.path.join(expected_dir, file)
						  for file in os.listdir(expected_dir)
						  if 'short' not in file]
		expected_files.sort()

		# Get config values
		genes_lists = [output_files.pop(0), expected_files.pop(0)]
		schemas_configs = [output_files.pop(0), expected_files.pop(0)]

		# Compare lists of loci
		assert fo.pickle_loader(genes_lists[0]).sort() == fo.pickle_loader(genes_lists[1]).sort()
		# Read config values
		# Ignore chewBBACA version value
		configs1 = fo.pickle_loader(schemas_configs[0])
		del configs1['chewBBACA_version']
		configs2 = fo.pickle_loader(schemas_configs[1])
		del configs2['chewBBACA_version']
		assert configs1 == configs2

		# Group test schema and expected schema files based on basename
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
	 (ta.CREATESCHEMA_TEST_EMPTY_FILES, ct.MISSING_FASTAS_EXCEPTION),
     (ta.CREATESCHEMA_TEST_ZERO_BYTES, ct.MISSING_FASTAS_EXCEPTION),
     (ta.CREATESCHEMA_INVALID_PTF_PATH, ct.INVALID_PTF_PATH),
     (ta.CREATESCHEMA_TEST_HEADER_ONLY, ct.CANNOT_PREDICT),
     (ta.CREATESCHEMA_TEST_INVALID_GENOME, ct.CANNOT_PREDICT)
    ],
	indirect=True
)
def test_createschema_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			# Check that the exit message includes expected message
			assert args_fixture[1] in e.value.code
