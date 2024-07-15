#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

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


@pytest.mark.parametrize(
		"args_fixture",
		[
		 (ta.PREPEXTERNALSCHEMA_TEST_VALID_INPUT, 'data/prep_data/expected_results_valid'),
		 (ta.PREPEXTERNALSCHEMA_TEST_EXTENSIONS, 'data/prep_data/expected_results_extension'),
		 (ta.PREPEXTERNALSCHEMA_TEST_GENE_LIST, 'data/prep_data/expected_results_extension')
		],
		indirect=True
)
def test_prep_valid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		# Get paths to output files
		output_dir = args_fixture[0][5]
		output_files = [os.path.join(output_dir, file)
						for file in os.listdir(output_dir)
						if 'short' != file]
		output_files.sort()

		# Get paths to files with expected results
		expected_dir = args_fixture[1]
		expected_files = [os.path.join(expected_dir, file)
						  for file in os.listdir(expected_dir)
						  if 'short' != file]
		expected_files.sort()

		# Get config files
		genes_lists = [output_files.pop(0), expected_files.pop(0)]
		schemas_configs = [output_files.pop(0), expected_files.pop(0)]

		# Compare configs
		assert fo.pickle_loader(genes_lists[0]).sort() == fo.pickle_loader(genes_lists[1]).sort()
		# Read config values
		# Ignore chewBBACA version value
		configs1 = fo.pickle_loader(schemas_configs[0])
		del configs1['chewBBACA_version']
		configs2 = fo.pickle_loader(schemas_configs[1])
		del configs2['chewBBACA_version']
		assert configs1 == configs2

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
		 (ta.PREPEXTERNALSCHEMA_TEST_EMPTY_DIR, ct.MISSING_FASTAS_EXCEPTION),
		 (ta.PREPEXTERNALSCHEMA_TEST_EMPTY_FILES, ct.MISSING_FASTAS_EXCEPTION),
		 (ta.PREPEXTERNALSCHEMA_TEST_ZERO_BYTES, ct.MISSING_FASTAS_EXCEPTION),
		 (ta.PREPEXTERNALSCHEMA_TEST_INVALID_PATH, ct.INVALID_INPUT_PATH),
		 (ta.PREPEXTERNALSCHEMA_TEST_BLANK_SPACE, ct.INPUTS_INCLUDE_BLANKS[:46]),
	 	 (ta.PREPEXTERNALSCHEMA_TEST_LONG_PREFIX, ct.INPUTS_LONG_PREFIX[:65]),
	 	 (ta.PREPEXTERNALSCHEMA_TEST_SAME_PREFIX, ct.INPUTS_SHARE_PREFIX[:56]),
	 	 (ta.PREPEXTERNALSCHEMA_TEST_PDB_CHAIN, ct.INPUTS_PDB_PREFIX[:70])
		],
		indirect=True
)
def test_prep_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		# Create empty dir for empty dir test
		if args_fixture[0][3] == 'empty_dir':
			# Create inside output directory that is deleted in teardown step
			empty_dirname = os.path.dirname(args_fixture[0][5])
			empty_dirpath = os.path.join(empty_dirname, args_fixture[0][3])
			os.mkdir(empty_dirpath)
			args_fixture[0][3] = empty_dirpath
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			assert args_fixture[1] in e.value.code
