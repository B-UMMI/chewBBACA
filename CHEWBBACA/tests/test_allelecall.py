#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify that the AlleleCall module works as expected.
"""


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
	 (ta.ALLELECALL_TEST_DEFAULT, 'data/allelecall_data/test_results/mode4'),
	 (ta.ALLELECALL_TEST_MODE1, 'data/allelecall_data/test_results/mode1'),
	 (ta.ALLELECALL_TEST_MODE2, 'data/allelecall_data/test_results/mode2'),
	 (ta.ALLELECALL_TEST_MODE3, 'data/allelecall_data/test_results/mode3'),
	 (ta.ALLELECALL_TEST_MODE4, 'data/allelecall_data/test_results/mode4'),
	 (ta.ALLELECALL_TEST_CDS_DEFAULT, 'data/allelecall_data/test_results/cds_input_mode4'),
	 (ta.ALLELECALL_TEST_CDS_MODE1, 'data/allelecall_data/test_results/cds_input_mode1'),
	 (ta.ALLELECALL_TEST_CDS_MODE2, 'data/allelecall_data/test_results/cds_input_mode2'),
	 (ta.ALLELECALL_TEST_CDS_MODE3, 'data/allelecall_data/test_results/cds_input_mode3'),
	 (ta.ALLELECALL_TEST_CDS_MODE4, 'data/allelecall_data/test_results/cds_input_mode4'),
	 (ta.ALLELECALL_TEST_GENOME_LIST, 'data/allelecall_data/test_results/mode4'),
	 (ta.ALLELECALL_TEST_LOCI_IDS_EXTENSION, 'data/allelecall_data/test_genes_list/test_genes_results'),
	 (ta.ALLELECALL_TEST_LOCI_IDS_NOEXTENSION, 'data/allelecall_data/test_genes_list/test_genes_results'),
	 (ta.ALLELECALL_TEST_LOCI_PATHS, 'data/allelecall_data/test_genes_list/test_genes_results'),
	],
	indirect=True # Pass parameters through args_fixture fixture
)
def test_allelecall_valid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		chewBBACA.main()

		# Get paths to output files
		output_dir = args_fixture[0][7]
		output_files = [os.path.join(output_dir, file)
						for file in os.listdir(output_dir)
						if file != 'logging_info.txt']

		# Get paths to files with expected results
		expected_dir = args_fixture[1]
		expected_files = [os.path.join(expected_dir, file)
						for file in os.listdir(expected_dir)]

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
	 (ta.ALLELECALL_TEST_EMPTY_DIR, ct.MISSING_FASTAS_EXCEPTION),
	 (ta.ALLELECALL_TEST_EMPTY_FILES, ct.MISSING_FASTAS_EXCEPTION),
	 (ta.ALLELECALL_TEST_ZERO_BYTES, ct.MISSING_FASTAS_EXCEPTION),
	 (ta.ALLELECALL_TEST_FAKE_PATH, ct.INVALID_INPUT_PATH),
	 (ta.ALLELECALL_TEST_BLANK_SPACE, ct.INPUTS_INCLUDE_BLANKS[:46]),
	 (ta.ALLELECALL_TEST_LONG_PREFIX, ct.INPUTS_LONG_PREFIX[:65]),
	 (ta.ALLELECALL_TEST_SAME_PREFIX, ct.INPUTS_SHARE_PREFIX[:56])
	],
	indirect=True
)
def test_allelecall_invalid_input(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		# Create empty dir for empty dir test
		if args_fixture[0][3] == 'empty_dir':
			# Create inside output directory that is deleted in teardown step
			empty_dirname = os.path.dirname(args_fixture[0][7])
			empty_dirpath = os.path.join(empty_dirname, args_fixture[0][3])
			os.mkdir(empty_dirpath)
			args_fixture[0][3] = empty_dirpath
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			assert args_fixture[1] in e.value.code


@pytest.mark.parametrize(
	"args_fixture",
	[
	 (ta.ALLELECALL_TEST_DEFAULT+['--bsr', '-1'], ct.INVALID_BSR),
	 (ta.ALLELECALL_TEST_DEFAULT+['--bsr', '1.1'], ct.INVALID_BSR),
	 (ta.ALLELECALL_TEST_DEFAULT+['--bsr', 'sus'], ct.INVALID_BSR_TYPE.format('sus')),
	 (ta.ALLELECALL_TEST_DEFAULT+['--l', '-1'], ct.INVALID_MINLEN),
	 (ta.ALLELECALL_TEST_DEFAULT+['--l', 'sus'], ct.INVALID_MINLEN_TYPE),
	 (ta.ALLELECALL_TEST_DEFAULT+['--st', '-1'], ct.INVALID_ST),
	 (ta.ALLELECALL_TEST_DEFAULT+['--st', '1.1'], ct.INVALID_ST),
	 (ta.ALLELECALL_TEST_DEFAULT+['--st', 'sus'], ct.INVALID_ST_TYPE)
	],
	indirect=True
)
def test_allelecall_invalid_args(monkeypatch, args_fixture):
	# Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			assert e.value.code == args_fixture[1]
