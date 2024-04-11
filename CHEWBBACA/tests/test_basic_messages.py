#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify the commands to select the main modules work.
"""


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
	yield args
	# Teardown step only runs after yield
	shutil.rmtree(tmp_path)


@pytest.mark.parametrize(
        "args_fixture",
        [
		 (ct.CHEWIE_TEST_V, 0),
		 (ct.CHEWIE_TEST_VERSION, 0),
		 (ct.CHEWIE_TEST_H, 0),
		 (ct.CHEWIE_TEST_HELP, 0),
		 (ct.CREATESCHEMA_TEST_H, 0),
		 (ct.CREATESCHEMA_TEST_HELP, 0),
		 (ct.ALLELECALL_TEST_H, 0),
		 (ct.ALLELECALL_TEST_HELP, 0),
		 (ct.SCHEMAEVALUATOR_TEST_H, 0),
		 (ct.SCHEMAEVALUATOR_TEST_HELP, 0),
		 (ct.ALLELECALL_EVALUATOR_TEST_H, 0),
		 (ct.ALLELECALL_EVALUATOR_TEST_HELP, 0),
		 (ct.EXTRACTCGMLST_TEST_H, 0),
		 (ct.EXTRACTCGMLST_TEST_HELP, 0),
		 (ct.REMOVEGENES_TEST_H, 0),
		 (ct.REMOVEGENES_TEST_HELP, 0),
		 (ct.PREPEXTERNALSCHEMA_TEST_H, 0),
		 (ct.PREPEXTERNALSCHEMA_TEST_HELP, 0),
		 (ct.JOINPROFILES_TEST_H, 0),
		 (ct.JOINPROFILES_TEST_HELP, 0),
		 (ct.UNIPROTFINDER_TEST_H, 0),
		 (ct.UNIPROTFINDER_TEST_HELP, 0),
		 (ct.DOWNLOADSCHEMA_TEST_H, 0),
		 (ct.DOWNLOADSCHEMA_TEST_HELP, 0),
		 (ct.LOADSCHEMA_TEST_H, 0),
		 (ct.LOADSCHEMA_TEST_HELP, 0),
		 (ct.SYNCSCHEMA_TEST_H, 0),
		 (ct.SYNCSCHEMA_TEST_HELP, 0),
		 (ct.NSSTATS_TEST_H, 0),
		 (ct.NSSTATS_TEST_HELP, 0),
		 (ct.FAKEMODULE_TEST_H, 1), # Test with module that does not exist
		 (ct.FAKEMODULE_TEST_HELP, 1)
        ],
		indirect=True
)
def test_help_codes(monkeypatch, args_fixture):
    # Add args to sys.argv
	with monkeypatch.context() as m:
		m.setattr(sys, 'argv', args_fixture[0])
		with pytest.raises(SystemExit) as e:
			chewBBACA.main()

			assert e.type == SystemExit
			assert e.value.code == args_fixture[1]
