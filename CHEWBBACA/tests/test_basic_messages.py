#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains tests to verify the commands to select the main modules work.
"""


import py
import os
import sys
import pytest
import filecmp
from unittest.mock import patch
#from contextlib import nullcontext as does_not_raise

from CHEWBBACA import chewBBACA


@pytest.mark.parametrize(
        'test_args, expected',
        [(['chewBBACA.py', '-h'], 0),
         (['chewBBACA.py', '--help'], 0),
         (['chewBBACA.py', '-v'], 0),
         (['chewBBACA.py', '--version'], 0),
         (['chewBBACA.py', 'CreateSchema', '-h'], 0),
         (['chewBBACA.py', 'CreateSchema', '--help'], 0),
         (['chewBBACA.py', 'AlleleCall', '-h'], 0),
         (['chewBBACA.py', 'AlleleCall', '--help'], 0),
         (['chewBBACA.py', 'SchemaEvaluator', '-h'], 0),
         (['chewBBACA.py', 'SchemaEvaluator', '--help'], 0),
		 (['chewBBACA.py', 'AlleleCallEvaluator', '-h'], 0),
         (['chewBBACA.py', 'AlleleCallEvaluator', '--help'], 0),
         (['chewBBACA.py', 'ExtractCgMLST', '-h'], 0),
         (['chewBBACA.py', 'ExtractCgMLST', '--help'], 0),
		 (['chewBBACA.py', 'RemoveGenes', '-h'], 0),
         (['chewBBACA.py', 'RemoveGenes', '--help'], 0),
		 (['chewBBACA.py', 'PrepExternalSchema', '-h'], 0),
         (['chewBBACA.py', 'PrepExternalSchema', '--help'], 0),
		 (['chewBBACA.py', 'JoinProfiles', '-h'], 0),
         (['chewBBACA.py', 'JoinProfiles', '--help'], 0),
		 (['chewBBACA.py', 'UniprotFinder', '-h'], 0),
         (['chewBBACA.py', 'UniprotFinder', '--help'], 0),
         (['chewBBACA.py', 'DownloadSchema', '-h'], 0),
         (['chewBBACA.py', 'DownloadSchema', '--help'], 0),
         (['chewBBACA.py', 'LoadSchema', '-h'], 0),
         (['chewBBACA.py', 'LoadSchema', '--help'], 0),
         (['chewBBACA.py', 'SyncSchema', '-h'], 0),
         (['chewBBACA.py', 'SyncSchema', '--help'], 0),
         (['chewBBACA.py', 'NSStats', '-h'], 0),
         (['chewBBACA.py', 'NSStats', '--help'], 0),
		 (['chewBBACA.py', 'FakeModule', '-h'], 1), # Test with module that does not exist
         (['chewBBACA.py', 'FakeModule', '--help'], 1),
         ])
def test_modules_exit_codes(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert e.value.code == expected 
