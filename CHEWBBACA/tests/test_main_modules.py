#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import py
import sys
import pytest
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
         (['chewBBACA.py', 'TestGenomeQuality', '-h'], 0),
         (['chewBBACA.py', 'TestGenomeQuality', '--help'], 0),
         (['chewBBACA.py', 'ExtractCgMLST', '-h'], 0),
         (['chewBBACA.py', 'ExtractCgMLST', '--help'], 0),
         (['chewBBACA.py', 'SchemaEvaluator', '-h'], 0),
         (['chewBBACA.py', 'SchemaEvaluator', '--help'], 0),
         (['chewBBACA.py', 'DownloadSchema', '-h'], 0),
         (['chewBBACA.py', 'DownloadSchema', '--help'], 0),
         (['chewBBACA.py', 'LoadSchema', '-h'], 0),
         (['chewBBACA.py', 'LoadSchema', '--help'], 0),
         (['chewBBACA.py', 'SyncSchema', '-h'], 0),
         (['chewBBACA.py', 'SyncSchema', '--help'], 0),
         (['chewBBACA.py', 'NSStats', '-h'], 0),
         (['chewBBACA.py', 'NSStats', '--help'], 0),
         ])
def test_modules_exit_codes(test_args, expected):
    with pytest.raises(SystemExit) as e:
        with patch.object(sys, 'argv', test_args):
            chewBBACA.main()

    assert e.type == SystemExit
    assert e.value.code == expected


#@pytest.mark.parametrize(
#        'test_args, expected',
#        [(['chewBBACA.py', 'AlleleCall', ], 0)
#         ])
#def test_allelecall(test_args, expected):
#    with patch.object(sys, 'argv', test_args):
#        capture = py.io.StdCapture()
#        chewBBACA.main()
#        stdout, stderr = capture.reset()
#
#    assert 'Select one of the following functions:' in stdout
