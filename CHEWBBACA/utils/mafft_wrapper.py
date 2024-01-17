#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with the execution
of the MAFFT software (https://mafft.cbrc.jp/alignment/software/).

Code documentation
------------------
"""


import os
import shutil
import subprocess

try:
    from utils import iterables_manipulation as im
except ModuleNotFoundError:
    from CHEWBBACA.utils import iterables_manipulation as im


def call_mafft(input_file, output_file):
    """Call MAFFT to generate an alignment.

    Parameters
    ----------
    input_file : str
        Path to a FASTA file with the sequences to align.

    Returns
    -------
    Path to the file with the computed MSA if successful, False otherwise.
    """
    mafft_cmd = [shutil.which('mafft'), '--thread', '1', '--treeout', '--retree', '1',
                 '--maxiterate', '0', input_file, '>', output_file]
    mafft_cmd = ' '.join(mafft_cmd)

    mafft_cmd = subprocess.Popen(mafft_cmd,
                                 shell=True,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
    mafft_cmd.wait()

    return os.path.exists(output_file)
