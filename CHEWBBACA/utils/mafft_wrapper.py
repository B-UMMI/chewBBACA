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
import subprocess

try:
    from utils import constants as ct
except ModuleNotFoundError:
    from CHEWBBACA.utils import constants as ct


def call_mafft(input_file, output_file):
    """Call MAFFT to compute a MSA.

    Parameters
    ----------
    input_file : str
        Path to a FASTA file with the sequences to align.
    output_file : str
        Path to the output file created by MAFFT.

    Returns
    -------
    output_file : str
        Path to the output file.
    outfile_exists : bool
        True if the output file was created, False otherwise.
    stdout : bytes
        MAFFT stdout.
    stderr : bytes
        MAFFT stderr.
    """
    mafft_cmd = [ct.MAFFT_ALIAS, '--thread', '1', '--treeout', '--retree', '1',
                 '--maxiterate', '0', input_file, '>', output_file]
    mafft_cmd = ' '.join(mafft_cmd)

    # Must include subprocess.PIPE to get stdout and stderr
    mafft_cmd = subprocess.Popen(mafft_cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True)
    stdout, stderr = mafft_cmd.communicate()

    outfile_exists = os.path.exists(output_file)

    return [output_file, outfile_exists, stdout, stderr]
