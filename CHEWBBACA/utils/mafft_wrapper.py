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
from Bio.Align.Applications import MafftCommandline


def call_mafft(genefile, output_directory):
    """Call MAFFT to generate an alignment.

    Parameters
    ----------
    genefile : str
        Path to a FASTA file with the sequences to align.

    Returns
    -------
    Path to the file with the computed MSA if successful, False otherwise.
    """
    try:
        mafft_cline = MafftCommandline(input=genefile,
                                       adjustdirection=True,
                                       treeout=True,
                                       thread=1,
                                       retree=1,
                                       maxiterate=0,
                                       )
        stdout, stderr = mafft_cline()
        outfile_basename = os.path.basename(genefile)
        outfile_basename = outfile_basename.replace('.fasta', '_aligned.fasta')
        outfile = os.path.join(output_directory, outfile_basename)
        with open(outfile, 'w') as handle:
            handle.write(stdout)

        return outfile
    except Exception as e:
        print(e)
        return False
