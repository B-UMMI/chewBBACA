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
	import global_constants as gct
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct


def call_mafft(input_file, output_file, custom_mafft_params=None):
	"""Call MAFFT to compute a MSA.

	Parameters
	----------
	input_file : str
		Path to a FASTA file with the sequences to align using MAFFT.
	output_file : str
		Path to the output file created by MAFFT with the MSA.
	custom_mafft_params : list, optional
		List of custom parameters to pass to MAFFT. If None, default
		parameters will be used. Default is None.

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
	if not custom_mafft_params:
		mafft_cmd = [gct.MAFFT_ALIAS] + gct.MAFFT_DEFAULT_PARAMETERS + [input_file, '>', output_file]
	else:
		mafft_cmd = [gct.MAFFT_ALIAS] + custom_mafft_params + [input_file, '>', output_file]
	# Join command list into single string
	mafft_cmd = ' '.join(mafft_cmd)

	# Must include subprocess.PIPE to get stdout and stderr
	mafft_cmd = subprocess.Popen(mafft_cmd,
								 stdout=subprocess.PIPE,
								 stderr=subprocess.PIPE,
								 shell=True)
	stdout, stderr = mafft_cmd.communicate()

	outfile_exists = os.path.exists(output_file)

	return [output_file, outfile_exists, stdout, stderr]
