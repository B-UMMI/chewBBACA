#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with the execution
of the FastTree software (http://www.microbesonline.org/fasttree/).

Code documentation
------------------
"""


import subprocess


def call_fasttree(alignment_file, tree_file):
	"""Compute a phylogenetic tree based on MSA data.

	Parameters
	----------
	alignment_file : str
		Path to a file with a MSA.
	tree_file : str
		Path to the output tree file in Newick format.
	"""
	proc = subprocess.Popen(['FastTree', '-fastest', '-nosupport',
							 '-noml', '-out', tree_file, alignment_file],
							stdout=subprocess.PIPE,
							stderr=subprocess.PIPE)

	# Read the stdout from FastTree
	stdout = proc.stdout.readlines()
	stderr = proc.stderr.readlines()

	return [stdout, stderr]
