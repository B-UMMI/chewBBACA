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
    """
    """
    proc = subprocess.Popen(['FastTree', '-fastest', '-nosupport',
                             '-noml', '-out', tree_file, alignment_file],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

    # Read the stdout from FastTree
    stdout = proc.stdout.readlines()
    stderr = proc.stderr.readlines()

    return [stdout, stderr]
