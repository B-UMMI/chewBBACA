#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:16:33 2019

@author: pcerqueira
"""

try:
    from createschema import PPanGen
except:
    from CHEWBBACA.createschema import PPanGen

from contextlib import nullcontext as does_not_raise
import pytest


@pytest.mark.parametrize(
        "test_input, expectation",
        [(5, pytest.raises(TypeError)), # Tests integer input
         ("Aello World", pytest.raises(KeyError)), # Tests string input
         ("AACTGGCATGCTATGCAT", does_not_raise()) # Tests valid input to confirm that no exception is raised
         ])
def test_reverseComplement_invalid_inputs(test_input, expectation):
    """Tests the behaviour of the reverseComplement function with unexpected inputs"""
    with expectation:
        PPanGen.reverseComplement(test_input)

@pytest.mark.parametrize(
        "test_input, expected",
        [(pytest.param("AAAAAACCCCTTTTGGCTTATCGX", "", marks=pytest.mark.xfail(raises=KeyError))), # Test ambiguous bases
        ])

def test_reverseComplement_invalid_sequences(test_input, expected):
    """Tests the behaviour of the reverseComplement function with different sequence inputs"""
    assert PPanGen.reverseComplement(test_input) == expected


