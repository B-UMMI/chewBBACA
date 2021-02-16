#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: rfm
"""


import pickle
from itertools import product

from Bio.SubsMat.MatrixInfo import blosum62

try:
    from utils import (constants as cnst,
                       str_utils as su)
except:
    from CHEWBBACA.utils import (constants as cnst,
                                 str_utils as su)


# list of aminoacids
aminoacids = set([k[0] for k in blosum62])

# create all possible 3-mers by generating permutations
# that allow repeated elements
aminoacids_permutations = list(product(aminoacids, repeat=3))

# get possible substitutions for all 3-mers positions
possible_substitutions = {}
for k in aminoacids_permutations:
    # identify neutral and positive substitutions for each position
    substitutions = []
    for p in k:
        left_substitutions = [s[1]
                              for s, score in blosum62.items()
                              if s[0] == p and score >= 0]
        right_substitutions = [s[0]
                               for s, score in blosum62.items()
                               if s[1] == p and score >= 0]
        substitutions.append(set(left_substitutions+right_substitutions))
    
    possible_substitutions[k] = substitutions


# create all possible 3-mers based on possible substitutions
matrix = {}
for k, v in possible_substitutions.items():
    vanilla_scores = [blosum62[(a, a)] for a in k]
    max_score = max(vanilla_scores)
    main_score = sum(vanilla_scores)
    mers3 = list(product(*v))
    # determine score for each 3-mer
    mers3_scores = []
    for m in mers3:
        # blosum62 dict only has entries below and including diagonal
        scores = []
        for a in range(len(m)):
            score = blosum62.get((k[a], m[a]), None)
            if score is None:
                score = blosum62.get((m[a], k[a]))
            scores.append(score)

        mers3_scores.append(sum(scores))

    # keep only 3-mers whose score is at least
    valid_scores = [mers3[i]
                    for i in range(len(mers3))
                    if mers3_scores[i] >= 11]  # BLASTp minimum score is 11
    matrix[k] = valid_scores


# create final matrix by joining all elements in tuples
spaced_kmers_matrix = {}
for k, v in matrix.items():
    kmer = ''.join(k)
    # key values as set to ensure fast lookup time
    spaced_kmers = set([''.join(s) for s in v])
    spaced_kmers_matrix[kmer] = spaced_kmers

# save pre-computed matrix
#matrix_file = cnst.SPACED_MATRIX
matrix_file = '/home/rfm/Desktop/rfm/Cloned_repos/chewBBACA/CHEWBBACA/utils/spaced_matrix'
with open(matrix_file, 'wb') as outfile:
    pickle.dump(spaced_kmers_matrix, outfile)






















