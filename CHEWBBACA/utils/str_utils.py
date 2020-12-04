#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


import re
import hashlib


def extract_subsequence(sequence, start, stop):
    """ Extract substring from string.

        Parameters
        ----------
        sequence : str
            Input string.
        start : int
            Substring start position in input string.
        stop : int
            Substring stop position in input string.

        Returns
        -------
        subsequence : str
            Substring extracted from input string.
    """

    subsequence = sequence[start:stop]

    return subsequence


def extract_single_cds(sequence, start, stop, strand):
    """ Extract coding sequence from contig.

        Parameters
        ----------
        sequence : str
            Contig that contains the coding sequence.
        start : int
            Coding sequence start position in contig.
        stop : int
            Coding sequence stop position in contig.
        strand : int
            Conding sequence orientation.

        Returns
        -------
        coding_sequence : str
            Coding sequence extracted from input contig
            in sense orientation.
    """

    coding_sequence = extract_subsequence(sequence, start, stop)

    if strand == 0:
        coding_sequence = reverse_complement(coding_sequence)

    return coding_sequence


def escape_special_characters(input_string):
    """ Escapes strings to use in regex.

        Parameters
        ----------
        input_string : str
            String containing characters to escape.

        Returns
        -------
        escaped_string : str
            Escaped string.
    """

    escaped_string = re.escape(input_string)

    return escaped_string


def replace_multiple_characters(input_string):
    """ Replaces multiple characters in a string.

        Parameters
        ----------
        input_string : str
            String with characters to be replaced.

        Returns
        -------
        replaced : str
            Input string without replaced characters.
    """

    replaced = input_string.replace("|", "_")\
                           .replace("_", "-")\
                           .replace("(", "")\
                           .replace(")", "")\
                           .replace("'", "")\
                           .replace("\"", "")\
                           .replace(":", "")

    return replaced


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.

        Returns
        -------
        reverse_complement_strand : str
            The reverse complement of the DNA sequence (lowercase
            is converted to uppercase).
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A',
                       'n': 'N', 'N': 'N'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    # default to 'N' if nucleotide is not in base_complement dictionary
    bases = [base_complement.get(base, 'N') for base in bases]

    complement_strand = ''.join(bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def reverse_str(string):
    """ Reverse character order in input string.

        Parameters
        ----------
        string : str
         String to be reversed.

        Returns
        -------
        revstr : str
            Reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def check_str_alphabet(string, alphabet):
    """ Determine if a string only contains characters from
        specified alphabet.

        Parameters
        ----------
        string : str
            Input string.
        alphabet : str
            String that has all characters from desired
            alphabet.

        Returns
        -------
        "True" if sequence only has characters from specified
        alphabet and string "ambiguous or invalid characters" if
        it any of its characters is not in the alphabet.
    """

    valid_chars = alphabet
    if all(n in valid_chars for n in string) is True:
        return True
    else:
        return 'ambiguous or invalid characters'


def check_str_multiple(string, number):
    """ Determine if length of input string is multiple of
        a specified number.

        Parameters
        ----------
        string : str
            Input string.
        number : int
            Length value should be a multiple of this number.

        Returns
        -------
        "True" if the length of the sequence is a multiple of the
        specified number and "sequence length is not a multiple of number"
        if condition is not satisfied.
    """

    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def hash_sequence(string):
    """ Compute SHA256 for an input string.

        Parameters
        ----------
        string : str
            Input string to hash.

        Returns
        -------
        sha256 : str
            String representation of the sha256 HASH object.
    """

    sha256 = hashlib.sha256(string.encode('utf-8')).hexdigest()

    return sha256


def sequence_kmerizer(sequence, k_value, offset=1, position=False):
    """ Decomposes a sequence into kmers.

        Parameters
        ----------
        sequence : str
            Sequence to divide into kmers.
        k_value : int
            Value for the size k of kmers.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        kmers : list
            List with the kmers determined for the input
            sequence. The list will contain strings if
            it is not specified that positions should be
            stored and tuples of kmer and start position
            if the position is stored.
    """

    if position is False:
        kmers = [sequence[i:i+k_value]
                 for i in range(0, len(sequence)-k_value+1, offset)]
    elif position is True:
        kmers = [(sequence[i:i+k_value], i)
                 for i in range(0, len(sequence)-k_value+1, offset)]

    return kmers


def determine_minimizers(sequence, adjacent_kmers, k_value, position=False):
    """ Determine the minimizers for a sequence based on
        lexicographical order. Skips windows that
        cannot have a minimizer based on the minimizer
        computed in the previous iteration.

        Parameters
        ----------
        sequence : str
            String representing the sequence.
        adjacent_kmers : int
            Window size value. Number of adjacent kmers per group.
        k_value : int
            Value of k for kmer size.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.

        Returns
        -------
        minimizers : list
            A list with the set of minimizers determined
            for the input sequence.
    """

    # break sequence into kmers
    kmers = sequence_kmerizer(sequence, k_value, position=position)

    i = 0
    previous = None
    sell = False
    minimizers = []
    # determine total number of windows
    last_window = (len(kmers)-adjacent_kmers)
    while i <= last_window:
        # get kmers in current window
        window = kmers[i:i+adjacent_kmers]
        # pick smallest kmer as minimizer
        minimizer = [min(window)]
        # get position in window of smallest minimizer
        minimizer_idx = window.index(minimizer[0])
        # sliding window that does not included last minimizer
        if previous is None:
            # simply store smallest minimizer
            minimizers.extend(minimizer)
        # sliding window includes last minimizer because we
        # skipped some sliding windows
        else:
            # check if minimizer is the same as the one picked
            # in the last window
            # Do not store minimizer if it is the same
            if minimizer[0] != previous:
                # get kmers smaller than last minimizer
                skipped = window[1:minimizer_idx]
                # determine if any of the smaller kmers is
                # the minimizer of a skipped window
                minimal = previous
                for m in skipped:
                    if m < minimal:
                        minimizer.append(m)
                        minimal = m
                minimizers.extend(minimizer)

        # slide by 1 if minimizer has index 0 in window
        if minimizer_idx == 0:
            i += 1
            previous = None
        # skip sliding windows based on minimizer position
        else:
            i += minimizer_idx
            # if adding minimizer index surpasses last window value we
            # might miss one last minimizer because it will fail the condition
            # find a better way to control this condition!
            if i > last_window and sell is False:
                i = last_window
                sell = True
            previous = minimizer[0]

    return minimizers
