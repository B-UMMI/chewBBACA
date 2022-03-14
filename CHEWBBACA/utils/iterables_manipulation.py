#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to manipulate iterables.

Code documentation
------------------
"""


import re
import shelve
import hashlib
import itertools
from copy import deepcopy

try:
    from utils import (file_operations as fo,
                       fasta_operations as fao)
except:
    from CHEWBBACA.utils import (file_operations as fo,
                                 fasta_operations as fao)


def match_regex(text, pattern):
    """ Extracts substrings that match a regex pattern.
    
    Parameters
    ----------
    text : str
        Input text to search the pattern in.
    pattern : str
        Regular expression pattern.

    Returns
    -------
    match : str
        Matched substring or NoneType if no matches
        were found.
    """

    match = re.search(pattern, text)
    if match is not None:
        match_span = match.span()
        match = text[match_span[0]:match_span[1]]

    return match


def join_list(lst, delimiter):
    """ Joins all elements in a list into a single string.

    Parameters
    ----------
    lst : list
        List with elements to be joined.
    delimiter : str
        Character used to join list elements.

    Returns
    -------
    joined_list : str
        A single string with all elements in the input
        list joined by the character chosen as link.
    """

    joined_list = delimiter.join(lst)

    return joined_list


def flatten_list(list_to_flatten):
    """ Flattens one level of a nested list.

    Parameters
    ----------
    list_to_flatten : list
        List with nested lists.

    Returns
    -------
    flattened_list : str
        Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list


def isListEmpty(input_list):
    """ Checks if a nested list is empty. """

    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False


def divide_list_into_n_chunks(list_to_divide, n):
    """ Divides a list into a defined number of sublists.

    Parameters
    ----------
    list_to_divide : list
        List to divide into sublists.
    n : int
        Number of sublists to create.

    Returns
    -------
    sublists : list
        List with the sublists created by dividing
        the input list.
    """

    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    # exclude lists that are empty due to small number of elements
    sublists = [i for i in sublists if len(i) > 0]

    return sublists


def merge_dictionaries(starting_dictionary, dictionaries_list, overwrite=False):
    """ Merges several dictionaries.

    Parameters
    ----------
    dictionaries_list : list
        A list with the dictionaries to merge.

    Returns
    -------
    merged_dicts : dict
        A dictionary resulting from merging
        all input dictionaries.
    """

    merged_dicts = starting_dictionary
    if overwrite is True:
        for d in dictionaries_list:
            merged_dicts = {**merged_dicts, **d}
    elif overwrite is False:
        for d in dictionaries_list:
            for k, v in d.items():
                if k in merged_dicts:
                    merged_dicts[k] = list(set.union(set(v), set(merged_dicts[k])))
                else:
                    merged_dicts[k] = v

    return merged_dicts


def invert_dictionary(dictionary):
    """ Inverts a dictionary (key:value to value:key).

    Parameters
    ----------
    dictionary : dict
        Dictionary to be inverted.

    Returns
    -------
    inverted_dict : dict
        Inverted dictionary.
    """

    inverted_dict = {value: key for key, value in dictionary.items()}

    return inverted_dict


def split_iterable(iterable, size):
    """ Splits a dictionary.

    Parameters
    ----------
    iterable : dict
        Dictionary to split.
    size : int
        Size of dictionaries created from the input
        dictionary.

    Returns
    -------
    chunks : list
        List with dictionaries of defined size
        resulting from splitting the input dictionary.
    """

    chunks = []
    it = iter(iterable)
    for i in range(0, len(iterable), size):
        chunks.append({k: iterable[k] for k in itertools.islice(it, size)})

    return chunks


def select_clusters(clusters, cluster_size):
    """ Determines clusters that contain a specified number
        of sequences.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        cluster_size : int
            Number of sequences in the clusters that will
            be selected.

        Returns
        -------
        clusters_ids : list
            List with cluster identifiers for clusters that
            contain a number of sequences equal to specified
            cluster size.
    """

    clusters_ids = [k for k, v in clusters.items() if len(v) == cluster_size]

    return clusters_ids


def remove_entries(dictionary, keys):
    """ Creates new dictionary without entries with
        specified keys.

        Parameters
        ----------
        dictionary : dict
            Input dictionary.
        keys : list
            List of keys for the entries that should
            not be included in the new dictionary.

        Returns
        -------
        new_dict : dict
            Dictionary without entries with keys in
            the input list.
    """

    new_dict = {k: v for k, v in dictionary.items() if k not in keys}

    return new_dict


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
    """ Extracts a coding sequence from another sequence.

        Parameters
        ----------
        sequence : str
            Contig that contains the coding sequence.
        start : int
            Coding sequence start position in contig.
        stop : int
            Coding sequence stop position in contig.
        strand : int
            Coding sequence orientation.

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


def replace_multiple_characters(input_string, replacements):
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

    for r in replacements:
        if r[0] in input_string:
            input_string = input_string.replace(*r)

    return input_string


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

    if all(n in alphabet for n in string) is True:
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
            String representation of the SHA256 HASH object
            in hexadecimal digits.
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
            Value for the size of kmers.
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


def determine_minimizers(sequence, adjacent_kmers, k_value, offset=1,
                         position=False):
    """ Determines minimizers for a sequence based on
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
            Value of k for the kmer size.
        offset : int
            Value to indicate offset of consecutive kmers.
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
    kmers = sequence_kmerizer(sequence, k_value,
                              offset=offset, position=position)

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


def decode_str(str_list, encoding):
    """ Decodes bytes objects in the input list and
        strips decoded strings from whitespaces and
        newlines.

        Parameters
        ----------
        str_list
            List with string or bytes objects to decode
            and strip of whitespaces and newlines.
        encoding : str
            Encoding codec to use.

        Returns
        -------
        decoded : list
            List with strings without whitespaces or
            newlines.
    """

    decoded = [m.decode(encoding).strip()
               if type(m) == bytes
               else m.strip()
               for m in str_list]

    return decoded


# sorted key parameter is None by default
def sort_iterable(data, sort_key=None, reverse=False):
    """ Sorts an iterable.

    Parameters
    ----------
    data : iter
        Iterable to sort.
    sort_key
        If provided, data will be sorted based
        on this function.
    reverse : bool
        If sorting order should be inverted.

    Returns
    -------
    sorted_data
        Iterable with sorted elements.
    """

    # sorted key parameter is None by default
    sorted_data = sorted(data, key=sort_key, reverse=reverse)

    return sorted_data


def add_prefix(ids, prefix):
    """ Adds a prefix to a set of identifiers.
        Identifiers are split by underscore and
        prefix is added to last element.

        Parameters
        ----------
        ids : iter
            Iterable with identifiers
            (e.g.: list, set, dictionary keys).
        prefix : str
            Prefix to add to all identifiers.

        Returns
        -------
        ids_map : dict
            Dictionary with input identifiers as
            keys and prefixed identifiers as values.
    """

    ids_map = {}
    for i in ids:
        new_id = '{0}_{1}'.format(prefix, i.split('_')[-1])
        ids_map[i] = new_id

    return ids_map


def filter_list(lst, remove):
    """ Removes elements from a list.

    Parameters
    ----------
    lst : list
        Input list.
    remove : list
        List of elements to remove from input list.

    Returns
    -------
    filtered_list : list
        List without the removed elements.
    """

    filtered_list = list(set(lst) - set(remove))

    return filtered_list


def find_missing(lst):
    """ Finds missing integers in list
        of consecutive integers.

        Parameters
        ----------
        lst : list
            List containing consecutive integers.

        Returns
        -------
        list
            Sorted list of missing integers.
    """

    start = lst[0]
    end = lst[-1]

    return sorted(set(range(start, end + 1)).difference(lst))


def kmer_index(fasta_file, word_size):
    """ Creates a kmer index based on a set
        of sequences in a Fasta file.

    Parameters
    ----------
    fasta_file : str
        Path to a Fasta file.
    word_size : int
        Value k for the kmer size.

    Returns
    -------
    kmers_mapping : dict
        Dictionary with kmers as keys and the
        list of sequence identifiers of the
        sequences that contain the kmers as
        values.
    seqs_kmers : dict
        Dictionary with sequence identifiers
        as keys and the set of distinct kmers
        for each sequence as values.
    """

    sequences = fao.sequence_generator(fasta_file)

    kmers_mapping = {}
    for record in sequences:
        seqid = record.id
        sequence = str(record.seq)
        minimizers = determine_minimizers(sequence, word_size,
                                          word_size, position=False)
        kmers = set(minimizers)

        # create dict with kmers as keys and list
        # of sequences with given kmers as values
        for kmer in kmers:
            kmers_mapping.setdefault(kmer, []).append(seqid)

    return kmers_mapping



def contained_terms(iterable, terms):
    """ Finds elements in an iterable that contain
        any term from a set of terms.

    Parameters
    ----------
    iterable
        An iterable such as a list with strings
        or a dictionary.
    terms : list
        Terms to search for.

    Returns
    -------
    matches : list
        List with the elements of the iterable
        that contain any of the searched terms.
    """

    matches = []
    for e in iterable:
        if any([term in e for term in terms]) is True:
            matches.append(e)

    return matches


def integer_mapping(values, inverse=False):
    """ Creates dictionary with the mapping between the
        elements in the input iterable and integer values
        assigned sequentially.

    Parameters
    ----------
    values : iter
        Input iterable.
    inverse : bool
        Invert mapping order, integers:values instead
        of values:integers.

    Returns
    -------
    mapping : dict
        Dictionary with the mapping between the
        elements in the input iterable and the
        result of applying the function to each
        value.
    """

    mapping = {}
    if inverse is False:
        mapping = {v: i+1 for i, v in enumerate(values)}
    elif inverse is True:
        mapping = {i+1: v for i, v in enumerate(values)}

    return mapping


def multiprocessing_inputs(inputs, common_args, function):
    """ Creates lists of inputs for the `map_async_parallelizer`
        function from the `multiprocessing_operations` module.

    Parameters
    ----------
    inputs : list
        List with the distinct inputs.
    common_args : list
        List with the common arguments/inputs.
    function : func
        Function that will be parallelized by the
        `map_async_parallelizer` function.

    Returns
    -------
    input_groups : list
        List with one sublist per distinct input. Each
        sublist also contains the common arguments and
        the function that will receive the arguments.
    """

    input_groups = deepcopy(inputs)
    # create a list for each distinct input
    for g in input_groups:
        # add the common arguments that will be used in each function call
        g.extend(common_args)
        # add the function that will be parallelized
        g.append(function)

    return input_groups


def aggregate_iterables(iterables):
    """ Aggregates elements with same index from
        several iterables.

    Parameters
    ----------
    iterables : list
        List of iterables to aggreagate.

    Returns
    -------
    aggregated_inputs : list
        List with a sublist per group of elements
        with same index that were aggregated.
    """

    aggregated_inputs = [list(i) for i in zip(*iterables)]

    return aggregated_inputs


def mapping_function(values, function, args):
    """ Creates a dictionary with the mapping between
        the elements in the input iterable and the
        result of applying a function to those values.

    Parameters
    ----------
    values : iter
        Input iterable.
    function : func
        Funtion to apply to the elements in the
        input iterable.
    args : list
        List of arguments to pass to the function.

    Returns
    -------
    mapping : dict
        Dictionary with the mapping between the
        elements in the input iterable and the
        result of applying the function to each
        value.
    """

    mapping = {v: function(v, *args) for v in values}

    return mapping


def polyline_encoding(number_list, precision=0):
    """ Uses the polyline algorithm to encode/compress
        a list of numbers.

    Parameters
    ----------
    number_list : list
        List with numbers to encode.
    precision : int
        Number of decimal places preserved by the encoding.

    Returns
    -------
    compressed_values : str
        A single string composed of ASCII characters
        that represents the compressed list of numbers
        to the desired level of precision.

    Notes
    -----
    This implementation is based on the `numcompress` package
    (https://github.com/amit1rrr/numcompress).
    """

    compressed_values = ''
    # store the precision value
    # to ensure proper character display, encoded values are summed with 63 (the ASCII character '?')
    compressed_values += chr(precision+63)
    previous_num = 0
    for number in number_list:
        # encode the difference between numbers
        difference = number - previous_num
        # multiply and round to get integer that preserves decimal places
        difference = int(round(difference*(10**precision)))
        # left bitwise shift 1 position (0b11111 --> 0b111110)
        # and invert the encoding if the original number was negative
        # the bit added by the shift will be inverted and that allows
        # to know that the number is negative
        difference = ~(difference << 1) if difference < 0 else difference << 1

        # process 5-bit chunks from right to left until we get a value that is smaller than 0b100000/0x20/32
        while difference >= 0x20:
            # use 0x1f (0b11111) as bitmask to get the smallest 5-bit chunk
            # and 0x20 is used as continuation bit to help decode
            compressed_values += (chr((0x20 | (difference & 0x1f)) + 63))
            # right bitwise shift to exclude the 5-bit chunk that was encoded
            difference >>= 5

        # encode the last 5-bit chunk
        compressed_values += (chr(difference + 63))
        # store number to subtract from next
        previous_num = number

    return compressed_values


def decompress_number(text, index):
    """ Decode a single number from a string created through
        the compression of a list of numbers with polyline encoding.

    Parameters
    ----------
    text : str
        String representing a compressed list of numbers.
    index : int
        Index of the first character to start decoding a number.

    Returns
    -------
    Index to start decoding the next number and the number decoded
    in the current function call.
    """

    number = 0
    bitwise_shift = 0

    while True:
        # subtract 63 and remove bit from OR with 0x20 if 5-bit chunk
        # is not the last to decode a number that was in the original list
        n = (ord(text[index]) - 63)
        index += 1
        # only continue if there is a continuation bit (0b100000)
        # e.g.: 0b11111 only has 5 bits, meaning it is the last chunk
        # that is necessary to decode the number
        if n >= 0x20:
            # subtract 0x20 (0b100000) to get the 5-bit chunk
            n -= 0x20
            # contruct the binary number with biwise shift to add each
            # 5-bit chunk to original position
            number = number | (n << bitwise_shift)
            # increment bitwise shift value for next 5-bit chunk
            bitwise_shift += 5
        else:
            break

    # add the last chunk, without continuation bit, to the leftmost position
    number = number | (n << bitwise_shift)

    # invert bits to get negative number if sign bit is 1
    # remove sign bit and keep only bits for decoded value
    return index, (~number >> 1) if (number & 1) != 0 else (number >> 1)


def polyline_decoding(text):
    """ Decode a list of integers compressed with polyline
        encoding.

    Parameters
    ----------
    text : str
        String representing a compressed list of numbers.

    Returns
    -------
    number_list : list
        List with the decoded numbers.
    """

    number_list = []
    index = last_num = 0
    # decode precision value
    precision = ord(text[index]) - 63
    index += 1
    while index < len(text):
        # decode a number and get index to start decoding next number
        index, difference = decompress_number(text, index)
        # add decoded difference to get next number in original list
        last_num += difference
        number_list.append(last_num)

    number_list = [round(item * (10 ** (-precision)), precision) for item in number_list]

    return number_list



