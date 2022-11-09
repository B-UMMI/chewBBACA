#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions to work with iterables.

Code documentation
------------------
"""


import re
import hashlib
import itertools

try:
    from utils import (constants as ct,
                       fasta_operations as fao)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 fasta_operations as fao)


def match_regex(text, pattern):
    """Extract substrings that match a regex pattern.

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
    """Join all elements in a list into a single string.

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
    """Flatten one level of a nested list.

    Parameters
    ----------
    list_to_flatten : list
        Nested list to flatten.

    Returns
    -------
    flattened_list : str
        Input list flattened by one level.
    """
    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list


def isListEmpty(input_list):
    """Check if a nested list is empty."""
    if isinstance(input_list, list):
        return all(map(isListEmpty, input_list)) if isinstance(input_list, list) else False


def divide_list_into_n_chunks(list_to_divide, n):
    """Divides a list into a defined number of sublists.

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
    sublists = [line for line in sublists if len(line) > 0]

    return sublists


def merge_dictionaries(dictionaries_list, overwrite=False):
    """Merge several dictionaries.

    Parameters
    ----------
    dictionaries_list : list
        A list with the dictionaries to merge.
    overwrite : bool
        True to overwrite values when there is a key
        collision, False otherwise.

    Returns
    -------
    merged_dicts : dict
        A dictionary resulting from merging
        all input dictionaries.
    """
    merged_dicts = dictionaries_list[0]
    if overwrite is True:
        for d in dictionaries_list[1:]:
            merged_dicts = {**merged_dicts, **d}
    elif overwrite is False:
        for d in dictionaries_list[1:]:
            for k, v in d.items():
                if k in merged_dicts:
                    merged_dicts[k] = list(set.union(set(v), set(merged_dicts[k])))
                else:
                    merged_dicts[k] = v

    return merged_dicts


def invert_dictionary(dictionary):
    """Invert a dictionary (key:value to value:key).

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
    """Split a dictionary.

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
    # create iterable from input object
    # need to convert or slicing will always return the first elements
    it = iter(iterable)
    for i in range(0, len(iterable), size):
        chunks.append({k: iterable[k] for k in itertools.islice(it, size)})

    return chunks


def select_keys(dictionary, size):
    """Select dictionary keys based on the number of values.

    Parameters
    ----------
    dictionary : dict
        Dictionary with key:value pairs to select.
    size : int
        Number of values used to select dictionary entries.

    Returns
    -------
    selected : list
        List with the dictionary keys that have a number
        of matching values equal or greater/lesser than
        `size`.
    """
    selected = [k for k, v in dictionary.items() if len(v) == size]

    return selected


def prune_dictionary(dictionary, keys):
    """Remove a set of keys from a dictionary.

    Parameters
    ----------
    dictionary : dict
        Input dictionary.
    keys : list
        List of keys to remove.

    Returns
    -------
    pruned_dict : dict
        Dictionary without the keys in the list
        of keys to remove.
    """
    pruned_dict = {k: v for k, v in dictionary.items() if k not in keys}

    return pruned_dict


def extract_subsequence(sequence, start, stop):
    """Extract a substring from a string.

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
    """Extract a coding sequence from a larger sequence.

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
        Coding sequence extracted from input contig.
    """
    coding_sequence = extract_subsequence(sequence, start, stop)

    if strand == 0:
        coding_sequence = reverse_complement(coding_sequence, ct.DNA_BASES)

    return coding_sequence


def escape_special_characters(input_string):
    """Escape strings to use in regex.

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
    """Replace multiple characters in a string.

    Parameters
    ----------
    input_string : str
        String with characters to be replaced.
    replacements : list
        List that contains sublists with two elements,
        the element to be replaced and the element it should
        be replaced by.

    Returns
    -------
    replaced : str
        Input string with replaced characters.
    """
    for r in replacements:
        if r[0] in input_string:
            input_string = input_string.replace(*r)

    return input_string


def reverse_complement(input_string, alphabet):
    """Determine the reverse complement of a string.

    Parameters
    ----------
    input_string : str
        Input string.
    alphabet : str
        String that contains the alphabet used in
        the input string (the reverse must be equal
        to the desired complement).

    Returns
    -------
    reverse_complement_string : str
        The reverse complement of the input string.
    """
    translation_table = str.maketrans(alphabet, alphabet[::-1])

    upper_string = input_string.upper()
    complement_string = upper_string.translate(translation_table)
    reverse_complement_string = reverse_str(complement_string)

    return reverse_complement_string


def reverse_str(input_string):
    """Reverse character order in a string.

    Parameters
    ----------
    input_string : str
     String to be reversed.

    Returns
    -------
    revstr : str
        Reverse of input string.
    """
    revstr = input_string[::-1]

    return revstr


def check_str_alphabet(input_string, alphabet):
    """Determine if a string only contains characters from specified alphabet.

    Parameters
    ----------
    input_string : str
        Input string.
    alphabet : str
        String that includes all characters in the
        alphabet.

    Returns
    -------
    True if sequence only has characters from specified
    alphabet, False otherwise.
    """
    alphabet_chars = set(alphabet)
    string_chars = set(input_string)

    diff = string_chars - alphabet_chars

    return len(diff) == 0


def check_str_multiple(input_string, number):
    """Determine if length of input string is a multiple of a specified number.

    Parameters
    ----------
    input_string : str
        Input string.
    number : int
        Length value should be a multiple of this number.

    Returns
    -------
    True if the length of the sequence is a multiple of the
    specified number, False otherwise.
    """
    return (len(input_string) % number) == 0


def hash_sequence(input_string, hash_type='sha256'):
    """Compute hash of an input string.

    Parameters
    ----------
    input_string : str
        Input string to hash.
    hash_type : str
        Hash type/function that will be used to compute the
        hash (any of the hash functions available in the
        hashlib module).

    Returns
    -------
    hashed_string : str
        String representation of the HASH object
        in hexadecimal digits.
    """
    # get hash function object from hashlib
    hashing_function = getattr(hashlib, hash_type)

    # default encoding is UTF-8
    hashed_string = hashing_function(input_string.encode()).hexdigest()

    return hashed_string


def string_kmerizer(input_string, k_value, offset=1, position=False):
    """Decompose a string into k-mers.

    Parameters
    ----------
    input_string : str
        String to divide into k-mers.
    k_value : int
        Value for the size of k-mers.
    offset : int
        Value to indicate offset of consecutive k-mers.
    position : bool
        If the start position of the k-mers in the string
        should be stored.

    Returns
    -------
    kmers : list
        List that contains the k-mers determined for the
        input string. The list will contain strings if
        it is not specified that positions should be
        stored and tuples of k-mer and start position
        if the position is stored.
    """
    if position is False:
        kmers = [input_string[i:i+k_value]
                 for i in range(0, len(input_string)-k_value+1, offset)]
    elif position is True:
        kmers = [(input_string[i:i+k_value], i)
                 for i in range(0, len(input_string)-k_value+1, offset)]

    return kmers


def determine_minimizers(input_string, adjacent_kmers, k_value, offset=1,
                         position=False):
    """Determine minimizers for a input string.

    Determines minimizers for a string based on lexicographical
    order. Skips windows that cannot have a minimizer based on
    the minimizer computed in the previous iteration.

    Parameters
    ----------
    input_string : str
        String representing the sequence.
    adjacent_kmers : int
        Window size value. Number of adjacent k-mers per group.
    k_value : int
        Value of k for the k-mer size.
    offset : int
        Value to indicate offset of consecutive k-mers.
    position : bool
        If the start position of the k-mers in the sequence
        should be stored.

    Returns
    -------
    minimizers : list
        A list with the set of minimizers determined
        for the input string.
    """
    # break string into k-mers
    kmers = string_kmerizer(input_string, k_value, offset, position)

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
    """Decode bytes objects.

    Decodes bytes objects in the input list and strips decoded
    strings from whitespaces and newlines.

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
    """Sort an iterable.

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


def filter_list(lst, remove):
    """Remove elements from a list.

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


def find_missing(integer_list):
    """Find the set of integers missing from a list to make it consecutive.

    Parameters
    ----------
    integer_list : list
        List containing integers.

    Returns
    -------
    missing_integers : list
        Sorted list of missing integers.
    """
    first = integer_list[0]
    last = integer_list[-1]
    missing_integers = sorted(set(range(first, last + 1)).difference(integer_list))

    return missing_integers


# add new parameter that accepts method to sample kmers
def kmer_index(fasta_file, word_size, fasta=True):
    """Create a k-mer index from a set of sequences in a FASTA file.

    Parameters
    ----------
    fasta_file : str
        Path to a Fasta file.
    word_size : int
        Value k for the k-mer size.
    fasta : bool
        True if input is a FASTA file, False if it is a dictionary
        with sequence identifiers as keys and sequences as values.

    Returns
    -------
    kmers_mapping : dict
        Dictionary with k-mers as keys and the
        list of sequence identifiers of the
        sequences that contain the kmers as
        values.
    seqs_kmers : dict
        Dictionary with sequence identifiers
        as keys and the set of distinct k-mers
        for each sequence as values.
    """
    if fasta is True:
        sequences = fao.sequence_generator(fasta_file)
    else:
        sequences = fasta_file

    kmers_mapping = {}
    for record in sequences:
        try:
            seqid = record
            sequence = sequences[record]
        except Exception as e:
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
    """Find elements in an iterable that contain any term from a set of terms.

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
    """Create mapping between list elements and consecutive integers.

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
        sequential integers.
    """
    mapping = {v: i+1 for i, v in enumerate(values)}
    if inverse is True:
        mapping = invert_dictionary(mapping)

    return mapping


def multiprocessing_inputs(inputs, common_args, function):
    """Create lists of inputs to process in parallel.

    Creates lists of inputs for the `map_async_parallelizer`
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
    input_groups = []
    # create a list for each distinct input
    for g in inputs:
        new_input = g + common_args + [function]
        input_groups.append(new_input)

    return input_groups


def aggregate_iterables(iterables):
    """Aggregate elements with same index from several iterables.

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
    """Create mapping between input elements and function result.

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
    """Use the polyline algorithm to encode/compress a list of numbers.

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
    # to ensure proper character display, encoded values are summed with 63
    # (the ASCII character '?')
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

        # process 5-bit chunks from right to left until we get a value that
        # is smaller than 0b100000/0x20/32
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
    """Decode a single number from a string created with polyline encoding.

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
    """Decode a list of integers compressed with polyline encoding.

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


def replace_list_values(input_list, replace_dict):
    """Replace values in list based on provided substitutions.

    Parameters
    ----------
    input_list : list
        List with values to substitute.
    replace_dict : dict
        Mapping between the values to substitute
        and the values to substitute by.

    Returns
    -------
    replaced_list : list
        List with substituted values (values
        that are not in `replace_dict` are kept
        unchanged).
    """
    replaced_list = [replace_dict.get(e, e) for e in input_list]

    return replaced_list


def replace_chars(column, missing_replace='0'):
    """Replace all non-numeric characters in a column with allele identifiers.

    Parameters
    ----------
    column : pandas.core.series.Series
        Pandas dataframe column.

    Returns
    -------
    replace_missing : pandas.core.series.Series
        Input column with cells that only contain
        numeric characters.
    """
    # remove 'INF-' from inferred alleles
    replace_inf = column.replace(to_replace='INF-',
                                 value='', regex=True)
    # replace '*' in novel alleles from schemas in Chewie-NS
    # before replacing missing data cases to avoid replacing '*' with '0'
    replace_inf = replace_inf.replace(to_replace='\*',
                                      value='', regex=True)
    # replace missing data with
    replace_missing = replace_inf.replace(to_replace='\D+.*',
                                          value=missing_replace, regex=True)

    return replace_missing
