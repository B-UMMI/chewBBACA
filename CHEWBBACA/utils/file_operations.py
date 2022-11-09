#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related with file operations,
such as read and write files, create and delete files,
manipulate file paths, compress files, verify file contents,
etc.

Code documentation
------------------
"""


import os
import sys
import csv
import math
import time
import gzip
import shutil
import pickle
import zipfile
import urllib.request
from itertools import islice
from multiprocessing import TimeoutError
from multiprocessing.pool import ThreadPool

import pandas as pd

try:
    from utils import (constants as ct,
                       iterables_manipulation as im)
except ModuleNotFoundError:
    from CHEWBBACA.utils import (constants as ct,
                                 iterables_manipulation as im)


def file_basename(file_path, file_extension=True, delimiter='.'):
    """Extract file basename from a path.

    Parameters
    ----------
    file_path : str
        Path to the file.
    file_extension : bool
        Specify if the basename should include the file
        extension.
    delimiter : str
        Delimiter used to split the basename to exclude
        the file extension.

    Returns
    -------
    basename : str
        File basename extracted from input path.
    """
    basename = os.path.basename(file_path)

    if file_extension is False:
        basename = basename.split(delimiter)[0]

    return basename


def get_locus_id(locus_path):
    """Extract the locus identifier from a path.

    Parameters
    ----------
    locus_path : str
        Path to the locus Fasta file.

    Returns
    -------
    locus_id : str
        Locus identifier without the .fasta extension.
    """
    locus_basename = file_basename(locus_path)
    locus_id = im.match_regex(locus_basename, ct.LOCUS_ID_PATTERN)

    return locus_id


def remove_files(files):
    """Delete a list of files.

    Parameters
    ----------
    files : list
        List with paths to the files to be deleted.
    """
    for f in files:
        if os.path.isfile(f) is True:
            os.remove(f)


def hash_file(file, hash_object, buffer_size=65536):
    """Compute hash based on the contents of a file.

    Parameters
    ----------
    file : str
        Path to a file.
    hash_object : _hashlib.HASH
        Hashlib object to update based on file
        contents.
    buffer_size : int
        Buffer size (amount of data, in KB, to
        read from file).

    Returns
    -------
    hash_str : str
        Hash computed from file contents.
    """
    updated_hash = hash_object

    with open(file, 'rb') as f:
        # read file in chunks
        while True:
            data = f.read(buffer_size)
            if not data:
                break
            updated_hash.update(data)

    hash_str = updated_hash.hexdigest()

    return hash_str


def filter_files(files, suffixes, reverse=False):
    """Filter files names based on a list of suffixes.

    Parameters
    ----------
    files : list
        A list with filenames or file paths.
    suffixes : list
        List with suffixes.
    reverse : bool
        False to select files that contain any of the suffixes.
        True to select files that do not contain any of the
        suffixes.

    Returns
    -------
    filtered : list
        List with files that passed filtering.
    """
    filtered = [file for file in files
                if any([True for suffix in suffixes if suffix in file])]

    if reverse is True:
        filtered = list(set(files)-set(filtered))

    return filtered


def create_directory(directory_path):
    """Create a diretory if it does not exist."""
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        return True
    else:
        return False


def join_paths(parent_path, child_paths):
    """Create path by joining a parent directory and a list of child paths."""
    joined_paths = os.path.join(parent_path, *child_paths)

    return joined_paths


def listdir_fullpath(directory_path, substring_filter=False):
    """Get the full path for all files in a directory.

    Parameters
    ----------
    directory_path : str
        Path to a directory.
    substring_filter : str
        Only list files that contain this substring.

    Returns
    -------
    file_list : list
        List containing the full path for every selected
        file in the input directory.
    """
    file_list = [os.path.join(directory_path, f)
                 for f in os.listdir(directory_path)]

    if substring_filter is not False:
        file_list = filter_files(file_list, [substring_filter])

    return file_list


def delete_directory(directory_path, max_retries=5):
    """Delete a directory."""
    for i in range(max_retries):
        # might fail to delete files due to permission issues
        try:
            shutil.rmtree(directory_path)
        # path does not exist
        # FileNotFoundError is a subclass of OSError
        # the latter will catch the former if checked first
        except FileNotFoundError:
            break
        # sleep before retry
        except OSError:
            time.sleep(i)

    # check if directory still exists
    exists = os.path.isdir(directory_path)
    if exists is True:
        print('Could not remove {0}'.format(directory_path))

    return exists


def read_lines(input_file, strip=True, num_lines=None):
    """Read lines in a file.

    Parameters
    ----------
    input_file : str
        Path to the input file.
    strip : bool
        Specify if lines should be stripped of leading
        and trailing white spaces and new line characters.

    Returns
    -------
    lines : list
        List with the lines read from the input file.
    """
    with open(input_file, 'r') as infile:
        if num_lines is None:
            lines = [line for line in infile.readlines()]
        else:
            lines = list(islice(infile, num_lines))

    if strip is True:
        lines = [line.strip() for line in lines]

    return lines


def pickle_dumper(content, output_file):
    """Use the Pickle module to serialize an object.

    Parameters
    ----------
    content : type
        Variable that refers to the object that will
        be serialized and written to the output file.
    output_file : str
        Path to the output file.
    """
    with open(output_file, 'wb') as poutfile:
        pickle.dump(content, poutfile)


def pickle_loader(input_file):
    """Use the Pickle module to de-serialize an object.

    Parameters
    ----------
    input_file : str
        Path to file with byte stream to be de-serialized.

    Returns
    -------
    content : type
        Variable that refers to the de-serialized
        object.
    """
    with open(input_file, 'rb') as pinfile:
        content = pickle.load(pinfile)

    return content


def file_zipper(input_file, zip_file):
    """Zip (compresses) a file.

    Parameters
    ----------
    input_file : str
        Path to the file that will be compressed.
    zip_file : str
        Path to the ZIP file that will be created.

    Returns
    -------
    zip_file : str
        Path to the ZIP file that was created by
        compressing the input file.
    """
    with zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(input_file, os.path.basename(input_file))

    return zip_file


def unzip_file(compressed_file, archive_type='.gz'):
    """Uncompresse a file.

    Parameters
    ----------
    compressed_file : str
        Path to the compressed file.
    archive_type : str
        Archive format.

    Returns
    -------
    uncompressed_file : str
        Path to the uncompressed file.
    """
    lines = []
    with gzip.open(compressed_file, 'rb') as f:
        for line in f:
            lines.append(line.decode())

    # save uncompressed contents
    uncompressed_file = compressed_file.rstrip('.gz')
    write_lines(lines, uncompressed_file, joiner='')

    return uncompressed_file


def download_file(file_url, outfile, max_tries=3):
    """Download a file.

    Parameters
    ----------
    file_url : str
        URL to download file.
    outfile : str
        Path to the downloaded file.
    max_tries : int
        Maximum number of retries if the download
        fails.
    """
    tries = 0
    downloaded = False
    while downloaded is False and tries <= max_tries:
        try:
            res = urllib.request.urlretrieve(file_url, outfile)
            if os.path.isfile(outfile) is True:
                downloaded = True
        except Exception as e:
            print(e)
            time.sleep(1)
        tries += 1

    return downloaded


def concatenate_files(files, output_file, header=None):
    """Concatenate files.

    Parameters
    ----------
    files : list
        List with the paths to the files to concatenate.
    output_file : str
        Path to the output file that will store the concatenation
        of input files.
    header : str or NoneType
        Specify a header that should be written as the first line
        in the output file.

    Returns
    -------
    output_file : str
        Path to the output file that was created with
        the concatenation of input files.
    """
    with open(output_file, 'w') as outfile:
        if header is not None:
            outfile.write(header)
        for file in files:
            with open(file, 'r') as infile:
                shutil.copyfileobj(infile, outfile)

    return output_file


def write_to_file(text, output_file, write_mode, end_char):
    """Write a single string to a file.

    Parameters
    ----------
    text : str
        A single string to write to the output file.
    output_file : str
        Path to the output file.
    write_mode : str
        Specify write mode ('w' creates file if it does not
        exist and truncates and over-writes existing file,
        'a' creates file if it does not exist and appends to
        the end of file if it exists.).
    end_char : str
        Character added to the end of the file.
    """
    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def write_lines(lines, output_file, joiner='\n', write_mode='w'):
    """Write a list of strings to a file.

    Parameters
    ----------
    lines : list
        List with the lines/strings to write to the output
        file.
    output_file : str
        Path to the output file.
    joiner : str
        Character used to join lines.
    write_mode : str
        Specify write mode ('w' creates file if it does not
        exist and truncates and over-writes existing file,
        'a' creates file if it does not exist and appends to
        the end of file if it exists).
    """
    joined_lines = im.join_list(lines, joiner)

    write_to_file(joined_lines, output_file, write_mode, '\n')


def read_tabular(input_file, delimiter='\t'):
    """Read a tabular (TSV) file.

    Parameters
    ----------
    input_file : str
        Path to a tabular file.
    delimiter : str
        Delimiter used to separate file fields.

    Returns
    -------
    lines : list
        A list with a sublist per line in the input file.
        Each sublist has the fields that were separated by
        the defined delimiter.
    """
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def input_timeout(prompt, timeout=30):
    """Add timeout feature when requesting user input.

    Parameters
    ----------
    prompt : str
        Message to print to stdout to request for user
        input.
    timeout : int
        Maximum number of seconds that the process will
        wait for input.

    Returns
    -------
    String with user input.

    Raises
    ------
    SystemExit
        - If there is no user input before timeout.
    """
    pool = ThreadPool(processes=1)
    answer = pool.apply_async(input, args=[prompt])

    try:
        return answer.get(timeout=timeout)
    except TimeoutError:
        sys.exit('Timed out.')


def is_file_empty(file_path):
    """Check if file is empty by confirming if its size is 0 bytes."""
    # Check if file exist and it is empty
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0


def is_file_empty_2(file_name):
    """Check if file is empty by confirming if its size is 0 bytes."""
    # Check if file exist and it is empty
    return os.path.isfile(file_name) and os.path.getsize(file_name) == 0


def is_file_empty_3(file_name):
    """Check if file is empty by reading first character in it."""
    # open ile in read mode
    with open(file_name, 'r') as read_obj:
        # read first character
        one_char = read_obj.read(1)
        # if not fetched then file is empty
        if not one_char:
            return True
        elif one_char == " ":
            return True
    return False


def create_short(schema_files, schema_dir):
    """Create the 'short' directory for a chewBBACA schema.

    Creates the 'short' directory inside a schema's main
    directory and copies the FASTA files with the representative
    alleles into the 'short' folder.

    Parameters
    ----------
    schema_files : list
        List with paths to FASTA files that contain the loci
        representative alleles.
    schema_dir : str
        Path to the schema's directory.

    Returns
    -------
    True on completion.
    """
    short_path = join_paths(schema_dir, ['short'])
    create_directory(short_path)

    for file in schema_files:
        short_file = join_paths(short_path, [file_basename(file)])
        short_file = short_file.replace('.fasta', '_short.fasta')
        shutil.copy(file, short_file)

    return True


def move_file(source, destination):
    """Move a file to specified destination."""
    shutil.move(source, destination)


def matching_lines(input_file, pattern):
    """Retrieve lines from a file that contain a pattern.

    Parameters
    ----------
    input_file : str
        Path to input file.
    pattern : str
        Pattern to match.

    Returns
    -------
    matched_lines : list
        List with lines that contain the pattern.
    """
    with open(input_file, 'r') as infile:
        matched_lines = [line for line in infile if pattern in line]

    return matched_lines


def transpose_matrix(input_file, output_directory):
    """Transpose a TSV file.

    Parameters
    ----------
    input_file : str
        Path to the input TSV file.
    output_directory : str
        Path to the directory to which intermediate files
        and the complete transposed file will be written.

    Returns
    -------
    transposed_file : str
        Path to the file with the transposed matrix.
        This file is created by concatenating all
        intermediate files.
    """
    intermediate_files = []
    with open(input_file, 'r') as infile:
        # get column identifiers
        columns = [c.strip() for c in (infile.__next__()).split('\t')]
        # divide into smaller sets to avoid loading complete file
        total_column_sets = math.ceil(len(columns)/500)
        column_sets = im.divide_list_into_n_chunks(columns, total_column_sets)
        # use Pandas to read columns sets and save transpose
        for i, c in enumerate(column_sets):
            # dtype=str or Pandas converts values into floats
            df = pd.read_csv(input_file, usecols=c, delimiter='\t', dtype=str)
            output_file = join_paths(output_directory, ['chunk{0}.tsv'.format(i)])
            # transpose columns
            df = df.T
            # do not save header that contains row indexes
            df.to_csv(output_file, sep='\t', header=False)
            intermediate_files.append(output_file)

    # concatenate all files with transposed lines
    transposed_file = input_file.replace('.tsv', '_transpose.tsv')
    concatenate_files(intermediate_files, transposed_file)

    # delete intermediate files
    remove_files(intermediate_files)

    return transposed_file
