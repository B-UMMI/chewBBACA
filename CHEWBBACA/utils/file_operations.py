#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to file operations,
such as read and write files, create and delete files,
manipulate file paths, compress files, verify file contents,
etc.

Code documentation
------------------
"""


import os
import sys
import csv
import shutil
import pickle
import hashlib
import zipfile
from multiprocessing import TimeoutError
from multiprocessing.pool import ThreadPool

try:
    from utils import (iterables_manipulation as im,
                       file_operations as fo)
except:
    from CHEWBBACA.utils import (iterables_manipulation as im,
                                 file_operations as fo)


def file_basename(file_path, suffix=True):
    """ Extract file basename from path.

        Parameters
        ----------
        file_path : str
            Path to the file.
        suffix : bool
            Specify if the basename should include the file
            extension.

        Returns
        -------
        basename : str
            File basename extracted from input path.
    """

    basename = os.path.basename(file_path)

    if suffix is False:
        basename = basename.split('.')[0]

    return basename


def remove_files(files):
    """ Deletes a list of files.

        Parameters
        ----------
        files : list
            List with paths to the files to be deleted.
    """

    for f in files:
        os.remove(f)


def hash_file(file, read_mode):
    """ Computes BLAKE2b hash based on the contents of a file.

        Parameters
        ----------
        file : str
            Path to a file.
        read_mode : str
            File read mode.

        Returns
        -------
        hash_str : str
            Hash computed from file contents.
    """

    with open(file, read_mode) as f:
        hash_obj = hashlib.blake2b()
        file_content = f.read()
        hash_obj.update(file_content)
        hash_str = hash_obj.hexdigest()

    return hash_str


def filter_files(files, suffixes, reverse=False):
    """ Filters files names based on a list of suffixes.

        Parameters
        ----------
        files : list
            A list with filenames of file paths.
        suffixes : list
            List with suffixes.
        reverse : bool
            True if files should be filtered out from
            input list. False to filter out files without
            any of the suffixes.

        Returns
        -------
        filtered : list
            List with files that passed filtering.
    """

    if reverse is False:
        filtered = [file for file in files
                    if any([True for suffix in suffixes if suffix in file])]
    elif reverse is True:
        filtered = [file for file in files
                    if not any([True for suffix in suffixes if suffix in file])]

    return filtered


def create_directory(directory_path):
    """ Creates a diretory if it does not exist."""

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)


def join_paths(parent_path, child_paths):
    """ Creates a new path by joining a parent directory
        and a list with child paths."""

    joined_paths = os.path.join(parent_path, *child_paths)

    return joined_paths


def listdir_fullpath(directory_path):
    """ Gets the full path for all files in a directory.

        Parameters
        ----------
        directory_path : str
            Path to a directory.

        Returns
        -------
        List containing the full path for every file
        in the input directory.
    """

    return [os.path.join(directory_path, f)
            for f in os.listdir(directory_path)]


def delete_directory(directory_path):
    """ Deletes a directory. """

    shutil.rmtree(directory_path)


def read_lines(input_file, strip=True):
    """ Reads lines in an input file and stores those lines
        in a list.

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
        if strip is True:
            lines = [file.strip() for file in infile.readlines()]
        else:
            lines = [file for file in infile.readlines()]

    return lines


def pickle_dumper(content, output_file):
    """ Use the Pickle module to serialize an object.

        Parameters
        ----------
        content : type
            Variable that refers to the object that will
            be serialized and written to the output file.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'wb') as po:
        pickle.dump(content, po)


def pickle_loader(input_file):
    """ Use the Pickle module to de-serialize an object.

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

    with open(input_file, 'rb') as pi:
        content = pickle.load(pi)

    return content


def file_zipper(input_file, zip_file):
    """ Zips (compresses) a file.

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


def concatenate_files(files, output_file, header=None):
    """ Concatenates the contents of a set of files.

        Parameters
        ----------
        files : list
            List with the paths to the files to concatenate.
        output_file : str
            Path to the output file that will store the
            concatenation of input files.
        header : str or NoneType
            Specify a header that should be written as the
            first line in the output file.

        Returns
        -------
        output_file : str
            Path to the output file that was created with
            the concatenation of input files.
    """

    with open(output_file, 'w') as of:
        if header is not None:
            of.write(header)
        for f in files:
            with open(f, 'r') as fd:
                shutil.copyfileobj(fd, of)

    return output_file


def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single string to a file.

        Parameters
        ----------
        text : str
            A single string to write to the output file.
        output_file : str
            Path to the output file.
        write_mode : str
            Write mode can be 'w', writes text and overwrites
            any text in file, or 'a', appends text to text
            already in file.
        end_char : str
            Character added to the end of the file.
    """

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


def write_lines(lines, output_file):
    """ Writes a list of strings to a file. The strings
        are joined with newlines before being written to
        file.

        Parameters
        ----------
        lines : list
            List with the lines/strings to write to the
            output file.
        output_file : str
            Path to the output file.
    """

    joined_lines = im.join_list(lines, '\n')

    write_to_file(joined_lines, output_file, 'a', '\n')


def read_tabular(input_file, delimiter='\t'):
    """ Read tabular file.

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


def write_list(lines, output_file):
    """ Writes list elements to file.

        Parameters
        ----------
        lines : list
            List with the lines that will be written to
            the output file.
        output_file : str
            Path to the output file.

        Returns
        -------
        Writes contents of input list to the output file
        (function does not add any character between lines).
    """

    with open(output_file, 'w') as file:
        file.writelines(lines)


def input_timeout(prompt, timeout=30):
    """ Adds timeout feature when requesting user input.

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
    except TimeoutError as e:
        sys.exit('Timed out.')


def is_file_empty(file_path):
    """ Check if file is empty by confirming if its size is 0 bytes. """

    # Check if file exist and it is empty
    return os.path.exists(file_path) and os.stat(file_path).st_size == 0


def is_file_empty_2(file_name):
    """ Check if file is empty by confirming if its size is 0 bytes. """

    # Check if file exist and it is empty
    return os.path.isfile(file_name) and os.path.getsize(file_name) == 0


def is_file_empty_3(file_name):
    """ Check if file is empty by reading first character in it. """

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
    """ Creates the 'short' directory for a schema.
        Creates the directory and copies schema files
        to the directory (should be used when the schema
        only has 1 sequence per gene/locus).

        Parameters
        ----------
        schema_files : list
            List with paths to all FASTA files in the schema.
        schema_dir : str
            Path to the schema's directory.

        Returns
        -------
        True on completion.
    """

    short_path = join_paths(schema_dir, ['short'])
    create_directory(short_path)

    for file in schema_files:
        short_file = join_paths(short_path, [fo.file_basename(file)])
        short_file = short_file.replace('.fasta', '_short.fasta')
        shutil.copy(file, short_file)

    return True
