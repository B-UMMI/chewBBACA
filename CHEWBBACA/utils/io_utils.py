#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


import sys
import csv
import shutil
import pickle
import zipfile
from multiprocessing import TimeoutError
from multiprocessing.pool import ThreadPool

try:
    from utils import list_utils as lu
except:
    from CHEWBBACA.utils import list_utils as lu


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

    joined_lines = lu.join_list(lines, '\n')

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


def input_timeout(prompt, timeout):
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
