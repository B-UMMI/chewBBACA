#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


import os
import shutil
import hashlib


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
    """ Creates a hash based on the contents of a file.

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
            input list or False otherwise.

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
    """ Gets the full path of the files from a directory.

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
    """ Deletes a directory.
    """

    shutil.rmtree(directory_path)
