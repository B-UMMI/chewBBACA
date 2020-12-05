#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


from itertools import islice


def merge_dictionaries(dictionaries_list):
    """ Merges several dictionaries into a single dictionary.

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

    merged_dicts = {}
    for d in dictionaries_list:
        merged_dicts = {**merged_dicts, **d}

    return merged_dicts


def invert_dictionary(dictionary):
    """ Inverts a dictionary (Keys become values and vice-versa).

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
        chunks.append({k: iterable[k] for k in islice(it, size)})

    return chunks


def select_clusters(clusters, cluster_size):
    """ Determines clusters that contain a specified number
        of sequences.

        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are clusters representatives as keys and
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
