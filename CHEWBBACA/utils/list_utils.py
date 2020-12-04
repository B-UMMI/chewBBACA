#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


def join_list(lst, link):
    """ Joins all elements in a list into a single string.

        Parameters
        ----------
        lst : list
            List with elements to be joined.
        link : str
            Character used to join list elements.

        Returns
        -------
        joined_list : str
            A single string with all elements in the input
            list joined by the character chosen as link.
    """

    joined_list = link.join(lst)

    return joined_list


def concatenate_list(str_list, join_char):
    """ Concatenates list elements with specified
        character between each original list element.

        Parameters
        ----------
        str_list : list
            List with strings that will be concatenated.
        join_char : str
            Character that will be used to join list
            elements.

        Returns
        -------
        ids_str : str
            String resulting from the concatenation of
            all strings in the input list.
    """

    concat = join_char.join(str_list)

    return concat


def extend_list(input_list, *elements):
    """
    """

    input_list.extend(elements)

    return input_list


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
