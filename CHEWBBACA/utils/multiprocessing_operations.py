#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This modules contains functions used to paralellize
function calls.

Code documentation
------------------
"""


import time
import traceback
from multiprocessing import Pool


def function_helper(input_args):
    """ Runs function by passing set of provided inputs and
        captures exceptions raised during function execution.

        Parameters
        ----------
        input_args : list
            List with function inputs and function object to call
            in the last index.

        Returns
        -------
        results : list
            List with the results returned by the function.
            If an exception is raised it returns a list with
            the name of the function and the exception traceback.
    """

    try:
        results = input_args[-1](*input_args[0:-1])
    except Exception as e:
        func_name = (input_args[-1]).__name__
        traceback_lines = traceback.format_exception(etype=type(e), value=e,
                                                     tb=e.__traceback__)
        traceback_text = ''.join(traceback_lines)
        print('Error on {0}:\n{1}\n'.format(func_name, traceback_text))
        results = [func_name, traceback_text]

    return results


def map_async_parallelizer(inputs, function, cpu, callback='extend',
                           chunksize=1, show_progress=False):
    """ Parallelizes function calls by creating several processes
        and distributing inputs.

        Parameters
        ----------
        inputs : list
            List with inputs to process.
        function
            Function to be parallelized.
        cpu : int
            Number of processes to create (based on the
            number of cores).
        callback : str
            Results can be appended, 'append', to the
            list that stores results or the list of results
            can be extended, 'extend'.
        chunksize : int
            Size of input chunks that will be passed to
            each process. The function will create groups
            of inputs with this number of elements.
        show_progress: bool
            True to show a progress bar with the percentage
            of inputs that have been processed, False
            otherwise.

        Returns
        -------
        results : list
            List with the results returned for each function
            call.
    """

    results = []
    pool = Pool(cpu)
    if callback == 'extend':
        rawr = pool.map_async(function, inputs,
                              callback=results.extend, chunksize=chunksize)
    elif callback == 'append':
        rawr = pool.map_async(function, inputs,
                              callback=results.append, chunksize=chunksize)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = progress_bar(rawr, len(inputs))

    rawr.wait()

    return results


def progress_bar(process, total, tickval=5, ticknum=20, completed=False):
    """ Creates and prints progress bar to stdout.

        Parameters
        ----------
        process : multiprocessing.pool.MapResult
            Multiprocessing object.
        total : int
            Total number of inputs that have to be processed.
        tickval : int
            Progress completion percentage value for each
            tick.
        ticknum : int
            Total number of ticks in progress bar.
        completed : bool
            Boolean indicating if process has completed.

        Returns
        -------
        completed : bool
            Boolean indicating if process has completed.
    """

    # check if process has finished
    if (process.ready()):
        # print full progress bar and satisfy stopping condition
        progress_bar = '[{0}] 100%'.format('='*ticknum)
        completed = True

    # check how many inputs have been processed
    remaining = process._number_left
    if remaining == total:
        # print empty progress bar
        progress_bar = '[{0}] 0%'.format(' '*ticknum)
    else:
        # print progress bar, incremented by 5%
        progress = int(100-(remaining/total)*100)
        progress_tick = progress//tickval
        progress_bar = '[{0}{1}] {2}%'.format('='*progress_tick,
                                              ' '*(ticknum-progress_tick),
                                              progress)

    print('\r', progress_bar, end='')
    time.sleep(0.5)

    return completed


def split_genes_by_core(inputs, cores, method):
    """ Creates balanced lists of loci to distribute per number
        of available cores. Loci lists can be created based
        on the number of sequence per locus (seqcount), the mean
        length of the sequences in each locus or the product of
        both values.

        Parameters
        ----------
        inputs : list
            List with one sublist per locus. Each sublist has
            a locus identifier, the total number of sequences
            and sequence mean legth for that locus.
        cores : int
            The number of loci groups that should be created.
            Based on the number of CPU cores that will be
            used to process the inputs.
        method : str
            "seqcount" to create loci lists based on the total
            number of sequences, "length" to split based
            on mean length of sequences and "seqcount+length" to
            split based on both criteria.

        Returns
        -------
        splitted_ids : list
            List with sublists that contain loci identifiers.
            Sublists are balanced based on the chosen method.
    """

    # initialize list with sublists to store inputs
    splitted_ids = [[] for cpu in range(cores)]
    # initialize list with chosen criterion values
    # for each sublist of inputs
    splitted_values = [0 for cpu in range(cores)]
    i = 0
    for locus in inputs:
        if method == 'seqcount':
            splitted_values[i] += locus[1]
        elif method == 'length':
            splitted_values[i] += locus[2]
        elif method == 'seqcount+length':
            splitted_values[i] += locus[1] * locus[2]
        splitted_ids[i].append(locus[0])
        # at the end of each iteration, choose the sublist
        # with lowest criterion value
        i = splitted_values.index(min(splitted_values))

    return splitted_ids
