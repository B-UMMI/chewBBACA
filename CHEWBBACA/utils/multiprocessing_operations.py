#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions used to paralellize
function calls.

Code documentation
------------------
"""


import time
import traceback
import multiprocessing.pool

from tqdm import tqdm
from loguru import logger

from CHEWBBACA.utils import iterables_manipulation as im


def function_helper(input_args):
	"""Run function with provided inputs and capture exceptions.

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
	except Exception as err:
		func_name = (input_args[-1]).__name__
		traceback_lines = traceback.format_exc()
		traceback_text = "".join(traceback_lines)
		logger.error("\nError on {func_name}:\n{traceback_text}\n")
		results = [func_name, traceback_text]

	return results


def map_async_parallelizer(inputs, function, cpu, chunksize=1, pool_type="pool"):
	"""Run function in parallel.

	Parameters
	----------
	inputs : list
		List with inputs to process.
	function : func
		Function to be parallelized.
	cpu : int
		Number of processes to create (based on the
		number of CPU cores).
	chunksize : int
		Size of input chunks that will be passed to
		each process. The function will create groups
		of inputs with this number of elements.
	pool_type : str
		The multiprocessing.pool object that will be used,
		Pool or ThreadPool.

	Returns
	-------
	results : list
		List with the results returned for each function
		call.
	"""
	if pool_type == "pool":
		multiprocessing_function = multiprocessing.pool.Pool
	# Gene prediction uses ThreadPool because Pyrodigal might hang with Pool
	elif pool_type == "threadpool":
		multiprocessing_function = multiprocessing.pool.ThreadPool

	# Use context manager to join and close pool automatically
	with multiprocessing_function(cpu) as pool:
		# Use tqdm to show progress bar
		# imap_unordered returns results as processes complete, not in the same order as inputs
		# Gather results into a list (e.g., if the function being executed in parallel returns a dict, it will be a list of dict)
		results = list(tqdm(pool.imap_unordered(function, inputs, chunksize=chunksize), total=len(inputs), desc="Processing"))

	logger.info(f"Finished running {function} in parallel. Processed {len(inputs)} items.")

	return results


def distribute_loci(inputs, cores, method):
	"""Create balanced lists of loci to efficiently parallelize function calls.

	Creates balanced lists of loci to distribute per number of
	available cores. Loci lists can be created based on the number
	of sequences per locus (seqcount), the mean length of the
	sequences (length) in each locus or the product of both values
	(seqcount+length).

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
			splitted_values[i] += locus[4]
		elif method == 'seqcount+length':
			splitted_values[i] += locus[1] * locus[4]
		splitted_ids[i].append(locus[0])
		# at the end of each iteration, choose the sublist
		# with lowest criterion value
		i = splitted_values.index(min(splitted_values))

	return splitted_ids


def parallelize_function(function, inputs, common_args=None, cpu_cores=1):
	"""Create list of inputs and parallelize function calls.

	Parameters
	----------
	function : func
		Function to be parallelized.
	inputs : list
		List of inputs to divide into sublists.
	common_args : list
		List of common arguments to add to each sublist.
		The common args are values that will be passed to the
		function.
	cpu_cores : int
		Number of CPU cores used to parallelize the function.

	Returns
	-------
	results : list
		List with the results returned by the function calls.
	"""
	# Create chunks to distribute per cores
	input_lists = im.divide_list_into_n_chunks(inputs, len(inputs))

	if common_args is None:
		common_args = []

	# Add common arguments to all sublists
	input_lists = im.multiprocessing_inputs(input_lists,
											common_args,
											function)

	results = map_async_parallelizer(input_lists,
									 function_helper,
									 cpu_cores)

	return results
