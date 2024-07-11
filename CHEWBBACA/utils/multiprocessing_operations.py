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

try:
	from utils import iterables_manipulation as im
except ModuleNotFoundError:
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
	except Exception as e:
		func_name = (input_args[-1]).__name__
		traceback_lines = traceback.format_exc()
		traceback_text = ''.join(traceback_lines)
		print('\nError on {0}:\n{1}\n'.format(func_name, traceback_text), flush=True)
		results = [func_name, traceback_text]

	return results


def map_async_parallelizer(inputs, function, cpu, callback='extend',
						   chunksize=1, show_progress=False, pool_type='pool'):
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
	callback : str
		Results can be appended, "append", to the
		list that stores results or the list of results
		can be extended, "extend".
	chunksize : int
		Size of input chunks that will be passed to
		each process. The function will create groups
		of inputs with this number of elements.
	show_progress : bool
		True to show a progress bar with the percentage
		of inputs that have been processed, False
		otherwise.
	pool_type : str
		The multiprocessing.pool object that will be used,
		Pool or ThreadPool.

	Returns
	-------
	results : list
		List with the results returned for each function
		call.
	"""
	if pool_type == 'pool':
		multiprocessing_function = multiprocessing.pool.Pool
	# Gene prediction uses ThreadPool because Pyrodigal might hang with Pool
	elif pool_type == 'threadpool':
		multiprocessing_function = multiprocessing.pool.ThreadPool

	results = []
	# Use context manager to join and close pool automatically
	with multiprocessing_function(cpu) as pool:
		if callback == 'extend':
			rawr = pool.map_async(function, inputs,
								  callback=results.extend, chunksize=chunksize)
		elif callback == 'append':
			rawr = pool.map_async(function, inputs,
								  callback=results.append, chunksize=chunksize)

		if show_progress is True:
			progress = None
			while progress != 100:
				progress = progress_bar(rawr._number_left, len(inputs), progress)

		rawr.wait()

	return results


def progress_bar(remaining, total, previous, tickval=5, ticknum=20):
	"""Create and print a progress bar to the stdout.

	Parameters
	----------
	remaining : int
		Number of remaining tasks to complete.
	total : int
		Total number of inputs that have to be processed.
	previous : int
		Percentage of tasks that had been completed in the
		previous function call.
	tickval : int
		Progress completion percentage value for each
		tick.
	ticknum : int
		Total number of ticks in progress bar.

	Returns
	-------
	completed : bool
		Boolean indicating if all inputs have been processed.
	"""
	# determine percentage of processed inputs
	progress = int(100-(remaining/total)*100)
	# only print if percentage has changed
	if progress != previous:
		progress_tick = progress//tickval
		progress_bar = '[{0}{1}] {2}%'.format('='*progress_tick,
											  ' '*(ticknum-progress_tick),
											  progress)
		print('\r', progress_bar, end='')

	time.sleep(0.1)

	return progress


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


def parallelize_function(function, inputs, common_args=None,
						 cpu_cores=1, show_progress=False):
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
	show_progress : bool
		True to show a progress bar for the percentage of
		inputs that have been processed, False otherwise.

	Returns
	-------
	results : list
		List with the results returned by the function calls.
	"""
	# create chunks to distribute per cores
	input_lists = im.divide_list_into_n_chunks(inputs, len(inputs))

	if common_args is None:
		common_args = []

	# add common arguments to all sublists
	input_lists = im.multiprocessing_inputs(input_lists,
											common_args,
											function)

	results = map_async_parallelizer(input_lists,
									 function_helper,
									 cpu_cores,
									 show_progress=show_progress)

	return results
