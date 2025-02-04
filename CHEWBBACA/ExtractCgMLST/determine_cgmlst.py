#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module determines the set of loci that constitute the core genome
based on a matrix with allelic profiles and a loci presence threshold.

Code documentation
------------------
"""


import os
import numpy as np
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   iterables_manipulation as im)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 file_operations as fo,
								 iterables_manipulation as im)


def binarize_matrix(column):
	"""Convert non-zero numeric values in a Pandas Series to 1.

	Parameters
	----------
	column : pandas.core.series.Series
		Pandas Series.

	Returns
	-------
	Numpy array corresponding to the input Pandas
	Series with non-zero numeric values converted to 1.
	"""
	numeric_column = pd.to_numeric(column)

	return np.int64(numeric_column > 0)


def remove_genomes(matrix, genomesToRemove):
	"""Remove rows from an allele calling matrix.

	Remove rows from a Pandas dataframe if the
	index identifier matches the identifier of
	a genome to remove.

	Parameters
	----------
	matrix : pandas.core.frame.DataFrame
		Pandas dataframe with allelic profiles.
		Each row has the allelic profile of a genome
		and each column has the allele identifiers
		determined for a locus.
	genomesToRemove : list
		List of genomes to remove.

	Returns
	-------
	pruned_matrix : pandas.core.frame.DataFrame
		Input dataframe without the rows whose
		index matched an identifier of a genome
		to remove.
	"""
	# determine row indexes that match any genome to remove
	to_remove_bool = matrix.index.isin(genomesToRemove)
	# create new matrix without rows that matched any genome to remove
	pruned_matrix = matrix.loc[~ to_remove_bool]

	return pruned_matrix


def remove_columns(matrix, columnsToRemove):
	"""Remove columns from a matrix.

	Parameters
	----------
	matrix : pandas.core.frame.DataFrame
		Pandas dataframe.
	columnsToRemove : list
		List of columns to remove. Must match
		the column names in the dataframe.
	"""
	matrix = matrix[matrix.columns[~matrix.columns.isin(columnsToRemove)]]

	return matrix


def above_threshold(column, column_length, threshold):
	"""Determine if the sum of a column values is equal or above a presence/absence threshold.

	Parameters
	----------
	column : pandas.core.series.Series
		Pandas Series with presence (1) and absence (0) values.
	column_length : int
		Total number of values.
	threshold : float
		Presence/absence threshold value.

	Returns
	-------
	bool
		True if the sum of the column values is equal or above the
		presence/absence threshold, False otherwise.
	"""
	return (np.sum(column) / column_length) >= threshold


def compute_cgMLST(matrix, sorted_genomes, threshold, step, compute_accessory):
	"""Compute the core genome based on loci presence/absence.

	Parameters
	----------
	matrix : pandas.core.frame.DataFrame
		Pandas dataframe with allelic profiles.
		Each row has the allelic profile of a genome
		and each column has the allele identifiers
		determined for a single gene.
	sorted_genomes : list
		List of genome identifiers sorted in order of
		decresing number of missing loci.
	threshold : float
		Loci presence/absence threshold to determine the core genome.
	step : int
		Number of genomes added to the core genome computation at
		each step.
	compute_accessory : bool
		Determine the accessory genome.

	Returns
	-------
	pruned_df : pandas.core.frame.DataFrame
		Dataframe with the cgMLST profiles for the last step value
		(the last step includes all genomes in the input matrix).
	cgMLST_size : dict
		Dictionary with the number of genomes used to compute the
		cgMLST as keys and the size of the core-genome as values.
	"""
	print('Determining core genome for locus presence/absence threshold of {0}...'.format(threshold))
	# Determine loci at or above threshold
	above = None
	cgMLST_size = {}
	below = None
	agMLST_size = {}
	for i in im.inclusive_range(1, len(sorted_genomes), step):
		# Get subdataframe for current genomes
		current_df = matrix.loc[sorted_genomes[:i]]
		pa_rows, _ = current_df.shape
		is_above_threshold = current_df.apply(above_threshold,
											  args=(pa_rows, threshold,))
		above = current_df.columns[is_above_threshold]
		cgMLST_size[pa_rows] = len(above)
		if compute_accessory:
			# Compute accessory genome
			below = current_df.columns[~is_above_threshold]
			agMLST_size[pa_rows] = len(below)
		print('\r', 'Computed for...{0} genomes.'.format(i), end='')

	# Return list of genes in cgMLST and cgMLST count per step iteration
	return [above, cgMLST_size, below, agMLST_size]


def compute_presence_absence(matrix, output_directory):
	"""Compute a presence-absence matrix.

	Parameters
	----------
	matrix : pandas.core.frame.DataFrame
		Pandas dataframe where zero values indicate absence.
	output_directory : str
		Path to the directory where the TSV file with
		the presence absence matrix will be stored.

	Returns
	-------
	presence_absence : pandas.core.frame.DataFrame
		Pandas dataframe with all non-zero values converted to 1.
	pa_outpath : str
		Path to the output TSV file that contains the
		presence-absence matrix.
	"""
	presence_absence = matrix.apply(binarize_matrix)

	pa_outpath = fo.join_paths(output_directory, [ct.PRESENCE_ABSENCE_BASENAME])
	presence_absence.to_csv(pa_outpath, sep='\t')

	return [presence_absence, pa_outpath]


def count_zeros(matrix, axis=1, column_names=None):
	"""Count zeros per dataframe row.

	Parameters
	----------
	matrix : pandas.core.frame.DataFrame
		Pandas dataframe with numeric values. Zero values
		correspond to missing data.

	Returns
	-------
	zeros_df : pandas.core.frame.DataFrame
		Dataframe with the number and percentage of zeros per
		row in the input matrix.
	"""
	nrows, ncols = matrix.shape
	non_zeros = matrix.apply(np.count_nonzero, axis=axis)

	zeros_data = [matrix.index if axis == 1 else matrix.columns,
				  (ncols-non_zeros) if axis == 1 else (nrows-non_zeros),
				  (1-(non_zeros/ncols))*100 if axis == 1 else (1-(non_zeros/nrows))*100]

	zeros_df = pd.DataFrame(list(zip(*zeros_data)),
							columns=column_names)

	return zeros_df


def main(input_file, output_directory, threshold, step,
		 compute_accessory, exclude_loci, exclude_genomes):
	"""Determine the cgMLST based on allele calling results.

	Parameters
	----------
	input_file : str
		Path to a TSV file containing allelic profiles.
	output_directory : str
		Path to the directory where the process will
		store output files.
	threshold : list
		Loci presence threshold values used to determine the core genome.
	step : int
		Number of genomes added to the core genome computation at
		each step.
	compute_accessory : bool
		Compute the results for the accessory genome. The accessory
		genome corresponds to all the loci not included in the core
		genome. The accessory genome is determined for each core
		genome threshold.
	exclude_loci : str
		Path to TXT file with a list of loci to exclude from
		the analysis.
	exclude_genomes : str
		Path to TXT file with a list of genomes to exclude from
		the analysis.
	"""
	fo.create_directory(output_directory)

	# Import allelic profiles
	profiles = pd.read_csv(input_file, header=0, index_col=0,
						 sep='\t', low_memory=False)

	# Get number of genomes and loci
	total_genomes, total_loci = profiles.shape
	print('Input file has {0} profiles for {1} '
		  'loci.'.format(total_genomes, total_loci))

	cgMLST_thresholds = sorted(threshold)
	print('Core genome thresholds: {0}'.format(', '.join(map(str, cgMLST_thresholds))))
	
	# Read lists of loci and genomes to exclude from the analysis
	genomes_to_remove = []
	if exclude_genomes:
		genomes_to_remove = fo.read_lines(exclude_genomes)
	print('{0} genomes to exclude.'.format(len(genomes_to_remove)))

	# Remove genomes
	if len(genomes_to_remove) > 0:
		profiles = remove_genomes(profiles, genomes_to_remove)
		print('Excluded {0} genomes.'.format(total_genomes - profiles.shape[0]))

	loci_to_remove = []
	if exclude_loci:
		loci_to_remove = fo.read_lines(exclude_loci)
	print('{0} loci to exclude.'.format(len(loci_to_remove)))

	# Remove loci
	if len(loci_to_remove) > 0:
		profiles = remove_columns(profiles, loci_to_remove)
		print('Excluded {0} loci.'.format(total_loci - profiles.shape[1]))

	total_genomes, total_loci = profiles.shape
	print('Processed file has {0} profiles for {1} loci.'.format(total_genomes, total_loci))

	# Mask special classifications and remove 'INF-' prefixes
	print('Masking profiles...')
	masked_profiles = profiles.apply(im.replace_chars)
	print('Masked {0} profiles.'.format(total_genomes))

	# Compute presence-absence matrix
	print('Computing presence-absence matrix...')
	pa_matrix, pa_outfile = compute_presence_absence(masked_profiles, output_directory)
	print('Presence-absence matrix saved to {0}'.format(pa_outfile))

	# Count number of special classifications per genome
	print('Computing missing data per genome...')
	genome_mdata_df = count_zeros(pa_matrix, column_names=ct.GENOMES_MISSING_COLUMNS)
	genome_mdata_stats = genome_mdata_df['missing'].describe()

	# Sort based on decreasing number of special classifications
	genome_mdata_df = genome_mdata_df.sort_values('missing', ascending=True)
	sorted_genomes = genome_mdata_df['FILE'].tolist()

	# Write TSV with special classifications statistics
	genome_mdata_path = os.path.join(output_directory, ct.GENOMES_MISSING_BASENAME)
	genome_mdata_df.to_csv(genome_mdata_path, sep='\t', index=False)
	print('Missing data per genome saved to {0}'.format(genome_mdata_path))

	# Count number of special classifications per locus
	print('Computing missing data per locus...')
	loci_mdata_df = count_zeros(pa_matrix, axis=0, column_names=ct.LOCI_MISSING_COLUMNS)
	loci_mdata_stats = loci_mdata_df['missing'].describe()

	# Sort based on decreasing number of special classifications
	loci_mdata_df = loci_mdata_df.sort_values('missing', ascending=True)
	sorted_loci = loci_mdata_df.index.tolist()

	# Write TSV with special classifications statistics
	loci_mdata_path = os.path.join(output_directory, ct.LOCI_MISSING_BASENAME)
	loci_mdata_df.to_csv(loci_mdata_path, sep='\t', index=False)
	print('Missing data per locus saved to {0}'.format(loci_mdata_path))

	trace_colors = ['#134130', '#4C825D', '#8CAE9E', '#8DC7DC',
					'#508CA7', '#1A5270', '#0E2A4D']
	# Compute the core genome for each threshold
	cgMLST_traces = []
	agMLST_traces = []
	for i, t in enumerate(cgMLST_thresholds):
		current_color = trace_colors[i] if i < len(trace_colors) else trace_colors[i%len(trace_colors)]
		cgMLST_results = compute_cgMLST(pa_matrix, sorted_genomes,
										t, step, compute_accessory)

		cgMLST_loci, cgMLST_counts, agMLST_loci, agMLST_counts = cgMLST_results

		print('\nCore genome for locus presence threshold of {0} composed '
			  'of {1}/{2} loci.'.format(t, len(cgMLST_loci), total_loci))

		# Write cgMLST matrix
		# Get subset from masked matrix
		cgMLST_matrix = masked_profiles[cgMLST_loci]
		cgMLST_path = os.path.join(output_directory, 'cgMLST{0}.tsv'.format(int(t*100)))
		cgMLST_matrix.to_csv(cgMLST_path, sep='\t')
		print('cgMLST profiles for threshold {0} saved to {1}'.format(t, cgMLST_path))

		# Write list of cgMLST loci
		cgMLST_loci_path = os.path.join(output_directory, 'cgMLSTschema{0}.txt'.format(int(t*100)))
		fo.write_lines(list(cgMLST_loci), cgMLST_loci_path)
		print('List of core loci for threshold {0} saved to {1}'.format(t, cgMLST_loci_path))

		# Create line plot for core genome values at each step value
		trace = go.Scattergl(x=list(cgMLST_counts.keys()),
							 y=list(cgMLST_counts.values()),
							 mode='lines',
							 line=dict(color=current_color),
							 name='cgMLST{0}'.format(int(t*100)),
							 hovertemplate=('%{y}'))
		cgMLST_traces.append(trace)

		if compute_accessory:
			# Write agMLST matrix
			# Get subset from masked matrix
			agMLST_matrix = masked_profiles[agMLST_loci]
			agMLST_path = os.path.join(output_directory, 'agMLST{0}.tsv'.format(int(t*100)))
			agMLST_matrix.to_csv(agMLST_path, sep='\t')
			print('agMLST profiles for threshold {0} saved to {1}'.format(t, agMLST_path))

			# Write list of agMLST loci
			agMLST_loci_path = os.path.join(output_directory, 'agMLSTschema{0}.txt'.format(int(t*100)))
			fo.write_lines(list(agMLST_loci), agMLST_loci_path)
			print('List of accessory loci for threshold {0} saved to {1}'.format(t, agMLST_loci_path))

			# Create line plot for accessory genome values at each step value
			trace = go.Scattergl(x=list(agMLST_counts.keys()),
								 y=list(agMLST_counts.values()),
								 mode='lines',
								 line=dict(dash='dash', color=current_color),
								 name='agMLST{0}'.format(int(t*100)),
								 hovertemplate=('%{y}'))
			agMLST_traces.append(trace)

	# Create trace with present loci
	present = [total_loci-i for i in genome_mdata_df['missing']]
	genomes_index = list(range(1, len(present)+1))
	miss_trace = go.Scattergl(x=genomes_index,
							  y=present,
							  mode='lines',
							  name='Present in genome',
							  line=dict(dash='dot', color='#000000'),
							  hovertemplate=('%{y}<br>'
											 'Genome: %{text}<br>'),
							  text=sorted_genomes)

	traces = cgMLST_traces + agMLST_traces + [miss_trace]

	fig = go.Figure(data=traces)
	fig.update_layout(title={'text': 'Number of loci in core genome',
							 'font_size': 30},
					  xaxis_title='Number of genomes',
					  yaxis_title='Number of loci',
					  template='simple_white',
					  hovermode='x')
	fig.update_xaxes(range=[0, len(sorted_genomes)],
					 tickfont=dict(size=18),
					 titlefont=dict(size=20),
					 showgrid=True)
	fig.update_yaxes(tickfont=dict(size=18),
					 titlefont=dict(size=20),
					 showgrid=True)
	output_html_basename = 'cgMLST.html' if not compute_accessory else 'cgMLST_plus_agMLST.html'
	output_html_path = os.path.join(output_directory, output_html_basename)
	plot(fig, filename=output_html_path, auto_open=False)
	print('Line plots for loci presence per threshold and step saved to {0}'.format(output_html_path))
