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
import plotly.colors as pc
import plotly.express as px
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
		Pandas Series object.

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
		Pandas Series object with presence (1) and absence (0) values.
	column_length : int
		Total number of values/rows.
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
		Whether to compute the accessory genome in addition to the core genome.

	Returns
	-------
	pruned_df : pandas.core.frame.DataFrame
		Dataframe with the cgMLST profiles for the last step value
		(the last step includes all genomes in the input matrix).
	cgMLST_size : dict
		Dictionary with the number of genomes used to compute the
		cgMLST as keys and the size of the core-genome as values.
	"""
	# Determine loci at or above threshold
	above = None
	cgMLST_size = {}
	below = None
	agMLST_size = {}
	for i in im.inclusive_range(1, len(sorted_genomes), step):
		# Get subdataframe for current genomes
		current_df = matrix.loc[sorted_genomes[:i]]
		pa_rows, _ = current_df.shape
		is_above_threshold = current_df.apply(above_threshold, args=(pa_rows, threshold))
		above = current_df.columns[is_above_threshold]
		cgMLST_size[pa_rows] = len(above)
		if compute_accessory:
			# Compute accessory genome
			below = current_df.columns[~is_above_threshold]
			agMLST_size[pa_rows] = len(below)

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


def compute_stats(matrix, axis=1, column_names=None):
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

	zeros_data = [matrix.index if axis == 1 else matrix.columns, # labels used in the first column
			   	  non_zeros, # Number of loci identified or samples with each locus
				  (non_zeros/ncols) if axis == 1 else (non_zeros/nrows)] # Percentage of loci identified or samples with each locus 

	zeros_df = pd.DataFrame(list(zip(*zeros_data)), columns=column_names)

	return zeros_df


def plot_loci_presence(loci_presence_data, output_directory):
	"""Plot the distribution of loci presence frequency percentages.

	Parameters
	----------
	loci_presence_data : pandas.core.frame.DataFrame
		Dataframe with the number and percentage of samples with each locus.
	output_directory : str
		Path to the directory where the HTML file with
		the plot will be stored.

	Returns
	-------
	output_html_path : str
		Path to the output HTML file that contains the
		plot with the distribution of loci presence frequency percentages.s
	"""
	fig = px.histogram(
		loci_presence_data,
		x='Sample presence proportion',
		marginal="box",
		nbins=200,
		color_discrete_sequence=['#225ea8'],
		hover_data=loci_presence_data.columns, # Set the hover data fr the points in the points above the histogram
	)

	fig.update_traces(boxpoints="all", jitter=1,
					fillcolor="rgba(0,0,0,0)",
					line_color="rgba(0,0,0,0)",
					marker=dict(color="#225ea8", size=3),
					selector=dict(type="box"))

	fig.update_layout(title={'text': 'Loci presence proportion', 'font_size': 30},
					template='simple_white',
					yaxis=dict(title=dict(text='Count', font=dict(size=20)),
									tickfont=dict(size=18), showgrid=True,
									domain=[0, 0.90]),
					xaxis=dict(title=dict(text='Locus presence proportion', font=dict(size=20)),
									tickfont=dict(size=18), showgrid=True),
					showlegend=False,
					yaxis2=dict(domain=[0.93, 1.0], # Reduce the vertical space for the trace with the invisible boxplot
					range=[-0.75, 0])
	)

	output_html_basename = 'loci_presence_proportion.html'
	output_html_path = os.path.join(output_directory, output_html_basename)
	plot(fig, filename=output_html_path, auto_open=False)

	return output_html_path


def main(input_file, output_directory, threshold, step,
		 compute_accessory, exclude_loci, exclude_genomes):
	"""Determine the cgMLST based on allele calling results.

	Parameters
	----------
	input_file : str
		Path to the TSV file that contains the allelic 
		profiles determined by the AlleleCall module.
	output_directory : str
		Path to the directory where the process will
		store output files.
	threshold : list
		Loci presence threshold values used to determine the set of core loci.
	step : int
		Number of genomes added to the core genome computation at
		each step.
	compute_accessory : bool
		Determine the set of accessory loci. The accessory genome 
		corresponds to all the loci not included in the core genome. 
		The accessory genome is determined for each core genome threshold.
	exclude_loci : str
		Path to TXT file with a list of loci to exclude from the analysis.
	exclude_genomes : str
		Path to TXT file with a list of genomes to exclude from the analysis.
	"""
	fo.create_directory(output_directory)

	# Import allelic profiles
	profiles = pd.read_csv(input_file, header=0, index_col=0,
						 sep='\t', low_memory=False)

	# Get number of genomes and loci
	total_genomes, total_loci = profiles.shape
	print(f'Input file includes profiles for {total_genomes} samples and {total_loci} loci.')

	cgMLST_thresholds = sorted(threshold)
	print(f'Will compute the set of core loci based on the following loci presence thresholds: {", ".join(map(str, cgMLST_thresholds))}')

	# Read list of samples to exclude from the analysis
	if exclude_genomes:
		genomes_to_remove = fo.read_lines(exclude_genomes)
		print(f'User provided list with {len(genomes_to_remove)} samples to exclude.')
		print(f'List of genomes to exclude: {", ".join(genomes_to_remove)}')
		profiles = remove_genomes(profiles, genomes_to_remove)
		print(f'Excluded {total_genomes - profiles.shape[0]} genomes.')

	if exclude_loci:
		loci_to_remove = fo.read_lines(exclude_loci)
		print(f'User provided list with {len(loci_to_remove)} loci to exclude.')
		profiles = remove_columns(profiles, loci_to_remove)
		print(f'Excluded {total_loci - profiles.shape[1]} loci.')

	if exclude_genomes or exclude_loci:
		total_genomes, total_loci = profiles.shape
		print(f'After excluding {len(genomes_to_remove)} sample profiles and {len(loci_to_remove)} '
			  f'loci columns, the profiles to analyze include {total_genomes} samples and {total_loci} loci.')

	# Mask special classifications and remove 'INF-' prefixes
	print('Masking profiles...')
	masked_profiles = profiles.apply(im.replace_chars)
	print(f'Masked {total_genomes} profiles.')

	# Compute presence-absence matrix
	print('Computing presence-absence matrix...')
	pa_matrix, pa_outfile = compute_presence_absence(masked_profiles, output_directory)
	print(f'Presence-absence matrix saved to {pa_outfile}')

	# Count number of special classifications per genome
	print('Computing the presence-absence statistics per sample...')
	genome_mdata_df = compute_stats(pa_matrix, column_names=ct.GENOMES_MISSING_COLUMNS)
	# Sort based on increasing number of special classifications
	genome_mdata_df = genome_mdata_df.sort_values('Loci presence count', ascending=False)
	sorted_genomes = genome_mdata_df['Sample'].tolist()
	# Round percentage values to 2 decimal places
	genome_mdata_df['Loci presence proportion'] = genome_mdata_df['Loci presence proportion'].round(3)

	# Count number of special classifications per locus
	print('Computing the presence-absence statistics per locus...')
	loci_mdata_df = compute_stats(pa_matrix, axis=0, column_names=ct.LOCI_MISSING_COLUMNS)
	# Sort based on increasing number of special classifications
	loci_mdata_df = loci_mdata_df.sort_values('Sample presence count', ascending=False)
	# Round percentage values to 2 decimal places
	loci_mdata_df['Sample presence proportion'] = loci_mdata_df['Sample presence proportion'].round(3)

	# Create custom color palette for line plots
	# Select the color at the midpoint of the colorscale if only one threshold is provided, otherwise sample colors from the colorscale
	# Providing a float value allows to get the color at a specific percentage mapping
	# Providing 1 when there is a single threshold would lead to a ZeroDivisionError because it assumes we want to sample at least 2 colors
	if len(cgMLST_thresholds) == 1:
		midpoint_color = pc.sample_colorscale("Cividis", 0.5)[0]
		custom_palette = [midpoint_color]
	else:
		custom_palette = pc.sample_colorscale("Cividis", len(cgMLST_thresholds))

	# Compute the core genome for each threshold
	cgMLST_traces = []
	agMLST_traces = []
	for i, t in enumerate(cgMLST_thresholds):
		print(f'Analyzing results for threshold {t}...')
		current_color = custom_palette[i]
		cgMLST_results = compute_cgMLST(pa_matrix, sorted_genomes, t, step, compute_accessory)
		cgMLST_loci, cgMLST_counts, agMLST_loci, agMLST_counts = cgMLST_results
		print(f'Based on the threshold of {t}, the core genome is composed of {len(cgMLST_loci)}/{total_loci} loci.')

		# Write cgMLST matrix
		# Get subset from masked matrix
		cgMLST_matrix = masked_profiles[cgMLST_loci]
		cgMLST_path = os.path.join(output_directory, 'cgMLST{0}.tsv'.format(int(t*100)))
		cgMLST_matrix.to_csv(cgMLST_path, sep='\t')
		print(f'cgMLST profiles for threshold {t} saved to {cgMLST_path}')

		# Write list of cgMLST loci
		cgMLST_loci_path = os.path.join(output_directory, 'cgMLSTschema{0}.txt'.format(int(t*100)))
		fo.write_lines(list(cgMLST_loci), cgMLST_loci_path)
		print(f'List of core loci for threshold {t} saved to {cgMLST_loci_path}')

		# Determine the percentage of cgMLST loci present in each genome
		# Get the values for each genome in the order of the sorted_genomes list, cast them to int32 and compute the percentage of core loci in each genome
		cgMLST_sample_pct = [np.sum(cgMLST_matrix.loc[genome].astype(np.int32) > 0) / len(cgMLST_loci) for genome in sorted_genomes]
		genome_mdata_df[f'Loci presence proportion (cgMLST{int(t*100)})'] = cgMLST_sample_pct
		genome_mdata_df[f'Loci presence proportion (cgMLST{int(t*100)})'] = genome_mdata_df[f'Loci presence proportion (cgMLST{int(t*100)})'].round(3)

		# Create line plot for core genome values at each step value
		trace = go.Scattergl(x=list(cgMLST_counts.keys()),
							 y=list(cgMLST_counts.values()),
							 mode='lines',
							 line=dict(color=current_color),
							 name='cgMLST{0}'.format(int(t*100)),
							 hovertemplate=('%{y}'))
		cgMLST_traces.append(trace)

		if compute_accessory:
			print(f'Based on the threshold of {t}, the accessory genome is composed of {len(agMLST_loci)}/{total_loci} loci.')
			# Write agMLST matrix
			# Get subset from masked matrix
			agMLST_matrix = masked_profiles[agMLST_loci]
			agMLST_path = os.path.join(output_directory, 'agMLST{0}.tsv'.format(int(t*100)))
			agMLST_matrix.to_csv(agMLST_path, sep='\t')
			print(f'agMLST profiles for threshold {t} saved to {agMLST_path}')

			# Write list of agMLST loci
			agMLST_loci_path = os.path.join(output_directory, 'agMLSTschema{0}.txt'.format(int(t*100)))
			fo.write_lines(list(agMLST_loci), agMLST_loci_path)
			print(f'List of accessory loci for threshold {t} saved to {agMLST_loci_path}')

			# Create line plot for accessory genome values at each step value
			trace = go.Scattergl(x=list(agMLST_counts.keys()),
								 y=list(agMLST_counts.values()),
								 mode='lines',
								 line=dict(dash='dash', color=current_color),
								 name='agMLST{0}'.format(int(t*100)),
								 hovertemplate=('%{y}'))
			agMLST_traces.append(trace)

	# Write TSV with presence-absence statistics per sample
	genome_mdata_path = os.path.join(output_directory, ct.GENOMES_MISSING_BASENAME)
	genome_mdata_df.to_csv(genome_mdata_path, sep='\t', index=False)
	print(f'Saved presence-absence statistics per sample to {genome_mdata_path}')

	# Write TSV with presence-absence statistics per locus
	loci_mdata_path = os.path.join(output_directory, ct.LOCI_MISSING_BASENAME)
	loci_mdata_df.to_csv(loci_mdata_path, sep='\t', index=False)
	print(f'Saved presence-absence statistics per locus to {loci_mdata_path}')

	# Plot loci presence frequency percentage
	print('Creating plot with the distribution of loci presence values...')
	output_html_path = plot_loci_presence(loci_mdata_df, output_directory)
	print(f'Saved plot to {output_html_path}')

	print('Creating line plots for each threshold and step value...')
	# Create trace with present loci
	present = genome_mdata_df['Loci presence count'].tolist()
	genomes_index = list(range(1, len(present)+1))
	miss_trace = go.Scattergl(x=genomes_index,
							  y=present,
							  mode='lines',
							  name='No. of loci per sample',
							  line=dict(dash='dot', color='#000000'),
							  hovertemplate=('%{y}<br>'
											 'Sample: %{text}<br>'),
							  text=sorted_genomes)

	traces = cgMLST_traces + agMLST_traces + [miss_trace]

	fig = go.Figure(data=traces)
	fig.update_layout(title={'text': 'Size of the core genome for each threshold and step value',
							 'font_size': 30},
					  template='simple_white',
					  hovermode='x',
					  yaxis=dict(title=dict(text='Number of loci', font=dict(size=20)),
                                      tickfont=dict(size=18),
                                      showgrid=True),
					  xaxis=dict(title=dict(text='Number of samples', font=dict(size=20)),
                                      tickfont=dict(size=18),
                                      showgrid=True,
                                      range=[0, len(sorted_genomes)])
	)

	output_html_basename = 'cgMLST.html' if not compute_accessory else 'cgMLST_plus_agMLST.html'
	output_html_path = os.path.join(output_directory, output_html_basename)
	plot(fig, filename=output_html_path, auto_open=False)
	print(f'Line plots saved to {output_html_path}')
