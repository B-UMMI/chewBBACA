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
import sys
import math

import numpy as np
import pandas as pd
import plotly.colors as pc
import plotly.express as px
from plotly.offline import plot
import plotly.graph_objects as go
from scipy.optimize import curve_fit

try:
	from utils import (constants as ct,
					   file_operations as fo,
					   iterables_manipulation as im,
					   multiprocessing_operations as mo)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (constants as ct,
								 file_operations as fo,
								 iterables_manipulation as im,
								 multiprocessing_operations as mo)


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


def compute_cgMLST(step, matrix, sorted_genomes, threshold, compute_accessory):
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
	# Get subdataframe for current genomes
	current_df = matrix.loc[sorted_genomes[:step]]
	pa_rows, _ = current_df.shape
	is_above_threshold = current_df.apply(above_threshold, args=(pa_rows, threshold))
	above = current_df.columns[is_above_threshold]
	cgMLST_size = (step, len(above))
	below = None
	agMLST_size = None
	if compute_accessory:
		# Compute accessory genome
		below = current_df.columns[~is_above_threshold]
		agMLST_size = (step, len(below))

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
		marginal="box", # Add a boxplot above the histogram
		color_discrete_sequence=['#225ea8'],
		hover_data=loci_presence_data.columns, # Set the hover data for the points above the histogram
	)

	# Set the bin size and the range for the histogram
	fig.update_traces(xbins=dict(size=0.011, start=0, end=1),
				   	  selector=dict(type="histogram"))

	fig.update_traces(boxpoints="all", # Show all the points
				    jitter=1,
					fillcolor="rgba(0,0,0,0)", # Hide the boxplot by making it transparent
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


def power_law_decay(x, k, alpha):
	"""Power law function for decay.

	Parameters
	----------
	x : array-like
		Input values.
	k : float
		Scaling factor.
	alpha : float
		Exponent parameter.

	Returns
	-------
	array-like
		Output values computed using the power law for decay.
	"""
	return k * (x**(-alpha))


def power_law_growth(x, k, alpha):
	"""Power law function for growth.

	Parameters
	----------
	x : array-like
		Input values.
	k : float
		Scaling factor.
	alpha : float
		Exponent parameter.

	Returns
	-------
	array-like
		Output values computed using the power law for growth.
	"""
	return k * (x ** alpha)


def perform_rarefaction_analysis(permutation_i, pa_matrix, sample_ids, n_samples, n_loci, threshold=None, analysis_type='core'):
	"""
	"""
	rarefaction_values = np.zeros((1, n_samples))
	shuffled_samples = np.random.choice(sample_ids, size=n_samples, replace=False)
	cumulative_presence = np.zeros(n_loci)
	for j, sample in enumerate(shuffled_samples):
		cumulative_presence += pa_matrix.loc[sample].values
		if analysis_type == 'core':
			rarefaction_values[0, j] = np.sum(cumulative_presence >= (threshold * (j + 1)))
		elif analysis_type == 'pangenome':
			rarefaction_values[0, j] = np.sum(cumulative_presence > 0)

	return rarefaction_values

	

	return X, Y_mean, Y_std, Y_predicted, alpha, alpha_error


def main(input_file, output_directory, threshold, step, compute_accessory, rarefaction_analysis,
		 permutation_number, permutation_samples, exclude_loci, exclude_genomes, cpu_cores):
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
	outdir_created = fo.create_directory(output_directory)
	if not outdir_created:
		sys.exit(ct.OUTPUT_DIRECTORY_EXISTS)

	# Import allelic profiles
	profiles = pd.read_csv(input_file, header=0, index_col=0, sep='\t', low_memory=False)

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

	# Use the fixed thresholds to compute the core genome for each threshold
	cgMLST_traces = []
	agMLST_traces = []
	for i, t in enumerate(cgMLST_thresholds):
		print(f'Analyzing results for threshold {t}...')
		# Use multiprocessing to determine set of core loci for each step in parallel
		step_ranges = [[step_i] for step_i in im.inclusive_range(1, len(sorted_genomes), step)]
		# Define the chunk size so that it creates groups of inputs instead of treating each input as a separate group which can be inefficient
		chunk_size = math.ceil(len(step_ranges)/cpu_cores)
		cgMLST_inputs = im.multiprocessing_inputs(step_ranges, [pa_matrix, sorted_genomes, t, compute_accessory], compute_cgMLST)
		cgMLST_results = mo.map_async_parallelizer(cgMLST_inputs, mo.function_helper, cpu_cores, show_progress=True, pool_type='threadpool', chunksize=chunk_size)

		# The set of core loci are determined on the last step value
		cgMLST_loci = cgMLST_results[-1][0]
		# Get number of core loci per step value
		cgMLST_counts = {r[1][0]: r[1][1] for r in cgMLST_results}
		print(f'\nBased on the threshold of {t}, the core genome is composed of {len(cgMLST_loci)}/{total_loci} loci.')

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

		current_color = custom_palette[i]
		# Create line plot for core genome values at each step value
		trace = go.Scattergl(x=list(cgMLST_counts.keys()),
							 y=list(cgMLST_counts.values()),
							 mode='lines',
							 line=dict(color=current_color),
							 name='cgMLST{0}'.format(int(t*100)),
							 hovertemplate=('%{y}'))
		cgMLST_traces.append(trace)

		if compute_accessory:
			# Get list of accessory loci
			agMLST_loci = cgMLST_results[-1][2]
			# Get size of accessory genome per step value
			agMLST_counts = {r[3][0]: r[3][1] for r in cgMLST_results}
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

	# Use rarefaction analysis to determine if the core stabilizes
	if rarefaction_analysis:
		if permutation_samples:
			n_samples = permutation_samples
		else:
			n_samples = total_genomes
		n_loci = total_loci

		n_permutations = permutation_number
		sample_ids = pa_matrix.index.tolist()
		rarefaction_traces = []
		analysis_type = 'core'
		for i, t in enumerate(cgMLST_thresholds):
			print(f'Performing rarefaction analysis to measure core genome stability for threshold {t}...')
			print(f'Using {n_permutations} permutations and {n_samples} samples for the rarefaction analysis.')

			chunk_size = math.ceil(n_permutations/cpu_cores)
			permutations = [[perm_i] for perm_i in list(range(n_permutations))]
			rarefaction_inputs = im.multiprocessing_inputs(permutations, [pa_matrix, sample_ids, n_samples, n_loci, t, analysis_type], perform_rarefaction_analysis)
			rarefaction_results = mo.map_async_parallelizer(rarefaction_inputs, mo.function_helper, cpu_cores, show_progress=True, pool_type='threadpool', chunksize=chunk_size)

			# Concatenate the results for all permutations
			rarefaction_perms = np.vstack(rarefaction_results)

			X = np.arange(1, n_samples + 1)
			Y_mean_core = np.mean(rarefaction_perms, axis=0)
			Y_std_core = np.std(rarefaction_perms, axis=0)
		
			# Use non-linear least squares to fit a power law function to the data
			fit_function = power_law_decay if analysis_type == 'core' else power_law_growth
			# Define the initial value for k as the value for the first empirical observation
			# Define the initial value for the exponent as 0.5 so that the initial value is within the valid interval and indicates there is some change
			# May need to set maxfev value to increase the number of times the parameter values are tweaked with very noisy data
			# Parameter value bounds are set to 0 as minimum and unbounded for upper value
			popt, pcov = curve_fit(fit_function, X, Y_mean_core, p0=[Y_mean_core[0], 0.5], bounds=[[0, 0], [np.inf, np.inf]])
			# Get optimized parameter values
			kappa, alpha_core = popt
			# Get values for the fitted curve
			Y_predicted = fit_function(X, kappa, alpha_core)
			# Get estimate for the variance of the fit for the exponent
			# The elements in the diagonal of pcov represent the variance of the parameters, which we can use to compute the standard error of the exponent
			# The non-diagonal elements of pcov represent the covariance between the parameters
			# Get only the standard error for exponent alpha
			_, alpha_core_error = np.sqrt(np.diag(pcov))

			# Print information about decay exponent α
			# α > 0 indicates that the core-genome is not completely stable, with higher values pointing to a more steep drop
			# Values close to 0 indicate the core-genome is stabilizing, while high values may be related to the strains not sharing a lot of the genes or the data being of low quality
			# High values can make it drop to 0 very quickly
			print(f'\nThe value of the decay exponent (α) for the core-genome determined based on a threshold of {t} is {alpha_core:.3f} ± {alpha_core_error:.3f}')

			# Create traces for empirical permutation averages and Power law curve
			current_color = custom_palette[i]
			# Trace for empirical permutation averages
			# Need to include "<extra></extra>" in the hovertemplate so that the trace name is not shown next to the hover label
			mean_trace_core = go.Scatter(x=X, y=Y_mean_core,
									   mode='lines',
									   name=f'Core-genome Permutation Means (t={t})',
									   line=dict(color=current_color),
									   hovertemplate=f"Core-genome Permutation Means (t={t})<br># samples: %{{x}}<br># loci: %{{y}}<extra></extra>")

			# Trace for the standard deviation values
			y_upper = Y_mean_core + Y_std_core
			y_lower = Y_mean_core - Y_std_core
			std_color = current_color.replace('rgb', 'rgba').replace(')', ', 0.35)')
			std_trace_core = go.Scatter(x=np.concatenate([X, X[::-1]]),
							   y=np.concatenate([y_upper, y_lower[::-1]]),
							   fill='toself',
							   fillcolor=std_color,
							   line=dict(color='rgba(255,255,255,0)'),
							   hoverinfo="skip",
							   showlegend=False
							   )

			# Trace for Power law curve
			pcore_trace = go.Scatter(x=X, y=Y_predicted,
							mode='lines',
							name=f"Core-genome Power Law Fit (t={t}, α={alpha_core:.3f})",
							line=dict(color=current_color, dash='dot'),
							hovertemplate=f"Power Law Fit (t={t}, α={alpha_core:.3f})<br># samples: %{{x}}<br># loci: %{{y}}<extra></extra>")

			rarefaction_traces.extend([mean_trace_core, std_trace_core, pcore_trace])

		print(f'Performing rarefaction analysis to evaluate pangenome stability...')
		print(f'Using {n_permutations} permutations and {n_samples} samples for the rarefaction analysis.')

		analysis_type = 'pangenome'
		chunk_size = math.ceil(n_permutations/cpu_cores)
		permutations = [[perm_i] for perm_i in list(range(n_permutations))]
		rarefaction_inputs = im.multiprocessing_inputs(permutations, [pa_matrix, sample_ids, n_samples, n_loci, t, analysis_type], perform_rarefaction_analysis)
		rarefaction_results = mo.map_async_parallelizer(rarefaction_inputs, mo.function_helper, cpu_cores, show_progress=True, pool_type='threadpool', chunksize=chunk_size)

		# Concatenate the results for all permutations
		rarefaction_perms = np.vstack(rarefaction_results)

		X = np.arange(1, n_samples + 1)
		Y_mean_pan = np.mean(rarefaction_perms, axis=0)
		Y_std_pan = np.std(rarefaction_perms, axis=0)

		# Use non-linear least squares to fit a power law function to the data
		fit_function = power_law_decay if analysis_type == 'core' else power_law_growth
		# Define the initial value for k as the value for the first empirical observation
		# Define the initial value for the exponent as 0.5 so that the initial value is within the valid interval and indicates there is some change
		# May need to set maxfev value to increase the number of times the parameter values are tweaked with very noisy data
		# Parameter value bounds are set to 0 as minimum and unbounded for upper value
		popt, pcov = curve_fit(fit_function, X, Y_mean_pan, p0=[Y_mean_pan[0], 0.5], bounds=[[0, 0], [np.inf, np.inf]])
		# Get optimized parameter values
		kappa, gamma_pan = popt
		# Get values for the fitted curve
		Y_pan_predicted = fit_function(X, kappa, gamma_pan)
		# Get estimate for the variance of the fit for the exponent
		# The elements in the diagonal of pcov represent the variance of the parameters, which we can use to compute the standard error of the exponent
		# The non-diagonal elements of pcov represent the covariance between the parameters
		# Get only the standard error for exponent alpha
		_, gamma_pan_error = np.sqrt(np.diag(pcov))

		# Print information about the growth exponent
		# Values above 0 mean that new genes keep being discovered as we add genomes
		# Values closer to 0 mean that fewer genes are discovered, with 0 meaning that the pangenome is closed
		print(f'\nThe growth exponent (y) for the pangenome is {gamma_pan:.3f} ± {gamma_pan_error:.3f}.')

		# Create traces for empirical permutation averages and Power law curve
		mean_trace_pan = go.Scatter(x=X, y=Y_mean_pan,
									 mode='lines',
									 name=f'Pangenome Permutation Means',
									 line=dict(color='black'),
									 hovertemplate=f"Pangenome Permutation Means<br># samples: %{{x}}<br># loci: %{{y}}<extra></extra>")

		# Trace for the standard deviation values
		y_upper = Y_mean_pan + Y_std_pan
		y_lower = Y_mean_pan - Y_std_pan
		std_trace_pan = go.Scatter(x=np.concatenate([X, X[::-1]]),
							y=np.concatenate([y_upper, y_lower[::-1]]),
							fill='toself',
							fillcolor='rgba(211, 211, 211, 0.35)',
							line=dict(color='rgba(255,255,255,0)'),
							hoverinfo="skip",
							showlegend=False
							)

		ppan_trace = go.Scatter(x=X, y=Y_pan_predicted,
						  mode='lines',
						  name=f"Pangenome Power Law Fit (y={gamma_pan:.3f})",
						  line=dict(color='black', dash='dot'),
						  hovertemplate=f"Power Law Fit (y={gamma_pan:.3f})<br># samples: %{{x}}<br># loci: %{{y}}<extra></extra>")

		rarefaction_traces.extend([mean_trace_pan, std_trace_pan, ppan_trace])

		fig = go.Figure()
		fig.add_traces(rarefaction_traces)

		fig.update_layout(title={'text': 'Core-genome and pangenome rarefaction analysis',
							 'font_size': 30},
							 template='simple_white',
							 hovermode='x',
							 xaxis=dict(title=dict(text='Number of samples', font=dict(size=20)),
                                      tickfont=dict(size=18),
                                      showgrid=True,
                                      range=[0, n_samples]),
							 yaxis=dict(title=dict(text='Number of loci', font=dict(size=20)),
                                      tickfont=dict(size=18),
                                      showgrid=True)
		)

		output_html_basename = 'rarefaction_analysis.html'
		output_html_path = os.path.join(output_directory, output_html_basename)
		plot(fig, filename=output_html_path, auto_open=False)
		print(f'Line plots for the rarefaction analysis saved to {output_html_path}')
