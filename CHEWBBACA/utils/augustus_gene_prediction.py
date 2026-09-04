#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions related to gene prediction with AUGUSTUS.

Code documentation
------------------
"""


import sys
import subprocess

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


def parse_augustus_gff(augustus_gff, genome_basename):
	"""
	"""
	protid = 1
	gene_info = {}
	with open(augustus_gff, 'r') as infile:
		current_cds_id = None
		reading_cds = False
		for line in infile:
			if line.startswith("# start gene"):
				gene_id = line.split("start gene ")[-1]
				contig_id = gene_id.rsplit(".", 1)[0]
				cds_id = f"{genome_basename}_{protid}"
				gene_info[cds_id] = [genome_basename, contig_id]
				current_cds_id = cds_id
				continue
			if "AUGUSTUS\tgene" in line:
				start, end, score, strand = line.split("\t")[3:7]
				strand = '1' if strand == '+' else '-1'
				gene_info[current_cds_id].extend([start, end, str(protid), strand, score])
				protid += 1
				continue
			if line.startswith("# coding sequence"):
				cds_start = (line.strip()).split("[")[1]
				gene_info[current_cds_id].append(cds_start.upper())
				reading_cds = True
				continue
			if reading_cds:
				if "]" not in line:
					cds_rest = (line.strip()).split("# ")[1]
					gene_info[current_cds_id][-1] += cds_rest.upper()
				else:
					cds_rest = (line.strip()).split("# ")[1][:-1]
					gene_info[current_cds_id][-1] += cds_rest.upper()
					# Add SHA256 hash
					sha256_hash = im.hash_sequence(gene_info[current_cds_id][-1], 'sha256')
					gene_info[current_cds_id].append(sha256_hash)
					reading_cds = False

	return gene_info


def predict_genome_genes(input_file, coordinates_outfile, output_directory, augustus_path, species, output_formats):
	"""Predict genes for sequences in a FASTA file.
	"""
	# Get genome unique identifier
	genome_basename = input_file[1]

	# Create output file paths
	output_file = fo.join_paths(output_directory, [f"{genome_basename}.gff"])
	# Add `--genemodel=complete` to predict only complete genes and does not identify incomplete genes in sequence boundaries
	# Add `--noInFrameStop=true` so it does not report CDSs with in-frame stop codons
	# Add `--gff3=true` to output in GFF format, including the stop codon
	# Include `--uniqueGeneId=true` to use unique gene names and avoid conflicting gene names when parsing results for multiple inputs
	augustus_cmd = [augustus_path, f"--species={species}", "--genemodel=complete", "--noInFrameStop=true",
				 	"--gff3=true", "--uniqueGeneId=true", "--codingseq=true", "--protein=false",
					f"--outfile={output_file}", input_file[0]]

	# AUGUSTUS outputs a warning to stderr when using `--gff3=true` which causes to process to hang indefinitely
	# To avoid this, redirect stderr to stdout and capture the output to still be able to check for errors in the output
	# The warning is the following:
	# "Warning: The gff3 standard requires that the stop codon is included in the CDS. Unless this is your intention, 
	# set stopCodonExcludedFromCDS to false in your species' configuration file or on the command line.""
	augustus_process = subprocess.Popen(augustus_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
	stdout, _ = augustus_process.communicate()

	## Need to implement logging and write captured stdout to the log file

	# Parse results in non-standard GFF3 format
	# A custom parser had to be created to extract the data and allow to create output files
	# equivalent to what is created when Pyrodigal is used to ensure consistent output formats
	gene_info = parse_augustus_gff(output_file, genome_basename)

	total_genome = len(gene_info)

	if "gff" not in output_formats:
		fo.remove_files([output_file])
	else:
		# Copy file to subfolder for GFF output format
		gff_outdir = fo.join_paths(output_directory, ["gff"])
		fo.copy_file(output_file, gff_outdir)

	output_files = [None, output_file] if "gff" in output_formats else [None]
	# Delete GFF file if it is not requested as an output format

	if total_genome > 0:
		# Create FASTA file with the predicted CDSs
		outfasta = fo.join_paths(output_directory, ["genes", f"{genome_basename}.fasta"])
		records = []
		for seqid, cds_data in gene_info.items():
			current_record = ct.FASTA_RECORD_TEMPLATE.format(seqid, cds_data[-2])
			records.append(current_record)
		fo.write_lines(records, outfasta)
		output_files[0] = outfasta

	# Append gene coordinates to process TSV file
	coordinate_data = []
	for seqid, cds_data in gene_info.items():
		coordinate_data.append([seqid]+cds_data[:7]+[cds_data[-1]])
	coordinate_outlines = ["\t".join(l) for l in coordinate_data]
	fo.write_lines(coordinate_outlines, coordinates_outfile, write_mode="a")

	# No detection of PLOT classes with AUGUSTUS because the command includes `--genemodel=complete`
	# To avoid identifying partial genes near contig boundaries
	close_to_tip = {genome_basename: {}}

	return [input_file, total_genome, close_to_tip, output_files]


def predict_genes(fasta_files, output_directory, gene_prediction_parameters, cpu_cores):
	"""
	"""
	augustus_inputs = im.divide_list_into_n_chunks(fasta_files, len(fasta_files))
	# Add path to TSV file to store CDS coordinate data per input subset
	for i in range(len(augustus_inputs)):
		augustus_inputs[i].append(fo.join_paths(output_directory, [f'gene_coordinates_{i}']))

	# Create subfolders to store each output format
	for output_format in gene_prediction_parameters["output_formats"]:
		format_outdir = fo.join_paths(output_directory, [output_format])
		fo.create_directory(format_outdir)

	# Add common arguments to all sublists
	augustus_exe = fo.join_paths(gene_prediction_parameters["augustus_path"], [ct.AUGUSTUS_ALIAS])
	common_args = [output_directory, augustus_exe,
				   gene_prediction_parameters["species"], gene_prediction_parameters["output_formats"]]
	augustus_inputs = im.multiprocessing_inputs(augustus_inputs, common_args, predict_genome_genes)

	# Run Pyrodigal to predict genes
	# Need to use ThreadPool. Pyrodigal might hang when using Pool
	augustus_results = mo.map_async_parallelizer(augustus_inputs,
												 mo.function_helper,
												 cpu_cores,
												 show_progress=True,
												 pool_type='threadpool')

	# Get number of inputs for which gene prediction failed
	# Inputs with 0 CDSs and inputs with error messages
	failed = {line[0][0]: line[1] for line in augustus_results if line[1] == 0}

	# Get number of CDSs predicted per valid input
	cds_counts = {line[0][1]: line[1] for line in augustus_results if isinstance(line[1], int) is True}

	# Get total number of CDSs predicted
	total_cds = sum(cds_counts.values())

	# Get paths to FASTA files with the extracted CDSs
	cds_fastas = [line[-1][0] for line in augustus_results if line[-1][0] is not None]

	# Merge TSV files with CDS coordinates
	coordinate_files = [i[1] for i in augustus_inputs]
	merged_coordinates = fo.join_paths(output_directory, ['gene_coordinates.tsv'])
	fo.concatenate_files(coordinate_files, merged_coordinates, header=ct.GENE_TABLE_HEADER)
	# Delete separate coordinate files
	fo.remove_files(coordinate_files)

	# Merge dictionaries with info about CDSs close to contig tips
	close_to_tip = [line[2] for line in augustus_results if len(line[2]) > 0]
	close_to_tip = im.merge_dictionaries(close_to_tip)

	return [failed, total_cds, cds_fastas, merged_coordinates, cds_counts, close_to_tip]
