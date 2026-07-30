#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains functions used to work with FASTA files.

Code documentation
------------------
"""


import warnings

from Bio import SeqIO
from Bio.SeqIO import FastaIO
# Import warning class to ignore warnings during FASTA file validation
from Bio import BiopythonDeprecationWarning

try:
	from utils import (file_operations as fo,
					   iterables_manipulation as im,
					   sequence_manipulation as sm,
					   constants as ct)
except ModuleNotFoundError:
	from CHEWBBACA.utils import (file_operations as fo,
								 iterables_manipulation as im,
								 sequence_manipulation as sm,
								 constants as ct)


def sequence_generator(input_file):
	"""Create a SeqRecord iterator.

	Parameters
	----------
	input_file : str
		Path to a Fasta file.

	Returns
	-------
	records : Bio.SeqIO.FastaIO.FastaIterator
		SeqRecord iterator.
	"""
	# Useful to create the generator
	# Need to exhaust the generator to avoid high memory usage
	records = SeqIO.parse(input_file, 'fasta')

	return records


def index_fasta(fasta_file, indexed_file=False):
	"""Create index to retrieve records from a FASTA file more efficiently.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.
	indexed_file : bool
		If True, use SeqIO.index_db to store the record information
		as a file on disk. Otherwise, use SeqIO.index to create
		the index in memory. A file on disk is more efficient for
		large FASTA files and can be used with multiprocessing by
		recalling this function inside each process.

	Returns
	-------
	fasta_index : Bio.File._IndexedSeqFileDict
		FASTA file index.
	"""
	if indexed_file:
		# Define the index file name
		index_file = fasta_file.replace('fasta', 'idx')
		# Create the index file on disk
		fasta_index = SeqIO.index_db(index_file, fasta_file, 'fasta')
	else:
		fasta_index = SeqIO.index(fasta_file, 'fasta')

	return fasta_index


def import_sequences(input_file):
	"""Import sequences from a FASTA file.

	Parameters
	----------
	input_file : str
		Path to a FASTA file.

	Returns
	-------
	records_dict : dict
		Dictionary with sequence identifiers as keys and
		sequences as values.
	"""
	records = sequence_generator(input_file)
	# Only want record identifier and sequence, no need to use SeqIO.to_dict
	records_dict = {rec.id: str(rec.seq.upper()) for rec in records}

	return records_dict


def count_sequences(fasta_file):
	"""Count the number of sequences in a FASTA file.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.

	Returns
	-------
	total_seqs : int
		Number of sequences in the FASTA file.
	"""
	with open(fasta_file, 'r') as infile:
		sequence_headers = [line for line in infile if line.startswith('>')]
		total_seqs = len(sequence_headers)

	return total_seqs


def write_records(records, output_file, mode='w'):
	"""Write FASTA records (BioPython SeqRecord) to a file.

	Parameters
	----------
	records : list
		List with BioPython SeqRecord objects.
	output_file : str
		Path to the output file.
	"""
	with open(output_file, 'w') as output_handle:
		fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
		fasta_out.write_file(records)


def integer_headers(input_fasta, output_fasta, start=1, limit=50000,
					prefix='', id_map=True, delete_input=False):
	"""Switch sequence headers in Fasta file by integer values.

	Parameters
	----------
	input_fasta : str
		Path to a FASTA file.
	output_fasta : str
		Path to the output file with modified headers.
	start : int
		Integer value for the first identifier.
	limit : int
		Maximum number of FASTA records to keep in
		memory.
	prefix : str
		Prefix added to the sequence identifiers.
	id_map : bool
		If True, the function returns a dictionary
		mapping the original sequence IDs to the
		new sequence IDs. If False, returns the
		number of sequences renamed.

	Returns
	-------
	ids_map : dict
		Dictionary with mapping between integer and original
		headers.
	"""
	seqs = []
	ids_map = {}
	exhausted = False
	seq_generator = sequence_generator(input_fasta)
	while exhausted is False:
		record = next(seq_generator, None)
		if record is not None:
			new_id = '{0}{1}'.format(prefix, start)
			ids_map[new_id] = record.id
			sequence = str(record.seq)
			new_rec = fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
									   [new_id, sequence])
			seqs.append(new_rec)
			start += 1
		elif record is None:
			exhausted = True

		if len(seqs) == limit or exhausted is True:
			fo.write_lines(seqs, output_fasta, write_mode='a')
			seqs = []

	# Remove input file if requested
	if delete_input is True:
		fo.remove_files([input_fasta])

	if id_map is True:
		return ids_map
	else:
		return (input_fasta, len(ids_map))


def fasta_str_record(record_template, record_data):
	"""Create the string representation of a FASTA record.

	Parameters
	----------
	record_template : str
		String template to construct the FASTA record.
	record_data : list
		List with the elements to add to the string.

	Returns
	-------
	record : str
		String representation of the FASTA record.
	"""
	record = record_template.format(*record_data)

	return record


def fasta_lines(template, records_data):
	"""Create a list with FASTA records.

	Parameters
	----------
	template : str
		String template to construct the FASTA record.
	records_data : list
		A list with one sublist per FASTA record.
		Each sublist contains the elements to insert
		inside the template placeholders.

	Returns
	-------
	seqs_lines : list
		A list with strings representing FASTA records.
	"""
	seqs_lines = [fasta_str_record(template, arg) for arg in records_data]

	return seqs_lines


def validate_fasta(file_path):
	"""Check if a file is a FASTA file.

	Parameters
	----------
	file_path : str
		Path to the file.

	Returns
	-------
	True if file is a valid FASTA, False otherwise.
	"""
	# Should return and empty list if file is empty or has any problem
	try:
		# Need to ignore warning if it is not a FASTA file or if it is but has comments or special characters at the start of the file
		# Biopython used to simply return an empty list for non-FASTA files and ignore comments or special characters, but Biopython 1.86 raises a deprecation warning
		# Need to use the warnings module to ignore warnings of category BiopythonDeprecationWarning
		with warnings.catch_warnings():
			warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)
			# Need to convert to list to exhaust generator
			# Otherwise it might increase memory usage
			records = list(sequence_generator(file_path))
	# Future Biopython versions will raise a ValueError and this try/except block was added to take that into consideration
	# This assumes that either it is not a FASTA file or has comments or special characters at the start
	except ValueError:
		records = []

	return any(records)


def filter_non_fasta(files):
	"""Select FASTA files from a list with file paths.

	Parameters
	----------
	files : list
		A list that contains file paths.

	Returns
	-------
	fasta_files : list
		List that contains paths to FASTA files.
	"""
	fasta_files = [file for file in files if validate_fasta(file) is True]

	invalid_files = list(set(files)-set(fasta_files))

	return fasta_files, invalid_files


def sequence_lengths(fasta_file, hashed=False):
	"""Determine length of sequences in a FASTA file.

	Read Fasta file and create dictionary with mapping
	between sequence identifiers and sequence lengths.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.
	hashed : bool
		If False, sequence headers are used as
		keys. If True, sequence hashes will be
		used as keys.

	Returns
	-------
	lengths : dict
		Dictionary with sequence identifiers as keys and
		sequence lengths as values.
	"""
	records = sequence_generator(fasta_file)
	if hashed is False:
		lengths = {rec.id: len(rec.seq) for rec in records}
	else:
		lengths = {im.hash_sequence(str(rec.seq)):
				   len(rec.seq) for rec in records}

	return lengths


def get_sequences_by_id(sequences, seqids, output_file, limit=50000):
	"""Retrieve sequences from indexed FASTA file or dictionary.

	Retrieves sequences based on sequence identifiers from a FASTA
	file indexed with SeqIO.index or from a dictionary with sequence
	identifiers as keys and sequences as values.

	Parameters
	----------
	sequences : dict or Bio.File._IndexedSeqFileDict
		Dictionary with seqids as keys and sequences as values or
		a FASTA file index created with BioPython.
	seqids : list
		List with the identifiers of the sequences that should be
		retrieved.
	output_file : str
		Path to the FASTA file to which selected sequences will
		be saved.
	limit : int
		Maximum number of sequences that will be kept in memory
		at a time (to avoid keeping huge datasets in memory).

	Returns
	-------
	total_selected : int
		Total number of records written to the output file.
	"""
	# using a generator
	if type(sequences) == dict:
		seqs = ((seqid, sequences[seqid]) for seqid in seqids)
	else:
		seqs = ((seqid, str(sequences[seqid].seq)) for seqid in seqids)

	records = []
	total_selected = 0
	exhausted = False
	while exhausted is False:
		record = next(seqs, None)
		if record is not None:
			record = fasta_str_record(ct.FASTA_RECORD_TEMPLATE,
									  [record[0], record[1]])
			records.append(record)
		else:
			exhausted = True

		# write records when it reaches the maximum number of records to
		# keep in memory or there are no records left to fetch
		if len(records) == limit or exhausted is True:
			fo.write_lines(records, output_file, write_mode='a')
			total_selected += len(records)
			records = []

	return total_selected


def split_seqcount(fasta_path, output_directory, max_seqs):
	"""Split a FASTA file based on a maximum number of sequences per file.

	Parameters
	----------
	fasta_path : str
		Path to a FASTA file.
	output_directory : str
		Path to the output directory.
	max_seqs : int
		Split FASTA file into files with a maximum number
		of sequences equal to this value.

	Returns
	-------
	split_files : list
		List with paths to the new files that were
		created by splitting the input FASTA file.
	"""
	file_count = 1
	exhausted = False
	current_recs = []
	split_files = []
	record_generator = sequence_generator(fasta_path)
	while exhausted is False:
		record = next(record_generator, None)
		if record is not None:
			current_recs.append(record)
		else:
			exhausted = True

		if len(current_recs) == max_seqs or exhausted is True:
			if len(current_recs) > 0:
				file_path = fo.join_paths(output_directory,
										  ['seqcount{0}.fasta'.format(file_count)])
				seqids = (rec.id for rec in current_recs)
				split_files.append([file_path, seqids])
				write_records(current_recs, file_path)
				current_recs = []
				file_count += 1

	return split_files


def split_seqlength(fasta_path, output_directory, length_cutoff):
	"""Split a FASTA file based on a sequence length threshold.

	Parameters
	----------
	fasta_path : str
		Path to a FASTA file.
	output_directory : str
		Path to the output directory.
	length_cutoff : int
		Sequence length threshold used to split the FASTA file.

	Returns
	-------
	List with two tuples: the first contains the path to a FASTA
	file with sequences equal or above the length cutoff and the
	list of sequence identifiers of those sequences or is None if
	there were no sequences above the cutoff. The second tuple
	contains the same data for the sequences below the cutoff or
	is None if there were no sequences below the cutoff.
	"""
	length_values = sequence_lengths(fasta_path)
	below_cutoff = [seqid for seqid, length in length_values.items()
					if length < length_cutoff]
	above_cutoff = list(set(length_values) - set(below_cutoff))

	fasta_index = index_fasta(fasta_path)

	above_data = [len(above_cutoff), above_cutoff]
	if len(above_cutoff) > 0:
		above_outfile = fo.join_paths(output_directory, ['above_cutoff.fasta'])
		above_count = get_sequences_by_id(fasta_index, above_cutoff, above_outfile)
		above_data.append(above_outfile)

	# File includes sequences shorter than cutoff value
	below_data = [len(below_cutoff), below_cutoff]
	if len(below_cutoff) > 0:
		below_outfile = fo.join_paths(output_directory, ['below_cutoff.fasta'])
		below_count = get_sequences_by_id(fasta_index, below_cutoff, below_outfile)
		below_data.append(below_outfile)

	return [above_data, below_data]


def fasta_stats(fasta_file):
	"""Determine the number of sequences in a FASTA file and length stats.

	Parameters
	----------
	fasta_file : str
		Path to a FASTA file.

	Returns
	-------
	fasta_file : str
		Path to the FASTA file.
	total_seqs: int
		Total number of records in the FASTA file.
	mean_length: float
		Mean sequence length.
	"""
	seq_lengths = sequence_lengths(fasta_file)
	min_length = min(seq_lengths.values())
	max_length = max(seq_lengths.values())
	mean_length = sum(seq_lengths.values())/len(seq_lengths)
	total_seqs = len(seq_lengths)

	return [fasta_file, total_seqs, min_length, max_length, mean_length]


def translate_fasta(input_fasta, output_directory, translation_table, write_dna=False):
	"""Translate DNA sequences in a FASTA file.

	Parameters
	----------
	input_fasta : str
		Path to the FASTA file that contains the DNA sequences
		to translate.
	output_directory : str
		Path to the output directory where the FASTA file with
		protein sequences will be written to.
	translation_table : int
		Genetic code used to translate DNA sequences.

	Returns
	-------
	input_fasta : str
		Path to the input FASTA file.
	protein_file : str
		Path to the FASTA file that contains the translated
		sequences.
	translated : int
		Number of sequences that were translated successfully.
	invalid : list
		List with one sublist for each sequence that could not
		be translated. Each sublist includes a sequence identifier
		and a exception message.
	"""
	# Import DNA sequences from the input FASTA file
	records = import_sequences(input_fasta)
	# Translate the DNA sequences
	translated_records = [[seqid,
						   sm.translate_dna(seq, translation_table, 0)]
						  for seqid, seq in records.items()]

	# Get the exceptions for the sequences that could not be translated
	invalid = [[rec[0], rec[1]]
			   for rec in translated_records
			   if type(rec[1]) == str]

	# Get the valid translated sequences
	valid_proteins = [[rec[0], str(rec[1][0][0])]
			 		  for rec in translated_records
			 		  if type(rec[1]) == list]

	translated = len(valid_proteins)

	# Do not attempt to create output FASTA files if none of the DNA sequences could be translated
	if len(valid_proteins) == 0:
		# Return NoneType for the output file paths
		return [None, None, translated, invalid]

	# Save protein sequences if there were any valid translations
	valid_protein_lines = fasta_lines(ct.FASTA_RECORD_TEMPLATE, valid_proteins)
	protein_file_basename = fo.file_basename(input_fasta, False) + '_protein.fasta'
	protein_file_path = fo.join_paths(output_directory, [protein_file_basename])
	fo.write_lines(valid_protein_lines, protein_file_path)

	# Write DNA sequences that could be translated to a FASTA file if requested
	# This is useful to get a FASTA file containing only the valid alleles
	if write_dna is True:
		# Identify DNA sequences that could be translated
		valid_dna = [[rec[0], str(rec[1][0][1])]
					 for rec in translated_records
					 if type(rec[1]) == list]
		valid_dna_lines = fasta_lines(ct.FASTA_RECORD_TEMPLATE, valid_dna)
		dna_file_basename = fo.file_basename(input_fasta)
		dna_file_path = fo.join_paths(output_directory, [dna_file_basename])
		fo.write_lines(valid_dna_lines, dna_file_path)
	else:
		dna_file_path = None

	return [dna_file_path, protein_file_path, translated, invalid]
