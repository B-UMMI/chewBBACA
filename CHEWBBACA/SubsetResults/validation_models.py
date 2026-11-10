#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the SubsetResults module.
"""


import sys
from typing import Annotated
from functools import partial

import pandas as pd
from pydantic import (BaseModel,
                      Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

import CHEWBBACA.SubsetResults.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
OutputDirectory = Annotated[str,
                            Field(description=ct.SUBSETRESULTS_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

ResultsFiles = Annotated[str,
                         Field(description=ct.SUBSETRESULTS_RESULTS_FILES_DESCRIPTION),
                         AfterValidator(partial(gval.input_paths_exist,
                                                function=fo.exists,
												exception_str=gct.MISSING_RESULTS)),
						 AfterValidator(partial(gval.contains_results,
						                        results_files_basenames=[gct.ALLELIC_PROFILES_BASENAME]))]

LociList = Annotated[str,
                     Field(description=ct.SUBSETRESULTS_LOCI_LIST_DESCRIPTION),
					 AfterValidator(partial(gval.input_paths_exist,
					                        function=fo.exists,
											exception_str=gct.MISSING_LOCI_LIST)),
					 AfterValidator(gval.input_is_fasta),
					 AfterValidator(fo.list_input_files)]

SampleList = Annotated[str,
                       Field(description=ct.SUBSETRESULTS_SAMPLE_LIST_DESCRIPTION),
                       AfterValidator(partial(gval.input_paths_exist, function=fo.exists, exception_str=ct.MISSING_EXCLUDE_GENOMES)),
					   AfterValidator(gval.input_is_fasta),
					   AfterValidator(fo.list_input_files)]

InvertLoci = Annotated[bool,
                        Field(default=False,
                              description=ct.SUBSETRESULTS_INVERT_LOCI_DESCRIPTION)]

InvertSample = Annotated[bool,
                        Field(default=False,
                              description=ct.SUBSETRESULTS_INVERT_SAMPLES_DESCRIPTION)]


# Main validation class
class SubsetResultsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	results_files: ResultsFiles
	loci_list: LociList
	sample_list: SampleList
	invert_loci: InvertLoci
	invert_samples: InvertSample

	@model_validator(mode="after")
	def check_lists(self):
		# Do not proceed if user provided neither a list of loci nor of samples
		if not self.loci_list and not self.sample_list:
			sys.exit(ct.SUBSETRESULTS_MISSING_LISTS_EXCEPTION)

		return self

	@model_validator(mode="after")
	def check_files(self):
		# Do not proceed if the input folder does not include any of the supported output files
		input_files = fo.listdir(self.results_files)
		if len(set(gct.SUPPORTED_RESULTS_FILES)-set(input_files)) == len(gct.SUPPORTED_RESULTS_FILES):
			sys.exit("Input directory does not include any of the supported output files that can be subsetted."
			         "Please ensure that the input directory includes at least one of the following files:\n"
					 f"{'\n'.join(gct.SUPPORTED_RESULTS_FILES)}")

		return self

	# Check if the loci IDs in the input list are all in the results
	@model_validator(mode="after")
	def check_loci(self):
		file_types = [gct.ALLELIC_PROFILES_BASENAME, gct.COORDINATES_PROFILES_BASENAME, gct.PRESENCE_ABSENCE_PROFILES_BASENAME, gct.LOCI_STATS_BASENAME]
		if self.loci_list:
			for fname in file_types:
				# Determine if file exists in the input directory
				file_path = fo.join_paths(self.results_files, [fname])
				if fo.is_file(file_path)[0]:
					if fname in file_types[:3]:
						results_loci = fo.read_lines(file_path, strip=True, num_lines=1)
						results_loci = results_loci[0].split("\t")[1:]
					elif fname in file_types[3:]:
						results_loci = pd.read_csv(file_path, usecols=["Locus"], sep="\t", dtype=str)
						results_loci = results_loci["Locus"].tolist()

					non_matching = set(self.loci_list) - set(results_loci)
					if len(non_matching) > 0:
						sys.exit(f"The following loci IDs in the input loci list do not match "
							     f"loci IDs in the results files:\n{'\n'.join(non_matching)}\n"
								 "Please ensure that all loci IDs in the input loci list match "
								 "loci IDs in the resutls files.")
			# It is possible that the input directory does not include any of the files containing the loci IDs
			if not fo.is_file(file_path)[0]:
				print(f"Could not get the list of loci IDs from the results files in the input directory "
				      "to validate the loci IDs in the input loci list.")

		return self

	# Check if the sample IDs in the input list are all in the results
	@model_validator(mode="after")
	def check_samples(self):
		file_types = [gct.ALLELIC_PROFILES_BASENAME, gct.COORDINATES_PROFILES_BASENAME, gct.PRESENCE_ABSENCE_PROFILES_BASENAME, gct.SAMPLE_STATISTICS_BASENAME]
		if self.sample_list:
			for fname in file_types:
				file_path = fo.join_paths(self.results_files, [fname])
				if fo.is_file(file_path)[0]:
					results_samples = pd.read_csv(file_path, usecols=["FILE"], sep="\t", dtype=str)
					results_samples = results_samples["FILE"].tolist()
					non_matching = set(self.sample_list) - set(results_samples)
					if len(non_matching) > 0:
						sys.exit(f"The following sample IDs in the input sample list do not match "
							  	 f"sample IDs in the results files:\n{'\n'.join(non_matching)}\n"
							  	 "Please ensure that all sample IDs in the input sample list match "
							  	 "sample IDs in the resutls files.")
			if not fo.is_file(file_path)[0]:
				print(f"Could not get the list of sample IDs from the results files in the input directory "
				      "to validate the sample IDs in the input loci list.")

		return self
