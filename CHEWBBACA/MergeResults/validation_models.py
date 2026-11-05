#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


import sys
from typing import Annotated
from functools import partial


import pandas as pd
from pydantic import (BaseModel,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

try:
	import global_constants as gct
	import global_validation as gval
	import MergeResults.constants as ct
	from utils import file_operations as fo
except ModuleNotFoundError:
	import CHEWBBACA.global_constants as gct
	import CHEWBBACA.global_validation as gval
	import CHEWBBACA.MergeResults.constants as ct
	from CHEWBBACA.utils import file_operations as fo


ResultsDirectories = Annotated[list,
								AfterValidator(partial(gval.input_paths_exist,
								                       function=fo.is_directory,
													   exception_str=ct.MERGERESULTS_INPUTFILE_EXCEPTION))]


class MergeResultsValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	results_directories: ResultsDirectories
	common: gval.BooleanFlag

	# Determine if user passed more than one input path
	@model_validator(mode="after")
	def check_multiple_input_folders(self):
		if len(self.results_directories) <= 1:
			sys.exit("Provided a single input path. Need multiple input paths to merge results files from different folders.")

		return self

	# Check if results have loci in common
	@model_validator(mode="after")
	def validate_loci(self):
		file_types = [gct.ALLELIC_PROFILES_BASENAME, gct.COORDINATES_PROFILES_BASENAME, gct.PRESENCE_ABSENCE_PROFILES_BASENAME, gct.LOCI_STATS_BASENAME]
		results_files = {fname: [] for fname in file_types}
		for folder in self.results_directories:
			# List files in folder
			current_files = fo.listdir(folder)
			for fname, files in results_files.items():
				if fname in current_files:
					results_files[fname].append(fo.join_paths(folder, [fname]))

		for fname, files in results_files.items():
			if len(files) != len(self.results_directories):
				print(f"Cannot determine the set of shared loci based on the {fname} "
				      "file type because it is not present in all input directories.")
			# All input directories contain the file type
			# Can use the file type to check if results have common loci
			else:
				headers = []
				for f in files:
					if fname in file_types[:3]:
						header = fo.read_lines(f, strip=True, num_lines=1)
						results_loci = header[0].split("\t")[1:]
					elif fname in file_types[3:]:
						results_loci = pd.read_csv(f, usecols=["Locus"], sep="\t", dtype=str)
						results_loci = results_loci["Locus"].tolist()
					headers.append(results_loci)
	
				loci = set(headers[0])
				for h in headers[1:]:
					# Determine common loci ordered based on headers in first file
					loci = loci.intersection(set(h))

				# Sort list of common loci
				loci = sorted(list(loci))

				# Exit if there are no common loci between all results
				if len(loci) == 0:
					sys.exit(ct.MERGERESULTS_NOCOMMONLOCI_EXCEPTION)
				# Exit if results are not for the same set of loci
				elif len(loci) != len(headers[0]) and not self.common:
					sys.exit(ct.MERGERESULTS_LOCIDIFFER_EXCEPTION)

				self.common = [self.common, loci, len(headers[0])]

				# No need to check based on other file type
				break

			# Exit if it was not possible to determine if results have loci IDs in common
			if type(self.common) == bool:
				sys.exit("Could not determine if the results files included in the input "
				         "directories share loci IDs. Will not attempt to merge the results."
						 "Please make sure that all input directories include at least one "
						 f"of the following file types: {'\n'.join(file_types)}")

		return self
