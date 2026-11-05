#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the ComputeCore module.
"""


from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

import CHEWBBACA.ComputeCore.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
from CHEWBBACA.utils import file_operations as fo


# Annotated Fields
LociPresenceThreshold = Annotated[list,
                            Field(default=ct.CGMLST_THRESHOLDS,
                                  ge=ct.CGMLST_THRESHOLD_MIN,
								  le=ct.CGMLST_THRESHOLD_MAX)]

SampleStep = Annotated[int,
                 Field(default=1)]

Permutation = Annotated[int,
                        Field(default=ct.PERMUTATION_NUMBER_DEFAULT)]

ExcludeLoci = Annotated[str,
                        AfterValidator(partial(gval.input_paths_exist,
                                               function=fo.exists,
											   exception_str=ct.MISSING_EXCLUDE_LOCI)),
						AfterValidator(gval.input_is_fasta),
						AfterValidator(fo.list_input_files)]

ExcludeGenomes = Annotated[str,
                           AfterValidator(partial(gval.input_paths_exist,
                                                  function=fo.exists,
												  exception_str=ct.MISSING_EXCLUDE_GENOMES)),
						   AfterValidator(gval.input_is_fasta),
						   AfterValidator(fo.list_input_files)]


# Main validation class
class ComputeCoreValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	# Need to validate the profiles (header starting with FILE followed by loci IDs, header line followed by profiles, etc?)
	profiles: gval.AllelicProfiles
	loci_presence_threshold: LociPresenceThreshold
	sample_step: SampleStep
	compute_accessory: gval.BooleanFlag
	rarefaction_analysis: gval.BooleanFlag
	permutation_number: Permutation
	permutation_samples: int | None = None
	cpu_cores: gval.CPUCores
	exclude_loci: ExcludeLoci | None = None
	exclude_genomes: ExcludeGenomes | None = None

	@model_validator(mode="after")
	def check_exclude_loci(self):
		if self.exclude_loci:
			# Read list of loci in the profiles file
			results_loci = fo.read_lines(self.profiles, strip=True, num_lines=1)
			results_loci = results_loci[0].split("\t")[1:]
			extra_loci = set(self.exclude_loci) - set(results_loci)
			if len(extra_loci) > 0:
				print(gct.INVALID_EXCLUDE_LOCI.format("\n".join(extra_loci)))

		return self

	@model_validator(mode="after")
	def check_exclude_genomes(self):
		if self.exclude_genomes:
			# Read list of samples in the profiles file
			results_samples = fo.extract_column(self.profiles)
			extra_samples = set(self.exclude_genomes) - set(results_samples)
			if len(extra_samples) > 0:
				print(gct.INVALID_EXCLUDE_SAMPLES.format("\n".join(extra_samples)))

		return self
