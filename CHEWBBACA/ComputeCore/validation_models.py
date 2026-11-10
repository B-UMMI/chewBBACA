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
# Should implement a validator for the profile filr format
AllelicProfiles = Annotated[str,
                            Field(description=ct.COMPUTECORE_PROFILES_DESCRIPTION),
                            AfterValidator(partial(gval.input_paths_exist,
                                                   function=fo.exists,
												   exception_str=gct.PROFILES_FILE_MISSING))]

OutputDirectory = Annotated[str,
                            Field(description=ct.COMPUTECORE_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

LociPresenceThreshold = Annotated[list,
                            Field(default=ct.CGMLST_THRESHOLDS,
								  description=ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_DESCRIPTION),
							AfterValidator(partial(gval.validate_value_bounds,
								            	   lower_limit=ct.CGMLST_THRESHOLD_MIN,
												   upper_limit=ct.CGMLST_THRESHOLD_MAX,
												   exception=gct.LOCI_PRESENCE_VALUE_OUT_OF_LIMITS))]

SampleStep = Annotated[int,
                 Field(default=1,
                       description=ct.COMPUTECORE_SAMPLE_STEP_DESCRIPTION)]

ComputeAccessory = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTECORE_COMPUTE_ACCESSORY_DESCRIPTION)]

RarefactionAnalysis = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTECORE_RAREFACTION_ANALYSIS_DESCRIPTION)]

PermutationNumber = Annotated[int,
                        Field(default=ct.PERMUTATION_NUMBER_DEFAULT,
                              description=ct.COMPUTECORE_PERMUTATION_NUMBER_DESCRIPTION)]

PermutationSamples = Annotated[int,
                        Field(default=ct.PERMUTATION_SAMPLES_DEFAULT,
                              description=ct.COMPUTECORE_PERMUTATION_SAMPLES_DESCRIPTION)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.COMPUTECORE_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

ExcludeLoci = Annotated[str,
                        Field(default=None, description=ct.COMPUTECORE_EXCLUDE_LOCI_DESCRIPTION),
                        AfterValidator(partial(gval.input_paths_exist,
                                               function=fo.exists,
											   exception_str=ct.MISSING_EXCLUDE_LOCI)),
						AfterValidator(gval.input_is_fasta),
						AfterValidator(fo.list_input_files)]

ExcludeGenomes = Annotated[str,
                           Field(default=None, description=ct.COMPUTECORE_EXCLUDE_GENOMES_DESCRIPTION),
                           AfterValidator(partial(gval.input_paths_exist,
                                                  function=fo.exists,
												  exception_str=ct.MISSING_EXCLUDE_GENOMES)),
						   AfterValidator(gval.input_is_fasta),
						   AfterValidator(fo.list_input_files)]


# Main validation class
class ComputeCoreValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	# Need to validate the profiles (header starting with FILE followed by loci IDs, header line followed by profiles, etc?)
	profiles: AllelicProfiles
	loci_presence_threshold: LociPresenceThreshold
	sample_step: SampleStep
	compute_accessory: ComputeAccessory
	rarefaction_analysis: RarefactionAnalysis
	permutation_number: PermutationNumber
	permutation_samples: PermutationSamples
	cpu_cores: CPUCores
	exclude_loci: ExcludeLoci
	exclude_genomes: ExcludeGenomes

	@model_validator(mode="after")
	def check_exclude_loci(self):
		if self.exclude_loci:
			# Read list of loci in the profiles file
			results_loci = fo.read_lines(self.profiles, strip=True, num_lines=1)
			results_loci = results_loci[0].split("\t")[1:]
			extra_loci = set(self.exclude_loci) - set(results_loci)
			if len(extra_loci) > 0:
				print(ct.INVALID_EXCLUDE_LOCI.format("\n".join(extra_loci)))

		return self

	@model_validator(mode="after")
	def check_exclude_genomes(self):
		if self.exclude_genomes:
			# Read list of samples in the profiles file
			results_samples = fo.extract_column(self.profiles)
			extra_samples = set(self.exclude_genomes) - set(results_samples)
			if len(extra_samples) > 0:
				print(ct.INVALID_EXCLUDE_SAMPLES.format("\n".join(extra_samples)))

		return self
