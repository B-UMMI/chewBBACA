#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the ComputeDistances module.
"""


from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  AfterValidator)

import CHEWBBACA.ComputeDistances.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
import CHEWBBACA.utils.file_operations as fo


# Annotated Fields
OutputDirectory = Annotated[str,
                            Field(description=ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

AllelicProfiles = Annotated[str,
                            Field(description=ct.COMPUTEDISTANCES_ALLELIC_PROFILES_DESCRIPTION),
                            AfterValidator(partial(gval.input_paths_exist,
                                                   function=fo.exists,
												   exception_str=gct.PROFILES_FILE_MISSING))]

DistanceComputationMethod = Annotated[str,
                   					  Field(default=ct.DEFAULT_DISTANCE_METHOD,
                                               description=ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_DESCRIPTION),
                   					  AfterValidator(partial(gval.validate_choice, valid_options=ct.DISTANCE_METHODS, exception_str=ct.INVALID_DISTANCE_METHOD))]

OutputFormat = Annotated[str,
                   Field(default=ct.DEFAULT_OUTPUT_FORMAT,
                         description=ct.COMPUTEDISTANCES_OUTPUT_FORMAT_DESCRIPTION),
				   AfterValidator(partial(gval.validate_choice, valid_options=ct.OUTPUT_FORMATS, exception_str=ct.INVALID_OUTPUT_FORMAT))]

ComputeSimilarity = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTEDISTANCES_COMPUTE_SIMILARITY_DESCRIPTION)]

NoMask = Annotated[bool,
                        Field(default=False,
                              description=ct.COMPUTEDISTANCES_NO_MASK_DESCRIPTION)]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.COMPUTEDISTANCES_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]


# Main validation class
class ComputeDistancesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: OutputDirectory
	allelic_profiles: AllelicProfiles
	distance_computation_method: DistanceComputationMethod
	output_format: OutputFormat
	compute_similarity: ComputeSimilarity
	no_mask: NoMask
	cpu_cores: CPUCores
