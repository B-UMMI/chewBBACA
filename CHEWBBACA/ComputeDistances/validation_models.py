#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------


"""


from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  model_validator)

try:
	import ComputeDistances.constants as ct
	import global_validation as gval
except ModuleNotFoundError:
	import CHEWBBACA.ComputeDistances.constants as ct
	import CHEWBBACA.global_validation as gval


DistanceComputationMethod = Annotated[str,
                   					  Field(default=ct.DEFAULT_DISTANCE_METHOD),
                   					  AfterValidator(partial(gval.validate_choice, valid_options=ct.DISTANCE_METHODS, exception_str=ct.INVALID_DISTANCE_METHOD))]

OutputFormat = Annotated[str,
                   Field(default=ct.DEFAULT_OUTPUT_FORMAT),
				   AfterValidator(partial(gval.validate_choice, valid_options=ct.OUTPUT_FORMATS, exception_str=ct.INVALID_OUTPUT_FORMAT))]


class ComputeDistancesValidator(BaseModel):
	model_config = ConfigDict(from_attributes=True)

	output_directory: gval.OutputDirectory
	allelic_profiles: gval.AllelicProfiles
	distance_computation_method: DistanceComputationMethod
	output_format: OutputFormat
	compute_similarity: gval.BooleanFlag
	no_mask: gval.BooleanFlag
	cpu_cores: gval.CPUCores
