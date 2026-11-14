#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file includes the Pydantic classes used to valdiate the argument values passed to the PredictGenes module.
"""


import sys
import shutil
from typing import Annotated
from functools import partial

from pydantic import (BaseModel,
					  Field,
					  ConfigDict,
					  AfterValidator,
					  BeforeValidator,
					  field_validator,
					  model_validator)

import CHEWBBACA.PredictGenes.constants as ct
import CHEWBBACA.global_constants as gct
import CHEWBBACA.global_validation as gval
import CHEWBBACA.utils.file_operations as fo
import CHEWBBACA.utils.iterables_manipulation as im
import CHEWBBACA.utils.augustus_gene_prediction as agp


def check_augustus(augustus_path, augustus_alias=ct.AUGUSTUS_ALIAS):
	"""Determine if AUGUSTUS executable is available in the system.

	Parameters
	----------
	augustus_path : str
		Path to the folder containing AUGUSTUS executable.
	augustus_alias : str
		Alias used to call AUGUSTUS.

	Returns
	-------
	augustus_path : str
		Path to the folder containing AUGUSTUS executable.

	Raises
	------
	SystemExit
		- If AUGUSTUS executable is not available in the system.
	"""
	# User did not provide a path, try to determine path to the AUGUSTUS executable using its alias
	if not augustus_path:
		augustus_path = fo.get_parent_directory(shutil.which(augustus_alias))

	if augustus_path is None or fo.exists(augustus_path) is False:
		sys.exit(ct.AUGUSTUS_MISSING)

	return augustus_path


def validate_augustus_species(species_id, augustus_path, augustus_alias=ct.AUGUSTUS_ALIAS):
	"""Validate if the provided AUGUSTUS species is supported by the installed version of AUGUSTUS.

	Parameters
	----------
	species_id : str
		Species ID provided by the user.
	augustus_path : str
		Path to the folder containing AUGUSTUS executable.
	augustus_alias : str
		Alias used to call AUGUSTUS.

	Returns
	-------
	species_id : str
		Species ID that will be used for AUGUSTUS gene prediction.

	Raises
	------
	SystemExit
		- If the provided species ID is not supported by the installed version of AUGUSTUS.
	"""
	augustus_exe = fo.join_paths(augustus_path, [augustus_alias])
	species_list = agp.get_augustus_species_list(augustus_exe)

	if species_id not in species_list.values():
		# Create string with list of options
		species_list_str = [": ".join([k, v]) for k, v in species_list.items()]
		species_list_str = "\n".join(species_list_str)
		sys.exit(ct.AUGUSTUS_INVALID_SPECIES.format(species_list_str))

	return species_id


# Annotated Fields
OutputDirectory = Annotated[str,
                            Field(description=ct.PREDICTGENES_OUTPUT_DIRECTORY_DESCRIPTION),
							AfterValidator(gval.parentdir_exists),
							AfterValidator(gval.create_output_directory)]

# Need to include list type to avoid PydanticSerializationUnexpectedValue
# Input is a str but AfterValidator converts it into a list of valid file paths
InputPath = Annotated[str | list,
                      Field(description=ct.PREDICTGENES_INPUT_PATH_DESCRIPTION),
					   AfterValidator(partial(gval.input_paths_exist,
					                          function=fo.exists,
											  exception_str=gct.MISSING_INPUT_ARG)),
					   AfterValidator(gval.input_is_fasta),
					   AfterValidator(fo.list_input_files),
					   AfterValidator(gval.filter_inputs_extension),
					   AfterValidator(gval.inputs_exist),
					   AfterValidator(gval.validate_inputs_fastas),
					   AfterValidator(gval.check_unique_prefixes),
					   AfterValidator(gval.check_blanks)]

TranslationTable = Annotated[int,
							 Field(default=gct.GENETIC_CODE_DEFAULT,
							       description=ct.PREDICTGENES_TRANSLATION_TABLE_DESCRIPTION),
							 AfterValidator(partial(gval.validate_choice,
							                        valid_options=gct.GENETIC_CODES,
													exception_str=gct.INVALID_GENETIC_CODE))]

PyrodigalTrainingFile = Annotated[str | None,
                                  Field(default=None)]

PyrodigalMode = Annotated[str,
						  Field(default=ct.PYRODIGAL_DEFAULT_MODE),
						  AfterValidator(partial(gval.validate_choice, valid_options=ct.PYRODIGAL_MODES, exception_str=ct.PREDICTGENES_INVALID_PYRODIGAL_MODE))]

PyrodigalOutputFormats = Annotated[list,
								   Field(default=ct.PYRODIGAL_DEFAULT_OUTFMTS),
								   BeforeValidator(im.ensure_list),
								   AfterValidator(partial(gval.validate_multiple_choices, valid_options=ct.PYRODIGAL_OUTFMTS, exception_str=ct.PYRODIGAL_INVALID_OUTFMT))]

PyrodigalMinimumConfidence = Annotated[float | None,
                                       Field(default=None,
                                             ge=ct.PYRODIGAL_MIN_CONFIDENCE,
											 le=ct.PYRODIGAL_MAX_CONFIDENCE)]

PyrodigalTrainingReference = Annotated[str | None,
                                       Field(default=None)]

PyrodigalJustTraining = Annotated[bool,
                                  Field(default=False)]

AugustusPath = Annotated[str | None,
                         Field(default=None,
                               validate_default=True),
						 AfterValidator(check_augustus)]

AugustusSpecies = Annotated[str,
                            Field(default=None)]

AugustusOutputFormat = Annotated[list,
								 Field(default=ct.AUGUSTUS_OUTFMT_DEFAULT),
								 BeforeValidator(im.ensure_list),
								 AfterValidator(partial(gval.validate_multiple_choices, valid_options=ct.AUGUSTUS_OUTFMTS, exception_str=ct.AUGUSTUS_INVALID_OUTFMT))]

CPUCores = Annotated[int,
					 Field(default=gct.CPU_CORES_DEFAULT,
					       ge=gct.CPU_CORES_MIN,
						   le=gval.get_system_cpu_cores(),
						   description=ct.PREDICTGENES_CPU_CORES_DESCRIPTION),
					 AfterValidator(gval.verify_cpu_usage)]

# Use separate models to validate Pyrodigal and AUGUSTUS arguments, since they have different fields and validation requirements
class PyrodigalArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	pyrodigal_training_file: PyrodigalTrainingFile
	pyrodigal_mode: PyrodigalMode
	pyrodigal_output_formats: PyrodigalOutputFormats
	pyrodigal_minimum_confidence: PyrodigalMinimumConfidence
	pyrodigal_training_reference: PyrodigalTrainingReference
	pyrodigal_just_training: PyrodigalJustTraining

	@field_validator("pyrodigal_training_file", "pyrodigal_training_reference")
	@classmethod
	def ptf_exists(cls, value):
		if value is not None:
			return gval.input_paths_exist(value, fo.exists, ct.INVALID_PTF_PATH)

		return value


# Create str with the defautl values for the gene prediction parameter=argument pairs
# Added here to avoid circular imports
def create_default_gene_prediction_options():
	default_gene_prediction_config = PyrodigalArgs.model_validate({}).model_dump()
	default_gene_prediction_options = []
	for parameter, argument in default_gene_prediction_config.items():
		if type(argument) == list:
			default_gene_prediction_options.append(f"{parameter}={','.join(argument)}")
		else:
			default_gene_prediction_options.append(f"{parameter}={argument}")

	return default_gene_prediction_options


class AugustusArgs(BaseModel):
	# Enforce strict field checking
	model_config = ConfigDict(extra="forbid")

	augustus_path: AugustusPath
	augustus_species: AugustusSpecies
	augustus_output_format: AugustusOutputFormat

	@model_validator(mode="after")
	def validate_species(self):
		self.augustus_species = validate_augustus_species(self.augustus_species, self.augustus_path)

		return self


GenePredictor = Annotated[str,
						  Field(default=ct.GENE_PREDICTOR_DEFAULT,
						        description=ct.PREDICTGENES_GENE_PREDICTOR_DESCRIPTION),
						  AfterValidator(partial(gval.validate_choice, valid_options=ct.GENE_PREDICTORS, exception_str=ct.PREDICTGENES_INVALID_GENE_PREDICTOR))]

# Validate if input is str, not if it is None
# Need to include dict type to avoid PydanticSerializationUnexpectedValue
# Input is a list but AfterValidator converts it into a dictionary
GenePredictionOptions = Annotated[list | dict,
                                  Field(default=create_default_gene_prediction_options(),
                                        validate_default=True,
                                        description=ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_DESCRIPTION),
                                  AfterValidator(gval.parse_parameter_string)]

# Field definition for the validated gene prediction arguments
ValidatedGenePredictionOptions = Annotated[PyrodigalArgs | AugustusArgs | None, Field(default=None)]


# Main validation class
class PredictGenesValidator(BaseModel):
	# Use `from_attributes=True` to inspect the attributes of the 
	# argparse.Namespace object directly without having to convert to a dictionary
	# Do not use `extra="forbid"` so that it accepts the arg with the module name without having to remove it explicitly
	model_config = ConfigDict(from_attributes=True)

	# Field order is preserved when serializing with model_dump()
	output_directory: OutputDirectory
	input_path: InputPath
	gene_predictor: GenePredictor
	gene_prediction_options: GenePredictionOptions
	validated_gene_prediction_options: ValidatedGenePredictionOptions
	translation_table: TranslationTable
	cpu_cores: CPUCores

	# Further validation for gene prediction arguments
	@model_validator(mode="after")
	def validate_gene_prediction_arguments(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options = PyrodigalArgs(**self.gene_prediction_options)
		elif self.gene_predictor == "augustus":
			self.validated_gene_prediction_options = AugustusArgs(**self.gene_prediction_options)

		return self

	# Need to get genetic code from training file
	@model_validator(mode="after")
	def get_ptf_genetic_code(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_file is not None:
			self.translation_table = gval.add_ptf_genetic_code(self.validated_gene_prediction_options.pyrodigal_training_file, self.translation_table)

		return self

	@model_validator(mode="after")
	def mode_and_ptf(self):
		if self.gene_predictor == "pyrodigal":
			self.validated_gene_prediction_options.pyrodigal_mode = gval.check_meta(self.validated_gene_prediction_options.pyrodigal_mode, self.validated_gene_prediction_options.pyrodigal_training_file)

		return self

	@model_validator(mode="after")
	def ptf_and_reference(self):
		if self.gene_predictor == "pyrodigal" and self.validated_gene_prediction_options.pyrodigal_training_reference:
			# Check if user provided path to Pyrodigal training file and to a training reference, which is not allowed
			gval.check_ptf_tref_conflict(self.validated_gene_prediction_options.pyrodigal_training_file, self.validated_gene_prediction_options.pyrodigal_training_reference)
			# Create training file based on training reference
			print(f'Creating Pyrodigal training file based on {self.validated_gene_prediction_options.pyrodigal_training_reference}...')
			self.validated_gene_prediction_options.pyrodigal_training_file = gval.create_training_file(self.validated_gene_prediction_options.pyrodigal_training_reference, self.output_directory, self.translation_table)
			print(f'Training file saved to {self.validated_gene_prediction_options.pyrodigal_training_file}')
			# Exit if user only requested to create training file
			if self.validated_gene_prediction_options.pyrodigal_just_training:
				sys.exit(gct.JUST_TRAINING)

		return self
