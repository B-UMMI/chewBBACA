#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This is the main script of the chewBBACA suite. It parses the options and 
arguments provided through the command line and calls the specified module.
"""


import sys
import argparse

try:
	import global_constants as gct
	import global_validation as gval
	from utils.process_datetime import process_timer
except ModuleNotFoundError:
	from CHEWBBACA import global_constants as gct
	from CHEWBBACA import global_validation as gval
	from CHEWBBACA.utils.process_datetime import process_timer


class ModifiedHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):

	# prog is the name of the program 'ex: chewBBACA.py'
	def __init__(self, prog, indent_increment=1, max_help_position=56, width=100):
		super().__init__(prog, indent_increment, max_help_position, width)

	# Override split lines method
	def _split_lines(self, text, width):
		lines = super()._split_lines(text, width) + ['']
		return lines

	def _format_action_invocation(self, action):
		if not action.option_strings:
			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			return metavar
		else:
			parts = []
			parts.extend(action.option_strings)
			parts_text = ', '.join(parts)
			return f'{parts_text}'


@process_timer
def run_predict_genes():
	"""Run the PredictGenes module to predict genes from a set of input genome assemblies."""

	try:
		from PredictGenes import constants as ct
		from PredictGenes import validation_models as vm
		from PredictGenes import predict_genes
	except ModuleNotFoundError:
		from CHEWBBACA.PredictGenes import constants as ct
		from CHEWBBACA.PredictGenes import validation_models as vm
		from CHEWBBACA.PredictGenes import predict_genes

	def msg(name=None):
		usage_msg = "chewBBACA.py PredictGenes --input-files <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="PredictGenes",
									 description="Predict genes from a set of input genome assemblies.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.PREDICTGENES_MODULE_NAME)}")

	parser.add_argument('PredictGenes', nargs='+', help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-files", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_INPUT_FILES_ARGNAME],
						required=False, dest=ct.PREDICTGENES_INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the input files or to a file "
							 "with a list of full paths to the input files, one per line. Input "
							 "files must be in FASTA.")

	parser.add_argument("-o", "--output-directory", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.PREDICTGENES_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store the "
							 "files with the predicted CDSs.")

	parser.add_argument("--gp", "--gene-predictor", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_GENE_PREDICTOR_ARGNAME],
					 	required=False, dest=ct.PREDICTGENES_GENE_PREDICTOR_ARGNAME,
						help="Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")

	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME],
					    nargs="+", required=False, dest=ct.PREDICTGENES_GENE_PREDICTION_OPTIONS_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When providing "
							 "genome assemblies in FASTA format, the list of arguments for each parameter "
							 "used to configure the gene prediction can be passed as the long format of "
							 "the parameter name followed by the argument value (e.g., pyrodigal-training"
							 "-file=/path/to/file).")

	parser.add_argument("--tt", "--translation-table", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_TRANSLATION_TABLE_ARGNAME],
						required=False, dest=ct.PREDICTGENES_TRANSLATION_TABLE_ARGNAME,
						help="The translation table, matching a specific genetic code, used for gene "
							 "prediction. This value is ignored if a Pyrodigal training file is passed to "
							 "`-gpo`, `--gene-prediction-options`. In that case, the translation table "
							 "specified in the Pyrodigal training file is used.")

	parser.add_argument("--cpu", "--cpu-cores", type=ct.PREDICTGENES_ARGUMENT_TYPES[ct.PREDICTGENES_CPU_CORES_ARGNAME],
						required=False, dest=ct.PREDICTGENES_CPU_CORES_ARGNAME,
						help="Number of CPU cores that will be used to run the process (chewie resets "
							 "to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores).")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.PredictGenesValidator.model_validate(clean_args)

	# Predict CDSs
	predict_genes.main(validated_args)


@process_timer
def run_create_schema():
	"""Run the CreateSchema module to create a schema seed."""

	try:
		from CreateSchema import constants as ct
		from CreateSchema import validation_models as vm
		from CreateSchema import create_schema
	except ModuleNotFoundError:
		from CHEWBBACA.CreateSchema import constants as ct
		from CHEWBBACA.CreateSchema import validation_models as vm
		from CHEWBBACA.CreateSchema import create_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py CreateSchema --input-files <path> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="CreateSchema",
									 description="Create a schema seed.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog="It is strongly advised to provide a training file to create a schema. "
											f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.CREATESCHEMA_MODULE_NAME)}")

	parser.add_argument("CreateSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-files", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_INPUT_FILES_ARGNAME],
						required=True, dest=ct.CREATESCHEMA_INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the input FASTA files or to a file"
							 " with a list of full paths to FASTA files, one per line.")

	parser.add_argument("-o", "--output-directory", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.CREATESCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help="Output directory where the process will store intermediate files and "
							 "create the schema's directory.")

	parser.add_argument("--sn", "--schema-name", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_SCHEMA_NAME_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_SCHEMA_NAME_ARGNAME,
						help="Name given to the schema folder.")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp "
							 "alignment and aligned sequences with a BSR >= than the defined value "
							 "are considered to be alleles of the same gene.")

	parser.add_argument("--ml", "--minimum-length", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value. Predicted coding sequences (CDSs) "
							 "shorter than this value are excluded.")

	parser.add_argument("--tt", "--translation-table", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_TRANSLATION_TABLE_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_TRANSLATION_TABLE_ARGNAME,
						help="The translation table, matching a specific genetic code, used for gene "
							 "prediction and to translate CDSs. This value is ignored if a valid Pyrodigal "
							 "training file is passed to `-gpo`, `--gene-prediction-options`. In that case, "
							 "the translation table specified in the Pyrodigal training file is used.")

	parser.add_argument("--st", "--size-threshold", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_SIZE_THRESHOLD_ARGNAME,
						help="CDS size variation threshold. This value is used to define how much allele "
							 "size can vary for each locus (e.g. setting this value to 0.2 will define a "
							 "threshold that limits allele size variation to +/-20%% the allele size mode "
							 "for each locus). This value is not used in the CreateSchema process, but is "
							 "added to the schema's config file to be used during allele calling to check "
							 "if the sequence size of the potential new alleles is within the defined size "
							 "threshold.")

	parser.add_argument("--gp", "--gene-predictor", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_GENE_PREDICTOR_ARGNAME],
					 	required=False, dest=ct.CREATESCHEMA_GENE_PREDICTOR_ARGNAME,
						help="Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")

	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME],
					    nargs="+", required=False, dest=ct.CREATESCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When "
							 "providing genome assemblies in FASTA format, the list of arguments "
							 "for each parameter used to configure the gene prediction can be "
							 "passed as the long format of the parameter name followed by the "
							 "argument value (e.g., pyrodigal-training-file=/path/to/file).")

	parser.add_argument("--co", "--clustering-options", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME],
					 	nargs="+", required=False, dest=ct.CREATESCHEMA_CLUSTERING_OPTIONS_ARGNAME,
						help="List of arguments passed to configure the clustering of predicted "
							 "coding sequences (CDSs). The list of arguments for each parameter "
							 "used to configure the clustering can be passed as the long format "
							 "of the parameter name followed by the argument value (e.g., word-size=).")

	parser.add_argument("--bp", "--blast-path", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_BLAST_PATH_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--ci", "--cds-input", action="store_true",
					 	required=False, dest=ct.CREATESCHEMA_CDS_INPUT_ARGNAME,
						help="If provided, chewBBACA skips the gene prediction step and "
							 "assumes the input FASTA files contain coding sequences.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.CREATESCHEMA_ARGUMENT_TYPES[ct.CREATESCHEMA_CPU_CORES_ARGNAME],
						required=False, dest=ct.CREATESCHEMA_CPU_CORES_ARGNAME,
						help="Number of CPU cores that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores).")

	parser.add_argument("--no-cds-renaming", action="store_true",
					 	required=False, dest=ct.CREATESCHEMA_NO_CDS_RENAMING_ARGNAME,
						help="Do not rename the sequence/CDS identifiers when using the `--ci` "
							 "option. Provide this parameter when the input FASTA files containing "
							 "CDSs were generated by the PredictGenes module or if you are sure that "
							 "the CDS identifiers conform to the format used by chewBBACA (i.e. the "
							 "input file basename and an integer joined by `_`).")

	parser.add_argument("--no-cleanup", action='store_true',
					 	required=False, dest=ct.CREATESCHEMA_NO_CLEANUP_ARGNAME,
						help="If provided, intermediate files generated during process execution "
							 "are not deleted at the end.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.CreateSchemaValidator.model_validate(clean_args)

	# Run the CreateSchema process
	nloci = create_schema.main(validated_args)
	print(f'Created schema seed with {nloci} loci.')


@process_timer
def run_allele_call():
	"""Run the AlleleCall module to perform allele calling."""

	try:
		from AlleleCall import constants as ct
		from AlleleCall import validation_models as vm
		from AlleleCall import allele_call
	except ModuleNotFoundError:
		from CHEWBBACA.AlleleCall import constants as ct
		from CHEWBBACA.AlleleCall import validation_models as vm
		from CHEWBBACA.AlleleCall import allele_call

	def msg(name=None):
		usage_msg = "chewBBACA.py AlleleCall --input-files <path> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="AlleleCall",
									 description="Determine the allelic profiles of a set of genomes.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog="It is strongly advised to perform allele calling with the schema's "
									 		"parameters to ensure the consistency of the results. "
											f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.ALLELECALL_MODULE_NAME)}")

	parser.add_argument("AlleleCall", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-files", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_INPUT_FILES_ARGNAME],
						required=True, dest=ct.ALLELECALL_INPUT_FILES_ARGNAME,
						help="Path to the directory that contains the input FASTA files or to "
							 "a file with a list of full paths to FASTA files, one per line.")

	parser.add_argument("-s", "--schema-directory", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.ALLELECALL_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema directory. The schema directory contains the "
							 "loci FASTA files and a folder named `short` that contains the FASTA "
							 "files with the loci representative alleles.")

	parser.add_argument("-o", "--output-directory", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.ALLELECALL_OUTPUT_DIRECTORY_ARGNAME,
						help="Output directory where the process will store intermediate files "
							 "and allele calling results (will create a subdirectory named "
							 "`results_<TIMESTAMP>` if the path passed by the user already exists).")

	parser.add_argument("--ll", "--loci-list", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_LOCI_LIST_ARGNAME],
						required=False, dest=ct.ALLELECALL_LOCI_LIST_ARGNAME,
						help="Path to a file with the list of schema loci to perform allele "
							 "calling for. The file must include the full paths to the loci FASTA "
							 "files or the loci IDs, one per line. The process will perform "
							 "allele calling only for the subset of loci provided in the file.")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_BLAST_SCORE_RATIO_ARGNAME],
						required=False, dest=ct.ALLELECALL_BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp "
							 "alignment and aligned sequences with a BSR >= than the defined value "
							 "are considered to be alleles of the same locus.")

	parser.add_argument("--ml", "--minimum-length", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.ALLELECALL_MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value. Predicted coding sequences (CDSs) "
							 "shorter than this value are excluded.")

	parser.add_argument("--tt", "--translation-table", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_TRANSLATION_TABLE_ARGNAME],
						required=False, dest=ct.ALLELECALL_TRANSLATION_TABLE_ARGNAME,
						help="Genetic code used to predict genes and to translate coding DNA "
							 "sequences (CDSs). This value will be ignored if a training file "
							 "is used.")

	parser.add_argument("--st", "--size-threshold", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.ALLELECALL_SIZE_THRESHOLD_ARGNAME,
						help="Coding sequence (CDS) size variation threshold. At the default "
							 "value of 0.2, CDSs with a size that deviates +-20 percent from "
							 "the locus length mode are classified as ASM/ALM.")

	parser.add_argument("-gp", "--gene-predictor", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_GENE_PREDICTOR_ARGNAME],
					 	required=False, dest=ct.ALLELECALL_GENE_PREDICTOR_ARGNAME,
						help="Specify which gene prediction software to use. Default is Pyrodigal "
							 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
							 "for prokaryotic and eukaryotic genomes.")	

	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME],
					 	nargs="+", required=False, dest=ct.ALLELECALL_GENE_PREDICTION_OPTIONS_ARGNAME,
						help="List of options passed to configure the gene prediction. When "
							 "providing genome assemblies in FASTA format, the list of arguments "
							 "for each parameter used to configure the gene prediction can be "
							 "passed as the long format of the parameter name followed by the "
							 "argument value (e.g., pyrodigal-training-file=/path/to/file).")

	parser.add_argument("--co", "--clustering-options", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_CLUSTERING_OPTIONS_ARGNAME],
					 	nargs="+", required=False, dest=ct.ALLELECALL_CLUSTERING_OPTIONS_ARGNAME,
						help="List of arguments passed to configure the clustering of predicted "
							 "coding sequences (CDSs). The list of arguments for each parameter "
							 "used to configure the clustering can be passed as the long format "
							 "of the parameter name followed by the argument value (e.g., word-size=).")

	parser.add_argument("--bp", "--blast-path", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_BLAST_PATH_ARGNAME],
						required=False, dest=ct.ALLELECALL_BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--ci", "--cds-input", action='store_true',
						required=False, dest=ct.ALLELECALL_CDS_INPUT_ARGNAME,
						help="If provided, chewBBACA skips the gene prediction step and "
							 "assumes the input FASTA files contain coding sequences "
							 "(one FASTA file per strain).")

	parser.add_argument("--ni", "--no-inferred", action="store_true",
						required=False, dest=ct.ALLELECALL_NO_INFERRED_ARGNAME,
						help="If provided, the process will not add the sequences of "
							 "inferred alleles (INF) to the schema. Allelic profiles "
							 "will still include the allele identifiers attributed to "
							 "the inferred alleles. Use this parameter if the schema "
							 "is being accessed by multiple processes/users simultaneously.")

	parser.add_argument("--oo", "--output-unclassified", action="store_true",
						required=False, dest=ct.ALLELECALL_OUTPUT_UNCLASSIFIED_ARGNAME,
						help="Create a Fasta file with the coding sequences (CDSs) that "
							 "were not classified.")

	parser.add_argument("--om", "--output-missing", action="store_true",
					 	required=False, dest=ct.ALLELECALL_OUTPUT_MISSING_ARGNAME,
						help="Create a Fasta file with coding sequences (CDSs) classified "
							 "as NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5 and LOTSC.")

	parser.add_argument("--on", "--output-novel", action='store_true',
						required=False, dest=ct.ALLELECALL_OUTPUT_NOVEL_ARGNAME,
						help="Create a Fasta file with the novel alleles inferred during "
							 "allele calling. The sequence headers include the locus and "
							 "allele identifiers attributed by chewBBACA based on the "
							 "allele calling results.")

	parser.add_argument("--omp", "--output-masked-profiles", action="store_true",
						required=False, dest=ct.ALLELECALL_OUTPUT_MASKED_ARGNAME,
						help="Create a TSV file with the masked allelic profiles. The "
							 "masking process removes the `INF-` prefix from inferred "
							 "alleles and substitutes all special classes (NIPH, NIPHEM, "
							 "ASM, ALM, PLOT3, PLOT5, LOTSC, PAMA) with `0`.")

	parser.add_argument("--fc", "--force-continue", action='store_true',
						required=False, dest=ct.ALLELECALL_FORCE_CONTINUE_ARGNAME,
						help="If provided, chewie will not warn users and ask for "
							 "permission to continue if any of the provided argument "
							 "values does not match the values in the config file.")

	parser.add_argument("--em", "--execution-mode", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_EXECUTION_MODE_ARGNAME],
					 	required=False, dest=ct.ALLELECALL_EXECUTION_MODE_ARGNAME,
						help="Execution mode (1: only exact matches at DNA level; 2: "
							 "exact matches at DNA and Protein level; 3: exact matches "
							 "and minimizer-based clustering to find similar alleles "
							 "based on BSR+0.1; 4: run the full process to find exact "
							 "matches and similar matches based on BSR value, including "
							 "the determination of new representative alleles to add to "
							 "the schema).")

	parser.add_argument("--cc", "--cpu-cores", type=ct.ALLELECALL_ARGUMENT_TYPES[ct.ALLELECALL_CPU_CORES_ARGNAME],
						required=False, dest=ct.ALLELECALL_CPU_CORES_ARGNAME,
						help="Number of CPU cores that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores).")

	parser.add_argument('--no-cds-renaming', action='store_true',
					 	required=False, dest=ct.ALLELECALL_NO_CDS_RENAMING_ARGNAME,
						help="Do not rename the sequence/CDS identifiers when using "
							 "the `--cds` option. Provide this parameter when the input "
							 "FASTA files containing CDSs were generated by the PredictGenes "
							 "module or if you are sure that the CDS identifiers conform "
							 "to the format used by chewBBACA (the input file basename "
							 "and an integer joined by `_`).")

	parser.add_argument("--no-cleanup", action='store_true',
						required=False, dest=ct.ALLELECALL_NO_CLEANUP_ARGNAME,
						help="If provided, intermediate files generated during process execution "
							 "are not removed at the end.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.AlleleCallValidator.model_validate(clean_args)

	# Run the AlleleCall process
	allele_call.main(validated_args)


@process_timer
def run_evaluate_schema():
	"""Run the SchemaEvaluator module to evaluate a typing schema."""

	try:
		from SchemaEvaluator import constants as ct
		from SchemaEvaluator import validation_models as vm
		from SchemaEvaluator import evaluate_schema
	except ModuleNotFoundError:
		from CHEWBBACA.SchemaEvaluator import constants as ct
		from CHEWBBACA.SchemaEvaluator import validation_models as vm
		from CHEWBBACA.SchemaEvaluator import evaluate_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py SchemaEvaluator --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="SchemaEvaluator",
									 description="Build an interactive report for schema evaluation.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog="The module can evaluate schemas created with chewBBACA or from external platforms. "
									 		f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.SCHEMAEVALUATOR_MODULE_NAME)}")

	parser.add_argument("SchemaEvaluator", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-s", "--schema-directory", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory.")

	parser.add_argument("-o", "--output-directory", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the report HTML files will be created.")

	parser.add_argument("--ll", "--loci-list", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_LOCI_LIST_ARGNAME],
						required=False, dest=ct.SCHEMAEVALUATOR_LOCI_LIST_ARGNAME,
						help="Path to a file with the list of loci in the schema that the process "
							 "should analyse (one per line, full paths or loci IDs).")

	parser.add_argument("--af", "--annotations-file", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME],
					 	required=False, dest=ct.SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME,
						help="Path to the TSV file with loci annotations to include in a datatable "
							 "component of the report.")

	parser.add_argument("--tt", "--translation-table", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME],
						required=False, dest=ct.SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME,
						help="Translation table used to translate CDSs.")

	parser.add_argument("--st", "--size-threshold", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME,
						help="Coding sequence (CDS) size variation threshold. The module identifies "
							 "the alleles with size that deviates from the locus length mode +/- this "
							 "size threshold value.")

	parser.add_argument("--ml", "--minimum-length", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value. The module identifies alleles shorter "
							 "than this value.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.SCHEMAEVALUATOR_ARGUMENT_TYPES[ct.SCHEMAEVALUATOR_CPU_CORES_ARGNAME],
						required=False, dest=ct.SCHEMAEVALUATOR_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--lr", "--loci-reports", action="store_true",
					 	required=False, dest=ct.SCHEMAEVALUATOR_LOCI_REPORTS_ARGNAME,
						help="Create a detailed report page for each locus. The locus report "
							 "includes components with relevant data and analysis results, such "
							 "as allele diversity charts, a MSA for the alignment of the distinct "
							 "translated alleles and a tree drawn with Phylocanvas based on the "
							 "MAFFT guide tree.")

	parser.add_argument("--light", action="store_true",
					 	required=False, dest=ct.SCHEMAEVALUATOR_LIGTH_ARGNAME,
						help="Skips MSA computation with MAFFT and does not add the Phylogenetic "
							 "Tree and MSA components to the loci reports.")

	parser.add_argument("--add-sequences", action="store_true",
					 	required=False, dest=ct.SCHEMAEVALUATOR_ADD_SEQUENCES_ARGNAME,
						help="Adds Code Editor components with the DNA and Protein sequences to "
							 "the loci reports. The Code Editor is in readonly mode (allows to "
							 "search for and copy text).")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.SchemaEvaluatorValidator.model_validate(clean_args)

	evaluate_schema.main(validated_args)


@process_timer
def run_evaluate_calls():
	"""Run the AlleleCallEvaluator module to evaluate allele calling results."""

	try:
		from AlleleCallEvaluator import constants as ct
		from AlleleCallEvaluator import validation_models as vm
		from AlleleCallEvaluator import evaluate_calls
	except ModuleNotFoundError:
		from CHEWBBACA.AlleleCallEvaluator import constants as ct
		from CHEWBBACA.AlleleCallEvaluator import validation_models as vm
		from CHEWBBACA.AlleleCallEvaluator import evaluate_calls

	def msg(name=None):
		usage_msg = "chewBBACA.py AlleleCallEvaluator --input-files <dir> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="AlleleCallEvaluator",
									 description="Build an interactive report for allele calling results evaluation.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.ALLELECALLEVALUATOR_MODULE_NAME)}")

	parser.add_argument("AlleleCallEvaluator", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-files", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME],
					 	required=True, dest=ct.ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME,
						help="Path to the directory that contains the allele calling results "
							 "generated by the AlleleCall module.")

	parser.add_argument("-s", "--schema-directory", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory.")

	parser.add_argument("-o", "--output-directory", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME],
					 	required=True, dest=ct.ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the module will store intermediate "
							 "files and create the report HTML files.")

	parser.add_argument("--af", "--annotations-file", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME],
					 	required=False, dest=ct.ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME,
						help="Path to the TSV file with loci annotations to include in a datatable "
							 "component of the report.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.ALLELECALLEVALUATOR_ARGUMENT_TYPES[ct.ALLELECALLEVALUATOR_CPU_CORES_ARGNAME],
						required=False, dest=ct.ALLELECALLEVALUATOR_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--npa", "--no-presence-absence", action="store_true",
					 	required=False, dest=ct.ALLELECALLEVALUATOR_NO_PRESENCE_ABSENCE_ARGNAME,
						help="Do not compute the presence-absence matrix.")

	parser.add_argument("--ndm", "--no-distance-matrix", action='store_true',
					 	required=False, dest=ct.ALLELECALLEVALUATOR_NO_DISTANCE_MATRIX_ARGNAME,
						help="Do not compute the distance matrix.")

	parser.add_argument("--nnj", "--no-neighbor-joining", action="store_true",
					 	required=False, dest=ct.ALLELECALLEVALUATOR_NO_NEIGHBOR_JOINING_ARGNAME,
						help="Do not compute the Neighbor-Joining tree.")

	parser.add_argument("--fcm", "--force-core-msa", action="store_true",
					 	required=False, dest=ct.ALLELECALLEVALUATOR_FORCE_CORE_MSA_ARGNAME,
						help="Compute the MSA of the core loci, even if `--nnj` is provided.")

	parser.add_argument("--light", action="store_true",
					 	required=False, dest=ct.ALLELECALLEVALUATOR_LIGTH_ARGNAME,
						help="Do not compute the presence-absence matrix, the distance matrix, and "
							 "the Neighbor-Joining tree.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.AlleleCallEvaluatorValidator.model_validate(clean_args)

	evaluate_calls.main(validated_args)


@process_timer
def run_compute_core():
	"""Run the ComputeCore module to determine the set of core loci based on allele calling results."""

	try:
		from ComputeCore import constants as ct
		from ComputeCore import validation_models as vm
		from ComputeCore import compute_core
	except ModuleNotFoundError:
		from CHEWBBACA.ComputeCore import constants as ct
		from CHEWBBACA.ComputeCore import validation_models as vm
		from CHEWBBACA.ComputeCore import compute_core

	def msg(name=None):
		usage_msg = "chewBBACA.py ComputeCore --input-file <file> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="ComputeCore",
									 description="Determine the set of core loci based on allele calling results.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.COMPUTECORE_MODULE_NAME)}")

	parser.add_argument("ComputeCore", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-p", "--profiles", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_PROFILES_ARGNAME],
						required=True, dest=ct.COMPUTECORE_PROFILES_ARGNAME,
						help="Path to the TSV file that contains the allelic profiles determined by "
							 "the AlleleCall module.")

	parser.add_argument("-o", "--output-directory", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.COMPUTECORE_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the directory where the process will store the output files.")

	parser.add_argument("--lpt", "--loci-presence-threshold", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME], nargs="+",
						required=False, dest=ct.COMPUTECORE_LOCI_PRESENCE_THRESHOLD_ARGNAME,
						help="Loci/genes that constitute the core genome must be in a proportion of "
							 "genomes that is at least equal to this value. Provide multiple values "
							 "to compute the core genome for multiple threshold values.")

	parser.add_argument("--ss", "--sample-step", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_SAMPLE_STEP_ARGNAME],
					 	required=False, dest=ct.COMPUTECORE_SAMPLE_STEP_ARGNAME,
						help="The allele calling results are processed iteratively to evaluate the "
							 "impact of adding subsets of the results in computing the core genome. "
							 "The step value controls the number of allelic profiles added in each "
							 "iteration until all profiles are included.")

	parser.add_argument("--ca", "--compute-accessory", action="store_true",
					 	required=False, dest=ct.COMPUTECORE_COMPUTE_ACCESSORY_ARGNAME,
						help="Determine the set of accessory loci. The accessory genome corresponds "
							 "to all the loci not included in the core genome. The accessory genome "
							 "is determined for each core genome threshold.")

	parser.add_argument("--ra", "--rarefaction-analysis", action="store_true",
					 	required=False, dest=ct.COMPUTECORE_RAREFACTION_ANALYSIS_ARGNAME,
						help="Perform rarefaction analysis to evaluate the stability of the core genome "
							 "are randomly selected and the number of core loci (for the core genome) or "
							 "total loci (for the pangenome) is computed for each subset from 1 to the "
							 "total number of samples. The average number of core loci (for the core genome) "
							 "or total loci (for the pangenome) is computed for each sample size subset "
							 "using the values of all permutations to plot the rarefaction curve. The "
							 "rarefaction analysis is performed for each threshold. A power law model "
							 "is fitted to the rarefaction curve to estimate the stability of the core genome "
							 "and the openness of the pangenome.")

	parser.add_argument("--pn", "--permutation-number", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_PERMUTATION_NUMBER_ARGNAME],
					 	required=False, dest=ct.COMPUTECORE_PERMUTATION_NUMBER_ARGNAME,
						help="Number of permutations for the rarefaction analysis. The rarefaction analysis "
							 "is repeated a number of times equal to the value provided to this parameter.")

	parser.add_argument("--ps", "--permutation-samples", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME],
					 	required=False, dest=ct.COMPUTECORE_PERMUTATION_SAMPLES_ARGNAME,
						help="Number of samples randomly selected for each permutation. All samples will be "
							 "used if this value is not provided.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_CPU_CORES_ARGNAME],
						required=False, dest=ct.COMPUTECORE_CPU_CORES_ARGNAME,
						help="Maximum number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads).")

	parser.add_argument("--el", "--exclude-loci", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_EXCLUDE_LOCI_ARGNAME],
						required=False, dest=ct.COMPUTECORE_EXCLUDE_LOCI_ARGNAME,
						help="Path to a file with a list of loci identifiers to exclude from the analysis "
							 "(one locus identifier per line).")

	parser.add_argument("--eg", "--exclude-genomes", type=ct.COMPUTECORE_ARGUMENT_TYPES[ct.COMPUTECORE_EXCLUDE_GENOMES_ARGNAME],
						required=False, dest=ct.COMPUTECORE_EXCLUDE_GENOMES_ARGNAME,
						help="Path to a file with a list of genome identifiers to exclude from the analysis "
							 "(one genome identifier per line).")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.ComputeCoreValidator.model_validate(clean_args)

	compute_core.main(validated_args)


@process_timer
def run_subset_results():
	"""Run the SubsetResults module to subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers."""

	try:
		from SubsetResults import constants as ct
		from SubsetResults import validation_models as vm
		from SubsetResults import subset_results
	except ModuleNotFoundError:
		from CHEWBBACA.SubsetResults import constants as ct
		from CHEWBBACA.SubsetResults import validation_models as vm
		from CHEWBBACA.SubsetResults import subset_results

	def msg(name=None):
		usage_msg = "chewBBACA.py SubsetResults --input-file <file> --loci-list <file> --samples-list <file> --output-file <file> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="SubsetResults",
									 description="Subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.SUBSETRESULTS_MODULE_NAME)}")

	parser.add_argument("SubsetResults", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-files", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_RESULTS_FILES_ARGNAME],
						required=True, dest=ct.SUBSETRESULTS_RESULTS_FILES_ARGNAME,
						help="Path to the directory containing the files to be subsetted.")

	parser.add_argument("-o", "--output-directory", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.SUBSETRESULTS_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--ll", "--loci-list", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_LOCI_LIST_ARGNAME],
						required=False, dest=ct.SUBSETRESULTS_LOCI_LIST_ARGNAME,
						help="Path to a TXT/TSV file containing a list of loci to select, one locus identifier "
							 "per line. If the file contains multiple columns, the loci identifiers must be in "
							 "the first column.")

	parser.add_argument("--sl", "--sample-list", type=ct.SUBSETRESULTS_ARGUMENT_TYPES[ct.SUBSETRESULTS_SAMPLE_LIST_ARGNAME],
						required=False, dest=ct.SUBSETRESULTS_SAMPLE_LIST_ARGNAME,
						help="Path to a TXT/TSV file containing a list of samples to select, one sample identifier "
							 "per line. If the file contains multiple columns, the sample identifiers must be in "
							 "the first column.")

	parser.add_argument("--invert-loci", action="store_true",
						required=False, dest=ct.SUBSETRESULTS_INVERT_LOCI_ARGNAME,
						help="If provided, the process will select the loci that are not in the input loci list.")

	parser.add_argument("--invert-samples", action="store_true",
						required=False, dest=ct.SUBSETRESULTS_INVERT_SAMPLES_ARGNAME,
						help="If provided, the process will select the samples that are not in the input samples list.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.SubsetResultsValidator.model_validate(clean_args)

	subset_results.main(validated_args)


@process_timer
def run_merge_results():
	"""Run the MergeResults module to merge results files created by chewBBACA."""

	try:
		from MergeResults import constants as ct
		from MergeResults import validation_models as vm
		from MergeResults import merge_results
	except ModuleNotFoundError:
		from CHEWBBACA.MergeResults import constants as ct
		from CHEWBBACA.MergeResults import validation_models as vm
		from CHEWBBACA.MergeResults import merge_results

	def msg(name=None):
		usage_msg = "chewBBACA.py MergeResults --results-directories <dir> <dir> ... --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="MergeResults",
									 description="Merge results files created by chewBBACA.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.MERGERESULTS_MODULE_NAME)}")

	parser.add_argument("MergeResults", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-r", "--results-directories", nargs="+", type=ct.MERGERESULTS_ARGUMENT_TYPES[ct.MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME],
						required=True, dest=ct.MERGERESULTS_RESULTS_DIRECTORIES_ARGNAME,
						help="Paths to the directories containing the results files created by chewBBACA. The "
							 "results must have been determined with the same schema and share all the loci "
							 "or a subset of the loci if using the --common parameter.")

	parser.add_argument("-o", "--output-directory", type=ct.MERGERESULTS_ARGUMENT_TYPES[ct.MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.MERGERESULTS_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--common", action="store_true",
						required=False, dest=ct.MERGERESULTS_COMMON_ARGNAME,
						help="Merge the results based on the subset of loci shared between all inputs.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.MergeResultsValidator.model_validate(clean_args)

	merge_results.main(validated_args)


@process_timer
def run_hash_profiles():
	"""Run the HashProfiles module to hash allelic profiles."""

	try:
		from HashProfiles import constants as ct
		from HashProfiles import validation_models as vm
		from HashProfiles import hash_profiles
	except ModuleNotFoundError:
		from CHEWBBACA.HashProfiles import constants as ct
		from CHEWBBACA.HashProfiles import validation_models as vm
		from CHEWBBACA.HashProfiles import hash_profiles

	def msg(name=None):
		usage_msg = "chewBBACA.py HashProfiles --input-file <file> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="HashProfiles",
									 description="Hash allelic profiles.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.HASHPROFILES_MODULE_NAME)}")

	parser.add_argument("HashProfiles", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-a", "--allelic-profiles", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_ALLELIC_PROFILES_ARGNAME],
						required=True, dest=ct.HASHPROFILES_ALLELIC_PROFILES_ARGNAME,
						help="Path to the TSV file that contains the allelic profiles determined by "
							 "the AlleleCall module.")

	parser.add_argument("-s", "--schema-directory", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.HASHPROFILES_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory to get the allele sequences and compute the hashes.")

	parser.add_argument("-o", "--output-directory", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.HASHPROFILES_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--ha", "--hashing-algo", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_HASHING_ALGO_ARGNAME],
					 	required=False, dest=ct.HASHPROFILES_HASHING_ALGO_ARGNAME,
						help="Hashing algorithm used to hash the allele sequences. The hashing algorithms "
							 "implemented in the hashlib and zlib Python libraries are supported.")

	parser.add_argument("--nc", "--nrow-chunk", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_NROW_CHUNK_ARGNAME],
					 	required=False, dest=ct.HASHPROFILES_NROW_CHUNK_ARGNAME,
						help="Divide the input file into chunks of this many rows to process larger files"
						 	 " more efficiently.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.HASHPROFILES_ARGUMENT_TYPES[ct.HASHPROFILES_CPU_CORES_ARGNAME],
						required=False, dest=ct.HASHPROFILES_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.HashProfilesValidator.model_validate(clean_args)

	hash_profiles.main(validated_args)


@process_timer
def run_get_alleles():
	"""Run the GetAlleles module to create FASTA files containing the alleles identified by the AlleleCall module."""

	try:
		from GetAlleles import constants as ct
		from GetAlleles import validation_models as vm
		from GetAlleles import get_alleles
	except ModuleNotFoundError:
		from CHEWBBACA.GetAlleles import constants as ct
		from CHEWBBACA.GetAlleles import validation_models as vm
		from CHEWBBACA.GetAlleles import get_alleles

	def msg(name=None):
		usage_msg = "chewBBACA.py GetAlleles --input-file <file> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="GetAlleles",
									 description="Create FASTA files containing the alleles identified by the AlleleCall module.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.GETALLELES_MODULE_NAME)}")

	parser.add_argument("GetAlleles", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-a", "--allelic-profiles", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_ALLELIC_PROFILES_ARGNAME],
						required=True, dest=ct.GETALLELES_ALLELIC_PROFILES_ARGNAME,
						help="Path to the TSV file containing the allelic profiles.")

	parser.add_argument("-s", "--schema-directory", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.GETALLELES_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema directory.")

	parser.add_argument("-ll", "--loci-list", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_LOCI_LIST_ARGNAME],
						required=False, dest=ct.GETALLELES_LOCI_LIST_ARGNAME,
						help="Path to a file with the list of genes/loci to create FASTA files for. The "
							 "file must include the identifiers of the loci, one per line, without the "
							 ".fasta extension.")

	parser.add_argument("-o", "--output-directory", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.GETALLELES_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_CPU_CORES_ARGNAME],
						required=False, dest=ct.GETALLELES_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

	parser.add_argument("--tt", "--translation-table", type=ct.GETALLELES_ARGUMENT_TYPES[ct.GETALLELES_TRANSLATION_TABLE_ARGNAME],
					 	required=False, dest=ct.GETALLELES_TRANSLATION_TABLE_ARGNAME,
						help="Genetic code used to translate coding DNA sequences (CDSs). If no value"
							 " is specified, the process tries to get the value stored in the schema "
							 "config file. If the schema does not include a config file, the process "
							 "uses the default translation table (11).")

	parser.add_argument("--ta", "--translate-alleles", action="store_true",
						required=False, dest=ct.GETALLELES_TRANSLATE_ALLELES_ARGNAME,
						help="Create FASTA files with the translated alleles.")

	parser.add_argument("--distinct", action="store_true",
						required=False, dest=ct.GETALLELES_DISTINCT_ARGNAME,
						help="Only include the distinct alleles in the output FASTA files.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.GetAllelesValidator.model_validate(clean_args)

	get_alleles.main(validated_args)


@process_timer
def run_adapt_schema():
	"""Run the AdaptSchema module to adapt a typing schema."""

	try:
		from AdaptSchema import constants as ct
		from AdaptSchema import validation_models as vm
		from AdaptSchema import adapt_schema
	except ModuleNotFoundError:
		from CHEWBBACA.AdaptSchema import constants as ct
		from CHEWBBACA.AdaptSchema import validation_models as vm
		from CHEWBBACA.AdaptSchema import adapt_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py AdaptSchema --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="AdaptSchema",
									 description="Adapt an external schema to be used with chewBBACA.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.ADAPTSCHEMA_MODULE_NAME)}")

	parser.add_argument("AdaptSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-s", "--schema-directory", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.ADAPTSCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the directory of the schema to adapt. The schema must contain one "
							 "FASTA file per gene/locus.")

	parser.add_argument("-o", "--output-directory", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.ADAPTSCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the adapted schema will be created.")

	parser.add_argument("--ll", "--loci-list", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_LOCI_LIST_ARGNAME],
						required=False, dest=ct.ADAPTSCHEMA_LOCI_LIST_ARGNAME,
						help="Path to a file with the list of loci in the schema that the process should "
							 "adapt (one per line, full paths or loci IDs).")

	parser.add_argument("-gp", "--gene-predictor", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME],
							required=False, dest=ct.ADAPTSCHEMA_GENE_PREDICTOR_ARGNAME,
							help="Specify which gene prediction software to use. Default is Pyrodigal "
								 "to predict genes from prokaryotic genomes. AUGUSTUS can predict genes "
								 "for prokaryotic and eukaryotic genomes.")
	
	parser.add_argument("--gpo", "--gene-prediction-options", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME],
						nargs="+", required=False, dest=ct.ADAPTSCHEMA_GENE_PREDICTION_OPTIONS_ARGNAME,
						help="List of arguments passed to configure the gene prediction. When providing "
							 "genome assemblies in FASTA format, the list of arguments for each parameter "
							 "used to configure the gene prediction can be passed as the long format of "
							 "the parameter name followed by the argument value (e.g., pyrodigal-training"
							 "-file=/path/to/file).")

	parser.add_argument("--co", "--clustering-options", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME],
						nargs="+", required=False, dest=ct.ADAPTSCHEMA_CLUSTERING_OPTIONS_ARGNAME,
						help="List of arguments passed to configure the clustering of predicted "
								"coding sequences (CDSs). The list of arguments for each parameter "
								"used to configure the clustering can be passed as the long format "
								"of the parameter name followed by the argument value (e.g., word-size=).")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME],
						required=False, dest=ct.ADAPTSCHEMA_BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio (BSR) value. The process selects representative alleles "
							 "for each locus based on this value. Representative alleles are selected "
							 "until all alleles in a locus align against one of the representatives "
							 "with a BSR >= than the specified value.")

	parser.add_argument("--ml", "--minimum-length", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME],
						required=False, dest=ct.ADAPTSCHEMA_MINIMUM_LENGTH_ARGNAME,
						help="Minimum sequence length value stored in the schema config file. The "
							 "schema adaptation process will only discard sequences smaller than "
							 "this value if the --size-filter parameter is provided.")

	parser.add_argument("--tt", "--translation-table", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME],
					 	required=False, dest=ct.ADAPTSCHEMA_TRANSLATION_TABLE_ARGNAME,
						help="Genetic code used for allele translation. This value is ignored if "
							 "a valid training file is passed to `--ptf`, `--training-file`.")

	parser.add_argument("--st", "--size-threshold", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME],
						required=False, dest=ct.ADAPTSCHEMA_SIZE_THRESHOLD_ARGNAME,
						help="Allele size variation threshold value stored in the schema config "
							 "file. The schema adaptation process will only discard alleles with "
							 "a size that deviates from the locus length mode +- the size theshold "
							 "value if the --size-filter parameter is provided.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_CPU_CORES_ARGNAME],
						required=False, dest=ct.ADAPTSCHEMA_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--bp", "--blast-path", type=ct.ADAPTSCHEMA_ARGUMENT_TYPES[ct.ADAPTSCHEMA_BLAST_PATH_ARGNAME],
						required=False, dest=ct.ADAPTSCHEMA_BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--size-filter", action="store_true",
						required=False, dest=ct.ADAPTSCHEMA_SIZE_FILTER_ARGNAME,
						help="Apply the minimum length and size threshold values to filter out "
							 "alleles during schema adaptation.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.AdaptSchemaValidator.model_validate(clean_args)

	adapt_schema.main(validated_args)


@process_timer
def run_annotate_schema():
	"""Run the AnnotateSchema module to annotate loci in a schema."""

	try:
		from AnnotateSchema import constants as ct
		from AnnotateSchema import validation_models as vm
		from AnnotateSchema import annotate_schema
	except ModuleNotFoundError:
		from CHEWBBACA.AnnotateSchema import constants as ct
		from CHEWBBACA.AnnotateSchema import validation_models as vm
		from CHEWBBACA.AnnotateSchema import annotate_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py AnnotateSchema --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="AnnotateSchema",
									 description="Retrieve annotations for loci in a schema.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.ANNOTATESCHEMA_MODULE_NAME)}")

	parser.add_argument("AnnotateSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-s", "--schema-directory", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.ANNOTATESCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema's directory.")

	parser.add_argument("-o", "--output-directory", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.ANNOTATESCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store intermediate "
							 "files and save the final TSV file with the loci annotations.")

	parser.add_argument("--ll", "--loci-list", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_LOCI_LIST_ARGNAME],
						required=False, dest=ct.ANNOTATESCHEMA_LOCI_LIST_ARGNAME,
						help="Path to a file with the list of loci in the schema that the process "
							 "should find annotations for (one per line, full paths or loci IDs).")

	parser.add_argument("--pt", "--protein-table", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME],
						required=False, dest=ct.ANNOTATESCHEMA_PROTEIN_TABLE_ARGNAME,
						help="Path to the TSV file with coding sequence (CDS) coordinate data, "
							 "`cds_coordinates.tsv`, created by the CreateSchema process.")

	parser.add_argument("--bsr", "--blast-score-ratio", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME],
					 	required=False, dest=ct.ANNOTATESCHEMA_BLAST_SCORE_RATIO_ARGNAME,
						help="BLAST Score Ratio value. The BSR is only used when taxa names are "
							 "provided to the --taxa parameter and local sequences are aligned "
							 "against reference proteomes downloaded from UniProt. Annotations "
							 "are selected based on a BSR >= than the specified value.")

	parser.add_argument("--tt", "--translation-table", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME],
					 	required=False, dest=ct.ANNOTATESCHEMA_TRANSLATION_TABLE_ARGNAME,
						help="Genetic code used to translate coding DNA sequences (CDSs). If no value"
							 " is specified, the process tries to get the value stored in the schema "
							 "config file. If the schema does not include a config file, the process "
							 "uses the default translation table (11).")

	parser.add_argument("--cc", "--cpu-cores", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_CPU_CORES_ARGNAME],
						required=False, dest=ct.ANNOTATESCHEMA_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")

	parser.add_argument("--ptax", "--proteome-taxa", nargs="+", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME],
						required=False, dest=ct.ANNOTATESCHEMA_PROTEOME_TAXA_ARGNAME,
						help="List of scientific names for a set of taxa. The process will download "
							 "reference proteomes from UniProt associated to taxa names that contain "
							 "any of the provided terms. The schema representative alleles are aligned "
							 "against the reference proteomes to assign annotations based on high-BSR "
							 "matches.")

	parser.add_argument("--pm", "--proteome-matches", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME],
					 	required=False, dest=ct.ANNOTATESCHEMA_PROTEOME_MATCHES_ARGNAME,
						help="Maximum number of proteome matches to report.")

	parser.add_argument("--bp", "--blast-path", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_BLAST_PATH_ARGNAME],
						required=False, default='', dest=ct.ANNOTATESCHEMA_BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--rapi", "--rest-api", action="store_true",
	                    required=False, dest=ct.ANNOTATESCHEMA_UNIPROT_REST_API_ARGNAME,
						help="Search for annotations using UniProt's Proteins API.")

	parser.add_argument("--rm", "--rest-api-matches", type=ct.ANNOTATESCHEMA_ARGUMENT_TYPES[ct.ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME],
					 	required=False, dest=ct.ANNOTATESCHEMA_UNIPROT_REST_API_MATCHES_ARGNAME,
						help="Maximum number of matches retrieved from UniProt's Proteins API to report.")

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, dest=ct.ANNOTATESCHEMA_NO_CLEANUP_ARGNAME,
						help="If provided, intermediate files generated during process execution "
							 "are not removed at the end.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.AnnotateSchemaValidator.model_validate(clean_args)

	annotate_schema.main(validated_args)


@process_timer
def run_compute_distances():
	"""Run the ComputeDistances module to compute pairwise distances based on allele calling results."""

	try:
		from ComputeDistances import constants as ct
		from ComputeDistances import validation_models as vm
		from ComputeDistances import compute_distances
	except ModuleNotFoundError:
		from CHEWBBACA.ComputeDistances import constants as ct
		from CHEWBBACA.ComputeDistances import validation_models as vm
		from CHEWBBACA.ComputeDistances import compute_distances

	def msg(name=None):
		usage_msg = "chewBBACA.py ComputeDistances --input-file <file> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="ComputeDistances",
									 description="Compute pairwise distances based on allele calling results.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.COMPUTEDISTANCES_MODULE_NAME)}")

	parser.add_argument("ComputeDistances", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-a", "--allelic-profiles", type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME],
						required=True, dest=ct.COMPUTEDISTANCES_ALLELIC_PROFILES_ARGNAME,
						help="Path to a TSV file containing allelic profiles determined by the AlleleCall "
							 "module.")

	parser.add_argument("-o", "--output-directory", type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.COMPUTEDISTANCES_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store intermediate "
							 "and final results.")

	parser.add_argument("--dcm", "--distance-computation-method", type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME],
					 	required=False, dest=ct.COMPUTEDISTANCES_DISTANCE_COMPUTATION_METHOD_ARGNAME,
						help="Distance method used to compute the distance matrix. The module supports "
							 "the hamming, jaccard, loci (number of loci not shared), and core (number of "
							 "different alleles for core loci) methods.")

	parser.add_argument('--of', '--output-format', type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME],
					 	required=False, dest=ct.COMPUTEDISTANCES_OUTPUT_FORMAT_ARGNAME,
						help="Output format for the distance matrix (upper_triangular, lower_triangular, "
							 "symmetric, table).")

	parser.add_argument("--cs", "--compute-similarity", action="store_true",
					 	required=False, dest=ct.COMPUTEDISTANCES_COMPUTE_SIMILARITY_ARGNAME,
						help="Compute similarity values instead of distance values.")

	parser.add_argument("--no-mask", action="store_true",
					 	required=False, dest=ct.COMPUTEDISTANCES_NO_MASK_ARGNAME,
						help="Do not mask missing data when computing the distance matrix. This option "
							 "is useful when the input profiles are already masked.")

	parser.add_argument('--cc', '--cpu-cores', type=ct.COMPUTEDISTANCES_ARGUMENT_TYPES[ct.COMPUTEDISTANCES_CPU_CORES_ARGNAME],
						required=False, default=1, dest=ct.COMPUTEDISTANCES_CPU_CORES_ARGNAME,
						help='Number of CPU cores/threads that will be '
							 'used to run the process (chewie resets to a '
							 'lower value if it is equal to or exceeds the total '
							 'number of available CPU cores/threads).')

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.ComputeDistancesValidator.model_validate(clean_args)

	compute_distances.main(validated_args)


@process_timer
def run_compute_msa():
	"""Run the ComputeMSA module to compute a Multiple Sequence Alignment based on allele calling results."""

	try:
		from ComputeMSA import constants as ct
		from ComputeMSA import validation_models as vm
		from ComputeMSA import compute_msa
	except ModuleNotFoundError:
		from CHEWBBACA.ComputeMSA import constants as ct
		from CHEWBBACA.ComputeMSA import validation_models as vm
		from CHEWBBACA.ComputeMSA import compute_msa

	def msg(name=None):
		usage_msg = "chewBBACA.py ComputeMSA --input-file <file> --schema-directory <dir> --output-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="ComputeMSA",
									 description="Compute a Multiple Sequence Alignment based on allele calling results.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.COMPUTEMSA_MODULE_NAME)}")

	parser.add_argument("ComputeMSA", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-i", "--input-path", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_INPUT_FILES_ARGNAME],
						required=True, dest=ct.COMPUTEMSA_INPUT_FILES_ARGNAME,
						help="Path to a TSV file containing allelic profiles or to a folder containing "
							 "FASTA files. If a TSV file containing allelic profiles is provided, it is "
							 "necessary to provide the path to the schema to the `--schema-directory` "
							 "parameter. The module will create a FASTA file with the alleles identified "
							 "in the samples for each schema locus and compute a MSA. The loci MSAs are "
							 "joined to create the complete MSA based on the allele calling results. If "
							 "a path to a folder is provided, the module computes a MSA for each FASTA "
							 "file in the folder, but will not attempt to join the MSAs as it does not "
							 "have the sample information (in this case, it is not necessary to pass the "
							 "schema path).")

	parser.add_argument("-o", "--output-directory", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.COMPUTEMSA_OUTPUT_DIRECTORY_ARGNAME,
						help="Path to the output directory where the process will store intermediate and "
							 "final results.")

	parser.add_argument("-s", "--schema-directory", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME],
						required=False, dest=ct.COMPUTEMSA_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the schema\'s directory. This parameter is only required if the "
							 "input is a TSV file with allelic profiles.")

	parser.add_argument("--tt", "--translation-table", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_TRANSLATION_TABLE_ARGNAME],
						required=False, dest=ct.COMPUTEMSA_TRANSLATION_TABLE_ARGNAME,
						help="Genetic code used for sequence translation.")

	parser.add_argument("--dm", "--dna-msa", action="store_true",
						required=False, dest=ct.COMPUTEMSA_DNA_MSA_ARGNAME,
						help="Converts the protein MSA back to DNA to create an additional output file "
							 "with the DNA MSA.")

	parser.add_argument("--ov", "--output-variable", action="store_true",
					 	required=False, dest=ct.COMPUTEMSA_OUTPUT_VARIABLE_ARGNAME,
						help="Output a reduced MSA including only the variable positions. If the "
							 "`--dna-msa` parameter is provided, the process will output a reduced MSA "
							 "for both the protein and DNA MSAs.")

	parser.add_argument("--gaps", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_GAPS_ARGNAME],
						required=False, dest=ct.COMPUTEMSA_GAPS_ARGNAME,
						help="How to treat gaps when determining the reduced MSA for the variable "
							 "positions. The default value, `exclude`, removes variable positions if "
							 "any of the aligned sequences contain a gap. The `ignore` option allows "
							 "to consider variable positions that include gaps in some sequences as "
							 "long as other sequences include variable non-gap characters. The character "
							 "used to represent gaps is `-`.")

	parser.add_argument("--ambiguous", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_AMBIGUOUS_ARGNAME],
						required=False, dest=ct.COMPUTEMSA_AMBIGUOUS_ARGNAME,
						help="How to treat ambiguous amino acids or nucleotides when determining the "
							 "reduced MSA for the variable positions. The default value, `exclude`, "
							 "removes variable positions if any of the aligned sequences contain an "
							 "ambiguous amino acid or nucleotide. The `ignore` option allows to consider "
							 "variable positions that include ambiguous amino acids or nucleotides in "
							 "some sequences as long as other sequences include variable non-ambiguous "
							 "characters. The characters interpreted as ambiguous amino acids are "
							 "[B, Z, X, J]. The characters interpreted as ambiguous nucleotides are"
                             " [R, Y, S, W, K, M, B, D, H, V, N].")

	parser.add_argument("--only-loci-msas", action="store_true",
						required=False, dest=ct.COMPUTEMSA_ONLY_LOCI_MSAS_ARGNAME,
						help="Do not compute the full MSA when the input file is a TSV file containing "
							 "allelic profiles (this is already the default when the input is a path to "
							 "a folder with FASTA files).")

	parser.add_argument("--cmo", "--custom-mafft-options", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME],
						required=False, dest=ct.COMPUTEMSA_CUSTOM_MAFFT_OPTIONS_ARGNAME,
						help="Custom parameters to pass to MAFFT when computing the loci MSAs. The "
							 "value must be a single string with all parameters enclosed in quotes "
							 "(e.g. '--retree 1 --maxiterate 0').")

	parser.add_argument("--pi", "--protein-input", action="store_true",
						required=False, dest=ct.COMPUTEMSA_PROTEIN_INPUT_ARGNAME,
						help="Input files contain protein sequences. This option skips sequence translation "
							 "and is only valid for cases when users provide a path to a directory containing "
							 "FASTA files.")

	parser.add_argument("--no-cleanup", action="store_true",
						required=False, dest=ct.COMPUTEMSA_NO_CLEANUP_ARGNAME,
						help="Keep intermediate files with locus/file MSAs and sample MSAs if input "
							 "is a TSV file containing allelic profiles.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.COMPUTEMSA_ARGUMENT_TYPES[ct.COMPUTEMSA_CPU_CORES_ARGNAME],
						required=False, dest=ct.COMPUTEMSA_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads).")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.ComputeMSAValidator.model_validate(clean_args)

	compute_msa.main(validated_args)


@process_timer
def run_download_schema():
	"""Run the DownloadSchema module to download a schema from Chewie-NS."""

	try:
		from CHEWBBACA_NS import constants as ct
		from CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA_NS import download_schema
	except ModuleNotFoundError:
		from CHEWBBACA.CHEWBBACA_NS import constants as ct
		from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA.CHEWBBACA_NS import download_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py DownloadSchema --species-id <id> --schema-id <id> --download-folder <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="DownloadSchema",
									 description="Download a schema from Chewie-NS.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.DOWNLOADSCHEMA_MODULE_NAME)}")

	parser.add_argument("DownloadSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-o", "--output-directory", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME],
						required=True, dest=ct.DOWNLOADSCHEMA_OUTPUT_DIRECTORY_ARGNAME,
						help="Output folder to which the schema will be saved.")

	parser.add_argument("--spi", "--species-id", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_SPECIES_ID_ARGNAME],
						required=True, dest=ct.DOWNLOADSCHEMA_SPECIES_ID_ARGNAME,
						help="The integer identifier or name of the species that the schema is "
							 "associated to in Chewie-NS.")

	parser.add_argument("--sci", "--schema-id", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME],
						required=True, dest=ct.DOWNLOADSCHEMA_SCHEMA_ID_ARGNAME,
						help="The URI, integer identifier or name of the schema to download from "
							 "Chewie-NS.")

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, dest=ct.DOWNLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help="The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https://"
							 "tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")

	parser.add_argument("--sd", "--schema-date", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME],
						required=False, dest=ct.DOWNLOADSCHEMA_SCHEMA_DATE_ARGNAME,
						help="Download schema with state from specified date. Must be in the format "
							 "`Y-m-dTH:M:S`.")

	parser.add_argument("--lv", "--latest-version", action="store_true",
						required=False, dest=ct.DOWNLOADSCHEMA_LATEST_VERSION_ARGNAME,
						help="If the compressed version that is available is not the latest, downloads "
							 "all loci FASTA files and constructs schema locally using the AdaptSchema module.")

	parser.add_argument("--bp", "--blast-path", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_BLAST_PATH_ARGNAME],
						required=False, dest=ct.DOWNLOADSCHEMA_BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.DOWNLOADSCHEMA_ARGUMENT_TYPES[ct.DOWNLOADSCHEMA_CPU_CORES_ARGNAME],
						required=False, dest=ct.DOWNLOADSCHEMA_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the total "
							 "number of available CPU cores/threads). This value is only used if it is "
							 "necessary to construct the schema locally.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.DownloadSchemaValidator.model_validate(clean_args)

	download_schema.main(validated_args)


@process_timer
def run_upload_schema():
	"""Run the UploadSchema module to upload a schema to Chewie-NS."""

	try:
		from CHEWBBACA_NS import constants as ct
		from CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA_NS import upload_schema
	except ModuleNotFoundError:
		from CHEWBBACA.CHEWBBACA_NS import constants as ct
		from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA.CHEWBBACA_NS import upload_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py UploadSchema --schema-directory <dir> --species-id <id> --schema-name <name> --loci-prefix <prefix> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="UploadSchema",
									 description="Upload a schema to Chewie-NS.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.UPLOADSCHEMA_MODULE_NAME)}")

	parser.add_argument("UploadSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-s", "--schema-directory", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.UPLOADSCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the directory of the schema to upload.")

	parser.add_argument("--spi", "--species-id", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_SPECIES_ID_ARGNAME],
						required=True, dest=ct.UPLOADSCHEMA_SPECIES_ID_ARGNAME,
						help="The integer identifier or name of the species that the schema will be "
							 "associated to in Chewie-NS.")

	parser.add_argument("--sn", "--schema-name", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_SCHEMA_NAME_ARGNAME],
						required=True, dest=ct.UPLOADSCHEMA_SCHEMA_NAME_ARGNAME,
						help="A brief and meaningful name that should help understand the type and "
							 "content of the schema.")

	parser.add_argument("--lp", "--loci-prefix", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_LOCI_PREFIX_ARGNAME],
						required=True, dest=ct.UPLOADSCHEMA_LOCI_PREFIX_ARGNAME,
						help="Prefix included in the name of each locus of the schema.")

	parser.add_argument("--df", "--description-file", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME],
						required=False, dest=ct.UPLOADSCHEMA_DESCRIPTION_FILE_ARGNAME,
						help="Path to a text file with a description about the schema. Markdown "
							 "syntax is supported in order to offer greater customizability of "
							 "the rendered description in the Frontend. Will default to the schema's "
							 "name if the user does not provide a valid path for a file.")

	parser.add_argument("--af", "--annotations-file", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME],
						required=False, dest=ct.UPLOADSCHEMA_ANNOTATIONS_FILE_ARGNAME,
						help="Path to a TSV file with loci annotations. The first column has loci "
							 "identifiers (w/o .fasta extension), the second has UniProt protein "
							 "names, the third has UniProt gene names, the fourth has UniProt URIs, "
							 "the fifth has user annotations, and the sixth has custom annotations.")

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, dest=ct.UPLOADSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help="The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https://"
							 "tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")

	parser.add_argument("--continue-upload", action="store_true",
						required=False, dest=ct.UPLOADSCHEMA_CONTINUE_UPLOAD_ARGNAME,
						help="Check if the schema upload was interrupted and attempt to continue upload.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.UPLOADSCHEMA_ARGUMENT_TYPES[ct.UPLOADSCHEMA_CPU_CORES_ARGNAME],
						required=False, dest=ct.UPLOADSCHEMA_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads). This value is used to accelerate the quality "
							 "control step that checks all alleles in the schema.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.UploadSchemaValidator.model_validate(clean_args)

	upload_schema.main(validated_args)


@process_timer
def run_synchronize_schema():
	"""Run the SyncSchema module to synchronize a local schema with the remote version in Chewie-NS."""

	try:
		from CHEWBBACA_NS import constants as ct
		from CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA_NS import synchronize_schema
	except ModuleNotFoundError:
		from CHEWBBACA.CHEWBBACA_NS import constants as ct
		from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA.CHEWBBACA_NS import synchronize_schema

	def msg(name=None):
		usage_msg = "chewBBACA.py SyncSchema --schema-directory <dir> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="SyncSchema",
									 description="Synchronize a schema with its remote version in Chewie-NS.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.SYNCSCHEMA_MODULE_NAME)}")

	parser.add_argument("SyncSchema", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("-s", "--schema-directory", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME],
						required=True, dest=ct.SYNCSCHEMA_SCHEMA_DIRECTORY_ARGNAME,
						help="Path to the directory with the schema to be synced.")

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, dest=ct.SYNCSCHEMA_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help="The base URL for the Chewie-NS instance. The default option will get the base "
							 "URL from the schema's URI. It is also possible to specify other options that "
							 "are available in chewBBACA\'s configs, such as: `main` will establish a connection"
							 " to `https://chewbbaca.online/`, `tutorial` to `https://tutorial.chewbbaca.online/` "
							 "and `local` to `http://127.0.0.1:5000/NS/api/` (localhost). Users may also provide "
							 "the IP address to other Chewie-NS instances.")

	parser.add_argument("--sa", "--submit-alleles", action="store_true",
						required=False, dest=ct.SYNCSCHEMA_SUBMIT_ALLELES_ARGNAME,
						help="If the process should identify new alleles in the local schema and send them "
							 "to the Chewie-NS instance. (only authorized users can submit new alleles).")

	parser.add_argument("--bp", "--blast-path", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_BLAST_PATH_ARGNAME],
						required=False, dest=ct.SYNCSCHEMA_BLAST_PATH_ARGNAME,
						help="Path to the directory that contains the BLAST executables.")

	parser.add_argument("--cc", "--cpu-cores", type=ct.SYNCSCHEMA_ARGUMENT_TYPES[ct.SYNCSCHEMA_CPU_CORES_ARGNAME],
						required=False, dest=ct.SYNCSCHEMA_CPU_CORES_ARGNAME,
						help="Number of CPU cores/threads that will be used to run the process (chewie "
							 "resets to a lower value if it is equal to or exceeds the total number of "
							 "available CPU cores/threads). This value is only used if the process retrieves "
							 "novel alleles from the remote schema and needs to redetermine the set of "
							 "representative alleles for the local schema.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.SynchronizeSchemaValidator.model_validate(clean_args)

	synchronize_schema.main(validated_args)


@process_timer
def run_stats_requests():
	"""Run the NSStats module to get information about schemas in Chewie-NS."""

	try:
		from CHEWBBACA_NS import constants as ct
		from CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA_NS import stats_requests
	except ModuleNotFoundError:
		from CHEWBBACA.CHEWBBACA_NS import constants as ct
		from CHEWBBACA.CHEWBBACA_NS import validation_models as vm
		from CHEWBBACA.CHEWBBACA_NS import stats_requests

	def msg(name=None):
		usage_msg = "chewBBACA.py NSStats --mode <mode> [options]"

		return usage_msg

	parser = argparse.ArgumentParser(prog="NSStats",
									 description="Retrieve basic information about the species and schemas in Chewie-NS.",
									 usage=msg(),
									 formatter_class=ModifiedHelpFormatter,
									 epilog=f"Module documentation available at {gct.DOCS_URL_TEMPLATE.format(ct.NSSTATS_MODULE_NAME)}")

	parser.add_argument("NSStats", nargs="+", help=argparse.SUPPRESS)

	parser.add_argument("--sm", "--stats-mode", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_STATS_MODE_ARGNAME],
						required=True, dest=ct.NSSTATS_STATS_MODE_ARGNAME,
						help="The process can retrieve the list of species (`species` option) in Chewie-NS "
							 "or the list of schemas for a species (`schemas` option).")

	parser.add_argument("--spi", "--species-id", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_SPECIES_ID_ARGNAME],
						required=False, dest=ct.NSSTATS_SPECIES_ID_ARGNAME,
						help="The integer identifier of a species in Chewie-NS.")

	parser.add_argument("--sci", "--schema-id", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_SCHEMA_ID_ARGNAME],
						required=False, dest=ct.NSSTATS_SCHEMA_ID_ARGNAME,
						help="The integer identifier of a schema in Chewie-NS.")

	parser.add_argument("--nsi", "--nomenclature-server-instance", type=ct.NSSTATS_ARGUMENT_TYPES[ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME],
						required=False, dest=ct.NSSTATS_NOMENCLATURE_SERVER_INSTANCE_ARGNAME,
						help="The base URL for the Chewie-NS instance. The default value, `main`, will "
							 "establish a connection to `https://chewbbaca.online/`, `tutorial` to `https:"
							 "//tutorial.chewbbaca.online/` and `local` to `http://127.0.0.1:5000/NS/api/` "
							 "(localhost). Users may also provide the IP address to other Chewie-NS instances.")

	args = parser.parse_args()
	# Strip out None values so Pydantic defaults take effect
	clean_args = {k: v for k, v in vars(args).items() if v is not None}
	# Use Pydantic model to validate argument values
	validated_args = vm.NSStatsValidator.model_validate(clean_args)

	stats_requests.main(validated_args)


def main():

	functions_info = {'PredictGenes': ['Predict genes from a set of input genome assemblies.',
									  run_predict_genes],
					  'CreateSchema': ["Create a gene-by-gene schema based on a set of FASTA files with genome assemblies or coding sequences.",
									   run_create_schema],
					  'AlleleCall': ['Determine the allelic profiles of a set of bacterial genomes based on a schema.',
									 run_allele_call],
					  'SchemaEvaluator': ['Build an interactive report for schema evaluation.',
										  run_evaluate_schema],
					  'AlleleCallEvaluator': ['Build an interactive report for allele calling results evaluation.',
											  run_evaluate_calls],
					  'ComputeCore': ['Determines the set of loci that constitute the core genome based on loci presence thresholds.',
										run_compute_core],
					  'SubsetResults': ['Subset the data in files created by chewBBACA based on a list of loci and/or sample identifiers.',
										run_subset_results],
					  'AdaptSchema': ['Adapt an external schema to be used with chewBBACA.',
											 run_adapt_schema],
					  'MergeResults': ['Merge results files created by chewBBACA.',
									   run_merge_results],
					  'HashProfiles': ['Hash allelic profiles.',
					  					run_hash_profiles],
					  'GetAlleles': ['Create FASTA files containing the alleles identified by the AlleleCall module.',
					  				 run_get_alleles],
					  'AnnotateSchema': ['Retrieve annotations for loci in a schema.',
										 run_annotate_schema],
					  'ComputeDistances': ['Compute pairwise distances based on allele calling results.',
										  run_compute_distances],
					  'ComputeMSA': ['Compute a Multiple Sequence Alignment based on allele calling results.',
									 run_compute_msa],
					  'DownloadSchema': ['Download a schema from Chewie-NS.',
										 run_download_schema],
					  'UploadSchema': ['Upload a schema to Chewie-NS.',
									 run_upload_schema],
					  'SyncSchema': ['Synchronize a schema with its remote version in Chewie-NS.',
									 run_synchronize_schema],
					  'NSStats': ['Retrieve basic information about the species and schemas in Chewie-NS.',
								  run_stats_requests]}

	print(f'chewBBACA version: {gct.__version__}')
	version_triggers = ['-v', '--v', '-version', '--version']
	if len(sys.argv) > 1 and sys.argv[1] in version_triggers:
		# Exit after printing version
		sys.exit(0)

	print(f'Authors: {gct.AUTHORS}')
	print(f'Github: {gct.REPOSITORY}')
	print(f'Documentation: {gct.DOCUMENTATION}')
	print(f'Contacts: {gct.CONTACTS}\n')

	# Display help message if selected process is not valid
	help_triggers = ['-h', '--h', '-help', '--help']
	if len(sys.argv) == 1 or sys.argv[1] not in functions_info or sys.argv[1] in help_triggers:
		exit_code = 0
		# Detect if user passed module name that does not exist
		if len(sys.argv) > 1 and sys.argv[1] not in help_triggers:
			print(f"No module named {sys.argv[1]}.\n")
			exit_code = 1
		print("USAGE: chewBBACA.py [module] -h, --help\n")
		print("Select one of the following modules:")
		for f in functions_info:
			print('{0}: {1}'.format(f, functions_info[f][0]))
		sys.exit(exit_code)

	# Check python version
	python_uptodate, python_version = gval.validate_python_version(gct.MIN_PYTHON)
	if not python_uptodate:
		sys.exit(gct.PYTHON_VERSION.format(python_version, gct.MIN_PYTHON))

	# Trigger module help message if no arguments are provided
	if len(sys.argv) == 2 and sys.argv[1] in functions_info:
		sys.argv.append("-h")

	process = sys.argv[1]
	functions_info[process][1]()


if __name__ == "__main__":

	main()
