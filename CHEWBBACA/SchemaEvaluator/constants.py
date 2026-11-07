#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the SchemaEvaluator module.
"""


# Parameter names
SCHEMAEVALUATOR_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
SCHEMAEVALUATOR_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
SCHEMAEVALUATOR_LOCI_LIST_ARGNAME = "loci_list"
SCHEMAEVALUATOR_ANNOTATIONS_FILE_ARGNAME = "annotations_file"
SCHEMAEVALUATOR_TRANSLATION_TABLE_ARGNAME = "translation_table"
SCHEMAEVALUATOR_SIZE_THRESHOLD_ARGNAME = "size_threshold"
SCHEMAEVALUATOR_MINIMUM_LENGTH_ARGNAME = "minimum_length"
SCHEMAEVALUATOR_CPU_CORES_ARGNAME = "cpu_cores"
SCHEMAEVALUATOR_LOCI_REPORTS_ARGNAME = "loci_reports"
SCHEMAEVALUATOR_LIGTH_ARGNAME = "light"
SCHEMAEVALUATOR_ADD_SEQUENCES_ARGNAME = "add_sequences"

# Argument types
SCHEMAEVALUATOR_ARGUMENT_TYPES = {
	"schema_directory": str,
	"output_directory": str,
	"loci_list": str,
	"annotations_file": str,
	"translation_table": int,
	"size_threshold": float,
	"minimum_length": int,
	"cpu_cores": int,
	"loci_reports": bool,
	"light": bool,
	"add_sequences": bool
}

# Parameter descriptions
SCHEMAEVALUATOR_SCHEMA_DIRECTORY_DESCRIPTION = ("Path to the schema's directory.")
SCHEMAEVALUATOR_OUTPUT_DIRECTORY_DESCRIPTION = ("Path to the output directory where the report HTML files will be created.")
SCHEMAEVALUATOR_LOCI_LIST_DESCRIPTION = ("Path to a file with the list of loci in the schema that the process "
							 "should analyse (one per line, full paths or loci IDs).")
SCHEMAEVALUATOR_ANNOTATIONS_FILE_DESCRIPTION = ("Path to the TSV file with loci annotations to include in a datatable "
							 "component of the report.")
SCHEMAEVALUATOR_TRANSLATION_TABLE_DESCRIPTION = ("Translation table used to translate CDSs.")
SCHEMAEVALUATOR_SIZE_THRESHOLD_DESCRIPTION = ("Coding sequence (CDS) size variation threshold. The module identifies "
							 "the alleles with size that deviates from the locus length mode +/- this "
							 "size threshold value.")
SCHEMAEVALUATOR_MINIMUM_LENGTH_DESCRIPTION = ("Minimum sequence length value. The module identifies alleles shorter "
							 "than this value.")
SCHEMAEVALUATOR_CPU_CORES_DESCRIPTION = ("Number of CPU cores/threads that will be used to run the process "
							 "(chewie resets to a lower value if it is equal to or exceeds the "
							 "total number of available CPU cores/threads).")
SCHEMAEVALUATOR_LOCI_REPORTS_DESCRIPTION = ("Create a detailed report page for each locus. The locus report "
							 "includes components with relevant data and analysis results, such "
							 "as allele diversity charts, a MSA for the alignment of the distinct "
							 "translated alleles and a tree drawn with Phylocanvas based on the "
							 "MAFFT guide tree.")
SCHEMAEVALUATOR_LIGTH_DESCRIPTION = ("Skips MSA computation with MAFFT and does not add the Phylogenetic "
							 "Tree and MSA components to the loci reports.")
SCHEMAEVALUATOR_ADD_SEQUENCES_DESCRIPTION = ("Adds Code Editor components with the DNA and Protein sequences to "
							 "the loci reports. The Code Editor is in readonly mode (allows to "
							 "search for and copy text).")

# Basename for the main report file
SCHEMA_REPORT_BASENAME = "schema_report.html"

# HTML template to create Schema Report
# Need to include "./" at start of file paths to work properly when referencing local files
SCHEMA_REPORT_HTML = ("""
<!DOCTYPE html>
<html lang="en">
<head>
	<meta charset="UTF-8" />
	<meta name="viewport" content="width=device-width, initial-scale=1.0" />
	<title>Schema Evaluator - React Edition</title>
</head>
<body style="background-color: #f6f6f6">
	<noscript> You need to enable JavaScript to run this app. </noscript>
	<div id="root"></div>
	<script> preComputedData = {0} </script>
	<script src="./report_bundle.js"></script>
</body>
</html>
""")

# HTML template to create Loci Reports
# need to include '.' at start of file paths to work properly when referencing local files
LOCUS_REPORT_HTML = ("""
<!DOCTYPE html>
<html lang="en">
	<head>
		<meta charset="UTF-8" />
		<meta name="viewport" content="width=device-width, initial-scale=1.0" />
		<title>Schema Evaluator - Individual Analysis</title>
	</head>
	<body style="background-color: #f6f6f6">
		<noscript> You need to enable JavaScript to run this app. </noscript>
		<div id="root"></div>
		<script src="https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.js"></script>
		<link type=text/css rel=stylesheet href=https://s3-eu-west-1.amazonaws.com/biojs/msa/latest/msa.css />
		<script> preComputedDataInd = {0} </script>
		<script src="./report_bundle.js"></script>
	</body>
</html>
""")

# Relative paths to the JS bundles used by SchemaEvaluator
# Main page
SCHEMA_EVALUATOR_SCHEMA_BUNDLE = "report_template_components/src/bundles/SchemaEvaluator/schema_report/report_bundle.js"
# Loci pages
SCHEMA_EVALUATOR_LOCI_BUNDLE = "report_template_components/src/bundles/SchemaEvaluator/loci_reports/report_bundle.js"

# Do not use list of strings as constants if the strings include formatting
# placeholders. Multiple references to the list of strings will have the same
# id and altering the strings with format will not change the list id. In
# multiprocessing it can reference the same list/id in different processess
# and use the latest changes to a string in the list/id when those changes
# might not refer to the current process (returning an incorrect value).

# Column headers for the Invalid Alleles table
INVALID_ALLELES_TABLE_HEADER = ["Allele ID", "Exception Category", "Exception Description"]

# List of possible translation exceptions
TRANSLATION_EXCEPTIONS = ["Extra in frame stop codon",
						  "is not a start codon",
						  "is not a stop codon",
						  "sequence length is not a multiple of 3",
						  "ambiguous or invalid characters"]

# Column headers/names for the table with data about distinct alleles
DISTINCT_ALLELES_TABLE_HEADER = ["Protein Allele ID", "Count", "List of Distinct Alleles"]

# Exceptions
# Table headers for the Schema report summary data table
SCHEMA_SUMMARY_TABLE_HEADER = ("Loci\tAlleles\tValid Alleles\tInvalid alleles\tIncomplete ORF\tAmbiguous Bases\t"
								"Missing Start/Stop Codon\tIn-frame Stop Codon\tAlleles < {0}bp\tAlleles below "
								"threshold\tAlleles above threshold")

# Table headers for the Loci Analysis Table in the Schema report
LOCI_ANALYSIS_TABLE_HEADER = ("Locus\tTotal Alleles\tValid Alleles\tInvalid Alleles\tProportion of Validated Alleles\tDistinct "
						 	  "Protein Alleles\tIncomplete ORF\tAmbiguous Bases\tMissing Start/Stop Codon\tIn-frame Stop Codon\tAlleles "
						 	  "< {0}bp\tAlleles below threshold\tAlleles above threshold\tMissing Allele IDs")

# Column headers for the Summary Table in the Loci reports
LOCI_REPORTS_SUMMARY_TABLE_HEADER = ("Locus\tTotal Alleles\tValid Alleles\tInvalid Alleles\tProportion of Validated Alleles\tDistinct "
				 					 "Protein Alleles\tIncomplete ORF\tAmbiguous Bases\tMissing Start/Stop Codon\tIn-frame Stop Codon\t"
									 "Alleles < {0}bp\tSize Range (bp)\tLength Median (bp)\tLength Mode (bp)\tAlleles below threshold "
									 "({1}bp)\tAlleles above threshold ({2}bp)\tMissing Allele IDs")
