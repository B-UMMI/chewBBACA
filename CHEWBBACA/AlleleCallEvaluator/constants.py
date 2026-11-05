#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This file stores constants used by the AlleleCallEvaluator module.
"""


# Parameter names
ALLELECALLEVALUATOR_RESULTS_FILES_ARGNAME = "results_files"
ALLELECALLEVALUATOR_SCHEMA_DIRECTORY_ARGNAME = "schema_directory"
ALLELECALLEVALUATOR_OUTPUT_DIRECTORY_ARGNAME = "output_directory"
ALLELECALLEVALUATOR_ANNOTATIONS_FILE_ARGNAME = "annotations_file"
ALLELECALLEVALUATOR_CPU_CORES_ARGNAME = "cpu_cores"
ALLELECALLEVALUATOR_NO_PRESENCE_ABSENCE_ARGNAME = "no_presence_absence"
ALLELECALLEVALUATOR_NO_DISTANCE_MATRIX_ARGNAME = "no_distance_matrix"
ALLELECALLEVALUATOR_NO_NEIGHBOR_JOINING_ARGNAME = "no_neighbor_joining"
ALLELECALLEVALUATOR_FORCE_CORE_MSA_ARGNAME = "force_core_msa"
ALLELECALLEVALUATOR_LIGTH_ARGNAME = "light"

# Argument types
ALLELECALLEVALUATOR_ARGUMENT_TYPES = {
	"results_files": str,
	"schema_directory": str,
	"output_directory": str,
	"annotations_file": str,
	"cpu_cores": int,
	"no_presence_absence": bool,
	"no_distance_matrix": bool,
	"no_neighbor_joining": bool,
	"force_core_msa": bool,
	"light": bool
}

# HTML template to create the main page of the HTML report
ALLELECALL_REPORT_HTML = ("""
<!DOCTYPE html>
<html lang="en">
<head>
	<meta charset="UTF-8" />
	<meta name="viewport" content="width=device-width, initial-scale=1.0" />
	<title>AlleleCall Report - React Edition</title>
</head>
<body style="background-color: #f6f6f6">
	<noscript> You need to enable JavaScript to run this app. </noscript>
	<div id="root"></div>
	<script> preComputedData = {0} </script>
	<script src="./report_bundle.js"></script>
</body>
</html>
""")

# Relative path to the main JS bundle
ALLELECALL_EVALUATOR_BUNDLE = "report_template_components/src/bundles/AlleleCallEvaluator/report_bundle.js"

# Basename for the output files created by the AlleleCallEvaluator module
DISTANCE_MATRIX_BASENAME = "distance_matrix.tsv"
CORE_MSA_BASENAME = "cgMLST_MSA.fasta"
CGMLST_PROFILES_BASENAME = "cgMLST_profiles.tsv"
ALLELECALL_REPORT_BASENAME = "allelecall_report.html"

# Table headers for the Sample Stats table in the HTML report
SAMPLE_STATS_HEADER = ["Sample", "Total Contigs", "Total CDSs", "Proportion of Classified CDSs",
                       "Identified Loci", "Proportion of Identified Loci", "Valid Classifications",
					   "Invalid Classifications"]

# Column headers for the Loci Stats table in the HTML report
LOCUS_STATS_HEADER = ["Locus", "Total CDSs", "Valid Classes", "Invalid Classes", "Proportion Samples"]

# Column headers for the Summary Stats table in the HTML report
SUMMARY_STATS_HEADER = ["Total Samples", "Total Loci", "Total CDSs", "Total CDSs Classified", "EXC",
                        "INF", "PLOT3", "PLOT5", "LOTSC", "NIPH", "NIPHEM", "ALM", "ASM", "PAMA", "LNF"]
