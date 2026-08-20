#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module contains variables used to test chewBBACA modules with pytest.

Code documentation
------------------
"""

# Test variables used by pytest

# Test version and help messages for all modules
CHEWIE_TEST_H = ['chewBBACA.py', '-h']
CHEWIE_TEST_HELP = ['chewBBACA.py', '--help']
CHEWIE_TEST_V = ['chewBBACA.py', '-v']
CHEWIE_TEST_VERSION = ['chewBBACA.py', '--version']
CREATESCHEMA_TEST_H = ['chewBBACA.py', 'CreateSchema', '-h']
CREATESCHEMA_TEST_HELP = ['chewBBACA.py', 'CreateSchema', '--help']
ALLELECALL_TEST_H = ['chewBBACA.py', 'AlleleCall', '-h']
ALLELECALL_TEST_HELP = ['chewBBACA.py', 'AlleleCall', '--help']
SCHEMAEVALUATOR_TEST_H = ['chewBBACA.py', 'SchemaEvaluator', '-h']
SCHEMAEVALUATOR_TEST_HELP = ['chewBBACA.py', 'SchemaEvaluator', '--help']
ALLELECALL_EVALUATOR_TEST_H = ['chewBBACA.py', 'AlleleCallEvaluator', '-h']
ALLELECALL_EVALUATOR_TEST_HELP = ['chewBBACA.py', 'AlleleCallEvaluator', '--help']
EXTRACTCGMLST_TEST_H = ['chewBBACA.py', 'ExtractCgMLST', '-h']
EXTRACTCGMLST_TEST_HELP = ['chewBBACA.py', 'ExtractCgMLST', '--help']
REMOVEGENES_TEST_H = ['chewBBACA.py', 'RemoveGenes', '-h']
REMOVEGENES_TEST_HELP = ['chewBBACA.py', 'RemoveGenes', '--help']
PREPEXTERNALSCHEMA_TEST_H = ['chewBBACA.py', 'PrepExternalSchema', '-h']
PREPEXTERNALSCHEMA_TEST_HELP = ['chewBBACA.py', 'PrepExternalSchema', '--help']
JOINPROFILES_TEST_H = ['chewBBACA.py', 'JoinProfiles', '-h']
JOINPROFILES_TEST_HELP = ['chewBBACA.py', 'JoinProfiles', '--help']
UNIPROTFINDER_TEST_H = ['chewBBACA.py', 'UniprotFinder', '-h']
UNIPROTFINDER_TEST_HELP = ['chewBBACA.py', 'UniprotFinder', '--help']
DOWNLOADSCHEMA_TEST_H = ['chewBBACA.py', 'DownloadSchema', '-h']
DOWNLOADSCHEMA_TEST_HELP = ['chewBBACA.py', 'DownloadSchema', '--help']
LOADSCHEMA_TEST_H = ['chewBBACA.py', 'LoadSchema', '-h']
LOADSCHEMA_TEST_HELP = ['chewBBACA.py', 'LoadSchema', '--help']
SYNCSCHEMA_TEST_H = ['chewBBACA.py', 'SyncSchema', '-h']
SYNCSCHEMA_TEST_HELP = ['chewBBACA.py', 'SyncSchema', '--help']
NSSTATS_TEST_H = ['chewBBACA.py', 'NSStats', '-h']
NSSTATS_TEST_HELP = ['chewBBACA.py', 'NSStats', '--help']
FAKEMODULE_TEST_H = ['chewBBACA.py', 'FakeModule', '-h']
FAKEMODULE_TEST_HELP = ['chewBBACA.py', 'FakeModule', '--help']

# AlleleCall
# AlleleCall template command
ALLELECALL_TEST_GENOME_TEMPLATE = ['chewBBACA.py', 'AlleleCall',
				 				   '-i', 'data/allelecall_data/test_genome',
				 				   '-g', 'data/allelecall_data/sagalactiae_schema',
				 				   '-o', 'allelecall_results']

# AlleleCall execution modes
ALLELECALL_TEST_DEFAULT = ALLELECALL_TEST_GENOME_TEMPLATE[:]
ALLELECALL_TEST_MODE1 = ALLELECALL_TEST_GENOME_TEMPLATE[:]+['--mode', '1']
ALLELECALL_TEST_MODE2 = ALLELECALL_TEST_GENOME_TEMPLATE[:]+['--mode', '2']
ALLELECALL_TEST_MODE3 = ALLELECALL_TEST_GENOME_TEMPLATE[:]+['--mode', '3']
ALLELECALL_TEST_MODE4 = ALLELECALL_TEST_GENOME_TEMPLATE[:]+['--mode', '4']

# AlleleCall --cds-input option
ALLELECALL_TEST_CDS_TEMPLATE = ['chewBBACA.py', 'AlleleCall',
				 				'-i', 'data/allelecall_data/test_cds_input',
				 				'-g', 'data/allelecall_data/sagalactiae_schema',
				 				'-o', 'allelecall_results',
								'--cds-input']

# AlleleCall --cds-input option
ALLELECALL_TEST_CDS_DEFAULT = ALLELECALL_TEST_CDS_TEMPLATE[:]
ALLELECALL_TEST_CDS_MODE1 = ALLELECALL_TEST_CDS_TEMPLATE[:]+['--mode', '1']
ALLELECALL_TEST_CDS_MODE2 = ALLELECALL_TEST_CDS_TEMPLATE[:]+['--mode', '2']
ALLELECALL_TEST_CDS_MODE3 = ALLELECALL_TEST_CDS_TEMPLATE[:]+['--mode', '3']
ALLELECALL_TEST_CDS_MODE4 = ALLELECALL_TEST_CDS_TEMPLATE[:]+['--mode', '4']

# AlleleCall with genome list
ALLELECALL_TEST_GENOME_LIST = ['chewBBACA.py', 'AlleleCall',
				 			   '-i', 'data/allelecall_data/test_genomes_list/test_genomes.txt',
				 			   '-g', 'data/allelecall_data/sagalactiae_schema',
				 			   '-o', 'allelecall_results']

# AlleleCall with loci list
ALLELECALL_TEST_LOCI_IDS_EXTENSION = ALLELECALL_TEST_GENOME_TEMPLATE[:]+ \
	['--gl', 'data/allelecall_data/test_genes_list/test_genes_extension.txt']
ALLELECALL_TEST_LOCI_IDS_NOEXTENSION = ALLELECALL_TEST_GENOME_TEMPLATE[:]+ \
	['--gl', 'data/allelecall_data/test_genes_list/test_genes_no_extension.txt']
ALLELECALL_TEST_LOCI_PATHS = ALLELECALL_TEST_GENOME_TEMPLATE[:]+ \
	['--gl', 'data/allelecall_data/test_genes_list/test_genes_path.txt']

# AlleleCall with empty directory
ALLELECALL_TEST_EMPTY_DIR = ['chewBBACA.py', 'AlleleCall',
		   					 '-i', 'empty_dir',
		   					 '-g', 'data/allelecall_data/sagalactiae_schema',
		   					 '-o', 'allelecall_results']

# AlleleCall with empty files
ALLELECALL_TEST_EMPTY_FILES = ['chewBBACA.py', 'AlleleCall',
		   					   '-i', 'data/createschema_data/genome_dir_with_empty_genomes',
		   					   '-g', 'data/allelecall_data/sagalactiae_schema',
		   					   '-o', 'allelecall_results']

# AlleleCall with files that contain no data
ALLELECALL_TEST_ZERO_BYTES = ['chewBBACA.py', 'AlleleCall',
		   					  '-i', 'data/createschema_data/zero_bytes_pair',
		   					  '-g', 'data/allelecall_data/sagalactiae_schema',
		   					  '-o', 'allelecall_results']

# AlleleCall with fake input path
ALLELECALL_TEST_FAKE_PATH = ['chewBBACA.py', 'AlleleCall',
		   					 '-i', 'this/path/aint/real',
		   					 '-g', 'data/allelecall_data/sagalactiae_schema',
		   					 '-o', 'allelecall_results']

# AlleleCall input name includes blank space
ALLELECALL_TEST_BLANK_SPACE = ['chewBBACA.py', 'AlleleCall',
		   					   '-i', 'data/allelecall_data/test_invalid_input_names/blank_spaces',
		   					   '-g', 'data/allelecall_data/sagalactiae_schema',
		   					   '-o', 'allelecall_results']

# AlleleCall input file has unique prefix longer than 30 chars
ALLELECALL_TEST_LONG_PREFIX = ['chewBBACA.py', 'AlleleCall',
		   					   '-i', 'data/allelecall_data/test_invalid_input_names/long_prefix',
		   					   '-g', 'data/allelecall_data/sagalactiae_schema',
		   					   '-o', 'allelecall_results']

# AlleleCall some input files have the same prefix
ALLELECALL_TEST_SAME_PREFIX = ['chewBBACA.py', 'AlleleCall',
		   					   '-i', 'data/allelecall_data/test_invalid_input_names/same_prefix',
		   					   '-g', 'data/allelecall_data/sagalactiae_schema',
		   					   '-o', 'allelecall_results']

# AlleleCall input prefix interpreted as PDB chain ID
ALLELECALL_TEST_PDB_CHAIN = ['chewBBACA.py', 'AlleleCall',
							 '-i', 'data/allelecall_data/test_invalid_input_names/pdb_prefix',
							 '-g', 'data/allelecall_data/sagalactiae_schema',
							 '-o', 'allelecall_results']

# CreateSchema template command
CREATESCHEMA_TEST_GENOME_TEMPLATE = ['chewBBACA.py', 'CreateSchema',
       								 '-i', 'data/createschema_data/mock_genome_dir',
       								 '-o', 'createschema_results',
       								 '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema with genome list
CREATESCHEMA_TEST_GENOME_LIST = ['chewBBACA.py', 'CreateSchema',
       							 '-i', 'data/createschema_data/mock_genome_list/mock_genomes.txt',
       							 '-o', 'createschema_results',
       							 '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema --cds-input option
CREATESCHEMA_TEST_CDS = ['chewBBACA.py', 'CreateSchema',
       					 '-i', 'data/createschema_data/mock_cds_dir',
       					 '-o', 'createschema_results',
       					 '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn',
	   					 '--cds-input']

# CreateSchema with empty files
CREATESCHEMA_TEST_EMPTY_FILES = ['chewBBACA.py', 'CreateSchema',
								 '-i', 'data/createschema_data/genome_dir_with_empty_genomes',
								 '-o', 'createschema_results',
								 '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema with files that contain no data
CREATESCHEMA_TEST_ZERO_BYTES = ['chewBBACA.py', 'CreateSchema',
								'-i', 'data/createschema_data/zero_bytes_pair',
								'-o', 'createschema_results',
								'--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema with invalid path for PTF
CREATESCHEMA_INVALID_PTF_PATH = ['chewBBACA.py', 'CreateSchema',
								 '-i', 'data/createschema_data/mock_schema_dir',
								 '-o', 'createschema_results',
								 '--ptf', 'path/does/not/exist']

# CreateSchema with FASTA files that only contain sequence header
CREATESCHEMA_TEST_HEADER_ONLY = ['chewBBACA.py', 'CreateSchema',
								 '-i', 'data/createschema_data/genome_dir_with_header_fasta_only',
								 '-o', 'createschema_results',
								 '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema with invalid genome
CREATESCHEMA_TEST_INVALID_GENOME = ['chewBBACA.py', 'CreateSchema',
							   		'-i', 'data/createschema_data/invalid_genome_dir',
							   		'-o', 'createschema_results',
							   		'--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema input name includes blank space
CREATESCHEMA_TEST_BLANK_SPACE = ['chewBBACA.py', 'CreateSchema',
		   					     '-i', 'data/allelecall_data/test_invalid_input_names/blank_spaces',
		   					     '-o', 'createschema_results',
							     '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema input file has unique prefix longer than 30 chars
CREATESCHEMA_TEST_LONG_PREFIX = ['chewBBACA.py', 'CreateSchema',
		   					     '-i', 'data/allelecall_data/test_invalid_input_names/long_prefix',
		   					     '-o', 'createschema_results',
							     '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema some input files have the same prefix
CREATESCHEMA_TEST_SAME_PREFIX = ['chewBBACA.py', 'CreateSchema',
		   					     '-i', 'data/allelecall_data/test_invalid_input_names/same_prefix',
		   					     '-o', 'createschema_results',
							     '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# CreateSchema input prefix interpreted as PDB chain ID
CREATESCHEMA_TEST_PDB_CHAIN = ['chewBBACA.py', 'CreateSchema',
							   '-i', 'data/allelecall_data/test_invalid_input_names/pdb_prefix',
							   '-o', 'createschema_results',
							   '--ptf', 'data/createschema_data/Streptococcus_agalactiae.trn']

# AlleleCallEvaluator
# AlleleCallEvaluator invalid path
ALLELECALL_EVALUATOR_INVALID_PATH = ['chewBBACA.py', 'AlleleCallEvaluator',
						 			 '-i', 'data/allelecall_data/fake_results',
                					 '-g', 'data/allelecall_data/sagalactiae_schema',
                					 '-o', 'results_report']

# AlleleCallEvaluator valid input
ALLELECALL_EVALUATOR_VALID = ['chewBBACA.py', 'AlleleCallEvaluator',
                			  '-i', 'data/allelecall_data/test_results/mode4',
                			  '-g', 'data/allelecall_data/sagalactiae_schema',
                			  '-o', 'results_report']

# PrepExternalSchema
# PrepExternalSchema empty directory
PREPEXTERNALSCHEMA_TEST_EMPTY_DIR = ['chewBBACA.py', 'PrepExternalSchema',
									 '-g', 'empty_dir',
									 '-o', 'adapted_schema']

# PrepExternalSchema empty input files
PREPEXTERNALSCHEMA_TEST_EMPTY_FILES = ['chewBBACA.py', 'PrepExternalSchema',
									   '-g', 'data/prep_data/empty_files',
									   '-o', 'adapted_schema']

# PrepExternalSchema files with no data
PREPEXTERNALSCHEMA_TEST_ZERO_BYTES = ['chewBBACA.py', 'PrepExternalSchema',
									  '-g', 'data/prep_data/zero_bytes_pair',
									  '-o', 'adapted_schema']

# PrepExternalSchema invalid input path
PREPEXTERNALSCHEMA_TEST_INVALID_PATH = ['chewBBACA.py', 'PrepExternalSchema',
										'-g', 'this/path/aint/real',
										'-o', 'adapted_schema']

# PrepExternalSchema valid input
PREPEXTERNALSCHEMA_TEST_VALID_INPUT = ['chewBBACA.py', 'PrepExternalSchema',
									   '-g', 'data/prep_data/valid_input',
									   '-o', 'preped_schema']

# PrepExternalSchema different file extensions
PREPEXTERNALSCHEMA_TEST_EXTENSIONS = ['chewBBACA.py', 'PrepExternalSchema',
									  '-g', 'data/prep_data/file_extensions',
									  '-o', 'preped_schema']

# PrepExternalSchema loci list
PREPEXTERNALSCHEMA_TEST_GENE_LIST = ['chewBBACA.py', 'PrepExternalSchema',
									 '-g', 'data/prep_data/file_extensions',
									 '-o', 'preped_schema',
									 '--gl', 'data/prep_data/test_genes_list/test_genes_extension.txt']

# PrepExternalSchema input name includes blank space
PREPEXTERNALSCHEMA_TEST_BLANK_SPACE = ['chewBBACA.py', 'PrepExternalSchema',
									 '-g', 'data/allelecall_data/test_invalid_input_names/blank_spaces',
									 '-o', 'preped_schema']

# PrepExternalSchema input file has unique prefix longer than 30 chars
PREPEXTERNALSCHEMA_TEST_LONG_PREFIX = ['chewBBACA.py', 'PrepExternalSchema',
									 '-g', 'data/allelecall_data/test_invalid_input_names/long_prefix',
									 '-o', 'preped_schema']

# PrepExternalSchema some input files have the same prefix
PREPEXTERNALSCHEMA_TEST_SAME_PREFIX = ['chewBBACA.py', 'PrepExternalSchema',
									 '-g', 'data/allelecall_data/test_invalid_input_names/same_prefix',
									 '-o', 'preped_schema']

# PrepExternalSchema input prefix interpreted as PDB chain ID
PREPEXTERNALSCHEMA_TEST_PDB_CHAIN = ['chewBBACA.py', 'PrepExternalSchema',
									 '-g', 'data/allelecall_data/test_invalid_input_names/pdb_prefix',
									 '-o', 'preped_schema']

# SchemaEvaluator
# SchemaEvaluator empty input files
SCHEMAEVALUATOR_TEST_EMPTY_FILES = ['chewBBACA.py', 'SchemaEvaluator',
									'-g', 'data/schemaevaluator_data/empty_files',
									'-o', 'schema_report']

# SchemaEvaluator files with no data
SCHEMAEVALUATOR_TEST_ZERO_BYTES = ['chewBBACA.py', 'SchemaEvaluator',
								   '-g', 'data/schemaevaluator_data/zero_bytes_pair',
								   '-o', 'schema_report']

# SchemaEvaluator invalid input path
SCHEMAEVALUATOR_TEST_FAKE_PATH = ['chewBBACA.py', 'SchemaEvaluator',
								  '-g', 'this/path/aint/real',
								  '-o', 'schema_report']

# SchemaEvaluator valid input
SCHEMAEVALUATOR_TEST_VALID_INPUT = ['chewBBACA.py', 'SchemaEvaluator',
                					'-g', 'data/schemaevaluator_data/test_schema',
									'-o', 'schema_report',
									'--loci-reports', '--add-sequences']

# SchemaEvaluator single FASTA with single allele
SCHEMAEVALUATOR_TEST_SINGLE_ALLELE = ['chewBBACA.py', 'SchemaEvaluator',
									  '-g', 'data/schemaevaluator_data/single_allele',
									  '-o', 'schema_report',
									  '--loci-reports', '--add-sequences']

# SchemaEvaluator single FASTA with single invalid allele
SCHEMAEVALUATOR_TEST_SINGLE_INVALID_ALLELE = ['chewBBACA.py', 'SchemaEvaluator',
											  '-g', 'data/schemaevaluator_data/single_invalid_allele',
											  '-o', 'schema_report',
											  '--loci-reports', '--add-sequences']

# SchemaEvaluator several invalid alleles
SCHEMAEVALUATOR_TEST_SEVERAL_INVALID_ALLELES =  ['chewBBACA.py', 'SchemaEvaluator',
                								 '-g', 'data/schemaevaluator_data/several_invalid_alleles',
												 '-o', 'schema_report',
												 '--loci-reports', '--add-sequences']

# SubsetResults
# SubsetResults valid input - loci list
SUBSETRESULTS_TEST_LOCI_LIST = ['chewBBACA.py', 'SubsetResults',
                				'-i', 'data/subsetresults_data/input_data/valid_input',
								'-o', 'test_results',
								'-l', 'data/subsetresults_data/loci.txt']

# SubsetResults valid input - sample list
SUBSETRESULTS_TEST_SAMPLE_LIST = ['chewBBACA.py', 'SubsetResults',
                				  '-i', 'data/subsetresults_data/input_data/valid_input',
								  '-o', 'test_results',
								  '-s', 'data/subsetresults_data/samples.txt']

# SubsetResults valid input - loci and sample lists
SUBSETRESULTS_TEST_LOCI_SAMPLE_LIST = ['chewBBACA.py', 'SubsetResults',
                				  	   '-i', 'data/subsetresults_data/input_data/valid_input',
								  	   '-o', 'test_results',
								  	   '-l', 'data/subsetresults_data/loci.txt',
									   '-s', 'data/subsetresults_data/samples.txt']

# SubsetResults no loci and sample lists
SUBSETRESULTS_TEST_MISSING_LISTS = ['chewBBACA.py', 'SubsetResults',
                				  	'-i', 'data/subsetresults_data/input_data/valid_input',
								  	'-o', 'test_results']

# SubsetResults missing all results files
SUBSETRESULTS_TEST_MISSING_FILES = ['chewBBACA.py', 'SubsetResults',
                				  	'-i', 'data/subsetresults_data/input_data/missing_files',
								  	'-o', 'test_results'
									'-l', 'data/subsetresults_data/loci.txt',
									'-s', 'data/subsetresults_data/samples.txt']

# SubsetResults missing profiles
SUBSETRESULTS_TEST_MISSING_PROFILES = ['chewBBACA.py', 'SubsetResults',
                				  	   '-i', 'data/subsetresults_data/input_data/missing_profiles',
								  	   '-o', 'test_results'
									   '-l', 'data/subsetresults_data/loci.txt',
									   '-s', 'data/subsetresults_data/samples.txt']

# SubsetResults subset loci not in dataset
SUBSETRESULTS_TEST_ABSENT_LOCI = ['chewBBACA.py', 'SubsetResults',
                				  '-i', 'data/subsetresults_data/input_data/missing_profiles',
								  '-o', 'test_results'
								  '-l', 'data/subsetresults_data/loci_absent.txt']

# SubsetResults subset samples not in dataset
SUBSETRESULTS_TEST_ABSENT_SAMPLES = ['chewBBACA.py', 'SubsetResults',
                				     '-i', 'data/subsetresults_data/input_data/missing_profiles',
								     '-o', 'test_results'
								     '-s', 'data/subsetresults_data/samples_absent.txt']

# MergeResults
# MergeResults valid input - default
MERGERESULTS_TEST_DEFAULT = ['chewBBACA.py', 'MergeResults',
							 '-i', 'data/mergeresults_data/input_data/default_input1', 'data/mergeresults_data/input_data/default_input2',
							 '-o', 'test_default_results']

# MergeResults valid input - common
MERGERESULTS_TEST_COMMON = ['chewBBACA.py', 'MergeResults',
							'-i', 'data/mergeresults_data/input_data/common_input1', 'data/mergeresults_data/input_data/common_input2',
							'-o', 'test_common_results',
							'--common']

# MergeResults invalid input - input file
MERGERESULTS_TEST_INPUT_FILE = ['chewBBACA.py', 'MergeResults',
								'-i', 'data/mergeresults_data/input_data/default_input1', 'data/mergeresults_data/input_data/input_file.txt',
								'-o', 'test_input_file_results']

# MergeResults invalid input - no common loci
MERGERESULTS_TEST_NOCOMMONLOCI = ['chewBBACA.py', 'MergeResults',
								  '-i', 'data/mergeresults_data/input_data/nocommon_input1', 'data/mergeresults_data/input_data/nocommon_input2',
								  '-o', 'test_nocommonloci_results']

# MergeResults invalid input - different sets of loci
MERGERESULTS_TEST_LOCIDIFFER = ['chewBBACA.py', 'MergeResults',
								'-i', 'data/mergeresults_data/input_data/locidiffer_input1', 'data/mergeresults_data/input_data/locidiffer_input2',
								'-o', 'test_locidiffer_results']

# HashProfiles
# HashProfiles valid input
HASHPROFILES_TEST_VALID = ['chewBBACA.py', 'HashProfiles',
						   '-i', 'data/hashprofiles_data/input_data/results_alleles.tsv',
						   '-g', 'data/hashprofiles_data/input_data/sagalactiae_schema',
						   '-o', 'test_hashprofiles_results']

# HashProfiles invalid hashign algorithm
HASHPROFILES_TEST_INVALID_HASH = ['chewBBACA.py', 'HashProfiles',
						   		  '-i', 'data/hashprofiles_data/input_data/results_alleles.tsv',
						   		  '-g', 'data/hashprofiles_data/input_data/sagalactiae_schema',
						   		  '-o', 'test_hashprofiles_results',
						   		  ' --hash-type', 'invalid']

# ComputeDistances
# ComputeDistances valid inputs

# ComputeDistances default
COMPUTEDISTANCES_TEST_DEFAULT = ['chewBBACA.py', 'ComputeDistances',
								 '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
								 '-o', 'data/computedistances_data/expected_results/default']

# ComputeDistances --outfmt symmetric
COMPUTEDISTANCES_TEST_DEFAULT_SYMMETRIC = ['chewBBACA.py', 'ComputeDistances',
										   '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
										   '-o', 'data/computedistances_data/expected_results/default',
										   '--outfmt', 'symmetric']

# ComputeDistances --similarity
COMPUTEDISTANCES_TEST_DEFAULT_SIMILARITY = ['chewBBACA.py', 'ComputeDistances',
											'-i', 'data/computedistances_data/input_data/results_alleles.tsv',
											'-o', 'data/computedistances_data/expected_results/default',
											'--similarity']

# ComputeDistances --outfmt symmetric and --similarity
COMPUTEDISTANCES_TEST_DEFAULT_SYMMETRIC_SIMILARITY = ['chewBBACA.py', 'ComputeDistances',
										   			  '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
										   			  '-o', 'data/computedistances_data/expected_results/default',
										   			  '--outfmt', 'symmetric',
										   			  '--similarity']

# ComputeDistances --outfmt table
COMPUTEDISTANCES_TEST_DEFAULT_TABLE = ['chewBBACA.py', 'ComputeDistances',
									   '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
									   '-o', 'data/computedistances_data/expected_results/default',
									   '--outfmt', 'table']

# ComputeDistances --m jaccard and --outfmt symmetric
COMPUTEDISTANCES_TEST_JACCARD_SYMMETRIC = ['chewBBACA.py', 'ComputeDistances',
										   '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
										   '-o', 'data/computedistances_data/expected_results/default',
										   '--m', 'jaccard',
										   '--outfmt', 'symmetric']

# ComputeDistances --m jaccard, --outfmt symmetric and --similarity
COMPUTEDISTANCES_TEST_JACCARD_SYMMETRIC_SIMILARITY = ['chewBBACA.py', 'ComputeDistances',
										   			  '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
										   			  '-o', 'data/computedistances_data/expected_results/default',
										   			  '--m', 'jaccard',
										   			  '--outfmt', 'symmetric',
										   			  '--similarity']

# ComputeDistances --m loci and --outfmt symmetric
COMPUTEDISTANCES_TEST_LOCI_SYMMETRIC = ['chewBBACA.py', 'ComputeDistances',
										'-i', 'data/computedistances_data/input_data/results_alleles.tsv',
										'-o', 'data/computedistances_data/expected_results/default',
										'--m', 'loci',
										'--outfmt', 'symmetric']

# ComputeDistances --m loci, --outfmt symmetric and --similarity
COMPUTEDISTANCES_TEST_LOCI_SYMMETRIC_SIMILARITY = ['chewBBACA.py', 'ComputeDistances',
												   '-i', 'data/computedistances_data/input_data/results_alleles.tsv',
												   '-o', 'data/computedistances_data/expected_results/default',
												   '--m', 'loci',
												   '--outfmt', 'symmetric',
												   '--similarity']
