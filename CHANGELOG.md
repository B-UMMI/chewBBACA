# Changelog

## 3.3.0 - 2023-10-11

- Added the AlleleCallEvaluator module. This module generates an interactive HTML report for the allele calling results. The report provides summary statistics to evaluate results per sample and per locus (with the possibility to provide a TSV file with loci annotations to include on a table). The report includes components to display a heatmap representing the loci presence-absence matrix, a heatmap representing the distance matrix based on allelic differences and a Neighbor-Joining tree based on the MSA of the core genome loci.

- Added [pyrodigal](https://github.com/althonos/pyrodigal) for gene prediction. This simplified the processing of the gene prediction results and reduced runtime.

- Fixed an issue where the AlleleCall module would try to create results files for excluded inputs.

- Fixed exception capturing during multiprocessing when using Python>=3.11.

- Fixed PLOT5/3 identification when coding sequences are in the reverse strand.

- Fixed computation of the representative self-scores when performing allele calling for a subset of the loci in a schema (would only compute the self-scores for the subset of loci if the 'self_scores' file had still not been created).

- Fixed issue related to the classification of single EXC/INF and single/multiple ASM/ALM (would classify some inputs as NIPH instead of EXC/INF).

- Fixed issue related to protein exact match classification when multiple pre-computed PROTEINtable files include the same protein hash.

- Changed the `-i`, `--input-files` parameter in the PrepExternalSchema and UniprotFinder modules to `-g`, `--schema-directory` and added the `--gl`, `--genes-list` parameter to enable adapting or annotating a subset of the loci in the schema.

## 3.2.0 - 2023-04-27

- New version of the SchemaEvaluator module. The updated version fixes several issues related with outdated dependencies that were leading to errors in the previous version. The new version also includes new features and components. Read the [docs page](https://chewbbaca.readthedocs.io/en/latest/user/modules/SchemaEvaluator.html) to know more about the latest version of the SchemaEvaluator module.
- Updated the link to the UniProt FTP used by the UniprotFinder module. 
- Added the `.fas` file extension to the list of file extensions accepted by chewBBACA. chewBBACA accepts genome assemblies and external schemas with FASTA files that use any of the following file extensions: `.fasta`, `.fna`, `.ffn`, `.fa` and `.fas`. The FASTA files created by chewBBACA use the `.fasta` extension.
- Fixed issue in the PrepExternalSchema module where it would only detect FASTA files if they ended with the `.fasta` extension.
- Added the `--size-filter` parameter to the PrepExternalSchema module to define if the adaptation process should filter out alleles based on the minimum length and size threshold values.
- Added the `--output-novel` parameter to the AlleleCall module. If this parameter is used, the AlleleCall module creates a FASTA file with the novel alleles inferred during the allele calling. This file is created even if the `--no-inferred` parameter is used and the novel alleles are not added to the schema.

## 3.1.2 - 2023-03-13

- Fixed issue related with sequence header format when FASTA files with coding sequences (CDSs) were provided to the CreateSchema and AlleleCall modules through the `--cds-input` parameter. chewBBACA expected the sequence headers to have the same format used to save the CDSs extracted during the gene prediction step, `<genomeID>-protein<cdsID>`. This would lead to errors if the FASTA files with CDSs provided by users had a different sequence header format. To fix this issue, chewBBACA determines the unique genome identifier for each input/genome and renames the sequence headers to match the `<genomeID>-protein<cdsID>` format (e.g: given a file named `GCF_000007125.1_ASM712v1_cds_from_genomic.fna`, the unique genome identifier is `GCF_000007125`, and the sequence headers are renamed to `GCF_000007125-protein1`, `GCF_000007125-protein2`, ..., `GCF_000007125-proteinN`). FASTA files with renamed headers are stored in the temporary directory created by chewBBACA. The sequence headers are renamed sequentially, so the integer after `protein` indicates the order of the sequences in the original FASTA file

## 3.1.1 - 2023-03-09

- Fixed issue in the ExtractCgMLST module. The cgMLST matrix only contained 1's and 0's because it was a subset of the presence-absence matrix. The module now uses the list of loci in the core genome to subset the masked matrix with the allele identifiers.

## 3.1.0 - 2023-01-12

- Updated the ExtractCgMLST module ([PR](https://github.com/B-UMMI/chewBBACA/pull/153)). It can now accept several loci presence threshold values and creates a HTML file with a line plot that displays the number of loci in the cgMLST per threshold value (default: compute for 0.95, 0.99 and 1 thresholds).
- Removed the TestGenomeQuality module (superseded by the new functionalities implemented in the ExtractCgMLST module).
- Bugfix for the determination of the BLASTp raw scores for the representative alleles ([PR](https://github.com/B-UMMI/chewBBACA/pull/156)).
- Changed the valued passed to BLASTp `-max_target_seqs` parameter during representative determination from 20 to the default value used by BLASTp (500). This leads to a slight improvement in classification accuracy.

## 3.0.0 - 2022-12-14

New implementation of the **AlleleCall** process. The new implementation was developed to reduce execution time, improve accuracy and provide more detailed results. It uses available computational resources more efficiently to allow for analyses with thousands of strains in a laptop. This new version is fully compatible with schemas created with previous versions.

### AlleleCall changes

- The new implementation avoids redundant comparisons through the identification of the set of distinct CDSs in the input files. The classification for a distinct CDS is propagated to classify all input genomes that contain the CDS.
- Implemented a clustering step based on minimizers to cluster the translated CDSs. This step complements the alignment-based strategy with BLASTp to increase computational efficiency and classification accuracy.
- The AlleleCall process has 4 execution modes (1: only exact matches at DNA level; 2: exact matches at DNA and Protein level; 3: exact matches and minimizer-based clustering to find similar alleles with BSR > 0.7; 4: runs the full process to find exact matches and all matches with BSR >= 0.6).
- Files with information about loci length modes (`loci_modes`) and the self-alignment raw score for the representative alleles (`short/self_scores`) are pre-computed and automatically updated (the process no longer creates and updates a file with the self-alignment raw score per locus).
- The process creates the `pre_computed` folder to store files with hash tables that are used to speedup exact matching and avoid running the step to translate the schema alleles in every run.
- Added the `--cds` parameter to accept FASTA files with CDSs (one FASTA file per genome) and skip gene prediction with Prodigal.
- Users can control the addition of novel alleles to the schema with the `--no-inferred` parameter.
- Added the `--output-unclassified` parameter to write a FASTA file (`unclassified_sequences.fasta`) with the distinct CDSs that were not classified in a run.
- Added the `--output-missing` parameter to write a FASTA file (`missing_classes.fasta`) and a TSV file with information about the classified sequences that led to a locus being classified as ASM, ALM, PLOT3, PLOT5, LOTSC, NIPH, NIPHEM and PAMA.
- Added the `--no-cleanup` parameter to keep the temporary folder with intermediate files created during a run.
- Removed the `--contained`, `--force-reset`, `--store-profiles` (to be reimplemented in a future release), `--json` and `--verbose` parameters.
- The `--force-continue` parameter no longer allows users to continue a run that was interrupted. This parameter is now used to ignore warnings and prompts about missing configuration files and the usage of multiple argument values per parameter.
- The allelic profiles in the `results_alleles.tsv` file can be hashed by providing the `--hash-profiles` parameter and a valid hash type as argument (hash algorithms available from the [hashlib](python) library and crc32 and adler32 from the [zlib](https://docs.python.org/3/library/zlib.html) library).
- The process creates a TSV file, `cds_coordinates.tsv`, with the genomic coordinates for all CDSs identified in the input files.
- The process creates a TSV file, `loci_summary_stats.tsv`, with summary statistics for loci classifications.
- The process no longer creates the `RepeatedLoci.txt` file. It now creates the `paralogous_counts.tsv` and `paralogous_loci.tsv` files with more detailed information about the loci identified as paralogous.
- The PLNF class is attributed in modes 1, 2 and 3 to indicate that a more thorough analysis might have found a match for the loci that were not found (LNF).
- CDSs that match several loci are classified as PAMA.
- Bugfix for PLOT3, PLOT5 and LOTSC classification types. LOTSC classification was not always attributed when a contig was smaller than the matched representative allele and some PLOT5 cases were classified as LOTSC. LOTSC cases counted as exact matches in the `results_statistics.tsv` file.

### Additional changes

- The UniprotFinder allows users to search for annotations through UniProt's SPARQL endpoint or based on matches against UniProt's reference proteomes or both.
- Bugfix for an issue in the UniprotFinder module that was leading to errors when the data returned by UniProt's SPARQL endpoint only contained one set of annotation terms.
- Bugfix for an issue in the UniprotFinder module that was preventing the annotations from being written to the output file.
- Bugfix for an issue in the [map_async_parallelizer](https://github.com/B-UMMI/chewBBACA/blob/d7572c085677319500546dbb4ed8eee69cc3d2c2/CHEWBBACA/utils/multiprocessing_operations.py#L51) function that led to high memory usage.
- Implemented and changed several functions in the modules included in the `utils` folder to optimize code reusability, reduce runtime and peak memory usage, especially for large schemas and datasets (these changes affect mostly the CreateSchema and AlleleCall modules).
- Updated function docstrings and added comments.

## 2.8.5 - 2021-07-05

- Updated **JoinProfiles** process to accept any number of files to join. Added `--common` flag to identify and join results for the set of loci common to all input files.
- Added line breaks at the end of file for **AlleleCall** outputs.
- Improved code readability and efficiency for the **RemoveGenes** process.

## 2.8.4 - 2021-04-27

- Fixed issue related with absence of UniprotFinder module from setup.py.

## 2.8.3 - 2021-04-23

Added the Sequence Logo component to the individual reports of SchemaEvaluator.

## 2.8.2 - 2021-04-22

- Bugfix to fix issue when users provided a text file with the paths to the input genome assemblies or genes to use in the CreateSchema or AlleleCall processes.

## 2.8.1 - 2021-04-08

- Bugfix to fix issue in the PrepExternalSchema process when BLAST >= v2.10 was used.

## 2.8.0 - 2021-03-30

New implementation of the **UniprotFinder** process. A new feature allows users to specify a set of taxa and determine annotation terms based on matches against UniProt's reference proteomes for those taxa.

## 2.7.0 - 2021-03-10

New implementation of the **CreateSchema** process. This new implementation significantly reduces execution time. It is designed to enable schema creation based on hundreds or thousands of assemblies on a laptop. The schemas generated by the new implementation are fully compatible with previous versions.

### Additional changes

- Improved detection of invalid inputs (inputs that do not contain coding sequences (CDSs), that contain invalid sequences/characters, empty files, etc).
- New parameter `--pm` allows users to set Prodigal's execution mode. The `single` mode is the default mode. Use the `meta` mode for input files that have less than 100kbp (e.g.: plasmids, viruses).
- `CreateSchema` accepts a single or several FASTA files with CDSs if the `--CDS` option is included in the command. This option skips the gene prediction step with Prodigal and creates a schema seed based on the CDSs in the input files.
- `AlleleCall` can automatically detect parameter values previously used with a schema. Users only need to provide values for the `-i`, `-g` and `-o` parameters.

## 2.6.0 - 2021-03-09

The **SchemaEvaluator** module was refactored due to unsupported dependencies that were not allowing the report generation. The style of the report is similar to the what can be found on [Chewie-NS](https://chewbbaca.online/).

More in-depth information can be found about the module on its [wiki page](https://github.com/B-UMMI/chewBBACA/wiki/4.-Schema-Evaluation).

## 2.5.6 - 2020-11-10

- Fixed BLAST version detection. BLAST versions greater than 2.9 (e.g.: 2.10) were not correctly detected. Thanks to [Eric DEVEAUD](https://github.com/EricDeveaud) for creating a pull request with a fix for this issue!
- Fixed issue that would lead to error when users provided a file with the list of genes to the AlleleCall process
- Added function that uses the subprocess module to run BLAST and capture warnings raised during normal execution of BLAST >= 2.10.

## 2.5.5 - 2020-09-15

- Removed Bio.Alphabet imports and generic_dna mentions. New version of the BioPython package was not compatible with the way chewBBACA was using those features.

## 2.5.4 - 2020-08-10

- Corrected problem related with versioning inconsistencies that would lead to errors during validation steps.

## 2.5.3 - 2020-08-10

- Organized argparsing for several processes and fixed argparsing issues introduced in version 2.5.0 for the `TestGenomeQuality`, `SchemaEvaluator` and `UniprotFinder` processes.

## 2.5.2 - 2020-08-08

- Removed overconservative constraints leading to compatibility issues between different chewBBACA versions.
- Implemented check to verify if Users have authorization to submit novel alleles during the SyncSchema process.
- Corrected argparsing for the `RemoveGenes` process.
- Corrected issue that would not let Users properly Sync Tutorial schemas.

## 2.5.1 - 2020-08-05

- Chewie verifies if Chewie-NS instance is available before starting processes.
- SyncSchema gets Chewie-NS base URL from configuration file.
- Schemas can be created without a Prodigal training file and without a size threshold value.
- LoadSchema and SyncSchema processes display info about loci and allele insertion progress.
- Corrected bug leading to errors during the AlleleCall process. Schema adaptation during the SyncSchema process may change loci representatives and pre-computed BSR values would become outdated.
- Corrected bug leading to error during the update of allelic profiles if new allele identifiers contained '*'.

## 2.5.0 - 2020-06-30

We've developed [Chewie-NS](https://chewbbaca.online/), a Nomenclature Server that is based on the [TypOn](https://jbiomedsem.biomedcentral.com/articles/10.1186/2041-1480-5-43) ontology and integrates with chewBBACA to provide access to gene-by-gene typing schemas and to allow a common and global allelic nomenclature to be maintained.

To allow all users to interact with the Chewie-NS, we've implemented the following set of modules:

- `LoadSchema`: enables upload of new schemas to the Chewie-NS.
- `DownloadSchema`: enables download of any schema from the Chewie-NS.
- `SyncSchema`: compares local schemas, previously downloaded from the Chewie-NS, with the remote versions in the Chewie-NS to download and add new alleles to local schemas, submit new alleles to update remote schemas and ensure that a common allele identifier nomenclature is maintained.
- `NSStats`:  retrieves basic information about species and schemas in the Chewie-NS.

The [documentation](https://chewie-ns.readthedocs.io/en/latest/) includes information about the integration with chewBBACA and how to run the new [LoadSchema](https://chewie-ns.readthedocs.io/en/latest/user/upload_api.html), [DownloadSchema](https://chewie-ns.readthedocs.io/en/latest/user/download_api.html), [SyncSchema](https://chewie-ns.readthedocs.io/en/latest/user/synchronize_api.html) and [NSStats](https://chewie-ns.readthedocs.io/en/latest/user/nsstats_api.html) processes.
The Chewie-NS [source code](https://github.com/B-UMMI/Nomenclature_Server_docker_compose) is freely available and deployment of local instances can be easily achieved through Docker Compose.

This version also includes other changes:

- The `AlleleCall` process will detect if a schema was created with previous chewBBACA versions and ask users if they wish to convert the schema to the latest version. The conversion process **will not alter** your schema files, it will simply add configuration files and copy the Prodigal training file to the schema's directory. You can force schema conversion with the `--fc` argument.
- It is now required that users provide a valid Prodigal training file to the `CreateSchema` and `PrepExternalSchema` processes. The training file will be included in the schema and can be automatically detected by the `AlleleCall` process.
- Schemas created with the `CreateSchema` process or adapted with the `PrepExternalSchema` retain information about parameters values (BLAST Score Ratio, Prodigal training file, genetic code, minimum sequence length and sequence size variation threshold) and users are advised to keep performing allele call with those parameters values to ensure consistent results and provide the possibility of schema upload to the Chewie-NS. The AlleleCall process detects if a user provides parameters values that differ from the original values and requests confirmation before proceeding (you may force execution with the `--fc` argument).
- The AlleleCall process creates a SQLite database in the schema's directory that is used to store the allelic profiles determined with that schema.
- Further optimizations in the `PrepExternalSchema` process.

## 2.1.0 - 2019-11-05

- New `PrepExternalSchema` implementation: algorithmic optimizations to improve speed and maintain memory efficiency. New output files with summary information about schema adaptation and excluded sequences and options to control the Blast Score Ratio, minimum sequence length and genetic code values passed to the process.

## chewBBACA released as a galaxy module!

Many Thanks to Stefano Morabito and Arnold Knijn (https://github.com/aknijn) for EURL VTEC in ISS, Rome! 
https://toolshed.g2.bx.psu.edu/repository?repository_id=88fd7663075eeae9&changeset_revision=093352878303

## 2.0.17 - 2019-02-10

- New alleles also have a timestamp added to the allele name.

## 2.0.16 - 2018-01-06

- Corrected bug from 2.0.15 when no Prodigal training file was provided.

## 2.0.15 - 2018-01-06

- Added Prodigal training files to the package. Available at `CHEWBBACA/prodigal_training_files`.
 
## 2.0.13 - 2018-09-18

- When using the function `PrepExternalSchema`, older behavior would remove any locus with a single translation error while the latest change (2.0.12) would not change the original source fasta, which would make the schema unusable. It is now enforced that the alleles that do not translate are removed from the fasta, be sure to backup your data before using this function.
 
## 2.0.11 - 2018-06-05

- Corrected bug when `-h` parameter was provided for allele call.
- New option for the schema creation. A schema can be created based on a single fasta file, skipping the gene prediction step. Use `--CDS` and provide a single FASTA file to the `-i` parameter.
 
## 2.0.10 - 2018-05-21

- cgMLST profile extraction function (ExtractCgMLST) more efficient (thanks Dillon Barker).
- New option for the allele call, size threshold previously hardcoded at 0.2 can now be changed using the `--st` option. Size threshold is important for the definition of ASM and ALM (alleles smaller/larger than mode).
 
## 2.0.9 - 2018-04-04

- Blast results during allele call are not saved as a file, instead are piped directly for processing.
- New option for the allele call, if genome fasta input is already a fasta of CDS use the `--CDS` option.
 
## 2.0.7 - 2018-02-23

- Corrected bug that prevented usage of latest blast version (>=2.7.0).
- Version flag can now be used `--version`.
- Instead of calling the main script `chewBBACA.py` you can now use `chewie` (if installed trought pip).
 
## 2.0.5 - 2018-02-18

- AlleleCall: `-i` parameter accepts a single fasta file now.
