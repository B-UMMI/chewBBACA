# Changelog

## 3.0.0 - 2022-08-25

New implementation of the **AlleleCall** process. The new implementation was developed to reduce execution time, improve accuracy and provide more detailed results. It uses available computational resources more efficiently to allow for analyses with thousands of strains in a laptop. The main steps of the new implementation are the following:

- Gene prediction with Prodigal followed by coding sequence (CDS) extraction to create FASTA files that contain all CDS extracted from the inputs. CDS coordinates are saved to pickle files.
- CDS deduplication to identify the distinct set of CDS and keep information about the inputs that contain each distinct CDS (hashtable with mapping between CDS SHA-256 and list of unique integer identifiers for the inputs that contain each CDS compressed with [polyline encoding](https://developers.google.com/maps/documentation/utilities/polylinealgorithm) adapted from [numcompress](https://github.com/amit1rrr/numcompress)).
- Exact match between distinct CDS and schema alleles at DNA level. Information about CDS classified at this stage is stored in pickle files that are updated throughout the process.
- Translation of distinct CDS that were not an exact match in the previous step. This step excludes CDS with ambiguous bases and CDS that are below the minimum length value defined.
- Translated CDS deduplication to identify the distinct set of translated CDS and keep information about the inputs that contain CDS that encode each distinct protein (hashtable with mapping between translated CDS SHA-256 and list of unique integer identifiers for the distinct CDS encoded with polyline encoding).
- Minimizer-based clustering (interior minimizers selected based on lexicographic order, k=5, w=5) to find translated CDS similar to schema representative alleles.
- Alignment, using BLASTp, of each cluster representative against all CDS added to its cluster to identify and classify matches with a BLAST Score Ratio (BSR) > 0.7.
- Alignment, using BLASTp, of each schema representative against the remaining CDS that were not classified in the previous steps to find matches 0.6 >= BSR <= 0.7 that are evaluated to determine new representative alleles.
- Assessment of the classifications for each locus to attribute allele identifiers and other classification types (ASM, ALM, PLOT3, PLOT5, LOTSC, NIPH, NIPHEM and LNF).

### Additional changes

- Users can control the addition of novel alleles to the schema with the `--no-inferred` parameter.
- The AlleleCall process has 4 execution modes (1: only exact matches at DNA level; 2: exact matches at DNA and Protein level; 3: exact matches and minimizer-based clustering to find similar alleles with BSR > 0.7; 4: runs the full process to find exact matches and all matches with BSR >= 0.6).
- The allelic profiles in the `results_alleles.tsv` file can be hashed by providing the `--hash-profiles` parameter and a valid hash type as argument (hash algorithms available from the [hashlib](python) library and crc32 and adler32 from the [zlib](https://docs.python.org/3/library/zlib.html) library).
- The UniprotFinder allows users to search for annotations through UniProt's SPARQL endpoint or based on matches against UniProt's reference proteomes or both.
- Bugfix for an issue in the UniprotFinder module that was leading to errors when the data returned by UniProt's SPARQL endpoint only contained one set of annotation terms).
- Bugfix for an issue in the [map_async_parallelizer](https://github.com/B-UMMI/chewBBACA/blob/d7572c085677319500546dbb4ed8eee69cc3d2c2/CHEWBBACA/utils/multiprocessing_operations.py#L51) function that led to high memory usage.
- Bugfix for PLOT3, PLOT5 and LOTSC classification types. LOTSC classification was not always attributed when a contig was smaller than the matched representative allele and some PLOT5 cases were classified as LOTSC. LOTSC cases counted as exact matches in the `results_statistics.tsv` file.

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

* Fixed BLAST version detection. BLAST versions greater than 2.9 (e.g.: 2.10) were not correctly detected. Thanks to [Eric DEVEAUD](https://github.com/EricDeveaud) for creating a pull request with a fix for this issue!
* Fixed issue that would lead to error when users provided a file with the list of genes to the AlleleCall process
* Added function that uses the subprocess module to run BLAST and capture warnings raised during normal execution of BLAST >= 2.10.

## 2.5.5 - 2020-09-15

* Removed Bio.Alphabet imports and generic_dna mentions. New version of the BioPython package was not compatible with the way chewBBACA was using those features.

## 2.5.4 - 2020-08-10

* Corrected problem related with versioning inconsistencies that would lead to errors during validation steps.

## 2.5.3 - 2020-08-10

*  Organized argparsing for several processes and fixed argparsing issues introduced in version 2.5.0 for the `TestGenomeQuality`, `SchemaEvaluator` and `UniprotFinder` processes.

## 2.5.2 - 2020-08-08

* Removed overconservative constraints leading to compatibility issues between different chewBBACA versions.
* Implemented check to verify if Users have authorization to submit novel alleles during the SyncSchema process.
* Corrected argparsing for the `RemoveGenes` process.
* Corrected issue that would not let Users properly Sync Tutorial schemas.

## 2.5.1 - 2020-08-05

* Chewie verifies if Chewie-NS instance is available before starting processes.
* SyncSchema gets Chewie-NS base URL from configuration file.
* Schemas can be created without a Prodigal training file and without a size threshold value.
* LoadSchema and SyncSchema processes display info about loci and allele insertion progress.
* Corrected bug leading to errors during the AlleleCall process. Schema adaptation during the SyncSchema process may change loci representatives and pre-computed BSR values would become outdated.
* Corrected bug leading to error during the update of allelic profiles if new allele identifiers contained '*'.

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

* New `PrepExternalSchema` implementation: Algorithmic optimizations to improve speed and maintain memory efficiency. New output files with summary information about schema adaptation and excluded sequences and options to control the Blast Score Ratio, minimum sequence length and genetic code values passed to the process.

## chewBBACA released as a galaxy module!

Many Thanks to Stefano Morabito and Arnold Knijn (https://github.com/aknijn) for EURL VTEC in ISS, Rome ! 
https://toolshed.g2.bx.psu.edu/repository?repository_id=88fd7663075eeae9&changeset_revision=093352878303

## 2.0.17 - 2019-02-10

* New alleles also have a timestamp added to the allele name.

## 2.0.16 - 2018-01-06

* Corrected bug from 2.0.15 when no prodigal training file provided.

## 2.0.15 - 2018-01-06

* Added prodigal training files to the package. They are now at CHEWBBACA/prodigal_training_files.
 
## 2.0.13 - 2018-09-18

* when using the function `PrepExternalSchema`, older behavior would remove any locus with a single translation error while the latest change(2.0.12) would not change the original source fasta, this would make the schema unusable. It is now enforced that the alleles that do not translate are removed from the fasta, be sure to backup your data before using this function.
 
## 2.0.11 - 2018-06-05

* corrected bug when -h on allele call
* new option for the schema creation. A schema can be created based on a single fasta file, jumping the prodigal gene prediction running. Use `--CDS` and provide a sinfle fasta file on the `-i` input.
 
## 2.0.10 - 2018-05-21

* cgMLST profile extraction function (ExtractCgMLST) more efficient (thanks Dillon Barker)
* new option for the allele call, size threshold previously hardcoded at 0.2 can now be changed using the `--st` option. Size threshold is important for the definition of ASM and ALM (alleles smaller/larger than mode).
 
## 2.0.9 - 2018-04-04

* blast results during allele call are not saved as a file, instead are piped directly for processing
* new option for the allele call, if genome fasta input is already a fasta of CDS use the `--CDS` option  
 
## 2.0.7 - 2018-02-23

* corrected bug that prevented usage of latest blast version (>=2.7.0)
* version flag can now be used `--version`
* instead of calling the main script `chewBBACA.py` you can now use `chewie` (if installed trought pip).
 
## 2.0.5 - 2018-02-18

* AlleleCall : -i option accepts a single fasta file now
