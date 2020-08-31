# Changelog

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

The [documentation](https://chewie-ns.readthedocs.io/en/latest/) includes information about the integration with chewBBACA and how to run the new [LoadSchema](https://chewie-ns.readthedocs.io/en/latest/user/upload_api.html), [DownloadSchema](https://chewie-ns.readthedocs.io/en/latest/user/download_api.html), [SyncSchema](https://chewie-ns.readthedocs.io/en/latest/user/synchronize_api.html) and [NSStats]() processes.
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