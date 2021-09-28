
[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/chewBBACA/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/chewbbaca)
[![chewBBACA](https://github.com/B-UMMI/chewBBACA/workflows/chewbbaca/badge.svg)](https://github.com/B-UMMI/chewBBACA/actions?query=workflow%3Achewbbaca)
[![License: GPL v3](https://img.shields.io/github/license/B-UMMI/chewBBACA)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:10.1099/mgen.0.000166](https://img.shields.io/badge/DOI-10.1099%2Fmgen.0.000166-blue)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166)

# chewBBACA

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be thought of as "Comprehensive and  Highly Efficient Workflow" 
but at this point still it needs a bit of work to make that claim so we just add "chew" to add extra coolness to the software name. BSR stands for 
BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2) 

chewBBACA is a comprehensive pipeline including a set of functions for the creation and validation of whole genome and core genome MultiLocus Sequence 
Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on Blast Score Ratio that can be run in multiprocessor 
settings and a set of functions to visualize and validate allele variation in the loci. chewBBACA performs the schema creation and allele calls on complete or draft genomes resulting from de novo assemblers.

chewBBACA has been published (version 2.0.5 at the time) in Microbial Genomics under the title:
**chewBBACA: A complete suite for gene-by-gene schema creation and strain identification**  - [Link to paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166) 

When using chewBBACA please use the following citation:

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)

## Contents

- [Important Notes](#important-notes).
- [Useful links](#useful-links).
- [Installation](#installation).
- [Quick Start](#quick-start).
- [Detailed Usage](#detailed-usage).
  1. [Whole Genome Multilocus Sequence Typing (wgMLST) schema creation](#i-whole-genome-multilocus-sequence-typing-wgmlst-schema-creation).
  2. [Allele call using a cg/wgMLST schema](#ii-allele-call-using-a-cgwgmlst-schema).
  3. [Determine annotations for loci in the schema](#iii-determine-annotations-for-loci-in-the-schema).
  4. [Evaluate wgMLST call quality per genome](#iv-evaluate-wgmlst-call-quality-per-genome).
  5. [Defining the cgMLST schema](#v-defining-the-cgmlst-schema).
  6. [Evaluate your schema](#vi-evaluate-your-schema).
  7. [Adapt an external schema](#vii-adapt-an-external-schema).
- [FAQ](#faq).
- [Citation](#citation).

## Important Notes

- chewBBACA only works with **python 3** (automatic testing for Python 3.7 and Python 3.8 with GitHub Actions).
- We strongly recommend that users install and use BLAST 2.9.0+ with chewBBACA, as chewBBACA's processes have been extensively tested with that version of BLAST.
- chewBBACA includes Prodigal training files for some species. You can consult the list of Prodigal training files that are readily available [here](https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files). We strongly recommend using the same Prodigal training file for schema creation and allele calling to ensure consistent results.
- chewBBACA defines an allele as a complete Coding DNA Sequence, with start and stop codons according 
 to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) identified using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). It will automatically exclude any allele for which the DNA sequence does not contain start or stop codons and for which the length is not multiple of three. 
- Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems, you should do a a quick conversion using for example [dos2unix](https://waterlan.home.xs4all.nl/dos2unix.html).

## Useful links

 - Check the [wiki pages](https://github.com/B-UMMI/chewBBACA/wiki) for a much more thorough chewBBACA walkthrough.
 - An extensive [tutorial repository](https://github.com/B-UMMI/chewBBACA_tutorial) is available as example on how to run an analysis pipeline using chewBBACA.
 - Look through the list of schemas in [Chewie-NS](https://chewbbaca.online/). Chewie-NS allows chewBBACA users to download and update cg/wgMLST schemas, allowing the easy sharing of results, while ensuring the reproducibility and consistency of these steps.
 - Use [BBACA gitter](https://gitter.im/BBACA/Lobby) if you have any pressing question. Chat can be faster and better than email for troubleshooting purposes.
 - A ready to use [docker image](https://hub.docker.com/r/ummidock/chewbbaca) automatically built from the latest version of chewBBACA in Ubuntu 16.04.
 - chewBBACA is available as a [Galaxy module](https://toolshed.g2.bx.psu.edu/repository?repository_id=88fd7663075eeae9&changeset_revision=093352878303) many thanks to Stefano Morabito and [Arnold Knijn](https://github.com/aknijn) for EURL VTEC in ISS, Rome!
 - Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to know more about the latest changes!

## Installation

Install using [conda](https://anaconda.org/bioconda/chewbbaca):

```
conda install -c bioconda chewbbaca
```

Install using [pip](https://pypi.org/project/chewBBACA/):

```
pip3 install chewbbaca
```

chewBBACA has the following dependencies:

Python dependencies (defined in the [requirements](https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/requirements.txt) file, should be automatically installed):
* numpy>=1.14.0
* scipy>=0.13.3
* biopython>=1.70
* plotly>=1.12.9
* SPARQLWrapper>=1.8.0
* requests>=2.2.1
* pandas>=0.22.0

Main dependencies:
* [BLAST 2.9.0+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/).
* [Prodigal 2.6.0](https://github.com/hyattpd/prodigal/releases/) or above

Other dependencies (for schema evaluation only):
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)

Installation through conda should take care of all dependencies. If you install through pip you will need to ensure that you have BLAST, Prodigal and MAFFT installed and added to the PATH.

## Quick Start

### Create a schema

Option 1 - Genome assemblies

Include all genome assemblies in a directory and adapt the following command:
```
chewBBACA.py CreateSchema -i InputAssemblies -o OutputSchemaFolder --ptf ProdigalTrainingFile
```

Option 2 - Adapt an external schema

Include all loci files in a directory and adapt the following command:
```
chewBBACA.py PrepExternalSchema -i ExternalSchemaFastaFiles -o OutputSchemaFolder --ptf ProdigalTrainingFile
```

### Perform allele calling

Determine the allelic profiles for genome assemblies:
```
chewBBACA.py AlleleCall -i InputAssemblies -g OutputSchemaFolder/SchemaName -o OutputFolderName
```

Use a subset of the loci in a schema:
```
chewBBACA.py AlleleCall -i InputAssemblies -g OutputSchemaFolder/SchemaName -o OutputFolderName --gl LociList.txt
```

**Important:**
- Genome assemblies and loci files from external schemas must be in FASTA format.
- We strongly advise users to provide a Prodigal training file and to keep using the same training file to ensure consistent results.
- Use the `--cpu` parameter to enable parallelization and considerably reduce execution time.
- The file passed to the `--gl` parameter must have one locus identifier per line (a locus identifier is the name of the FASTA file that contains the locus alleles).

## Detailed Usage

### i. Whole Genome Multilocus Sequence Typing (wgMLST) schema creation

Create a schema seed based on a set of FASTA files with genome assemblies or coding sequences.

Basic usage:

```
chewBBACA.py CreateSchema -i /path/to/InputAssemblies -o /path/to/OutputFolderName --n SchemaName --ptf /path/to/ProdigalTrainingFile --cpu 4
```

**Parameters:**

`-i`, `--input-files` Path to the directory that contains the input FASTA
     files. Alternatively, a single file with a list of
     paths to FASTA files, one per line.

`-o`, `--output-directory` Output directory where the process will store
     intermediate files and create the schema's
     directory.

`--n`, `--schema-name` (Optional) Name given to the folder that will store the schema
      files (default: schema_seed).

`--ptf`, `--training-file` (Optional) Path to the Prodigal training file. We strongly
        advise users to provide a Prodigal training file and to keep
        using the same training file to ensure consistent results
        (default: None).

`--bsr`, `--blast-score-ratio` (Optional) BLAST Score Ratio value. Sequences with alignments
        with a BSR value equal to or greater than this
        value will be considered as sequences from the same
        gene (default: 0.6).

`--l`, `--minimum-length` (Optional) Minimum sequence length value. Coding sequences
      shorter than this value are excluded (default: 201).

`--t`, `--translation-table` (Optional) Genetic code used to predict genes and to translate
      coding sequences (default: 11).

`--st`, `--size-threshold` (Optional) CDS size variation threshold. Added to the schema's
       config file and used to identify alleles with a
       length value that deviates from the locus length
       mode during the allele calling process (default: 0.2).

`--cpu`, `--cpu-cores` (Optional) Number of CPU cores that will be used to run the
        CreateSchema process (will be redefined to a lower
        value if it is equal to or exceeds the total number
        of available CPU cores)(default: 1).

`--b`, `--blast-path` (Optional) Path to the BLAST executables (default: assumes BLAST
	executables were added to PATH).

`--pm`, `--prodigal-mode` (Optional) Prodigal running mode (default: single).

`--CDS` (Optional) If provided, input is a single or several FASTA
        files with coding sequences (default: False).
	
`--no-cleanup` (Optional) If provided, intermediate files generated during process
	execution are not removed at the end (default: False).

**Outputs:**

```
OutputFolderName
├── SchemaName
│   ├── short
│   │   ├── GenomeID_proteinN_short.fasta
│   │   ├── ...
│   │   └── GenomeID_proteinN_short.fasta
│   ├── GenomeID_proteinN.fasta
│   ├── ...
│   ├── GenomeID_proteinN.fasta
│   └── Training_file.trn
├── invalid_alleles.txt
└── cds_info.tsv
```

One fasta file per distinct gene identified in the schema creation process in the `OutputFolderName/SchemaName` directory. The name attributed to each fasta file in the schema is based on the genome of origin of the first allele identified for that gene and on the order of gene prediction (e.g.: `GCA-000167715-protein12.fasta`, first allele for the gene was identified in an assembly with the prefix `GCA-000167715` and the gene was the 12th gene predicted by Prodigal in that assembly). The `OutputFolderName/SchemaName` directory also contains a directory named `short` that includes fasta files with the representative sequences for each locus. The training file passed to create the schema is also included in `OutputFolderName/SchemaName` and will be automatically detected during the allele calling process. A file with the locations of the identified genes in each genome passed to create the schema, `cds_info.tsv`, and a file with the list of alleles predicted by Prodigal that were excluded in the subsequent steps , `invalid_alleles.txt`, are included in `OutputFolderName`.

--------------

### ii. Allele call using a cg/wgMLST schema 

Perform allele calling to determine the allelic profiles of a set of samples in FASTA format. The
process identifies new alleles, assigns an integer identifier to those alleles and adds them to the
schema.

Basic usage:

```
chewBBACA.py AlleleCall -i /path/to/InputAssemblies -g /path/to/SchemaName -o /path/to/OutputFolderName --cpu 4
```

**Parameters:** 

`-i`, `--input-files` Path to the directory with the genome FASTA files or to a file
     with a list of paths to the FASTA files, one per line.

`-g`, `--schema-directory` Path to the schema directory with the genes FASTA files.  

`-o`, `--output-directory` Output directory where the allele calling results will be stored.

`--ptf`, `--training-file` (Optional) Path to the Prodigal training file. Default is to
        get training file from the schema's directory (default: searches for a training 
	file in the schema's directory).

`--gl`, `--genes-list` (Optional) Path to a file with the list of genes in the schema
       that the process should identify alleles for (default: False).

`--bsr`, `--blast-score-ratio` (Optional) BLAST Score Ratio value. Sequences with alignments
	with a BSR value equal to or greater than this value will be considered
	as sequences from the same gene (default: uses value defined in schema config).

`--l`, `--minimum-length` (Optional) Minimum sequence length accepted for a coding sequence
	to be included in the schema (default: uses value defined in schema config).

`--t`, `--translation-table` (Optional) Genetic code used to predict genes and to translate coding
	sequences. Must match the genetic code used to create the training file
	(default: uses value defined in schema config).

`--st`, `--size-threshold` (Optional) CDS size variation threshold. If set to a value
	of 0.2, alleles with size variation +-20 percent will be classified as
	ASM/ALM (default: uses value defined in schema config).

`--cpu`, `--cpu-cores` (Optional) Number of CPU cores/threads that will be used to
        run the CreateSchema process (will be redefined to
        a lower value if it is equal to or exceeds the
        total number of available CPU cores/threads)(default: 1).

`--b`, `--blast-path` (Optional) Path to the BLASTp executables. Use this option if chewBBACA cannot find
     BLASTp executables or if you want to use anoter BLAST installation that is not
     the one added to the PATH (default: assumes BLAST executables were added to PATH).

`--pm` (Optional) Prodigal running mode (default: single).

`--fc` (Optional) Continue the previous allele calling process if it was
       interrupted (default: False).

`--fr` (Optional) Force process reset even if there are temporary
       files from a previous process that was interrupted (default: False).

`--db`, `--store-profiles` (Optional) If the profiles in the output matrix should be
	stored in a local SQLite database. The SQLite database is stored in the "profiles_database"
	folder inside the schema's directory (default: False).

By default, the AlleleCall process uses the Prodigal training file included in the schema's directory and it is not necessary to pass a training file to the `--ptf` argument. If a text file with a list of gene identifiers, one per line, is passed to the `--gl` parameter, the process will only perform allele calling for the genes in the list.


**Outputs**:

```
OutputFolderName
└── results_datestamp
    ├── results_statistics.tsv
    ├── results_contigsInfo.tsv
    ├── results_alleles.tsv
    ├── RepeatedLoci.txt
    └── logging_info.txt
```

The `results_statistics.tsv` file contains the total number of exact matches (EXC), inferred new alleles (INF), loci not found (LNF), loci on contig tips (PLOT), non-informative paralogous hits (NIPH), alleles larger than locus length mode (ALM) and alleles smaller than locus length mode (ASM) classifications attributed for each genome.

The `results_contigsInfo.tsv` file contains the loci positions in the genomes analyzed. The first column contains the name of the genome files used in the allele calling and the other columns (with loci names in the headers) the locus position information or the classification attributed by chewBBACA if it was not an exact match or inferred allele.

The `results_alleles.tsv` file contains the allelic profiles determined for the input samples. The first column has the identifiers of the genome assemblies for which the allele call was performed. The remaining columns contain the allele call data for loci present in the schema, with the column headers being the locus identifiers.

The `RepeatedLoci.txt` file provides information about homologous loci detection. This output is useful to identify loci in the schema that are highly similar and loci that have a high number of CDS hits that are not exact matches or new inferred alleles.

Please visit the Wiki section about [Allele Calling](https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling) if you want to know more about the types of classifications attributed by chewBBACA and the output files.

--------------

### iii. Determine annotations for loci in the schema

The UniprotFinder process can be used to retrieve annotations for the loci in the schema through requests to [UniProt's SPARQL endpoint](http://sparql.uniprot.org/sparql) and through alignment against the reference proteomes for a set of taxa.

Basic usage:

```
chewBBACA.py UniprotFinder -i /path/to/SchemaName -o /path/to/OutputFolderName -t /path/to/cds_info.tsv --taxa "Species Name" --cpu 4
```

**Parameters:**

`-i`, `--input-files` Path to the schema's directory or to a file with a list of
     paths to loci FASTA files, one per line.

`-o`, `--output-directory` Output directory where the process will store
     intermediate files and save the final TSV file with the
     annotations.

`-t`, `--protein-table` (Optional) Path to the "cds_info.tsv" file created by the
     CreateSchema process (default: None).

`--bsr` (Optional) BLAST Score Ratio value. This value is only used when a
        taxon/taxa is provided and local sequences are aligned
        against reference proteomes (default: 0.6).

`--cpu`, `--cpu-cores` (Optional) Number of CPU cores used to run the process (default: 1).

`--taxa` (Optional) List of scientific names for a set of taxa. The process
         will search for and download reference proteomes with
         terms that match any of the provided taxa (default: None).

`--pm` (Optional) Maximum number of proteome matches to report (default: 1).

**Outputs:**

The `/path/to/OutputFolderName` directory contains a TSV file, `schema_annotations.tsv`, with the information found for each locus. The process will always search for annotations through UniProt's SPARQL endpoint, reporting the product name and UniProt URL for local loci with an exact match in UniProt's database. If the `cds_info.tsv` file is passed to the `-t` parameter, the output file will also include the information in that file. The `--taxa` parameter receives a set of taxa names and searches for reference proteomes that match the provided terms. The reference proteomes are downloaded and the process aligns schema representative sequences against the reference proteomes to include additional information in the `schema_annotations.tsv` file based on matches against the sequences in the reference proteomes.

--------------

### iv. Evaluate wgMLST call quality per genome

Basic usage:

```
chewBBACA.py TestGenomeQuality -i /path/to/AlleleCall/results/results_alleles.tsv -n 12 -t 200 -s 5 -o /path/to/OutputFolderName
```

`-i` Path to file with a matrix of allelic profiles (i.e. results_alleles.tsv).

`-n` Maximum number of iterations. Each iteration removes a set of genomes over the
     defined threshold (-t) and recalculates loci presence percentages.

`-t` Maximum threshold. This threshold represents the maximum number of missing loci
     allowed, for each genome independently, before removing the genome.

`-s` Step to add to each threshold (suggested 5).

`-o` Path to the output directory that will store output files.

The output is a HTML file with a plot with all thresholds and a `removedGenomes.txt` file with
information about which genomes were removed per threshold when it reaches a stable point
(no more genomes are removed).

Example of an output can be seen [here](http://im.fm.ul.pt/chewBBACA/GenomeQual/GenomeQualityPlot_all_genomes.html).
The example uses an original set of 714 genomes and a scheme consisting of 3266 loci with `-n 12`, `-t 300` and `-s 5`
passed to arguments.

--------------

### v. Defining the cgMLST schema

Determine the set of loci that constitute the core genome based on a threshold.

Basic usage:

```
chewBBACA.py ExtractCgMLST -i /path/to/AlleleCall/results/results_alleles.tsv -o /path/to/OutputFolderName
```
	
`-i`, `--input-file` Path to input file containing a matrix with allelic profiles.

`-o`, `--output-directory` Path to the directory where the process will store output files.

`--t`, `--threshold` (Optional) Genes that constitute the core genome must be in a
      proportion of genomes that is at least equal to this value.
      (e.g 0.95 to get a matrix with the loci that are present in at 
      least 95% of the genomes) (default: 1).

`--r`, `--genes2remove` (Optional) Path to file with a list of genes/columns to remove 
      from the matrix (one gene identifier per line, e.g. the list of
      genes listed in the RepeatedLoci.txt file created by the AlleleCall
      process) (default: False).

`--g`, `--genomes2remove` (Optional) Path to file with a list of genomes/rows to remove from the
      matrix (one genome identifier per line, e.g. list of genomes to be 
      removed based on the results from the TestGenomeQuality process) (default: False).

**Note:** The matrix with allelic profiles created by the ExtractCgMLST
          process can be imported into [**PHYLOViZ**](https://online.phyloviz.net/index)
	  to visualize and explore typing results.

--------------

### vi. Evaluate your schema

Evaluate the number of alelles and allele size variation for the loci in a schema or for a set of
selected loci. Provide information about problematic alleles per locus and individual pages for each
locus with a plot with allele size, a Neighbor Joining tree based on a multiple sequence alignment
(MSA) and a visualization of the MSA.
 
See an example [here](https://saureus-report.herokuapp.com/)

Basic usage:

```
chewBBACA.py SchemaEvaluator -i /path/to/SchemaName -o /path/to/OutputFolderName --cpu 4
```
	
`-i`, `--input-files` Path to the schema's directory or path to a file containing the
     paths to the FASTA files of the loci that will be evaluated, one 
     per line.

`-o`, `--output` Path to the output directory where the report HTML
     files will be generated.

`-a`, `--annotations` (Optional) Path to the TSV table created by the UniprotFinder
     process.
     
`--ta`, `--translation-table` (Optional) Genetic code used to translate coding sequences
       (default: [11](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1)).

`--th`, `--threshold` (Optional) Allele size variation threshold. If an allele has a
       size within the interval of the locus mode -/+ the
       threshold, it will be considered a conserved allele (default: 0.05).

`--ml`, `--minimum-length` (Optional) Minimum sequence length accepted for a coding
       sequence to be included in the schema.

`--cpu` (Optional) Number of CPU cores to use to run the process (default: 1).

`--light` (Optional) Skips the indepth analysis of the individual schema
          loci, including MAFFT (default: False).

`--no-cleanup` Stops the removal of intermediate files created during the
	report generation.

Please consult the [SchemaEvaluator's wiki page](https://github.com/B-UMMI/chewBBACA/wiki/4.-Schema-Evaluation) for more information.

--------------

### vii. Adapt an external schema

The PrepExternalSchema process enables the adaptation of external schemas so that it is possible to use those schemas with chewBBACA. An external schema may be a set of sequences from any number of genes that have been selected for a particular study or it may be a schema that has already been defined and is available for download from some well known databases, such as [Ridom](https://www.cgmlst.org/ncs), [BIGSdb](https://pubmlst.org/) and [Enterobase](http://enterobase.warwick.ac.uk/).

Basic usage:

```
chewBBACA.py PrepExternalSchema -i /path/to/ExternalSchemaFastaFiles -o /path/to/OutputFolderName --ptf /path/to/ProdigalTrainingFile --cpu 4
```

`-i`, `--input-files` Path to the folder containing the FASTA files, one FASTA
	file per gene/locus (alternatively, a file with a list of paths can be
	given).

`-o`, `--output-directory` The directory where the output files will be saved
	(will create the directory if it does not exist).

`--ptf`, `--training-file` Path to the Prodigal training file that will be
	included in the adapted schema (default: None).

`--bsr`, `--blast-score-ratio` The BLAST Score Ratio value that will be used to
	adapt the external schema (default: 0.6).

`--l`, `--minimum-length` Minimum sequence length accepted. Sequences with a
	length value smaller than the value passed to this argument will be
	discarded (default: 0).

`--t`, `--translation-table` Genetic code to use for CDS translation.
	Must match the genetic code used to create the training file (default: 11).

`--st`, `--size-threshold` CDS size variation threshold. At the default value
	of 0.2, alleles with size variation +-20 percent when compared to
	the representative will not be included in the final schema (default: 0.2).

`--cpu`, `--cpu-cores` The number of CPU cores to use (default: 1).

--------------

## FAQ

### Q: What strains should I use for [Schema Creation](#i-whole-genome-multilocus-sequence-typing-wgmlst-schema-creation)?
A: The set of genome assemblies used for schema creation should be carefully selected to avoid the inclusion of spurious loci in the schema seed that result from low quality assemblies (e.g.: genome assemblies resulting from low quality sequencing data, highly fragmented genome assemblies, genome assemblies with many frameshifted proteins, genome length too large or too small). A set of high quality genome assemblies, ideally complete genomes, that capture the diversity of the species or lineage of interest should result in a good schema seed.

### Q: [Allele calling](#ii-allele-call-using-a-cgwgmlst-schema) just crashed, do I need to start over?  
A: chewBBACA should allow you to continue where you stopped, just re-run the same command and you should be prompted to continue the allele call or use the flag `--fc`. If the process keeps crashing with the same set of inputs, it is very likely that one or more of the inputs is misformatted or has a format that is incompatible with chewBBACA. Please ensure that your input files are genome assemblies or coding sequences in FASTA format. Consider opening an [issue](https://github.com/B-UMMI/chewBBACA/issues) to report the problem, we will do our best to help solve the issue and user feedback is very important for the continued improvement of chewBBACA.

### Q: I ran all the steps and my cgMLST loci size is smaller than traditional MLST, does this even work?  
A: In order to have a robust definition of a cgMLST schema for a given bacterial species, a set of representative strains of the diversity of a given species should be selected. Furthermore, since cgMLST schema definition is based on pre-defined thresholds, only when a sufficient number of strains have been analyzed can the cgMLST schema be considered stable. This number will always depend on the population structure and diversity of the species in question. cgMLST schemas are defined as the set of loci that are present in all strains under analysis, but defining a smaller loci presence threshold, such as 95%, might be necessary to include very frequent genes or ubiquitous genes that are not present in some strains due to sequencing/assembly limitations.
The quality of the genome assemblies is also an important factor that can impact profoundly the MLST schema definition (e.g.: genome assemblies with a high number of missing loci, contaminated genome assemblies, misclassified genome assemblies). The [TestGenomeQuality](#iv-evaluate-wgmlst-call-quality-per-genome) process can be used to identify genome assemblies that are responsible for a considerable loss of loci. You can pass a list of genomes to remove to the [ExtractCgMLST](#v-defining-the-cgmlst-schema) process to exclude those genomes from the analysis that determines the set of loci that constitute the core-genome. Identifying and removing low quality genome assemblies can significantly improve the determination of the core-genome.

### Q: Can I use a schema from an external source?  
A: Yes. The [PrepExternalSchema](#vii-adapt-an-external-schema) process enables the adaptation of external schemas so that it is possible to use those schemas with chewBBACA. An external schema may be a set of sequences from any number of genes that have been selected for a particular study or it may be a schema that has already been defined and is available for download from some well known databases, such as [Ridom cgMLST](http://www.cgmlst.org/ncs), [BIGSdb](https://pubmlst.org/) and [Enterobase](http://enterobase.warwick.ac.uk/).

### Q: Which species already have a training file?  
A: At the moment:
 - *Acinetobacter baumannii*
 - *Campylobacter jejuni*
 - *Enterococcus faecium*
 - *Escherichia coli*
 - *Haemophilus influenzae*
 - *Legionella pneumophila*
 - *Listeria monocytogenes*
 - *Salmonella enterica enteritidis*
 - *Staphylococcus aureus*
 - *Staphylococcus haemolyticus*
 - *Streptococcus agalactiae*
 - *Streptococcus canis*
 - *Streptococcus dysgalactiae*
 - *Streptococcus equi*
 - *Streptococcus pneumoniae*
 - *Streptococcus pyogenes*
 - *Yersinia enterocolitica*

get them [here](https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files).

### Q: My favorite species has no training file. What can I do?
A: You can propose a new one to be added to the repository or create your own training 
files. To create a training file make sure you have Prodigal installed and run the following command:

```
prodigal -i myGoldStandardGenome.fna -t myTrainedFile.trn -p single
```

## Citation

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
