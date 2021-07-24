
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

## IMPORTANT

- chewBBACA only works with **python 3** (automatic testing for Python 3.7 and Python 3.8 with GitHub Actions).
- We strongly recommend that users install and use BLAST 2.9.0+ with chewBBACA, as chewBBACA's processes have been extensively tested with that version of BLAST.
- chewBBACA includes Prodigal training files for some species. You can consult the list of Prodigal training files that are readily available [here](https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files). We strongly recommend using the same Prodigal training file for schema creation and allele calling to ensure consistent results.
- chewBBACA defines an allele as a complete Coding DNA Sequence, with start and stop codons according 
 to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) identified using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). It will automatically exclude any allele for which the DNA sequence does not contain start or stop codons and for which the length is not multiple of three. 
- Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems, you should do a a quick conversion using for example [dos2unix](https://waterlan.home.xs4all.nl/dos2unix.html).

## Contents

- [Useful links](#useful-links).
- [Installation](#installation).
- [Detailed Usage](#detailed-usage).
  1. [Whole Genome Multilocus Sequence Typing (wgMLST) schema creation]().
  2. [Allele call using a cg/wgMLST schema]().
  3. [Determine annotations for loci in the schema]().
  4. [Evaluate wgMLST call quality per genome]().
  5. [Defining the cgMLST schema]().
  6. [Evaluate your schema]().
  7. [Adapt an external schema](#vii-adapt-an-external-schema).
- [FAQ](#faq).
- [Citation](#citation).

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
* [mafft](https://mafft.cbrc.jp/alignment/software/)

Installation through conda should take care of all dependencies. If you install through pip you will need to ensure that you have BLAST and Prodigal installed and added to the PATH.

## Detailed Usage

### i. Whole Genome Multilocus Sequence Typing (wgMLST) schema creation

Create a schema seed based on a set of FASTA files with genome assemblies or coding sequences.

Basic usage:

```
chewBBACA.py CreateSchema -i /path/to/InputAssemblies -o /path/to/OutputFolderName --n SchemaName --ptf /path/to/ProdigalTrainingFile --cpu 4
```

**Parameters**

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

One fasta file per distinct gene identified in the schema creation process in the `/path/to/OutputFolderName/SchemaName` directory. The name attributed to each fasta file in the schema is based on the genome of origin of the first allele of that gene and on the order of gene prediction (e.g.: `GCA-000167715-protein12.fasta`, first allele for the gene was identified in an assembly with the prefix `GCA-000167715` and the gene was the 12th gene predicted by Prodigal in that assembly). The CreateSchema process also creates a file, "cds_info.tsv", in `/path/to/OutputFolderName/` with the locations of the identified genes in each genome passed to create the schema.


### ii. Allele call using a cg/wgMLST schema 

Perform allele calling to determine the allelic profiles of a set of samples in FASTA format. The
process identifies new alleles, assigns an integer identifier to those alleles and adds them to the
schema.

Basic usage:

```
chewBBACA.py AlleleCall -i /path/to/InputAssemblies -g /path/to/SchemaName -o /path/to/OutputFolderName --cpu 4
```

**Parameters** 

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

By default, the AlleleCall process uses the Prodigal training file included in the schema's directory and it is not necessary to pass a training file to the `--ptf` argument. If a text file with a list of gene identifiers is passed to the `--gl` parameter, the process will only perform allele calling for the genes in the list.


**Outputs**:

The following files will be saved into the `/path/to/OutputFolderName/results_datestamp` directory:

```
results_statistics.tsv
results_contigsInfo.tsv
results_alleles.tsv
logging_info.txt
RepeatedLoci.txt
```

### iii. Determine annotations for loci in the schema

The UniprotFinder process can be used to retrieve annotations for the loci in the schema through requests to [UniProt's SPARQL endpoint](http://sparql.uniprot.org/sparql) and through alignment against the reference proteomes for a set of taxa.

Basic usage:

```
chewBBACA.py UniprotFinder -i /path/to/SchemaName -o /path/to/OutputFolderName -t /path/to/cds_info.tsv --taxa "Species Name" --cpu 4
```

**Parameters**

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

...Add info about performing allele calling only for the loci in the cgMLST...

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


## FAQ

### Q: Step 2 is taking hours, will it ever end?  
A: Depending on the variability of the strains used to create the schema and the number 
of CPUs you have selected, the computing time used will vary. The more variable the strains, the more BLAST 
comparisons will be made, meaning more time will be needed for finishing the analysis.

### Q: Step 3 just crashed at 99% after 2 days running, do I need to start over?  
A: chewBBACA should allow you to continue where you stopped, just re-run the same command and you should be prompted to continue the allele call or use the flag `--fc`.

### Q: I ran all the steps and my cgMLST loci size is smaller than traditional MLST, does this even work?  
A: You probably forgot to eliminate from the analysis genomes responsible for a considerable loss of loci. 
Try to run again step 4, remove some of those genomes and check if the cgMLST loci number rises.

### Q: Can I use a schema from an external source?
A: Yes. Be sure to have a single fasta for each locus and use the `PrepExternalSchema` process.

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
