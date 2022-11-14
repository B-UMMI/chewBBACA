
[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/chewBBACA/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/chewbbaca)
[![Conda](https://img.shields.io/conda/dn/bioconda/chewbbaca?color=green)](https://anaconda.org/bioconda/chewbbaca)
[![chewBBACA](https://github.com/B-UMMI/chewBBACA/workflows/chewbbaca/badge.svg)](https://github.com/B-UMMI/chewBBACA/actions?query=workflow%3Achewbbaca)
[![Documentation Status](https://readthedocs.org/projects/chewbbaca/badge/?version=latest)](https://chewbbaca.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/github/license/B-UMMI/chewBBACA)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:10.1099/mgen.0.000166](https://img.shields.io/badge/DOI-10.1099%2Fmgen.0.000166-blue)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166)

# chewBBACA

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be thought of as "Comprehensive and  Highly Efficient Workflow" 
but at this point still it needs a bit of work to make that claim so we just add "chew" to add extra coolness to the software name. BSR stands for 
BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2) 

chewBBACA is a comprehensive pipeline including a set of functions for the creation and validation of whole genome and core genome MultiLocus Sequence 
Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on Blast Score Ratio that can be run in multiprocessor 
settings and a set of functions to visualize and validate allele variation in the loci. chewBBACA performs the schema creation and allele calls on complete or draft genomes resulting from de novo assemblers.

## Citation

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

- chewBBACA only works with **python 3** (automatic testing for Python 3.8 and Python 3.9 with GitHub Actions).
- We strongly recommend that users install and use **BLAST 2.9.0+** with chewBBACA, as chewBBACA's processes have been extensively tested with that version of BLAST.
- chewBBACA includes Prodigal training files for some species. You can consult the list of Prodigal training files that are readily available [here](https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files). We strongly recommend using the same Prodigal training file for schema creation and allele calling to ensure consistent results.
- chewBBACA defines an allele as a complete Coding DNA Sequence, with start and stop codons according 
 to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) identified using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). It will automatically exclude any allele for which the DNA sequence does not contain start or stop codons and for which the length is not multiple of three. Alleles that contain ambiguous bases are also excluded.
- Make sure that your FASTA files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems, you should do a quick conversion using for example [dos2unix](https://waterlan.home.xs4all.nl/dos2unix.html).

## Useful links

 - Check the [wiki pages](https://github.com/B-UMMI/chewBBACA/wiki) for a much more thorough chewBBACA walkthrough.
 - An extensive [tutorial repository](https://github.com/B-UMMI/chewBBACA_tutorial) is available as example on how to run an analysis pipeline using chewBBACA.
 - Look through the list of schemas in [Chewie-NS](https://chewbbaca.online/). Chewie-NS allows chewBBACA users to download and update cg/wgMLST schemas, allowing the easy sharing of results, while ensuring the reproducibility and consistency of these steps.
 - Use [BBACA gitter](https://gitter.im/BBACA/Lobby) if you have any pressing question. Chat can be faster and better than email for troubleshooting purposes.
 - A ready to use [docker image](https://hub.docker.com/r/ummidock/chewbbaca) automatically built from the latest version of chewBBACA in Ubuntu 16.04.
 - chewBBACA is available as a [Galaxy module](https://toolshed.g2.bx.psu.edu/repository?repository_id=88fd7663075eeae9&changeset_revision=093352878303) many thanks to Stefano Morabito and [Arnold Knijn](https://github.com/aknijn) for EURL VTEC in ISS, Rome!
 - Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to know more about the latest changes!

## Installation

Install the latest released version using [conda](https://anaconda.org/bioconda/chewbbaca):

```
conda create -c bioconda -c conda-forge -n chewie chewbbaca=2.8.5 blast=2.9
```

Install using [pip](https://pypi.org/project/chewBBACA/):

```
pip3 install chewbbaca
```

chewBBACA has the following dependencies:

Python dependencies (defined in the [requirements](https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/requirements.txt) file, should be automatically installed when using conda or pip):
* numpy>=1.14.0
* scipy>=0.13.3
* biopython>=1.70
* plotly>=1.12.9
* SPARQLWrapper>=1.8.0
* requests>=2.2.1
* pandas>=0.22.0

Main dependencies:
* [BLAST 2.9.0+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/)
* [Prodigal 2.6.0](https://github.com/hyattpd/prodigal/releases/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/) (for schema evaluation only)

Installation through conda should take care of all dependencies. If you install through pip you will need to ensure that you have BLAST, Prodigal and MAFFT installed and added to the PATH.

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
