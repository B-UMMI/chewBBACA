
[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/chewBBACA/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/chewbbaca)
[![Conda](https://img.shields.io/conda/dn/bioconda/chewbbaca?color=green)](https://anaconda.org/bioconda/chewbbaca)
[![chewBBACA](https://github.com/B-UMMI/chewBBACA/workflows/chewbbaca/badge.svg)](https://github.com/B-UMMI/chewBBACA/actions?query=workflow%3Achewbbaca)
[![Documentation Status](https://readthedocs.org/projects/chewbbaca/badge/?version=latest)](https://chewbbaca.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/github/license/B-UMMI/chewBBACA)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:10.1099/mgen.0.000166](https://img.shields.io/badge/DOI-10.1099%2Fmgen.0.000166-blue)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166)

# chewBBACA

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be thought of as "Comprehensive and  Highly Efficient Workflow" 
but at this point it still needs a bit of work to make that claim, so we just add "chew" to add extra coolness to the software name. BSR stands for 
BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2) 

chewBBACA is a comprehensive pipeline including a set of functions for the creation and validation of whole genome and core genome MultiLocus Sequence 
Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on Blast Score Ratio that can be run in multiprocessor 
settings and a set of functions to visualize and validate allele variation in the loci. chewBBACA performs the schema creation and allele calls on complete or draft genomes.

### Check the [documentation](https://chewbbaca.readthedocs.io/en/latest/index.html) for implementation details and guidance on using chewBBACA.

## News

## 3.3.3 - 2024-02-23

- Fixed warning related with BLASTp `--seqidlist` parameter. For BLAST>=2.9, the TXT file with the sequence IDs is converted to binary format with `blastdb_aliastool`.

- The `Bio.Application` modules are deprecated and might be removed from future Biopython versions. Modified the function that calls MAFFT so that it uses the subprocess module instead of `Bio.Align.Applications.MafftCommandline`. Changed the Biopython version requirement to >=1.79.

- Added a `pyproject.toml` configuration file and simplified the instructions in `setup.py`. The use of `setup.py` as a command line tool is deprecated and the `pyproject.toml` configuration file allows to install and build packages through the recommended method.

- Updated the Dockerfile to install chewBBACA with `python3 -m pip install .` instead of the deprecated `python setup.py install` command.

- Removed FASTA header integer conversion before running BLASTp. This was done to avoid a warning from BLAST related to sequence header length exceeding 50 characters.

- The seqids and coordinates of the CDSs closest to contig tips are stored in a dictionary during gene prediction to simplify LOTSC and PLOT5/3 determination (in many cases this reduces runtime by ~20%).

- Limited the number of values stored in memory while creating the `results_contigsInfo.tsv` and `results_alleles.tsv` output files to reduce memory usage.

- Adding data to the FASTA and TSV files for the missing classes per locus instead of storing the complete per input data to reduce memory usage.

- The data for novel alleles is saved to files to reduce memory usage.

- Fixed the in-frame stop codon count values displayed in the reports created by the SchemaEvaluator module.

- The `UniprotFinder` module now exits cleanly if the output directory already exists.

- Improved info printed to the stdout by the CreateSchema and AlleleCall modules, added comments, and changed variable names to better match data being stored.

Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to learn about the latest changes.

## Citation

When using chewBBACA, please use the following citation:

> Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
