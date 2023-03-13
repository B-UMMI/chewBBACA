
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
settings and a set of functions to visualize and validate allele variation in the loci. chewBBACA performs the schema creation and allele calls on complete or draft genomes.

### Check the [documentation](https://chewbbaca.readthedocs.io/en/latest/index.html) for implementation details and guidance on using chewBBACA.

## News

## 3.1.2 - 2023-03-13

- Fixed issue related with sequence header format when FASTA files with coding sequences were provided to the CreateSchema and AlleleCall modules through the `--cds-input` parameter.

Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to know more about the latest changes.

## Citation

When using chewBBACA please use the following citation:

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
