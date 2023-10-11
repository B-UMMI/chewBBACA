
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

## 3.3.0 - 2023-10-11

- Added the AlleleCallEvaluator module. This module generates an interactive HTML report for the allele calling
results. The report provides summary statistics to evaluate results per sample and per locus (with the possibility
to provide a TSV file with loci annotations to include on a table). The report includes components to display a
heatmap representing the loci presence-absence matrix, a heatmap representing the distance matrix based on allelic
differences and a Neighbor-Joining tree based on the MSA of the core genome loci.

- Added [pyrodigal](https://github.com/althonos/pyrodigal) for gene prediction. This simplified the processing of the gene prediction results and reduced runtime.

- Fixed an issue where the AlleleCall module would try to create results files for excluded inputs.

- Fixed exception capturing during multiprocessing when using Python>=3.11.

- Fixed PLOT5/3 identification when coding sequences are in the reverse strand.

- Fixed computation of the representative self-scores when performing allele calling for a subset of the loci in a schema (would only compute the self-scores for the subset of loci if the 'self_scores' file had still not been created).

- Fixed issue related to the classification of single EXC/INF and single/multiple ASM/ALM (would classify some inputs as NIPH instead of EXC/INF).

- Fixed issue related to protein exact match classification when multiple pre-computed PROTEINtable files include the same protein hash.

- Changed the `-i`, `--input-files` parameter in the PrepExternalSchema and UniprotFinder modules to `-g`, `--schema-directory` and added the `--gl`, `--genes-list` parameter to enable adapting or annotating a subset of the loci in the schema.

Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to learn about the latest changes.

## Citation

When using chewBBACA, please use the following citation:

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
