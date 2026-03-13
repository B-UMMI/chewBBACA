
[![PyPI](https://img.shields.io/badge/Install%20with-PyPI-blue)](https://pypi.org/project/chewBBACA/#description)
[![Bioconda](https://img.shields.io/badge/Install%20with-bioconda-green)](https://anaconda.org/bioconda/chewbbaca)
[![Conda](https://img.shields.io/conda/dn/bioconda/chewbbaca?color=green)](https://anaconda.org/bioconda/chewbbaca)
[![chewBBACA](https://github.com/B-UMMI/chewBBACA/workflows/chewbbaca/badge.svg)](https://github.com/B-UMMI/chewBBACA/actions?query=workflow%3Achewbbaca)
[![Documentation Status](https://readthedocs.org/projects/chewbbaca/badge/?version=latest)](https://chewbbaca.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/github/license/B-UMMI/chewBBACA)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI:10.1099/mgen.0.000166](https://img.shields.io/badge/DOI-10.1099%2Fmgen.0.000166-blue)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166)

# chewBBACA

**chewBBACA** is a software suite for the creation and evaluation of core genome and whole genome MultiLocus Sequence 
Typing (cg/wgMLST) schemas and results. The "BBACA" stands for "BSR-Based Allele Calling Algorithm". BSR stands for 
BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2). The "chew" part adds extra coolness to the name and could be thought of as "Comprehensive and Highly Efficient Workflow". chewBBACA allows to define the target loci in a schema based on multiple genomes (e.g. define target loci based on the distinct loci identified in a dataset of high-quality genomes for a species or lineage of interest) and performs allele calling to determine the allelic profiles of bacterial strains, easily scaling to thousands of genomes with modest computational resources. chewBBACA includes functionalities to annotate the schema loci, compute the set of loci that constitute the core genome for a given dataset, and generate interactive reports for schema and allele calling results evaluation to enable an intuitive analysis of the results in surveillance and outbreak detection settings or population studies. Pre-defined cg/wgMLST schemas can be downloaded from [Chewie-NS ](https://chewbbaca.online/) or adapted from other cg/wgMLST platforms.

### Check the [documentation](https://chewbbaca.readthedocs.io/en/latest/index.html) for implementation details and guidance on using chewBBACA.

## News

## 3.6.0 - 2026-03-

This version adds multiple modules that perform operations that were previously integrated into other modules and could not be used separately. The new modules are the following:

- **PredictGenes**: this module predicts genes from genome assemblies in FASTA format. In past versions, gene prediction was integrated into the CreateSchema and AlleleCall modules and could not be used separately to perform gene prediction and keep the FASTA files with the predicted CDSs. With this new module, users can predict genes for genomes of interest and provide the FASTA files for subsequent analyses with the CreateSchema and AlleleCall modules (through the `--cds-input` option) or to perform other analyses with external software. Performing gene prediciton with the PredictGenes module and storing the FASTA files also eliminates the need to perform gene prediction each time that dataset is analysed, reducing the runtime of the analyses.

- **ComputeDistances**: this module computes pairwise distances based on the allelic profiles determined by the AlleleCall module. Previously, pairwise ditance computation was integrated into the AlleleCallEvalautor module. By creating a separte module for pairwise distance computation, users do not have to run the AlleleCallEvaluator module if all they want is a matrix with the pariwise distances. Creating a separate module also allowed to implement more functionalities, such as selecting between different types of distances (hamming, jaccard and loci (number of loci not shared)) and output file formats (upper- and lower-triangular matrices, symmetrice and table), and to provide the option to compute pairwise similarity instead of distance.

- **HashProfiles**: the `--hash-profiles` option integrated into the AlleleCall module was isolated to implement the HashProfiles module. The HashProfiles module allows users to hash the allelic profiles by providing the `results_alleles.tsv` file created by the AlleleCall module. Previously, users had to run the AlleleCall module to get the hashed profiles, but with the HashProfiles module users can hash the allelic profiles from the latest or past runs.

- **SubsetResults**: the SubsetResults module enables users to subset results based on loci and sample lists. The SubsetResults module superseeds the RemoveGenes module while also providing more options, such as excluding the list of loci and/or samples provided to keep only remaining loci and/or samples.

- **MergeResults**: the MergeResults allows to merge allele calling results from different runs, superseeding the JoinProfiles module. It can merge results that were generated based on the same schema.

Additional changes:

- Added option to determine the set of loci that constitute the accessory genome to the ExtractCgMLST module.

Check our [Changelog](https://github.com/B-UMMI/chewBBACA/blob/master/CHANGELOG.md) to learn about the latest changes.

## Citation

When using chewBBACA, please use the following citation:

> Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)
