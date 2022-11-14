
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

Check the [documentations](https://chewbbaca.readthedocs.io/en/latest/index.html) for implementation details and guidance on using chewBBACA.

## Citation

When using chewBBACA please use the following citation:

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)

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
