Overview
========

About chewBBACA
---------------

**chewBBACA** is a software suite for the creation and evaluation of core genome and whole genome MultiLocus Sequence 
Typing (cg/wgMLST) schemas and results. The "BBACA" stands for "BSR-Based Allele Calling Algorithm". BSR stands for 
BLAST Score Ratio as proposed by `Rasko DA et al. <http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2>`_. 
The "chew" part adds extra coolness to the name and could be thought of as "Comprehensive and Highly Efficient Workflow". 
chewBBACA allows to define the target loci in a schema based on multiple genomes (e.g. define target loci based on the distinct 
loci identified in a dataset of high-quality genomes for a species or lineage of interest) and performs allele calling to determine 
the allelic profiles of bacterial strains, easily scaling to thousands of genomes with modest computational resources. chewBBACA 
includes functionalities to annotate the schema loci, compute the set of loci that constitute the core genome based on the results
for a given dataset, and generate interactive reports for schema and allele calling results evaluation to enable an intuitive
analysis of the results in surveillance and outbreak detection settings or population studies. Pre-defined cg/wgMLST schemas can
be downloaded from  `Chewie-NS <https://chewbbaca.online/>`_ or adapted from other cg/wgMLST platforms.

The main processes available in chewBBACA are represented in the following workflow:

.. image:: /_static/images/Overview.png
   :width: 1000px
   :align: center

- Steps labelled 1: Schema creation from genome assemblies or coding sequences in FASTA format.
- Steps labelled 2: Adaptation of external schemas for usage with chewBBACA.
- Steps labelled 3: Upload, download, and synchronize schemas from `Chewie-NS <https://chewbbaca.online/>`_.
- Steps labelled 4: Perform allele calling to determine the allelic profiles of strains of interest.
- Steps labelled 5: Determine the set of loci that constitute the core genome based on allele calling results.
- Steps labelled 6: Annotate schema loci based on UniProt data.
- Steps labelled 7: Evaluate schemas and explore loci diversity through an interactive report.
- Steps labelled 8: Evaluate allele calling results through an interactive report.

Citation
--------

If you use chewBBACA, please cite the following publication:

  Mamede R, Vila-Cerqueira P, Carriço JA, Ramirez M. 2026. chewBBACA 3: lowering the barrier for scalable and detailed whole- and core-genome multilocus sequence typing. Genome Med 18:51. `https://doi.org/10.1186/s13073-026-01625-x <https://doi.org/10.1186/s13073-026-01625-x>`_

The supplementary material for the publication is available on `Zenodo <https://doi.org/10.5281/zenodo.14637858>`_.

Other relevant citations:

- chewBBACA's first publication (if using a version < 2.1.0):


  Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. `https://doi.org/10.1099/mgen.0.000166 <https://doi.org/10.1099/mgen.0.000166>`_


- Chewie-NS (if using any of the schemas deposited on Chewie-NS):


  Mamede R, Vila-Cerqueira P, Silva M, Carriço JA, Ramirez M. 2021. Chewie Nomenclature Server (chewie-NS): a deployable nomenclature server for easy sharing of core and whole genome MLST schemas. Nucleic Acids Res 49:D660–D666. `https://doi.org/10.1093/nar/gkaa889 <https://doi.org/10.1093/nar/gkaa889>`_

Please consider also citing the awesome software dependencies used by chewBBACA:

- General dependencies:
  
  - `BLAST <https://blast.ncbi.nlm.nih.gov/doc/blast-help/references.html>`_.
  - `Prodigal <https://link.springer.com/article/10.1186/1471-2105-11-119>`_ (if using chewBBACA<3.3.0).
  - `MAFFT <https://academic.oup.com/mbe/article/30/4/772/1073398>`_.
  - `FastTree <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490>`_.

- Python dependencies:
  
  - `NumPy <https://numpy.org/citing-numpy/>`_.
  - `SciPy <https://scipy.org/citing-scipy/>`_.
  - `Biopython <https://academic.oup.com/bioinformatics/article/25/11/1422/330687>`_.
  - `Plotly <https://github.com/plotly/plotly.py/blob/main/CITATION.cff>`_.
  - `SPARQLWrapper <https://github.com/RDFLib/sparqlwrapper>`_.
  - `Requests <https://requests.readthedocs.io/en/latest/>`_.
  - `pandas <https://pandas.pydata.org/about/citing.html>`_.
  - `Pyrodigal <https://joss.theoj.org/papers/10.21105/joss.04296>`_ (if using chewBBACA>=3.3.0).

- JavaScript dependencies:

  - `Material UI <https://www.npmjs.com/package/@mui/material>`_.
  - `MUI-Datatables <https://www.npmjs.com/package/mui-datatables>`_.
  - `MSA Viewer <https://www.npmjs.com/package/@jlab-contrib/msa>`_.
  - `Phylocanvas.gl <https://www.npmjs.com/package/@phylocanvas/phylocanvas.gl>`_.
  - `react-plotly.js <https://www.npmjs.com/package/react-plotly.js>`_.
  - `react-scroll <https://www.npmjs.com/package/react-scroll>`_.
  - `Monaco Editor for React <https://www.npmjs.com/package/@monaco-editor/react>`_.
  - `react-markdown <https://www.npmjs.com/package/react-markdown>`_.
  - `remark-gfm <https://www.npmjs.com/package/remark-gfm>`_.

Licensing
---------

This project is licensed under the `GPLv3 license 
<https://github.com/B-UMMI/Nomenclature_Server_docker_compose/blob/master/LICENSE>`_.
The source code of chewBBACA is available on `GitHub <https://github.com/B-UMMI/chewBBACA>`_.

Funding
-------

- **INNUENDO** project co-funded by the European Food Safety Authority (EFSA), grant agreement
  GP/EFSA/AFSCO/2015/01/CT2 ("New approaches in identifying and characterizing microbial and
  chemical hazards"). The conclusions, findings, and opinions expressed in this review paper
  reflect only the view of the authors and not the official position of the European Food Safety
  Authority (EFSA).
- **ONEIDA** project (LISBOA-01-0145-FEDER-016417) co-funded by FEEI - "Fundos Europeus Estruturais
  e de Investimento" from "Programa Operacional Regional Lisboa 2020" FCT - "Fundação para a
  Ciência e a Tecnologia".
- **BacGenTrack** (TUBITAK/0004/2014) [FCT/ Scientific and Technological Research Council of Turkey
  (Türkiye Bilimsel ve Teknolojik Araşrrma Kurumu, TÜBİTAK)]

.. image:: /_static/images/chewie_funding.png
   :width: 700px
   :align: center
