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
- Steps labelled 4: Perform allele calling to determine allelic profiles.
- Steps labelled 5: Compute the core genome based on allele calling results.
- Steps labelled 6: Annotate schema loci based on UniProt data.
- Steps labelled 7: Evaluate schemas and explore loci diversity through an interactive report.
- Steps labelled 8: Evaluate allele calling results through an interactive report.

Citation
--------

**chewBBACA 1** was published in Microbial Genomics under the title:
**chewBBACA: A complete suite for gene-by-gene schema creation and strain identification** - `Link to paper 
<http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166>`_. 

**chewBBACA 2** added the integration with Chewie-NS, which was described in the following article:
**Chewie Nomenclature Server (chewie-NS): a deployable nomenclature server for easy sharing of core and whole genome MLST schemas** - `Link to paper <https://academic.oup.com/nar/article/49/D1/D660/5929238>`_.

When using chewBBACA please cite at least one of the papers above.

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
