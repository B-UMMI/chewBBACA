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
includes functionalities to annotate the schema loci, compute the set of loci that constitute the core genome for a given dataset, 
and generate interactive reports for schema and allele calling results evaluation to enable an intuitive analysis of the results 
in surveillance and outbreak detection settings or population studies. Pre-defined cg/wgMLST schemas can be downloaded from 
`Chewie-NS <https://chewbbaca.online/>`_ or adapted from other cg/wgMLST platforms.

The general workflow of chewBBACA is represented in the following image:

.. image:: /_static/images/Overview.png
   :width: 1200px
   :align: center

Citation
--------

chewBBACA has been published (version 2.0.5 at the time) in Microbial Genomics under the title:
**chewBBACA: A complete suite for gene-by-gene schema creation and strain identification** - `Link to paper 
<http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166>`_. 

When using chewBBACA please use the following citation:

::

  Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. doi:10.1099/mgen.0.000166

Licensing
---------

This project is licensed under the `GPLv3 license 
<https://github.com/B-UMMI/Nomenclature_Server_docker_compose/blob/master/LICENSE>`_.
The source code of chewBBACA is available at `<https://github.com/B-UMMI/chewBBACA>`_.

Funding
-------

- INNUENDO project co-funded by the European Food Safety Authority (EFSA), grant agreement
  GP/EFSA/AFSCO/2015/01/CT2 ("New approaches in identifying and characterizing microbial and
  chemical hazards"). The conclusions, findings, and opinions expressed in this review paper
  reflect only the view of the authors and not the official position of the European Food Safety
  Authority (EFSA).
- ONEIDA project (LISBOA-01-0145-FEDER-016417) co-funded by FEEI - "Fundos Europeus Estruturais
  e de Investimento" from "Programa Operacional Regional Lisboa 2020" FCT - "Fundação para a
  Ciência e a Tecnologia".
- BacGenTrack (TUBITAK/0004/2014) [FCT/ Scientific and Technological Research Council of Turkey
  (Türkiye Bilimsel ve Teknolojik Araşrrma Kurumu, TÜBİTAK)]

.. image:: http://i.imgur.com/XhvagNV.png
   :width: 500px
   :align: center
