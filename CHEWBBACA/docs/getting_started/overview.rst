Overview
========

About chewBBACA
---------------

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be
thought of as "Comprehensive and  Highly Efficient Workflow" but at this point still it
needs a bit of work to make that claim so we just add "chew" to add extra coolness to
the software name. BSR stands for BLAST Score Ratio as proposed by `Rasko DA et al 
<http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2>`_. 

chewBBACA is a comprehensive pipeline including a set of functions for the creation and
validation of whole genome and core genome MultiLocus Sequence Typing (wg/cgMLST) schemas,
providing an allele calling algorithm based on Blast Score Ratio that can be run in multiprocessor 
settings and a set of functions to visualize and validate allele variation in the loci.
chewBBACA performs the schema creation and allele calls on complete or draft genomes resulting
from de novo assemblers.

.. important:: chewBBACA only works with **python 3** (automatic testing for Python 3.8 and Python 3.9
               with GitHub Actions).

.. warning:: We strongly recommend that users install and use **BLAST 2.9.0+** with chewBBACA, as
             chewBBACA's processes have been extensively tested with that version of BLAST.

The general workflow of chewBBACA is represented in the following image:

.. image:: http://i.imgur.com/aqNsv7i.png

chewBBACA is developed by the `Molecular Microbiology and Infection Unit (UMMI) 
<http://im.fm.ul.pt>`_ at the 
`Instituto de Medicina Molecular João Lobo Antunes 
<https://imm.medicina.ulisboa.pt/>`_ and `Faculdade de Medicina, Universidade de Lisboa 
<https://www.medicina.ulisboa.pt/>`_.

Citation
--------

chewBBACA has been published (version 2.0.5 at the time) in Microbial Genomics under the title:
**chewBBACA: A complete suite for gene-by-gene schema creation and strain identification** - `Link to paper 
<http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166>`_. 

When using chewBBACA please use the following citation:

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. [doi:10.1099/mgen.0.000166](doi:10.1099/mgen.0.000166)

Licensing
---------

This project is licensed under the `GPLv3 license 
<https://github.com/B-UMMI/Nomenclature_Server_docker_compose/blob/master/LICENSE>`_.
The source code of chewBBACA is available at `<https://github.com/B-UMMI/chewBBACA>`_.

Contacts
--------

For any issues contact the development team at imm-bioinfo@medicina.ulisboa.pt or 
open an issue at the chewBBACA Github `repository <https://github.com/B-UMMI/chewBBACA>`_.

Funding
-------

- INNUENDO project (https://www.innuendoweb.org) co-funded by the European Food Safety
  Authority (EFSA), grant agreement GP/EFSA/AFSCO/2015/01/CT2 (“New approaches in
  identifying and characterizing microbial and chemical hazards”).” The conclusions,
  findings, and opinions expressed in this review paper reflect only the view of the
  authors and not the official position of the European Food Safety Authority (EFSA).
- ONEIDA project (http://www.itqb.unl.pt/oneida) (LISBOA-01-0145-FEDER-016417) co-funded
  by FEEI - "Fundos Europeus Estruturais e de Investimento" from "Programa Operacional
  Regional Lisboa 2020" FCT - "Fundação para a Ciência e a Tecnologia"
- BacGenTrack (http://darwin.phyloviz.net/dokuwiki/doku.php?id=pip:tubitak) (TUBITAK/0004/2014)
  [FCT/ Scientific and Technological Research Council of Turkey (Türkiye Bilimsel ve Teknolojik
  Araşrrma Kurumu, TÜBİTAK)]

.. image:: http://i.imgur.com/XhvagNV.png
