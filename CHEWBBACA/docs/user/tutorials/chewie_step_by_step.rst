Step-by-step Tutorial
=====================

Objective
:::::::::

The objective of this tutorial is to illustrate a complete workflow for creating a wgMLST and a
cgMLST schema for a collection of 712 *Streptococcus agalactiae* genomes (32 complete genomes
and 680 draft genome assemblies deposited on the NCBI databases) by providing step-by-step
instructions and displaying the obtained outputs.

Please start by going through the following steps:

- Install chewBBACA. Check :doc:`Installation </user/getting_started/installation>` for instructions
  on how to install chewBBACA.
- Download the ZIP file with the datasets and expected results for the tutorial `here <https://zenodo.org/records/10694715>`_.
- Uncompress the ZIP file (this will create a folder named ``chewBBACA_tutorial`` that has all
  the necessary data).

.. important::
	The paths and commands used in this tutorial assume that the working directory is the top-level
	directory of the folder **chewBBACA_tutorial** created after uncompressing the ZIP file.
	The commands should be modified if they are executed from a different working directory.
	The expected results for each section were included in the ``expected_results`` folder
	for reference (each subfolder has the name of one of the sections).

Metadata about the NCBI genomes used in this tutorial is available on the TSV file ``genomes/GBS_NCBI_metadata.tsv``.

chewBBACA includes Prodigal training files for several species, including for
*Streptococcus agalactiae*. You can check the list of available training files
`here <https://github.com/B-UMMI/chewBBACA/raw/master/CHEWBBACA/prodigal_training_files/>`_. We
have included the training file for *Streptococcus agalactiae*,
``Streptococcus_agalactiae.trn``, in the tutorial data.

.. note::
	The execution times reported in this tutorial were obtained for a Lenovo Legion Pro i5 16
	(i7-13700HX, using 6 CPU threads). The execution time will vary depending on the specifications
	of the computer used to perform the analyses and on the value passed to the ``--cpu`` parameter.

Schema creation
:::::::::::::::

We will start by creating a wgMLST schema based on 32 *Streptococcus agalactiae* complete
genomes (32 genomes with a level of assembly classified as complete genome or chromossome)
available at the NCBI. Uncompress the ``genomes/sagalactiae_32_complete_genomes.zip`` file
to create the folder ``genomes/sagalactiae_32_complete_genomes``. To create the wgMLST schema,
run the following command:  

::

	chewBBACA.py CreateSchema -i genomes/sagalactiae_32_complete_genomes/ -o tutorial_schema --ptf Streptococcus_agalactiae.trn --cpu 6

The schema seed will be available at ``tutorial_schema/schema_seed``. We passed the value ``6`` to
the ``--cpu`` parameter to use 6 CPU cores, but you should pass a value based on the
specifications of your machine. In our system, the process took 17 seconds to complete
resulting on a wgMLST schema with 3,127 loci. At this point the schema is defined as a set of
loci each with a single representative allele.

Allele calling
::::::::::::::

The next step is to perform allele calling with the wgMLST schema created in the previous step
for the 32 complete genomes. The allele call step determines the allelic profiles of the
analyzed strains, identifying known and novel alleles in the analyzed genomes. Novel alleles
are assigned an allele identifier and added to the schema. To perform allele call, run the
following command:

::

	chewBBACA.py AlleleCall -i genomes/sagalactiae_32_complete_genomes -g tutorial_schema/schema_seed -o results32_wgMLST --cpu 6

The allele call used the default BLAST Score Ratio (BSR) value of ``0.6`` and took 30s to complete. The allele call identified 14,703
novel alleles and added those alleles to the schema, increasing the number of alleles in the schema from 3,127 to 17,830.

Paralog detection
:::::::::::::::::

The next step in the analysis is to determine if some of the loci can be considered paralogs
based on the result of the wgMLST allele calling. The *AlleleCall* module returns a list of
Paralogous genes in the ``paralogous_counts.tsv`` file that can be found on the
``results32_wgMLST`` folder. The ``paralogous_counts.tsv`` file contains a set
of 12 loci that were identified as possible paralogs. These loci should be removed from the schema
due to the potential uncertainty in allele assignment. To remove the set of 12 paralogous loci
from the allele calling results, run the following command:

::

	chewBBACA.py RemoveGenes -i results32_wgMLST/results_alleles.tsv -g results32_wgMLST/paralogous_counts.tsv -o results32_wgMLST/results_alleles_NoParalogs.tsv

This will remove the columns matching the 12 paralogous loci from the allele calling results and
save the allelic profiles into the ``results32_wgMLST/results_alleles_NoParalogs.tsv`` file (the new file contains
allelic profiles with 3,115 loci).

cgMLST schema determination
:::::::::::::::::::::::::::

We can now determine the set of loci in the core genome based on the allele calling results.
The set of loci in the core genome is determined based on a threshold of loci presence in the
analysed genomes. We can run the *ExtractCgMLST* module to determine the set of loci in
the core genome for the loci presence thresholds of 95%, 99% and 100%.

::

	chewBBACA.py ExtractCgMLST -i results32_wgMLST/results_alleles_NoParalogs.tsv -o results32_wgMLST/cgMLST

The *ExtractCgMLST* module creates a HTML file (available at ``results32_wgMLST/cgMLST/cgMLST.html``) with
an interactive line plot that displays the number of loci in the cgMLST per threshold value and the number
of loci in each genome added to the analysis.

.. image:: /_static/images/cgMLST_tutorial.png
   :width: 900px
   :align: center

.. note::
	The ExtractCgMLST module converts/masks all non-integer classifications in the profile matrix to ``0``
	and removes all the ``INF-`` prefixes.

We selected the threshold of 95% to account for loci that may not be identified due to sequencing
coverage and assembly problems. The list with the 1,271 loci in the core genome at 95% is in the
``results32_wgMLST/cgMLST/cgMLSTschema95.txt`` file. This file can be passed
to the ``--gl`` parameter of the AlleleCall process to perform allele calling only for the set of
genes that constitute the core genome.

Allele call for 680 *Streptococcus agalactiae* assemblies
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

682 assemblies of *Streptococcus agalactiae* available on NCBI were downloaded (03-08-2016) and
analyzed with `MLST <https://github.com/tseemann/mlst>`_ in order to exclude possibly mislabeled
samples as *Streptococcus agalactiae*. Out of the 682 genomes, 2 (GCA_000323065.2_ASM32306v2 and
GCA_001017915.1_ASM101791v1) were detected as being of a different species/contamination and
were removed from the analysis. Uncompress the ``genomes/sagalactiae_680_draft_genomes.zip`` file to create a
folder named ``sagalactiae_680_draft_genomes``.

Allele call was performed on the *bona fide* *Streptococcus agalactiae* **680 genomes** using the
**1,271 loci** that constitute the core genome at 95%.

::

	chewBBACA.py AlleleCall -i genomes/sagalactiae_680_draft_genomes/ -g tutorial_schema/schema_seed --gl results32_wgMLST/cgMLST/cgMLSTschema95.txt -o results680_cgMLST --cpu 6

The process took 1m32s to complete and added 23,767 novel alleles to the schema. Paralog detection found no paralogous loci.

Redetermination of the cgMLST
:::::::::::::::::::::::::::::

We can now concatenate the cgMLST results for the 32 complete genomes with the cgMLST results
for the 680 genomes to have all the results in a single file. To concatenate the allelic profiles
of both analyses run the following command:

::

	chewBBACA.py JoinProfiles -p results32_wgMLST/cgMLST/cgMLST95.tsv results680_cgMLST/results_alleles.tsv -o cgMLST_712.tsv

We also redetermined the cgMLST based on the allele calling results for this more diverse set of
strains:

::

	chewBBACA.py ExtractCgMLST -i cgMLST_712.tsv -o cgMLST_712

The number of loci present in 95% of genomes based on the 712 assemblies is 1,194, a slight decrease
from the number of loci present in 95% of the 32 genomes used for schema creation.

Evaluate genome quality
:::::::::::::::::::::::

One important factor that was not evaluated, and that can greatly affect the cgMLST determination,
is the quality of the genome assemblies. Since the quality of the used assemblies was not confirmed,
it is possible that some of the assemblies included were of low quality. A general analysis of the
assemblies (available at ``genomes/sagalactiae_assembly_stats.tsv``) shows a N50 variation that ranges from 8,055
to over 2.2M, while the number of contigs ranges between 1 and 553. These results made us suspect
that the quality of the genomes could have affected the allele call results and consequently caused
a significant drop in the number of loci that constitute the cgMLST. We defined a set of minimum quality
criteria to select high quality genome assemblies, that are the following:

- Less than 150 contigs.
- Genome size between 1,674,000 and 2,512,000 bases (defined according to the species genome size values provided by the NCBI on 16-12-2022 and available `here <https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt.gz>`_).
- Less than 1,000 N bases.
- Less than 5% missing loci from the cgMLST (64 loci).

We identified 65 genome assemblies that did not meet the minimum quality criteria and 2 genomes that the NCBI excluded from RefSeq,
`GCA_000221325.2 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000221325.2>`_ and
`GCA_000427055.1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000427055.1>`_, due to ``genome length too large`` and
``many frameshifted proteins``, respectively (the list of excluded genome assemblies is available at
``expected_results/Evaluate_genome_quality/excluded_genomes.txt``).
We used the following command to recompute the cgMLST:

::

	chewBBACA.py ExtractCgMLST -i cgMLST_712.tsv -o cgMLST_645 --g expected_results/Evaluate_genome_quality/excluded_genomes.txt

The determined cgMLST at 95% includes 1,248 loci, an additional 54 loci (~+4% of the previously defined cgMLST).

Minimum Spanning Tree
:::::::::::::::::::::

You can upload the file ``cgMLST_645/cgMLST95.tsv`` and any of the associated metadata to `PHYLOViZ Online <https://online.phyloviz.net>`_
to visualize a Minimum Spanning Tree and perform various dataset operations that allow you to explore and analyse the results generated
during this tutorial. PHYLOViZ Online considers all classifications when computing the distances between samples. This means that classifications
such as ``ASM`` and ``LNF`` will be treated in the same way as valid allele identifiers. If you want to define that your profile matrix includes missing
data that should not count for cgMLST analysis, you can upload the profile matrix created by the ExtractCgMLST module, which converts/masks
all non-integer classifications to ``0`` and removes ``INF-`` prefixes, to `PHYLOViZ Online 2 <https://online2.phyloviz.net/index>`_,
which allows users to specify that a profile matrix includes missing data. Please note that PHYLOViZ Online 2 is not in production phase,
and you might run into some feature bugs.
