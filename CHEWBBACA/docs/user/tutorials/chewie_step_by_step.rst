Step-by-step Tutorial
=====================

Objective
:::::::::

The objective of this tutorial is to illustrate a complete workflow for creating a wgMLST and a
cgMLST schema for a collection of 714 *Streptococcus agalactiae* genomes (32 complete genomes
and 682 draft genome assemblies deposited on the NCBI databases) by providing step-by-step
instructions and displaying the obtained outputs.

Please start by going through the following steps:

- Install chewBBACA. Check :doc:`Installation </user/getting_started/installation>` for instructions
  on how to install chewBBACA.
- Download the ZIP file with the datasets and expected results for the tutorial `here <https://zenodo.org/record/7129300#.Yzb26_fTUUE>`_.
- Uncompress the ZIP file (this will create a folder named ``chewBBACA_tutorial`` that has all
  the necessary data).

.. important::
	The paths and commands used in this tutorial assume that the working directory is the top-level
	directory of the folder **chewBBACA_tutorial** created after uncompressing the ZIP file.
	The commands should be modified if they are executed from a different working directory.
	The expected results for each section were included in the ``expected_results`` folder
	for reference (each subfolder has the name of one of the sections).

Metadata about the NCBI genomes used in this tutorial is on the TSV file ``genomes/GBS_NCBI_metadata.tsv``.

chewBBACA includes Prodigal training files for several species, including for
*Streptococcus agalactiae*. You can check the list of available training files
`here <https://github.com/B-UMMI/chewBBACA/raw/master/CHEWBBACA/prodigal_training_files/>`_. We
have included the training file for *Streptococcus agalactiae*,
``Streptococcus_agalactiae.trn``, in the tutorial data.

.. note::
	The execution times reported in this tutorial were obtained for a DELL XPS13 (10th
	Generation Intel® Core™ i7-10710U Processor - 12MB Cache, up to 4.7 GHz, using 6 cores).
	Using a computer with less powerful specifications can increase the duration
	of the analyses.  

Schema creation
:::::::::::::::

We will start by creating a wgMLST schema based on 32 *Streptococcus agalactiae* complete
genomes (32 genomes with a level of assembly classified as complete genome or chromossome)
available at NCBI. Uncompress the ``genomes/GBS_32_complete_genomes.zip`` file
to create the folder ``genomes/complete_genomes``. To create the wgMLST schema, run the following command:  

::

	chewBBACA.py CreateSchema -i genomes/complete_genomes/ -o tutorial_schema --ptf Streptococcus_agalactiae.trn --cpu 6

The schema seed will be available at ``tutorial_schema/schema_seed``. We passed the value ``6`` to
the ``--cpu`` parameter to use 6 CPU cores, but you should pass a value based on the
specifications of your machine. In our system, the process took 56 seconds to complete
resulting on a wgMLST schema with 3128 loci. At this point the schema is defined as a set of
loci each with a single representative allele.

Allele calling
::::::::::::::

The next step is to perform allele calling with the wgMLST schema created in the previous step
for the **32** complete genomes. The allele call step determines the allelic profiles of the
analyzed strains, identifying known and novel alleles in the analyzed genomes. Novel alleles
are assigned an allele identifier and added to the schema. To perform allele call, run the
following command:

::

	chewBBACA.py AlleleCall -i genomes/complete_genomes/ -g tutorial_schema/schema_seed -o results32_wgMLST --cpu 6

The allele call used the default BSR threshold of ``0.6`` and took approximately 17
minutes to complete (an average of 32 seconds per genome). The allele call identified 14,720
novel alleles and added those alleles to the schema, increasing the number of alleles in the
schema from 3,128 to 17,848.

Paralog detection
:::::::::::::::::

The next step in the analysis is to determine if some of the loci can be considered paralogs
based on the result of the wgMLST allele calling. The *AlleleCall* module returns a list of
Paralogous genes in the ``paralogous_counts.tsv`` file that can be found on the
``results32_wgMLST/results_<datestamp>`` folder. The ``paralogous_counts.tsv`` file contains a set
of 20 loci that were identified as possible paralogs. These loci should be removed from the schema
due to the potential uncertainty in allele assignment. To remove the set of 20 paralogous loci
from the allele calling results, run the following command:

::

	chewBBACA.py RemoveGenes -i results32_wgMLST/results_<datestamp>/results_alleles.tsv -g results32_wgMLST/results_<datestamp>/paralogous_counts.tsv -o results32_wgMLST/results_<datestamp>/results_alleles_NoParalogs.tsv

This will remove the columns matching the 20 paralogous loci from the allele calling results and
save the allelic profiles into the ``results_alleles_NoParalogs.tsv`` file (the new file contains
allelic profiles with 3108 loci).

cgMLST schema determination
:::::::::::::::::::::::::::

We can now determine the set of loci in the core genome based on the allele calling results.
The set of loci in the core genome is determined based on a threshold of loci presence in the
analysed genomes. We can run the *TestGenomeQuality* module to determine the impact of several
threshold values on the number of loci in the core genome.

::

	chewBBACA.py TestGenomeQuality -i results32_wgMLST/results_<datestamp>/results_alleles_NoParalogs.tsv -n 13 -t 200 -s 5 -o results32_wgMLST/results_<datestamp>/genome_quality_32

The process will automatically open a HTML file with the following plot:

.. image:: https://i.imgur.com/uf3Hygd.png
	:width: 1400px
	:align: center

Interactive version `here <http://im.fm.ul.pt/chewBBACA/GenomeQual/GenomeQualityPlot_complete_genomes.html>`_.

A set of **1136 loci** were found to be present in all the analyzed complete genomes, while
**1267 loci** were present in at least 95%. For further analysis only the **1267** loci present
in at least 95% of the complete genomes will be used. We selected that threshold value to account
for loci that may not be identified due to sequencing coverage and assembly problems.

We can run the *ExtraCgMLST* module to quickly determine the set of loci in the core genome at 95%.

::

	chewBBACA.py ExtractCgMLST -i results32_wgMLST/results_<datestamp>/results_alleles_NoParalogs.tsv -o results32_wgMLST/results_<datestamp>/cgMLST_95 --t 0.95

The list with the 1267 loci in the core genome at 95% is in the
``results32_wgMLST/results_<datestamp>/cgMLST_95/cgMLSTschema.txt`` file. This file can be passed
to the ``--gl`` parameter of the AlleleCall process to perform allele calling only for the set of
genes that constitute the core genome.

Allele call for 682 *Streptococcus agalactiae* assemblies
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

682 assemblies of *Streptococcus agalactiae* available on NCBI were downloaded (03-08-2016) and
analyzed with `MLST <https://github.com/tseemann/mlst>`_ in order to exclude possibly mislabeled
samples as *Streptococcus agalactiae*. Out of the 682 genomes, 2 (GCA_000323065.2_ASM32306v2 and
GCA_001017915.1_ASM101791v1) were detected as being of a different species/contamination and
were removed from the analysis. Uncompress the ``genomes/GBS_680_genomes.zip`` file to create a
folder named ``GBS_Aug2016``.

Allele call was performed on the *bona fide* *Streptococcus agalactiae* **680 genomes** using the
**1267 loci** that constitute the core genome at 95%. Paralog detection found no paralog loci.

::

	chewBBACA.py AlleleCall -i genomes/GBS_Aug2016/ -g tutorial_schema/schema_seed --gl results32_wgMLST/results_<datestamp>/cgMLST_95/cgMLSTschema.txt -o results680_cgMLST --cpu 6

The process took approximately 39 minutes to complete (an average of 3.4 secs per genome).

Evaluate genome quality
:::::::::::::::::::::::

We can now concatenate the cgMLST results for the 32 complete genomes with the cgMLST results
for the 680 genomes to have all the results in a single file. To concatenate the allelic profiles
of both analyses run the following command:

::

	chewBBACA.py JoinProfiles -p results32_wgMLST/results_<datestamp>/cgMLST_95/cgMLST.tsv results680_cgMLST/results_<datestamp>/results_alleles.tsv -o cgMLST_all.tsv


The concatenated file was analyzed in order to assess the cgMLST allele quality attribution for
all the genomes.

::

	chewBBACA.py TestGenomeQuality -i cgMLST_all.tsv -n 13 -t 300 -s 5

.. image:: https://i.imgur.com/m1OSycz.png
	:width: 1400px
	:align: center

Interactive version `here <http://im.fm.ul.pt/chewBBACA/GenomeQual/GenomeQualityPlot_all_genomes.html>_`.

While the number of loci present in 95% of genomes remains virtually constant at around **1200**
loci, considering all or most of the genomes (90%<x≤100%) the number of loci present is lower and
presents some variation when specific genomes are removed from the analysis.

We selected the results at the threshold of 25 for further analysis. Although this selection is
somewhat arbitrary, when moving to a lower threshold there is step increase in the number of loci
present in 95% and 99% of genomes that could represent the exclusion of a more divergent clade
from the analysis. Furthermore, at the threshold of 25 there is an acceptable number of loci
present in all considered genomes (650 genomes/440 loci), which we felt would afford a good
discriminatory power.

The genomes that were removed at each threshold are indicated in the file
``expected_results/Evaluate_genome_quality/removedGenomes.txt`` and a file
``expected_results/Evaluate_genome_quality/removedGenomes_25.txt`` was created
with only the genomes removed at the 25 threshold.

The following command creates a folder named ``cgMLST_25`` and saves the cgMLST schema selected
at the chosen threshold to the file ``cgMLST_25/cgMLST.tsv``.

::

	chewBBACA.py ExtractCgMLST -i cgMLST_all.tsv -o cgMLST_25 --g expected_results/Evaluate_genome_quality/removedGenomes_25.txt

Minimum Spanning Tree
:::::::::::::::::::::

The file ``cgMLST_25/cgMLST.tsv`` was uploaded to `Phyloviz online <https://online.phyloviz.net>`_
and can be accessed `here <https://online.phyloviz.net/main/dataset/share/cfab1610a3ca3a80cf9c139e436ce741fc5fa29dcc5aeb3988025491d7194044fc73f5284eafad8356322fb0e29e50d6e06d5808ae369a2b37d1ece96e4e716d8d7eeb5c85a5a30c5d3d63bf014643013fa981108bd5bfbacf0a145ab41656a9a67c489b878cb0aa9f2de534ee81b201e198>`_.

Genome Quality analysis
:::::::::::::::::::::::

Since the quality of the used assemblies was not confirmed, it is possible that some of the
assemblies included were of low quality. A general analysis of the assemblies show a N50
variation that ranges from 8,055 to over 2.2M, while the number of contigs ranges between
1 and 553. These results made us suspect that the quality of the genomes could have affected
the allele call results and consequently caused a significant drop in the number of loci
detected as present in all genomes.  

As stated previously, to obtain the cgMLST schema, some genomes (n=62) had to be removed since
they were extremes cases of missing data. In order to assess the possible reason for their poor
allele call performance, two plots were built. The removed genomes were then highlighted and
dashed lines were drawn linking the values for the same genomes.

The first plot represents the total number of bp in contigs with a size >10 kbp and the N50 of
the assemblies, sorted by decreasing values.

.. image:: http://i.imgur.com/I0fNqtd.png
	:width: 1400px
	:align: center

The second plot represents the total number of contigs and the number of contigs >10kbp.

.. image:: http://i.imgur.com/fabxi0Z.png
	:width: 1400px
	:align: center

Interactive version `here <http://im.fm.ul.pt/chewBBACA/GenomeQual/AssemblyStatsStack.html>`_

At first sight, most of the removed genomes (56/62) were located on the lower range of N50 and
bp in contigs >10 kbp (fig.3) and the higher number of contigs (fig.4).

The 5 genomes that were outside this pattern were individually checked:

- `GCA_000186445.1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000186445.1>`_ - 21 contigs
  but only 1 is above 10k (Scaffold with lot of Ns, 134 real contigs)
- `GCA_000221325.2 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000221325.2>`_ - NCBI curated
  it out of RefSeq because it had a genome length too large
- `GCA_000427055.1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000427055.1>`_ - NCBI curated
  it out of RefSeq because it had many frameshifted proteins
- `GCA_000289455.1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000289455.1>`_ - No ST found.
  We concluded the assembly has a problem but we have not yet identified it.
- `GCA_000288835.1 <https://www.ncbi.nlm.nih.gov/assembly/GCA_000288835.1>`_ - NCBI curated
  it out of RefSeq because it had many frameshifted proteins.
