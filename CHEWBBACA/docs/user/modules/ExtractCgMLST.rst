ExtractCgMLST - Determine the set of loci that constitute the core genome
==========================================================================

Requirements to define a core genome MLST (cgMLST) schema
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cgMLST schemas are defined as the set of loci that are present in all strains under analysis
or, due to sequencing/assembly limitations, >95% of strains analyzed. In order to have a
robust definition of a cgMLST schema for a given bacterial species, a set of representative
strains of the diversity of a given species should be selected. Furthermore, since cgMLST
schema definition is based on pre-defined thresholds, only when a sufficient number of strains
have been analyzed can the cgMLST schema be considered stable. This number will always depend
on the population structure and diversity of the species in question, with non-recombinant
monomorphic species possibly requiring a smaller number of strais to define cgMLST schemas
than panmictic highly recombinogenic species that are prone to have large numbers of accessory
genes and mobile genetic elements. It is also important to refer that the same strategy
described here can be used to defined lineage specific schemas for more detailed analysis
within a given bacterial lineage. Also, by definition, all the loci that are not considered
core genome, can be classified as being part of an accessory genome MLST (agMLST) schema.

Determine the loci that constitute the cgMLST
:::::::::::::::::::::::::::::::::::::::::::::

Determine the set of loci that constitute the core genome based on loci presence thresholds.

Basic Usage
-----------

::

	chewBBACA.py ExtractCgMLST -i /path/to/AlleleCallResultsFolder/results_alleles.tsv -o /path/to/OutputFolder

Parameters
----------

::

    -i, --input-file       (Required) Path to the TSV file that contains the allelic profiles determined
                           by the AlleleCall module.

    -o, --output-directory (Required) Path to the directory where the process will store the output
                           files.

    --t, --threshold       (Optional) Genes that constitute the core genome must be in a proportion
                           of genomes that is at least equal to this value. Users can provide multiple
                           values to compute the core genome for multiple threshold values (default:
                           [0.95, 0.99, 1]).

    --s, --step            (Optional) The allele calling results are processed iteratively to evaluate
                           the impact of adding subsets of the results in computing the core genome. The
                           step value controls the number of profiles added in each iteration until all
                           profiles are included (default: 1).

    --r, --genes2remove    (Optional) Path to a file with a list of gene IDs to exclude from the
                           analysis (one gene identifier per line) (default: False).

    --g, --genomes2remove  (Optional) Path to a file with a list of genome IDs to exclude from
                           the analysis (one genome identifier per line) (default: False).

Outputs
-------

The output folder contains 3 files:

- ``Presence_Abscence.tsv`` - allele presence and absence matrix (1 or 0, respectively) for
  all the loci found in the ``-i`` file (excluding the loci and genomes that were flagged
  to be excluded).
- ``mdata_stats.tsv`` - total number and percentage of loci missing from each genome.
- ``cgMLST<threshold>.tsv`` - a file for each specified threshold that contains the matrix with
  the allelic profiles for the cgMLST (already excluding the list of loci and list of genomes
  passed to the ``--r`` and ``--g`` parameters, respectively).
- ``cgMLSTschema<threshold>.txt`` - a file for each specified threshold that contains the list of
  loci that constitute the cgMLST schema.
- ``cgMLST.html`` - HTML file with a line plot for the number of loci in the cgMLST per threshold.
  Also includes a black line with the number of loci present in each genome that is added to the
  analysis.

.. important::
	The ExtractCgMLST module converts/masks all non-integer classifications in the profile matrix to ``0``
	and removes all the ``INF-`` prefixes.

Example of the plot created by the ExtractCgMLST module based on the allelic profiles for 680
*Streptococcus agalactiae* genomes:

.. image:: /_static/images/cgMLST_docs.png
   :width: 900px
   :align: center

.. important::
	The ``cgMLSTschema<threshold>.txt`` file can be passed to the ``--gl`` parameter of the *AlleleCall*
	module to perform allele calling only for the loci in the cgMLST schema.

.. note::
	The matrix with allelic profiles created by the *ExtractCgMLST* process can be imported
	into `PHYLOViZ <https://online.phyloviz.net/index>`_ to visualize and explore typing results.

Workflow of the ExtractCgMLST module
::::::::::::::::::::::::::::::::::::

.. image:: /_static/images/ExtractCgMLST.png
   :width: 1000px
   :align: center
