TestGenomeQuality - Analyze your allele call output to refine schemas
=====================================================================

The quality of draft genomes can impact profoundly the MLST schema definition. Therefore, chewBBACA offers a way to test the impact
of each genome on the number of loci selected for inclusion in the final schema, based on the number of missing loci when considering
the original schema being tested.

The algorithm description is the following:

	1) For a given genome, let *nl* be the number of loci that are not present (neither EXC nor INF, see the :doc:`AlleleCall section </user/modules/AlleleCall>`
	   for more detail in loci classification) but are found in at least a 95% of the genomes in analysis;
	2) For a specific exclusion threshold (*et*, maximum number of absent loci), genomes that have *nl > et* are excluded from the analysis in the next iteration.
	   If no genomes are removed, do not proceed to step 3;
	3) Return to step 1. The total number of genomes to be used in step 1 will be recalculated leaving out the genomes excluded in step 2. The *et* value
	   (starting at 0) is incremented by a *s* value previously defined by the user.

Basic Usage
-----------

::

	chewBBACA.py TestGenomeQuality -i /path/to/AlleleCall/results/results_alleles.tsv -o /path/to/OutputFolderName -n 12 -t 200 -s 5

Parameters
----------

::

    -i, --input-files      (Required) Path to file with a matrix of allelic profiles (default: None).

    -o, --output-directory (Required) Path to the output directory that will store output files (default: None).

    -n, --max-iteration    (Required) Maximum number of iterations. Each iteration removes a set of genomes over
	                       the defined threshold (-t) and recalculates loci presence percentages (default: None).

    -t, --max-threshold    (Required) This threshold represents the maximum number of missing loci allowed, for
	                       each genome independently, before removing the genome (default: None).

    --s --step             (Required) Step to add to each threshold (default: 5).

Outputs
-------

The outputs consist of an interactive plot and a ``.../outputDirectory/removedGenomes.txt`` text file. The latter output contains the
list of genomes removed in each exclusion threshold (*et*) for loci that must be in at least the 95% of the genomes in analysis.
An example of the interactive plot output can be seen `here <http://im.fm.ul.pt/chewBBACA/GenomeQual/GenomeQualityPlot_all_genomes.html>`_.
This plot shows the number of genomes used with each exclusion threshold (*et*) considering loci that must be in ≥ 95% of the genomes under
analysis, together with the number of loci present in 100%, ≥99.5%, ≥99% and ≥95% of the genomes considered in each iteration. This
example uses a set of 712 genomes and a schema of 1264 loci. The -n 12 -t 300 -s 5 parameters were used. Based on the plot and specific
knowledge of the species under analysis, the user can decide which exclusion threshold value (*et*) to use for genome exclusion in the
cgMLST schema refinement. This might be important since poor quality genome assemblies/annotations can increase the false negatives and
introduce artificial loci.
