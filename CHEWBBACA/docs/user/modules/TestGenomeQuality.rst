TestGenomeQuality - Analyze your allele call output to refine schemas
=====================================================================

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

The output is a HTML file with a plot with all thresholds and a ``removedGenomes.txt`` file with
information about which genomes were removed per threshold when it reaches a stable point (no more genomes are removed).

Example of an output can be seen [here](http://im.fm.ul.pt/chewBBACA/GenomeQual/GenomeQualityPlot_all_genomes.html).
The example uses a set of 714 genomes and a schema consisting of 3266 loci with ``-n 12 -t 300 -s 5`` passed to arguments.
