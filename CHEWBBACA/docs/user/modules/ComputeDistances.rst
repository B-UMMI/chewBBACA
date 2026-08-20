ComputeDistances - Compute pairwise distances based on allele calling results
=============================================================================

The ComputeDistances module computes pairwise distances based on allele calling results.

Basic Usage
-----------

::

	chewBBACA.py ComputeDistances -i /path/to/AlleleCallResultsFolder/results_alleles.tsv -o /path/to/OutputFolder

Parameters
----------

::

    -i, --input-file          (Required) Path to a TSV file containing allelic profiles determined by the AlleleCall module.

    -o, --output-directory    (Required) Path to the output directory where the process will store intermediate and final results.

    --m, --method             (Optional) Distance method used to compute the distance matrix. The module supports the hamming, 
                              jaccard and loci (number of loci not shared) methods (default: hamming).

    --outfmt, --output-format (Optional) Output format for the distance matrix (upper_triangular, lower_triangular, symmetric, table) 
                              (default: upper_triangular).

    --no-mask                 (Optional) Do not mask missing data when computing the distance matrix. This option is useful when the 
                              input profiles are already masked (default: False).

    --similarity              (Optional) Compute similarity values instead of distance values (default: False).

    --cpu, --cpu-cores        (Optional) Number of CPU cores/threads that will be used to run the process (chewie
                              resets to a lower value if it is equal to or exceeds the total number of available 
                              CPU cores/threads) (default: 1).

Outputs
-------

The process creates a TSV file containing the masked input profiles (all special classes assigned by chewBBACA are replaced by ``0`` and ``INF-`` prefixes are removed) and a TSV file containing the distance matrix or a table with the pairwise distances/similarities between all samples.
