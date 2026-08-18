MergeResults - Merge results files created by chewBBACA
=======================================================

The MergeResults module merges results files created by chewBBACA based on the same or a common set of loci.

Basic Usage
-----------

::

	chewBBACA.py MergeResults -i /path/to/InputFolder1 /path/to/InputFolder2 -o /path/to/OutputFolder

Parameters
----------

::

    -i, --input-directory  (Required) Paths to the directories containing the results files created by chewBBACA.
                           The results must have been determined with the same schema and share all the loci or 
                           a subset of the loci if using the --common parameter.

    -o, --output-directory (Required) Path to the output directory.

    --common               (Optional) Merge the results based on the subset of loci shared between all inputs (default: False).

.. important::
  The results to merge must have been determined with the same schema so that the loci identifiers are the same in all input directories. If the results were determined with different schemas, the module will raise an exception because the loci identifiers will not match.

.. note::
  If the ``--common`` parameter is used, the module will identify the set of loci shared between the results in all input directories and merge the results based on that set of loci. The module calls the ``SubsetResults`` module to subset the results in each input directory based on the set of shared loci before merging the results.

Outputs
-------

The process lists the files in the input directories and groups them by file type, then it merges the files matching each file type and saves the merged results to the output directory. The module can merge the following files:

- ``results_alleles.tsv`` (**This file is mandatory, it has to be included in the input directory**): concatenates the allelic profile matrices from the input directories and saves the merged matrix to the output directory.
- ``results_contigsInfo.tsv``: concatenates the allele coordinates matrices from the input directories and saves the merged matrix to the output directory.
- ``results_statistics.tsv``: concatenates all files with samples statistics from the input directories and saves the merged TSV file to the output directory.
- ``loci_summary_stats.tsv``: sums the values in the loci summary statistics files from the input directories and saves the merged TSV file to the output directory.
- ``invalid_cds.txt``: concatenates the files with the list of invalid CDSs from the input directories and saves the merged list to the output directory.
- ``cds_coordinates.tsv``:  concatenates all files with CDS coordinates from the input directories and saves the merged TSV file to the output directory.
- ``missing_classes.fasta``: concatenates the FASTA files in all input directories and saves the merged FASTA file to the output directory.
- ``missing_classes.tsv``: concatenates all TSV files with information about the CDSs assigned special classifications and saves the merged TSV file to the output directory.
- ``unclassified_sequences.fasta``: concatenates the FASTA files in all input directories and saves the merged FASTA file to the output directory.
- ``paralogous_counts.tsv``: sums the values for all loci for which oaralogous matches were identified and saves the merged TSV file to the output directory.
- ``paralogous_loci.tsv``: concatenates all TSV files with information about paralogous matches and saves the merged TSV file to the output directory.
- ``novel_alleles.fasta``: concatenates the FASTA files in all input directories and saves the merged FASTA file to the output directory.
- ``presence_absence.tsv``: concatenates the presence-absence matrices from the input directories and saves the merged matrix to the output directory.

.. important::
  The process identifies each file type based on its basename. It cannot identify the files to merge if the filenames are different than the ones in the list above.
