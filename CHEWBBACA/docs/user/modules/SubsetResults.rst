SubsetResults - Subset allele calling results
=============================================

The SubsetResults module subsets the data in files created by chewBBACA based on a list of loci and/or sample identifiers.

Basic Usage
-----------

::

	chewBBACA.py SubsetResults -i /path/to/InputFolder -o /path/to/OutputFolder -l /path/to/loci_list.txt -s /path/to/sample_list.txt

Parameters
----------

::

    -i, --input-directory  (Required) Path to the directory containing the files to be subsetted.

    -o, --output-directory (Required) Path to the output directory.

    -l, --loci-list        (Optional) Path to a TXT/TSV file containing a list of loci to 
                           select, one locus identifier per line. If the file contains multiple 
                           columns, the loci identifiers must be in the first column (default: False).

    -s, --sample-list      (Optional) Path to a TXT/TSV file containing a list of samples to select, 
                           one sample identifier per line. If the file contains multiple columns, 
						   the sample identifiers must be in the first column (default: None).

    --inverse-loci         (Optional) If provided, the process will select the loci that are not in 
                           the input loci list (default: False).

    --inverse-samples      (Optional) If provided, the process will select the samples that are not in 
                           the input samples list(default: False).

.. important::
  The loci and sample identifiers included in the input lists must match identifiers used in the files to subset. The process warns users and exits if any of the identifers in the input lists is not in the files to subset. The process gets the complete lists of loci and sample identifiers from the ``results_alleles.tsv`` file, meaning that the ``results_alleles.tsv`` file must always be included in the input directory.

.. important::
  The module was designed to subset files with results for a single dataset. It will not work properly if the input directory includes files with results for multiple datasets.

.. note::
  You can provide the path to a folder created by the AlleleCall module or copy the files to be subsetted to a new folder and provide the path to that folder as input.

Outputs
-------

The process lists the files in the input directory and subsets the data in each file based on the provided loci and/or sample lists, saving the subsetted results in the specified output directory. The module can subset the following files:

- ``results_alleles.tsv`` (**This file is mandatory, it has to be included in the input directory**): subsets the profile matrix to extract and save the columns corresponding to the selected loci and/or rows corresponding to the selected samples.
- ``results_contigsInfo.tsv``: subsets the allele coordinates matrix to extract and save the columns corresponding to the selected loci and/or rows corresponding to the selected samples.
- ``results_statistics.tsv``: recomputes the sample summary statistics based on the subsetted profile matrix.
- ``loci_summary_stats.tsv``: recomputes the loci summary statistics based on the subsetted profile matrix.
- ``invalid_cds.txt``: identifies and extracts the rows matching CDSs identified in samples in the input sample list.
- ``cds_coordinates.tsv``:  subsets the table to extract and save the rows containing data for the subsets of loci and/or samples in the input lists.
- ``missing_classes.fasta``: identifies the records matching the loci and samples in the input lists and saves them to a new FASTA file.
- ``missing_classes.tsv``: subsets the table to extract and save the rows containing data for the subsets of loci and/or samples in the input lists.
- ``unclassified_sequences.fasta``: identifies the records matching the sample identifiers in the input sample list and saves them to a new FASTA file.
- ``paralogous_counts.tsv``: recomputes the number of times each schema locus is detected multiple times (NIPH and NIPHEM classifications) in the samples included in the dataset based on the subsetted profile matrix.
- ``paralogous_loci.tsv``: subsets the table to extract and save the rows containing data for the subsets of loci and/or samples in the input lists.
- ``novel_alleles.fasta``: identifies the records matching the loci identifiers in the input loci list and saves them to a new FASTA file.
- ``presence_absence.tsv``: subsets the presence-absence matrix to extract and save the columns corresponding to the selected loci and/or rows corresponding to the selected samples.

.. important::
  The process identifies each file type based on its basename. It cannot identify the files to subset if the filenames are different than the ones in the list above.
