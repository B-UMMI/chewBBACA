JoinProfiles - Join allele calling results from different runs
===============================================================

Basic Usage
-----------

::

	chewBBACA.py JoinProfiles -p <profiles1> <profiles2> -o <output_file>

Parameters
----------

::

    -p, --profiles    (Required) Paths to the files containing allelic profiles determined by
                      the AlleleCall module. It is possible to provide any number of files.
                      The results must have been determined with the same schema and share
                      all the loci or a subset of the loci if using the --common parameter.

    -o, --output-file (Required) Path to the output file.

    --common          (Optional) Merge the results based on the subset of loci shared between
                      all files. (default: False).

.. important::
	It is necessary to pass the ``--common`` argument if the input files do not have the same
	set of loci (this option creates a new file only with the set of loci shared between all
	input files).

Outputs
-------

The process creates a TSV file with the allelic profiles of all samples in the input files.
