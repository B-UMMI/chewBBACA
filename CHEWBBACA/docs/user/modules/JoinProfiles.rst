JoinProfiles - Join allele calling results from different runs
===============================================================

Basic Usage
-----------

::

	chewBBACA.py JoinProfiles -p <profiles1> <profiles2> -o <output_file>

Parameters
----------

::

    -p, --profiles    (Required) Path to files containing the results from the AlleleCall process.

    -o, --output-file (Required) Path to the output file.

    --common          (Optional) Create file with profiles for the set of common loci (default:
                      False).

.. important::
	It is necessary to pass the ``--common`` argument if the input files do not have the same
	set of loci (this option creates a new file only with the set of loci shared between all
	input files).

Outputs
-------

The process creates a TSV file with the allelic profiles of all samples in the input files.
