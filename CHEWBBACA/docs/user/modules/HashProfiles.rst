HashProfiles - Hash allelic profiles
====================================

The HashProfiles module hashes allelic profiles determined with chewBBACA, allowing users operating under stricter data privacy policies to share and compare their results without having to share allele data such as allele identifiers and sequences.

Basic Usage
-----------

::

	chewBBACA.py HashProfiles -i /path/to/InputFile -g /path/to/SchemaFolder -o /path/to/OutputFolder

Parameters
----------

::

    -i, --input-directory  (Required) Path to the TSV file that contains the allelic profiles determined by the
                           AlleleCall module.

    -g, --schema-directory (Required) Path to the schema's directory to get the allele sequences and compute the
                           hashes.

    -o, --output-directory (Required) Path to the output directory.

    --hash-type            (Optional) Hashing algorithm used to hash the profiles. The hashing algorithms 
                           implemented in the hashlib and zlib Python libraries are supported.

    --nrows                (Optional) Divide the input file into chunks of this many rows to process larger files 
                           more efficiently.

    --cpu, --cpu-cores     (Optional) Number of CPU cores/threads that will be used to run the process (chewie
                           resets to a lower value if it is equal to or exceeds the total number of
                           available CPU cores/threads).

.. important::
  The HashProfiles module computes the hashes based on the allele sequences matching the allele identifiers in the profiles. The schema passed to the process must include all the alleles referenced in the profiles. When using the ``--no-inferred`` option for allele calling, newly identified alleles are not added to the schema, and it will not be possible to compute the hashes for those alleles.

Outputs
-------

The process creates a single TSV file in the output directory with the same name as the input file, but with the suffix ``_hashed``. The output file contains the same columns as the input file, but with the allele identifiers replaced by their corresponding hashes.
