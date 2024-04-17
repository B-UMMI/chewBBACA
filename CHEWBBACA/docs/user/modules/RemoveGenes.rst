RemoveGenes - Remove a set of loci from allele calling results
==============================================================

Basic Usage
-----------

::

	chewBBACA.py RemoveGenes -i /path/to/ProfilesFile -g /path/to/GenesListFile -o path/to/OutputFile

Parameters
----------

::

    -i, --input-file    (Required) Path to a TSV file with allelic profiles determined by the AlleleCall module.

    -g, --genes-list    (Required) Path to a file with a list of genes to remove, one identifier per line.

    -o, --output-file   (Required) Path to the output file.

    --inverse           (Optional) If provided, the genes included in the list will be kept, and all other genes
                        will be removed (default: False).

Outputs
-------

The process creates a TSV file without the columns matching the genes in the input genes list or
only with the columns matching those genes if the ``--inverse`` parameter was provided.
