UniprotFinder - Retrieve annotations for loci in a schema
=========================================================

The UniprotFinder module can be used to retrieve annotations for the loci in a schema
through requests to `UniProt's SPARQL endpoint <http://sparql.uniprot.org/sparql>`_ and through
alignment with BLASTp against UniProt's reference proteomes for a set of taxa.

Basic Usage
-----------

::

	chewBBACA.py UniprotFinder -g /path/to/SchemaFolder -o /path/to/OutputFolder -t /path/to/cds_coordinates.tsv --taxa "Species Name" --cpu 4

Parameters
----------

::

    -g, --schema-directory (Required) Path to the schema's directory.

    -o, --output-directory (Required) Output directory where the process will store intermediate
                           files and save the final TSV file with the loci annotations.

    --gl, --genes-list     (Optional) Path to a file with the list of loci in the schema that the process should
                           find annotations for (one per line, full paths or loci IDs) (default: False).

    -t, --protein-table    (Optional) Path to the TSV file with coding sequence (CDS) coordinate data,
                           "cds_coordinates.tsv", created by the CreateSchema process. (default: None).

    --bsr                  (Optional) BLAST Score Ratio value. The BSR is only used when taxa names are provided
                           to the --taxa parameter and local sequences are aligned against reference
                           proteomes downloaded from UniProt. Annotations are selected based on a BSR
                           >= than the specified value (default: 0.6).

    --cpu, --cpu-cores     (Optional) Number of CPU cores/threads that will be used to run the process
                           (chewie resets to a lower value if it is equal to or exceeds the total number
                           of available CPU cores/threads) (default: 1).

    --taxa                 (Optional) List of scientific names for a set of taxa. The process will download
                           reference proteomes from UniProt associated to taxa names that contain any
                           of the provided terms. The schema representative alleles are aligned
                           against the reference proteomes to assign annotations based on high-BSR matches
                           (default: None).

    --pm                   (Optional) Maximum number of proteome matches to report (default: 1).

    --no-sparql            (Optional) Do not search for annotations through the UniProt SPARQL endpoint.

    --no-cleanup           (Optional) If provided, intermediate files generated during process execution are not
                           removed at the end (default: False).

    --b, --blast-path      (Optional) Path to the directory that contains the BLAST executables.

Outputs
-------

The process writes a TSV file, ``schema_annotations.tsv``, with the annotations found for each
locus in the directory passed to the ``-o`` parameter. By default, the process searches for
annotations through UniProt's SPARQL endpoint, reporting the product name and UniProt URL for
local loci with an exact match in UniProt's database. If the ``cds_coordinates.tsv`` file is
passed to the ``-t`` parameter, the output file will also include the loci coordinates. The
``--taxa`` parameter receives a set of taxa names and searches for reference proteomes that match
the provided terms. The reference proteomes are downloaded and the process aligns the loci
representative alleles against the reference proteomes to include the product and gene name
in the reference proteomes.
