UniprotFinder - Retrieve annotations for loci in a schema
=========================================================

The UniprotFinder module can be used to retrieve annotations for the loci in a schema
through requests to `UniProt's SPARQL endpoint <http://sparql.uniprot.org/sparql>`_ and through
alignment with BLASTp against UniProt's reference proteomes for a set of taxa.

Basic Usage
-----------

::

	chewBBACA.py UniprotFinder -g /path/to/SchemaName -o /path/to/OutputFolderName -t /path/to/cds_coordinates.tsv --taxa "Species Name" --cpu 4

Parameters
----------

::

    -g, --schema-directory (Required) Path to the schema directory. The schema directory contains the loci FASTA
                           files and a folder named "short" that contains the FASTA files with the loci representative
                           alleles.

    -o, --output-directory (Required) Output directory where the process will store intermediate
                           files and save the final TSV file with the annotations.

    --gl, --genes-list     (Optional) Path to a file that contains the list of full paths to the loci FASTA files or
                           the loci IDs, one per line, the process should find annotations for (default: False).

    -t, --protein-table    (Optional) Path to the "cds_coordinates.tsv" file created by the CreateSchema
                           process (default: None).

    --bsr                  (Optional) BLAST Score Ratio value. This value is only used when a taxon/taxa
                           is provided and local sequences are aligned against reference proteomes
                           (default: 0.6).

    --cpu, --cpu-cores     (Optional) Number of CPU cores used to run the process (default: 1).

    --taxa                 (Optional) List of scientific names for a set of taxa. The process will
                           search for and download reference proteomes with terms that match any of
                           the provided taxa (default: None).

    --pm                   (Optional) Maximum number of proteome matches to report (default: 1).

    --no-sparql            (Optional) If provided, the process will not search for annotations 
                           through UniProt's SPARQL endpoint.

    --b, --blast-path      (Optional) Path to the BLAST executables. Use this option if chewBBACA
                           cannot find the BLASTp and makeblastdb executables or if you want to
                           use anoter BLAST installation that is not the one added to the PATH
                           (default: assumes BLAST executables were added to PATH).

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
