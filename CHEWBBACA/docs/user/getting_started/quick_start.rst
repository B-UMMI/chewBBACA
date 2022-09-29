Quick Start
===========

The following sections contain information about how to create a schema from scratch based on
a collection of genome assemblies or adapt an external schema for usage with chewBBACA, and
perform allele calling to determine the allelic profiles of a set of strains of interest.

Create a schema
:::::::::::::::

Option 1 - Genome assemblies
............................

Include all genome assemblies (complete or draft genome assemblies in FASTA format) in a directory
and adapt the following template command:

::
	
	chewBBACA.py CreateSchema -i InputAssemblies -o OutputSchemaFolder --ptf ProdigalTrainingFile

.. important::
	The CreateSchema module creates a schema seed with one representative allele per locus in the
	schema. To include more allele variants in the schema, we recommend starting by performing
	allele calling with the set of genome assemblies used for schema creation.


Option 2 - Adapt an external schema
...................................

Include all loci files (one FASTA file per locus, each file contains all alleles for a specific
locus) in a directory and adapt the following template command:

::

	chewBBACA.py PrepExternalSchema -i ExternalSchemaFastaFiles -o OutputSchemaFolder --ptf ProdigalTrainingFile

.. important::
	External schemas need to be processed to filter out sequences that do not meet a set of
	criteria applied to create every chewBBACA schema. This process might remove alleles or
	complete loci from the schemas. For more information see the page about the
	:doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module.

Perform allele calling
::::::::::::::::::::::

Determine the allelic profiles for genome assemblies:

::

	chewBBACA.py AlleleCall -i InputAssemblies -g OutputSchemaFolder/SchemaName -o OutputFolderName


Use a subset of the loci in a schema:

::

	chewBBACA.py AlleleCall -i InputAssemblies -g OutputSchemaFolder/SchemaName -o OutputFolderName --gl LociList.txt

.. important::
	The file passed to the `--gl` parameter must have one locus identifier per line (a locus
	identifier is the name of the FASTA file that contains the locus alleles).

.. important::
	- Genome assemblies and loci files from external schemas must be in FASTA format.
	- We strongly advise users to provide a Prodigal training file and to keep using the same training file to ensure consistent results (the training file used for schema creation is added to the schema's directory and automatically detected at start of allele calling without the need to pass it to the `--ptf` parameter).
	- Use the `--cpu` parameter to enable parallelization and considerably reduce execution time.
