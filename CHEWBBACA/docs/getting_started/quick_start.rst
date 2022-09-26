Quick Start
===========

Create a schema
:::::::::::::::

Option 1 - Genome assemblies

Include all genome assemblies (complete or draft genome assemblies in FASTA format) in a directory and adapt the following template command:

::
	
	chewBBACA.py CreateSchema -i InputAssemblies -o OutputSchemaFolder --ptf ProdigalTrainingFile


Option 2 - Adapt an external schema

Include all loci files (one FASTA file per locus, each file contains all alleles for a specific locus) in a directory and adapt the following template command:

::

	chewBBACA.py PrepExternalSchema -i ExternalSchemaFastaFiles -o OutputSchemaFolder --ptf ProdigalTrainingFile


Perform allele calling
::::::::::::::::::::::

Determine the allelic profiles for genome assemblies:

::

	chewBBACA.py AlleleCall -i InputAssemblies -g OutputSchemaFolder/SchemaName -o OutputFolderName


Use a subset of the loci in a schema:

::

	chewBBACA.py AlleleCall -i InputAssemblies -g OutputSchemaFolder/SchemaName -o OutputFolderName --gl LociList.txt


.. important::
	- Genome assemblies and loci files from external schemas must be in FASTA format.
	- We strongly advise users to provide a Prodigal training file and to keep using the same training file to ensure consistent results (the training file used for schema creation is added to the schema's directory and automatically detected at start of allele calling without the need to pass it to the `--ptf` parameter).
	- Use the `--cpu` parameter to enable parallelization and considerably reduce execution time.
	- The file passed to the `--gl` parameter must have one locus identifier per line (a locus identifier is the name of the FASTA file that contains the locus alleles).
