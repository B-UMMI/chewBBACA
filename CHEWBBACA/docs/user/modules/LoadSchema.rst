LoadSchema - Upload a schema to Chewie-NS
=========================================

The *LoadSchema* module enables the upload of schemas to Chewie-NS.

.. important::
    **You need to be a registered user in Chewie-NS with Contributor privileges to be able to upload schemas!**

    If you have registered at the `chewie-NS public website <https://chewbbaca.online/auth>`_
    and want to contribute new schemas please contact us via e-mail at: imm-bioinfo@medicina.ulisboa.pt

To upload a schema to Chewie-NS it is required to provide:

- The **path** to the local schema.

  - The schema must have been created with chewBBACA v2.5.0 or a later version. If your schema was created with
    an older version, please adapt the schema with the ``PrepExternalSchema`` process or run the 
    ``AlleleCall`` process to convert the schema to the latest version.

.. warning:: **Only schemas that have been used with the same valid
             value per parameter can be uploaded (this restriction applies
             to the BLAST Score Ratio, Prodigal training file, minimum 
             sequence length, genetic code and sequence size variation 
             threshold parameters).**
             
             Invalid or multiple values used in different ``AlleleCall`` runs
             for a single parameter can lead to inconsistent results; thus,
             it is strongly advised to always perform allele calling with
             the same set of parameters and refrain from altering the initial
             set of parameter values defined in the schema creation or
             adaptation processes.

- The **species ID** or **scientific name** of the species that the schema will be associated to.
  
  - There are 3 ways to know the **species ID** of a species: 1) you can consult the `Overview <https://chewbbaca.online/stats>`_ 
    table in the Chewie-NS website; 2) you can use the 
    `NSStats <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/CHEWBBACA_NS/stats_requests.py>`_ 
    process in the  chewBBACA suite to directly obtain information about the species and schemas in Chewie-NS or; 3) you can 
    query the ``/species/list`` API endpoint through  `Swagger <https://chewbbaca.online/api/NS/api/docs>`_ or a simple curl 
    command (``e.g.: curl -X GET "https://chewbbaca.online/NS/api/species/list" 
    -H  "accept: application/json"``).
  - e.g.: ``9`` or ``Escherichia coli``.

- A **name** for the schema.

  - The name should be short and concise. The name **must be unique** among the set of names for 
    schemas of the same species (this means that the using the same name of an existing schema will lead to an error)
    and should not include spaces.
  - e.g.: ``Project_cgMLST``, ``SRA_wgMLST``, ``Organization_cgMLST`` .

- A **prefix** for the loci identifiers to facilitate the identification of the schema they belong to.

  - You may use the name of the schema as prefix to ensure prefix uniqueness for the loci
    of a schema.

Users may provide a description about the schema. The file with the description 
will be sent to the Chewie-NS and displayed in the schema's page in the Chewie-NS website. Markdown syntax is 
supported in order to allow greater customizability of the rendered description.
For more information on the Markdown specification accepted by Chewie-NS please visit the
`Github Flavored Markdown Specification page <https://github.github.com/gfm/>`_.

Sample description

::

    # Whole-genome MLST schema for *Species name*

    This schema was created with [chewBBACA 2.5.0](https://github.com/B-UMMI/chewBBACA).

    ## Schema creation and validation

    A total of 100 *Species name* genomes were used to create this wgMLST schema.
    (Add more information that might be relevant to reproduce the process.
    For instance the assembly pipeline used or links to relevant external data repositories)

    ## Dataset

    All raw reads from SRA annotated as *Species name* were assembled into 100 genomes.
    (Add more information about dataset creation/collection. Include date of data download.
    A list of accession numbers or links to relevant external data repositories maybe useful)

    ## Citations

    (Add any relevant citations)

    For more information please access [external page](https://external/page)
    (If there is any external source with more information, link it here)

The process queries UniProt's SPARQL endpoint to retrieve annotations for the loci 
in the schema. The user that uploads the schema can provide a TSV file with annotations for some or all 
loci in the schema. The file with annotations must have the following structure:

- First column: locus identifier (name of locus file without ``.fasta`` extension).
- Second column: user annotation (name commonly attributed by the user).
- Third column: custom annotation (another term that the user might want to attribute).

However, no headers are necessary.

.. rst-class:: align-center

  +----------+--------------+---------------------------------------------------+
  | locus_1  |     dnaA     |  Chromosomal replication  initiator protein DnaA  |
  +----------+--------------+---------------------------------------------------+
  | locus_2  |     dnaG     |                    DNA primase                    |
  +----------+--------------+---------------------------------------------------+
  | locus_3  |              |            RNA-directed DNA polymerase            |
  +----------+--------------+---------------------------------------------------+
  | locus_4  |     pbp      |                                                   |
  +----------+--------------+---------------------------------------------------+

It is not necessary to provide both annotation types for each locus nor for every locus.
If no information is provided N/A will be automatically shown in the locus details page in
Chewie-NS.

Basic Usage
-----------

To upload a schema for *Escherichia coli*, we could run one of the following commands:

- Providing the species ID:

::

	$ chewBBACA.py LoadSchema -i path/to/schema/to/be/sent -sp 9 -sn cgMLST_95 -lp cgMLST_95

- Providing the species name:

::

	$ chewBBACA.py LoadSchema -i path/to/schema/to/be/sent -sp "Escherichia coli" -sn cgMLST_95 -lp cgMLST_95

To upload a schema and provide a description and annotations:

::

    $ chewBBACA.py LoadSchema -i path/to/schema/to/be/sent -sp 9 -sn cgMLST_95 -lp cgMLST_95 --df description.txt --a annotations.tsv

To continue an upload that was interrupted or that aborted, we should provide the command used in 
the process that failed and add the ``--continue_up`` argument

::
	
    $ chewBBACA.py LoadSchema -i path/to/schema/to/be/sent -sp 9 -sn cgMLST_95 -lp cgMLST_95 --continue_up

.. important:: **If you cannot complete schema upload or if the information in the
                 website is incorrect or missing, please contact us via e-mail:**
                 imm-bioinfo@medicina.ulisboa.pt

Parameters
----------

::

    -i, --schema-directory      (Required) Path to the directory of the schema to upload (default: None).

    -sp, --species-id           (Required) The integer identifier or name of the species that the
                                schema will be associated to in the NS (default: None).

    -sn, --schema-name          (Required) A brief and meaningful name that should help understand
                                the type and content of the schema (default: None).

    -lp, --loci-prefix          (Required) Prefix included in the name of each locus of the schema
                                (default: None).

    --df, --description-file    (Optional) Path to a text file with a description about the schema.
                                Markdown syntax is supported in order to offer greater customizability
                                of the rendered description in the Frontend. Will default to the schema's
                                name if the user does not provide a valid path for a file (default: None).

    --a, --annotations          (Optional) Path to a TSV file with loci annotations. The first column
                                has loci identifiers (w/o .fasta extension), the second has user
                                annotations and the third has custom annotations (default: None).

    --cpu, --cpu-cores          (Optional) Number of CPU cores that will be used in the Schema
                                Pre-processing step (default: 1).

    --ns, --nomenclature-server (Optional) The base URL for the Nomenclature Server. The default value,
                                "main", will establish a connection to "https://chewbbaca.online/",
                                "tutorial" to "https://tutorial.chewbbaca.online/" and "local" to
                                "http://127.0.0.1:5000/NS/api/" (localhost). Users may also provide
                                the IP address to other Chewie-NS instances (default: main).

    --continue_up               (Optional) If the process should check if the schema upload was
                                interrupted and try to finish it (default: False).
