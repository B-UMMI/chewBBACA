SyncSchema - Synchronize a schema with its remote version in Chewie-NS
======================================================================

The *SyncSchema* module allows users to synchronize local schemas, previously downloaded from
Chewie-NS, with their remote versions. All chewBBACA users can synchronize schemas to **get the
latest alleles added to Chewie-NS and to ensure that a common allele identifier nomenclature is
maintained for the alleles that are common between the local instance and the remote schemas**.
We also provide the option to **submit novel alleles**, that were identified locally and are not
present in Chewie-NS.

.. important::
    Only authorized users can submit new local alleles to update remote schemas, although all users
    can download schemas and novel alleles from the Chewie-NS public server. Please send a request to
    imm-bioinfo@medicina.ulisboa.pt if you wish to submit novel alleles.

To synchronize a local schema with its remote Chewie-NS public server version it is only necessary to provide the path
to the schema directory. The simplicity of the process is ensured by a configuration file,
present in all schemas downloaded from Chewie-NS, that contains the identifier of the
schema in Chewie-NS and the last modification date of the schema.

Configuration file content

::

    ['2020-06-30T19:10:37.466104', 'http://chewbbaca.online/NS/api/species/9/schemas/1']

Novel alleles identified locally are added to the schema with a ``*`` preceding their integer
identifier to indicate these are temporary designations. In the example below alleles 4 to 7
were detected locally but were not present in the database that had been retrieved from the
Chewie-NS public server.

Local FASTA file example

::

    >prefix-018550_1
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA
    >prefix-018550_2
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA
    >prefix-018550_3
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA
    >prefix-018550_*4
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA
    >prefix-018550_*5
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA
    >prefix-018550_*6
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA
    >prefix-018550_*7
    ATGAGCAAGCCTAATGTTGTTCAGTTAAATAATCAATATATTAACGATGAGAATCTAAAAAAACGTTACGAAGCTGAGGAGTTACGCTAA

On synchronizing, novel alleles retrieved from Chewie-NS are compared to novel local alleles and
the process reassigns allele identifiers to ensure that the alleles common to local and remote
schemas have the same identifiers. Local alleles that are not in Chewie-NS are shifted to the last
positions in the FASTA files and keep a ``*`` in the identifier. If the user wants to submit those
alleles (only available to authorized users) it should add the ``--submit`` flag and the necessary
data will be collected and uploaded to Chewie-NS. Chewie-NS will return the identifiers assigned to
the submitted alleles and the local process will remove the '*' from the submitted alleles and
assign the permanent identifiers. If the process retrieves new alleles from Chewie-NS, it will
redetermine representative sequences **ONLY** for the loci in the local schema that were altered
by the synchronization process.

.. important::
    It is strongly advised that users adjust the value of the ``--cpu`` argument in order to
    accelerate the determination of representative sequences.

Basic Usage
-----------

If we want to synchronize a local schema we only need to provide the path to the directory that
contains the schema:

::

    $ chewBBACA.py SyncSchema -i path/to/SchemaFolder

The ``--submit`` argument allows users to submit novel alleles in their local schema to the
Chewie-NS public server:

::

    $ chewBBACA.py SyncSchema -i path/to/SchemaFolder --submit

Parameters
----------

::

    -sc, --schema-directory     (Required) Path to the directory with the schema to be synced (default: None).

    --cpu, --cpu-cores          (Optional) Number of CPU cores/threads that will be used to run the process
                                (chewie resets to a lower value if it is equal to or exceeds the total
                                number of available CPU cores/threads). This value is only used if the
                                process retrieves novel alleles from the remote schema and needs to
                                redetermine the set of representative alleles for the local schema (default: 1).

    --ns, --nomenclature-server (Optional) The base URL for the Chewie-NS instance. The default
                                option will get the base URL from the schema's URI. It is also
                                possible to specify other options that are available in chewBBACA's
                                configs, such as: "main" will establish a connection to "https://chewbbaca.online/",
                                "tutorial" to "https://tutorial.chewbbaca.online/" and "local" to
                                "http://127.0.0.1:5000/NS/api/" (localhost). Users may also provide
                                the IP address to other Chewie-NS instances (default: None).

    --b, --blast-path           Path to the directory that contains the BLAST executables (default:None).
                                                   
    --submit                    If the process should identify new alleles in the local schema and
                                send them to the Chewie-NS instance. (only authorized users can submit new alleles)
                                (default: False).

Workflow of the LoadSchema module
:::::::::::::::::::::::::::::::::

.. image:: /_static/images/SyncSchema.png
   :width: 1200px
   :align: center
