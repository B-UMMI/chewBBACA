Integration with Chewie-NS
==========================

We prepared a `chewie-NS tutorial instance <https://tutorial.chewbbaca.online/>`_ that provides a
sandbox-style environment where anyone can test functionalities. This tutorial instance enables
testing with simple cases that were specially designed to demonstrate how users can leverage the
`chewBBACA modules <https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/CHEWBBACA_NS>`_
specifically designed to interact with chewie-NS to download, upload and synchronize schemas.
The tutorial also provides clues for the exploration of data uploaded to the chewie-NS tutorial instance in its 
`website <https://tutorial.chewbbaca.online/>`_ and a better understanding of how it is
possible to interact with the API.

.. important::
  This tutorial includes sample commands that you'll have to run to upload,
  download and synchronize schemas. To connect to the correct chewie-NS instance, you
  should always check if the sample commands include ``--ns tutorial`` at the
  end and also include that parameter at the end of your commands. If you do
  not do this, the processes will try to connect to the main chewie-NS instance (it
  may download data if the identifiers you have chosen match any record but will
  not upload any data to the main instance).

The tutorial is divided into 9 steps:

- Getting started - tutorial datasets: brief description of the general structure of the tutorial
  datasets that are available.
- chewBBACA installation: links to external resources that explain how to install the chewBBACA
  suite and its dependencies.
- Uploading the tutorial schema: example on how to upload a tutorial schema to the chewie-NS 
  tutorial instance.
- Downloading the schema: download the schema uploaded to the chewie-NS tutorial instance in the 
  previous step.
- Local analysis - subset1: perform allele call to populate the downloaded schema with new alleles.
- Schema synchronization: submit novel alleles that were identified in the local analysis to
  chewie-NS to update the remote schema.
- Getting schema snapshot: download the uploaded schema in a state prior to the synchronization
  process.
- Local analysis - subset2: perform allele call with the second set of genomes.
- Schema synchronization - conflicting identifiers: synchronize local schema and remote schema to reassign local allele 
  identifiers that are conflicting with identifiers in chewie-NS and submit novel alleles 
  inferred from the second subset of genomes.


Getting started - tutorial datasets
:::::::::::::::::::::::::::::::::::

You can start by downloading any of the archives with tutorial datasets that are available
at the chewie-NS tutorial `GitHub repository <https://github.com/B-UMMI/Chewie-NS_tutorial>`_.

Currently, there are datasets for the following species:

- *Streptococcus agalactiae*
  (`sagalactiae.zip <https://github.com/B-UMMI/Chewie-NS_tutorial/blob/master/tutorial_data/sagalactiae_tutorial.zip?raw=true>`_)

In this tutorial we will provide step-by-step instructions that use the
*Streptococcus agalactiae* tutorial dataset, but the procedure is valid for any dataset that
we make available.

You will have to extract the contents in the archive. Tutorial datasets have the following
directory structure

::

    sagalactiae_tutorial
    ├── sagalactiae_genomes
    │   ├── subset1
    │   │   └── ...
    │   └── subset2
    │       └── ...
    ├── sagalactiae_schema
    │   ├── short
    │   │   ├── sagalactiae_protein1_short.fasta
    │   │   ├── ...
    │   │   └── sagalactiae_protein10_short.fasta
    │   ├── sagalactiae_protein1.fasta
    │   ├── ...
    │   ├── sagalactiae_protein10.fasta
    │   └── Streptococcus_agalactiae.trn
    ├── sagalactiae_annotations.tsv
    └── sagalactiae_description.md

The ``subset1`` and ``subset2`` directories inside the ``sagalactiae_genomes`` directory contain two
sets of genomes that will be used at different steps of the tutorial.

The ``sagalactiae_schema`` directory contains a schema with 10 loci. The ``short`` directory contained
in the schema's directory has the set of FASTA files with representative sequences for each locus in the
schema. The ``Streptococcus_agalactiae.trn`` file is the Prodigal training file used to predict coding
sequences from input genomes.

The ``sagalactiae_annotations.tsv`` file contains User and Custom annotations for the loci in the schema
and the ``sagalactiae_description.md`` file is a sample description for the schema. The Custom annotation
field in the annotations file and the description file support `markdown syntax <https://github.github.com/gfm/>`_. You may change the
contents of the files before uploading the schema if you wish.

.. important:: The sample commands provided in this tutorial include relative paths that assume that the 
               working directory is the root of the tutorial directory (the topmost level in the directory 
               structure of the tutorial dataset, ``sagalactiae_tutorial``). It is strongly advised 
               that you run the commands from that directory to ensure that you can use the commands exactly 
               as provided and obtain the same results.

chewBBACA installation
::::::::::::::::::::::

By taking advantage of chewie-NS’ API, chewBBACA is capable of handling not only the schema creation,
but also its upload, synchronization and download. The set of modules to interact with chewie-NS
included in the chewBBACA suite provide a simple and automatic solution for the main tasks
that users will want to perform.

You can install `chewBBACA <https://github.com/B-UMMI/chewBBACA>`_ through 
`conda <https://anaconda.org/bioconda/chewbbaca>`_ or `pip <https://pypi.org/project/chewBBACA/>`_.
chewBBACA has dependencies that will not be included if you install it through pip. If you install
through pip you need to ensure that you have `Prodigal <https://github.com/hyattpd/Prodigal>`_ 
and `BLAST <https://www.ncbi.nlm.nih.gov/books/NBK279671/>`_ installed and added to PATH. 
Please visit `chewBBACA's <https://github.com/B-UMMI/chewBBACA>`_ homepage for detailed 
instructions on how to install it.

.. important:: We recommend using BLAST 2.9.0. Older versions might not include functionalities
               used in the latest version of chewBBACA.


Uploading the tutorial schema
:::::::::::::::::::::::::::::

To upload schemas to the main instance of chewie-NS it is necessary to have Contributor privileges, but
in the chewie-NS tutorial instance schema upload is available to anyone that wishes to test it.
Before uploading the schema please visit the :doc:`upload_api` documentation page to learn more about the
whole process.

.. important:: The name attributed to the schema needs to be unique. You will not be able to upload
               a new schema if the schema's name has already been attributed to a schema that is
               available in chewie-NS.

To upload the schema included in the *Streptococcus agalactiae* dataset, you can run the following command 
(do not forget to set you working directory to the topmost level of the directory structure of the tutorial 
dataset and to include ``--ns tutorial`` at the end):

::

    $ chewBBACA.py LoadSchema -i sagalactiae_schema/ -sp 1 -sn tut -lp tut --df sagalactiae_description.md --a sagalactiae_annotations.tsv --ns tutorial

    ==========================
      chewBBACA - LoadSchema
    ==========================

    -- User Permissions --
    User id: 
    User role: 
    Authorized: True

    -- Parameters Validation --
    Local schema: sagalactiae_schema
    Schema's species: Streptococcus agalactiae (id=1)  ------> Species ID <------
    Number of loci: 10
    Number of alleles: 10

    Verifying schema configs...
      bsr: 0.6
      translation_table: 11
      minimum_locus_length: 201
      chewBBACA_version: 2.5.0
      size_threshold: 0.2
      word_size: None
      cluster_sim: None
      representative_filter: None
      intraCluster_filter: None
    All configurations successfully validated.

    New schema name: "tut" 
    Schema description: sagalactiae_description.md

    -- Schema Pre-processing --
    Determining data to upload...
      Loci to create and associate with species and schema: 10
      Loci without the full set of alleles: 10

    Translating sequences based on schema configs...
      Found a total of 0 invalid alleles.

    Loci missing UniProt annotation: 10
    Creating SPARQL queries to search UniProt for annotations...
    Searching for annotations on UniProt...
    Searched annotations for 10/10 loci
    User provided valid annotations for 10 loci.

    -- Schema Upload --
    Created schema with name tut (id=1).  ------> Schema ID <------

    Loci data:
      Collecting loci data...
      Sending data to the NS...
        Inserted 10 loci; Linked 10 to species; Linked 10 to schema.
      The NS completed the insertion of 10 loci.

    Alleles data:
      Collecting alleles data...
      Compressing files with alleles data...
      Sending alleles data to the NS...
        Sent data for alleles of 10 loci.

    Uploading Prodigal training file...
    Provided training file is already in the NS.

    The NS has received the data and will insert the alleles into the database.
    Schema will be available for download as soon as the process has completed.
    Schema information will also be available on the NS website.

    Removing intermediate files...


We have included the command and the information that the process prints to the standard output.
It is important to know the unique identifier that chewie-NS attributed to the schema you 
have uploaded (the lines with the schema and species identifiers are highlighted in the
standard output).
When the `LoadSchema` process finishes, chewie-NS will insert the data that was sent 
into its database and unlock the schema to make it available for download. You can find
the schema you have uploaded listed in the ``Schemas Overview`` page for the species 
(`Schemas Overview page for *Streptococcus agalactiae* <https://tutorial.chewbbaca.online/species/1>`_).

.. important:: Schemas that are uploaded to the chewie-NS tutorial instance are deleted after 48h.

Downloading the schema
::::::::::::::::::::::

In order to use a schema you have uploaded to chewie-NS, you will have to download it.

To know more about the ``DownloadSchema`` process, please visit the :doc:`</user/modules/DownloadSchema>` page
in the documentation.

To download the schema you have uploaded, please run the following command:

.. important:: Substitute the species and schema ID values, ``-sp`` and ``-sc``, by the values that 
               serve to identify the schema you have uploaded.

::

    $ chewBBACA.py DownloadSchema -sp 1 -sc 1 -o sagalactiae_ns --ns tutorial

    ==============================
      chewBBACA - DownloadSchema
    ==============================

    Schema id: 1
    Schema name: tut
    Schema's species: Streptococcus agalactiae (id=1)

    Downloading compressed version...
    Decompressing schema...
    Schema is now available at: sagalactiae_ns/sagalactiae_tut

The process will download a ready-to-use schema to the output directory you have specified.
The loci and alleles included in the schema are the same that were in the original schema,
but chewie-NS has attributed new identifiers that will help to unmistakably identify
those loci and alleles and facilitate results comparison for anyone that is using the same
schema.

Local analysis - subset1
:::::::::::::::::::::::::::

You can use the schema you have downloaded to perform allele call and determine the allelic
profiles of a set of genomes. Allele calling is performed locally and privately, without the
need to provide any data or private information. You can learn more about the ``AlleleCall``
process in its :doc:`page </user/modules/AlleleCall>`.

If you open any FASTA file in the schema that you have downloaded, you will find sequences
that have the following header structure:

::

    $ cat tut-00000001.fasta

    >tut-00000001_1
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...

Headers start with the loci prefix (``tut``) followed by the loci integer identifier (``00000001``)
and end with the allele identifier (``1``).

To perform allele call and determine the allelic profiles of the genomes in the subset1, run
the following command:

::

    $ chewBBACA.py AlleleCall -i sagalactiae_genomes/subset1/ -g sagalactiae_ns/sagalactiae_tut/ -o subset1_results 

    ==========================
      chewBBACA - AlleleCall
    ==========================

    Prodigal training file: Streptococcus_agalactiae.trn
    Number of CPU cores: 1

    Checking dependencies...
    Blast installation...True
    Prodigal installation...True
    Blast version meets minimum requirements (>=2.5.0).

    Checking if genome files exist...
    Checking if gene files exist...

    Starting Prodigal at: ...
    done prodigal run on:GCA_000012705.1_ASM1270v1_genomic.fna
    done prodigal run on:GCA_000007265.1_ASM726v1_genomic.fna
    done prodigal run on:GCA_000302475.2_ASM30247v2_genomic.fna
    done prodigal run on:GCA_000196055.1_ASM19605v1_genomic.fna
    done prodigal run on:GCA_000299135.1_ASM29913v1_genomic.fna
    done prodigal run on:GCA_000427035.1_09mas018883_genomic.fna
    done prodigal run on:GCA_000427055.1_ILRI112_genomic.fna
    done prodigal run on:GCA_000427075.1_ILRI005_genomic.fna
    done prodigal run on:GCA_000599965.1_ASM59996v1_genomic.fna
    done prodigal run on:GCA_000689235.1_GBCO_p1_genomic.fna
    done prodigal run on:GCA_000730255.1_ASM73025v1_genomic.fna
    done prodigal run on:GCA_000730215.2_ASM73021v2_genomic.fna
    Finishing Prodigal at: ...

    Checking if Prodigal created all the necessary files...
    All files were created.

    Translating genomes...
    Creating Blast databases for all genomes...

    Starting Allele Calling at: ...
    Processing tut-00000002.fasta. Start ... Locus 9 of 10. Done 90%.
    Finished Allele Calling at: ...

    Wrapping up the results...
    ##################################################
    12 genomes used for 10 loci

    Used a BSR of: 0.6

    17 exact matches found out of 120

    14.17 percent of exact matches
    ##################################################

    Writing output files...

    ------------------------------------------------------------------------------------------
    Genome                                      EXC    INF    LNF   PLOT   NIPH    ALM    ASM 
    ------------------------------------------------------------------------------------------
    GCA_000007265.1_ASM726v1_genomic.fna         1      5      4      0      0      0      0  
    GCA_000012705.1_ASM1270v1_genomic.fna        1      4      5      0      0      0      0  
    GCA_000196055.1_ASM19605v1_genomic.fna       1      5      4      0      0      0      0  
    GCA_000299135.1_ASM29913v1_genomic.fna       4      1      4      0      0      0      1  
    GCA_000302475.2_ASM30247v2_genomic.fna       0      5      5      0      0      0      0  
    GCA_000427035.1_09mas018883_genomic.fna      2      3      5      0      0      0      0  
    GCA_000427055.1_ILRI112_genomic.fna          1      4      4      0      0      0      1  
    GCA_000427075.1_ILRI005_genomic.fna          1      5      4      0      0      0      0  
    GCA_000599965.1_ASM59996v1_genomic.fna       0      5      5      0      0      0      0  
    GCA_000689235.1_GBCO_p1_genomic.fna          0      5      5      0      0      0      0  
    GCA_000730215.2_ASM73021v2_genomic.fna       3      3      4      0      0      0      0  
    GCA_000730255.1_ASM73025v1_genomic.fna       3      2      4      0      0      0      1  
    ------------------------------------------------------------------------------------------

    Checking the existence of paralog genes...
    Detected number of paralog loci: 0

    Creating SQLite database to store profiles...done.
    Inserted 10 loci into database.

    Sending allelic profiles to SQLite database...done.
    Inserted 12 profiles (12 total, 12 total unique).

The ``AlleleCall`` process will print a table with the summary of the results to the standard
output. You can see a  :doc:`detailed description </user/modules/AlleleCall>` 
of each category but for the purpose of this tutorial, the ``INF`` cases are the most relevant. The alleles
that received this classification correspond to new alleles that have been inferred during the 
process and were added to the schema FASTA files. If we inspect the same file that we looked into
before the allele calling, you will notice that new alleles have been added to that file.

::

    $ cat tut-00000001.fasta

    >tut-00000001_1
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000007265.1-ASM726v1-genomic.fna_07/08/2020T21:00:08_*2
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000012705.1-ASM1270v1-genomic.fna_07/08/2020T21:00:08_*3
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000196055.1-ASM19605v1-genomic.fna_07/08/2020T21:00:08_*4
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000302475.2-ASM30247v2-genomic.fna_07/08/2020T21:00:08_*5
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000427055.1-ILRI112-genomic.fna_07/08/2020T21:00:08_*6
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000427075.1-ILRI005-genomic.fna_07/08/2020T21:00:08_*7
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000599965.1-ASM59996v1-genomic.fna_07/08/2020T21:00:08_*8
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_S_GCA-000689235.1-GBCO-p1-genomic.fna_07/08/2020T21:00:08_*9
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...

New alleles added to loci files that belong to a schema that was downloaded from chewie-NS will
include a ``*`` before the allele number in the end of the sequence identifier (e.g.: ``*4``). The ``*`` serves to indicate that
the alleles were identified locally and that it has not been verified if those alleles exist in
chewie-NS and, if they exist, what was the identifier that chewie-NS attributed.

Schema synchronization
::::::::::::::::::::::

To verify if newly identified alleles exist in chewie-NS, and submit those alleles if they are
not in chewie-NS, we will need to run the ``SyncSchema`` process. This process will retrieve
alleles added to the remote schema in chewie-NS since the last time we synchronized the local
and remote schemas and offer the option to submit novel alleles that have been identified in
local analyses and are not in chewie-NS. To learn more about the ``SyncSchema`` process, please
visit the :doc:`</user/modules/SyncSchema>` page.

Running the ``SyncSchema`` process is fairly simple. To retrieve new alleles added to the remote
schema since the last synchronization process, we only need to provide the path to the directory
with the schema files. We also want to submit any novel alleles that our local schema might have,
so we include the ``--submit`` argument (there is no need to include ``--ns tutorial`` because
the ``SyncSchema`` process automatically detects what is the chewie-NS instance the schema was 
downloaded from).

::

    $ chewBBACA.py SyncSchema -sc sagalactiae_ns/sagalactiae_tut/ --submit

    ==========================
      chewBBACA - SyncSchema
    ==========================

    Schema id: 1
    Schema name: tut
    Schema's species: Streptococcus agalactiae (id=1)
    Last synced: 2020-08-07T22:46:52.406869

    Remote schema was last modified on: 2020-08-07T22:46:52.406869

    Retrieving alleles added to remote schema after 2020-08-07T22:46:52.406869...
    Retrieved 0 alleles for 0 loci.
    Local schema has 47 novel alleles for 7 loci.
    Collecting data and creating files to submit local alleles...
    Sending and inserting new alleles...
        Sent data for alleles of 7/7 loci.
        Inserted 47 alleles.
    The Chewie-NS inserted 47 new alleles and detected 0 repeated alleles.

    Updating local allele identifiers...
    Updated 12 profiles.

    Adapting schema in the following directory:
    sagalactiae_ns/temp
    Prodigal training file:

    Number of cores: 1
    BLAST Score Ratio: 0.6
    Translation table: 11
    Minimum accepted sequence length: 201
    Size threshold: None
    Number of genes to adapt: 7

    Determining the total number of alleles and allele mean length per gene...

    Adapting 7 genes...

    [==========] 100%

    Number of invalid genes: 0
    Number of invalid alleles: 0

    Successfully adapted 7/7 genes present in the external schema.

    Received 0 new alleles for 7 loci and sent 47 for 7 loci. 

Since the schema has not been modified since the upload date, the synchronization process 
will not retrieve alleles from chewie-NS. Our local schema includes alleles that are not in chewie-NS
and the synchronization process will send those alleles to chewie-NS, waiting for the insertion 
process to finish and return the set of identifiers that were attributed to the novel alleles.
The ``SyncSchema`` process will reassign allele identifiers to local alleles based on the 
identifiers attributed by chewie-NS and re-determine representative sequences for the loci
that were altered. The schema had not been altered since its upload and chewie-NS attributed
the same allele identifiers that were already being used in the local schema. Thus, the sequence
headers will be shortened and the synchronization process will simply remove the ``*`` from the 
headers. The file structure will be changed to the following:

::

    $ cat tut-00000001.fasta

    >tut-00000001_1
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_2 <----- *2
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_3 <----- *3
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_4 <----- *4
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_5 <----- *5
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_6 <----- *6
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_7 <----- *7
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_8 <----- *8
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_9 <----- *9
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...

We have included the mapping between the new identifiers and the old identifiers with ``*`` to 
highlight allele identifiers reassignments, e.g.: ``>tut-00000001_2 <----- *2``. This mapping 
serves to clearly indicate the changes made by the SyncShcema process during this tutorial and 
is not added to the FASTA files.

Getting schema snapshot
:::::::::::::::::::::::

To demonstrate a synchronization process that will need to perform more complicated reassignments
to ensure that local and remote schemas share the same identifiers, we will start by using a 
feature that allows users to download a snapshot of any schema. Quickly consult the ``Schemas Overview``
table and copy the ``Last Change Date`` of the schema that you have uploaded. We will subtract 2 minutes 
from that date and slightly modify the date format so that it matches the ``yyyy-mm-ddThh:mm:ss`` format 
(if the ``Last Change Date`` is ``2020-08-07T22:49:52``, the date that you should include in the command 
is ``2020-08-07T22:47:52`` as indicated in the example shown).

.. important:: Substitute the data and time below with the time you have calculated, do NOT simply copy 
               the command!

A sample command would be:

::

    $ chewBBACA.py DownloadSchema -sp 1 -sc 1 -o sagalactiae_snapshot --ns tutorial --d "2020-08-07T22:47:52"

    ==============================
      chewBBACA - DownloadSchema
    ==============================

    Schema id: 1
    Schema name: tut
    Schema's species: Streptococcus agalactiae (id=1)

    Downloading schema FASTA files...
    Number of loci to download: 10
    Downloading schema files...
    Downloaded: 10/10
    Downloaded and wrote FASTA files for 10/10 loci
    Failed download for 0 loci.

    Adapting schema in the following directory:
    sagalactiae_snapshot
    Prodigal training file:
    sagalactiae_snapshot/Streptococcus_agalactiae.trn
    Number of cores: 1
    BLAST Score Ratio: 0.6
    Translation table: 11
    Minimum accepted sequence length: 201
    Size threshold: None
    Number of genes to adapt: 10

    Determining the total number of alleles and allele mean length per gene...

    Adapting 10 genes...

    [==========] 100%

    Number of invalid genes: 0
    Number of invalid alleles: 0

    Successfully adapted 10/10 genes present in the external schema.

    Schema is now available at: sagalactiae_snapshot/sagalactiae_tut

This will download all FASTA files for all loci in the schema and construct the schema locally.
Since we have requested for the schema in a state prior to its ``Last Change Date``, we will
retrieve a schema that does not include all alleles in the latest version of the remote schema
and is outdated.

Local analysis - subset2
:::::::::::::::::::::::::::

We will perform allele call with the genomes in subset2 to demonstrate how the ``SyncSchema``
process would behave if the remote schema had already been modified by another user and the
sequences and allele identifiers in our local schema and in the remote schema did not fully
match.

::

    $ chewBBACA.py AlleleCall -i sagalactiae_genomes/subset2/ -g sagalactiae_snapshot/sagalactiae_tut/ -o subset2_results 

    ...

    ------------------------------------------------------------------------------------------
    Genome                                      EXC    INF    LNF   PLOT   NIPH    ALM    ASM 
    ------------------------------------------------------------------------------------------
    GCA_000782855.1_ASM78285v1_genomic.fna       0      6      4      0      0      0      0  
    GCA_000831105.1_ASM83110v1_genomic.fna       1      6      3      0      0      0      0  
    GCA_000831125.1_ASM83112v1_genomic.fna       4      3      3      0      0      0      0  
    GCA_000831145.1_ASM83114v1_genomic.fna       5      2      3      0      0      0      0  
    GCA_001026925.1_ASM102692v1_genomic.fna      0      6      3      0      0      0      1  
    GCA_001190865.1_ASM119086v1_genomic.fna      0      5      5      0      0      0      0  
    GCA_001190885.1_ASM119088v1_genomic.fna      1      5      4      0      0      0      0  
    GCA_001266635.1_ASM126663v1_genomic.fna      1      5      4      0      0      0      0  
    GCA_001275545.2_ASM127554v2_genomic.fna      3      3      4      0      0      0      0  
    GCA_001448985.1_ASM144898v1_genomic.fna      5      2      3      0      0      0      0  
    GCA_001655175.1_ASM165517v1_genomic.fna      0      5      5      0      0      0      0  
    GCA_001683515.1_ASM168351v1_genomic.fna      5      1      4      0      0      0      0  
    ------------------------------------------------------------------------------------------

    ...

Once again, we verify that the ``AlleleCall`` process inferred some alleles during its execution
and that those alleles have been added to the local schema. Since we have used a different set of
genomes we do not know if the set of alleles that were added to the schema are in the remote schema,
nor if the alleles that are common to both schemas have been attributed the same identifiers (in this
case they have not and it is very unlikely that different sets of genomes will lead to the same results
and schema modifications).

Schema synchronization - conflicting identifiers
::::::::::::::::::::::::::::::::::::::::::::::::

In the final step we will synchronize our schema with the remote schema. This process will retrieve
alleles that are in the remote schema and add them to our schema with the identifier they have in
chewie-NS. The alleles that are not in chewie-NS will be shifted to the end of the FASTA files and
assigned sequential identifiers with ``*`` at the end and in the same order as they were added to the 
schema. This ensures that there are no conflicts between remote and strictly local identifiers. Local 
alleles with ``*`` in their identifiers will be sent to chewie-NS and inserted into the schema's database. 
The ``SyncSchema`` process wil receive the identifiers attributed by chewie-NS and assign them to 
the local sequences that still had no global identifier, ensuring that all alleles have the correct 
identifier and that there is a common and global nomenclature.

To perform this last synchronization, execute:

::

    $ chewBBACA.py SyncSchema -sc sagalactiae_snapshot/sagalactiae_tut/ --submit

    ...

    Received 47 new alleles for 7 loci and sent 33 for 7 loci.

    ...

The synchronization process will retrieve 47 alleles that were inferred from subset1
and send 33 local alleles that were inferred from subset2. Identifier reassignmnent results
in the following file structure:

::

    $ cat tut-00000001.fasta

    >tut-00000001_1
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_2 <----- *7
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_3 <----- *5
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_4
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_5 <----- *9
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_6 <----- *8
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_7
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_8 <----- *6
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_9
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_10 <----- *2
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_11 <----- *3
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...
    >tut-00000001_12 <----- *4
    ATGTTTAAAGGTAATAAGAAGTTGAATAGTTCTAAATTAGGTGATTACACACCACTTGAATTTGGTTCT...

-------------------------------------------------------------------------------------------

Reading the documentation and completing the tutorial should provide a good overview of
how chewie-NS works and how you can interact with it through the chewBBACA suite. You can 
head to `chewie-NS' main instance website <https://chewbbaca.online/>`_ to explore available schema data for several species
and download data through the website or using the chewBBACA modules that were used in the 
tutorial. Schema upload and allele submission during synchronization in chewie-NS' main instance
are only possible to authorized users. If you want to submit data or provide
any type of feedback, please contact us through imm-bioinfo@medicina.ulisboa.pt.
