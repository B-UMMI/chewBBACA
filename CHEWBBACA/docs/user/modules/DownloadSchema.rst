DownloadSchema -  Download a schema from Chewie-NS
==================================================

Compressed versions of every schema in the Chewie-NS are available for download through the 
`Chewie-NS public website <https://chewbbaca.online/>`_ or through the ``/species/{species_id}/schemas/{schema_id}/zip``
API endpoint in `Swagger <https://chewbbaca.online/api/NS/api/docs>`_ or with a simple curl command
(``e.g.: curl -X GET "https://chewbbaca.online/NS/api/species/9/schemas/1/zip?request_type=download" -H  "accept: application/json"``).
You can also take advantage of the integration with the `chewBBACA suite <https://github.com/B-UMMI/chewBBACA>`_ and use the 
`download_schema.py <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/CHEWBBACA_NS/down_schema.py>`_ script.

.. note:: Compressed versions are ZIP archives that contain ready-to-use schemas. Simply extract
          and you can start performing allele calls using the `chewBBACA suite <https://github.com/B-UMMI/chewBBACA>`_.

.. important:: Chewie-NS generates new compressed versions of each schema every 24h, if the
               schemas were updated since the last compression date. This means that the compressed
               version is not always the latest. If that is the case, the integration with
               chewBBACA allows to quickly update your local version with  the latest information
               by using the process described at :docs:`synchronize_api`.

To download a schema with the chewBBACA suite (``chewBBACA.py DownloadSchema``), it is necessary to provide:

- The **ID of the species** that the schema is associated with and the **ID of the schema**.

  - To know the **ID of a species** you can consult the `Overview <https://chewbbaca.online/stats>`_ 
    table in the Chewie-NS public website, query the ``/species/list`` API endpoint through  
    `Swagger <https://chewbbaca.online/api/NS/api/docs>`_, or a simple curl command 
    (``e.g.: curl -X GET "https://chewbbaca.online/NS/api/species/list" -H  "accept: application/json"``).
    To know the **ID of the schema** you want to download, you can click the ``SCHEMA DETAILS`` button 
    in the `Overview <https://chewbbaca.online/stats>`_ table to get a list with all the schemas and their 
    **IDs** for a given species, query the ``/species/{species_id}`` API endpoint through  
    `Swagger <https://chewbbaca.online/api/NS/api/docs>`_, or a simple curl command 
    (``e.g.: curl -X GET "https://chewbbaca.online/NS/api/species/1" -H  "accept: application/json"``).
    Alternatively, you can use the `NSStats <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/CHEWBBACA_NS/stats_requests.py>`_ 
    process in the  chewBBACA suite to get information about species and schemas in the Chewie-NS
  - e.g.: species ID = ``9`` and schema ID = ``1``.

- **Path to the output directory** that will store the schema.

  - If the directory does not exist, the process will create it (will not create
    parent directories that do not exist). If the directory exists, it **must be empty**
    or the process will exit without downloading the schema.

The chewBBACA suite provides the option to download a schema snapshot
at a given date. The date should be in the format ``yyyy-mm-ddThh:mm:ss`` (e.g.: ``2020-06-30T19:10:37``).
It also allows users to request the latest version of a schema (``--latest``), if the compressed version that
is available compressed is outdated. An alternative and more efficient approach that can be applied to get the latest 
version of the schema is to download the compressed version available and run the 
`SyncSchema <https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/CHEWBBACA_NS/sync_schema.py>`_ 
process to retrieve the alleles that were added to the schema after the creation of the compressed file.

.. note:: The DownloadSchema process will download the compressed version that is available
          by default. If the provided date matches the date of the latest compressed version,
          it will download the compressed version, otherwise it will download the FASTA files
          and construct the schema locally.

.. important:: It is strongly advised that users adjust the value of the ``--cpu`` argument
               if they antecipate that the process will have to construct the schema locally.
               Schema adaptation is relatively fast but will greatly benefit if it can distribute
               work to several CPU cores.

Basic Usage
-----------

To download a schema of *Escherichia coli* we need to provide the ID of the species and the ID of the schema that we want to download

::

    $ chewBBACA.py DownloadSchema -sp 9 -sc 1 -o path/to/download/folder

To download a snapshot of the schema at a given date

::

    $ chewBBACA.py DownloadSchema -sp 9 -sc 1 -o path/to/download/folder --date 2020-06-30T19:10:37

To retrieve the latest version of the schema

::

    $ chewBBACA.py DownloadSchema -sp 9 -sc 1 -o path/to/download/folder --latest 

Parameters
----------

::

    $ chewBBACA.py DownloadSchema -h

    ...
                                
      -sp SPECIES_ID            The integer identifier or name of the species that
                                the schema is associated to in the NS. (default:
                                None)
                                
      -sc SCHEMA_ID             The URI, integer identifier or name of the schema
                                to download from the NS. (default: None)
                                
      -o DOWNLOAD_FOLDER        Output folder to which the schema will be saved.
                                (default: None)
                                
      --cpu CPU_CORES           Number of CPU cores that will be passed to the
                                PrepExternalSchema process to determine
                                representatives and create the final schema.
                                (default: 1)
                                
      --ns NOMENCLATURE_SERVER  The base URL for the Nomenclature Server. The
                                default value, "main", will establish a connection
                                to "https://chewbbaca.online/", "tutorial" to
                                "https://tutorial.chewbbaca.online/" and "local"
                                to "http://127.0.0.1:5000/NS/api/" (localhost).
                                Users may also provide the IP address to other
                                Chewie-NS instances. (default: main)
                                
      --d DATE                  Download schema with state from specified date.
                                Must be in the format "Y-m-dTH:M:S". (default:
                                None)
                                
      --latest                  If the compressed version that is available is not
                                the latest, downloads all loci and constructs
                                schema locally. (default: False)
