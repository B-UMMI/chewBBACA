NSStats -  Retrieve basic information about the species and schemas in Chewie-NS
================================================================================

The *NSStats* module enables the retrieval of information from the Chewie-NS server. Its main
objective is to provide information about the list of species and schemas in Chewie-NS, so that
users can quickly identify a schema of interest and download it.

Basic Usage
-----------

Retrieve the list of species and the total number of schemas, loci and alleles per species:

::

    $ chewBBACA.py NSStats -m species

    ------------------------------------------------------------------------------
    Species                             id       #schemas     #loci      #alleles 
    ------------------------------------------------------------------------------
    Mycobacterium tuberculosis          1           1          2891       748546  
    Clostridioides difficile            2           1          2147       44556   
    Salmonella enterica                 4           2         11558      7556832  
    Acinetobacter baumannii             5           1          2390       701466  
    Arcobacter butzleri                 7           4         10592       29956   
    Campylobacter jejuni                8           1          2794       298415  
    Escherichia coli                    9           1          7601      2958460  
    Yersinia enterocolitica             10          1          6344       126383  
    ------------------------------------------------------------------------------

Retrieve the list of schemas for a species and the total number of loci and alleles per schema:

::

    $ chewBBACA.py NSStats -m schemas --sp 4

    Salmonella enterica (id=4)
    ------------------------------------------------------------------
    Schema_name                         id        #loci      #alleles 
    ------------------------------------------------------------------
    cgMLSTEnterobase                    1          3000      4726816  
    INNUENDO_cgMLST                     2          8558      2830016  
    ------------------------------------------------------------------

Retrieve property values for a schema:

::

    $ chewBBACA.py NSStats -m schemas --sp 4 --sc 2

    -------------------------------------
    Salmonella enterica - INNUENDO_cgMLST
    -------------------------------------

    ID: 2
    Created by: chewie
    Total loci: 8558
    Total alleles: 2830016
    BLAST Score Ratio: 0.6
    chewBBACA version: 2.1.0
    Genetic code: 11
    Minimum sequence length: 0
    Sequence length variation threshold: None
    Clustering word size: None
    Clustering similarity: None
    Representative similarity filter: None
    Intracluster similarity filter: None
    Creation date: 2020-06-04T23:23:24.099616
    Last modified: 2020-06-04T23:23:24.099616

Parameters
----------

::

    -m, --mode                  (Required) The process can retrieve the list of species ("species"
                                option) in the Chewie-NS or the list of schemas for a species
                                ("schemas" option).

    --sp, --species-id          (Optional) The integer identifier of a species in Chewie-NS
                                (default: None).

    --sc, --schema-id           (Optional) The integer identifier of a schema in Chewie- NS
                                (default: None).

    --ns, --nomenclature-server (Optional) The base URL for the Chewie-NS instance. The default
                                value, "main", will establish a connection to "https://chewbbaca.online/",
                                "tutorial" to "https://tutorial.chewbbaca.online/" and "local" to
                                "http://127.0.0.1:5000/NS/api/" (localhost). Users may also provide
                                the IP address to other Chewie-NS instances (default: main).
