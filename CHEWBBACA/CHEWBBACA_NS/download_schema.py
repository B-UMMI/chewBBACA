#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------

This module enables the download of schemas from a Chewie-NS instance.
The process enables the download of ZIP archives that contain ready-to-use
versions of any schema in Chewie-NS. It also allows users to download
any schema with the structure it had at a specific time point. It is also
possible to download the latest version of the schema through requests to
the Chewie-NS API, if the compressed version that is available does not
match the latest version of the schema. An alternative approach that can
be applied to get the latest version of the schema, if the compressed version
does not provide it, is to download the compressed version that is available
and run the SyncSchema process to retrieve the alleles that were added to the
schema after the compression date.

Code documentation
------------------
"""


import os
import shutil
import requests
import concurrent.futures
from itertools import repeat
from urllib3.exceptions import InsecureRequestWarning

from loguru import logger

import CHEWBBACA.global_constants as gct
import CHEWBBACA.CHEWBBACA_NS.constants as ct
import CHEWBBACA.CHEWBBACA_NS.validation_models as vm
from CHEWBBACA.AdaptSchema import adapt_schema
from CHEWBBACA.AdaptSchema.validation_models import AdaptSchemaValidator
from CHEWBBACA.utils import (
    core_functions as cf,
    file_operations as fo,
    process_datetime as pd,
    chewiens_requests as cr)


# Suppress only the single warning from urllib3 needed.
requests.packages.urllib3.disable_warnings(category=InsecureRequestWarning)


def check_compressed(schema_uri, headers_get):
    """Determine if there is a compressed version of a schema.

    Parameters
    ----------
    schema_uri : str
        The URI of the schema in the Chewie-NS.
    headers_get : dict
        HTTP headers for GET requests.

    Returns
    -------
    list
        A list with the following elements:

        - The URI for the compressed version of the schema (str).
        - The timestamp of the compressed version. Indicates the
          last modification date of the schema at time of compression.
    """
    zip_uri, zip_response = cr.simple_get_request(schema_uri, headers_get, ["zip"], parameters={"request_type": "check"})
    zip_info = zip_response.json()
    if "zip" in zip_info:
        zip_file = zip_info["zip"][0]
        zip_date = zip_file.split("_")[-1].split(".zip")[0]
    else:
        zip_date = None

    return [zip_uri, zip_date]


def download_date(user_date, zip_date, latest, insertion_date,
                  modification_date):
    """Determine the date that will be used as timepoint for schema download.

    Parameters
    ----------
    user_date : str or None
        Date value that was received by :py:func:`main`.
        If the user provided a value, it should be a date
        in the format %Y-%m-%dT%H:%M:%S or %Y-%m-%dT%H:%M:%S.%f.
        If the user did not provide a value for the argument, it
        will be None by default.
    zip_date : str or None
        Date value of the compressed version of the schema that
        is available for download. Will be in the format
        %Y-%m-%dT%H:%M:%S or %Y-%m-%dT%H:%M:%S.%f if there is
        a compressed version. Will be None if there is no compressed
        version.
    latest : bool
        Boolean value that indicates if the user wants to download
        the latest version of the schema.
    insertion_date : str
        Date on which the initial schema upload/insertion was
        completed.
    modification_date : str
        Date on which the schema was last modified.

    Returns
    -------
    schema_date : str
        Date selected as timepoint. Schema will be downloaded with
        the data structure it had on this date.

    Raises
    ------
    SystemExit
        - If the value received for the `user_date` is not a
          valid date in the format %Y-%m-%dT%H:%M:%S or
          %Y-%m-%dT%H:%M:%S.%f.
        - If the value received for the `user_date` is a date
          that is prior to `insertion_date` or greater than the
          `modification_date`.
    """
    # user did not provide a value for 'date' or 'latest' arguments
    if user_date is None and latest is False:
        if zip_date is not None:
            schema_date = zip_date
        # there is no compressed version
        else:
            schema_date = modification_date
    # user wants latest schema version
    elif latest is True:
        schema_date = modification_date
    # user wants schema at a particular time point
    elif user_date is not None:
        # get schema insertion and last modification date
        insertion_date_obj = pd.datetime_obj(insertion_date, "%Y-%m-%dT%H:%M:%S.%f")
        modification_date_obj = pd.datetime_obj(modification_date, "%Y-%m-%dT%H:%M:%S.%f")

        if user_date is False:
            logger.warning("Provided date is invalid. Please provide a date in the format '%Y-%m-%dT%H:%M:%S' or '%Y-%m-%dT%H:%M:%S.%f'")
            raise vm.DownloadSchemaInvalidDateFormat
        if user_date >= insertion_date_obj and user_date <= modification_date_obj:
            schema_date = pd.datetime_str(user_date, "%Y-%m-%dT%H:%M:%S.%f")
        elif user_date < insertion_date_obj:
            logger.warning(f"Provided date is prior to the date of schema insertion. Please provide a date later than schema insertion date ({insertion_date}).")
            raise vm.DownloadSchemaInvalidDate
        elif user_date > modification_date_obj:
            logger.warning(f"Provided date is greater than the last modification date. Please provide a date equal or prior to the schema modification date ({insertion_date}).")
            raise vm.DownloadSchemaInvalidDate

    return schema_date


def build_fasta(locus_id, locus_info, download_folder):
    """Write DNA sequences from a response object into a FASTA file.

    Parameters
    ----------
    locus_id : str
        Name of the locus in the Chewie-NS.
    locus_info : dict
        Response object with the list of DNA
        sequences of the locus.
    download_folder : str
        Path to the directory where the FASTA
        file with the locus sequences will be
        created.

    Returns
    -------
    locus_file : str
        Path to the FASTA file with the locus
        sequences.
    """
    locus_name = locus_id.rstrip(".fasta")
    locus_file = os.path.join(download_folder, locus_id+".fasta")
    ns_data = locus_info.json()["Fasta"]

    # extract identifiers and DNA sequences from response
    locus_alleles = []
    for allele in ns_data:
        allele_id = int(allele["allele_id"]["value"])
        allele_seq = allele["nucSeq"]["value"]
        locus_alleles.append((allele_id, allele_seq))

    # create FASTA records
    records = [">{0}_{1}\n{2}".format(locus_name, a[0], a[1])
               for a in locus_alleles]

    # write sequences to FASTA file
    with open(locus_file, "w") as lf:
        concat_records = "\n".join(records)
        lf.write(concat_records)

    return locus_file


def get_fasta_seqs(url, headers_get, schema_date):
    """Retrieve the DNA sequences of a locus in Chewie-NS.

    Parameters
    ----------
    url : str
        Endpoint URL to make the request.
    headers_get : dict
        HTTP headers for GET requests.
    schema_date : str
        The function will only retrieve alleles
        that were inserted up to this date.

    Returns
    -------
    tuple
        Tuple with the following elements:
        - URI of the locus.
        - Response object with the DNA sequences
          that were downloaded.
    """
    payload = {"date": schema_date}
    tries = 0
    max_tries = 3
    downloaded = False
    while downloaded is False:
        res = cr.simple_get_request(url, headers_get, [], payload, False, 180)[1]
        tries += 1
        if res.status_code in [200, 201] or tries == max_tries:
            downloaded = True

    return (url.rstrip("/fasta"), res)


def schema_loci(schema_uri, headers_get):
    """Retrieve the list of loci for a schema.

    Parameters
    ----------
    schema_uri : str
        The URI of the schema in the Chewie-NS.
    headers_get : dict
        HTTP headers for GET requests.

    Returns
    -------
    loci : dict
        A dictionary with loci URIs as keys and
        loci names as values.
    """
    # get the list of loci
    loci_url, loci_res = cr.simple_get_request(schema_uri, headers_get, ["loci"])
    loci_res = loci_res.json()["Loci"]

    # locus URI to locus name
    loci = {}
    for locus in loci_res:
        loci[str(locus["locus"]["value"])] = locus["name"]["value"]

    return loci


def download_fastas(loci, download_folder, headers_get, schema_date):
    """Download and write FASTA files for the loci of a schema in Chewie-NS.

    Parameters
    ----------
    loci : dict
        A dictionary with loci URIs as keys and
        loci names as values.
    download_folder : str
        Path to the directory where the FASTA files
        will be created.
    headers_get : dict
        HTTP headers for GET requests.
    schema_date : str
        The function will only retrieve alleles
        that were inserted up to this date.

    Returns
    -------
    ns_files : list
        List with the paths to the schema's FASTA
        files that were created.
    """
    # Total number of loci
    total_loci = len(loci)
    logger.info(f"Number of loci to download: {total_loci}")

    # Build the list of urls to get
    fasta_urls = [cr.make_url(locus, "fasta") for locus in loci]

    # Multithread the requests
    logger.info("Downloading schema files...")
    total = 0
    failed = []
    downloaded = 0
    ns_files = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        for result in executor.map(get_fasta_seqs, fasta_urls, repeat(headers_get), repeat(schema_date)):
            locus_id = loci[result[0]]
            locus_info = result[1]
            if locus_info.status_code in [200, 201]:
                locus_file = build_fasta(locus_id, locus_info, download_folder)
                ns_files.append(locus_file)
                downloaded += 1
            else:
                failed.append(locus_id)
            total += 1
            logger.info("Downloaded: {0}/{1}".format(downloaded, total_loci))

    logger.info(f"Downloaded and wrote FASTA files for {downloaded}/{total} loci.")
    logger.info(f"Failed download for {failed} loci.")
    if len(failed) > 0:
        logger.error(f"Failed download for following loci: {','.join(failed)}\n"
                     "Please download files for failed loci through the API or retry full schema download.")
        raise vm.DownloadSchemaFailedLociFastaDownload

    return ns_files


def download_compressed(zip_uri, schema_filename, download_folder, headers_get):
    """Download and extract ZIP archive with a ready-to-use schema.

    Parameters
    ----------
    zip_uri : str
        Endpoint URL to make the request to download
        the compressed schema.
    schema_filename : str
        Name assigned to the local schema.
    download_folder : str
        Path to the directory to which the ZIP archive
        will be saved.
    headers_get : dict
        HTTP headers for GET requests.

    Returns
    -------
    schema_path : str
        ZIP archive contents will be extracted to this
        directory.
    """
    schema_path = os.path.join(download_folder, schema_filename)
    fo.create_directory(schema_path)

    # Download ZIP archive
    zip_name = f"{schema_filename}.zip"
    url, zip_response = cr.simple_get_request(zip_uri, headers_get, parameters={"request_type": "download"})
    zip_path = os.path.join(schema_path, zip_name)
    open(zip_path, "wb").write(zip_response.content)
    # Uncompress
    logger.info("Decompressing schema...")
    shutil.unpack_archive(zip_path, extract_dir=schema_path)
    # Delete ZIP
    os.remove(zip_path)

    return schema_path


def download_ptf(download_folder, schema_id, species_id, species_name, headers_get, base_url):
    """Download the Prodigal training file for a schema.

    Parameters
    ----------
    download_folder : str
        Path to the directory to which the Prodigal
        training file should be saved.
    schema_id : str
        The identifier of the schema in the Chewie-NS.
    species_id : str
        The identifier of the schema's species in the
        Chewie-NS.
    species_name : str
        Scientific name of the schema species.
    headers_get : dict
        HTTP headers for GET requests.
    base_url : str
        Base URL of the Chewie Nomenclature server.

    Returns
    -------
    ptf_file : str
        Path to the Prodigal training file.
    """
    _, ptf_response = cr.simple_get_request(base_url, headers_get, ["species", species_id, "schemas", schema_id, "ptf"])

    ptf_file = os.path.join(download_folder, "{0}.trn".format(species_name.replace(" ", "_")))

    open(ptf_file, "wb").write(ptf_response.content)

    return ptf_file


def main(config):
    # GET request headers
    headers_get = gct.HEADERS_GET_JSON

    logger.info(f"Schema id: {config.schema_id[0]}")
    logger.info(f"Schema name: {config.schema_id[2]["name"]}")
    logger.info(f"Schema's species: {config.species_id[1]} (id={config.species_id[0]})")

    # Assign name to local schema
    schema_basename = "_".join(config.species_id[1].split(" "))
    schema_basename = f"{schema_basename}_{config.schema_id[2]["name"]}"

    # Check if schema is locked
    lock_status = config.schema_id[2]["Schema_lock"]
    if lock_status != "Unlocked":
        logger.warning("Schema is locked. This may be because it is being uploaded, updated or compressed."
                       " Please try again later and contact the Administrator if the schema stays locked for "
                       "a long period of time.")
        raise vm.SchemaLocked

    # Get zip information
    zip_uri, zip_date = check_compressed(config.schema_id[1], headers_get)

    schema_date = download_date(config.schema_date, zip_date, config.latest_version,
                                config.schema_id[2]["dateEntered"],
                                config.schema_id[2]["last_modified"])

    if schema_date == zip_date:
        logger.info("Downloading compressed version...")
        # Chewie-NS does not add clustering parameters to config, change that
        schema_path = download_compressed(zip_uri, schema_basename, config.output_directory, headers_get)
    else:
        logger.info("Downloading schema FASTA files...")
        # Download schema FASTA files
        loci = schema_loci(config.schema_id[1], headers_get)
        # Create folder to store FASTA files
        fastas_folder = fo.join_paths(config.output_directory, ["loci_fastas"])
        fo.create_directory(fastas_folder)
        ns_files = download_fastas(loci, fastas_folder, headers_get, schema_date)

        # Use PrepExternalSchema to determine representatives
        logger.info("Calling the AdaptSchema module to determine a set of representative alleles and create the schema structure...")
        adapt_schema_outdir = fo.join_paths(config.output_directory, [schema_basename])
        adapt_schema_config = AdaptSchemaValidator(adapt_schema_outdir,
                                                   fastas_folder,
                                                   blast_score_ratio=float(config.schema_id[2]["bsr"]),
                                                   minimum_length=int(config.schema_id[2]["minimum_length"]),
                                                   adaptation_minimum_length=0,
                                                   translation_table=int(config.schema_id[2]["translation_table"]),
                                                   size_threshold=float(config.schema_id[2]["size_threshold"]),
                                                   adaptation_size_threshold=0,
                                                   cpu_cores=config.cpu_cores,
                                                   blast_path=config.blast_path)

        adapt_schema.main(adapt_schema_config)
        schema_path = fo.join_paths(adapt_schema_outdir, ["adapted_schema"])

        # Download Prodigal training file
        ptf_hash = config.schema_id[2]["prodigal_training_file"]

        ptf_file = download_ptf(schema_path, config.schema_id[0],
                                config.species_id[0], config.species_id[1], headers_get,
                                config.nomenclature_server_instance)

        # Delete folder with FASTA files downloaded from Chewie-NS
        fo.delete_directory(fastas_folder)

        # Write schema config file
        config.schema_id[2]["ptf_path"] = ptf_hash
        config.schema_id[2]["window_size"] = config.schema_id[2].get("window_size", None)
        schema_config = cf.write_schema_config(config.schema_id[2],
                                               config.schema_id[2]["chewBBACA_version"],
                                               schema_path)

        # Create ns_config file
        ns_config = os.path.join(schema_path, ".ns_config")
        download_info = [schema_date, config.schema_id[1]]
        if not os.path.exists(ns_config):
            fo.pickle_dumper(download_info, ns_config)

        genes_list_file = cf.write_gene_list(schema_path)

    logger.info(f"Schema is now available at: {schema_path}")
