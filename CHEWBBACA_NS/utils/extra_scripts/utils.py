#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 2 15:37:42 2019

@author: pcerqueira
"""

import os
import time
import json
import requests
import itertools
import multiprocessing
from urllib.parse import urlparse, urlencode, urlsplit, parse_qs

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna



def verify_cpu_usage(cpuToUse):
    """ Verify the cpu usage for chewBBACA
    """
    total_cpu = multiprocessing.cpu_count()

    # do not allow a value of cpuToUse greater than the number of cores/threads
    if cpuToUse > total_cpu:
    	print("Warning! You have provided a CPU core count value that exceeds the number of cores in your machine!")
    	print("Setting a different value for the CPU core count...")
    	# define a value that is safe according to the number of available cores/threads
    	if total_cpu > 2:
    		cpuToUse = total_cpu - 2
    	elif total_cpu == 2:
    		cpuToUse = 1
    	print("CPU core count value set to: ", cpuToUse)
    
    if cpuToUse > total_cpu - 2:
        print("Warning! You have provided a CPU core count value that is close to the maximum core count of your machine (" \
        	+ str(cpuToUse) + '/' + str(total_cpu) + "). This may affect your system responsiveness.")

    return cpuToUse


def flatten_list(list_to_flatten):
    """Flattens one level of a nested list

        Args:
            list_to_flatten (list)

        Returns:
            flattened list

        Example:

            >>> flatten_list([[[1,2],[3,4]]])
            [[1, 2], [3, 4]]

    """

    return list(itertools.chain(*list_to_flatten))


def translateSeq(DNASeq):

    seq = DNASeq
    tableid = 11
    try:
        myseq = Seq(seq)
        protseq = Seq.translate(myseq, table=tableid, cds=True)
    except:
        try:
            seq = reverseComplement(seq)
            myseq = Seq(seq)
            protseq = Seq.translate(myseq, table=tableid, cds=True)
        except:
            try:
                seq = seq[::-1]
                myseq = Seq(seq)
                protseq = Seq.translate(myseq, table=tableid, cds=True)
            except:
                try:
                    seq = seq[::-1]
                    seq = reverseComplement(seq)
                    myseq = Seq(seq)
                    protseq = Seq.translate(myseq, table=tableid, cds=True)
                except Exception as e:
                    # print("translation error")
                    # print(e)
                    raise

    return myseq


def reverseComplement(strDNA):

    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    strDNArevC = ''
    for l in strDNA:
        strDNArevC += basecomplement[l]

    return strDNArevC[::-1]


def is_fasta(filename: str):
    """ Checks if a file is a FASTA file.

        Args:
            filename (str): the full path to the FASTA file

        Returns:
            True if FASTA file,
            False otherwise
    
    """
    
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        
        # returns True if FASTA file, False otherwise
        return any(fasta)
    


def is_url(url):
    """ Checks if a url is valid
    
        Args: 
        url (str): the url to be checked

        Returns:
        True if url is valid, False otherwise.
    
    """
    
    try:
        
        result = urlparse(url)
        return all([result.scheme, result.netloc, result.path])
    
    except:
        return False



def check_if_list_or_folder(folder_or_list: str):
    """ Checks if the input is a file or a directory.

        Args: 
            folder_or_list (str): the full path to the file or directory

        Returns:
            list_files (str) if folder_or_list is a path to a file,
            list_files (list) if folder_or_list is a path to a directory,
            Raises Exception otherwise
    """
    
    # check if input argument is a file or a directory
    if os.path.isfile(folder_or_list):
        list_files = folder_or_list
    
    elif os.path.isdir(folder_or_list):
        
        fasta_files = []
            
        for genome in os.listdir(folder_or_list):
                
            genepath = os.path.join(folder_or_list, genome)
            
            # do not continue if genepath is a dir
            if os.path.isdir(genepath):
                continue
            
            # check if file is a FASTA file
            if is_fasta(genepath):
                fasta_files.append(os.path.abspath(genepath))
        
        # if there are FASTA files
        if fasta_files:
            # store full paths to FASTA files
            with open("listGenes.txt", "w") as f:
                for genome in fasta_files:
                    f.write(genome + "\n")
        else:
            raise Exception("There were no FASTA files in the given directory. Please provide a directory \
                            with FASTA files or a file with the list of full paths to the FASTA files.")

        list_files = "listGenes.txt"
    
    else:
        raise Exception("Input argument is not a valid directory or file with a list of paths. \
                        Please provide a valid input, either a folder with FASTA files or a file with \
                        the list of full paths to FASTA files (one per line).")

    return list_files



def make_url(base_url , *res, **params):
    """ Creates a url. 
    
        Args: 
            base_url (str): the base url
            res (str): endpoint(s) to add to the base url
            params (str): addtional parameters (WIP)

        Returns:
            url (str) with the provided parameters.
            Otherwise, returns base_url.

    """
    
    url = base_url
    
    # Check if the url is valid
    if is_url(url):
        
        if url[-1] == "/":
            url = url[:-1]
    
        # Add the endpoints
        for r in res:
    #        url = '{}/{}'.format(url, r)
            url = f'{url}/{r}'
        
        # Add params if they are provided
        if params:
    #        url = '{}?{}'.format(url, urllib.urlencode(params))
            url = f'{url}?{urlencode(params)}'
        
        return url
    
    else:
        return "An invalid URL was provided."
    

def get_sequence_from_url(url):
    """
    """
    
    seq = parse_qs(urlsplit(url).query)["sequence"][0]
    
    return seq


def login_user_to_NS(server_url, email, password):
    """ Logs a user in Nomenclature Server
    
        Args:
            server_url (str): url of Nomeclature Server API
            email (str): email of the user in NS
            password (str): password of the user in NS
            
        Returns:
            token (str): authorization token to perform requests to NS
    """
    
    auth_params = {}
    auth_params["email"] = email 
    auth_params["password"] = password
    
    auth_headers = {}
    auth_headers["Content-Type"] = "application/json"
    auth_headers["accepts"] = "application/json"
    
    auth_url = make_url(server_url, "auth", "login")
    
    auth_r = requests.post(auth_url, data=json.dumps(auth_params), headers=auth_headers)
    
    auth_result = auth_r.json() 
    
    token = auth_result["X-API-KEY"]
    
    return token
    

def send_data(sparql_query, url_send_local_virtuoso, virtuoso_user, virtuoso_pass):
    """ Sends data to virtuoso.
    
        Args: 
            sparql_query (str): the sparql query to use
            url_send_local_virtuoso (str): the url for the local Virtuoso folder
            virtuoso_user (str): the Virtuoso user
            virtuoso_pass (str): the Virtuoso password

        Returns:
            r (Response) of the Virtuoso server
    """
    
    url = url_send_local_virtuoso
    headers = {'content-type': 'application/sparql-query'}
    r = requests.post(url, data=sparql_query, headers=headers, auth=requests.auth.HTTPBasicAuth(virtuoso_user, virtuoso_pass))

	#sometimes virtuoso returns 405 God knows why ¯\_(ツ)_/¯ retry in 2 sec
    if r.status_code >201:
        time.sleep(2)
        r = requests.post(url, data=sparql_query, headers=headers, auth=requests.auth.HTTPBasicAuth(virtuoso_user, virtuoso_pass))
        
    return r


def send_post(loci_uri, sequence, token, noCDSCheck):
    """ """

    params = {}
    params['sequence'] = sequence

    if not noCDSCheck:
        params['enforceCDS'] = "False"
    
    headers = {'X-API-KEY': token,
               'Content-type': 'application/json',
               'accept': 'application/json'}
    
#    url = loci_uri + "/alleles"
    
    url = make_url(loci_uri, "alleles")

    req_success = False
    sleepfactor = 4
    while not req_success:
        try:
            r = requests.post(url, data=json.dumps(params), headers=headers, timeout=30)
            
            if r.status_code == 418:
                print("Sequence is already attributed to a loci/allele")
            
            elif r.status_code > 201:
                print(r)
                print("failed sending sequence, retrying in seconds "
                      + str(sleepfactor))
                time.sleep(sleepfactor)
                sleepfactor = sleepfactor * 2
            else:
                req_success = True
        except:
            time.sleep(sleepfactor)
            sleepfactor = sleepfactor * 2
            pass

    req_code = int(r.status_code)
    # allele_url=((r.content).decode("utf-8")).replace('"', '').strip()

    return req_code


def send_sequence(token, sequence, loci_uri, noCDSCheck):
    """ """

    req_success = False
    sleepfactor = 4
    while not req_success:

        reqCode = send_post(loci_uri, sequence, token, noCDSCheck)
        if reqCode > 201:
            print("failed, retrying in seconds "+str(sleepfactor))
            time.sleep(sleepfactor)
            sleepfactor = sleepfactor * 2

        else:
            req_success = True
	
	# if reqCode==401:
		# print ("Token is not valid")
	# elif reqCode>201:
		#
		# try:
			#~ allele_url,reqCode=send_post(loci_uri,sequence,token)
		# except:
			#~ print ("Server returned code "+str(reqCode))
			#~ print(loci_uri)
	# else:
		# new_allele_id=str(int(allele_url.split("/")[-1]))
        
    return reqCode


def process_locus(gene, token, loci_url, auxBar, noCDSCheck):
    """ """

    for allele in SeqIO.parse(gene, "fasta", generic_dna):

        sequence = (str(allele.seq)).upper()
        try:
            sequence = translateSeq(sequence)
            sequence = str(sequence)

        except:
            continue

        #reqCode = send_sequence(token, sequence, loci_url, noCDSCheck)
        reqCode = send_post(loci_url, sequence, token, noCDSCheck)
        
#    if reqCode == 418:
#        print(gene)

    if gene in auxBar:
        auxlen = len(auxBar)
        index = auxBar.index(gene)
        print("[" + "=" * index + ">" +
                " " * (auxlen - index) +
                "] Sending alleles " +
                    str(int((index / auxlen) * 100)) + "%")

    return reqCode


#def process_locus2(sequence, token, loci_url, auxBar, noCDSCheck):
#    """ """
#
#    #reqCode = send_sequence(token, sequence, loci_url, noCDSCheck)
#    reqCode = send_post(loci_url, sequence, token, noCDSCheck)
#        
#
#
##    if gene in auxBar:
##        auxlen = len(auxBar)
##        index = auxBar.index(gene)
##        print("[" + "=" * index + ">" +
##                " " * (auxlen - index) +
##                "] Sending alleles " +
##                    str(int((index / auxlen) * 100)) + "%")
#
#    return reqCode

 
