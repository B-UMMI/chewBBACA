#!/usr/bin/env python3

import itertools
import multiprocessing
import os
import re
import shutil

from Bio import SeqIO


def verify_cpu_usage(cpu_to_use: int) -> int:
    """ Verify the cpu usage for chewBBACA.

        Args:
            cpu_to_use (int): the number of cpu provided to chewBBACA

        Returns:
            cpu_to_use (int): the number of cpu to use after verification

        Example:

            >>> verify_cpu_usage(6)
            6
    """
    total_cpu = multiprocessing.cpu_count()

    # do not allow a value of cpuToUse greater than the number of cores/threads
    if cpu_to_use > total_cpu:
    	print("Warning! You have provided a CPU core count value that exceeds the number of cores in your machine!")
    	print("Setting a different value for the CPU core count...")
    	# define a value that is safe according to the number of available cores/threads
    	if total_cpu > 2:
    		cpu_to_use = total_cpu - 2
    	elif total_cpu == 2:
    		cpu_to_use = 1
    	print("CPU core count value set to: ", cpu_to_use)
    
    elif cpu_to_use < total_cpu and cpu_to_use > total_cpu - 2:
        print("Warning! You have provided a CPU core count value that is close to the maximum core count of your machine (" \
        	+ str(cpu_to_use) + '/' + str(total_cpu) + "). This may affect your system responsiveness.")

    return cpu_to_use


def check_prodigal_output_files(path_to_temp: str, list_of_genomes):
	""" Checks if Prodigal created ORF files 
		equal to the number of genome files provided.

        Args:
            path_to_temp (str): the full path to the 'temp' directory created by chewBBACA.
            list_of_genomes (list): list containing the full path to the input genomes.

        Returns:
            prints a message if Prodigal created all the necessary files, 
            otherwise raises a ValueError 

	"""
	
	# list ORF files created in the "temp" directory
	listOfORFCreated = []
	for orffile in os.listdir(path_to_temp):
		if orffile.endswith("_ORF.txt"):
			listOfORFCreated.append(orffile)

	# raise exception if the number of ORF files is not equal to the number of genome files provided
	if len(list_of_genomes) > len(listOfORFCreated):
		message = "Missing some ORF files from Prodigal run. Missing " + str((len(list_of_genomes)) - (len(listOfORFCreated))) \
		+ " ORF files out of " + str(len(list_of_genomes)) + " genome files provided."
		# remove "temp" directory
		shutil.rmtree(path_to_temp)
		raise ValueError(message)
	else:
		print("Prodigal created all the necessary files.")
	

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
            with open("listGenomes2Call.txt", "w") as f:
                for genome in fasta_files:
                    f.write(genome + "\n")
        else:
            raise Exception("There were no FASTA files in the given directory. Please provide a directory \
                            with FASTA files or a file with the list of full paths to the FASTA files.")

        list_files = "listGenomes2Call.txt"
    
    else:
        raise Exception("Input argument is not a valid directory or file with a list of paths. \
                        Please provide a valid input, either a folder with FASTA files or a file with \
                        the list of full paths to FASTA files (one per line).")

    return list_files


def escape_special_characters(a_string: str):
    """ Escapes strings to use in regex

        Args:
            a_string (str): string containing characters to escape

        Returns:
            escaped (str): escaped string
    """

    escaped = re.escape(a_string)

    return escaped

def replace_multiple_characters(namefile: str):
    """ Replaces multiple characters in a string

        Args:
            namefile (str): string containing the name of the contig 
            with characters to replace

        Returns:
            replaced (str): string containing the name of the contig 
            without characters to replace
    """

    replaced = namefile.replace("|", "_")\
                       .replace("_", "-")\
                       .replace("(", "")\
                       .replace(")", "")\
                       .replace("'", "")\
                       .replace("\"", "")\
                       .replace(":", "")

    return replaced


def listdir_fullpath(path: str):
    """ Gets the full path of the files from a directory

        Args:
            path (str): full path to a directory

        Returns:
            list containing the full path of every file contained in the input directory
    
    """
    
    return [os.path.join(path, f) for f in os.listdir(path)]


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


def invert_dictionary(dictionary):
    """ Inverts a dictionary. Keys become values and vice-versa

        Args:
            dictionary (dict)

        Returns:
            inverted (dict): inverted dictionary
        
        Example:

            >>> inverted_dictionary({key:value})
            {value:key}
    """
    
    inverted = {value:key for key, value in dictionary.items()}
    
    return inverted

def threads_for_blast(files_to_blast, cpu_to_apply):
    """ Define the number of threads for BLAST

        Args:
            files_to_blast (list): list containing the full path to the files to BLAST
            cpu_to_apply (int): number of cpu to use

        Returns:
            blast_threads (list): list contaning the number of threads to use for each file
            proc (int): Number of processes to use in multiprocessing
    """
    
    # define number of processes and available cores for each BLAST
    if len(files_to_blast) >= cpu_to_apply:
        blast_threads = [1 for protogenome in files_to_blast]
        proc = cpu_to_apply
    
    elif cpu_to_apply % len(files_to_blast) == 0:
        blast_threads = [int(cpu_to_apply / len(files_to_blast)) for protogenome in files_to_blast]
        proc = len(blast_threads)
    
    elif cpu_to_apply % len(files_to_blast) == 1:
        blast_threads = [2] + [1 for protogenome in range(0,len(files_to_blast)-1)]
        proc = len(blast_threads)
    
    elif cpu_to_apply % len(files_to_blast) > 1:
        base_cpu = int(cpu_to_apply / len(files_to_blast))
        blast_threads = [base_cpu for protogenome in range(0,len(files_to_blast))]
        extra_cpu = cpu_to_apply - sum(blast_threads)
        i = 0
        while extra_cpu > 0:
            blast_threads[i] += 1
            extra_cpu -= 1
            i += 1
        proc = len(blast_threads)

    return blast_threads, proc


