#!/usr/bin/env python3

import itertools
import multiprocessing
import os
import re
import shutil

from Bio import SeqIO


def check_prodigal_output_files(prodigal_path, fasta_files, genomes_dir, prodigal_results, genomes_identifiers, parent_dir):
    """ Checks if Prodigal created ORF files
        equal to the number of genome files provided.

        Args:
            path_to_temp (str): the full path to the 'temp' directory created by chewBBACA.
            list_of_genomes (list): list containing the full path to the input genomes.
    """

    no_cds = [l for l in prodigal_results if l[1] == 0]
    errors = [l for l in prodigal_results if isinstance(l[1], str) is True]
    failed = no_cds + errors

    if len(failed) > 0:
        outfile = os.path.join(parent_dir, 'prodigal_fails.tsv')
        with open(outfile, 'w') as pf:
            lines = ['{0}\t{1}'.format(l[0], l[1]) for l in failed]
            pf.writelines(lines)

        # remove failed genomes from paths
        for f in failed:
            file_path = os.path.join(genomes_dir, '{0}.fasta'.format(f[0]))
            fasta_files.remove(file_path)
            genomes_identifiers.remove(f[0])

    return [fasta_files, genomes_identifiers]


def is_fasta(filename):
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


