#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR
    
    Rafael Mamede
    github: @rfm-targa

DESCRIPTION





"""


import os
import csv
import pickle
import itertools

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna

import new_utils
import init_schema_4_bbaca


def genome_id(fasta_path):
    """ Get genome/assembly identifier from full path.
    
        Args: 
            fasta_path (str): the full path to the FASTA file.
        
        Returns: 
            file (str): the last string resulting from the splitting of the path by 
            "/" and stripped of the file extension suffix.
        
        Example:
            
            >>> get_id("/home/user/wd/file.fasta")
            file
        
    """
    
    # possible FASTA file suffixes
    extension_format = ['.fasta','.fna','.ffn']
    
    file = fasta_path.split('/')[-1]
    
    for extension in extension_format:
        if extension in file:
            file = file.rstrip(extension)
    
    file = file.split('.')[0]
            
    return file
    

# Import contigs for each genome
def import_contigs(fasta_path):
    """ Imports contigs from a FASTA file.
    
        Args: 
            fasta_path (str): full path to the FASTA file.
        
        Returns: 
            dictionary that has contigs ids as keys and contigs DNA 
            sequences as values.
            
        Example:
            
            >>> import_contigs("/home/user/wd/file.fasta")
            {contig_1:'TCGAACCACCGACCTCACGCTTATCAGG...',
            contig_2:'ATAAATGGCGCGAGACGGAATCGAACCGCCGA...'
            ...}
            
    """
    
    contigs_dict = {}
    # use BioPython to read FASTA file and get each contig sequence
    for contig in SeqIO.parse(fasta_path, 'fasta', generic_dna):
        # seq object has to be converted to string
        sequence = str(contig.seq.upper())
        contig_id = contig.id
        
        # add contig id as key and DNA sequence as value
        contigs_dict[contig_id] = sequence
    
    return contigs_dict


def extract_coding_sequences(reading_frames, contigs, starting_id):
    """ Extracts CDSs from contigs based on the start codon and stop codon 
        positions determined by Prodigal.
        
        Args:
            orf_file_path (str): full path to the ORF file derived from Prodigal.
            contigs (dict): a dictionary with contigs ids as keys and contigs 
            sequences as values.
            starting_protid (int): integer identifier to give to the first CDS
            extracted and that will be incremented to serve as identifier for
            subsequent CDSs.
            genome_id (str): id of the genome or assembly.
            
        Returns:
            cds_dict (dict): dictionary with CDSs ids as keys and CDSs DNA 
            sequences as values.
            cdss_lines (list): list of lists where each sublist has information 
            about the CDS.
            cdss_contigs (dict): dictionary with CDSs ids as keys and contigs 
            ids as values.
            protid (int): last extracted CDS id + 1. integer value that might 
            serve as starting id if the function will be used to extract CDSs 
            from more genomes/assemblies.
        
        Example:
            >>> extract_cdss("/home/user/wd/file_ORF.txt", contigs, 1, "genome1")
            {1: 'GTGACACCAAAACCAGTACCGGATAAA...',
            2: 'TTAGTCTAATTCTATTTGGAGAAATTTA...',
            ...},
            [['genome1','contig_1','102','528','1'],
            ['genome1','contig_1','1861','2110','2'],
            ...],
            {1: 'contig_1',
            2: 'contig_1,
            ...},
            1914
            
    """
    
    # load binary file with a list of lists for each contig
    # each sublist has a start codon and stop codon positions in the contig
    with open(reading_frames, 'rb') as orf_file:
        start_stop_codons = pickle.load(orf_file)

    protid = starting_id
    coding_sequences = {}
    coding_sequences_info = []
    # for each contig
    for contig_id, frames in start_stop_codons.items():
        # for each start and stop codon in that contig
        for coding_sequence in frames:
            start_codon = coding_sequence[0]
            stop_codon = coding_sequence[1]
            strand = coding_sequence[2]
            # extract CDS sequence
            cds_sequence = contigs[contig_id][start_codon:stop_codon].upper()
            if strand == 0:
                cds_sequence = reverse_complement(cds_sequence)

            # store CDS with unique id
            coding_sequences[protid] = cds_sequence

            # store CDS information
            coding_sequences_info.append([contig_id, str(start_codon),
                                          str(stop_codon), str(protid),
                                          str(strand)])

            # increment the CDS id by 1 so that it can be an unique identifier
            # for the next CDS
            protid += 1

    return [coding_sequences, coding_sequences_info]


def reverse_str(string):
    """ Reverse character order in input string.

        Args:
            string (str): string to be reversed.

        Returns:
            revstr (str): reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Args:
            strDNA (str): string representing a DNA sequence.

        Returns:
            revC_dna (str): the reverse complement of the DNA sequence, without
            lowercase letters.

        Example:
            >>> reverse_complement('ATCGgcaNn')
            'NNTGCCGAT'
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    complement_bases = []
    for base in bases:
        if base in base_complement:
            complement_bases.append(base_complement[base])
        else:
            complement_bases.append(base.upper())

    complement_strand = ''.join(complement_bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.

        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Tries to translate an input DNA sequence in specified orientation
        and stores exceptions when the input sequence cannot be translated.

        Args:
            dna_sequence (str): string representing a DNA sequence.
            method (str): a string specifying the way the sequence will
            be oriented to attempt translation.
            table_id (int): translation table identifier.

        Returns:
            List with following elements if translation is successful:
                protseq (str): string representing the translated DNA sequence.
                myseq (str): string representing the DNA sequence in the
                orientation used to translate it.
            Otherwise, returns string derived from captured exception.
    """

    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def check_str_alphabet(string, alphabet):
    """ Determine if a string only has characters from specified
        alphabet.

        Args:
            string (str): input string.
            alphabet (str): string that has all characters from desired
            alphabet.

        Returns:
            "True" if sequence only has characters from specified
            alphabet and string "ambiguous or invalid characters" if
            it any of its characters is not in the alphabet.
    """

    valid_chars = alphabet
    if all(n in valid_chars for n in string) is True:
        return True
    else:
        return 'ambiguous or invalid characters'


def check_str_multiple(string, number):
    """ Determine if length of input string is multiple of
        a specified number.

        Args:
            string (str): input string.
            number (int): integer that will be used to check if sequence
            length is multiple of.

        Returns:
            "True" if the length of the sequence is a multiple of the
            specified number and "sequence length is not a multiple of number"
            if condition is not satisfied.
    """

    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def translate_dna(dna_sequence, table_id):
    """ Checks if sequence is valid and attempts to translate it,
        calling several functions to ensure that the sequence only has
        'ACTG', is multiple of 3 and that it can be translated in any of 4
        different orientations. Stores exceptions so that it is possible to
        understand the sequence could not be translated.

        Args:
            dna_sequence (str):
            table_id (int):

        Returns:
            If the sequence can be translated,
            a list with following elements:
                sequence (list): a list with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
                coding_strand (str): the strand orientation that had could be
                translated.
            Otherwise:
                exception_str (str): a string containing the exceptions that
                determined that the sequence could not be translated.
    """

    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the string is DNA, without ambiguous bases
    valid_dna = check_str_alphabet(original_seq, 'ACTG')
    if valid_dna is not True:
        return valid_dna

    # check if sequence is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return valid_length

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        sequence, exception_collector = retranslate(original_seq,
                                                    translating_methods[i],
                                                    table_id, strands[i],
                                                    exception_collector)

        i += 1
        if i == len(strands) or isinstance(sequence, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(sequence, list):
        return [sequence, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str


def retranslate(sequence, method, table_id, strands, exception_collector):
    """ Sends sequence for translation and collects exceptions when
        the sequence cannot be translated.

        Args:
            sequence (str): string representing a DNA sequence.
            method (str): a string specifying the sequence orientation
            that should be used to attempt translation.
            table_id (int): translation table identifier.
            strands (list): list with 4 different orientations that can
            be checked.
            exception_collector (list): list used to store all exceptions
            arising from translation attempts.

        Returns:
            A list with following elements, if the sequence can be translated:
                translated_seq (list): a list with the protein sequence and
                with the DNA sequence in the orientation used for translation.
                exception_collector (list): a list with the exceptions that are
                captured when the sequence could not be translated.
            Otherwise:
                translated_seq (str): a string with the exception/reason why
                the sequence could not be translated.
                exception_collector (list): list with all exception that have
                been captured during translation attempts of the current
                sequence.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def translate_coding_sequences(sequences_file, valid_seqs, dna_valid_file, protein_valid_file, table_id):
    """ Translates CDSs into protein sequences.

        Args:
            cdss (dict): a dictionary with CDSs ids as keys and CDSs DNA 
            sequences as values.

        Returns:
            prots (dict): a dictionary with CDSs/proteins ids as keys and protein
            sequences as values.
            trans_state (dict): a dictionary with the CDSs/proteins ids as keys 
            and the DNA strand that coded for those proteins.
            ambiguous (dict): a dictionary with CDSs/proteins ids as keys and 
            CDSs DNA sequences that had ambiguous bases as values.
            untranslatable (dict): a dictionary with CDSs/proteins ids as keys 
            and CDSs DNA sequences that could not be translated.

        Example:

            >>> translate_cds(coding_sequences)
            {1: 'MTPKPVPDKDKYDPTG',
            2: 'MSPLGMIKDEGLFNTELD',
            ...},
            {1: 'sense',
            2: 'antisense',
            ...},
            {3: 'AAATTTNTGCATGA',
            4: 'ATGATTTTTGCBTGA',
            ...}
    """

    dna_lines = []
    prot_lines = []
    invalid_alleles = []
    dna_seqs_file = dna_valid_file
    protein_seqs_file = protein_valid_file
    total_seqs = 0

    cds_index = SeqIO.index(sequences_file, 'fasta')

    for seqid in valid_seqs:
        sequence = str(cds_index.get(seqid).seq)

        translation = translate_dna(sequence, table_id)

        if isinstance(translation, list):
            dna_lines.append('>{0}'.format(seqid))
            dna_lines.append(translation[0][1])
            prot_lines.append('>{0}'.format(seqid))
            prot_lines.append(str(translation[0][0]))
            total_seqs += 1
        # if returned value is a string, translation failed and
        # string contains exceptions
        elif isinstance(translation, str):
            invalid_alleles.append([seqid, translation])

    with open(dna_seqs_file, 'w') as dna_out:
        dna_lines = '\n'.join(dna_lines)
        dna_out.writelines(dna_lines)

    with open(protein_seqs_file, 'w') as protein_out:
        prot_lines = '\n'.join(prot_lines)
        protein_out.writelines(prot_lines)

    return [invalid_alleles, total_seqs]


def determine_repeated(sequences_file):
    """
    """

    seqs_dict = {}
    for entry in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        sequence = str(entry.seq.upper())
        seqid = entry.id

        if sequence not in seqs_dict:
            seqs_dict[sequence] = [seqid]
        elif sequence in seqs_dict:
            seqs_dict[sequence].append(seqid)

    # get sequences with more than one seqid
    repeated_seqs = {seq: seqids for seq, seqids in seqs_dict.items() if len(seqids) > 1}
    repeated_seqs = [seqids[1:] for seq, seqids in repeated_seqs.items()]
    repeated_seqs = list(itertools.chain.from_iterable(repeated_seqs))

    return repeated_seqs


def determine_small(sequences_file, minimum_length):
    """ Find protein sequences that are shorter than desired length.

        Args:
            prots (dict): a dictionary with protein ids as keys and protein 
            sequences as values.
            min_len (int): Proteins with a number of amino acids lower than 
            this value are considered small.

        Returns:
            small_proteins (dict): a dictionary with the ids of small proteins 
            as keys and their amino acid sequence as values.

        Example:

            >>> small_prots(prots, 67)
            {5: 'MTPKPVDKDKYD',
            19: 'MTLNEMVGYVISAHHGMYDFCYCSDDAE',
            ...}
    """

    small_seqs = []
    for entry in SeqIO.parse(sequences_file, 'fasta'):
        # seq object has to be converted to string
        sequence = str(entry.seq.upper())
        seqid = entry.id

        if len(sequence) < minimum_length:
            small_seqs.append(seqid)

    return small_seqs


def get_sequences_by_id(sequences_file, seqids, out_file):
    """
    """
    
    # index FASTA file, much faster and efficient with large files
    sequences_index = SeqIO.index(sequences_file, 'fasta')

    selected = []
    for seqid in seqids:
        identifier = sequences_index[seqid].id
        header = '>{0}'.format(identifier)
        sequence = str(sequences_index[seqid].seq)
        selected.append(header)
        selected.append(sequence)

    with open(out_file, 'w') as out:
        lines = '\n'.join(selected)
        out.write(lines)


def load_cluster_info(cluster_file, out_file):
    """
    """

    clusters = import_clusters(cluster_file)

    # create mapping between representative sequences and subjects in clusters
    clusters = clusters_dicts(clusters)

    with open(out_file, 'wb') as out:
        pickle.dump(clusters, out)


def create_protein_dict(protein_file, out_file):
    """
    """

    protein_lines = {}
    for seq in SeqIO.parse(protein_file, 'fasta'):
        # seq object has to be converted to string
        seq_id = seq.id

        sequence = str(seq.seq.upper())

        protein_lines[seq_id] = sequence

    with open(out_file, 'wb') as out:
        pickle.dump(protein_lines, out)


def create_fasta_lines(sequences, genome_id):
    """

        Args:

        Returns:

        Example:
    """

    lines = []
    for seqid in sequences:
        header = '>' + genome_id + '-protein' + str(seqid)
        sequence = sequences[seqid]

        lines.append(header)
        lines.append(sequence)

    return lines


def write_fasta(fasta_lines, output_file):
    """
    """

    joined_lines = '\n'.join(fasta_lines)

    with open(output_file, 'a') as file:
        file.write(joined_lines+'\n')

    return 'Wrote FASTA sequences to ' + output_file


def write_protein_table(file_name, genome_id, cds_info):
    """
    """

    total_proteins = 0
    with open(file_name, 'a') as file:
        for protein_list in cds_info:
            line = [genome_id] + protein_list
            line = '\t'.join(line)+'\n'
            file.writelines(line)

            total_proteins += len(line)

    return 'Wrote information about ' + str(total_proteins) + ' proteins to ' + file_name


def read_blast_tabular(blast_tabular_file):
    """ Read a file with BLAST results in tabular format

        Args: 
            blast_tabular_file (str): path to output file of BLAST.

        Returns:
            blasting_results (list): a list with a sublist per line in the input
            file.
    """

    with open(blast_tabular_file, 'r') as blastout:
        blasting_results = []
        reader = csv.reader(blastout, delimiter='\t')
        for row in reader:
            blasting_results.append(row)

    return blasting_results


def remove_same_locus_alleles(genes, genes_to_remove, proteinFIlePath, output_file, sizethresh):
    """
    """

    outputFIlePath = output_file
    concatenatedFile = ''

    # Create directory for final protogenome file
    if not proteinFIlePath and outputFIlePath and not os.path.exists(outputFIlePath):
        os.makedirs(outputFIlePath)

    for contig in SeqIO.parse(genes, "fasta", IUPAC.unambiguous_dna):
        name2 = contig.id

        # print name2
        if name2 not in genes_to_remove:
            if int(len(contig.seq)) >= sizethresh:

                concatenatedFile += ">" + contig.id + " \n" + str(contig.seq.upper()) + "\n" 

    if proteinFIlePath and outputFIlePath:
        with open(outputFIlePath, "w") as f:
            f.write(concatenatedFile)


def build_schema(last_file, output_file):
    """
    """

    outputFIlePath = output_file
    if not os.path.exists(output_file):
        os.mkdir(output_file)

    listfiles = []

    total_genes = 0

    for contig in SeqIO.parse(last_file, "fasta", IUPAC.unambiguous_dna):
        namefile = contig.name
        namefile_replaced = new_utils.replace_multiple_characters(namefile)

        if outputFIlePath:
            newFile = os.path.join(outputFIlePath, namefile_replaced + ".fasta")
            listfiles.append(newFile)
            with open(newFile, "w") as f:
                f.write(">" + namefile_replaced + "_1\n" + str(contig.seq).upper() + "\n")

                total_genes += 1

    if outputFIlePath:
        init_schema_4_bbaca.get_Short(listfiles)
        print('\nTotal of {0} loci that constitute the schema.'.format(total_genes))


def cd_hit_blast(inputs):
    """
    """

    blast_inputs = inputs[0:-2]
    blastp_path = inputs[-2]
    output_directory = inputs[-1]

    for cluster in blast_inputs:

        cluster_id = cluster[0]
        cluster_file = os.path.join(output_directory, '{0}.clstr'.format(cluster_id))
        with open(cluster_file, 'rb')as clstr:
            cluster_ids = pickle.load(clstr)

        blast_db = cluster[1]

        ids_str = '\n'.join(cluster_ids)
        ids_file = os.path.join(output_directory, '{0}_ids.txt'.format(cluster_id))
        with open(ids_file, 'w') as ids:
            ids.write(ids_str)

        # import protein sequences
        with open(cluster[2], 'rb') as proteins:
            protein_sequences = pickle.load(proteins)

        fasta_file = os.path.join(output_directory, '{0}_protein.fasta'.format(cluster_id))
        with open(fasta_file, "w") as output:
            seqs_lines = []
            for s in cluster_ids:
                header = '>{0}\n'.format(s)
                sequence = '{0}\n'.format(protein_sequences[s])
                seqs_lines.append(header)
                seqs_lines.append(sequence)

            output.writelines(seqs_lines)

        blast_output = os.path.join(output_directory, '{0}_blast_out.tsv'.format(cluster_id))
        blast_command = ('{0} -db {1} -query {2} -out {3} -outfmt "6 qseqid sseqid score" '
                         '-max_hsps 1 -num_threads {4} -evalue 0.001 -seqidlist {5}'.format(blastp_path, blast_db,
                                                                                            fasta_file, blast_output,
                                                                                            1, ids_file))

        os.system(blast_command)


def prune_clusters(clstr_file, cutoff):
    """
    """

    with open(clstr_file, 'r') as clstrs:
        lines = clstrs.readlines()

    new_lines = []
    excluded_seqids = []
    for l in lines:
        if '>Cluster' in l:
            new_lines.append(l)
        elif '*' in l:
            new_lines.append(l)
        else:
            current_line = l.split(' ')
            percentage = current_line[-1].strip()
            percentage = float(percentage.strip('%'))

            if percentage < cutoff:
                new_lines.append(l)
            else:
                seqid = l.split('...')[0].split('>')[1]
                excluded_seqids.append(seqid)

    os.remove(clstr_file)
    with open(clstr_file, 'w') as new_file:
        new_file.writelines(new_lines)

    return excluded_seqids


def hierarchical_clustering(protein_fasta, output_directory, threads, cd_hit_sim, cd_hit_word, cutoff_sim):
    """
    """

    # cluster file based on defined similarity
    sim_percentage = int(cd_hit_sim*100)
    print('Clustering at {0}%...'.format(sim_percentage))
    clstr_file = 'cd_hit_{0}perc_n{1}'.format(sim_percentage, cd_hit_word)
    clstr_path = '{0}/{1}'.format(output_directory, clstr_file)
    cd_hit_command = ('cd-hit -i {0} -o {1} -c {2} -n {3} '
                             '-d 200 -T {4} -sc 1 -sf 1 -g 1 >/dev/null 2>&1'.format(protein_fasta, clstr_path,
                                                                                     cd_hit_sim, cd_hit_word,
                                                                                     threads))
    os.system(cd_hit_command)

    # read clusters and remove sequences that share a cutoff_sim equal or greater than
    # the defined value. This avoids BLASTing sequences that are from the same gene
    pruned_seqids = prune_clusters(clstr_path, cutoff_sim)

    # sort clusters by cluster size
    sort_clusters_command = ('clstr_sort_by.pl < {0}/{1}.clstr '
                             '> {0}/cd_hit_{2}perc_n{3}_sorted.clstr'.format(output_directory, clstr_file,
                                                                             sim_percentage, cd_hit_word))
    os.system(sort_clusters_command)

    print('Finished clustering sequences!')

    return pruned_seqids


def import_clusters(clstr_file):
    """
    """

    with open(clstr_file, 'r') as file:
        clusters = {}
        # read all lines from clstr file
        lines = file.readlines()
        # remove new line chars
        lines = [line.strip() for line in lines]
        # variable to signal singletons
        singleton = False
        l = 0
        # increment to check the number of sequences per cluster
        # if a new cluster identifier appears when this variable is 1, the last
        # cluster was a singleton
        last_seen = 0
        # clusters are ordered by descending number of sequences
        # when the first singleton is found, stop iterating over the lines
        while singleton != True:
            new_line = lines[l]
            # if the new line is a cluster identifier line
            if 'Cluster' in new_line:
                # check if the last cluster was a singleton
                if last_seen == 1:
                    singleton = True
                    # remove last cluster, it was a singleton
                    clusters.pop(cluster_id)
                # create a new dictionary entry for the new cluster
                else:
                    cluster_id = new_line.split(' ')[-1]
                    clusters[cluster_id] = []
                    # reset variable value
                    last_seen = 0
            # if the new line is not a cluster identifier, keep adding the sequences
            # in the cluster to the dictionary entry
            else:
                clusters[cluster_id].append(new_line)
                last_seen += 1

            # increment variable with line number
            l += 1

            # if we reach the end of the file without finding singletons
            # stop iterating
            if l == len(lines):
                singleton = True

    return clusters


def clusters_dicts(clusters):
    """
    """

    rep_sub_mapping = {}
    # for each cluster (cluster identifiers are integers)
    for i in clusters:
        # get all sequences in that cluster
        cluster_lines = clusters[i]
        # list to store sequences identifiers in the cluster
        all_subs = []
        for line in range(len(cluster_lines)):
            # representative sequences have '*'
            if '*' in cluster_lines[line]:
                # split in order to get only the identifier
                rep_id = cluster_lines[line].split('>')[1]
                rep_id = rep_id.split('...')[0]
                # create dictionary entry with represetative identifier as key
                rep_sub_mapping[rep_id] = []
                # add representative identifier to list
                all_subs.append(rep_id)
            else:
                # add to the cluster list every other sequence
                sub_id = cluster_lines[line].split('>')[1]
                sub_id = sub_id.split('...')[0]
                all_subs.append(sub_id)

        # add the list of sequences identifiers as value
        rep_sub_mapping[rep_id] += all_subs

    return rep_sub_mapping


def blast_inputs(clusters_dict_file, blastdb_path, proteins_dict_file, output_directory):
    """
    """

    # import clusters dict
    with open(clusters_dict_file, 'rb') as infile:
        clusters = pickle.load(infile)

    ids_to_blast = []
    for i in clusters:

        cluster_file = os.path.join(output_directory, '{0}.clstr'.format(i))
        with open(cluster_file, 'wb') as out_clstr:
            pickle.dump(clusters[i], out_clstr)

        blast_input = [i, blastdb_path, proteins_dict_file]
        ids_to_blast.append(blast_input)

    return ids_to_blast


def split_blast_inputs_by_core(blast_inputs, threads, blast_files_dir):
    """
    """

    splitted_ids = [[] for cpu in range(threads)]
    splitted_values = [[] for cpu in range(threads)]
    cluster_sums = [0] * threads
    i = 0
    for cluster in blast_inputs:
        cluster_file = os.path.join(blast_files_dir, '{0}.clstr'.format(cluster[0]))
        with open(cluster_file, 'rb') as infile:
            cluster_seqs = pickle.load(infile)
        splitted_values[i].append(len(cluster_seqs))
        splitted_ids[i].append(cluster)
        cluster_sums[i] += len(cluster_seqs)
        i = cluster_sums.index(min(cluster_sums))

    return splitted_ids


def divide_list_into_n_chunks(list_to_divide, n):
    """
    """

    sublists = []
    d, r = divmod(len(list_to_divide), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list_to_divide[si:si+(d+1 if i < r else d)])

    return sublists


def apply_bsr(inputs):
    """
    """

    blast_file = inputs[0]
    blast_results = read_blast_tabular(blast_file)
    self_scores = {r[0]: r[2] for r in blast_results if r[0] == r[1]}
    blast_results = [r for r in blast_results if r[0] != r[1]]

    fasta_file = inputs[1]
    lengths = {}
    for k in self_scores:
        record = fasta_file.get(k)
        sequence = str(record.seq)
        lengths[k] = len(sequence)

    bsr = inputs[2]

    excluded_alleles = []
    for res in blast_results:

        query = res[0]
        hit = res[1]
        score = res[2]

        if query not in excluded_alleles:
            self_blast_score = self_scores[query]

            query_length = lengths[query]
            hit_length = lengths[hit]
            blast_score_ratio = float(score) / float(self_blast_score)

            # BSR has to be greater than threshold, just as in the original function
            if blast_score_ratio >= bsr and hit not in excluded_alleles:

                if hit_length > query_length and query not in excluded_alleles:
                    excluded_alleles.append(query)

                elif hit_length <= query_length:
                    excluded_alleles.append(hit)

    return excluded_alleles
