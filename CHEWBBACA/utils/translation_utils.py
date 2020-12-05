#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


DESCRIPTION

"""


from Bio.Seq import Seq

try:
    from utils import str_utils as su
except:
    from CHEWBBACA.utils import str_utils as su


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Parameters
        ----------
        dna_str : str
            String representing a DNA sequence.
        table_id : int
            Translation table identifier.

        Returns
        -------
        protseq : Bio.Seq.Seq
            Protein sequence created by translating the
            input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Attempts to translate an input DNA sequence in specified
        orientation and stores exceptions when the input sequence
        cannot be translated.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.
        method : str
            Sequence orientation to attempt translation.
        table_id : int
            Translation table identifier.

        Returns
        -------
        If the sequence can be translated:
            protseq : Bio.Seq.Seq
                Translated DNA sequence.
            myseq : str
                String representing the DNA sequence in the
                orientation used to translate it.
        Otherwise, returns string with the description of the
        exception that was raised.
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
            myseq = su.reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = su.reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = su.reverse_str(myseq)
            myseq = su.reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def translate_dna(dna_sequence, table_id, min_len):
    """ Checks if sequence is valid and attempts to translate
        it, calling several functions to ensure that the sequence
        only has 'ACTG', is multiple of 3 and that it can be
        translated in any of 4 different orientations. Stores
        exceptions so that it is possible to understand why the
        sequence could not be translated.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.
        table_id : int
            Translation table identifier.

        Returns
        -------
        If the sequence can be translated:
            sequence : list
                List with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
            coding_strand : str
                The sequence orientation that codes for the
                protein.
        Otherwise:
            exception_str : str
                A string containing the exceptions that
                explain why the the sequence could not be
                translated.
    """

    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the string is DNA, without ambiguous bases
    valid_dna = su.check_str_alphabet(original_seq, 'ACTG')
    if valid_dna is not True:
        return valid_dna

    # check if sequence is multiple of three
    valid_length = su.check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return valid_length

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

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

        Parameters
        ----------
        sequence : str
            String representing a DNA sequence.
        method : str
            Sequence orientation to attempt translation.
        table_id : int
            Translation table identifier.
        strands : list
            List with maximum of 4 different orientations
            to attempt translation, 'original', 'revcomp',
            'rev' and 'revrevcomp'.
        exception_collector : list
            List used to store all exceptions arising from
            translation attempts.

        Returns
        -------
        If the sequence can be translated:
            translated_seq : list
                List with the protein sequence and with the
                DNA sequence in the orientation used for translation.
            exception_collector : list
                List with the exceptions that were captured when the
                sequence could not be translated.
        Otherwise:
            translated_seq : str
                String with the description of the last exception
                captured because sequence could not be translated.
            exception_collector : list
                List with all exception that have been captured
                for translation attempts.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]
