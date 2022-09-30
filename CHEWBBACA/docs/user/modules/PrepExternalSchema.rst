PrepExternalSchema - Adapt an external schema to be used with chewBBACA
=======================================================================

The PrepExternalSchema module enables the adaptation of external schemas so that it is possible
to use those schemas with chewBBACA. An external schema may be a set of sequences from any number
of genes that have been selected for a particular study or it may be a schema that has already
been defined and is available for download from some well known databases, such as:

- `Ridom cgMLST <http://www.cgmlst.org/ncs>`_.
- `BIGSdb <https://pubmlst.org/>`_.
- `Enterobase <http://enterobase.warwick.ac.uk/>`_.

External schemas are defined with specific parameters that might differ from the parameters and
conditions enforced by chewBBACA. Therefore, these external schemas need to be processed to
filter out sequences that do not meet a set of criteria applied to create every chewBBACA schema.
Every sequence that is included in the final schema has to represent a complete coding sequence
(the first and last codons must be valid start and stop codons, the sequence length must be a
multiple of 3 and cannot contain in-frame stop codons) and contain no invalid or ambiguous
characters (sequences must be composed of ATCG only).

By default, the process will adapt the external schema based on a BLAST Score Ratio (BSR) value of
``0.6``, it will accept sequences of any length and will use the genetic code ``11`` (Bacteria and
Archaea) to translate sequences. These options can be changed by passing different values to
the ``--bsr``, ``--l`` and ``--t`` arguments. The process runs relatively fast with the default value
for the ``--cpu`` argument, but it will complete considerably faster if it can use several CPU cores
to evaluate several loci in parallel.

For each gene in the external schema, and assuming the default BSR value, the process will:

- Exclude sequences with invalid or ambiguous characters.
- Exclude sequences with length value that is not a multiple of 3.
- Try to translate sequences and exclude sequences that cannot be translated in any possible
  orientation due to invalid start and/or stop codons or in-frame stop codons.
- Select the longest (or one of the longest) sequence as the first representative for that gene;
- Use BLASTp to align the representative against all sequences that were not excluded.
- Calculate the BSR value for each alignment.
- If all BSR values are greater than 0.7, the current representative is considered appropriate
  to capture the gene sequence diversity when performing allele calling.
- Otherwise, an additional representative has to be chosen in order to find a suitable set of
  representatives for the gene. The new representative will be the longest sequence from the
  set of non-representative sequences that had a BSR value in the interval [0.6,0.7] (in this
  BSR value interval, aligned sequences are still considered to be alleles of the same gene but
  display a degree of dissimilarity that can contribute to an increase of the sensitivity
  compared to the utilization of only one of those sequences as representative). If there is
  no alignment with a BSR value in the interval [0.6,0.7], the next representative will be the
  longest (or one of the longest) sequence from the set of sequences that had an alignment with
  a BSR<0.6.
- The process will keep expanding the set of representatives until we have a set of
  representatives that when aligned against all alleles of the gene, guarantee that each allele
  has at least one alignment with a BSR>0.7.

After determining the representative sequences, the process writes the FASTA file with all valid
sequences to the adapted schema directory and the FASTA file with only the representatives to
the *short* directory inside the adapted schema directory.

Basic Usage
-----------

To adapt an external schema, the FASTA files with the allele sequences for each gene of the
schema have to be in a single directory. Alternatively, a file containing the list of paths
to the schema files may also be provided.

::

	chewBBACA.py PrepExternalSchema -i /path/to/ExternalSchema -o /path/to/OutputFolderName --ptf /path/to/ProdigalTrainingFile --cpu 4

Parameters
----------

::

    -i, --input-files          (Required) Path to the folder containing the FASTA files, one FASTA file
                               per gene/locus (alternatively, a file with a list of paths can be given).

    -o, --output-directory     (Required) The directory where the output files will be saved (will
                               create the directory if it does not exist).

    --ptf, --training-file     (Optional) Path to the Prodigal training file that will be included in
                               the adapted schema (default: None).

    --bsr, --blast-score-ratio (Optional) The BLAST Score Ratio value that will be used to adapt
                               the external schema (default: 0.6).

    --l, --minimum-length      (Optional) Minimum sequence length accepted. Sequences with a length
                               value smaller than the value passed to this argument will be discarded
                               (default: 0).

    --t, --translation-table   (Optional) Genetic code to use for CDS translation. Must match the
                               genetic code used to create the training file (default: 11).

    --st, --size-threshold     (Optional) CDS size variation threshold. At the default value of
                               0.2, alleles with size variation +-20 percent when compared to the
                               representative will not be included in the final schema (default:
                               0.2).

    --cpu, --cpu-cores         (Optional) The number of CPU cores to use (default: 1).

    --b, --blast-path          (Optional) Path to the BLAST executables. Use this option if
                               chewBBACA cannot find the BLASTp and makeblastdb executables or if
                               you want to use anoter BLAST installation that is not the one added
                               to the PATH (default: assumes BLAST executables were added to PATH).

Outputs
-------

- ``<adapted_schema>_invalid_alleles.txt`` - contains the identifiers of the alleles that were
  excluded and the reason for the exclusion of each allele.
- ``<adapted_schema>_invalid_genes.txt`` - contains the list of genes that had no valid alleles.
- ``<adapted_schema>_summary_stats.tsv`` - contains summary statistics for each gene. Number of
  alleles in the external schema, number of valid alleles included in the adapted schema and
  number of representatives.

.. note::
	For most genes, only one or a few sequences need to be chosen as representatives to
	represent the gene sequence diversity. Nevertheless, some genes will have a high number
	of representatives. This is more common for small genes, where a small number of
	differences has a big impact on the alignment score, for genes with repetitive or low
	complexity regions that are masked by BLAST and lead to lower alignment scores between
	highly similar sequences, and for genes that have inversions, deletions or insertions
	that can lead to several High-scoring Segment Pairs (HSPs), none of which have a score
	sufficiently high to identify both sequences as belonging to the same gene.
