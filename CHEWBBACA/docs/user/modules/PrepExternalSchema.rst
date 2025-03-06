PrepExternalSchema - Adapt an external schema to be used with chewBBACA
=======================================================================

The PrepExternalSchema module enables the adaptation of external schemas so that it is possible
to use those schemas with chewBBACA. An external schema may be a set of sequences from any number
of genes that have been selected for a particular study or it may be a schema that has already
been defined and is available for download from some well known databases, such as:

- `Ridom cgMLST <http://www.cgmlst.org/ncs>`_.
- `BIGSdb <https://pubmlst.org/>`_.
- `BIGSdb-Pasteur <https://bigsdb.pasteur.fr/>`_.
- `Enterobase <http://enterobase.warwick.ac.uk/>`_.

External schemas are defined with specific parameters that might differ from the parameters and
conditions enforced by chewBBACA. Therefore, these external schemas need to be processed to
filter out sequences that do not meet a set of criteria applied to create every chewBBACA schema.
Every sequence that is included in the final schema has to represent a complete coding sequence
(the first and last codons must be valid start and stop codons, the sequence length must be a
multiple of 3 and cannot contain in-frame stop codons) and contain no invalid or ambiguous
characters (sequences must be composed of ATCG only).

Basic Usage
-----------

To adapt an external schema, the FASTA files with the allele sequences for each gene of the
schema have to be in a single directory. Alternatively, a file containing the list of paths
to the schema files may also be provided.

::

	chewBBACA.py PrepExternalSchema -g /path/to/ExternalSchemaFolder -o /path/to/OutputFolder --ptf /path/to/ProdigalTrainingFile --cpu 4

Parameters
----------

::

    -g, --schema-directory     (Required) Path to the directory of the schema to adapt. The schema must contain
                               one FASTA file per gene/locus.

    -o, --output-directory     (Required) Path to the output directory where the adapted schema will be created.

    --gl, --genes-list         (Optional) Path to a file with the list of loci in the schema that the process
                               should adapt (one per line, full paths or loci IDs) (default: False).

    --ptf, --training-file     (Optional) Path to the Prodigal training file that will be included in the
                               directory of the adapted schema (default: None).

    --bsr, --blast-score-ratio (Optional) BLAST Score Ratio (BSR) value. The process selects representative
                               alleles for each locus based on this value. Representative alleles are
                               selected until all alleles in a locus align against one of the representatives
                               with a BSR >= than the specified value (default: 0.6).

    --l, --minimum-length      (Optional) Minimum sequence length value stored in the schema config file. The
                               schema adaptation process will only discard sequences smaller than this value
                               if the --size-filter parameter is provided (default: 0).

    --t, --translation-table   (Optional) Genetic code used for allele translation (default: 11).

    --st, --size-threshold     (Optional) Allele size variation threshold value stored in the schema config file.
                               The schema adaptation process will only discard alleles with a size that deviates
                               from the locus length mode +- the size theshold value if the --size-filter parameter
                               is provided(default: 0.2).

    --cpu, --cpu-cores         (Optional) TNumber of CPU cores/threads that will be used to run the process
                               (chewie resets to a lower value if it is equal to or exceeds the total number
                               of available CPU cores/threads) (default: 1).

    --b, --blast-path          (Optional) Path to the directory that contains the BLAST executables
                               (default: assumes BLAST executables were added to PATH).

    --size-filter              (Optional) Apply the minimum length and size threshold values to
                               filter out alleles during schema adaptation (default: False).

Outputs
-------

- ``<adapted_schema>_invalid_alleles.txt`` - contains the identifiers of the alleles that were
  excluded and the reason for the exclusion of each allele.
- ``<adapted_schema>_invalid_genes.txt`` - contains the list of genes that had no valid alleles, one gene identifier per line.
- ``<adapted_schema>_summary_stats.tsv`` - contains summary statistics for each gene (number of
  alleles in the external schema, number of valid alleles included in the adapted schema and
  number of representative alleles chosen by chewBBACA).

.. note::
	For most genes, only one or a few sequences need to be chosen as representatives to
	represent the gene sequence diversity. Nevertheless, some genes will have a high number
	of representatives. This is more common for small genes, where a small number of
	differences has a big impact on the alignment score, for genes with repetitive or low
	complexity regions that might be masked by BLAST and lead to lower alignment scores between
	highly similar sequences, and for genes that have inversions, deletions or insertions
	that can lead to several High-scoring Segment Pairs (HSPs), none of which have a score
	sufficiently high to identify both sequences as belonging to the same gene.

Workflow of the PrepExternalSchema module
:::::::::::::::::::::::::::::::::::::::::

.. image:: /_static/images/PrepExternalSchema.png
   :width: 1200px
   :align: center

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
