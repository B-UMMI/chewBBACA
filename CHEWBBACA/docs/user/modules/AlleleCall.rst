AlleleCall -  Determine the allelic profiles of a set of genomes
================================================================

What is an allele?
::::::::::::::::::

In Biology, an allele is a specific sequence variant that occurs at a given locus.
However, given a DNA sequence, the assignment of a putative allele to a locus can be
confounded by several factors:

- Quality of the sequence assembly (influenced by several aspects, such as the sequencing
  method, the assembler used, etc);
- If the alleles must correspond to coding sequences (CDSs);
- Presence of possibly homologous loci (this situation can result in a wrong allele assignment
  to a given locus given the difficulty in distinguishing closely related homologs).

Therefore, in gene-by-gene methods, the definition of an allele is determined by the sequence
similarity search method and all the parameters used to decide if an allele can be identified
as a *bona fide* allele of a given locus.

In chewBBACA, an allele needs to be a CDS defined by `Prodigal <https://github.com/hyattpd/Prodigal>`_.
To ensure reproducibility of the CDS prediction, the same Prodigal training file for
each bacterial species should be used and provided as input. 

.. important::
	Please read the `Prodigal wiki <https://github.com/hyattpd/prodigal/wiki>`_ for more
	information about the requirements to create a training file.

.. image:: http://i.imgur.com/H9pKjHQ.png
	:width: 1400px
	:align: center

A 100% DNA identity comparison is performed first on all the genome CDSs against
each locus allele database. If an exact match is found an allele identification is attributed
to the CDS. The next step translates the remaining unclassified CDSs and searches for
exact matches at protein level against the schema. Exact matches are classified as new inferred alleles (INF).
If not a BLAST Score Ratio `(BSR) <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2>`_
approach is used to identify the allele. To improve computational efficiency, chewBBACA
processes each locus in the schema separately for homology search, parallelising the jobs
using the specified number of CPUs. Therefore, for each locus, the short list containing the
more divergent alleles of each locus is queried against the BLASTp database. The BSR is
calculated for each hit and based on these results and a size validation step, the locus
is either considered not found (*LNF*) or a new allele for the locus is inferred. After
running each genome, the loci database is updated with the new alleles found and, whenever
required, a locus short list is also updated with a new divergent allele.
There are also other options of allele classification (*NIPH*, *ASM*, *ALM* and *PLOT*)
that are described later on in this document. All thresholds for considering or excluding
alleles can be altered. For BSR we considered a value of 0.6 which is related to a DNA
identity of 70%-80%. For size exclusion (*ASM* and *ALM* classifications), the allele
length mode (of the updated database) +/- 20% considered as default was based in empirical
observations of manually curated allele variation for several species.

Perform allele calling
::::::::::::::::::::::

Having defined a cgMLST or wgMLST schema, we can proceed to use it to call alleles on target
genomes. chewBBACA's *Allele Calling* algorithm can also use schemas and alleles from existing
databases such as `PubMLST <https://pubmlst.org/>`_, `EnteroBase <https://enterobase.warwick.ac.uk/>`_
or `RIDOM(tm) cgMLST schemas <http://www.cgmlst.org/ncs>`_. External schemas can be adapted for
usage with chewBBACA through the :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module.

.. important::
	Although the use of a training file is optional, it is highly recommended to ensure consistent
	results.

Basic Usage
-----------

::

	chewBBACA.py AlleleCall -i /path/to/InputAssemblies -g /path/to/SchemaName -o /path/to/OutputFolderName --cpu 4

Parameters
----------

::

    -i, --input-files           (Required) Path to the directory with the genome FASTA files or to a file with
                                a list of paths to the FASTA files, one per line.

    -g, --schema-directory      (Required) Path to the schema directory with the loci FASTA files.  

    -o, --output-directory      (Required) Output directory where the allele calling results will be stored
                                (will create a subdirectory named "results_\<TIMESTAMP\>" if the path passed
                                by the user already exists).

    --ptf, --training-file      (Optional) Path to the Prodigal training file. Default is to get training
                                file from the schema's directory (default: searches for a training file in
                                the schema's directory).

    --gl, --genes-list          (Optional) Path to a file with the list of genes in the schema that the process
                                should identify alleles for (default: False).

    --bsr, --blast-score-ratio  (Optional) BLAST Score Ratio value. Sequences with alignments with a BSR
                                value equal to or greater than this value will be considered as sequences
                                from the same gene (default: uses value defined in the schema config file).

    --l, --minimum-length       (Optional) Minimum sequence length accepted for a coding sequence to be included
                                in the schema (default: uses value defined in schema config file. Default value
                                added to the config file is 0).

    --t, --translation-table    (Optional) Genetic code used to predict genes and to translate coding sequences.
                                Must match the genetic code used to create the training file (default: uses value
                                defined in schema config).

    --st, --size-threshold      (Optional) CDS size variation threshold. If set to a value of 0.2, alleles with
                                size variation +-20 percent will be classified as ASM/ALM (default: uses value
                                defined in schema config).

    --cpu, --cpu-cores          (Optional) Number of CPU cores that will be used to run the AlleleCall process
                                (will be redefined to a lower value if it is equal to or exceeds the total number
                                of available CPU cores/threads)(default: 1).

    --b, --blast-path           (Optional) Path to the BLAST executables. Use this option if chewBBACA cannot find
                                the BLASTp and makeblastdb executables or if you want to use anoter BLAST installation
                                that is not the one added to the PATH (default: assumes BLAST executables were added
                                to PATH).

    --pm, --prodigal-mode       (Optional) Prodigal running mode (default: single).

    --cds, --cds-input          (Optional) Input files contain coding sequences (one Fasta file per strain). Skips
                                gene prediction with Prodigal (default: False).

    --no-inferred               (Optional) If provided, the process will not add the sequences of inferred alleles
                                to the schema (default: False).

    --output-unclassified       (Optional) Create a Fasta file with unclassified coding sequences (default: False).

    --output-missing            (Optional) Create a Fasta file with coding sequences classified as NIPH, NIPHEM,
                                ASM, ALM, PLOT3, PLOT5 and LOTSC (default: False).

    --no-cleanup                (Optional) If provided, intermediate files generated during process execution are
                                not removed at the end (default: False).

    --hash-profile              (Optional) Create TSV file with hashed allelic profiles. Profiles can be hashed
                                with any of the hash algorithms implemented in the hashlib and zlib libraries
                                (default: None).

    --force-continue            (Optional) If provided, chewBBACA will add config files with default parameter
                                values to schemas that are missing those files and will also proceed if any of
                                the argument values does not match the value in the config files. Otherwise, it
                                will prompt users for the parameter values to add to the config files and for
                                permission to proceed if the argument values differ from the ones in the config
                                files (default: False).

    --mode                      (Optional) Execution mode (1: only exact matches at DNA level; 2: exact matches
                                at DNA and Protein level; 3: exact matches and minimizer-based clustering to find
                                similar alleles based on BSR+0.1; 4: runs the full process to find exact matches
                                and similar matches based on BSR value) (default: 4).

.. important::
	By default, the *AlleleCall* process uses the Prodigal training file included in the schema's
	directory and it is not necessary to pass a training file to the ``--ptf`` parameter.

.. note::
	If a text file with a list of gene identifiers, one per line, is passed to the ``--gl``
	parameter, the process will only perform allele calling for the genes in that list.

Outputs
-------

::

	OutputFolderName
	├── cds_coordinates.tsv
	├── invalid_cds.txt
	├── loci_summary_stats.tsv
	├── results_statistics.tsv
	├── results_contigsInfo.tsv
	├── results_alleles.tsv
	├── paralogous_counts.tsv
	├── paralogous_loci.tsv
	└── logging_info.txt


- The ``cds_coordinates.tsv`` file contains the coordinates (genome unique identifier, contig
  identifier, start position, stop position, protein identifier attributed by chewBBACA and coding
  strand) of the coding sequences identified in each genome.

- The ``invalid_cds.txt`` file contains the list of alleles predicted by Prodigal that were
  excluded based on the minimum sequence size value and presence of ambiguous bases.

- The ``loci_summary_stats.tsv`` file contains the counts for each classification type (*EXC*,
  *INF*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM*, *ALM*, *ASM*, *LNF*) and the total number
  of classified CDS (non-*LNF*) per locus.

- The ``results_statistics.tsv`` file contains the total number of exact matches (*EXC*), inferred
  new alleles (*INF*), loci on contig tips (*PLOT3*/*PLOT5*), loci identified on contigs smaller than
  the matched schema representative (*LOTSC*), non-informative paralogous hits (*NIPH*/*NIPHEM*),
  alleles larger than locus length mode (*ALM*), alleles smaller than locus length mode (*ASM*)
  and loci not found (*LNF*) classifications attributed for each genome.

+--------------+-----+------+-------+-------+-------+------+--------+-----+-----+-----+
| FILE         | EXC | INF  | PLOT3 | PLOT5 | LOTSC | NIPH | NIPHEM | ALM | ASM | LNF |
+==============+=====+======+=======+=======+=======+======+========+=====+=====+=====+
| SAMD00008628 | 14  | 1722 | 0     | 0     | 0     |    8 |      0 |   1 |   2 |   1 |
+--------------+-----+------+-------+-------+-------+------+--------+-----+-----+-----+
| SAMD00053744 | 600 | 1138 | 0     | 0     | 0     | 4    | 4      | 1   | 1   | 0   |
+--------------+-----+------+-------+-------+-------+------+--------+-----+-----+-----+

The column headers stand for:

- *EXC* - alleles which have exact matches (100% DNA identity) with previously identified
  alleles.
- *INF* - inferred new alleles that had no exact match in the schema but are highly
  similar to loci in the schema. The *INF-* prefix in the allele identifier indicates that
  such allele was newly inferred in that genome, and the number following the prefix is the
  allele identifier attributed to such allele. Inferred alleles are added to the FASTA file of the locus they
  share high similarity with.
- *LNF* - loci not found. No alleles were found for the number of loci in the schema shown.
  This means that, for those loci, there were no BLAST hits or they were not within the BSR
  threshold for allele assignment.
- *PLOT3/PLOT5* - possible loci on the tip of the query genome contigs (see image below). A locus
  is classified as *PLOT* when the CDS of the query genome has a BLAST hit with a known larger
  allele that covers the CDS sequence entirely and the unaligned regions of the larger allele
  exceed one of the query genome contigs ends (a locus can be classified as *PLOT5* or *PLOT3*
  depending on whether the CDS in the genome under analysis matching the schema locus is located
  in the 5' end or 3' end (respectively) of the contig). This could be an artifact caused by
  genome fragmentation resulting in a shorter CDS prediction by Prodigal. To avoid locus
  misclassification, loci in such situations are classified as *PLOT*.

.. image:: http://i.imgur.com/41oONeS.png
	:width: 700px
	:align: center

- *LOTSC* - A locus is classified as *LOTSC* when the contig of the query genome is smaller
  than the matched allele.
- *NIPH* - non-informative paralogous hit (see image below). When ≥2 CDSs in the query
  genome match one locus in the schema with a BSR > 0.6, that locus is classified as *NIPH*.
  This suggests that such locus can have paralogous (or orthologous) loci in the query genome
  and should be removed from the analysis due to the potential uncertainty in allele assignment
  (for example, due to the presence of multiple copies of the same mobile genetic element (MGE)
  or as a consequence of gene duplication followed by pseudogenization). A high number of *NIPH*
  may also indicate a poorly assembled genome due to a high number of smaller contigs which
  result in partial CDS predictions. These partial CDSs may contain conserved domains that
  match multiple loci.
- *NIPHEM* - similar to the *NIPH* classification, but specifically
  referring to exact matches. Whenever several CDSs from the same genome match a single or
  multiple alleles of the same locus with 100% DNA similarity during the first DNA sequence
  comparison, the *NIPHEM* tag is attributed.

.. image:: http://i.imgur.com/4VQtejr.png
	:width: 700px
	:align: center

- *ALM* - alleles 20% larger than the length mode of the distribution of the matched
  loci (CDS length > (locus length mode + locus length mode * 0.2)) (see image below).
  This determination is based on the currently identified set of alleles for a given locus.
  It is important to remember that, although infrequently, the mode may change as more
  alleles for a given locus are called and added to a schema.
- *ASM* - similar to *ALM* but for alleles 20% smaller than the length mode distribution
  of the matched loci (CDS length < (locus length mode - locus length mode * 0.2)). As with
  *ALMs* it is important to remember that, although infrequently, the mode may change as
  more alleles for a given locus are called and added to a schema.

.. image:: http://i.imgur.com/l1MDyEz.png
	:width: 700px
	:align: center

.. note::
	The *ALM* and *ASM* classifications impose a limit on size variation since for the
	majority of loci the allele lengths are quite conserved. However, some loci can have larger
	variation in allele length and those should be manually curated.

The statistics file also helps the user to identify bad quality draft genomes among the
analyzed genomes since with a proper schema most identified loci should be exact matches
or inferred alleles. A high number of *PLOT*, *ASM*, *ALM* and/or *NIPH* usually indicates
bad quality or contaminated assemblies.

- The ``results_contigsInfo.tsv`` file contains the loci coordinates in the genomes analyzed. The
  first column contains the identifier of the genome used in the allele calling and the other
  columns (with loci names in the headers) the locus coordinate information or the classification
  attributed by chewBBACA if it was not an exact match or inferred allele.

+--------------+--------------------------+-------------------------+-----+
| FILE         | locus1                   | locus2                  | ... |
+==============+==========================+=========================+=====+
| SAMD00008628 | contig2&162560-161414&0  |             LNF         | ... |
+--------------+--------------------------+-------------------------+-----+
| SAMD00053744 | contig4&268254-269400&1  | contig3&272738-274082&1 | ... |
+--------------+--------------------------+-------------------------+-----+

Example for the ``SAMD00008628`` genome:

	- locus1 with ``contig2&161414-162560&0`` information was found in this genome. It is located
	  in (``&`` character is the field delimiter):

	    - the sequence with identifier ``contig2``.
	    - between 161,414 bp and 162,560 bp (reported as ``162560-161414`` because the CDS is encoded
	      in the reverse strand). These nucleotide positions are inclusive positions and include the
	      stop codon as well.
	    - in the reverse strand (represented by a ``0`` signal). ``1`` means that the CDS is encoded
	      in the direct strand.

	- locus2 was not found (*LNF*).

- The ``results_alleles.tsv`` file contains the allelic profiles determined for the input samples.
  The first column has the identifiers of the genome assemblies for which the allele call was
  performed. The remaining columns contain the allele call data for loci present in the schema,
  with the column headers being the locus identifiers.

+--------------+--------+--------+--------+--------+--------+-----+
| FILE         | locus1 | locus2 | locus3 | locus4 | locus5 | ... |
+==============+========+========+========+========+========+=====+
| SAMD00008628 | INF-2  | 1      | 3      | ASM    | PLOT3  | ... |
+--------------+--------+--------+--------+--------+--------+-----+
| SAMD00053744 | 10     | 1      | 3      | ALM    | PLOT5  | ... |
+--------------+--------+--------+--------+--------+--------+-----+

.. note::
	The allelic profile output can be transformed and imported into
	`PHYLOViZ <http://www.phyloviz.net/>`_ to generate and visualize a Minimum Spanning
	Tree.

.. important::
	The *ExtractCgMLST* module was designed to determine the set of loci that
	constitute the core genome based on a given threshold, but it can also be used to
	convert the TSV file with allelic profiles into a suitable format that can be imported
	into PHYLOViZ. To convert an allelic profile output simply run the *ExtractCgMLST* module
	with a threshold value, ``--t``, of ``0``.

- The ``paralogous_counts.tsv`` file contains the list of paralogous loci and the number of times
  those loci matched a CDS that was also similar to other loci in the schema.

- The ``paralogous_loci.tsv`` file contains the sets of paralogous loci identified per genome
  (genome identifier, identifiers of the paralogous loci and the coordinates of the CDS that
  is similar to the group of paralogous loci).

.. image:: http://i.imgur.com/guExrGx.png
	:width: 700px
	:align: center

- The ``logging_info.txt`` contains summary information about the allele calling process.

- If the ``--output-unclassified`` parameter is provided, the process will create a FASTA file
  with the DNA sequences of the distinct CDSs that were not classified.

- If the ``--output-missing`` parameter is provided, the process will create a FASTA file with
  the DNA sequences of the CDSs classified as *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM*, *ALM*
  and *ASM*.

- If the ``--hash-profiles`` parameter is provided, the process will use the provided hash
  algorithm to create a TSV file with hashed profiles (each allele identifier is substituted
  by the hash of the DNA sequence).
