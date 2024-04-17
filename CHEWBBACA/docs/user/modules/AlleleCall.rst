AlleleCall -  Determine the allelic profiles of a set of genomes
================================================================

What is an allele?
::::::::::::::::::

In Biology, an allele is a specific sequence variant that occurs at a given locus.
However, given a DNA sequence, the assignment of a putative allele to a locus is
influenced by several factors:

- Quality of the sequence assembly (influenced by several aspects, such as the sequencing
  method, the assembler used, etc);
- If the alleles must correspond to Coding DNA Sequences (CDSs) and open reading frames (ORFs);
- Presence of possibly homologous loci (this situation can result in an allele assignment
  to a possibly wrong locus given the difficulty in distinguishing closely related homologs).

Therefore, in gene-by-gene methods, the definition of an allele is determined by the sequence
similarity search method and all the parameters used to decide if an allele can be identified
as a *bona fide* allele of a given locus.

In chewBBACA, by default, an allele needs to be a CDS defined by `Prodigal <https://github.com/hyattpd/Prodigal>`_
(for chewBBACA <= 3.2.0) or `Pyrodigal <https://github.com/althonos/pyrodigal>`_ (for chewBBACA >= 3.3.0).
To ensure reproducibility of the CDS prediction, the same Prodigal training file for each bacterial species should
be used and provided as input. Users can also provide input files with CDSs, in which case the gene prediction step
with Prodigal or Pyrodigal will be skipped.

.. important::
  chewBBACA >=3.3.0 uses `Pyrodigal <https://github.com/althonos/pyrodigal>`_ for gene
  prediction. Pyrodigal is a Python module that provides bindings to Prodigal, including
  several bug fixes and performance improvements.

.. important::
  Please read the `Prodigal wiki <https://github.com/hyattpd/prodigal/wiki>`_ for more
  information about the requirements to create a training file.

The allele calling algorithm has the following main steps:

- Gene predictipon with Prodigal followed by CDS extraction to create FASTA files
  that contain all CDSs extracted from the inputs (There is also the option to provide FASTA files
  with CDSs and the ``--cds`` parameter to skip the gene prediction step with Prodigal).

- Identification of the distinct CDSs (chewBBACA stores information about the distinct CDSs and the
  genomes that contain those CDSs in a hashtable with the mapping between CDS SHA-256 and list of unique
  integer identifiers for the inputs that contain each CDS compressed with `polyline encoding <https://developers.google.com/maps/documentation/utilities/polylinealgorithm>`_
  adapted from `numcompress <https://github.com/amit1rrr/numcompress>`_).

- Exact matches at DNA level between the alleles of each locus in the schema and the CDSs identified
  in the input files (Information about the exact matches found for each locus is saved to
  classification files, one per locus in the schema. The classification files are updated throughout
  the process with information about the matches and classifications at each step. If ``--mode`` equals 1,
  skips the next classification steps and goes directly to the last step).

- Translation of distinct CDSs that were not an exact match in the previous step (This step identifies
  and excludes CDSs that contain ambiguous bases and with length below the theshold defined by the ``--l``
  parameter).

- Protein deduplication to identify the distinct set of proteins and keep information about the inputs that
  contain CDSs that encode each distinct protein (Hashtable with mapping between protein SHA-256 and list of
  unique integer identifiers for the distinct CDSs encoded with polyline encoding).

- Exact matches at protein level between the translated alleles of each locus in the schema and the
  translated CDSs identified in the input files (Classification files are updated. If ``--mode`` equals 2,
  skips the next classification steps and goes directly to the last step).

- Minimizer-based clustering. The distinct proteins are sorted in order of decreasing length and
  clustered based on the percentage of shared distinct minimizers (Default >= 20%, interior minimizers
  selected based on lexicographic order, k=5, w=5) with the schema representative alleles.

- Alignment, using BLASTp, of each cluster representative against all proteins added to its cluster to
  identify and classify matches with a BLAST Score Ratio (BSR) >0.7 (Classification files are updated.
  If ``--mode`` equals 3, skips the next step).

- Align the schema representatives against the distinct proteins that have not been classified. Compute the
  BSR for each alignment and classify proteins with alignments with BSR >=0.6. For each locus, determine new
  representative alleles from the set of proteins that had a BSR between 0.6 and 0.7. Realign the new
  representatives against the remaining unclassified proteins. Keep iterating until no proteins are classified
  for any locus (classification files are updated in each iteration).

- Assign novel allele identifiers, add novel alleles to the schema and write the output files based on the
  information stored in the classification files. If ``--no-inferred`` is passed, this step is bypassed and
  no novel alleles are added to the database. This results in keeping a stable database to make allele calls
  against. When in doubt about the quality of the genomes one is performing allele call on, this parameter
  will prevent entering into the database alleles that may result from poor quality data rather than actual
  alleles present in the bacterial population.

Perform allele calling
::::::::::::::::::::::

Having defined a cgMLST or wgMLST schema, we can proceed to use it to call alleles on target
genomes. chewBBACA's allele calling algorithm can also use schemas and alleles from existing
databases such as `Ridom cgMLST <http://www.cgmlst.org/ncs>`_, `BIGSdb <https://pubmlst.org/>`_,
`BIGSdb-Pasteur <https://bigsdb.pasteur.fr/>`_, `Enterobase <http://enterobase.warwick.ac.uk/>`_ and
`Chewie-NS <https://chewbbaca.online/>`_. External schemas can be adapted for
usage with chewBBACA through the :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module.

.. warning::
  If you installed chewBBACA v3 and want to use a schema created with chewBBACA v2, please use the
  :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module to convert the schema to a format
  fully compatible with chewBBACA v3.

.. warning::
  If you want to run chewBBACA in an environment with multiple processes/users accessing the same schema,
  please read the section about concurrent allele calling at the end of this page.

.. important::
  Although the use of a training file is optional, it is highly recommended to ensure consistent
  results.

Basic Usage
-----------

::

	chewBBACA.py AlleleCall -i /path/to/InputAssembliesFolder -g /path/to/SchemaFolder -o /path/to/OutputFolder --cpu 4

Parameters
----------

::

    -i, --input-files           (Required) Path to the directory that contains the input FASTA files or to a file
                                with a list of full paths to FASTA files, one per line.

    -g, --schema-directory      (Required) Path to the schema directory. The schema directory contains the loci
                                FASTA files and a folder named "short" that contains the FASTA files with the
                                loci representative alleles.

    -o, --output-directory      (Required) Output directory where the process will store intermediate files and
                                allele calling results (will create a subdirectory named "results_<TIMESTAMP>"
                                if the path passed by the user already exists).

    --ptf, --training-file      (Optional) Path to the Prodigal training file used by Pyrodigal to predict genes.
                                Default is to use the training file included in the schema's directory (default: None)

    --gl, --genes-list          (Optional) Path to a file with the list of genes/loci to perform allele calling.
                                The file must include the full paths to the loci FASTA files or the loci IDs, one
                                per line. The process will perform allele calling only for the subset of genes
                                provided in the file (default: False).

    --bsr, --blast-score-ratio  (Optional) BLAST Score Ratio (BSR) value. The BSR is computed for each BLASTp
                                alignment and aligned sequences with a BSR >= than the defined value are
                                considered to be alleles of the same gene (default: uses value defined in the
                                schema config file).

    --l, --minimum-length       (Optional) Minimum sequence length value. Predicted coding sequences (CDSs)
                                shorter than this value are excluded (default value added to the config file is 0).

    --t, --translation-table    (Optional) Genetic code used to predict genes and to translate coding sequences.
                                Must match the genetic code used to create the training file (default: uses value
                                defined in schema config).

    --st, --size-threshold      (Optional) Coding sequence (CDS) size variation threshold. At the default value of
                                0.2, CDSs with a size that deviates +-20 percent from the locus length mode are
                                classified as ASM/ALM (default: uses value defined in schema config).

    --cpu, --cpu-cores          (Optional) Number of CPU cores that will be used to run the process (chewie resets
                                to a lower value if it is equal to or exceeds the total number of available CPU cores)
                                (default: 1).

    --b, --blast-path           (Optional) Path to the directory that contains the BLAST executables. (default: assumes
                                BLAST executables were added to PATH).

    --pm, --prodigal-mode       (Optional) Prodigal running mode ("single" for finished genomes, reasonable
                                quality draft genomes and big viruses. "meta" for metagenomes, low quality
                                draft genomes, small viruses, and small plasmids) (default: single).

    --cds, --cds-input          (Optional) If provided, chewBBACA skips the gene prediction step and assumes the
                                input FASTA files contain coding sequences (one FASTA file per strain) (default: False).

    --no-inferred               (Optional) If provided, the process will not add the sequences of inferred alleles
                                (INF) to the schema. Allelic profiles will still include the allele identifiers
                                attributed to the inferred alleles. Use this parameter if the schema is being
                                accessed by multiple processes/users simultaneously (default: False).

    --output-unclassified       (Optional) Create a Fasta file with the coding sequences (CDSs) that were not
                                classified (default: False).

    --output-missing            (Optional) Create a Fasta file with coding sequences classified as NIPH, NIPHEM,
                                ASM, ALM, PLOT3, PLOT5 and LOTSC (default: False).

    --output-novel              (Optional) Create Fasta file with the novel alleles inferred during allele calling.
                                The sequence headers include the locus and allele identifiers attributed by
                                chewBBACA based on the allele calling results (default: False).

    --no-cleanup                (Optional) If provided, intermediate files generated during process execution are
                                not removed at the end (default: False).

    --hash-profile              (Optional) Create a TSV file with hashed allelic profiles. Profiles can be hashed
                                with any of the hashing algorithms implemented in the hashlib and zlib Python libraries
                                (default: None).

    --force-continue            (Optional) If provided, chewie will not warn users and ask for permission to
                                continue if any of the provided argument values does not match the values in the
                                config file (default: False).

    --mode                      (Optional) Execution mode (1: only exact matches at DNA level; 2: exact matches
                                at DNA and Protein level; 3: exact matches and minimizer-based clustering to find
                                similar alleles based on BSR+0.1; 4: runs the full process to find exact matches
                                and similar matches based on BSR value, including the determination of new
                                representative alleles to add to the schema) (default: 4).

.. important::
	By default, the *AlleleCall* module uses the Prodigal training file included in the schema's
	directory and it is not necessary to pass a training file to the ``--ptf`` parameter.

.. important::
  If you provide the ``--cds-input`` parameter, chewBBACA assumes that the input FASTA files contain
  CDSs and skips the gene prediction step. To avoid issues related to the
  format of the sequence headers, chewBBACA renames the sequence headers based on the unique basename
  prefix determined for each input file and on the order of the CDSs (e.g.: CDSs
  inside a file named ``GCF_000007125.1_ASM712v1_cds_from_genomic.fna`` are renamed to
  ``GCF_000007125-protein1``, ``GCF_000007125-protein2``, ..., ``GCF_000007125-proteinN``).

.. note::
  If a text file that contains a list of full paths to loci FASTA files or loci IDs, one per line, is
  passed to the ``--gl`` parameter, the process will only perform allele calling for the loci in that list.

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
  identifier, start position, stop position, protein identifier attributed by chewBBACA, and coding
  strand (chewBBACA<=3.2.0 assigns 1 to the forward strand and 0 to the reverse strand and
  chewBBACA>=3.3.0 assigns 1 and -1 to the forward and reverse strands, respectively)) of the CDSs
  identified in each genome. 

- The ``invalid_cds.txt`` file contains the list of alleles predicted by Prodigal that were
  excluded based on the minimum sequence size value and presence of ambiguous bases.

- The ``loci_summary_stats.tsv`` file contains the classification type counts (*EXC*,
  *INF*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM*, *ALM*, *ASM*, *PAMA*, *LNF*) and the total number
  of classified CDSs (non-*LNF*) per locus.

- The ``results_statistics.tsv`` file contains the classification type counts (*EXC*,
  *INF*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM*, *ALM*, *ASM*, *PAMA*, *LNF*), the total number
  of invalid CDSs, the total number of classified CDSs (non-*LNF*) and the total number of predicted
  CDSs per genome.

The column headers stand for:

- *EXC* - EXaCt matches (100% DNA identity) with previously identified alleles.
- *INF* - INFerred new alleles that had no exact match in the schema but are highly
  similar to loci in the schema. The *INF-* prefix in the allele identifier indicates that
  such allele was newly inferred in that genome, and the number following the prefix is the
  allele identifier attributed to such allele. Inferred alleles are added to the FASTA file of the locus they
  share high similarity with.
- *LNF* - Locus Not Found. No alleles were found for the number of loci in the schema shown.
  This means that, for those loci, there were no BLAST hits or they were not within the BSR
  threshold for allele assignment.
- *PLNF* - Probable Locus Not Found. Attributed when a locus is not found during execution modes 1, 2 and 3.
  Those modes do not perform the complete analysis, that is only performed in mode 4 (default), and the
  distinct classification indicates that a more thorough analysis might have found a match for the loci
  that were not found.
- *PLOT3/PLOT5* - Possible Locus On the Tip of the query genome contigs (see image below). A locus
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
- *NIPH* - Non-Informative Paralogous Hit (see image below). When ≥2 CDSs in the query
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
- *PAMA* - PAralogous MAtch. Attributed to CDSs that are highly similar to more than one locus.
  This type of classification allows the identification of groups of similar loci in the
  schema that are classified as paralogous loci and listed in the ``paralogous_counts.tsv`` and
  ``paralogous_loci.tsv`` files.

.. image:: http://i.imgur.com/4VQtejr.png
	:width: 700px
	:align: center

- *ALM* - Alleles 20% Larger than the length Mode of the distribution of the matched
  loci (CDS length > (locus length mode + locus length mode * 0.2)) (see image below).
  This determination is based on the currently identified set of alleles for a given locus.
  It is important to remember that, although infrequently, the mode may change as more
  alleles for a given locus are called and added to a schema.
- *ASM* - similar to *ALM* but for Alleles 20% Smaller than the length Mode distribution
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
	into PHYLOViZ by substituting all non-numeric classifications by ``0``. To convert an
	allelic profile output simply run the *ExtractCgMLST* module with a threshold
	value, ``--t``, of ``0``.

- The ``paralogous_counts.tsv`` file contains the list of paralogous loci and the number of times
  those loci matched a CDS that was also similar to other loci in the schema.

- The ``paralogous_loci.tsv`` file contains the sets of paralogous loci identified per genome
  (genome identifier, identifiers of the paralogous loci and the coordinates of the CDS that
  is similar to the group of paralogous loci).

.. image:: http://i.imgur.com/guExrGx.png
	:width: 700px
	:align: center

- The ``logging_info.txt`` contains summary information about the allele calling process.

- If the ``--output-unclassified`` parameter is provided, the process will create a FASTA file, ``unclassified_sequences.fasta``,
  with the DNA sequences of the distinct CDSs that were not classified.

- If the ``--output-missing`` parameter is provided, the process will create a FASTA file, ``missing_classes.fasta``, and a
  TSV file with information about the classified sequences that led to a locus being classified
  as *ASM*, *ALM*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM* and *PAMA*.

- If the ``--hash-profiles`` parameter is provided, the process will use the provided hash
  algorithm to create a TSV file, ``results_alleles_hashed.tsv``, with hashed profiles (each allele identifier is substituted
  by the hash of the DNA sequence).

Identify genetic clusters
:::::::::::::::::::::::::

We recommend that you use `ReporTree <https://github.com/insapathogenomics/ReporTree>`_ to identify genetic clusters
based on the allelic profiles (contained in the ``results_alleles.tsv`` output file) determined by chewBBACA. ReporTree
includes functionalities to identify genetic clusters at any distance threshold level(s), obtain summary reports with relevant
statistics computed based on sample metadata, identify regions of cluster stability, etc. Cluster nomenclature can be maintained
and updated in subsequent analyses, which is especially useful in surveillance-oriented workflows. Check the
`publication <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1>`_ and the
`GitHub repository <https://github.com/insapathogenomics/ReporTree>`_ to know more about ReporTree.

Concurrent allele calling
:::::::::::::::::::::::::

In its default mode, mode 4, and in the execution modes 2 and 3, the AlleleCall module updates the schema with the
novel alleles inferred during the allele calling. This is incompatible with concurrent access to the same schema.
If you run chewBBACA in an environment with multiple processes/users accessing the same schema, please use the
``--no-inferred`` parameter. By providing this parameter, chewBBACA will still identify novel alleles but will not
update the schema files with the information about those novel alleles. When you create a new schema, adapt an external
schema or download a schema from Chewie-NS, you must perform a single allele calling before using the schema for
concurrent allele calling. You can use a single genome assembly; it's only essential to generate the pre-computed data
that chewBBACA uses to speed up the allele calling. After that, multiple users can concurrently perform allele calling
based on the same schema if they pass the ``--no-inferred`` parameter. chewBBACA will still identify novel alleles and
include them in the final results, but those alleles will not be added to the schema, and the pre-computed files will
not be updated. If you ever want to add new alleles to the schema, you'll have to perform allele calling without the
``--no-inferred`` parameter and ensure that there's only one process working with the schema while it is updated.

.. warning::
	The schema will most likely become corrupted and unusable if you attempt to run multiple concurrent processes
	with the same schema without providing the ``--no-inferred`` parameter.
