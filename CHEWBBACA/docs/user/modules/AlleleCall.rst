AlleleCall -  Determine the allelic profiles of a set of genomes
================================================================

What is an allele?
::::::::::::::::::

In Biology, an allele is a specific sequence variant that occurs at a given locus. However, given a DNA sequence, the assignment of a putative allele to a locus is influenced by several factors, such as:

- Quality of the sequence assembly (influenced by several aspects, such as the sequencing method, the assembler used, etc);
- If the alleles must correspond to Coding DNA Sequences (CDSs) and open reading frames (ORFs);
- Presence of possibly homologous loci (this situation can result in an allele assignment to a possibly wrong locus given the difficulty in distinguishing closely related homologs).

Therefore, in gene-by-gene methods, the definition of an allele is determined by the sequence similarity search method and all the parameters used to decide if an allele can be identified as a *bona fide* allele of a given locus.

In chewBBACA, by default, an allele needs to be a CDS defined by `Prodigal <https://github.com/hyattpd/Prodigal>`_ (for chewBBACA <= 3.2.0) or `Pyrodigal <https://github.com/althonos/pyrodigal>`_ (for chewBBACA >= 3.3.0). To ensure reproducibility of the CDS prediction, the same Prodigal training file for each bacterial species should be used and provided as input. Users can also provide input files with CDSs, in which case the gene prediction step with Prodigal or Pyrodigal will be skipped.

.. important::
  The `Prodigal wiki <https://github.com/hyattpd/prodigal/wiki/Gene-Prediction-Modes#training-mode>`_ includes information about how to create a training file. chewBBACA >=3.3.0 uses `Pyrodigal <https://github.com/althonos/pyrodigal>`_ for gene prediction. Pyrodigal is a Python module that provides bindings to Prodigal, including several bug fixes and performance improvements.

Perform allele calling
::::::::::::::::::::::

Having defined a cgMLST or wgMLST schema, we can proceed to use it to call alleles on target genomes. chewBBACA's allele calling algorithm can also use schemas from existing databases such as `Ridom cgMLST <http://www.cgmlst.org/ncs>`_, `BIGSdb <https://pubmlst.org/>`_, `BIGSdb-Pasteur <https://bigsdb.pasteur.fr/>`_, `Enterobase <http://enterobase.warwick.ac.uk/>`_ and `Chewie-NS <https://chewbbaca.online/>`_. External schemas can be adapted for usage with chewBBACA through the :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module.

.. warning::
  If you installed chewBBACA v3 and want to use a schema created with chewBBACA v2, please use the :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module to convert the schema to a format fully compatible with chewBBACA v3.

.. warning::
  If you want to run chewBBACA in an environment with multiple processes/users accessing the same schema, please read the section about concurrent allele calling at the end of this page.

.. important::
  Although the use of a training file is optional, it is highly recommended to ensure consistent results.

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

    --output-masked             (Optional) Create a TSV file with the masked allelic profiles. The masking
                                process removes the ``INF-` prefix from inferred alleles and substitutes all
                                special classes (NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5, LOTSC, PAMA) with `0` (default: False).

    --no-cleanup                (Optional) If provided, intermediate files generated during process execution are
                                not removed at the end (default: False).

    --hash-profiles             (Optional) Create a TSV file with hashed allelic profiles. Profiles can be hashed
                                using any of the hash algorithms implemented in the hashlib and zlib Python libraries
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
	By default, the *AlleleCall* module uses the Prodigal training file included in the schema's directory and it is not necessary to pass a training file to the ``--ptf`` parameter.

.. important::
  If you provide the ``--cds-input`` parameter, chewBBACA assumes that the input FASTA files contain CDSs and skips the gene prediction step. To avoid issues related to the format of the sequence headers, chewBBACA renames the sequence headers based on the unique basename prefix determined for each input file and on the order of the CDSs (e.g.: CDSs inside a file named ``GCF_000007125.1_ASM712v1_cds_from_genomic.fna`` are renamed to ``GCF_000007125-protein1``, ``GCF_000007125-protein2``, ..., ``GCF_000007125-proteinN``).

.. note::
  If a text file that contains a list of full paths to loci FASTA files or loci IDs, one per line, is passed to the ``--gl`` parameter, the process will only perform allele calling for the loci in that list.

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


- The ``cds_coordinates.tsv`` file contains the coordinates (genome unique identifier, contig identifier, start position, stop position, protein identifier attributed by chewBBACA, and coding strand (chewBBACA<=3.2.0 assigns 1 to the forward strand and 0 to the reverse strand and chewBBACA>=3.3.0 assigns 1 and -1 to the forward and reverse strands, respectively)) of the CDSs identified in each genome. 

- The ``invalid_cds.txt`` file contains the list of alleles predicted by Prodigal that were excluded based on the minimum sequence size value and presence of ambiguous bases.

- The ``loci_summary_stats.tsv`` file contains the classification type counts (*EXC*, *INF*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM*, *ALM*, *ASM*, *PAMA*, *LNF*) and the total number of classified CDSs (non-*LNF*) per locus.

- The ``results_statistics.tsv`` file contains the classification type counts (*EXC*, *INF*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM*, *ALM*, *ASM*, *PAMA*, *LNF*), the total number of invalid CDSs, the total number of classified CDSs (non-*LNF*) and the total number of predicted CDSs per genome.

	The column headers stand for:

	- **EXC** - EXaCt match (100% DNA identity) with previously identified alleles.
	- **INF** - INFerred new alleles that had no exact match in the schema but are highly similar to loci in the schema. In the allelic profiles (included in the ``results_alleles.tsv`` output file), the *INF-* prefix in the allele identifier indicates that the allele was newly inferred and the number following the prefix is the allele identifier attributed to the allele. Inferred alleles are added to the FASTA file of the locus they share high similarity with (unless the ``--no-inferred`` parameter is used).
	- **LNF** - Locus Not Found. The number of schema loci for which it was not possible to identify a similar allele in the input genome. This means that, for those loci, there were no BLAST hits or they were not within the BSR threshold for allele assignment.
	- **PLNF** - Probable Locus Not Found. Attributed when a locus is not found during execution modes 1, 2 and 3. Those modes do not perform the complete analysis, which is only performed in mode 4 (default), and the distinct classification indicates that a more thorough analysis may have found a match for the loci that were not found.
	- **PLOT3/PLOT5** - Possible Locus On the Tip of a contig (see image below). A CDS is classified as *PLOT5* or *PLOT3* if it is close to the contig 5’- or 3’-end and if the unaligned portion of the matched representative allele exceeds the contig end. This could be an artifact caused by genome fragmentation resulting in the prediction of a shorter CDS. To avoid locus misclassification, loci in such situations are classified as *PLOT*.
	- **LOTSC** - A CDS is classified as *LOTSC* if the matched representative allele is bigger than the contig containing the CDS.

|

	.. image:: /_static/images/PLOT5_PLOT3_LOTSC.png
		:width: 700px
		:align: center

|

	- **NIPH/NIPHEM** - The *NIPH* and *NIPHEM* classifications are assigned when multiple CDSs from the same genome match the same schema locus. *NIPH* (Non-Informative Paralogous Hit) is assigned when multiple CDSs from the same genome match a single locus. *NIPHEM* (Non-Informative Paralogous Hit Exact Match) is assigned when multiple CDSs from the same genome are exact matches to alleles of a single locus. These classifications suggest that a locus can have paralogous (or orthologous) loci in the query genome and should be removed from the analysis due to the potential uncertainty in allele assignment (for example, due to the presence of multiple copies of the same mobile genetic element (MGE) or as a consequence of gene duplication followed by pseudogenization). A high number of these classifications may also indicate a poorly assembled genome due to a high number of smaller contigs which result in partial CDS predictions.

|

	.. image:: /_static/images/NIPH_NIPHEM.png
		:width: 700px
		:align: center

|

	- **PAMA** - PAralogous MAtch. Attributed to CDSs that are highly similar to more than one locus. This type of classification allows the identification of groups of similar loci in the schema that are classified as paralogous loci and listed in the ``paralogous_counts.tsv`` and ``paralogous_loci.tsv`` output files.

|

	.. image:: /_static/images/PAMA.png
		:width: 700px
		:align: center

|

	- **ASM/ALM** - The *ASM* (Allele Smaller than Mode) and *ALM* (Allele Larger than Mode) classifications are assigned when the size of a CDS that matches a schema locus is below or above the locus size variation interval, respectively. The default behaviour is to assign these classifications to alleles that are 20% shorter or longer than the locus allele size mode. It is important to remember that, although infrequently, the mode may change as more alleles for a given locus are called and added to a schema. The *ALM* and *ASM* classifications impose a limit on allele size variation since for the majority of loci the allele lengths are quite conserved. However, some loci can have larger variation in allele length and those should be manually curated.

|

	.. image:: /_static/images/ASM_ALM.png
		:width: 700px
		:align: center

|

.. important::
	A high number of *PLOT3/PLOT5*, *LOTSC*, *NIPH/NIPHEM*, *PAMA*, and *ASM/ALM* classifications assigned to a genome usually indicates that the genome is of poor quality (e.g. highly fragmented or contaminated).

- The ``results_contigsInfo.tsv`` file contains the loci coordinates in the genomes analyzed. The first column contains the identifiers of the input genomes and the other columns (with loci names in the headers) the locus coordinate information or the classification attributed by chewBBACA if it was not an exact match or inferred allele.

	.. rst-class:: align-center

		+--------------+---------------------------+-------------------------+-----+
		| FILE         | locus1                    | locus2                  | ... |
		+==============+===========================+=========================+=====+
		| SAMD00008628 | contig2&162560-161414&-1  |             LNF         | ... |
		+--------------+---------------------------+-------------------------+-----+
		| SAMD00053744 | contig4&268254-269400&1   | contig3&272738-274082&1 | ... |
		+--------------+---------------------------+-------------------------+-----+

	- Example for the ``SAMD00008628`` genome:

		- locus1 with ``contig2&161414-162560&0`` information was found in this genome. It is located in (``&`` character is the field delimiter):

			- the sequence with identifier ``contig2``.
			- between 161,414 bp and 162,560 bp (reported as ``162560-161414`` because the CDS is encoded in the reverse strand). These nucleotide positions are inclusive positions and include the stop codon as well.
			- in the reverse strand (represented by a ``-1`` signal). ``1`` means that the CDS is encoded in the forward strand.

		- locus2 was not found (*LNF*).

- The ``results_alleles.tsv`` file contains the allelic profiles determined for the input samples. The first column has the identifiers of the genome assemblies. The remaining columns contain the allele identifiers or other classification labels for loci present in the schema, with the column headers being the locus identifiers.

.. rst-class:: align-center

	+--------------+--------+--------+--------+--------+--------+-----+
	| FILE         | locus1 | locus2 | locus3 | locus4 | locus5 | ... |
	+==============+========+========+========+========+========+=====+
	| SAMD00008628 | INF-2  | 1      | 3      | ASM    | PLOT3  | ... |
	+--------------+--------+--------+--------+--------+--------+-----+
	| SAMD00053744 | 10     | 1      | 3      | ALM    | PLOT5  | ... |
	+--------------+--------+--------+--------+--------+--------+-----+

.. note::
	The *ExtractCgMLST* module was designed to determine the set of loci that constitute the core genome based on a given threshold, but it can also be used to convert the TSV file with allelic profiles into a suitable format that can be imported into `PHYLOViZ <http://www.phyloviz.net/>`_ by substituting all non-numeric classifications by ``0``. To convert an allelic profile output simply run the *ExtractCgMLST* module with a threshold value, ``--t``, of ``0``.

- The ``paralogous_counts.tsv`` file contains the list of paralogous loci and the number of times those loci matched a CDS that was also similar to other loci in the schema.

- The ``paralogous_loci.tsv`` file contains the sets of paralogous loci identified per genome (genome identifier, identifiers of the paralogous loci and the coordinates of the CDS that is similar to the group of paralogous loci).

- The ``logging_info.txt`` contains summary information about the allele calling process.

- If the ``--output-unclassified`` parameter is provided, the process will create a FASTA file, ``unclassified_sequences.fasta``, with the DNA sequences of the distinct CDSs that were not classified.

- If the ``--output-missing`` parameter is provided, the process will create a FASTA file, ``missing_classes.fasta``, and a TSV file with information about the classified sequences that led to a locus being classified as *ASM*, *ALM*, *PLOT3*, *PLOT5*, *LOTSC*, *NIPH*, *NIPHEM* and *PAMA*.

- If the ``--output-novel`` parameter is provided, the process will create a FASTA file, ``novel_alleles.fasta``, with the DNA sequences of the novel alleles inferred during allele calling.

- If the ``--output-masked`` parameter is provided, the process will create a TSV file, ``results_alleles_masked.tsv``, with the masked allelic profiles. The masking process removes the ``INF-`` prefix from inferred alleles and substitutes all special classes (NIPH, NIPHEM, ASM, ALM, PLOT3, PLOT5, LOTSC, PAMA) with ``0``.

- If the ``--hash-profiles`` parameter is provided, the process will use the provided hash algorithm to create a TSV file, ``results_alleles_hashed.tsv``, with hashed profiles (each allele identifier is substituted by the hash of its DNA sequence). chewBBACA can use any of the hash algorithms implemented in the `zlib <https://docs.python.org/3/library/zlib.html>`_ and `hashlib <https://docs.python.org/3/library/hashlib.html>`_ libraries. Please check the documentation of those libraries to know which hash algorithms are available and select one (e.g., if using the `crc32 <https://docs.python.org/3/library/zlib.html#zlib.crc32>`_ algorithm implemented in zlib, include ``--hash-profiles crc32`` in your command).

Workflow of the AlleleCall module
:::::::::::::::::::::::::::::::::

.. image:: /_static/images/AlleleCall.png
   :width: 1000px
   :align: center

The AlleleCall module determines the allelic profiles for strains of interest. Brief description of the workflow:

	- The process accepts FASTA files with genome assemblies or CDSs (``--cds`` parameter). If genome assemblies are given, the process starts by predicting CDSs for each genome using Pyrodigal.

	- The CDSs identified in the input files are deduplicated (chewBBACA stores information about the distinct CDSs and the genomes that contain those CDSs in a hashtable that maps CDS SHA-256 hashes to the list of unique integer identifiers for the input genomes that contain each CDS compressed with `polyline encoding <https://developers.google.com/maps/documentation/utilities/polylinealgorithm>`_ adapted from `numcompress <https://github.com/amit1rrr/numcompress>`_) and compared against the schema alleles to find and classify exact matches at the DNA level (information about the matches found for each locus is saved to classification files, one per locus in the schema. The classification files are updated throughout the process with information about the matches and classifications at each step). If the process runs in mode 1, the results are evaluated to write the output files and exit.

	- Otherwise, the CDSs that do not match any schema alleles at the DNA level are translated (CDS translation identifies and excludes CDSs that contain ambiguous bases and with length below the theshold defined by the ``--l`` parameter) and matched against the translated schema alleles to find exact matches at the protein level (hashtable mapping distinct protein SHA-256 hashes to the list of unique integer identifiers for the distinct CDSs encoded with polyline encoding). If the process runs in mode 2, the results are evaluated to write the output files, add new alleles to the schema and exit.

	- Otherwise, the distinct translated CDSs not classified through exact matching are compared against the schema representative alleles through minimizer-based clustering (minimizers selected based on lexicographic order, k=5, w=5) to identify CDSs that share a proportion of minimizers ≥ 0.2 with the representative alleles.

	- Each cluster's representative allele is aligned against the clustered CDSs with BLASTp to classify CDSs based on the defined BLAST Score Ratio (BSR) value (the default is 0.6) plus 0.1. At this point, if the process runs in mode 3, the results are evaluated to write the output files, add new alleles to the schema and exit.

	- Otherwise, the representative alleles are aligned against the remaining unclassified CDSs to classify them based on the defined BSR value and identify new representative alleles whose BSR is not above the defined BSR value plus 0.1. If the process finds new representative alleles, it aligns them against the unclassified CDSs to find new matches. This process repeats until no new representative alleles are identified.

	- When no new representative alleles are found, the process evaluates the results to create the output files, add new alleles to the schema, and exit

Identify genetic clusters
:::::::::::::::::::::::::

We recommend that you use `ReporTree <https://github.com/insapathogenomics/ReporTree>`_ to identify genetic clusters based on the allelic profiles (contained in the ``results_alleles.tsv`` output file) determined by chewBBACA. ReporTree includes functionalities to identify genetic clusters at any distance threshold level(s), obtain summary reports with relevant statistics computed based on sample metadata, identify regions of cluster stability, etc. Cluster nomenclature can be maintained and updated in subsequent analyses, which is especially useful in surveillance-oriented workflows. Check the `publication <https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-023-01196-1>`_ and the `GitHub repository <https://github.com/insapathogenomics/ReporTree>`_ to know more about ReporTree.

Concurrent allele calling
:::::::::::::::::::::::::

In its default mode, mode 4, and in the execution modes 2 and 3, the AlleleCall module updates the schema with the novel alleles inferred during the allele calling. This is incompatible with concurrent access to the same schema. If you run chewBBACA in an environment with multiple processes/users accessing the same schema, please use the ``--no-inferred`` parameter. By providing this parameter, chewBBACA will still identify novel alleles but will not update the schema files with the information about those novel alleles. **When you create a new schema, adapt an external schema or download a schema from Chewie-NS, you must perform a single allele calling before using the schema for concurrent allele calling**. You can use a single genome assembly; it's only essential to generate the pre-computed data that chewBBACA uses to speed up the allele calling. After that, multiple users can concurrently perform allele calling based on the same schema if they pass the ``--no-inferred`` parameter. chewBBACA will still identify novel alleles and include them in the final results, but those alleles will not be added to the schema, and the pre-computed files will not be updated. If you ever want to add new alleles to the schema, you'll have to perform allele calling without the ``--no-inferred`` parameter and ensure that there's only one process working with the schema while it is updated.

.. warning::
	The schema will most likely become corrupted and unusable if you attempt to run multiple concurrent processes with the same schema without providing the ``--no-inferred`` parameter.
