CreateSchema - Create a gene-by-gene schema
===========================================

What is a Schema and how is it defined
::::::::::::::::::::::::::::::::::::::

A Schema is a pre-defined set of loci that is used in MLST analyses. Traditional MLST schemas
relied in 7 loci that were internal fragments of housekeeping genes and each locus was defined
by its amplification by a pair of primers yielding a fragment of a defined size.

In genomic analyses, schemas are a set of loci that are:

- Present in the majority of strains for core genome (cg) MLST schemas, typically a threshold
  of presence in 95% of the strains is used in schema creation. The assumption is that in each
  strain up to 5% of loci may not be identified due to sequencing coverage problems, assembly
  problems or other issues related to the use of draft genome assemblies.

- Present in at least one of the analyzed strains in the schema creation for pan genome/whole 
  genome (pg/wg) MLST schemas.

- Present in less than 95% of the strains for accessory genome (ag) MLST schemas.

It is important to consider that these definitions are always operational in nature, in the sense
that the analyses are performed on a limited number of strains representing part of the biological
diversity of a given species or genus and are always dependent on the definition of thresholds.  

In most cg/wg/pg/ag MLST schemas, contrary to MLST schemas, each locus corresponds to a coding sequence
(CDS). However, depending on the allele calling algorithm, the alleles called for a given locus can be
CDSs or best matches to existing CDSs without enforcing the need for the identified allele to be a CDS.  

In **chewBBACA**, schemas are composed of loci defined by CDSs and, by default, all the called alleles of a given
locus are CDSs as defined by `Prodigal <https://github.com/hyattpd/Prodigal>`_ (it is also possible to provide
FASTA files with CDSs and the ``--cds`` parameter to skip the gene prediction step with Prodigal).
The use of Prodigal, instead of simply ensuring the presence of start and stop codons, adds an extra layer
of confidence in identifying the most probable CDS for each allele. Because of this approach there may
be variability in the size of the alleles identified by **chewBBACA** and by default a threshold of +/-20%
of the mode of the size of the alleles of a given locus is used to identify a locus as present.

Create a wgMLST schema
::::::::::::::::::::::

Given a set of genome assemblies in FASTA format, chewBBACA offers the option to create a new schema by defining
the distinct loci present in the genomes.

The schema creation algorithm has the following main steps:

- Gene predictipon with Prodigal followed by coding sequence (CDS) extraction to create FASTA files
  that contain all CDSs extracted from the inputs. (there is also the option to provide FASTA files
  with CDSs and the ``--cds`` parameter to skip the gene prediction step with Prodigal).

- Identification of the distinct CDSs (chewBBACA stores information about the distinct CDSs and the
  genomes that contain those CDSs in a hashtable with the mapping between CDS SHA-256 and list of unique
  integer identifiers for the inputs that contain each CDS compressed with `polyline encoding <https://developers.google.com/maps/documentation/utilities/polylinealgorithm>`_
  adapted from `numcompress <https://github.com/amit1rrr/numcompress>`_).

- Exclusion of the CDSs smaller than the value passed to the ``--l`` parameter (default: 201).

- Translation of distinct CDSs that were not an exact match in the previous step (This step identifies
  and excludes CDSs that contain ambiguous bases).

- Protein deduplication to identify the distinct set of proteins and keep information about the inputs that
  contain CDSs that encode each distinct protein (hashtable with mapping between protein SHA-256 and list of
  unique integer identifiers for the distinct CDSs encoded with polyline encoding).

- Minimizer-based clustering. The distinct proteins are sorted in order of decreasing length and
  clustered based on the percentage of shared distinct minimizers (default >= 20%, interior minimizers
  selected based on lexicographic order, k=5, w=5). The first protein is chosen as representative of
  the first cluster and a new cluster is defined each time a protein cannot be added to any of the
  previously defined clusters based on the percentage of minimizers shared with the cluster repsentatives.

- Exclude proteins that share >=90% minimizers with cluster representatives (we assume that these
  sequences represent alleles for the same gene and only keep one representative per gene).

- Exclude proteins that share >=90% minimizers with other proteins in the same cluster (a cluster
  might include sequences from multiple genes and we want to keep only one representative sequence
  per gene).

- Align proteins in each cluster with BLASTp to select a set of representative proteins per cluster
  based on the BLAST Score Ratio (BSR) computed for each alignment.

- Align the selected representatives for all clusters with BLASTp to identify and exclude representative
  proteins that are highly similar (default: BSR >= 0.6) to other representative proteins. The remaining
  set of proteins is not considered highly similar based on the clustering or alignment approach and
  constitutes the schema seed.

- Create the schema seed directory structure with one FASTA file per representative CDS (proteins are converted
  back into DNA). The schema seed can be used to perform allele calling.

.. image::

Basic Usage
-----------

::

	$ chewBBACA.py CreateSchema -i /path/to/InputAssemblies -o /path/to/OutputFolderName --n SchemaName --ptf /path/to/ProdigalTrainingFile --cpu 4

.. important::
	You should adjust the value passed to the ``--cpu`` parameter based on the specifications of
	your machine. chewBBACA will automatically adjust the value if it matches or exceeds the number
	of available CPU cores.

.. important::
	The use of a prodigal training file for schema creation is highly recommended.

Parameters
----------

::

    -i, --input-files           (Required) Path to the directory that contains the input FASTA files.
                                Alternatively, a file with a list of full paths to FASTA files, one
                                per line.

    -o, --output-directory      (Required) Output directory where the process will store intermediate
                                files and create the schema's directory.

    --n, --schema-name          (Optional) Name given to the folder that will store the schema files
                                (default: schema_seed).

    --ptf, --training-file      (Optional) Path to the Prodigal training file. We strongly advise users
                                to provide a Prodigal training file and to keep using the same training
                                file to ensure consistent results (default: None).

    --bsr, --blast-score-ratio  (Optional) BLAST Score Ratio value. Sequences with alignments with a BSR
                                value equal to or greater than this value will be considered as sequences
                                from the same gene (default: 0.6).

    --l, --minimum-length       (Optional) Minimum sequence length value. Coding sequences shorter than
                                this value are excluded (default: 201).

    --t, --translation-table    (Optional) Genetic code used to predict genes and to translate coding
                                sequences (default: 11).

    --st, --size-threshold      (Optional) CDS size variation threshold. Added to the schema's config
                                file and used to identify alleles with a length value that deviates
                                from the locus length mode during the allele calling process (default: 0.2).

    --cpu, --cpu-cores          (Optional) Number of CPU cores that will be used to run the CreateSchema
                                process (will be redefined to a lower value if it is equal to or exceeds
                                the total number of available CPU cores)(default: 1).

    --b, --blast-path           (Optional) Path to the BLAST executables (default: assumes BLAST executables
                                were added to PATH).

    --pm, --prodigal-mode       (Optional) Prodigal running mode ("single" for finished genomes, reasonable
                                quality draft genomes and big viruses. "meta" for metagenomes, low quality
                                draft genomes, small viruses, and small plasmids) (default: single).

    --cds, --cds-input          (Optional) If provided, input is a single or several FASTA files with coding
                                sequences (one per input genome, default: False).
		
    --no-cleanup                (Optional) If provided, intermediate files generated during process execution
                                are not removed at the end (default: False).

.. important::
  If you provide the ``--cds-input`` parameter, chewBBACA assumes that the input FASTA files contain
  coding sequences and skips the gene prediction step with Prodigal. To avoid issues related with the
  format of the sequence headers, chewBBACA renames the sequence headers based on the unique basename
  prefix determined for each input file and on the order of the coding sequences (e.g.: coding sequences
  inside a file named ``GCF_000007125.1_ASM712v1_cds_from_genomic.fna`` are renamed to
  ``GCF_000007125-protein1``, ``GCF_000007125-protein2``, ..., ``GCF_000007125-proteinN``).

Outputs
-------

::

	OutputFolderName
	├── SchemaName
	│   ├── short
	│   │   ├── GenomeID_proteinN_short.fasta
	│   │   ├── ...
	│   │   └── GenomeID_proteinN_short.fasta
	│   ├── GenomeID_proteinN.fasta
	│   ├── ...
	│   ├── GenomeID_proteinN.fasta
	│   └── Training_file.trn
	├── invalid_cds.txt
	└── cds_coordinates.tsv

- One FASTA file per distinct gene identified in the schema creation process in the
  ``OutputFolderName/SchemaName`` directory. The name attributed to each FASTA file in
  the schema is based on the genome of origin of the representative allele chosen for that
  gene and on the order of gene prediction (e.g.: ``GCA-000167715-protein12.fasta``,
  first allele for the gene was identified in a genome assembly with the prefix ``GCA-000167715``
  and the gene was the 12th gene predicted by Prodigal in that assembly).

- The ``OutputFolderName/SchemaName`` directory also contains a directory named ``short`` that
  includes FASTA files with the representative alleles for each locus.

- The training file passed to create the schema is also included in ``OutputFolderName/SchemaName``
  and will be automatically detected during the allele calling process.

- The ``cds_coordinates.tsv`` file contains the coordinates (genome unique identifier, contig
  identifier, start position, stop position, protein identifier attributed by chewBBACA, and coding
  strand (chewBBACA<=3.2.0 assigns 1 to the forward strand and 0 to the reverse strand and
  chewBBACA>=3.3.0 assigns 1 and -1 to the forward and reverse strands, respectively)) of the CDSs
  identified in each genome. 

- The ``invalid_cds.txt`` file contains the list of alleles predicted by Prodigal that were
  excluded based on the minimum sequence size value and presence of ambiguous bases.
