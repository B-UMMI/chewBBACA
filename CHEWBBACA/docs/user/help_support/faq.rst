FAQ
===

What strains should I use for :doc:`Schema Creation </user/modules/CreateSchema>`?
..................................................................................
The set of genome assemblies used for schema creation should be carefully selected to avoid
the inclusion of spurious loci in the schema seed that result from low quality assemblies
(e.g.: genome assemblies resulting from low quality sequencing data, highly fragmented genome
assemblies, genome assemblies with many frameshifted proteins, genome length too large or too
small). A set of high quality genome assemblies, ideally complete genomes, that capture the
diversity of the species or lineage of interest should result in a good schema seed.

Should I worry about file naming best practices for input files?
................................................................
Yes, definitely. Ideally, input FASTA files should have short names without blank spaces or
special characters (e.g.: ``!@#?$^*()+``). It is also important to ensure each input file has a
unique basename prefix (everything before the first ``.`` in the basename). chewBBACA uses the
basename prefix as a unique identifier during process execution and to label the final results.
chewBBACA will warn you and exit if any files share the same basename prefix (e.g.: 
``GCF_002209245.1_ASM220924v1_genomic.fna`` and ``GCF_002209245.2_ASM220924v2_genomic.fna``).
chewBBACA accepts input files with the following file extensions: ``.fasta``, ``.fna``, ``.ffn``,
``.fa`` and ``.fas``.

I ran all the steps and my cgMLST loci size is smaller than traditional MLST...does this even work?
...................................................................................................
In order to have a robust definition of a cgMLST schema for a given bacterial species, a set
of representative strains of the diversity of a given species should be selected. Furthermore,
since cgMLST schema definition is based on pre-defined thresholds, only when a sufficient number
of strains have been analyzed can the cgMLST schema be considered stable. This number will always
depend on the population structure and diversity of the species in question. cgMLST schemas are
defined as the set of loci that are present in all strains under analysis, but defining a smaller
loci presence threshold, such as 95%, might be necessary to include very frequent genes or
ubiquitous genes that are not present in some strains due to sequencing/assembly limitations.
The quality of the genome assemblies is also an important factor that can impact profoundly
the MLST schema definition (e.g.: genome assemblies with a high number of missing loci,
contaminated genome assemblies, misclassified genome assemblies).
The results of the :doc:`ExtractCgMLST </user/modules/ExtractCgMLST>` module can help in
identifying genome assemblies that are responsible for a considerable loss of loci. You can
provide a file with the list of genome assemblies that you've identified as being of low quality
to exclude them from the cgMLST analysis. Identifying and removing low quality genome assemblies
can significantly improve the determination of the core-genome.

Can I use a schema from an external source?
...........................................
Yes. The :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module enables the adaptation
of external schemas so that it is possible to use those schemas with chewBBACA. An external
schema may be a set of sequences from any number of genes that have been selected for a particular
study or it may be a schema that has already been defined and is available for download from
some well known databases, such as `Ridom cgMLST <http://www.cgmlst.org/ncs>`_,
`BIGSdb <https://pubmlst.org/>`_, `BIGSdb-Pasteur <https://bigsdb.pasteur.fr/>`_ and `Enterobase <http://enterobase.warwick.ac.uk/>`_.
You can also download schemas compatible with chewBBACA from `Chewie-NS <https://chewbbaca.online/>`_.
Chewie-NS allows chewBBACA users to download and update cg/wgMLST schemas, allowing the easy sharing of
results, while ensuring the reproducibility and consistency of these steps.

Is it possible to perform allele calling without adding novel alleles to the schema?
....................................................................................
Yes, you can use the ``--no-inferred`` parameter to specify that you do not want to add novel alleles
identified during allele calling. This can be useful to avoid contaminating the schema when you
are not sure if the genome assemblies are from the same taxon as the schema was created for or if they might be
of low quality (e.g. contaminated, highly fragmented, high number of ambiguous bases, etc.).

I would like to inspect the sequences that were not classified during allele calling. Is there an output with information about those sequences?
................................................................................................................................................
chewBBACA can create a FASTA file, ``unclassified_sequences.fasta``, with the coding sequences that were not classified if you pass the ``--output-unclassified``
parameter.

Is there any information about the sequences classified as non-*EXC* or non-*INF*?
..............................................................................
chewBBACA can create a FASTA file, ``missing_classes.fasta``, and a TSV file, ``missing_classes.tsv``, with the list of coding sequences classified as
*ASM/ALM/PLOT3/PLOT5/LOTSC/NIPH/NIPHEM* if you pass the ``--output-missing`` parameter.

Can I get the list of new alleles identified during allele calling?
...................................................................
The ``results_alleles.tsv`` file includes the ``INF-`` prefix in the identifiers of newly inferred alleles. chewBBACA can also create a FASTA file, ``novel_alleles.fasta``,
with the sequences of all new alleles if you pass the ``--output-novel`` parameter.

Are the schemas created with chewBBACA v2 compatible with chewBBACA v3?
.......................................................................
To use a schema created with chewBBACA v2 in chewBBACA v3, you need to use the
:doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module to convert the schema to a format
fully compatible with chewBBACA v3. The adaptation process removes the files terminating in ``bsr.txt``
from the ``short`` directory and reformats the sequence headers.

Which species already have a training file?
...........................................
At the moment:

- *Acinetobacter baumannii*
- *Campylobacter jejuni*
- *Enterococcus faecium*
- *Escherichia coli*
- *Haemophilus influenzae*
- *Legionella pneumophila*
- *Listeria monocytogenes*
- *Salmonella enterica enteritidis*
- *Staphylococcus aureus*
- *Staphylococcus haemolyticus*
- *Streptococcus agalactiae*
- *Streptococcus canis*
- *Streptococcus dysgalactiae*
- *Streptococcus equi*
- *Streptococcus pneumoniae*
- *Streptococcus pyogenes*
- *Yersinia enterocolitica*

get them `here <https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files>`_.

My favorite species has no training file. What can I do?
........................................................
You can propose a new one to be added to the repository or create your own training files.
To create a training file make sure you have Prodigal installed and run the following command:

::

	prodigal -i myGoldStandardGenome.fna -t myTrainedFile.trn -p single

How should I cite chewBBACA?
............................
If you use chewBBACA, please cite:

Silva M, Machado MP, Silva DN, Rossi M, Moran-Gilad J, Santos S, Ramirez M, Carriço JA. 2018. chewBBACA: A complete suite for gene-by-gene schema creation and strain identification. Microb Genom 4:000166. doi:10.1099/mgen.0.000166
