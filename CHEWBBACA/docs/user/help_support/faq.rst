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
The :doc:`TestGenomeQuality </user/modules/TestGenomeQuality>` process can be used to
identify genome assemblies that are responsible for a considerable loss of loci. You can pass
a list of genomes to remove to the :doc:`ExtractCgMLST </user/modules/ExtractCgMLST>` process to
exclude those genomes from the analysis that determines the set of loci that constitute the
core-genome. Identifying and removing low quality genome assemblies can significantly improve
the determination of the core-genome.

Can I use a schema from an external source?
...........................................
Yes. The :doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` process enables the adaptation
of external schemas so that it is possible to use those schemas with chewBBACA. An external
schema may be a set of sequences from any number of genes that have been selected for a particular
study or it may be a schema that has already been defined and is available for download from
some well known databases, such as `Ridom cgMLST <http://www.cgmlst.org/ncs>`_,
`BIGSdb <https://pubmlst.org/>`_ and `Enterobase <http://enterobase.warwick.ac.uk/>`_.

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
