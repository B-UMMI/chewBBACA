SchemaEvaluator - Build a html report to better navigate/visualize your schema
==============================================================================

Evaluate the number of alelles and allele size variation for the loci in a schema or for a set
of selected loci. Provide information about problematic alleles per locus and individual pages
for each locus with a plot with allele size, a Neighbor Joining tree based on a multiple sequence
alignment (MSA) and a visualization of the MSA.

Basic Usage
:::::::::::

::

	chewBBACA.py SchemaEvaluator -g /path/to/SchemaName -o /path/to/OutputFolderName --cpu 4 --loci-reports

Parameters
::::::::::

::

    -g, --schema-directory      (Required) Path to the schema's directory (default: None).

    -o, --output-directory      (Required) Path to the output directory where the report HTML
                                files will be generated (default: None).

    --gl, --genes-list          (Optional) Path to a file with the list of genes in the schema
                                that the process should analyse (one per line) (default: False).

    -a, --annotations           (Optional) Path to the TSV file created by the UniprotFinder module (
                                default: None).

    --ta, --translation-table   (Optional) Genetic code used to translate coding sequences (default: None).

    --st, --size-threshold      (Optional) Allele size variation threshold. If an allele has a size
                                 within the interval of the locus mode -/+ the threshold (default: None).

    --ml, --minimum-length      (Optional) Minimum sequence length accepted for a coding sequence to
                                be included in the schema (default: None).

    --cpu, --cpu-cores          (Optional) Number of CPU cores/threads that will be used to run the
                                process (will be redefined to a lower value if it is equal to or
                                exceeds the totalnumber of available CPU cores/threads) (default: 1).

    --loci-reports              (Optional) If the process should create an individual report for each
                                locus (default: False).

    --light                     (Optional) If determining loci reports, skips MSA computation with
                                MAFFT and does not add the Phylogenetic Tree and MSA components
                                (default: False).

    --add-sequences             (Optional) Adds Code Editor with the DNA and Protein sequences to
                                loci reports. The Code Editor is in readonly mode (allows to search
                                for and copy text) (default: False).

Outputs
:::::::

- A main report HTML file that contains information about the schema.
- Various HTML files containing a detailed report about each locus.

.. note::
	If the ``--no-cleanup`` flag is used, then the following intermediate files will not be removed:

		- The pre-computed data in JSON format contain the data necessary to generate the plots and
		  tables.

		- The output files from MAFFT: multiple sequence alignments in FASTA format and the
		  Neighbour-Joining tree in Newick format.

		- Various JSON files containing exceptions encountered during the analysis such as: short
		  sequences, invalid start/stop codons and sequence lengths not being a multiple of 3.

Main report components
----------------------

The first component gives a small introduction that details the type of information contained in
the main and individual locus reports.

.. image:: https://user-images.githubusercontent.com/33930141/116111815-e87da200-a6ae-11eb-95d4-ee4a74a71e96.png
	:width: 1400px
	:align: center

Schema Summary Statistics
.........................

The second component displays summary statistics about the schema such as:

- chewBBACA version used to create it.
- BLAST Score Ratio (BSR) used to create it.
- Total no. of Loci.
- Total no. of Alleles.
- Total no. of Alleles not multiple of 3.
- Total no. of Alleles w/ >1 stop codons.
- Total no. of Alleles wo/ Start/Stop Codon.
- Total no. of Alleles shorter than ``--ml``, the minimum sequence length (in no. of nucleotides).

.. image:: https://user-images.githubusercontent.com/33930141/116112126-30042e00-a6af-11eb-9647-bba82ce433eb.png
	:width: 1400px
	:align: center

Loci with high variability
..........................

This analysis calculates the mode size per locus and using that value -/+ a threshold
(0.05 default) considers an allele "conserved" if it falls within the sequence length interval.
The user is given the choice of threshold and the choice to consider if a locus is classified
as having "high length variability" if 1 allele is outside the threshold (default) or to be
less stringent and classify a locus as having "high length variability" if >1 of the alleles
is outside the threshold.

.. image:: https://user-images.githubusercontent.com/33930141/116112200-414d3a80-a6af-11eb-83a5-bbaa37ca0c87.png
	:width: 1400px
	:align: center

Loci with only one allele
.........................

The module detects loci that have a single allele, allowing the users to quickly identify possible
problematic loci.

.. image:: https://user-images.githubusercontent.com/33930141/116112246-4ad6a280-a6af-11eb-92e8-9087d0d3d2ef.png
	:width: 1400px
	:align: center

In both tables, clicking on the locus name will open the individual report HTML for that locus.

Loci shorter than the minimum sequence length threshold
.......................................................

This table displays the loci that are shorter than the value passed to the ``--ml`` parameter.

.. image:: https://user-images.githubusercontent.com/33930141/116112665-abfe7600-a6af-11eb-81a6-2c930f7afbb2.png
	:width: 1400px
	:align: center

Schema Evaluation
.................

The third component contains 4 panels with summary charts displaying relevant information about
the schema. The panel is presented in the same way as in Chewie-NS.

- Panel A displays the distribution of loci by number of alleles.

.. image:: https://user-images.githubusercontent.com/33930141/102388113-37148480-3fc9-11eb-9dc4-963837eb8663.png
	:width: 1400px
	:align: center

- Panel B displays the distribution of loci by allele mode size.

.. image:: https://user-images.githubusercontent.com/33930141/105173595-294aa580-5b19-11eb-8b40-69223e760084.png
	:width: 1400px
	:align: center

- Panel C contains a representation of summary statistics (minimum allele size in blue, maximum
  allele size in orange and median size in green).

.. image:: https://user-images.githubusercontent.com/33930141/102388587-e0f41100-3fc9-11eb-840a-09ed0437839e.png
	:width: 1400px
	:align: center

- Panel D displays box plots of locus size distribution.

.. image:: https://user-images.githubusercontent.com/33930141/102388782-20baf880-3fca-11eb-9e88-1dba1b73dab1.png
	:width: 1400px
	:align: center

Loci Analysis
.............

The final component of the report presents a stacked bar chart and a table. In this component the
alleles of each locus are checked for their integrity as CDSs. The table includes the
*Uniprot Annotation*, the product name found through UniProt's SPARQL endpoint, and the
*Proteome Product*, the product name attributed based on high similarity to proteins included
in UniProt's reference proteomes. In addition, the *Missing Allele IDs* column presents the IDs
of alleles that are missing in the initial list of each locus and the *Total Invalid Alleles*
and *Valid Alleles* columns present the sum of invalid alleles and the total no. of valid alleles,
respectively.

The stacked bar chart presents, per locus, and sorted by the total number of alleles, the number
of alleles per locus. The alleles are divided into 5 classes:

	a) more than one stop codon (green);
	b) allele length not a multiple of 3 (orange);
	c) missing start or stop codon (red);
	d) alleles shorter than the ``--ml`` minimum length (purple);
	e) the number of alleles which are valid CDSs (blue).

.. note::
	In order to identify the *Missing Allele IDs*, the module expects the headers of the input
	FASTA files to have the locus identifier followed by the allele integer identifier
	(e.g.: >lmo_1) or simply the allele integer identifier (e.g.: >1).

.. image:: https://user-images.githubusercontent.com/33930141/116113169-27f8be00-a6b0-11eb-99a4-a03e8e8fedc7.png
	:width: 1400px
	:align: center

.. image:: https://user-images.githubusercontent.com/33930141/105173895-9b22ef00-5b19-11eb-9013-9db6835d2704.png
	:width: 1400px
	:align: center

Individual Report Components
----------------------------

Clicking on a point (locus) on Panel C or Panel D or on the name of the locus on the Loci
Analysis table will open a new page containing a detailed report about the selected locus.

Locus Individual Analysis
.........................

The first component presents a panel with 2 charts:

- A histogram summarizing the size distribution of the alleles (frequency of binned sizes).

- A scatter plot representing the actual sizes of each allele ordered by allele number.

.. note::
	The red line represents the minimum sequence value, ``--ml``, minus a size variation threshold
	of 20% (the default value for the size variation threshold used by the AlleleCall module).
	Alleles shorter than this value are below the size variation threshold. The yellow area
	represents the values that are within the size threshold.

.. image:: https://user-images.githubusercontent.com/33930141/116114802-9d18c300-a6b1-11eb-90d5-5b86a721b095.png
	:width: 1400px
	:align: center

.. image:: https://user-images.githubusercontent.com/33930141/116114827-a3a73a80-a6b1-11eb-8a69-d9f53ef8aa19.png
	:width: 1400px
	:align: center

Locus Information
.................

The second component presents a table containing the CDS analysis of the selected locus. It also
presents 4 new columns, in comparison with the table in the *Loci Analysis* of the main report,
with information on the:

- Number of alleles.
- Size Range, in nucleotides (nt).
- Allele median size (nt).
- Allele mode size (nt).

.. image:: https://user-images.githubusercontent.com/33930141/105175131-6b74e680-5b1b-11eb-845f-5121c91cf5be.png
	:width: 1400px
	:align: center

Exceptions
..........

The third component displays a table containing the list of alleles that are considered exceptions
based on the parameters used to evaluate the schema.

.. image:: https://user-images.githubusercontent.com/33930141/105175517-f524b400-5b1b-11eb-9554-e2094d4c1639.png
	:width: 1400px
	:align: center

NJ Tree
.......

The fourth component displays a Neighbor Joining tree built by ClustalW2 based on the
`MAFFT <https://mafft.cbrc.jp/alignment/software/>`_ alignment. The tree visualization
is produced using `Phylocanvas <http://phylocanvas.org/>`_.

.. image:: https://user-images.githubusercontent.com/33930141/105175900-6c5a4800-5b1c-11eb-98c3-f8e4beb15d6b.png
	:width: 1400px
	:align: center

Sequence Logo
.............

The fifth component displays a sequence logo obtained from the multiple sequence alignment
produced by `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_.

The *Change mode to frequency/information_content* button allows users to change how letter
heights are computed.

.. image:: https://user-images.githubusercontent.com/33930141/116115456-51b2e480-a6b2-11eb-88ad-747d542f9e98.png
	:width: 1400px
	:align: center

Multiple Sequence Analysis
..........................

The final component of the individual report presents the multiple sequence alignment produced by
`MAFFT <https://mafft.cbrc.jp/alignment/software/>`_. In order to visualize a different region of
the alignment, hover over the alignment until the hand cursor appears and then drag the alignment
to check the remaining rows and columns.

.. image:: https://user-images.githubusercontent.com/33930141/105175977-885de980-5b1c-11eb-86ad-b68b13f09cb0.png
	:width: 1400px
	:align: center
