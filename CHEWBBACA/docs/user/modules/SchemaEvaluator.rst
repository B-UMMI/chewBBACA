SchemaEvaluator - Build an interactive report for schema evaluation
===================================================================

The SchemaEvaluator module allows users to generate an interactive HTML report to better explore
the structure and evaluate typing schemas. The report includes data tables and charts with detailed
information about the allele size variation and integrity for each locus in the schema. The module
can be used to analyse schemas created with chewBBACA and external schemas from platforms such as
`Ridom cgMLST <http://www.cgmlst.org/ncs>`_, `BIGSdb <https://pubmlst.org/>`_,
`BIGSdb-Pasteur <https://bigsdb.pasteur.fr/>`_ and `Enterobase <http://enterobase.warwick.ac.uk/>`_.

Basic Usage
:::::::::::

::

	chewBBACA.py SchemaEvaluator -g /path/to/SchemaDirectory -o /path/to/OutputFolderName --cpu 4 --loci-reports

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

::

	OutputFolderName
	├── loci_reports
	│   ├── locus1.html
	│   ├── ...
	│   ├── locusN.html
	│   └── loci_bundle.js
	├── schema_report.html
	└── schema_bundle.js

- A HTML report, ``schema_report.html``, that contains the following components:

  - A table with summary data about the schema.
  - A tab panel with charts for the distribution of the number of alleles and allele size per locus.
  - If a TSV file with annotations is provided to the ``--annotations`` parameter, the schema report
    will also include a table with the provided annotations. Otherwise, it will display a warning informing that
    no annotations were provided.
  - A table with the results of the allele integrity and diversity analysis per locus.

- A HTML report per locus that contains the following components:

  - A table with summary data about the locus.
  - A tab panel with charts for the locus allele diversity.
  - A table with the total number and list of alleles that encode each distinct protein.
  - A table with the list of alleles that are not complete coding sequences and/or that are
    considered size outliers based on the minimum length and size threshold values used for
    schema evaluation.
  - A MSA component with the alignment for the distinct proteins.
  - A tree drawn with Phylocanvas based on the Neighbor-Joining (NJ) tree created by MAFFT.
  - If the ``--add-sequences`` parameter is provided, the locus report will also include a
    code editor with the allele DNA sequences and a code editor with the distinct protein
    sequences.

- Two JavaScript bundle files. The ``schema_bundle.js`` is used by the schema report, and the ``loci_bundle.js``,
  located inside the ``loci_reports`` folder is used by the loci reports.

.. note::
  You need to provide the ``--loci-reports`` parameter if you want a detailed report per locus.
  Running the SchemaEvaluator module with its default parameter values will only generate the Schema
  Report.

.. warning::
  The JS bundles are necessary to visualize the HTML reports. Do not delete these files. You should
  not move or delete any of the files in the output folder. If you want to share the report files,
  simply compress the output folder and share the compressed archive. The receiver can simply uncompress
  the archive and open the HTML files in a browser to visualize the report.

Schema Report Components
------------------------

The first component gives a small introduction that details the type of information contained in
the schema report.

.. image:: /_static/images/schema_report_description.png
   :width: 1400px
   :align: center

Schema Summary Data
...................

The second component is a table with summary statistics about the schema such as:

- Total no. of loci in the schema/evaluated.
- Total no. of alleles.
- Total no. of valid alleles.
- Total no. of invalid alleles.
- Total no. of incomplete alleles (sequence size not multiple of 3).
- Total number of alleles that contain ambiguous bases.
- Total no. of alleles missing the Start and/or Stop codons.
- Total no. of alleles with in-frame stop codons.
- Total no. of alleles shorter than ``--ml``, the minimum sequence length (in no. of nucleotides).
- Total no. of alleles below the locus sequence size threshold.
- Total no. of alleles above the locus sequence size threshold.

.. image:: /_static/images/schema_report_summary.png
   :width: 1400px
   :align: center

Loci Statistics
...............

The third component contains 4 panels with summary charts displaying relevant information about
the schema. The panel is presented in the same way as in Chewie-NS.

- Panel A displays the distribution of loci by number of alleles.

.. image:: /_static/images/schema_report_panelA.png
   :width: 1400px
   :align: center

- Panel B displays the distribution of loci by allele mode size.

.. image:: /_static/images/schema_report_panelB.png
   :width: 1400px
   :align: center

- Panel C contains a representation of summary statistics (minimum allele size in blue, maximum
  allele size in orange and median size in green).

.. image:: /_static/images/schema_report_panelC.png
   :width: 1400px
   :align: center

- Panel D displays box plots of locus size distribution.

.. image:: /_static/images/schema_report_panelD.png
   :width: 1400px
   :align: center

Loci annotations
................

If a TSV file with loci annotations is provided, the fourth component of the schema report is a table
with the list of annotations provided.

.. image:: /_static/images/schema_report_annotations.png
   :width: 1400px
   :align: center

Allele Analysis
...............

The final component of the report presents a table. In this component the alleles of each locus are
checked for their integrity as CDSs. In addition, the *Missing Allele IDs* column presents the number
o fIDs of alleles that are missing in the initial list of each locus.

.. note::
	In order to identify the *Missing Allele IDs*, the module expects the headers of the input
	FASTA files to have the locus identifier followed by the allele integer identifier
	(e.g.: >lmo_1) or simply the allele integer identifier (e.g.: >1).

.. image:: /_static/images/schema_report_allele_analysis.png
   :width: 1400px
   :align: center

.. note::
	If the ``--loci-reports`` parameter was provided, clicking on a point (locus) on Panel C or
	Panel D or on the name of the locus on the Allele Analysis table will open a new page containing
	a detailed report about the selected locus.

Locus Report Components
-----------------------

The first component gives a small introduction that details the type of information contained in
the locus report.

.. image:: /_static/images/loci_reports_description.png
   :width: 1400px
   :align: center

Locus Summary Data
..................

The second component is a table with summary statistics about the locus such as:

- Locus identifier.
- Total no. of alleles.
- Total no. of valid alleles.
- Total no. of invalid alleles.
- Proportion of validated alleles.
- Distinct protein alleles.
- Total no. of incomplete alleles (sequence size not multiple of 3).
- Total number of alleles that contain ambiguous bases.
- Total no. of alleles missing the Start and/or Stop codons.
- Total no. of alleles with in-frame stop codons.
- Total no. of alleles shorter than ``--ml``, the minimum sequence length (in no. of nucleotides).
- Allele length range.
- Allele length median.
- Allele length mode.
- Total no. of alleles below the locus sequence size threshold.
- Total no. of alleles above the locus sequence size threshold.
- Number of missing allele IDs.

.. image:: /_static/images/loci_reports_summary.png
   :width: 1400px
   :align: center

Locus Annotation Data
.....................

.. image:: /_static/images/loci_reports_annotations.png
   :width: 1400px
   :align: center

Locus Size Plots
................

The next component presents a panel with 3 charts:

- A histogram summarizing the size distribution of the alleles (frequency of binned sizes).

.. image:: /_static/images/loci_reports_allele_size_counts.png
   :width: 1400px
   :align: center

- A scatter plot representing the size of each allele ordered by allele number.

.. image:: /_static/images/loci_reports_allele_size.png
   :width: 1400px
   :align: center

- A bar chart with the number of distinct alleles that encode each distinct protein.

.. image:: /_static/images/loci_reports_protein_alleles.png
   :width: 1400px
   :align: center

.. note::
	The red line represents the minimum sequence value, ``--ml``, minus a size variation threshold
	of 20% (the default value for the size variation threshold used by the AlleleCall module).
	Alleles shorter than this value are below the size variation threshold. The yellow area
	represents the values that are within the size threshold.

Distinct Protein Alleles
........................

The fith component presents a table with the list of distinct protein alleles and the list of
distinct alleles that encode for each protein alleles.

.. image:: /_static/images/loci_reports_protein_table.png
   :width: 1400px
   :align: center

Invalid Alleles and Size Outliers
.................................

The sixth component presents a table with the list of alleles that are invalid and/or that are
size outliers based on the minimum length and size threshold values.

.. image:: /_static/images/loci_reports_invalid_alleles.png
   :width: 1400px
   :align: center

Multiple Sequence Alignment
...........................

The seventh component of the locus report presents the multiple sequence alignment produced by
`MAFFT <https://mafft.cbrc.jp/alignment/software/>`_.

.. image:: /_static/images/loci_reports_msa.png
   :width: 1400px
   :align: center

Neighbor-Joining Tree
.....................

The next component displays a Neighbor-Joining tree based on the
`MAFFT <https://mafft.cbrc.jp/alignment/software/>`_ alignment. The tree visualization
is produced using `Phylocanvas.gl <https://www.npmjs.com/package/@phylocanvas/phylocanvas.gl>`_.

.. image:: /_static/images/loci_reports_nj.png
   :width: 1400px
   :align: center

DNA sequences and Protein sequences
...................................

If the ``--add-sequences`` parameter was provided, the report will include two Monaco Code Editor components.

.. image:: /_static/images/loci_reports_dna_editor.png
   :width: 1400px
   :align: center

.. image:: /_static/images/loci_reports_protein_editor.png
   :width: 1400px
   :align: center
