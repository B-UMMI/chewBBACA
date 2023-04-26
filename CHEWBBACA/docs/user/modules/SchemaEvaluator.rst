SchemaEvaluator - Build an interactive report for schema evaluation
===================================================================

The SchemaEvaluator module allows users to generate an interactive HTML report to better explore
the structure and evaluate typing schemas. The report includes data tables and charts with detailed
information about the allele size variation and allele integrity for each locus in the schema. The module
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
  - A table with the locus annotations, if the user provided annotations for the locus. Otherwise, it will
    display a warning informing that no annotations were provided.
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
  located inside the ``loci_reports`` folder, is used by the loci reports.

.. note::
  You need to provide the ``--loci-reports`` parameter if you want a detailed report per locus.
  Running the SchemaEvaluator module with its default parameter values will only generate the Schema
  Report.

.. warning::
  The JS bundles are necessary to visualize the HTML reports. Do not delete these files. You should
  not move or delete any of the files in the output folder. If you want to share the report, simply
  compress the output folder and share the compressed archive. The receiver can simply uncompress
  the archive and open the HTML files in a browser to visualize the report.

Schema Report Components
------------------------

The first component gives a small introduction that details the type of information contained in
each component of the schema report.

.. image:: /_static/images/schema_report_description.png
   :width: 1400px
   :align: center

The two alerts on top of the expandable component provide information about the parameter values
used to create and evaluate the schema, respectively. The SchemaEvaluator can only determine the
parameter values used for schema creation if the schema was created with chewBBACA.

Schema Summary Data
...................

The second component is a table with summary statistics about the schema, such as:

  - **Loci**: Total number of loci that were evaluated.
  - **Alleles**: Total number of alleles.
  - **Valid Alleles**: Total number of valid alleles. An allele is considered valid if its sequence size is a multiple
    of 3, if it has a single start and stop codon and if it contains no ambiguous bases.
  - **Invalid Alleles**: Total number of invalid alleles. The value in this column is the sum of the values in the ``Incomplete ORF``,
    ``Ambiguous Bases``, ``Missing Start/Stop Codon`` and ``In-frame Stop Codon`` columns.
  - **Incomplete ORF**: Total number of incomplete alleles (sequence size not multiple of 3).
  - **Ambiguous Bases**: Total number of alleles that contain ambiguous bases (non-ACTG characters).
  - **Missing Start/Stop Codon**: Total number of alleles missing the Start and/or Stop codons.
  - **In-frame Stop Codon**: Total number of alleles with in-frame stop codons.
  - **Alleles <bp**: Total number of alleles shorter than ``--ml``, the minimum sequence length value used
    for schema evaluation (in number of nucleotides).
  - **Alleles below threshold**: Total number of alleles below the locus sequence size bot threshold. This threshold identifies
    alleles with a sequence size that is -20% of the allele size mode.
  - **Alleles above threshold**: Total number of alleles above the locus sequence size top threshold. This threshold identifies
    alleles with a sequence size that is +20% of the allele size mode.

.. image:: /_static/images/schema_report_summary.png
   :width: 1400px
   :align: center

Loci Statistics
...............

The third component contains 4 panels with charts displaying relevant information about
the distribution of the number of alleles and allele size variation per evaluated locus.

- Panel A, ``Total Alleles``, displays the distribution of loci by number of alleles.

.. image:: /_static/images/schema_report_panelA.png
   :width: 1400px
   :align: center

- Panel B, ``Allele Mode Size``, displays the distribution of loci by allele mode size.

.. image:: /_static/images/schema_report_panelB.png
   :width: 1400px
   :align: center

- Panel C, ``Locus Statistics``, displays a scatter chart with points for the minimum allele size (blue), maximum allele
  size (orange) and median allele size (green) per locus.

.. image:: /_static/images/schema_report_panelC.png
   :width: 1400px
   :align: center

- Panel D, ``Allele Size Variation``, displays box plots for the locus size distribution. The range slider
  beneath the xaxis line can be used to redefine the boxplots that are visible in the plot area.

.. image:: /_static/images/schema_report_panelD.png
   :width: 1400px
   :align: center

.. note::
  If you have provided the ``--loci-reports`` parameter, the points in Panel C and the
  boxplots in Panel D are clickable and will open the detailed report of the selected locus.

Loci annotations
................

If a TSV file with loci annotations is provided, the fourth component of the schema report is a table
with the list of annotations. Otherwise, it will display a warning informing that no annotations
were provided.

.. image:: /_static/images/schema_report_annotations.png
   :width: 1400px
   :align: center

If you have provided the ``--loci-reports`` parameter, the loci identifiers in the first column will
link to the loci report pages. If a column name includes ``URL``, the SchemaEvaluator module assumes
that the values in that column are URLs and creates links to the web pages.

.. important::
  The first column in the TSV file with annotations must be named ``Locus`` and contain the identifiers
  of the loci (the basename of the locus FASTA file without the ``.fasta`` extension).

You can use the :doc:`UniprotFinder </user/modules/UniprotFinder>` module to annotate the loci in a schema
created with chewBBACA. If you want to annotate an external schema, you can adapt it with the
:doc:`PrepExternalSchema </user/modules/PrepExternalSchema>` module followed by annotation with the
:doc:`UniprotFinder </user/modules/UniprotFinder>` module.

Allele Analysis
...............

The final component of the schema report presents a table with the results of the allele integrity and
diversity analysis per locus. The table includes values per locus for most column categories in the
``Schema Summary Data`` table. It also includes the following additional columns:

  - **Proportion of Validated Alleles**: the proportion of the total alleles in the locus FASTA file that
    were considered valid.
  - **Distinct Protein Alleles**: the number of distinct protein alleles encoded by all alleles.
  - **Missing Allele IDs**: the number of allele identifiers that are missing, assuming that allele identifiers
    in the FASTA file should be sequential.

.. note::
	In order to identify the *Missing Allele IDs*, the module expects the headers of the input
	FASTA files to have the locus identifier followed by the allele integer identifier
	(e.g.: >lmo_1) or simply the allele integer identifier (e.g.: >1).

.. image:: /_static/images/schema_report_allele_analysis.png
   :width: 1400px
   :align: center

Locus Report Components
-----------------------

The first component gives a small introduction that details the type of information contained in
the locus report.

.. image:: /_static/images/loci_reports_description.png
   :width: 1400px
   :align: center

Locus Summary Data
..................

The second component is a table that includes the values for the locus presented in the ``Allele Analysis``
table and also includes the additional values:

- **Size Range (bp)**: the allele size range (minimum-maximum).
- **Length Median (bp)**: the allele median size.
- **Length Mode (bp)**: the allele mode size.

.. image:: /_static/images/loci_reports_summary.png
   :width: 1400px
   :align: center

Locus Annotation Data
.....................

The third component is a table with the annotations provided for the locus. An alert will be displayed if there
are no annotations for the locus.

.. image:: /_static/images/loci_reports_annotations.png
   :width: 1400px
   :align: center

Locus Size Plots
................

The fourth component contains 3 panels with charts displaying relevant information about the distribution
of allele sizes, the sequence size per allele and the diversity of distinct proteins.

- Panel A, ``Allele Size Counts``, display a histogram summarizing the size distribution of the alleles (frequency
  of binned sizes).

.. image:: /_static/images/loci_reports_allele_size_counts.png
   :width: 1400px
   :align: center

.. note::
	The bar corresponding to the allele size mode is colored in green.

- Panel B, ``Allele Size``, displays a scatter chart representing the size of each allele ordered by allele identifier.

.. image:: /_static/images/loci_reports_allele_size.png
   :width: 1400px
   :align: center

.. note::
	The points corresponding to valid alleles are colored in blue and points for invalid alleles are colored in grey.

- Panel C, ``Alleles Per Protein``, displays a bar chart with the number of distinct alleles that encode each
  distinct protein.

.. image:: /_static/images/loci_reports_protein_alleles.png
   :width: 1400px
   :align: center

.. note::
   In Panels A and B, the ``Show Thresholds`` switch can be toggled to adjust the axes limits to show the
   bot and top allele size thresholds (with the default parameter values, the thresholds are defined based
   on a -/+20% size variation from the allele size mode).

Distinct Protein Alleles
........................

The fith component presents a table with the list of distinct protein alleles and the list of
distinct alleles that encode for each protein allele. The identifiers of the protein alleles
are selected based on the first distinct allele that encodes for the protein.

.. image:: /_static/images/loci_reports_protein_table.png
   :width: 1400px
   :align: center

Invalid Alleles and Size Outliers
.................................

The sixth component presents a table with the list of alleles that are invalid and/or that are considered size
outliers based on the minimum length and size threshold values used for schema evaluation. The ``Exception Category``
is defined based on the first exception captured for each allele. The list of all exceptions captured for each allele
is displayed in the ``Exception Description`` column.

.. image:: /_static/images/loci_reports_invalid_alleles.png
   :width: 1400px
   :align: center

Multiple Sequence Alignment
...........................

The seventh component of the locus report presents the protein multiple sequence alignment (MSA) produced by
`MAFFT <https://mafft.cbrc.jp/alignment/software/>`_. The MSA only includes the distinct proteins encoded by
the valid alleles.

.. image:: /_static/images/loci_reports_msa.png
   :width: 1400px
   :align: center

Neighbor-Joining Tree
.....................

The eigth component displays the guide tree created by `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_
alignment. The tree visualization is produced using `Phylocanvas.gl <https://www.npmjs.com/package/@phylocanvas/phylocanvas.gl>`_.


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
