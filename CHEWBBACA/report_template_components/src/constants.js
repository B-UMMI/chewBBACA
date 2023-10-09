// Constants used in the reports

export const schemaReportHeaderMessage = `
### Schema Report

The schema report includes the following components:

- **Schema Summary Data**

  - Total number of loci and totals for the **Allele Analysis** categories.

- **Loci statistics**

  - Tab Panel with the following panels:
    - Total Alleles: distribution of the number of alleles for all loci.
    - Allele Mode Size: distribution of the allele mode size (most frequent allele size) for all loci.
    - Locus Statistics: total number of alleles and allele minimum, maximum and median size per locus.
    - Allele Size Variation: boxplots for the allele size variation per locus.  

	If you've provided the *--loci-reports* parameter to create the individual loci reports, the points in
	the "Locus Statistics" and the boxplots in the "Allele Size Variation" plots are clickable and will open
	the report for the corresponding locus.

- **Loci Annotations**
	- Annotations provided for each locus (added from the TSV file passed to the *--a* parameter).

	If the TSV file provided to the *--a* parameter included a column named "UniProt_URL" with the URLs to the UniProt records that matched each locus, the column values will link to the pages of those records.  
	If you've provided the *--loci-reports* parameter, the loci identifiers in the first column will link to the loci individual reports.
  
- **Allele Analysis**
	- Results of the allele analysis per locus. The results are organized into the following categories:
		- Total alleles: number of alleles in the locus FASTA file.
		- Valid alleles: number of alleles considered valid based on the configuration values.
		- Incomplete ORF: number of alleles with a size value not multiple of 3.
		- Missing Start/Stop Codon: number of alleles missing the start and/or stop codon.
		- In-frame Stop Codon: number of alleles with in-frame stop codons.
		- Alleles < minLen: number of alleles with smaller than the minimum sequence length defined in the configuration values.
		- Alleles below threshold: number of alleles smaller than the allele length mode bot threshold.
		- Alleles above threshold: number of alleles bigger than the allele length mode top threshold.

	If you've provided the *--loci-reports* parameter, the loci identifiers in the first column will link to the loci individual reports.

This report was created with the [React](https://react.dev/) library and uses the following JavaScript packages:

  - [Material UI](https://www.npmjs.com/package/@mui/material) for most of the components.
  - [MUI-Datatables](https://www.npmjs.com/package/mui-datatables) for datatables components.
  - [MSA Viewer](https://www.npmjs.com/package/@jlab-contrib/msa) for the Multiple Sequence Alignmnet.
  - [Phylocanvas.gl](https://www.npmjs.com/package/@phylocanvas/phylocanvas.gl) for the NJ Tree.
  - [react-plotly.js](https://www.npmjs.com/package/react-plotly.js) for the charts.
  - [react-scroll](https://www.npmjs.com/package/react-scroll) for animating vertical scrolling in the NJ Tree component.
  - [Monaco Editor for React](https://www.npmjs.com/package/@monaco-editor/react) for the read-only code editor that displays the DNA and Protein sequences.
  - [react-markdown](https://www.npmjs.com/package/react-markdown) and [remark-gfm](https://www.npmjs.com/package/remark-gfm) to render markdown.

You can find more information in the SchemaEvaluator module [documentation](https://chewbbaca.readthedocs.io/en/latest/user/modules/SchemaEvaluator.html).  
If you have any feature requests or issues to report, please visit chewBBACA's [GitHub repository](https://github.com/B-UMMI/chewBBACA).
`;

export const locusReportHeaderMessage = `
### Locus Report

The locus report includes the following components:

- **Locus Summary Data**

  - Locus identifier and totals for the **Allele Analysis** categories.
- **Locus Annotation Data**

  - Annotations provided for the locus (added from the TSV file passed to the *--a* parameter).
- **Locus Size plots**

  - Tab Panel with the following panels:
    - Allele Size Counts: distribution of the allele size values for the locus.
    - Allele Size: sequence size per allele.
- **Neighbor-Joining Tree**

  - A tree drawn with Phylocanvas based on the Neighbor-Joining (NJ) tree created by MAFFT.
- **Multiple Sequence Alignment**

  - Visualization of the MSA computed by MAFFT.
- **DNA sequences**

  - Text Editor in read-only mode with the alleles in the FASTA file.
- **Protein sequences**

  - Text Editor in read-only mode with the translated alleles that were considered valid based on the configuration values.

This report was created with the [React](https://react.dev/) library and uses the following JavaScript packages:

- [Material UI](https://www.npmjs.com/package/@mui/material) for most of the components.
- [MUI-Datatables](https://www.npmjs.com/package/mui-datatables) for datatables components.
- [MSA Viewer](https://www.npmjs.com/package/@jlab-contrib/msa) for the Multiple Sequence Alignmnet.
- [Phylocanvas.gl](https://www.npmjs.com/package/@phylocanvas/phylocanvas.gl) for the NJ Tree.
- [react-plotly.js](https://www.npmjs.com/package/react-plotly.js) for the charts.
- [react-scroll](https://www.npmjs.com/package/react-scroll) for animating vertical scrolling in the NJ Tree component.
- [Monaco Editor for React](https://www.npmjs.com/package/@monaco-editor/react) for the read-only code editor that displays the DNA and Protein sequences.
- [react-markdown](https://www.npmjs.com/package/react-markdown) and [remark-gfm](https://www.npmjs.com/package/remark-gfm) to render markdown.

You can find more information in the SchemaEvaluator module [documentation](https://chewbbaca.readthedocs.io/en/latest/user/modules/SchemaEvaluator.html).  
If you have any feature requests or issues to report, please visit chewBBACA's [GitHub repository](https://github.com/B-UMMI/chewBBACA).
`;

export const alleleCallReportHeaderMessage = `
### Allele Call Report

The allele call report includes the following components:

- **Results Summary Data**

  - The total number of samples, total number of loci, total number of coding sequences (CDSs) identified in the samples, total number of CDSs classified and totals per classification type.
- **Classification Counts**
  - Tab Panel with the following panels:
    - Counts per Sample: displays the stacked bar charts for the sample classification type count.
    - Counts per Locus: displays the stacked bar charts for the loci classification type counts.

	The samples/loci are sorted in order of decreasing number of valid (EXC+INF) classifications. The plot area
	will display at most data for 300 samples/loci. You can click the left/right arrows to view the previous/next
	300 samples/loci and the double left/right arrows to view the data for the first/last 300 samples/loci. The
	component includes a slider to select the range of sample/loci bars that are visible.
- **Detailed Statistics**
  - Tab Panel with the following panels:

    - Sample Stats: detailed sample statistics.
	- Loci Stats: detailed loci statistics.

	The dropdown menu below the tables allows the selection of a single column to generate a histogram for the values in the selected column.
- **Loci Annotations**
  - Annotations provided for each locus (added from the TSV file passed to the *--a* parameter).

    If the TSV file provided to the *--a* parameter included a column named "UniProt_URL" with the URLs to the UniProt records that matched each locus, the column values will link to the pages of those records.  
    If you've provided the *--loci-reports* parameter, the loci identifiers in the first column will link to the loci individual reports.
- **Loci Presence-Absence**

  - Heatmap representing the loci presence-absence matrix for all samples in the dataset. Blue cells (z=1) correspond to
    loci presence and grey cells (z=0) to loci absence. The "Select Sample" dropdown menu enables the selection of a single
	sample to display its heatmap on top of the main heatmap. The "Select Locus" dropdown menu enables the selection of a
	single locus to display its heatmap on the right of the main heatmap. 
- **Allelic distances**

  - Heatmap representing the allelic distance matrix for all samples in the dataset. The allelic distances are computed based
    on the cgMLST profiles (loci that are not present in all samples are not included in the computations). The "Select Sample"
	dropdown menu enables the selection of a single sample to display its heatmap on top of the main heatmap. The menu after
	the heatmap enables the selection of a single sample and of a distance threshold to display a table with the list of samples
	at a distance equal or smaller than the specified distance value.
- **Core-genome Neighbor-Joining Tree**

  - A tree drawn with Phylocanvas based on the Neighbor-Joining (NJ) tree computed by FastTree. The tree is computed based on the MSA for the set of loci that constitute the core-genome.

This report was created with the [React](https://react.dev/) library and uses the following JavaScript packages:

- [Material UI](https://www.npmjs.com/package/@mui/material) for most of the components.
- [MUI-Datatables](https://www.npmjs.com/package/mui-datatables) for datatables components.
- [MSA Viewer](https://www.npmjs.com/package/@jlab-contrib/msa) for the Multiple Sequence Alignmnet.
- [Phylocanvas.gl](https://www.npmjs.com/package/@phylocanvas/phylocanvas.gl) for the NJ Tree.
- [react-plotly.js](https://www.npmjs.com/package/react-plotly.js) for the charts.
- [react-scroll](https://www.npmjs.com/package/react-scroll) for animating vertical scrolling in the NJ Tree component.
- [Monaco Editor for React](https://www.npmjs.com/package/@monaco-editor/react) for the read-only code editor that displays the DNA and Protein sequences.
- [react-markdown](https://www.npmjs.com/package/react-markdown) and [remark-gfm](https://www.npmjs.com/package/remark-gfm) to render markdown.

You can find more information in the AlleleCallEvaluator module [documentation]().  
If you have any feature requests or issues to report, please visit chewBBACA's [GitHub repository](https://github.com/B-UMMI/chewBBACA).
`;

export let schemaReportAlertMessages = [
	["Each point (locus) is clickable and will open a page with more details about it.", "info"],
	["Each boxplot (locus) is clickable and will open a page with more details about it.", "info"],
	["The loci annotations are not displayed because the user did not provide a file with annotations.", "warning"],
];

export let locusReportAlertMessages = [
	["List of invalid alleles is not displayed because all alleles are considered valid.", "warning"],
	["The locus annotations are not displayed because the user did not provide a file with annotations or there were no annotations for this locus.", "warning"],
	["The red lines represent the bot and top allele length thresholds. The bar corresponding to the length mode is colored in green.", "info"],
	["The red lines represent the bot and top allele length thresholds.", "info"],
	["Valid alleles are colored in blue and invalid alleles are colored in grey.", "info"],
	["The protein allele identifiers were assigned based on the first allele in the schema FASTA file that encoded each distinct protein.", "info"],
	["Nodes are labeled with the allele identifiers. A tree might not be displayed if the distance between all alleles is 0.", "info"],
	["Use the mouse wheel to Zoom In/Out. In the Hierarchical and Rectangular Tree Types, holding the Ctrl key while using the mouse wheel allows to resize the branches.", "info"],
	["The NJ tree and MSA components are not displayed because this locus does not have valid alleles or has only one valid allele or the --light flag was provided.", "warning"],
	["Displaying the MSA for the %VALUE% distinct proteins.", "info"],
	["Click on a row/column label to select a row/column. Hold the Ctrl key for multiselection. The 'Row Selection' option exports the Full MSA if no rows are selected.", "info"],
	["Displaying the %VALUE% sequences in the schema FASTA file.", "info"],
	["Only displaying the %VALUE% alleles that were considered valid.", "info"],
	["This locus has %VALUE% alleles that have not been synchronized with the remote schema in Chewie-NS (%VALUE2%).", "warning"],
	["The Bar chart is not displayed because this locus does not have valid alleles.", "warning"],
	["The table for the distinct protein alleles is not displayed because this locus does not have valid alleles.", "warning"],
];

export let alleleCallReportAlertMessages = [
	["Click the left/right arrows to view the data for the previous/next 300 %VALUE%. Click the double left/right arrows to view the data for the first/last 300 %VALUE%.", "info"],
	["Click the left/right arrows to view the data for the previous/next 300 %VALUE%. Click the double left/right arrows to view the data for the first/last 300 %VALUE%.", "info"],
	["The loci annotations are not displayed because the user did not provide a file with annotations.", "warning"],
	["Select a sample/locus to display the heatmap for that sample/locus.", "info"],
	["Select a sample to display the heatmap for that sample.", "info"],
	["Select a sample and a distance value to display a table with the list of samples with a distance to the selected sample equal or smaller than the specified distance value.", "info"],
	["The distances were computed by determining the number of allelic differences from the set of %VALUE% core loci.", "info"],
	["The tree was computed based on the MSA for the set of %VALUE% core loci.", "info"],
];

export const msaColorSchemes = [
	["buried_index", "Buried"],
	["cinema", "Cinema"],
	["clustal", "Clustal"],
	["clustal2", "Clustal2"],
	["helix_propensity", "Helix propensity"],
	["hydro", "Hydro"],
	["lesk", "Lesk"],
	["mae", "Mae"],
	["nucleotide", "Nucleotide"],
	["purine_pyrimidine", "Purine Pyrimidine"],
	["strand_propensity", "Strand propensity"],
	["taylor", "Taylor"],
	["turn_propensity", "Turn propensity"],
	["zappo", "Zappo"],
];

export const msaDownloadOptions = ["Full MSA", "Row Selection", "Image (PNG)"];

export const treeDownloadOptions = ["SVG", "Newick", "JSON"];

export const globalTableOptions = {
	responsive: "simple",
	selectableRowsHeader: false,
	selectableRows: "none",
	selectableRowsOnClick: false,
	print: false,
	download: true,
	filter: false,
	search: false,
	viewColumns: true,
	pagination: true,
};

export const sampleTableHiddenColumns = ["EXC", "INF", "PLOT3", "PLOT5", "LOTSC", "NIPH", "NIPHEM", "ALM", "ASM", "PAMA", "LNF"];

export const lociTableHiddenColumns = ["EXC", "INF", "PLOT3", "PLOT5", "LOTSC", "NIPH", "NIPHEM", "ALM", "ASM", "PAMA", "LNF"];
