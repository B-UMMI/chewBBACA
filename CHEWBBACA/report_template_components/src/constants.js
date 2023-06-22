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
