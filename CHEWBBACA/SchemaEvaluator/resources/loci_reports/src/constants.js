// Constants used in other modules

export const headerMessage = `
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

export let alertMessages = [
	["List of invalid alleles is not displayed because all alleles are considered valid.", "warning"],
	["The locus annotations are not displayed because the user did not provide a file with annotations or there were no annotations for this locus.", "warning"],
	["The red lines represent the bot and top allele length thresholds. The bar corresponding to the length mode is colored in green.", "info"],
	["The red lines represent the bot and top allele length thresholds.", "info"],
	["Valid alleles are colored in blue and invalid alleles are colored in grey.", "info"],
	["Nodes are labeled with the allele identifiers. A tree might not be displayed if the distance between all alleles is 0.", "info"],
	["Use the mouse wheel to Zoom In/Out. In the Hierarchical and Rectangular Tree Types, holding the Ctrl key while using the mouse wheel allows to resize the branches.", "info"],
	["The NJ tree and MSA components are not displayed because this locus does not have valid alleles or has only one valid allele or the --light flag was provided.", "warning"],
	["Displaying the MSA for the %VALUE% alleles that were considered valid.", "info"],
	["Click on a row/column label to select a row/column. Hold the Ctrl key for multiselection. The 'Row Selection' option exports the Full MSA if no rows are selected.", "info"],
	["Displaying the %VALUE% alleles in the schema FASTA file.", "info"],
	["Only displaying the %VALUE% alleles that were considered valid.", "info"],
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