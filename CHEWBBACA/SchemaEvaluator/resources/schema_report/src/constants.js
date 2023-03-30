// Constants used in other modules

export const headerMessage = `
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

export let alertMessages = [
	["Each point (locus) is clickable and will open a page with more details about it.", "info"],
	["Each boxplot (locus) is clickable and will open a page with more details about it.", "info"],
	["The loci annotations are not displayed because the user did not provide a file with annotations.", "warning"],
];

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
