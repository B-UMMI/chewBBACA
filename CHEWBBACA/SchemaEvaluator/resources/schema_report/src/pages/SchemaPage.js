import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import TabPanelMUI from '../components/TabPanelMUI';
import AccordionMUI from '../components/AccordionMUI';

import Stack from '@mui/material/Stack';
import Alert from '@mui/material/Alert';
import Typography from '@mui/material/Typography';

import React from 'react'
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';


const SchemaPage = () => {

	const markdown = `
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

	// get pre-computed data
	const data = window.preComputedData;

	// data for Summary Data table
	const summaryData = data.summaryData;
	const summaryTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: "schema_summary.tsv",
			separator: "\t",
			filterOptions: {
				useDisplayedColumnsOnly: true,
				useDisplayedRowsOnly: true
			}
		},
		filter: false,
		search: false,
		viewColumns: true,
		pagination: false,
		draggableColumns: {
			enabled: true
		}
	};
	// Component for Summary table
	const summaryTable = <DataTable
						  tableData={summaryData} 
						  tableTitle="Schema Summary Data" 
						  tableOptions={summaryTableOptions}
						 >
						 </DataTable>

	// data for Panel A (Total Alleles Distribution)
	const xDataPanelA = data.total_alleles;
	const yDataPanelA = data.loci;
	const plotDataPanelA = [
		{x: xDataPanelA,
		 y: yDataPanelA,
		 type: "histogram",
		 name: "Total alleles",
		 marker: {
			 color: "#0570b0",
			 line: {
				 color: "#a6bddb",
				 width: 1
			 }
		 }
	    }
	];
	const layoutPanelA = {
		title: {
			text: "Distribution of the number of alleles"
		},
		xaxis: {
			title: { text: "Number of Alleles" },
			showgrid: false,
			zeroline: false,
			showline: true,
			ticks: "outside"
		},
		yaxis: {
			title: { text: "Number of Loci" },
			showgrid: true,
			zeroline: false,
			showline: true,
			ticks: "outside"
			//dtick: 100
		},
		bargroupgap: 0.05
	};
	const configPanelA = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'TotalAlleles',
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with total alleles distribution
	const TotalAllelesHistogram = (
		<PlotlyPlot
		 key="TotalAllelesHistogram"
		 plotData={plotDataPanelA}
		 layoutProps={layoutPanelA}
		 configOptions={configPanelA}
		>
		</PlotlyPlot>
	);

	// data for Panel B
	const xDataPanelB = data.mode;
	const yDataPanelB = data.loci;
	const plotDataPanelB = [
		{x: xDataPanelB,
		 y: yDataPanelB,
		 type: "histogram",
		 name: "Allele Mode Size",
		 marker: {
			color: "#0570b0",
			line: {
				color: "#a6bddb",
				width: 1
			}
		}
	    }
	];
	const layoutPanelB = {
		title: {
			text: "Distribution of the allele mode size"
		},
		xaxis: {
			title: { text: "Allele Mode Size" },
			showgrid: false,
			zeroline: false,
			showline: true,
			ticks: "outside"
		},
		yaxis: {
			title: { text: "Number of Loci" },
			zeroline: false,
			showline: true,
			ticks: "outside"
		},
		bargroupgap: 0.05
	};
	const configPanelB = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'AlleleModeSize',
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with allele mode size distribution
	const AlleleModeHistogram = (
		<PlotlyPlot
		 key="AlleleModeHistogram"
		 plotData={plotDataPanelB}
		 layoutProps={layoutPanelB}
		 configOptions={configPanelB}
		>
		</PlotlyPlot>
	);

	// data for Panel C
	const xDataMinPanelC = data.min
	const xDataMaxPanelC = data.max
	const xDataMedianPanelC = data.median

	const plotDataPanelC = [
		{
		  x: xDataMinPanelC,
		  y: xDataPanelA,
		  type: "scatter",
		  name: "Minimum",
		  mode: "markers",
		  marker: {
			opacity: 0.7,
			size: 4,
		  },
		  hovertemplate: "<b>ID</b>: %{text}",
		  text: yDataPanelB,
		},
		{
		  x: xDataMaxPanelC,
		  y: xDataPanelA,
		  type: "scatter",
		  name: "Maximum",
		  mode: "markers",
		  marker: {
			opacity: 0.7,
			size: 4,
		  },
		  hovertemplate: "<b>ID</b>: %{text}",
		  text: yDataPanelB,
		},
		{
		  x: xDataMedianPanelC,
		  y: xDataPanelA,
		  type: "scatter",
		  name: "Median",
		  mode: "markers",
		  marker: {
			opacity: 0.7,
			size: 4,
		  },
		  hovertemplate: "<b>ID</b>: %{text}",
		  text: yDataPanelB,
		}
	];

	// Event handler to open locus page on point click
	const clickPlotPointHandler = (pointID) => {
		// Create anchor element and associate link to locus HTML report
		const anchor=document.createElement("a");
		anchor.href=`./loci_reports/${pointID}.html`
		anchor.target="_blank";
		anchor.rel="noopener noreferrer";
		anchor.click();
	};

	const layoutPanelC = {
		title: {
			text: "Locus Length Statistics"
		},
		xaxis: {
			title: { text: "Allele size (bp)" },
			showgrid: true,
			zeroline: false,
			showline: true,
			ticks: "outside"
		},
		yaxis: {
			title: { text: "Number of alleles" },
			zeroline: false,
			showline: true,
			ticks: "outside"
		}
	};
	const configPanelC = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'LocusStatistics',
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};

	let scatterClickHandler = undefined;
	if (data.lociReports === 1) {
		// Define handler to get data from click in scatter point
		scatterClickHandler = (event) => {
			return clickPlotPointHandler(event.points[0].text)
		}
	};

	// Component for Plotly Scatter with loci statistics
	const LociStatsScatter = (
		<PlotlyPlot
		 key="LociStatsScatter"
		 plotData={plotDataPanelC}
		 layoutProps={layoutPanelC}
		 configOptions={configPanelC}
		 onClick={scatterClickHandler} // Get locus ID associated to the point that was clicked
		>
		</PlotlyPlot>
	);

	// data for Panel D
	const q1PanelD = data.q1
	const q3PanelD = data.q3

	const plotDataPanelD = [
		{
		type: "box",
		name: "Locus Size Variation",
		x: data.loci,
		q1: q1PanelD,
		median: xDataMedianPanelC,
		q3: q3PanelD,
		lowerfence: xDataMinPanelC,
		upperfence: xDataMaxPanelC,
		showlegend: false,
		marker: {color: "#0570b0"}
	  }
	];

	// Define xaxis maximum value to show
	let xaxisMax = 100;
	const totalLoci = summaryData[1].rows[0][0];
	// Show all boxplots if total number of loci < 100
	if (totalLoci <= xaxisMax) {
		xaxisMax = totalLoci;
	}
	const layoutPanelD = {
		title: {
			text: "Locus Size Variation"
		},
		xaxis: {
			title: { text: "Loci" },
			showgrid: false,
			zeroline: false,
			showline: true,
			rangeslider: {
				thickness: 0.15, 
			},
			showticklabels: false,
			range: [0, xaxisMax]
		},
		yaxis: {
			title: { text: "Allele Size" },
			zeroline: false,
			showline: true,
			ticks: "outside"
		},
		hovermode: 'x unified'
	};

	const configPanelD = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'AlleleSizeVariation',
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};

	let boxplotClickHandler = undefined;
	if (data.lociReports === 1) {
		// Define handler to get data from click in boxplot
		boxplotClickHandler = (event) => {
			return clickPlotPointHandler(event.points[0].x)
		}
	};

	// Component for Plotly Boxplots with loci allele size variation
	const LociSizeBoxplots = (
		<PlotlyPlot
		 key="LociSizeBoxplots"
		 plotData={plotDataPanelD}
		 layoutProps={layoutPanelD}
		 configOptions={configPanelD}
		 onClick={boxplotClickHandler} // Get locus ID associated to the point that was clicked
		>
		</PlotlyPlot>
	);

	// Create Alert to tell users that the points/boxplots are clickable
	// and link to the loci reports
	let scatterClickAlert = undefined;
	let boxplotClickAlert = undefined;
	if (data.lociReports === 1) {
		scatterClickAlert = (
			<Alert key="clickScatterAlert" severity="info">
				<Typography sx={{ fontSize: 14 }}>
					Each point (locus) is clickable and will
					open a page with more details about it.
				</Typography>
			</Alert>
		)
		boxplotClickAlert = (
			<Alert key="clickBoxplotAlert" severity="info">
				<Typography sx={{ fontSize: 14 }}>
					Each boxplot (locus) is clickable and will
					open a page with more details about it.
				</Typography>
			</Alert>
		)
	};

	const TabPanelTitles = [
		"Total Alleles", 
		"Allele Mode Size", 
		"Locus Statistics", 
		"Allele Size Variation"
	];
	const TabPanelData = [
		[TotalAllelesHistogram],
		[AlleleModeHistogram],
		[scatterClickAlert, LociStatsScatter],
		[boxplotClickAlert, LociSizeBoxplots]
	];

	const locusColumnFormatting = { 
		"Locus": {
			customBodyRender: (value) => (
				<a
		  			href={`./loci_reports/${value}.html`}
		  			target={"_blank"}
		  			rel={"noopener noreferrer"}
				>
		  			{value}
				</a>
	  		)
		}
	};

	const uniprotLinkColumnFormatting = {
		"UniProt_URL": {
			customBodyRender: (value) => (
				<a 
					href={value} 
					target={"_blank"} 
					rel={"noopener noreferrer"}
				>
					{value}
				</a>
			)
		}
	};

	// data for Annotations table
	const annotationsData = data.annotations;
	const annotationsTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: "loci_annotations.tsv",
			separator: "\t",
			filterOptions: {
				useDisplayedColumnsOnly: true,
				useDisplayedRowsOnly: true
			}
		},
		filter: true,
		filterType: 'multiselect',
		search: true,
		viewColumns: true,
		pagination: true,
		rowsPerPage: 10,
		rowsPerPageOptions: [10, 20, 40, 80, 100],
		jumpToPage: true,
		draggableColumns: {
			enabled: true
		}
	};

	// Select conditional formatting function to apply
	// Only add href to loci names if the loci reports were generated
	let LociAnnotationsFormatting = undefined;
	if (data.lociReports === 1) {
		LociAnnotationsFormatting = {...locusColumnFormatting, ...uniprotLinkColumnFormatting};
	} else {
		LociAnnotationsFormatting = {...uniprotLinkColumnFormatting};
	}
	// DataTable component to display loci annotations
	// Only define component if there is annotation data
	let LociAnnotationsTable = undefined;
	if (annotationsData[0].columns.length > 0) {
		LociAnnotationsTable = (
			<DataTable 
			 tableData={annotationsData} 
			 tableTitle="Loci Annotations" 
			 tableOptions={annotationsTableOptions}
			 tableConditionalFormatting={LociAnnotationsFormatting}
			>
			</DataTable>
		);
	}

	let LociAnnotationsAlert = undefined;
	if (LociAnnotationsTable === undefined) {
		LociAnnotationsAlert = (
			<Alert severity="warning">
				<Typography sx={{ fontSize: 14 }}>
					The loci annotations are not displayed because the
					user did not provide a file with annotations.
				</Typography>
			</Alert>
		)
	};

	// data for allele analysis
	const analysisData = data.analysis;
	const analysisTableOptions = {
		responsive: "simple",
		selectableRows: "none",
		selectableRowsHeader: false,
		selectableRowsHideCheckboxes: true,
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: "allele_analysis.tsv",
			separator: "\t",
			filterOptions: {
				useDisplayedColumnsOnly: true,
				useDisplayedRowsOnly: true
			}
		},
		filter: true,
		filterType: 'multiselect',
		search: true,
		viewColumns: true,
		pagination: true,
		rowsPerPage: 10,
		rowsPerPageOptions: [10, 20, 40, 80, 100],
		jumpToPage: true,
		draggableColumns: {
			enabled: true
		},
	};

	// Select conditional formatting function to apply
	// Only add href to loci names if the loci reports were generated
	let LociAnalysisFormatting = undefined;
	if (data.lociReports === 1) {
		LociAnalysisFormatting = {...locusColumnFormatting};
	}
	// DataTable component to display loci integrity analysis
	const LociAnalysisTable = (
		<DataTable 
		 tableData={analysisData} 
		 tableTitle="Allele Analysis" 
		 tableOptions={analysisTableOptions}
		 tableConditionalFormatting={LociAnalysisFormatting}
		>
		</DataTable>
	);

	const HeaderSummary = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Schema Report Description
		</Typography>
	);

	// Need to change the root component attribute to avoid error related with
	// passing something to <p>
	const HeaderDescription = (
		<Typography component={'div'} style={{ lineHeight: "18px" }} > 
			<ReactMarkdown 
				children={markdown} 
				remarkPlugins={[remarkGfm]}
			>
			</ReactMarkdown>
		</Typography>
	);

	// Get message about config used for schema creation
	const creationConfigMessage = data.creationConfig;
	const creationAlert = (
		<Alert severity="info">
			<Typography sx={{ fontSize: 14 }}>
				{creationConfigMessage}
			</Typography>
		</Alert>
	);

	// Get message about config used for schema evaluation
	const evaluationConfigMessage = data.evaluationConfig;
	const evaluationAlert = (
		<Alert severity="info">
			<Typography sx={{ fontSize: 14 }}>
				{evaluationConfigMessage}
			</Typography>
		</Alert>
	);

	return (
		<div style={{ marginTop: "10px", marginBottom: "10px" }}>
			<div style={{ marginBottom: "5px" }}>
				<Stack sx={{ width: '100%' }} spacing={0.5}>
					{creationAlert}
					{evaluationAlert}
				</Stack>
			</div>
			<div>
				<AccordionMUI 
					summaryText={HeaderSummary}
					detailsData={HeaderDescription} 
					expanded={false}
				>
				</AccordionMUI>
			</div>
			<div style={{ marginTop: "20px"}}>
				{summaryTable}
			</div>
			<div style={{ marginTop: "20px"}}>
				<TabPanelMUI
					ContentTitles={TabPanelTitles}
					ContentData={TabPanelData}
				>
				</TabPanelMUI>
			</div>
			<div style={{ marginTop: "20px"}}>
				<div style={{ marginBottom: "10px" }}>
					{LociAnnotationsAlert}
				</div>
				{LociAnnotationsTable}
			</div>
			<div style={{ marginTop: "20px"}}>
				{LociAnalysisTable}
			</div>
		</div>
	  );
};

export default SchemaPage;
