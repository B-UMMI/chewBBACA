import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import TabPanelMUI from '../components/TabPanelMUI';
import AccordionMUI from '../components/AccordionMUI';

import Typography from '@mui/material/Typography';
import Stack from '@mui/material/Stack';
import Alert from '@mui/material/Alert';

import React from 'react'
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

const SchemaPage = () => {

	const markdown = `
  ## Schema Evaluator
  
  The main report includes the following components:
  - Summary Data
	- Table with summary statistics for the schema.
  - Loci statistics
	- Tab Panel with four panels the users can alternate between.
		- Total Alleles: histogram for the number of alleles for all loci.
		- Allele Mode Size: histogram for the allele mode size (most frequent allele size) for all loci.
		- Locus Statistics: total number of alleles and allele minimum, maximum and median size per locus.
		- Allele Size Variation: boxplots for the allele size variation per locus.

		If you've provided the **--loci-reports** parameter to create the individual loci reports, the points in the Locus Statistics and the boxplots in the Allele Size Variation plots are clickable and will open the report for the corresponding locus.
  - Loci Annotations
	- Table that includes the annotations passed to the **--a** parameter.

	Adds links to Uniprot pages if the file includes a column named "UniPort_URL" with UniProt URLs.  
	If you've provided the **--loci-reports** parameter, the loci identifiers in the first column will link to the loci individual reports.
  - CDS Analysis
	- Table with the results of the CDS analysis for each locus.
		- Total alleles: number of sequences in the locus FASTA file.
		- Valid alleles: number of complete coding sequences.
		- Incomplete ORF:
		- Missing Start/Stop Codon:
		- In-frame Stop Codon:
		- Alleles < minLen:
		- Alleles below threshold:
		- Alleles above threshold:
	
		If you've provided the **--loci-reports** parameter, the loci identifiers in the first column will link to the loci individual reports.

  The loci individual reports include the following components:
  - Locus Summary Data
    - Table with summary statistics and CDS analysis for the locus.
  - Locus Annotation Data
    - Table that includes the annotations for the locus.
  - Locus Size plots
    - Tab Panel with two panels the users can alternate between.
		- Allele Size Counts: histogram for the allele size values for the locus.
		- Allele Size: scatter plot for the sequence size per allele.
  - Phylogenetic Tree
    - A tree drawn with Phylocanvas based on the Neighbor Joining tree created by MAFFT.
  - Multiple Sequence Alignment
    - Visualization of the MSA computed by MAFFT.
  - DNA sequences
    - Text Editor in read-only mode with the allele DNA sequences.
  - Protein sequences
    - Text Editor in read-only mode with the translated DNA sequences.
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
			filename: "schema:summary.tsv",
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
						  tableTitle="Summary Data" 
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
		 plotData={plotDataPanelD}
		 layoutProps={layoutPanelD}
		 configOptions={configPanelD}
		 onClick={boxplotClickHandler} // Get locus ID associated to the point that was clicked
		>
		</PlotlyPlot>
	);

	const TabPanelTitles = [
		"Total Alleles", 
		"Allele Mode Size", 
		"Locus Statistics", 
		"Allele Size Variation"
	];
	const TabPanelData = [
		TotalAllelesHistogram,
		AlleleModeHistogram,
		LociStatsScatter,
		LociSizeBoxplots
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
			<Alert severity="info" variant="outlined">
				<Typography sx={{ fontSize: 14 }}>
					Loci annotations were not provided.
				</Typography>
			</Alert>
		)
	};

	// data for CDS analysis
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
			filename: "cds_analysis.tsv",
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
		 tableTitle="CDS Analysis" 
		 tableOptions={analysisTableOptions}
		 tableConditionalFormatting={LociAnalysisFormatting}
		>
		</DataTable>
	);

	const HeaderSummary = (
		<Typography 
			sx={{ 
				color: '#bb7944', 
				fontSize: 20 
			}}
		>
			{'Report Description'}
		</Typography>
	);

	// Need to change the root component attribute to avoid error related with
	// passing something to <p>
	const HeaderDescription = (
		<Typography component={'div'}> 
			<ReactMarkdown 
				children={markdown} 
				remarkPlugins={[remarkGfm]}>
			</ReactMarkdown>
		</Typography>
	);

	// Get message about config used for schema creation
	const creationConfigMessage = data.creationConfig;
	const creationAlert = (
		<Alert severity="info" variant="outlined">
			<Typography sx={{ fontSize: 14 }}>
				{creationConfigMessage}
			</Typography>
		</Alert>
	);

	// Get message about config used for schema evaluation
	const evaluationConfigMessage = data.evaluationConfig;
	const evaluationAlert = (
		<Alert severity="info" variant="outlined">
			<Typography sx={{ fontSize: 14 }}>
				{evaluationConfigMessage}
			</Typography>
		</Alert>
	);

	// Create Alert to tell users that the points/boxplots are clickable
	// and link to the loci reports
	let clickAlert = undefined;
	if (data.lociReports === 1) {
		clickAlert = (
			<Alert severity="info" variant="outlined">
				<Typography sx={{ fontSize: 14 }}>
					In the Locus Statistics and Allele Size 
					Variation plots, each point (locus) is 
					clickable and will open a page with more 
					details about it.
				</Typography>
			</Alert>
		)
	};

	return (
		<div style={{ marginTop: "10px", marginBottom: "10px" }}>
			<div>
				<AccordionMUI 
					summaryText={HeaderSummary}
					detailsText={HeaderDescription} 
					expanded={false} >
				</AccordionMUI>
			</div>
			<div style={{ marginTop: "30px"}}>
				<div style={{ marginBottom: "10px" }}>
					<Stack sx={{ width: '100%' }} spacing={0.5}>
						{creationAlert}
						{evaluationAlert}
					</Stack>
				</div>
				{summaryTable}
			</div>
			<div style={{ marginTop: "30px"}}>
				<div style={{ marginBottom: "10px" }}>
					{clickAlert}
				</div>
				<TabPanelMUI
					ContentTitles={TabPanelTitles}
					ContentData={TabPanelData}
				>
				</TabPanelMUI>
			</div>
			<div style={{ marginTop: "30px"}}>
				<div style={{ marginBottom: "10px" }}>
					{LociAnnotationsAlert}
				</div>
				{LociAnnotationsTable}
			</div>
			<div style={{ marginTop: "30px"}}>
				{LociAnalysisTable}
			</div>
		</div>
	  );
};

export default SchemaPage;
