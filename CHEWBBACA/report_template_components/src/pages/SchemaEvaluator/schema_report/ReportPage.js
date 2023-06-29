import React from 'react'

import { schemaReportHeaderMessage, schemaReportAlertMessages, globalTableOptions } from '../constants';

import classes from './ReportPage.css';

import AlertMUI from '../components/AlertMUI';
import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import TabPanelMUI from '../components/TabPanelMUI';
import AccordionMUI from '../components/AccordionMUI';

import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';

import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';


const SchemaPage = () => {

	// get pre-computed data
	const data = window.preComputedData;

	// data for Summary Data table
	const summaryData = data.summaryData;

	// Replace values in Alert Messages data
	schemaReportAlertMessages.push([data.creationConfig, "info"]);
	schemaReportAlertMessages.push([data.evaluationConfig, "info"]);
	if ((data.nsSchema.length)>0) {
		schemaReportAlertMessages.push([data.nsSchema, "warning"])
	}

	// Create all Alert components
	const alertMessagesComponents = schemaReportAlertMessages.map((alert, index) => {
		return (
			<AlertMUI
				key={index}
				severity={alert[1]}
				fontSize={14}
			>
				{alert[0]}
			</AlertMUI>
		)
	})

	const summaryTableCustomOption = {
		downloadOptions: {
			filename: "schema_summary.tsv",
			separator: "\t",
			filterOptions: {
				useDisplayedColumnsOnly: true,
				useDisplayedRowsOnly: true
			}
		},
		pagination: false,
		draggableColumns: {
			enabled: true
		}
	}
	const summaryTableOptions = {
		...globalTableOptions,
		...summaryTableCustomOption,
	};

	// Component for Summary table
	const summaryTable = (
		<DataTable
			tableData={summaryData} 
			tableTitle="Schema Summary Data" 
			tableOptions={summaryTableOptions}
		>
		</DataTable>
	);

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
		scatterClickAlert = alertMessagesComponents[0];
		boxplotClickAlert = alertMessagesComponents[1];
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
		"URL": {
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
	const annotationsTableCustomOptions = {
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
		rowsPerPage: 10,
		rowsPerPageOptions: [10, 20, 40, 80, 100],
		jumpToPage: true,
		draggableColumns: {
			enabled: true
		}
	};
	const annotationsTableOptions = {
		...globalTableOptions,
		...annotationsTableCustomOptions,
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

	let LociAnnotationsAlert = alertMessagesComponents[2];

	// data for allele analysis
	const analysisData = data.analysis;
	const analysisTableCustomOptions = {
		selectableRowsHideCheckboxes: true,
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
		rowsPerPage: 10,
		rowsPerPageOptions: [10, 20, 40, 80, 100],
		jumpToPage: true,
		draggableColumns: {
			enabled: true
		},
	};
	const analysisTableOptions = {
		...globalTableOptions,
		...analysisTableCustomOptions,
	}

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
				children={schemaReportHeaderMessage} 
				remarkPlugins={[remarkGfm]}
			>
			</ReactMarkdown>
		</Typography>
	);

	// Get message about config used for schema creation
	const creationAlert = alertMessagesComponents[3];

	// Get message about config used for schema evaluation
	const evaluationAlert = alertMessagesComponents[4];

	// Get Alert if 
	let nsAlert = undefined;
	if ((data.nsSchema.length)>0) {
		nsAlert = alertMessagesComponents[5];
	}

	return (
		<div className={classes.mainDiv}>
			<div style={{ marginBottom: "5px" }}>
				<Stack sx={{ width: '100%' }} spacing={0.5}>
					{creationAlert}
					{evaluationAlert}
					{nsAlert}
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
			<div className={classes.secondaryDiv}>
				{summaryTable}
			</div>
			<div className={classes.secondaryDiv}>
				<TabPanelMUI
					ContentTitles={TabPanelTitles}
					ContentData={TabPanelData}
				>
				</TabPanelMUI>
			</div>
			<div className={classes.secondaryDiv}>
				<div style={{ marginBottom: "10px" }}>
					{LociAnnotationsTable ? undefined: LociAnnotationsAlert}
				</div>
				{LociAnnotationsTable}
			</div>
			<div className={classes.secondaryDiv}>
				{LociAnalysisTable}
			</div>
		</div>
	  );
};

export default SchemaPage;
