import React from 'react'
import { useState } from 'react';

import { headerMessage, alertMessages, globalTableOptions } from '../constants';

import classes from './ReportPage.css';

import Box from '@mui/material/Box';
import Grid from '@mui/material/Unstable_Grid2';
import Button from '@mui/material/Button';
import ButtonGroup from '@mui/material/ButtonGroup';
import AlertMUI from '../components/AlertMUI';
import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import TabPanelMUI from '../components/TabPanelMUI';
import AccordionMUI from '../components/AccordionMUI';
import InputLabel from '@mui/material/InputLabel';
import MenuItem from '@mui/material/MenuItem';
import FormControl from '@mui/material/FormControl';
import Select, { SelectChangeEvent } from '@mui/material/Select';

import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';

import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';


const ReportPage = () => {

	// Define bar plot step
	const [barStep, setBarStep] = useState(300); 

	// get pre-computed data
	const data = window.preComputedData;

	// data for Summary Data table
	const summaryData = data.summaryData;

	const summaryTableCustomOption = {
		downloadOptions: {
			filename: "report_summary.tsv",
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
			tableTitle="Report Summary Data" 
			tableOptions={summaryTableOptions}
		>
		</DataTable>
	);

	// Data for table with sample stats
	const sampleStats = data.sample_stats

	// Custom option for Sample stats Table
	const sampleTableCustomOptions = {
		downloadOptions: {
			filename: "sample_stats.tsv",
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
	const sampleTableOptions = {
		...globalTableOptions,
		...sampleTableCustomOptions,
	}
	// Component for Sample stats table
	const sampleTable = (
		<DataTable
			tableData={sampleStats} 
			tableTitle="Sample Stats" 
			tableOptions={sampleTableOptions}
		>
		</DataTable>
	);

	// Data for table with loci stats
	const lociStats = data.loci_stats

	// Custom option for Loci stats Table
	const lociTableCustomOptions = {
		downloadOptions: {
			filename: "loci_stats.tsv",
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
	const lociTableOptions = {
		...globalTableOptions,
		...lociTableCustomOptions,
	}
	// Component for Loci stats table
	const lociTable = (
		<DataTable
			tableData={lociStats} 
			tableTitle="Loci Stats" 
			tableOptions={lociTableOptions}
		>
		</DataTable>
	);

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

	// Component for histogram with sample stats
	const [sampleColumn, setSampleColumn] = useState(1);

	const handleSampleColumnChange = (event) => {
		setSampleColumn(event.target.value);
	}

	// Get column name
	let sampleColumnName = sampleStats[0].columns[sampleColumn];

	let xSampleData = [];
	let ySampleData = [];
	for (let item of sampleStats[1].rows) {
		xSampleData.push(item[sampleColumn])
		ySampleData.push(item[0])
	}

	const plotSample = [
		{x: xSampleData,
		 y: ySampleData,
		 type: "histogram",
		 name: `Distribution of ${sampleColumnName}`,
		 marker: {
			 color: "#0570b0",
			 line: {
				 color: "#a6bddb",
				 width: 1
			 }
		 }
	    }
	];
	const layoutSample = {
		title: {
			text: `Distribution of ${sampleColumnName}`
		},
		xaxis: {
			title: { text: sampleColumnName },
			showgrid: false,
			zeroline: false,
			showline: true,
			showticklabels: false,
		},
		yaxis: {
			title: { text: "Count" },
			showgrid: true,
			zeroline: false,
			showline: true,
			ticks: "outside"
			//dtick: 100
		},
		bargroupgap: 0,
		bargap: 0,
	};
	const configSample = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: sampleColumnName,
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with total alleles distribution
	const sampleHistogram = (
		<PlotlyPlot
			key="sampleHistogram"
			plotData={plotSample}
			layoutProps={layoutSample}
			configOptions={configSample}
		>
		</PlotlyPlot>
	);

	// Create Select Menu to select stat to plot
	const sampleStatsMenuOptions = (sampleStats[0].columns.slice(1,)).map((item, index) => {
		return (
			<MenuItem value={index+1}>{item}</MenuItem>
		)
	})

	const sampleStatsMenu = (
		<Box sx={{ minWidth: 120 }}>
			<FormControl fullWidth>
				<InputLabel id="sample-menu">Select Stat</InputLabel>
				<Select
					labelId="sample-select"
					id="sample-select"
					value={sampleColumn}
					label={sampleColumnName}
					onChange={handleSampleColumnChange}
				>
					{sampleStatsMenuOptions}
				</Select>
			</FormControl>
		</Box>
	)

	// Component for histogram with loci stats
	const [lociColumn, setLociColumn] = useState(1);

	// Get column name
	let lociColumnName = lociStats[0].columns[lociColumn];

	const handleLociColumnChange = (event) => {
		setLociColumn(event.target.value);
	}

	let xLociData = [];
	let yLociData = [];
	for (let item of lociStats[1].rows) {
		xLociData.push(item[lociColumn])
		yLociData.push(item[0])
	}

	const plotLoci = [
		{x: xLociData,
		 y: yLociData,
		 type: "histogram",
		 name: `Distribution of ${lociColumnName}`,
		 marker: {
			 color: "#0570b0",
			 line: {
				 color: "#a6bddb",
				 width: 1
			 }
		 }
	    }
	];
	const layoutLoci = {
		title: {
			text: `Distribution of ${lociColumnName}`
		},
		xaxis: {
			title: { text: lociColumnName },
			showgrid: false,
			zeroline: false,
			showline: true,
			showticklabels: false,
		},
		yaxis: {
			title: { text: "Count" },
			showgrid: true,
			zeroline: false,
			showline: true,
			ticks: "outside"
			//dtick: 100
		},
		bargroupgap: 0,
		bargap: 0,
	};
	const configLoci = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: lociColumnName,
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with total alleles distribution
	const lociHistogram = (
		<PlotlyPlot
			key="lociHistogram"
			plotData={plotLoci}
			layoutProps={layoutLoci}
			configOptions={configLoci}
		>
		</PlotlyPlot>
	);

	// Create Select Menu to select stat to plot
	const lociStatsMenuOptions = (lociStats[0].columns.slice(1,)).map((item, index) => {
		return (
			<MenuItem value={index+1}>{item}</MenuItem>
		)
	})

	const lociStatsMenu = (
		<Box sx={{ minWidth: 120 }}>
			<FormControl fullWidth>
				<InputLabel id="loci-menu">Select Stat</InputLabel>
				<Select
					labelId="loci-select"
					id="loci-select"
					value={lociColumn}
					label={lociColumnName}
					onChange={handleLociColumnChange}
				>
					{lociStatsMenuOptions}
				</Select>
			</FormControl>
		</Box>
	)

	// Bar charts with class counts per sample/loci
	const classes = ['EXC', 'INF', 'PLOT3', 'PLOT5',
		'LOTSC', 'NIPH', 'NIPHEM', 'ALM',
		'ASM', 'PAMA', 'LNF'];
	const classesColors = ['#1a9850', '#a6d96a', '#d6604d', '#b2182b',
		'#e08214', '#8073ac', '#542788', '#fee090',
		'#fdb462', '#80b1d3', '#878787'];

	// data for sample classification counts stacked barplots
	const sampleIDs = data.sample_ids;
	const sampleClassCounts = data.sample_data;

	const [sampleBarMin, setSampleBarMin] = useState(0);
	const [sampleBarMax, setSampleBarMax] = useState(sampleIDs.length > barStep ? barStep : sampleIDs.length);

	const sampleSliderMax = Math.min(sampleIDs.length, barStep/2);

	let sampleTraces = [];
	for (let i = 0; i < classes.length; i++) {
		let classTrace = {
			x: sampleIDs.slice(sampleBarMin, sampleBarMax),
			y: sampleClassCounts[i].slice(sampleBarMin, sampleBarMax),
			type: 'bar',
			name: classes[i],
			marker: {
				color: classesColors[i],
			}
		}
		sampleTraces.push(classTrace);
	}

	const sampleTracesLayout = {
		title: {
			text: "Class Counts Per Sample"
		},
		xaxis: {
			title: { text: "Sample" },
			showgrid: false,
			zeroline: false,
			showline: true,
			showticks: false,
			showticklabels: false,
			rangeslider: {
				thickness: 0.15, 
			},
			range: [-0.5, sampleSliderMax]
		},
		yaxis: {
			title: { text: "Class Counts" },
			showgrid: true,
			zeroline: false,
			showline: true,
			showticks: false,
			//dtick: 100
		},
		//bargap: 0,
		bargroupgap: 0,
		barmode: 'stack',
		hovermode: 'x unified',
	};
	const sampleTracesConfig = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'SampleCounts',
			 height: 500,
			 width: 1000,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with total alleles distribution
	const sampleCountsStackedBar = (
		<PlotlyPlot
			key="TotalAllelesHistogram"
			plotData={sampleTraces}
			layoutProps={sampleTracesLayout}
			configOptions={sampleTracesConfig}
		>
		</PlotlyPlot>
	);

	const handleSampleBarButtonClick = (event) => {
		if (event.target.id === 'previous') {
			setSampleBarMin(Math.max(sampleBarMin-barStep, 0))
			setSampleBarMax(Math.max(sampleBarMax-barStep, 0+barStep))
		} else if (event.target.id === 'next') {
			setSampleBarMax(Math.min(sampleBarMax+barStep, sampleIDs.length))
			setSampleBarMin(Math.min(sampleBarMin+barStep, sampleIDs.length-barStep))
		}
	}

	// ButtonGroup component to select sample range
	const sampleButtonGroup = (
		<Box
			key="sample-counts-menu"
			sx={{
				display: 'flex',
				flexDirection: 'row',
				flexWrap: 'wrap',
				justifyContent: 'right',
				alignItems: 'right',
			}}
		>
			<ButtonGroup variant="contained" aria-label="outlined primary button group">
				<Button id="previous" onClick={handleSampleBarButtonClick}>Previous {barStep}</Button>
				<Button id="next" onClick={handleSampleBarButtonClick}>Next {barStep}</Button>
			</ButtonGroup>
		</Box>
	)

	// data for sample classification counts stacked barplots
	const lociIDs = data.loci_ids;
	const lociClassCounts = data.loci_data;

	const [lociBarMin, setLociBarMin] = useState(0);
	const [lociBarMax, setLociBarMax] = useState(lociIDs.length > barStep ? barStep : lociIDs.length);

	const lociSliderMax = Math.min(lociIDs.length, barStep/2);

	let lociTraces = [];
	for (let i = 0; i < classes.length; i++) {
		let classTrace = {
			x: lociIDs.slice(lociBarMin, lociBarMax),
			y: lociClassCounts[i].slice(lociBarMin, lociBarMax),
			type: 'bar',
			name: classes[i],
			marker: {
				color: classesColors[i],
			}
		}
		lociTraces.push(classTrace);
	}

	const lociTracesLayout = {
		title: {
			text: "Class Counts Per Locus"
		},
		xaxis: {
			title: { text: "Sample" },
			showgrid: false,
			zeroline: false,
			showline: true,
			showticks: false,
			showticklabels: false,
			rangeslider: {
				thickness: 0.15, 
			},
			range: [-0.5, lociSliderMax]
		},
		yaxis: {
			title: { text: "Class Counts" },
			showgrid: true,
			zeroline: false,
			showline: true,
			showticks: false,
			//dtick: 100
		},
		//bargap: 0,
		bargroupgap: 0,
		barmode: 'stack',
		hovermode: 'x unified',
	};
	const lociTracesConfig = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'LociCounts',
			 height: 500,
			 width: 1000,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with total alleles distribution
	const lociCountsStackedBar = (
		<PlotlyPlot
			key="TotalAllelesHistogram"
			plotData={lociTraces}
			layoutProps={lociTracesLayout}
			configOptions={lociTracesConfig}
		>
		</PlotlyPlot>
	);

	const handleLociBarButtonClick = (event) => {
		if (event.target.id === 'previous') {
			setLociBarMin(Math.max(lociBarMin-barStep, 0))
			setLociBarMax(Math.max(lociBarMax-barStep, 0+barStep))
		} else if (event.target.id === 'next') {
			setLociBarMax(Math.min(lociBarMax+barStep, lociIDs.length))
			setLociBarMin(Math.min(lociBarMin+barStep, lociIDs.length-barStep))
		}
	}

	// ButtonGroup component to select loci range
	const lociButtonGroup = (
		<Box
			key="loci-counts-menu"
			sx={{
				display: 'flex',
				flexDirection: 'row',
				flexWrap: 'wrap',
				justifyContent: 'right',
				alignItems: 'right',
			}}
		>
			<ButtonGroup variant="contained" aria-label="outlined primary button group">
				<Button id="previous" onClick={handleLociBarButtonClick}>Previous {barStep}</Button>
				<Button id="next" onClick={handleLociBarButtonClick}>Next {barStep}</Button>
			</ButtonGroup>
		</Box>
	)

	// Tab Panel for classification Counts
	const TabPanelTitles = [
		"Counts Per Sample", 
		"Counts Per Locus", 
	];
	const TabPanelData = [
		[sampleButtonGroup, sampleCountsStackedBar],
		[lociButtonGroup, lociCountsStackedBar],
	];

	// Grid components
	const sampleStatsGrid = (
		<Box sx={{ flexGrow: 1 }}>
			<Grid container spacing={2}>
				<Grid xs={7}>
					{sampleTable}
				</Grid>
				<Grid xs={5}>
					{sampleStatsMenu}
					{sampleHistogram}
				</Grid>
			</Grid>
		</Box>
	)

	const lociStatsGrid = (
		<Box sx={{ flexGrow: 1 }}>
			<Grid container spacing={2}>
				<Grid xs={7}>
					{lociTable}
				</Grid>
				<Grid xs={5}>
					{lociStatsMenu}
					{lociHistogram}
				</Grid>
			</Grid>
		</Box>
	)

	// Tab Panel for Sample and Loci stats
	const StatsTabPanelTitles = [
		"Sample Stats", 
		"Loci Stats", 
	];
	const StatsTabPanelData = [
		[sampleStatsGrid],
		[lociStatsGrid],
	];


	return (
		<div className={classes.mainDiv}>
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
				<TabPanelMUI
					ContentTitles={StatsTabPanelTitles}
					ContentData={StatsTabPanelData}
				>
				</TabPanelMUI>
			</div>
			<div className={classes.secondaryDiv}>
				{LociAnnotationsTable}
			</div>
		</div>
	  );
};

export default ReportPage;
