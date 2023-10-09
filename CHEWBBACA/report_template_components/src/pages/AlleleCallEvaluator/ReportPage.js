import React from 'react'
import { useState } from 'react';

import Box from '@mui/material/Box';
import Select from '@mui/material/Select';
import AlertMUI from '../components/AlertMUI';
import MenuItem from '@mui/material/MenuItem';
import Grid from '@mui/material/Unstable_Grid2';
import DataTable from '../components/DataTable';
import TextField from '@mui/material/TextField';
import Typography from '@mui/material/Typography';
import IconButton from '@mui/material/IconButton';
import PlotlyPlot from '../components/PlotlyPlot';
import InputLabel from '@mui/material/InputLabel';
import FormControl from '@mui/material/FormControl';
import TabPanelMUI from '../components/TabPanelMUI';
import AccordionMUI from '../components/AccordionMUI';
import Autocomplete from '@mui/material/Autocomplete';
import KeyboardArrowLeftIcon from '@mui/icons-material/KeyboardArrowLeft';
import KeyboardDoubleArrowLeftIcon from '@mui/icons-material/KeyboardDoubleArrowLeft';
import KeyboardArrowRightIcon from '@mui/icons-material/KeyboardArrowRight';
import KeyboardDoubleArrowRightIcon from '@mui/icons-material/KeyboardDoubleArrowRight';

import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

import classes from './ReportPage.css';

import PhylogeneticTree from "../components/PhylogeneticTree";
import { alleleCallReportHeaderMessage, alleleCallReportAlertMessages, globalTableOptions, sampleTableHiddenColumns, lociTableHiddenColumns } from '../constants';


const ReportPage = () => {

	// get pre-computed data
	const data = window.preComputedData;

	// Define function to apply to column with locus identifiers
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

	// Define function to apply when column name includes "URL"
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

	// Create Alert components
	// Replace values in Alert Messages data
	alleleCallReportAlertMessages[0][0] = alleleCallReportAlertMessages[0][0].replaceAll("%VALUE%", 'samples')
	alleleCallReportAlertMessages[1][0] = alleleCallReportAlertMessages[1][0].replaceAll("%VALUE%", 'loci')
	// Get size of cgMLST
	const cgMLSTSize = data.cgMLST_size[0].cgMLST100
	alleleCallReportAlertMessages[6][0] = alleleCallReportAlertMessages[6][0].replaceAll("%VALUE%", cgMLSTSize)
	alleleCallReportAlertMessages[7][0] = alleleCallReportAlertMessages[7][0].replaceAll("%VALUE%", cgMLSTSize)

	// Create all Alert components
	const alertMessagesComponents = alleleCallReportAlertMessages.map((alert, index) => {
		return (
			<AlertMUI
				key={index}
				ind={index}
				severity={alert[1]}
				fontSize={14}
			>
				{alert[0]}
			</AlertMUI>
		)
	})

	// Create component to display report description
	const HeaderSummary = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Allele Call Report Description
		</Typography>
	);

	// Need to change the root component attribute to avoid error related with
	// passing something to <p>
	const HeaderDescription = (
		<Typography component={'div'} style={{ lineHeight: "18px" }} > 
			<ReactMarkdown children={alleleCallReportHeaderMessage} remarkPlugins={[remarkGfm]} />
		</Typography>
	);

	const HeaderComponent = (
		<AccordionMUI
			summaryText={HeaderSummary}
			detailsData={HeaderDescription} 
			expanded={false}
		>
		</AccordionMUI>
	);

	// Create component to display Summary Table
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
			tableTitle="Results Summary Data" 
			tableOptions={summaryTableOptions}
		>
		</DataTable>
	);

	// Create components to display stacked bar plots for sample and locus classification counts
	const classificationTypes = ['EXC', 'INF', 'PLOT3', 'PLOT5',
		'LOTSC', 'NIPH', 'NIPHEM', 'ALM',
		'ASM', 'PAMA', 'LNF'];
	const typesColors = ['#1a9850', '#a6d96a', '#d6604d', '#b2182b',
		'#e08214', '#8073ac', '#542788', '#fee090',
		'#fdb462', '#80b1d3', '#878787'];

	// Define bar plot step
	const [barStep, setBarStep] = useState(300); 

	// Data for sample classification counts stacked barplots
	const sampleIDs = data.sample_ids;
	const sampleClassCounts = data.sample_data;

	const [sampleBarIndex, setSampleBarIndex] = useState(0);

	let sampleBarMin = 0;
	let sampleBarMax = sampleIDs.length > barStep ? barStep : sampleIDs.length;
	// Define array with barplot ranges
	// Define number of ranges
	const maxRanges = Math.ceil(sampleIDs.length / barStep);
	let barplotRanges = [];
	for (let i = 0; i < maxRanges; i++) {
		barplotRanges.push([sampleBarMin, sampleBarMax])
		sampleBarMin = Math.min(sampleBarMin+barStep, sampleIDs.length-barStep);
		sampleBarMax = Math.min(sampleBarMax+barStep, sampleIDs.length)
	}

	const sampleSliderMax = Math.min(sampleIDs.length, barStep/2);

	let sampleTraces = [];
	for (let i = 0; i < classificationTypes.length; i++) {
		let classTrace = {
			x: sampleIDs.slice(barplotRanges[sampleBarIndex][0], barplotRanges[sampleBarIndex][1]),
			y: sampleClassCounts[i].slice(barplotRanges[sampleBarIndex][0], barplotRanges[sampleBarIndex][1]),
			type: 'bar',
			name: classificationTypes[i],
			marker: {
				color: typesColors[i],
			}
		}
		sampleTraces.push(classTrace);
	}

	const sampleTracesLayout = {
		title: {
			text: ""
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
		margin: {
			t: 15
		}
	};

	const sampleTracesConfig = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'SampleCounts',
			 height: 500,
			 width: 1500,
			 scale: 1
		}
	};

	// Component for Plotly bar chart with classification counts per locus
	const sampleCountsStackedBar = (
		<PlotlyPlot
			key="sampleStackedCounts"
			plotData={sampleTraces}
			layoutProps={sampleTracesLayout}
			configOptions={sampleTracesConfig}
		>
		</PlotlyPlot>
	);

	const handleSampleBarButtonClick = (event) => {
		if (event.target.id === 'minusSample') {
			setSampleBarIndex(Math.max(sampleBarIndex-1, 0))
		} else if (event.target.id === 'plusSample') {
			setSampleBarIndex(Math.min(sampleBarIndex+1, maxRanges-1))
		}
	}

	const handleSampleBarFast = (event) => {
		if (event.target.id === 'toSampleStart') {
			setSampleBarIndex(0)
		} else if (event.target.id === 'toSampleEnd') {
			setSampleBarIndex(maxRanges-1)
		}
	}

	// Component to select sample range
	const sampleButtonGroup = (
		<Box display="flex" alignItems="center" justifyContent="center">
			<IconButton onClick={handleSampleBarFast} size="small">
				<KeyboardDoubleArrowLeftIcon id="toSampleStart" fontSize="large"/>
			</IconButton>
			<IconButton onClick={handleSampleBarButtonClick} size="small">
				<KeyboardArrowLeftIcon id="minusSample" fontSize="large"/>
			</IconButton>
			<TextField
				id="sampleRange"
				size="small"
				InputProps={{ readOnly: true }}
				sx={{ width: 120, textAlign: "center" }}
				label={`${barplotRanges[sampleBarIndex][0]}-${barplotRanges[sampleBarIndex][1]}`}
				disabled
				InputLabelProps={{
					style: {
					  width: '100%',
					  color: 'black'
					}
				}}
			>
			</TextField>
			<IconButton onClick={handleSampleBarButtonClick} size="small">
				<KeyboardArrowRightIcon id="plusSample" fontSize="large"/>
			</IconButton>
			<IconButton onClick={handleSampleBarFast} size="small">
				<KeyboardDoubleArrowRightIcon id="toSampleEnd" fontSize="large"/>
			</IconButton>
		</Box>
	)

	// Data for loci classification counts stacked barplots
	const lociIDs = data.loci_ids;
	const lociClassCounts = data.loci_data;

	const [lociBarIndex, setLociBarIndex] = useState(0);

	let lociBarMin = 0;
	let lociBarMax = lociIDs.length > barStep ? barStep : lociIDs.length;
	// Define array with barplot ranges
	// Define number of ranges
	const lociMaxRanges = Math.ceil(lociIDs.length / barStep);
	let lociBarplotRanges = [];
	for (let i = 0; i < lociMaxRanges; i++) {
		lociBarplotRanges.push([lociBarMin, lociBarMax])
		lociBarMin = Math.min(lociBarMin+barStep, (lociIDs.length)-barStep);
		lociBarMax = Math.min(lociBarMax+barStep, lociIDs.length)
	}

	const lociSliderMax = Math.min(lociIDs.length, barStep/2);

	let lociTraces = [];
	for (let i = 0; i < classificationTypes.length; i++) {
		let classTrace = {
			x: lociIDs.slice(lociBarplotRanges[lociBarIndex][0], lociBarplotRanges[lociBarIndex][1]),
			y: lociClassCounts[i].slice(lociBarplotRanges[lociBarIndex][0], lociBarplotRanges[lociBarIndex][1]),
			type: 'bar',
			name: classificationTypes[i],
			marker: {
				color: typesColors[i],
			}
		}
		lociTraces.push(classTrace);
	}

	const lociTracesLayout = {
		title: {
			text: ""
		},
		xaxis: {
			title: { text: "Locus" },
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
		margin: {
			t:15
		}
	};

	const lociTracesConfig = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: 'LociCounts',
			 height: 500,
			 width: 1500,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with total alleles distribution
	const lociCountsStackedBar = (
		<PlotlyPlot
			key="lociStackedCounts"
			plotData={lociTraces}
			layoutProps={lociTracesLayout}
			configOptions={lociTracesConfig}
		>
		</PlotlyPlot>
	);

	const handleLociBarButtonClick = (event) => {
		if (event.target.id === 'minusLoci') {
			setLociBarIndex(Math.max(lociBarIndex-1, 0))
		} else if (event.target.id === 'plusLoci') {
			setLociBarIndex(Math.min(lociBarIndex+1, lociMaxRanges-1))
		}
	}

	const handleLociBarFast = (event) => {
		if (event.target.id === 'toLociStart') {
			setLociBarIndex(0)
		} else if (event.target.id === 'toLociEnd') {
			setLociBarIndex(lociMaxRanges-1)
		}
	}

	// ButtonGroup component to select loci range
	const lociButtonGroup = (
		<Box display="flex" alignItems="center" justifyContent="center">
			<IconButton onClick={handleLociBarFast} size="small">
				<KeyboardDoubleArrowLeftIcon id="toLociStart" fontSize="large"/>
			</IconButton>
			<IconButton onClick={handleLociBarButtonClick} size="small">
				<KeyboardArrowLeftIcon id="minusLoci" fontSize="large"/>
			</IconButton>
			<TextField
				id="lociRange"
				size="small"
				InputProps={{ readOnly: true }}
				sx={{ width: 120, textAlign: "center" }}
				label={`${lociBarplotRanges[lociBarIndex][0]}-${lociBarplotRanges[lociBarIndex][1]}`}
				disabled
				InputLabelProps={{
					style: {
					  width: '100%',
					  color: 'black'
					}
				}}
			>
			</TextField>
			<IconButton onClick={handleLociBarButtonClick} size="small">
				<KeyboardArrowRightIcon id="plusLoci" fontSize="large"/>
			</IconButton>
			<IconButton onClick={handleLociBarFast} size="small">
				<KeyboardDoubleArrowRightIcon id="toLociEnd" fontSize="large"/>
			</IconButton>
		</Box>
	)

	// Tab Panel for classification Counts
	const TabPanelTitles = [
		"Counts Per Sample", 
		"Counts Per Locus", 
	];

	const samplesTabPanelData = (
		<Box key="sample-box">
			<Box key="sample-info" sx={{ p: 1 }}>
				{alertMessagesComponents[0]}
			</Box>
			<Box key="sample-menu" sx={{ p: 1 }}>
				{sampleButtonGroup}
			</Box>
			<Box key="sample-data" sx={{ p: 1 }}>
				{sampleCountsStackedBar}
			</Box>
		</Box>
	)

	const lociTabPanelData = (
		<Box key="loci-box">
			<Box key="loci-info" sx={{ p: 1 }}>
				{alertMessagesComponents[1]}
			</Box>
			<Box key="loci-menu" sx={{ p: 1 }}>
				{lociButtonGroup}
			</Box>
			<Box key="loci-data" sx={{ p: 1 }}>
				{lociCountsStackedBar}
			</Box>
		</Box>
	)

	const TabPanelData = [
		[samplesTabPanelData],
		[lociTabPanelData],
	];

	const stackedBarsTabPanel = (
		<TabPanelMUI
			ContentTitles={TabPanelTitles}
			ContentData={TabPanelData}
		>
		</TabPanelMUI>
	)

	// Create components to display tables with sample and locus stats
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
		},
		elevation: 0
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
			hiddenColumns={sampleTableHiddenColumns}
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
		},
		elevation: 0
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
			hiddenColumns={lociTableHiddenColumns}
		>
		</DataTable>
	);

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
			<MenuItem key={index} value={index+1}>{item}</MenuItem>
		)
	})

	const sampleStatsMenu = (
		<Box sx={{ minWidth: 120, marginTop: "20px" }}>
			<FormControl fullWidth size='small'>
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
			<MenuItem key={index} value={index+1}>{item}</MenuItem>
		)
	})

	const lociStatsMenu = (
		<Box sx={{ minWidth: 120, marginTop: "20px" }}>
			<FormControl fullWidth size='small'>
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

	// Grid components
	const sampleStatsGrid = (
		<div key="sample-stats">
			{sampleTable}
			<div>
				{sampleStatsMenu}
				{sampleHistogram}
			</div>
		</div>
	)

	const lociStatsGrid = (
		<div key="loci-stats">
			{lociTable}
			<div>
				{lociStatsMenu}
				{lociHistogram}
			</div>
		</div>
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

	const statsTabPanel = (
		<TabPanelMUI
			ContentTitles={StatsTabPanelTitles}
			ContentData={StatsTabPanelData}
		>
		</TabPanelMUI>
	)

	// Create component to display table with loci annotations
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

	// Create component to display Loci Presence-Absence Heatmap
	const paRows = data.presence_absence[0].rows;
	const paLociIDs = data.presence_absence[1].loci_ids;
	const paSampleIDs = data.presence_absence[2].sample_ids;

	const colorscaleValue = [
		[0, '#f7f7f7'],
		[1, '#053061']
	];

	let paHeatmap = undefined;
	if (paRows.length > 0) {
		const paData = [
			{
				z: paRows,
				x: paLociIDs,
				y: paSampleIDs,
				type: 'heatmap',
				showscale: false,
				colorscale: colorscaleValue,
				xgap: 0.05,
				ygap: 0.05,
				// Necessary to force binary colorscale when values are all the same
				zmax: 1,
				zmin: 0
			}
		];

		const paLayout = {
			title: {
				text: ''
			},
			xaxis: {
				title: { text: "Locus" },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			yaxis: {
				title: { text: "Sample" },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			height: 700,
			margin: {
				t: 20,
				r: 20
			},
		}

		const paConfig = {
			toImageButtonOptions: 
				{format: 'svg',
				filename: 'presence_absence',
				height: 700,
				scale: 1
			}
		};

		paHeatmap = (
			<PlotlyPlot
				key="paHeatmap"
				plotData={paData}
				layoutProps={paLayout}
				configOptions={paConfig}
			>
			</PlotlyPlot>
		);
	}

	// Component to show the Presence-Absence heatmap for a single sample
	const [selectedSample, setSelectedSample] = useState(undefined);

	const handleSampleSelect = (event, value) => {
		setSelectedSample(value)
	}

	let paSampleHeatmap = undefined;
	if (paRows.length > 0 && paSampleIDs.includes(selectedSample)) {
		let selectedSampleIndex = paSampleIDs.indexOf(selectedSample)
		let paSampleData = [
			{
				z: [paRows[selectedSampleIndex]],
				x: paLociIDs,
				y: [paSampleIDs[selectedSampleIndex]],
				type: 'heatmap',
				showscale: false,
				colorscale: colorscaleValue,
				xgap: 0.05,
				ygap: 0.05,
				// Necessary to force binary colorscale when values are all the same
				zmax: 1,
				zmin: 0
			}
		];

		const paSampleLayout = {
			title: {
				text: ''
			},
			xaxis: {
				title: { text: '' },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			yaxis: {
				title: { text: '' },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			height: 80,
			margin: {
				t: 20,
				b: 0,
				r: 20
			}
		}

		const paSampleConfig = {
			toImageButtonOptions: 
				{format: 'svg',
				 filename: `pa_${selectedSample}`,
				 height: 80,
				 scale: 1
				},
			displayModeBar: false
		};

		paSampleHeatmap = (
			<PlotlyPlot
				key="paSampleHeatmap"
				plotData={paSampleData}
				layoutProps={paSampleLayout}
				configOptions={paSampleConfig}
			>
			</PlotlyPlot>
		);
	}

	// Create Menu to select single sample
	let sampleSelectMenu = undefined;
	if (paRows.length > 0) {
		sampleSelectMenu = (
			<Box key="select-sample" sx={{ p: 1 }}>
				<Autocomplete
					disablePortal
					id="sample-select"
					options={paSampleIDs}
					sx={{ width: 300 }}
					size='small'
					renderInput={(params) => <TextField {...params} label="Select Sample" />}
					onInputChange={handleSampleSelect}
				/>
			</Box>
		)
	}

	// Create component to display heatmap for a single locus
	const [selectedLocus, setSelectedLocus] = useState(undefined);

	const handleLocusSelect = (event, value) => {
		setSelectedLocus(value)
	}

	let paLocusHeatmap = undefined;
	if (paRows.length > 0 && paLociIDs.includes(selectedLocus)) {
		let selectedLocusIndex = paLociIDs.indexOf(selectedLocus)

		// Get locus values for all samples
		let locusData = []
		for (let item of paRows) {
			locusData.push([item[selectedLocusIndex]])
		}

		let paLocusData = [
			{
				z: locusData,
				x: [paLociIDs[selectedLocusIndex]],
				y: paSampleIDs,
				type: 'heatmap',
				showscale: false,
				colorscale: colorscaleValue,
				xgap: 0.05,
				ygap: 0.05,
				// Necessary to force binary colorscale when values are all the same
				zmax: 1,
				zmin: 0
			}
		];

		const paLocusLayout = {
			title: {
				text: ''
			},
			xaxis: {
				title: { text: '' },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			yaxis: {
				title: { text: '' },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			height: 700,
			width: 80,
			margin: {
				t: 20,
				l: 0,
				r: 20
			}
		}

		const paLocusConfig = {
			toImageButtonOptions: 
				{format: 'svg',
				filename: `pa_${selectedLocus}`,
				height: 700,
				scale: 1
			},
			displayModeBar: false
		};

		paLocusHeatmap = (
			<PlotlyPlot
				key="paLocusHeatmap"
				plotData={paLocusData}
				layoutProps={paLocusLayout}
				configOptions={paLocusConfig}
			>
			</PlotlyPlot>
		);
	}

	// Create Menu to select single locus
	let locusSelectMenu = undefined;
	if (paRows.length > 0) {
		locusSelectMenu = (
			<Box key="select-locus" sx={{ p: 1 }}>
				<Autocomplete
					disablePortal
					id="locus-select"
					options={paLociIDs}
					sx={{ width: 300 }}
					size='small'
					renderInput={(params) => <TextField {...params} label="Select Locus" />}
					onInputChange={handleLocusSelect}
				/>
			</Box>
		)
	}

	let paComponent = undefined;
	if (paRows.length > 0) {
		const paTitle = (
			<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
				Loci Presence-Absence
			</Typography>
		);

		// Menu component
		const paMenu = (
			<Box
				key="pa-options-menu"
				sx={{
					display: 'flex',
					flexDirection: 'row',
					flexWrap: 'wrap',
					justifyContent: 'center',
					alignItems: 'center',
					alignContent: 'space-around',
				}}
			>
				{sampleSelectMenu}
				{locusSelectMenu}
			</Box>
		);

		const paHeatmaps = (
			<Grid
				key="pa-heatmaps"
				container
				spacing={0}
			>
				<Grid key="sample-heatmap" xs={11}>
					{paSampleHeatmap}
				</Grid>
				<Grid key="sample-fill" xs={1}>
				</Grid>
				<Grid key="main-heatmap" xs={11}>
					{paHeatmap}
				</Grid>
				<Grid key="locus-heatmap" xs={1}>
					{paLocusHeatmap}
				</Grid>
			</Grid>
		);

		paComponent = (
			<AccordionMUI
				summaryText={paTitle}
				detailsData={[paMenu, paHeatmaps]}
				expanded={false}
				alerts={[alertMessagesComponents[3]]}
			>
			</AccordionMUI>
		);
	}

	// Create component to display Distance Matrix Heatmap
	const dmRows = data.distance_matrix[0].rows;
	const dmSampleIDs = data.distance_matrix[1].sample_ids;

	let dmHeatmap = undefined;
	if (dmRows.length > 0) {
		const dmData = [
			{
				z: dmRows,
				x: dmSampleIDs,
				y: dmSampleIDs,
				type: 'heatmap',
				showscale: true,
				colorscale: 'Viridis',
				xgap: 0.05,
				ygap: 0.05
			}
		];

		const dmLayout = {
			title: {
				text: ''
			},
			xaxis: {
				title: { text: "Sample" },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			yaxis: {
				title: { text: "Sample" },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			height: 700,
			margin: {
				t: 20,
			}
		}

		const dmConfig = {
			toImageButtonOptions: 
				{format: 'svg',
				filename: 'distance_matrix',
				height: 700,
				scale: 1
			}
		};

		dmHeatmap = (
			<PlotlyPlot
				key="dmHeatmap"
				plotData={dmData}
				layoutProps={dmLayout}
				configOptions={dmConfig}
			>
			</PlotlyPlot>
		);
	}

	// Component to show a heatmap for a single sample
	const [selectedDmSample, setSelectedDmSample] = useState(undefined);

	const handleDmSampleSelect = (event, value) => {
		setSelectedDmSample(value)
	}

	let dmSampleHeatmap = undefined;
	if (dmRows.length > 0 && dmSampleIDs.includes(selectedDmSample)) {
		let selectedDmSampleIndex = dmSampleIDs.indexOf(selectedDmSample)
		let dmSampleData = [
			{
				z: [dmRows[selectedDmSampleIndex]],
				x: dmSampleIDs,
				y: [dmSampleIDs[selectedDmSampleIndex]],
				type: 'heatmap',
				showscale: false,
				colorscale: 'Viridis',
				xgap: 0.05,
				ygap: 0.05,
			}
		];

		const dmSampleLayout = {
			title: {
				text: ''
			},
			xaxis: {
				title: { text: '' },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			yaxis: {
				title: { text: '' },
				ticks: '',
				showticklabels: false,
				showticks: false,
			},
			height: 80,
			margin: {
				t: 20,
				b: 0,
			}
		}

		const dmSampleConfig = {
			toImageButtonOptions: 
				{format: 'svg',
				filename: `dm_${selectedDmSample}`,
				height: 80,
				scale: 1
			},
			displayModeBar: false
		};

		dmSampleHeatmap = (
			<PlotlyPlot
				key="dmSampleHeatmap"
				plotData={dmSampleData}
				layoutProps={dmSampleLayout}
				configOptions={dmSampleConfig}
			>
			</PlotlyPlot>
		);
	}

	// Create Menu to select single sample
	let dmMenu = undefined;
	if (dmRows.length > 0) {
		const sampleDmSelectMenu = (
			<Autocomplete
				key="dm-sample-menu"
				disablePortal
				id="sample-select"
				options={dmSampleIDs}
				sx={{ width: 300 }}
				size='small'
				renderInput={(params) => <TextField {...params} label="Select Sample" />}
				onInputChange={handleDmSampleSelect}
			/>
		)

		// Upper Menu component
		dmMenu = (
			<Box
				key="dm-options-menu"
				sx={{
					display: 'flex',
					flexDirection: 'row',
					flexWrap: 'wrap',
					justifyContent: 'center',
					alignItems: 'center',
					alignContent: 'space-around'
				}}
			>
				{sampleDmSelectMenu}
			</Box>
		);
	}

	// Bottom Menu component to find closest strains
	const [sampleDistanceSelect, setSampleDistanceSelect] = useState(undefined);
	const [closestSamples, setClosestSamples] = useState(undefined);
	// Event handler used to change sample to compute distances
	const handleSampleDistanceSelect = (event, value) => {
		setSampleDistanceSelect(value)
	};

	// Component to select distance threshold
	// State used to set distance threshold value
	const [distanceValue, setDistanceValue] = useState(10);
	// Event handler used to change distance threshold value
	const handleDistanceValue = (value) => {
		setDistanceValue(+value)
	};

	let sampleClosestStrains = undefined;
	let sampleDistanceValue = undefined;
	if (dmRows.length > 0) {
	// Component to select sample
		sampleClosestStrains = (
			<Autocomplete
				disablePortal
				id="sample-distance-select"
				options={dmSampleIDs}
				sx={{ width: 300 }}
				size='small'
				renderInput={(params) => <TextField {...params} label="Select Sample" />}
				onInputChange={handleSampleDistanceSelect}
			/>
		)

		sampleDistanceValue = (
			<TextField
				id="distance-value"
				size="small"
				sx={{ height: 40, width: 100 }}
				label="Max Distance"
				variant="outlined"
				defaultValue={distanceValue}
				onKeyDown={(e) => {
					if (e.key === 'Enter') {
						e.preventDefault();
						handleDistanceValue(e.target.value)
					}
				}}
			/>
		)
	}

	// Find closest samples
	let selectedSamples = [{'columns': ['Sample', 'Distance']}, {'rows': []}];
	if (dmSampleIDs.includes(sampleDistanceSelect)) {
		// Find closest strains based on distance threshold
		const sampleIndex = dmSampleIDs.indexOf(sampleDistanceSelect)
		let sampleRow = dmRows[sampleIndex]
		for (let i = 0; i < dmSampleIDs.length; i++) {
			if (dmSampleIDs[i] !== sampleDistanceSelect) {
				if (sampleRow[i] <= distanceValue) {
					selectedSamples[1].rows.push([dmSampleIDs[i], sampleRow[i]])
				}
			}
		}
	}

	// Create table to display list of closest samples
	let distanceTable = undefined;
	if (sampleDistanceSelect) {
		const distanceTableCustomOptions = {
			downloadOptions: {
				filename: "closest_samples.tsv",
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
			},
			elevation: 0
		};
		const distanceTableOptions = {
			...globalTableOptions,
			...distanceTableCustomOptions,
		};
	
		// Component for Distance table
		distanceTable = (
			<DataTable
				tableData={selectedSamples} 
				tableTitle="Closest Samples" 
				tableOptions={distanceTableOptions}
			>
			</DataTable>
		);
	}

	let dmComponent = undefined;
	if (dmRows.length > 0) {
		const dmBottomMenu = (
			<Box key="dm-bot-menu">
				{alertMessagesComponents[5]}
				<Box
					key="dm-bottom-menu"
					sx={{
						display: 'flex',
						flexDirection: 'row',
						flexWrap: 'wrap',
						justifyContent: 'center',
						alignItems: 'center',
						alignContent: 'space-around',
						p: 1
					}}
				>
					<Box key="dm-sample" sx={{ p: 1 }}>
						{sampleClosestStrains}
					</Box>
					<Box key="dm-distance" sx={{ p: 1 }}>
						{sampleDistanceValue}
					</Box>
				</Box>
			</Box>
		);

		const dmHeatmaps = (
			<Grid
				key="dm-heatmaps"
				container
				spacing={0.5}
			>
				<Grid key="dm-sample-heatmap" xs={12}>
					{dmSampleHeatmap}
				</Grid>
				<Grid key="dm-main-heatmap" xs={12}>
					{dmHeatmap}
				</Grid>
			</Grid>
		);

		const dmTitle = (
			<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
				Allelic Distances
			</Typography>
		);

		dmComponent = (
			<AccordionMUI
				summaryText={dmTitle}
				detailsData={[dmMenu, dmHeatmaps, dmBottomMenu, distanceTable]}
				expanded={false}
				alerts={[alertMessagesComponents[4], alertMessagesComponents[6]]}
			>
			</AccordionMUI>
		);
	}

	// Create component to display NJ tree
	// Get data for Phylocanvas tree
	const phyloData = data.cgMLST_tree.phylo_data;

	// Define title for tree component
	const phylogeneticElementTitle = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Core-genome Neighbor-Joining Tree
		</Typography>
	);

	let phylogeneticElementTree = undefined;
	if (phyloData.length > 0) {
		phylogeneticElementTree = (
			<Box key="tree-box" sx={{ p: 1 }}>
				<PhylogeneticTree
					treeSource={phyloData}
					validIDs={sampleIDs}
				>
				</PhylogeneticTree>
			</Box>
		)
	};

	// Create component to display tree
	let PhylogeneticElement = undefined;
	if (phylogeneticElementTree !== undefined) {
		PhylogeneticElement = (
			<AccordionMUI
				summaryText={phylogeneticElementTitle}
				detailsData={[phylogeneticElementTree]}
				expanded={false}
				alerts={[alertMessagesComponents[7]]}
			>
			</AccordionMUI>
		);
	}

	return (
		<div className={classes.mainDiv}>
			<div>
				{HeaderComponent}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{summaryTable}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{stackedBarsTabPanel}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{statsTabPanel}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{LociAnnotationsTable ? undefined : alertMessagesComponents[2]}
				{LociAnnotationsTable}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{paComponent}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{dmComponent}
			</div>
			<div className={classes.secondaryDiv}>
			{/* <div style={{marginTop: "20px"}}> */}
				{PhylogeneticElement}
			</div>
		</div>
	  );
};

export default ReportPage;
