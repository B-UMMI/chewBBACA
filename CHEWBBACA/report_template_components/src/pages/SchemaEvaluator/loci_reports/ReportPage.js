import { useState } from 'react';

// Monaco code editor (options example at https://monaco-react.surenatoyan.com/)
import Editor from "@monaco-editor/react"

// Components to render markdown
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

import classes from './ReportPage.css';

// Material-UI components
import Box from '@mui/material/Box';
import Stack from '@mui/material/Stack';
import Switch from '@mui/material/Switch';
import Typography from '@mui/material/Typography';
import FormControlLabel from '@mui/material/FormControlLabel';

// Import Header Markdown text and array with alert messages
import { locusReportHeaderMessage, locusReportAlertMessages, globalTableOptions } from '../constants';

// Custom components
import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import AccordionMUI from '../components/AccordionMUI';
import TabPanelMUI from '../components/TabPanelMUI';
import AlertMUI from '../components/AlertMUI';
import MSA from '../components/MSA';
import PhylogeneticTree from "../components/PhylogeneticTree";


const LocusPage = () => {

	// Get pre-computed data
	const data = window.preComputedDataInd;
	const validIDs = data.validIDs;

	// Data for Summary Data table
	const summaryData = data.summaryData;
	const locusName = summaryData[1].rows[0][0]

	// Replace values in Alert Messages data
	locusReportAlertMessages[9][0] = locusReportAlertMessages[9][0].replace("%VALUE%", data.distinctAlleles[1].rows.length)
	locusReportAlertMessages[11][0] = locusReportAlertMessages[11][0].replace("%VALUE%", summaryData[1].rows[0][1])
	locusReportAlertMessages[12][0] = locusReportAlertMessages[12][0].replace("%VALUE%", summaryData[1].rows[0][2])
	// Create warning for loci that have novel alleles not submitted to Chewie-NS
	locusReportAlertMessages[13][0] = locusReportAlertMessages[13][0].replace("%VALUE%", data.nsAlleles.length)
	if (data.nsAlleles.length > 0) {
		const nsIDs = data.nsAlleles.join(', ')
		locusReportAlertMessages[13][0] = locusReportAlertMessages[13][0].replace("%VALUE2%", nsIDs)
	}

	// Create all Alert components
	const alertMessagesComponents = locusReportAlertMessages.map((alert, index) => {
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

	// Custom option for Summary Table
	const summaryTableCustomOptions = {
		downloadOptions: {
			filename: `${locusName}_summary.tsv`,
			separator: "\t"
		},
		pagination: false
	}
	const summaryTableOptions = {
		...globalTableOptions,
		...summaryTableCustomOptions,
	}
	// Component for Summary table
	const summaryTable = (
		<DataTable
			tableData={summaryData} 
			tableTitle="Locus Summary Data" 
			tableOptions={summaryTableOptions}
		>
		</DataTable>
	);

	// Data for exceptions table
	const exceptionsData = data.invalidAlleles;
	// Custom options for Exceptions Table
	const exceptionsTableCustomOptions = {
		downloadOptions: {
			filename: `${locusName}_invalidAlleles.tsv`,
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
	}
	const exceptionsTableOptions = {
		...globalTableOptions,
		...exceptionsTableCustomOptions,
	};
	// Component for Exceptions table
	let exceptionsTable = undefined;
	if (exceptionsData[1].rows.length > 0) {
		exceptionsTable = (
			<DataTable
				tableData={exceptionsData} 
				tableTitle="Invalid Alleles and Size Outliers" 
				tableOptions={exceptionsTableOptions}
			>
			</DataTable>
		)
	}

	// Data for Annotations Data table
	const annotationsData = data.annotations;
	// Custom options for Locus Annotations table
	const annotationsTableCustomOptions = {
		downloadOptions: {
			filename: `${locusName}_annotations.tsv`,
			separator: "\t"
		},
		pagination: false,
	}
	const annotationsTableOptions = {
		...globalTableOptions,
		...annotationsTableCustomOptions,
	}

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
	const LociAnnotationsFormatting = {...uniprotLinkColumnFormatting};

	// Component for Annotations table
	let annotationsTable = undefined;
	if (annotationsData[0].columns.length > 0) {
		annotationsTable = (
			<DataTable
				tableData={annotationsData} 
				tableTitle="Locus Annotation Data" 
				tableOptions={annotationsTableOptions}
				tableConditionalFormatting={LociAnnotationsFormatting}
			>
			</DataTable>
		)
	};

	// State for Plotly plots
	// State used to control threshold shape display in Bar plot
	const [showBarThreshold, setShowBarThreshold] = useState(false);
	// Event handler used to change showBarThreshold value
	const handleShowBarThreshold = (event) => {
		setShowBarThreshold(event.target.checked)
	};

	// State used to control threshold shape display in Scatter plot
	const [showScatterThreshold, setShowScatterThreshold] = useState(false);
	// Event handler used to change showScatterThreshold value
	const handleShowScatterThreshold = (event) => {
		setShowScatterThreshold(event.target.checked)
	};

	// Get bot and top thresholds to add shapes to plots
	const botThreshold = data.botThreshold;
	const topThreshold = data.topThreshold;

	// data for Panel A (Sequence Size Distribution)
	const xDataPanelA = data.counts[0];
	const yDataPanelA = data.counts[1];

	// Get mode value and its index in the array of xvalues for the bar plot
	const modeValue = summaryData[1].rows[0][13];
	const modeIndex = xDataPanelA.indexOf(modeValue);

	// Change color for mode bar
	let colorArray = Array(xDataPanelA.length).fill("#0570b0");
	colorArray[modeIndex] = "#41ab5d";

	const plotDataPanelA = [
		{x: xDataPanelA,
		 y: yDataPanelA,
		 type: "bar",
		 name: locusName,
		 marker: {
			 color: colorArray,
			 line: {
				 color: "#a6bddb",
				 width: 1
			 }
		 }
	    }
	];

	// Determine minimum allele length to define plot range
	const lengthMin = Math.min(...xDataPanelA);

	// Determine maximum allele length to define plot range
	const lengthMax = Math.max(...xDataPanelA);

	let xaxisBarMin = showBarThreshold ? botThreshold-5 : lengthMin
	let xaxisBarMax = showBarThreshold ? topThreshold+5 : lengthMax

	const layoutPanelA = {
		title: {
            text: locusName,
        },
        xaxis: {
        	title: { text: "Sequence Size (bp)" },
        	showgrid: false,
			zeroline: false,
			showline: true,
			ticks: "outside",
			range: [xaxisBarMin-5, xaxisBarMax+5]
        },
        yaxis: {
        	title: { text: "Number of Alleles" },
			zeroline: false,
			showline: true,
			ticks: "outside"
        },
		bargroupgap: 0.05,
		hovermode: "x",
		shapes: [
            {
              line: { color: "red", width: 1 },
              type: "line",
              x0: botThreshold,
              x1: botThreshold,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
            {
              line: { color: "red", width: 1 },
              type: "line",
              x0: topThreshold,
              x1: topThreshold,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
            {
              fillcolor: "red",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: botThreshold,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
            {
              fillcolor: "red",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: topThreshold,
              x1: topThreshold+botThreshold,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
          ],
	};
	const configPanelA = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: `${locusName}_AlleleSizeCounts`,
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with allele size distribution
	const AlleleSizeBar = (
		<PlotlyPlot
			key="AlleleSizeBar"
			plotData={plotDataPanelA}
			layoutProps={layoutPanelA}
			configOptions={configPanelA}
		>
		</PlotlyPlot>
	);

	let xaxisScatterMin = showScatterThreshold ? botThreshold-5 : lengthMin
	let xaxisScatterMax = showScatterThreshold ? topThreshold+5 : lengthMax

	// Create array with point colors
	let colors = []
	for (const id of data.ids) {
		// Valid identifiers are strings
		if (validIDs.includes(String(id))) {
			colors.push("#0570b0")
		} else {
			colors.push("#969696")
		}
	}

	// data for Panel B (Allele Size)
	const xDataPanelB = data.ids;
	const yDataPanelB = data.lengths;
	const plotDataPanelB = [
		{x: xDataPanelB,
		 y: yDataPanelB,
		 type: "scatter",
		 name: locusName,
		 mode: "markers",
		 marker: {
			color: colors,
		}
	    }
	];
	const layoutPanelB = {
		title: {
            text: locusName,
        },
        xaxis: {
        	title: { text: "Allele ID" },
        	showgrid: true,
			zeroline: false,
			showline: true,
			ticks: "outside"
        },
        yaxis: {
        	title: { text: "Sequence Size (bp)" },
			zeroline: false,
			showline: true,
			ticks: "outside",
			range: [xaxisScatterMin-5, xaxisScatterMax+5]
        },
		bargroupgap: 0.05,
		shapes: [
            {
              line: { color: "red", width: 1 },
              type: "line",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: botThreshold,
              y1: botThreshold,
              yref: "y",
            },
            {
              line: { color: "red", width: 1 },
              type: "line",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: topThreshold,
              y1: topThreshold,
              yref: "y",
            },
            {
              fillcolor: "red",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: 0,
              y1: botThreshold,
              yref: "y",
            },
            {
              fillcolor: "red",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: topThreshold,
              y1: topThreshold+botThreshold,
              yref: "y",
            },
          ],
	};
	const configPanelB = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: `${locusName}_AlleleSizes`,
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	const AlleleSizeScatter = (
		<PlotlyPlot
			plotData={plotDataPanelB}
			layoutProps={layoutPanelB}
			configOptions={configPanelB}
			key="AlleleSizeScatter"
		>
		</PlotlyPlot>
	);

	const dataPanelC = data.distinctAlleles[1].rows;
	let xDataPanelC = [];
	let yDataPanelC = [];
	for (const d of dataPanelC) {
		xDataPanelC.push(String(d[0]))
		yDataPanelC.push(d[1].length)
	}

	const xDataTickvals = Array(xDataPanelC.length).fill().map((v,i)=>i)

	const plotDataPanelC = [
		{x: xDataTickvals,
		 y: yDataPanelC,
		 type: "bar",
		 name: locusName,
		 marker: {
			 color: "#0570b0",
			 line: {
				 color: "#a6bddb",
				 width: 1
			 }
		 }
	    }
	];

	const layoutPanelC = {
		title: {
            text: locusName,
        },
        xaxis: {
        	title: { text: "Protein Allele ID" },
        	showgrid: false,
			zeroline: false,
			showline: true,
			showticklabels: false,
			ticktext: xDataPanelC,
			tickvals: xDataTickvals,
        },
        yaxis: {
        	title: { text: "Number of Distinct Alleles" },
			zeroline: false,
			showline: true,
			ticks: "outside",
        },
		bargroupgap: 0.05,
		hovermode: "x",
	};
	const configPanelC = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: `${locusName}_DistinctAlleles`,
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};
	// Component for Plotly Histogram with allele size distribution
	const DistinctAllelesBar = (
		<PlotlyPlot
			key="DistinctAlleles"
			plotData={plotDataPanelC}
			layoutProps={layoutPanelC}
			configOptions={configPanelC}
		>
		</PlotlyPlot>
	);

	// Components to resize plots to show length thresholds
	// Need to include checked prop or switch will change to disabled state when alternating tabs
	const lengthThresholdCheckBar = (
		<Box key="length-threshold-box-bar" sx={{ p: 1 }}>
			<FormControlLabel 
				key="length-threshold-form-bar"
				control={<Switch key="length-threshold-switch-bar" checked={showBarThreshold} size="medium" />} 
				label="Show Thresholds" 
				labelPlacement="start"
				onChange={handleShowBarThreshold}
				sx={{ height: 40 }}
			/>
		</Box>
	);

	const lengthThresholdCheckScatter = (
		<Box key="length-threshold-box-scatter" sx={{ p: 1 }}>
			<FormControlLabel 
				key="length-threshold-form-scatter"
				control={<Switch key="length-threshold-switch-scatter" checked={showScatterThreshold} size="medium" />} 
				label="Show Thresholds" 
				labelPlacement="start"
				onChange={handleShowScatterThreshold}
				sx={{ height: 40 }}
			/>
		</Box>
	);

	const scatterAlerts = (
		<Stack key="scatterStack" sx={{ width: '100%' }} spacing={0.5}>
			{alertMessagesComponents[3]}
			{alertMessagesComponents[4]}
		</Stack>
	);

	let panelCData = []
	if (data.distinctAlleles[1].rows.length > 0) {
		panelCData = [alertMessagesComponents[5], DistinctAllelesBar]
	} else {
		panelCData = [alertMessagesComponents[14]]
	}

	const AlleleSizePanelTitles = ["Allele Size Counts", "Allele Size", "Alleles Per Protein"];
	const AlleleSizePanelsData = [
		[alertMessagesComponents[2], lengthThresholdCheckBar, AlleleSizeBar],
		[scatterAlerts, lengthThresholdCheckScatter, AlleleSizeScatter],
		panelCData
	];

	const AlleleSizeTabs = (
		<TabPanelMUI
			ContentTitles={AlleleSizePanelTitles}
			ContentData={AlleleSizePanelsData}
		>
		</TabPanelMUI>
	);

	// get data for Phylocanvas tree
	const phyloData = data.phylo.phylo_data;

	// Define title for tree component
	const phylogeneticElementTitle = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Neighbor-Joining Tree
		</Typography>
	);

	let phylogeneticElementTree = undefined;
	if (phyloData.length > 0) {
		phylogeneticElementTree = (
			<Box key="tree-box" sx={{ p: 1 }}>
				<PhylogeneticTree
					treeSource={phyloData}
					validIDs={xDataPanelC}
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
				alerts={[alertMessagesComponents[6], alertMessagesComponents[7]]}
			>
			</AccordionMUI>
		);
	}

	// get data for MSA
	const msaData = data.msa.sequences;

	// Define title for MSA component
	const MSAComponentTitle = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Multiple Sequence Alignment
		</Typography>
	);

	let msaElement = undefined;
	if (msaData.length > 0) {
		// pass MSA component constructor instead of instance
		// this allows to get constructor and pass props in component that receives constructor
		msaElement = (
			<Box key="msa-element">
				<MSA
					MSADAta={msaData}
				>
				</MSA>
			</Box>
		);
	}

	// create component for MSA
	let MSAComponent = undefined;
	if (msaElement !== undefined) {
		MSAComponent = (
			<AccordionMUI
				summaryText={MSAComponentTitle}
				detailsData={[msaElement]}
				expanded={false}
				alerts={[alertMessagesComponents[9], alertMessagesComponents[10]]}
			>
			</AccordionMUI>
		);
	}

	// get DNA sequences
	const dnaSequences = data.dna.sequences
	let DNAEditor = undefined;
	if (dnaSequences.length > 0) {
		const dnaText = dnaSequences.map((seq) => {
			const seqid = seq.name;
			const sequence = seq.sequence;
			const sequenceStr = `>${seqid}\n${sequence}\n`
			return sequenceStr
		})

		const joinedDNA = dnaText.join('');
		// create DNA Editor component
		DNAEditor = (
			<Editor
				height="40vh"
				options={{"readOnly": true, "wordWrap": "on"}}
				defaultValue={`${joinedDNA}`}
			>
			</Editor>
		);

		const DNAEditorTitle = (
			<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
				DNA sequences
			</Typography>
		);

		DNAEditor = (
			<AccordionMUI
				summaryText={DNAEditorTitle}
				detailsData={DNAEditor}
				expanded={false}
				alerts={[alertMessagesComponents[11]]}
			>
			</AccordionMUI>
		);
	}

	// get Protein sequences
	const proteinSequences = data.protein.sequences
	let ProteinEditor = undefined;
	if (proteinSequences.length > 0) {
		const proteinText = proteinSequences.map((seq) => {
			const seqid = seq.name;
			const sequence = seq.sequence;
			const sequenceStr = `>${seqid}\n${sequence}\n`
			return sequenceStr
		})
		const joinedProtein = proteinText.join('');
		// create DNA Editor component
		ProteinEditor = (
			<Editor
				height="40vh"
				options={{"readOnly": true, "wordWrap": "on"}}
				defaultValue={`${joinedProtein}`}
			>
			</Editor>
		);

		const ProteinEditorTitle = (
			<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
				Protein sequences
			</Typography>
		);

		ProteinEditor = (
			<AccordionMUI
				summaryText={ProteinEditorTitle}
				detailsData={ProteinEditor}
				expanded={false}
				alerts={[alertMessagesComponents[12]]}
			>
			</AccordionMUI>
		);
	}

	const HeaderSummary = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Locus Report Description
		</Typography>
	);

	// Need to change the root component attribute to avoid error related with
	// passing something to <p>
	const HeaderDescription = (
		<Typography component={'div'} style={{ lineHeight: "18px" }} > 
			<ReactMarkdown children={locusReportHeaderMessage} remarkPlugins={[remarkGfm]} />
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

	// Data for distinct alleles table
	const distinctData = data.distinctAlleles;
	const distinctColumns = distinctData[0];
	const distinctRowsData = distinctData[1].rows;
	let distinctRows = {"rows": []}
	for (const d of distinctRowsData) {
		distinctRows.rows.push([d[0], d[1].length, d[1].join(', ')])
	}

	// Custom options for Exceptions Table
	const distinctTableCustomOptions = {
		downloadOptions: {
			filename: `${locusName}_distinctAlleles.tsv`,
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
	}
	const distinctTableOptions = {
		...globalTableOptions,
		...distinctTableCustomOptions,
	};
	// Component for distinct alleles table
	const distinctTableData = [distinctColumns, distinctRows]
	const distinctTable = (
		<DataTable
			tableData={distinctTableData} 
			tableTitle="Distinct Protein Alleles" 
			tableOptions={distinctTableOptions}
		>
		</DataTable>
	)

	return (
		<div className={classes.mainDiv}>
			<div>
				{HeaderComponent}
			</div>
			<div className={classes.secondaryDiv}>
				{summaryTable}
			</div>
			<div className={classes.secondaryDiv}>
				{annotationsTable ? undefined : alertMessagesComponents[1]}
				{annotationsTable}
			</div>
			<div className={classes.secondaryDiv}>
				{data.nsAlleles.length > 0 ? alertMessagesComponents[13] : undefined}
				{AlleleSizeTabs}
			</div>
			<div className={classes.secondaryDiv}>
				{distinctRowsData.length > 0 ? distinctTable: alertMessagesComponents[15]}
			</div>
			<div className={classes.secondaryDiv}>
				{exceptionsTable ? undefined : alertMessagesComponents[0]}
				{exceptionsTable}
			</div>
			<div className={classes.secondaryDiv}>
				{MSAComponent}
			</div>
			<div className={classes.secondaryDiv}>
				{phylogeneticElementTree ? undefined : alertMessagesComponents[8]}
				{PhylogeneticElement}
			</div>
			<div className={classes.secondaryDiv}>
				{DNAEditor}
			</div>
			<div className={classes.secondaryDiv}>
				{ProteinEditor}
			</div>
		</div>
	);
};

export default LocusPage;
