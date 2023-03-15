import { useState } from 'react';

import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import Resized from '../components/Resized';
import AccordionMUI from '../components/AccordionMUI';
import TabPanelMUI from '../components/TabPanelMUI';
import MSA from '../components/MSA';

// Material-UI components
import Box from '@mui/material/Box';
import Alert from '@mui/material/Alert';
import Switch from '@mui/material/Switch';
import Select from '@mui/material/Select';
import MenuItem from '@mui/material/MenuItem';
import InputLabel from '@mui/material/InputLabel';
import FormControl from '@mui/material/FormControl';
import FormControlLabel from '@mui/material/FormControlLabel';

// Phylocanvas
import { TreeTypes } from "@phylocanvas/phylocanvas.gl";
import PhylogeneticTree from "../components/PhylogeneticTree";

// 
import { Element } from 'react-scroll';

// Monaco code editor (options example at https://monaco-react.surenatoyan.com/)
import Editor from "@monaco-editor/react";
import { Typography } from '@mui/material';


const LocusPage = () => {

	// Get pre-computed data
	const data = window.preComputedDataInd;

	// Data for Summary Data table
	const summaryData = data.summaryData;
	const locusName = summaryData[1].rows[0][0]
	const summaryTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: `${locusName}_summary.tsv`,
			separator: "\t"
		},
		filter: false,
		search: false,
		viewColumns: true,
		pagination: false,
	};
	// Component for Summary table
	const summaryTable = <DataTable
						  tableData={summaryData} 
						  tableTitle="Locus Summary Data" 
						  tableOptions={summaryTableOptions}
						 >
						 </DataTable>

	// Data for exceptions table
	const exceptionsData = data.invalidAlleles;
	const exceptionsTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
		selectableRowsOnClick: false,
		print: false,
		download: true,
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
		viewColumns: true,
		pagination: true,
		rowsPerPage: 10,
		rowsPerPageOptions: [10, 20, 40, 80, 100],
		jumpToPage: true,
		draggableColumns: {
			enabled: true
		}
	};
	// Component for Exceptions table
	let exceptionsTable = undefined;
	if (exceptionsData[1].rows.length > 0) {
		exceptionsTable = <DataTable
						   tableData={exceptionsData} 
						   tableTitle="Invalid Alleles" 
						   tableOptions={exceptionsTableOptions}
						  >
						  </DataTable>
	}

	let ExceptionsAlert = undefined;
	if (exceptionsTable === undefined) {
		ExceptionsAlert = (
			<Alert severity="warning">
				<Typography sx={{ fontSize: 14 }}>
					List of invalid alleles is not displayed because all
					alleles are considered valid.
				</Typography>
			</Alert>
		)
	};

	// Data for Annotations Data table
	const annotationsData = data.annotations;
	const annotationsTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: `${locusName}_annotations.tsv`,
			separator: "\t"
		},
		filter: false,
		search: false,
		viewColumns: true,
		pagination: false,
	}

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
	const LociAnnotationsFormatting = {...uniprotLinkColumnFormatting};

	// Component for Annotations table
	let annotationsTable = undefined;
	if (annotationsData[0].columns.length > 0) {
		annotationsTable = <DataTable
						  tableData={annotationsData} 
						  tableTitle="Locus Annotation Data" 
						  tableOptions={annotationsTableOptions}
						  tableConditionalFormatting={LociAnnotationsFormatting}
						 >
						 </DataTable>
	};

	let LociAnnotationsAlert = undefined;
	if (annotationsTable === undefined) {
		LociAnnotationsAlert = (
			<Alert severity="warning">
				<Typography sx={{ fontSize: 14 }}>
					The locus annotations are not displayed because the
					user did not provide a file with annotations or there were
					no annotations for this locus.
				</Typography>
			</Alert>
		)
	};

	// Get bot and top thresholds to add shapes to plots
	const botThreshold = data.botThreshold;
	const topThreshold = data.topThreshold;

	// data for Panel A (Sequence Size Distribution)
	const xDataPanelA = data.counts[0];
	const yDataPanelA = data.counts[1];

	// Get mode value and its index in the array of xvalues for the bar plot
	const modeValue = summaryData[1].rows[0][9];
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
	let xaxisMin = lengthMin;
	if (lengthMin > botThreshold) {
		xaxisMin = botThreshold;
	};

	// Determine maximum allele length to define plot range
	const lengthMax = Math.max(...xDataPanelA);
	let xaxisMax = lengthMax;
	if (lengthMax < topThreshold) {
		xaxisMax = topThreshold;
	};

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
			range: [xaxisMin-20, xaxisMax+20]
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
		 plotData={plotDataPanelA}
		 layoutProps={layoutPanelA}
		 configOptions={configPanelA}
		 key="AlleleSizeBar"
		>
		</PlotlyPlot>
	);

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
			color: "#0570b0",
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
			range: [xaxisMin-20, xaxisMax+20]
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

	// Alert to explain red lines added for bot and top length thresholds
	const lengthThresholdsAlert = (
		<Alert severity="info" key="lengthThresholdsAlert">
			<Typography sx={{ fontSize: 14 }}>
				The red lines represent the bot and top allele length thresholds.
				The bar corresponding to the length mode is colored in green.
			</Typography>
		</Alert>
	)

	const AlleleSizePanelTitles = ["Allele Size Counts", "Allele Size"];
	const AlleleSizePanelsData = [
		[lengthThresholdsAlert, AlleleSizeBar],
		[lengthThresholdsAlert, AlleleSizeScatter]];

	const AlleleSizeTabs = (
		<TabPanelMUI
			ContentTitles={AlleleSizePanelTitles}
			ContentData={AlleleSizePanelsData}
		>
		</TabPanelMUI>
	);

	// get data for Phylocanvas tree
	const phyloData = data.phylo.phylo_data;

	const [treeSource, setTreeSource] = useState(phyloData);
	const [treeType, setTreeType] = useState(TreeTypes.Circular);
	const [showLabels, setShowLabels] = useState(true);

	const handleTreeTypeSelect = (event) => {
		setTreeType(event.target.value);
	};

	const handleTreeLabelsCheck = (event) => {
		setShowLabels(event.target.checked);
	};

	const phyloTreeMenu = (
		<div key="tree-options-menu" style={{display: 'flex', justifyContent: 'space-between'}}>
			<Box key="tree-type-select-box">
				<FormControl size="small">
					<InputLabel id="tree-type-select">Tree Type</InputLabel>
					<Select 
						labelId="tree-type-select"
						id="tree-type"
						value={treeType}
						label="Tree Type"
						onChange={handleTreeTypeSelect}
					>
						<MenuItem value={TreeTypes.Circular}>Circular</MenuItem>
						<MenuItem value={TreeTypes.Rectangular}>Rectangular</MenuItem>
						<MenuItem value={TreeTypes.Diagonal}>Diagonal</MenuItem>
						<MenuItem value={TreeTypes.Hierarchical}>Hierarchical</MenuItem>
						<MenuItem value={TreeTypes.Radial}>Radial</MenuItem>
					</Select>
				</FormControl>
			</Box>
			<FormControlLabel 
				key="tree-labels-form"
				control={<Switch key="tree-labels-switch" defaultChecked size="medium" />} 
				label="Leaf Labels" 
				labelPlacement="top"
				onChange={handleTreeLabelsCheck}
			/>
		</div>
	);

	// Define title for tree component
	const phylogeneticElementTitle = (
		<Typography sx={{ 
			color: '#bb7944', 
			fontSize: 20 
			}}
		>
			Neighbor-Joining Tree
		</Typography>

	);

	let phylogeneticElementTree = undefined;
	if (phyloData.length > 0) {
		phylogeneticElementTree = (
			<Box key="tree-box" sx={{ p: 3 }}>
				<Element 
				 name="phyloTree" 
				 className="element" 
				 id="containerElement"
				 style={{
					 position: 'relative',
					 height: '750px',
					 overflow: 'scroll',
					 marginBottom: '0px'
					 }}
				>
					<PhylogeneticTree
						treeSource={treeSource}
						treeType={treeType}
						showLabels={showLabels}
					>
					</PhylogeneticTree>
				</Element>
			</Box>
		);
	}

	// Alert to explain leaf node labels
	let alertLeafNodes = undefined;
	if (phylogeneticElementTree !== undefined) {
		alertLeafNodes = (
			<Alert severity="info" key="alertLeafNodes">
				<Typography sx={{ fontSize: 14 }}>
					Leaves are labeled with the allele identifiers.
					A tree might not be displayed if the distance between all alleles is 0.
				</Typography>
			</Alert>
		);
	}

	// Alert rendered when there is no Phylogenetic Tree and MSA
	let alertPhyloMSA = undefined;
	if (phylogeneticElementTree === undefined) {
		alertPhyloMSA = (
			<Alert variant="outlined" severity="warning">
				<Typography sx={{ fontSize: 14 }}>
					The NJ tree and MSA components are not displayed because this locus 
					does not have valid alleles or has only one valid allele or
					the --light flag was provided.
				</Typography>
			</Alert>
		);
	}

	// Create component to display tree
	let PhylogeneticElement = undefined;
	if (phylogeneticElementTree !== undefined) {
		PhylogeneticElement = (
			<AccordionMUI
				summaryText={phylogeneticElementTitle}
				detailsData={[phyloTreeMenu, phylogeneticElementTree]}
				expanded={false}
				alerts={[alertLeafNodes]}
			>
			</AccordionMUI>
		);
	}

	// get data for MSA
	const msaData = data.msa.sequences;

	const [colorScheme, setColorScheme] = useState("clustal");

	const handleColorChange = (event) => {
		setColorScheme(event.target.value)
	};

	const msaMenu = (
		<div key="msa-options-menu" style={{display: 'flex', justifyContent: 'space-between'}}>
			<Box key="msa-color-select-box">
				<FormControl size="small">
					<InputLabel id="msa-color-select">Color Scheme</InputLabel>
					<Select 
						labelId="msa-color-select"
						id="msa-color"
						value={colorScheme}
						label="Color Scheme"
						onChange={handleColorChange}
					>
						<MenuItem value={"buried_index"}>Buried</MenuItem>
						<MenuItem value={"cinema"}>Cinema</MenuItem>
						<MenuItem value={"clustal"}>Clustal</MenuItem>
						<MenuItem value={"clustal2"}>Clustal2</MenuItem>
						<MenuItem value={"helix_propensity"}>Helix propensity</MenuItem>
						<MenuItem value={"hydro"}>Hydro</MenuItem>
						<MenuItem value={"lesk"}>Lesk</MenuItem>
						<MenuItem value={"mae"}>Mae</MenuItem>
						<MenuItem value={"nucleotide"}>Nucleotide</MenuItem>
						<MenuItem value={"purine_pyrimidine"}>Purine Pyrimidine</MenuItem>
						<MenuItem value={"strand_propensity"}>Strand propensity</MenuItem>
						<MenuItem value={"taylor"}>Taylor</MenuItem>
						<MenuItem value={"turn_propensity"}></MenuItem>
						<MenuItem value={"zappo"}>Zappo</MenuItem>
					</Select>
				</FormControl>
			</Box>
		</div>
	);

	// Define title for MSA component
	const MSAComponentTitle = (
		<Typography sx={{ 
			color: '#bb7944', 
			fontSize: 20 
			}}
		>
			Multiple Sequence Alignment
		</Typography>
	);

	let msaElement = undefined;
	if (msaData.length > 0) {
		// pass MSA component constructor instead of instance
		// this allows to get constructor and pass props in component that receives constructor
		msaElement = (
			<Box key="msa-element" sx={{ p: 3 }}>
				<Resized
					divID="MSA"
					component={MSA}
					colorScheme={colorScheme}
				>
				</Resized>
			</Box>
		);
	}

	// Alert for MSA
	let alertMSA = undefined;
	if (msaElement !== undefined) {
		alertMSA = (
			<Alert key="alertMSA" severity="info">
				<Typography sx={{ fontSize: 14 }}>
					Displaying the MSA for the {summaryData[1].rows[0][2]} alleles
					that were considered valid.
				</Typography>
			</Alert>
		);
	}

	// create component for MSA
	let MSAComponent = undefined;
	if (msaElement !== undefined) {
		MSAComponent = (
			<AccordionMUI
				summaryText={MSAComponentTitle}
				detailsData={[msaMenu, msaElement]}
				expanded={false}
				alerts={[alertMSA]}
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
			<Typography sx={{ 
				color: '#bb7944', 
				fontSize: 20 
				}}
			>
				DNA sequences
			</Typography>

		);

		const dnaEditorAlert = (
			<Alert severity="info" key="dnaEditorAlert">
				<Typography sx={{ fontSize: 14 }}>
					Displaying the {summaryData[1].rows[0][1]} alleles
					in the schema FASTA file.
				</Typography>
			</Alert>
		)

		DNAEditor = (
			<AccordionMUI
			 summaryText={DNAEditorTitle}
			 detailsData={DNAEditor}
			 expanded={false}
			 alerts={[dnaEditorAlert]}
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
			<Typography sx={{ 
				color: '#bb7944', 
				fontSize: 20 
				}}
			>
				Protein sequences
			</Typography>

		);

		const proteinEditorAlert = (
			<Alert severity="info" key="proteinEditorAlert">
				<Typography sx={{ fontSize: 14 }}>
					Only displaying the {summaryData[1].rows[0][2]} alleles
					that were considered valid.
				</Typography>
			</Alert>
		)

		ProteinEditor = (
			<AccordionMUI
			 summaryText={ProteinEditorTitle}
			 detailsData={ProteinEditor}
			 expanded={false}
			 alerts={[proteinEditorAlert]}
			>
			</AccordionMUI>
		);
	}

	return (
		<div style={{ marginTop: "10px", marginBottom: "10px" }}>
			<div style={{ marginTop: "10px" }}>
				{summaryTable}
			</div>
			<div style={{ marginTop: "20px" }}>
				{LociAnnotationsAlert}
				{annotationsTable}
			</div>
			<div style={{ marginTop: "20px"}}>
				{AlleleSizeTabs}
			</div>
			<div style={{ marginTop: "20px"}}>
				{ExceptionsAlert}
				{exceptionsTable}
			</div>
			<div style={{ marginTop: "20px" }}>
				{alertPhyloMSA}
				{PhylogeneticElement}
			</div>
			<div style={{ marginTop: "20px" }}>
				{MSAComponent}
			</div>
			<div style={{ marginTop: "20px"}}>
				{DNAEditor}
			</div>
			<div style={{ marginTop: "20px"}}>
				{ProteinEditor}
			</div>
		</div>
	  );
};

export default LocusPage; 
