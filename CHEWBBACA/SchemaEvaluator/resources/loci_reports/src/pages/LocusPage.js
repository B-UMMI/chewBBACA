import { useState } from 'react';

import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import AccordionMUI from '../components/AccordionMUI';
import TabPanelMUI from '../components/TabPanelMUI';
import MSA from '../components/MSA';

// Material-UI components
import Box from '@mui/material/Box';
import Stack from '@mui/material/Stack';
import Alert from '@mui/material/Alert';
import Switch from '@mui/material/Switch';
import Select from '@mui/material/Select';
import MenuItem from '@mui/material/MenuItem';
import TextField from '@mui/material/TextField';
import InputLabel from '@mui/material/InputLabel';
import IconButton from '@mui/material/IconButton';
import FormControl from '@mui/material/FormControl';
import Autocomplete from '@mui/material/Autocomplete';
import FormControlLabel from '@mui/material/FormControlLabel';
import FileDownload from '@mui/icons-material/FileDownload';

// Phylocanvas
import { TreeTypes, Shapes } from "@phylocanvas/phylocanvas.gl";
import PhylogeneticTree from "../components/PhylogeneticTree";

// 
import { Element } from 'react-scroll';

// Monaco code editor (options example at https://monaco-react.surenatoyan.com/)
import Editor from "@monaco-editor/react";
import { Typography } from '@mui/material';

import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';


const LocusPage = () => {

	const markdown = `
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
    - A tree drawn with Phylocanvas based on the Neighbor Joining tree created by MAFFT.
  - **Multiple Sequence Alignment**
    - Visualization of the MSA computed by MAFFT.
  - **DNA sequences**
    - Text Editor in read-only mode with the alleles in the FASTA file.
  - **Protein sequences**
    - Text Editor in read-only mode with the translated alleles that were considered valid based on the configuration values.

  You can find more information in the SchemaEvaluator module [documentation](https://chewbbaca.readthedocs.io/en/latest/user/modules/SchemaEvaluator.html).  
  If you have any feature requests or issues to report, please visit chewBBACA's [GitHub repository](https://github.com/B-UMMI/chewBBACA).
  `;

	// Get pre-computed data
	const data = window.preComputedDataInd;
	const validIDs = data.validIDs;

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

	// Determine maximum allele length to define plot range
	const lengthMax = Math.max(...xDataPanelA);

	const [showBarThreshold, setShowBarThreshold] = useState(false);
	const [showScatterThreshold, setShowScatterThreshold] = useState(false);

	const handleShowBarThreshold = (event) => {
		setShowBarThreshold(event.target.checked)
	};

	const handleShowScatterThreshold = (event) => {
		setShowScatterThreshold(event.target.checked)
	};

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
		 plotData={plotDataPanelA}
		 layoutProps={layoutPanelA}
		 configOptions={configPanelA}
		 key="AlleleSizeBar"
		>
		</PlotlyPlot>
	);

	let xaxisScatterMin = showScatterThreshold ? botThreshold-5 : lengthMin
	let xaxisScatterMax = showScatterThreshold ? topThreshold+5 : lengthMax

	// Create array with point colors
	let colors = []
	for (const id of data.ids) {
		if (validIDs.includes(id)) {
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

	// Alerts to explain red lines added for bot and top length thresholds
	const lengthThresholdsBarAlert = (
		<Alert severity="info" key="lengthThresholdsAlert">
			<Typography sx={{ fontSize: 14 }}>
				The red lines represent the bot and top allele length thresholds.
				The bar corresponding to the length mode is colored in green.
			</Typography>
		</Alert>
	)

	const lengthThresholdsScatterAlert = (
		<Alert severity="info" key="lengthThresholdsAlert">
			<Typography sx={{ fontSize: 14 }}>
				The red lines represent the bot and top allele length thresholds.
			</Typography>
		</Alert>
	)

	const validAlleleColorScatterAlert = (
		<Alert severity="info" key="validAlleleColorAlert">
			<Typography sx={{ fontSize: 14 }}>
				Valid alleles are colored in blue and invalid alleles are colored in grey.
			</Typography>
		</Alert>
	)

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
			{lengthThresholdsScatterAlert}
			{validAlleleColorScatterAlert}
		</Stack>
	);

	const AlleleSizePanelTitles = ["Allele Size Counts", "Allele Size"];
	const AlleleSizePanelsData = [
		[lengthThresholdsBarAlert, lengthThresholdCheckBar, AlleleSizeBar],
		[scatterAlerts, lengthThresholdCheckScatter, AlleleSizeScatter]];

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
	const [showBranchLengths, setShowBranchLengths] = useState(false);
	const [selectedIds, setSelectedIds] = useState([]);
	const [treeNodeSize, setTreeNodeSize] = useState(16);
	const [treeFontSize, setTreeFontSize] = useState(16);
	const [treeNodeShape, setTreeNodeShape] = useState(Shapes.Circle);

	const handleTreeNodeShape = (event) => {
		setTreeNodeShape(event.target.value)
	}

	const handleTreeFontSize = (value) => {
		setTreeFontSize(+value)
	};

	const handleTreeNodeSize = (value) => {
		setTreeNodeSize(+value)
	};

	const handleTreeTypeSelect = (event) => {
		setTreeType(event.target.value);
	};

	const handleTreeLabelsCheck = (event) => {
		setShowLabels(event.target.checked);
	};

	const handleBranchLabelsCheck = (event) => {
		setShowBranchLengths(event.target.checked);
	};

	const handleChange = (event, newValue) => {
		setSelectedIds(newValue);
	};

	// Piece of state to trigger tree download
	const [download, setDownload] = useState(false);
	const [downloadType, setDownloadType] = useState("SVG");

	const handleClick = (event) => {
		setDownload(!download);
	}

	const handleExport = (event) => {
		setDownloadType(event.target.value);
	};

	const treeDownloadOptions = ["SVG", "Newick", "JSON"]
	const treeDownloadMenuOptions = treeDownloadOptions.map((format) => {
		return <MenuItem key={format} value={format}>{format}</MenuItem>
	})

	const treeDownloadMenu = (
		<Box key="tree-export-select" sx={{ p: 1 }}>
			<FormControl>
				<InputLabel id="tree-export-select">File Format</InputLabel>
				<Select 
					labelId="tree-export-select"
					id="tree-export"
					value={downloadType}
					label="Export format"
					onChange={handleExport}
					sx={{ height: 40, width: 110 }}
				>
					{treeDownloadMenuOptions}
				</Select>
			</FormControl>
			<IconButton
				onClick={handleClick}
				sx={{
					height: 40,
					width: 30,
					borderRadius: "4px",
					color: "rgb(255, 255, 255)",
					backgroundColor: "rgb(25, 118, 210)"
				}}
			>
				<FileDownload></FileDownload>
			</IconButton>
		</Box>
	);

	// Component to select tree nodes
	const treeIdsSelect = (
		<Box key="select-tree-id" sx={{ p: 1 }}>
			<Autocomplete
				limitTags={3}
				multiple
				disablePortal
				id="select-tree-id"
				options={validIDs}
				size="small"
				sx={{ height: 40, width: 300 }}
				renderInput={(params) => <TextField {...params} label="Select Ids" />}
				onChange={handleChange}
				clearOnEscape
			/>
		</Box>
	);

	const treeShapeOptions = Object.entries(TreeTypes).map((type) => {
		return <MenuItem key={type[0]} value={type[1]}>{type[0]}</MenuItem>
	})

	const treeNodeShapeOptions = Object.entries(Shapes).map((shape) => {
		return <MenuItem key={shape[0]} value={shape[1]}>{shape[0]}</MenuItem>
	})

	const phyloTreeMenu = (
		<Box
			key="tree-options-menu"
			sx={{
				display: 'flex',
				flexDirection: 'row',
				flexWrap: 'wrap',
				justifyContent: 'center',
				alignItems: 'center',
				alignContent: 'space-around'
			}}
		>
			<Box key="tree-type-select-box" sx={{ p: 1 }}>
				<FormControl>
					<InputLabel id="tree-type-select">Tree Type</InputLabel>
					<Select 
						labelId="tree-type-select"
						id="tree-type"
						value={treeType}
						label="Tree Type"
						onChange={handleTreeTypeSelect}
						sx={{ height: 40, width: 140 }}
					>
						{treeShapeOptions}
					</Select>
				</FormControl>
			</Box>
			<Box key="tree-node-shape-select" sx={{ p: 1 }}>
				<FormControl>
					<InputLabel id="tree-node-shape">Node Shape</InputLabel>
					<Select 
						labelId="tree-node-shape"
						id="node-shape"
						value={treeNodeShape}
						label="Node Shape"
						onChange={handleTreeNodeShape}
						sx={{ height: 40, width: 230 }}
					>
						{treeNodeShapeOptions}
					</Select>
				</FormControl>
			</Box>
			<Box
				component="form"
				noValidate
				autoComplete="off"
				sx={{ p: 1 }}
			>
				<TextField
					id="node-size"
					size="small"
					sx={{ height: 40, width: 80 }}
					label="Node Size"
					variant="outlined"
					defaultValue={treeNodeSize}
					onKeyPress={(e) => {
						if (e.key === 'Enter') {
							e.preventDefault();
							handleTreeNodeSize(e.target.value)
						}
					}}
				/>
			</Box>
			<Box
				component="form"
				noValidate
				autoComplete="off"
				sx={{ p: 1 }}
			>
				<TextField
					id="tree-font-size"
					size="small"
					sx={{ height: 40, width: 80 }}
					label="Font Size"
					variant="outlined"
					defaultValue={treeFontSize}
					onKeyPress={(e) => {
						if (e.key === 'Enter') {
							e.preventDefault();
							handleTreeFontSize(e.target.value)
						}
					}}
				/>
			</Box>
			{treeIdsSelect}
			{treeDownloadMenu}
			<Box key="node-labels-form" sx={{ p: 1 }}>
				<FormControlLabel 
					key="node-labels-form"
					control={<Switch key="node-labels-switch" defaultChecked size="medium" />} 
					label="Node Labels" 
					labelPlacement="start"
					onChange={handleTreeLabelsCheck}
					sx={{ height: 40 }}
				/>
			</Box>
			<Box key="branch-labels-form" sx={{ p: 1 }}>
				<FormControlLabel 
					key="branch-labels-form"
					control={<Switch key="branch-labels-switch" size="medium" />} 
					label="Branch Lengths" 
					labelPlacement="start"
					onChange={handleBranchLabelsCheck}
					sx={{ height: 40 }}
				/>
			</Box>
		</Box>
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
			<Box key="tree-box" sx={{ p: 1 }}>
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
						showBranchLengths={showBranchLengths}
						selectedIds={selectedIds}
						download={download}
						downloadType={downloadType}
						treeNodeSize={treeNodeSize}
						fontSize={treeFontSize}
						nodeShape={treeNodeShape}
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
					Nodes are labeled with the allele identifiers.
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

	const [colorScheme, setColorScheme] = useState("lesk");

	const handleColorChange = (event) => {
		setColorScheme(event.target.value)
	};

	const colorSchemes = [
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

	const colorSchemesOptions = colorSchemes.map((color) => {
		return <MenuItem key={color[0]} value={color[0]}>{color[1]}</MenuItem>
	})

	const [msaExport, setMSAExport] = useState(false);
	const handleMSAExport = (event) => {
		setMSAExport(!msaExport)
	};

	const [downloadMSAType, setDownloadMSAType] = useState("Full MSA");
	const handleMSATypeExport = (event) => {
		setDownloadMSAType(event.target.value);
	};

	const msaDownloadOptions = ["Full MSA", "Row Selection", "Image (PNG)"]
	const msaDownloadMenuOptions = msaDownloadOptions.map((format) => {
		return <MenuItem key={format} value={format}>{format}</MenuItem>
	})

	const [showConservation, setShowConservation] = useState(false);
	const handleMSAConservation = (event) => {
		setShowConservation(event.target.checked);
	};

	const [showSeqLogo, setShowSeqLogo] = useState(false);
	const handleSeqLogo = (event) => {
		setShowSeqLogo(event.target.checked);
	};

	// MSA component does not handle well simultaneous row and column selection
	// const [msaRow, setMSARow] = useState(undefined);
	// const handleMSARow = (value) => {
	// 	setMSARow(value);
	// };

	const [msaColumn, setMSAColumn] = useState(undefined);
	const handleMSAColumn = (value) => {
		setMSAColumn(value-1)
	};

	const [searchMotif, setSearchMotif] = useState(undefined);
	const handleSearchMotif = (value) => {
		setSearchMotif(value)
	};

	const exportMSA = (
		<Box key="msa-export-select" sx={{ p: 1 }}>
			<FormControl>
				<InputLabel id="msa-export-select">Export</InputLabel>
				<Select 
					labelId="msa-export-select"
					id="msa-export"
					value={downloadMSAType}
					label="Export"
					onChange={handleMSATypeExport}
					sx={{ height: 40, width: 160 }}
				>
					{msaDownloadMenuOptions}
				</Select>
			</FormControl>
			<IconButton
				onClick={handleMSAExport}
				sx={{
					height: 40,
					width: 30,
					borderRadius: "4px",
					color: "rgb(255, 255, 255)",
					backgroundColor: "rgb(25, 118, 210)"
				}}
			>
				<FileDownload></FileDownload>
			</IconButton>
		</Box>
	);

	const msaMenu = (
		<Box
			key="tree-options-menu"
			sx={{
				display: 'flex',
				flexDirection: 'row',
				flexWrap: 'wrap',
				justifyContent: 'center',
				alignItems: 'center',
				alignContent: 'space-around'
			}}
		>
			<Box key="msa-color-select-box" sx={{ p: 1 }}>
				<FormControl>
					<InputLabel id="msa-color-select">Color Scheme</InputLabel>
					<Select 
						labelId="msa-color-select"
						id="msa-color"
						value={colorScheme}
						label="Color Scheme"
						onChange={handleColorChange}
						sx={{ height: 40 , width: 190 }}
					>
						{colorSchemesOptions}
					</Select>
				</FormControl>
			</Box>
			{/* <Box
				component="form"
				noValidate
				autoComplete="off"
				sx={{ p: 1 }}
			>
				<TextField
					id="jump-to-row"
					size="small"
					sx={{ height: 40, width: 90 }}
					label="Row"
					variant="outlined"
					defaultValue={msaRow}
					onKeyPress={(e) => {
						if (e.key === 'Enter') {
							e.preventDefault();
							handleMSARow(e.target.value)
						}
					}}
				/>
			</Box> */}
			<Box
				component="form"
				noValidate
				autoComplete="off"
				sx={{ p: 1 }}
			>
				<TextField
					id="jump-to-col"
					size="small"
					sx={{ height: 40, width: 90 }}
					label="Column"
					variant="outlined"
					defaultValue={msaColumn}
					onKeyPress={(e) => {
						if (e.key === 'Enter') {
							e.preventDefault();
							handleMSAColumn(e.target.value)
						}
					}}
				/>
			</Box>
			<Box
				component="form"
				noValidate
				autoComplete="off"
				sx={{ p: 1 }}
			>
				<TextField
					id="search-motif"
					size="small"
					sx={{ height: 40, width: 130 }}
					label="Search Motif"
					variant="outlined"
					defaultValue={searchMotif}
					onKeyPress={(e) => {
						if (e.key === 'Enter') {
							e.preventDefault();
							handleSearchMotif(e.target.value)
						}
					}}
				/>
			</Box>
			<Box key="msa-export-all" sx={{ p: 1 }}>
				{exportMSA}
			</Box>
			<Box key="msa-conserv-switch" sx={{ p: 1 }}>
				<FormControlLabel 
					key="msa-conserv-switch"
					control={<Switch key="msa-conserv-switch" size="medium" />} 
					label="Conservation"
					labelPlacement="start"
					onChange={handleMSAConservation}
					sx={{ height: 40 }}
				/>
			</Box>
			<Box key="msa-seq-logo" sx={{ p: 1 }}>
				<FormControlLabel 
					key="msa-seq-logo"
					control={<Switch key="msa-seq-logo" size="medium" />} 
					label="Seq. Logo" 
					labelPlacement="start"
					onChange={handleSeqLogo}
					sx={{ height: 40 }}
				/>
			</Box>
		</Box>
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
			<Box key="msa-element">
				<MSA
					MSADAta={msaData}
					colorScheme={colorScheme}
					msaExport={msaExport}
					downloadMSAType={downloadMSAType}
					conservation={showConservation}
					seqLogo={showSeqLogo}
					// msaRow={msaRow}
					msaColumn={msaColumn}
					searchMotif={searchMotif}
				>
				</MSA>
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

	const HeaderSummary = (
		<Typography sx={{ color: '#bb7944', fontSize: 20 }}>
			Locus Report Description
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

	return (
		<div style={{ marginTop: "10px", marginBottom: "10px" }}>
			<div>
				<AccordionMUI 
					summaryText={HeaderSummary}
					detailsData={HeaderDescription} 
					expanded={false}
				>
				</AccordionMUI>
			</div>
			<div style={{ marginTop: "20px" }}>
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
				{MSAComponent}
			</div>
			<div style={{ marginTop: "20px" }}>
				{alertPhyloMSA}
				{PhylogeneticElement}
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
