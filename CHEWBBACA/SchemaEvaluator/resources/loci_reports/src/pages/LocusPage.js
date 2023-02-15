import { useState, useRef, useEffect, useLayoutEffect } from 'react';
import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import Box from '@mui/material/Box';
import Paper from '@mui/material/Paper';
import Tabs from '@mui/material/Tabs';
import Tab from '@mui/material/Tab';
//import Container from '@mui/material/Container';
// MSAViewer
//import MSAViewer from "react-msa-viewer";

// Phylocanvas
import PhylogeneticTree from "../components/PhylogeneticTree";

import { Element } from 'react-scroll';

// Material-UI ExpansionPanel components
import Accordion from '@mui/material/Accordion';
import AccordionSummary from '@mui/material/AccordionSummary';
import AccordionDetails from '@mui/material/AccordionDetails';
import Typography from '@mui/material/Typography'; 
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import Divider from '@mui/material/Divider';

// import Roboto font
// import '@fontsource/roboto/300.css';


function TabPanel(props) {
	const { children, value, index, ...other } = props;
  
	return (
	  <div
		role="tabpanel"
		hidden={value !== index}
		id={`simple-tabpanel-${index}`}
		aria-labelledby={`simple-tab-${index}`}
		{...other}
	  >
		{value === index && (
		  <Box sx={{ p: 3 }}>
			{children}
		  </Box>
		)}
	  </div>
	);
};


function a11yProps(index) {
	return {
		id: `simple-tab-${index}`,
		'aria-controls': `simple-tabpanel-${index}`,
	};
};


const LocusPage = () => {
	const [panel, setPanel] = useState(0);

	const handleChange = (event, newValue) => {
		// where is this value coming from?
		setPanel(newValue);
	};

	// get pre-computed data
	const data = window.preComputedDataInd;

	// data for Summary Data table
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
			filename: "schema:summary.tsv",
			separator: "\t"
		},
		filter: false,
		search: false,
		viewColumns: true,
		pagination: false,
	};

	// data for Panel A
	const xDataPanelA = data.lengths;
	const yDataPanelA = data.lengths;
	const plotDataPanelA = [
		{x: xDataPanelA,
		 y: yDataPanelA,
		 type: "histogram",
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
	const layoutPanelA = {bargroupgap: 0.05};
	const configPanelA = {
		toImageButtonOptions: 
			{format: 'svg',
			 filename: `${locusName}_AlleleSizes`,
			 height: 500,
			 width: 700,
			 scale: 1
		}
	};

	// data for Panel B
	const xDataPanelB = data.ids;
	const yDataPanelB = xDataPanelA;
	const plotDataPanelB = [
		{x: xDataPanelB,
		 y: yDataPanelB,
		 type: "scatter",
		 name: "Distribution of allele mode sizes per gene",
		 mode: "markers",
		 marker: {
			color: "#0570b0",
		}
	    }
	];
	const layoutPanelB = {
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

	// get data for Phylocanvas tree
	const phyloData = data.phylo.phylo_data
	console.log('LocusPage');

	// get DNA sequences
	const dnaSequences = data.dna.sequences
	const dnaText = dnaSequences.map((seq) => {
		const seqid = seq.name;
		const sequence = seq.sequence;
		const sequenceStr = `>${seqid}\n${sequence}\n`
		return <Typography id={seqid} style={{fontSize: 12}}>{sequenceStr}</Typography>
	})

	// get Protein sequences
	const proteinSequences = data.protein.sequences
	const proteinText = proteinSequences.map((seq) => {
		const seqid = seq.name;
		const sequence = seq.sequence;
		const sequenceStr = `>${seqid}\n${sequence}\n`
		return <Typography id={seqid} style={{fontSize: 12}}>{sequenceStr}</Typography>
	})


	return (
		<div style={{ marginTop: "40px"}}>
			<div style={{ marginTop: "40px"}}>
				<DataTable 
					tableData={summaryData} 
					tableTitle="Summary Data" 
					tableOptions={summaryTableOptions}
				>
				</DataTable>
			</div>
			<div style={{ marginTop: "40px"}}>
				<Box sx={{ width: "100%" }}>
					<Paper elevation={3}>
						<Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
							<Tabs 
								value={panel} 
								onChange={handleChange} 
								aria-label="basic tabs example" 
								variant="scrollable"
								scrollButtons={true}
								allowScrollButtonsMobile
							>
								<Tab label="Total Alleles" wrapped {...a11yProps(0)} />
								<Tab label="Allele Mode Size" wrapped {...a11yProps(1)} />
							</Tabs>
						</Box>
						<TabPanel value={panel} index={0}>
							<PlotlyPlot 
								plotData={plotDataPanelA}
								plotTitle={locusName}
								xaxisTitle="Sequence Size (bp)"
								yaxisTitle="Number of Alleles"
								layoutProps={layoutPanelA}
								configOptions={configPanelA}
							>
							</PlotlyPlot>
						</TabPanel>
						<TabPanel value={panel} index={1}>
							<PlotlyPlot 
								plotData={plotDataPanelB}
								plotTitle="Distribution of allele mode sizes"
								xaxisTitle="Allele Mode Size"
								yaxisTitle="Number of Loci"
								layoutProps={layoutPanelB}
								configOptions={configPanelB}
							>
							</PlotlyPlot>
						</TabPanel>
					</Paper>
				</Box>
			</div>
			<div style={{ marginTop: "40px" }}>
				<Accordion defaultExpanded>
					<AccordionSummary
						expandIcon={<ExpandMoreIcon />}
						aria-controls="panella-content"
						id="panella-header"
					>
						<Typography>Phylocanvas Tree</Typography>
					</AccordionSummary>
					<Divider />
					<AccordionDetails >
						<TabPanel>
							<Element 
								name="phyloTree" 
								className="element" 
								id="containerElement" style={{
								position: 'relative',
								height: '750px',
								overflow: 'scroll',
								marginBottom: '0px'
							}}>
								<div style={{ position: "relative", margin: "auto" }}>
									<div id="demo" style={{ position: "absolute", margin: "auto"}}>
										<PhylogeneticTree
											source={phyloData}
											treeWidth={600}
											treeHeight={700}
											showLabels
											showLeafLabels
											interactive
										>
										</PhylogeneticTree>
									</div>
								</div>
							</Element>
						</TabPanel>
					</AccordionDetails>
				</Accordion>
				<Accordion style={{ marginTop: "40px"}} defaultExpanded={false}>
					<AccordionSummary
						expandIcon={<ExpandMoreIcon />}
						aria-controls="panella-content"
						id="panella-header"
					>
						<Typography>DNA sequences</Typography>
					</AccordionSummary>
					<Divider />
					<AccordionDetails style={{overflowWrap: 'break-word'}}>
						{dnaText}
					</AccordionDetails>
				</Accordion>
				<Accordion defaultExpanded={false}>
					<AccordionSummary
						expandIcon={<ExpandMoreIcon />}
						aria-controls="panella-content"
						id="panella-header"
					>
						<Typography>Protein sequences</Typography>
					</AccordionSummary>
					<Divider />
					<AccordionDetails style={{overflowWrap: 'break-word'}}>
						{proteinText}
					</AccordionDetails>
				</Accordion>
			</div>
		</div>
	  );
};

export default LocusPage; 
