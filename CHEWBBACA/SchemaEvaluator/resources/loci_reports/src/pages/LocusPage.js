import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import Resized from '../components/Resized';
import AccordionMUI from '../components/AccordionMUI';
import TabPanelMUI from '../components/TabPanelMUI';
import MSA from '../components/MSA';

// Material-UI components
import Box from '@mui/material/Box';

// Phylocanvas
import PhylogeneticTree from "../components/PhylogeneticTree";

// 
import { Element } from 'react-scroll';

// Monaco code editor (options example at https://monaco-react.surenatoyan.com/)
import Editor from "@monaco-editor/react";


const LocusPage = () => {

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
			filename: "schema_summary.tsv",
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
						  tableTitle="Summary Data" 
						  tableOptions={summaryTableOptions}
						 >
						 </DataTable>

	// data for Panel A (Sequence Size Distribution)
	const xDataPanelA = data.counts[0];
	const yDataPanelA = data.counts[1];
	const plotDataPanelA = [
		{x: xDataPanelA,
		 y: yDataPanelA,
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
	const layoutPanelA = {
		title: {
            text: locusName,
        },
        xaxis: {
        	title: { text: "Sequence Size (bp)" },
        	showgrid: true,
        },
        yaxis: {
        	title: { text: "Number of Alleles" },
        },
		bargroupgap: 0.05
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
        },
        yaxis: {
        	title: { text: "Sequence Size (bp)" },
        },
		bargroupgap: 0.05
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
		>
		</PlotlyPlot>
	);

	const AlleleSizePanelTitles = ["Allele Size Counts", "Allele Size"];
	const AlleleSizePanelsData = [AlleleSizeBar, AlleleSizeScatter];

	// get data for Phylocanvas tree
	const phyloData = data.phylo.phylo_data;
	let PhylogeneticElement = undefined;
	if (phyloData.length > 0) {
		PhylogeneticElement = (
			<Box sx={{ p: 3 }}>
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
					<div id="demo" style={{ margin: "auto" }}>
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
				</Element>
			</Box>
		);

		PhylogeneticElement = (
			<AccordionMUI
				summaryText="Phylogenetic Tree"
				detailsData={PhylogeneticElement}
				expanded={false}
			>
			</AccordionMUI>
		);
	}

	// create component for MSA
	let MSAComponent = undefined;
	if (data.msa.sequences.length > 0) {
		// pass MSA component constructor instead of instance
		// this allows to get constructor and pass props in component that receives constructor
		MSAComponent = (
			<Box sx={{ p: 3 }}>
				<Resized
					divID="MSA"
					component={MSA}
				>
				</Resized>
			</Box>
		);

		MSAComponent = (
			<AccordionMUI
				summaryText="Multiple Sequence Alignment"
				detailsData={MSAComponent}
				expanded={false}
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

		DNAEditor = (
			<AccordionMUI
			 summaryText="DNA sequences"
			 detailsData={DNAEditor}
			 expanded={false}
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

		ProteinEditor = (
			<AccordionMUI
			 summaryText="Protein sequences"
			 detailsData={ProteinEditor}
			 expanded={false}
			>
			</AccordionMUI>
		);
	}

	return (
		<div style={{ marginTop: "40px" }}>
			<div style={{ marginTop: "40px" }}>
				{summaryTable}
			</div>
			<div style={{ marginTop: "40px"}}>
				<TabPanelMUI
					ContentTitles={AlleleSizePanelTitles}
					ContentData={AlleleSizePanelsData}
				>
				</TabPanelMUI>
			</div>
			<div style={{ marginTop: "40px" }}>
				{PhylogeneticElement ? PhylogeneticElement : PhylogeneticElement}
				{MSAComponent ? MSAComponent : MSAComponent}
			</div>
			<div style={{ marginTop: "40px"}}>
				{DNAEditor ? DNAEditor : DNAEditor}
				{ProteinEditor ? ProteinEditor : ProteinEditor}
			</div>
		</div>
	  );
};

export default LocusPage; 
