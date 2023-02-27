import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import Resized from '../components/Resized';
import AccordionMUI from '../components/AccordionMUI';
import TabsMUI from '../components/TabPanelMUI';
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
	const xDataPanelA = data.lengths;
	const yDataPanelA = data.ids;
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
	// Component for Plotly Histogram with allele size distribution
	const AlleleSizeHistogram = (
		<PlotlyPlot
		 plotData={plotDataPanelA}
		 plotTitle={locusName}
		 xaxisTitle="Sequence Size (bp)"
		 yaxisTitle="Number of Alleles"
		 layoutProps={layoutPanelA}
		 configOptions={configPanelA}
		>
		</PlotlyPlot>
	);

	// data for Panel B (Allele Mode Size, )
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
	const AlleleSizeScatter = (
		<PlotlyPlot 
		 plotData={plotDataPanelB}
		 plotTitle={locusName}
		 xaxisTitle="Allele ID"
		 yaxisTitle="Sequence Size (bp)"
		 layoutProps={layoutPanelB}
		 configOptions={configPanelB}
		>
		</PlotlyPlot>
	);

	const AlleleSizePanelTitles = ["Allele Size Distribution", "Allele Size"];
	const AlleleSizePanelsData = [AlleleSizeHistogram, AlleleSizeScatter];

	// get data for Phylocanvas tree
	const phyloData = data.phylo.phylo_data;

	const PhylogeneticElement = (
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

	// create component for MSA
	// pass MSA component constructor instead of instance
	// this allows to get constructor and pass props in component that receives constructor
	const MSAComponent = (
		<Box sx={{ p: 3 }}>
			<Resized
				divID="MSA"
				component={MSA}
			>
			</Resized>
		</Box>
	);

	// get DNA sequences
	const dnaSequences = data.dna.sequences
	const dnaText = dnaSequences.map((seq) => {
		const seqid = seq.name;
		const sequence = seq.sequence;
		const sequenceStr = `>${seqid}\n${sequence}\n`
		return sequenceStr
	})

	const joinedDNA = dnaText.join('');
	// create DNA Editor component
	const DNAEditor = (
		<Editor
		 height="40vh"
		 options={{"readOnly": true, "wordWrap": "on"}}
		 defaultValue={`${joinedDNA}`}
		>
		</Editor>
	);

	// get Protein sequences
	const proteinSequences = data.protein.sequences
	const proteinText = proteinSequences.map((seq) => {
		const seqid = seq.name;
		const sequence = seq.sequence;
		const sequenceStr = `>${seqid}\n${sequence}\n`
		return sequenceStr
	})

	const joinedProtein = proteinText.join('');
	// create DNA Editor component
	const ProteinEditor = (
		<Editor
		 height="40vh"
		 options={{"readOnly": true, "wordWrap": "on"}}
		 defaultValue={`${joinedProtein}`}
		>
		</Editor>
	);

	return (
		<div style={{ marginTop: "40px" }}>
			<div style={{ marginTop: "40px" }}>
				{summaryTable}
			</div>
			<div style={{ marginTop: "40px"}}>
				<TabsMUI
					ContentTitles={AlleleSizePanelTitles}
					ContentData={AlleleSizePanelsData}
				>
				</TabsMUI>
			</div>
			<div style={{ marginTop: "40px" }}>
				<AccordionMUI
					summaryText="Phylogenetic Tree"
					detailsData={PhylogeneticElement}
					expanded={false}
				>
				</AccordionMUI>
				<AccordionMUI
					summaryText="Multiple Sequence Alignment"
					detailsData={MSAComponent}
					expanded={false}
				>
				</AccordionMUI>
			</div>
			<div style={{ marginTop: "40px"}}>
				<AccordionMUI
					summaryText="DNA sequences"
					detailsData={DNAEditor}
					expanded={false}
				>
				</AccordionMUI>
				<AccordionMUI
					summaryText="Protein sequences"
					detailsData={ProteinEditor}
					expanded={false}
				>
				</AccordionMUI>
			</div>
		</div>
	  );
};

export default LocusPage; 
