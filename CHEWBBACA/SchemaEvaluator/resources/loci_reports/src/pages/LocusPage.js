import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import Resized from '../components/Resized';
import AccordionMUI from '../components/AccordionMUI';
import TabPanelMUI from '../components/TabPanelMUI';
import MSA from '../components/MSA';

// Material-UI components
import Box from '@mui/material/Box';
import Alert from '@mui/material/Alert';

// Phylocanvas
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
			filename: "locus_summary.tsv",
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
			filename: "locus_annotations.tsv",
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
			<Alert severity="info" variant="outlined">
				<Typography sx={{ fontSize: 14 }}>
					Loci annotations were not provided.
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
		>
		</PlotlyPlot>
	);

	const AlleleSizePanelTitles = ["Allele Size Counts", "Allele Size"];
	const AlleleSizePanelsData = [AlleleSizeBar, AlleleSizeScatter];

	// get data for Phylocanvas tree
	const phyloData = data.phylo.phylo_data;
	let PhylogeneticElement = undefined;
	if (phyloData.length > 0) {
		const phylogeneticElementTree = (
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

		const phylogeneticElementTitle = (
			<Typography sx={{ 
				color: '#bb7944', 
				fontSize: 20 
				}}
			>
				Phylogenetic Tree
			</Typography>

		);

		PhylogeneticElement = (
			<AccordionMUI
				summaryText={phylogeneticElementTitle}
				detailsData={phylogeneticElementTree}
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

		const MSAComponentTitle = (
			<Typography sx={{ 
				color: '#bb7944', 
				fontSize: 20 
				}}
			>
				Multiple Sequence Alignment
			</Typography>

		);

		MSAComponent = (
			<AccordionMUI
				summaryText={MSAComponentTitle}
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

		const DNAEditorTitle = (
			<Typography sx={{ 
				color: '#bb7944', 
				fontSize: 20 
				}}
			>
				DNA sequences
			</Typography>

		);

		DNAEditor = (
			<AccordionMUI
			 summaryText={DNAEditorTitle}
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

		const ProteinEditorTitle = (
			<Typography sx={{ 
				color: '#bb7944', 
				fontSize: 20 
				}}
			>
				Protein sequences
			</Typography>

		);

		ProteinEditor = (
			<AccordionMUI
			 summaryText={ProteinEditorTitle}
			 detailsData={ProteinEditor}
			 expanded={false}
			>
			</AccordionMUI>
		);
	}

	// Alert to explain red lines added for bot and top length thresholds
	const lengthThresholdsAlert = (
		<Alert severity="info" variant="outlined">
			<Typography sx={{ fontSize: 14 }}>
				The red lines represent the bottom and top allele length thresholds.
			</Typography>
		</Alert>
	)

	// Alert rendered when there is no Phylogenetic Tree and MSA
	const alertPhyloMSA = (
		<Alert variant="outlined" severity="warning">
        	<Typography>
              The NJ tree and MSA were not generated because this locus 
			  does not have valid alleles or has only one valid allele or
			  the --light parameter was provided.
            </Typography>
    	</Alert>
	);

	// Alert to explain leaf node labels
	const alertLeafNodes = (
		<Alert variant="outlined" severity="info">
            <Typography>
                Leaf labels correspond to the allele identifiers.
				A tree might not be displayed if the distance between all alleles is 0.
            </Typography>
        </Alert>
	);

	return (
		<div style={{ marginTop: "10px", marginBottom: "10px" }}>
			<div style={{ marginTop: "10px" }}>
				{summaryTable}
			</div>
			<div style={{ marginTop: "30px" }}>
				<div style={{ marginBottom: "10px" }}>
					{LociAnnotationsAlert}
				</div>
				{annotationsTable}
			</div>
			<div style={{ marginTop: "30px"}}>
				<div style={{ marginBottom: "10px" }}>
					{lengthThresholdsAlert}
				</div>
				<TabPanelMUI
					ContentTitles={AlleleSizePanelTitles}
					ContentData={AlleleSizePanelsData}
				>
				</TabPanelMUI>
			</div>
			<div style={{ marginTop: "30px" }}>
				<div style={{ marginBottom: "10px" }}>
					{PhylogeneticElement ? alertLeafNodes : undefined}
					{PhylogeneticElement ? undefined : alertPhyloMSA}
				</div>
				{PhylogeneticElement}
				{MSAComponent}
			</div>
			<div style={{ marginTop: "30px"}}>
				{DNAEditor}
				{ProteinEditor}
			</div>
		</div>
	  );
};

export default LocusPage; 
