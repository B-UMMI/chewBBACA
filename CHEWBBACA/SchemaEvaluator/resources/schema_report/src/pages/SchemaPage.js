import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import TabPanelMUI from '../components/TabPanelMUI';


const SchemaPage = () => {

	// get pre-computed data
	const data = window.preComputedData;

	// data for Summary Data table
	const summaryData = data.summaryData;
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
		draggableColumns: {
			enabled: true
		}
	};
	// Component for Summary table
	const summaryTable = <DataTable
						  tableData={summaryData} 
						  tableTitle="Summary Data" 
						  tableOptions={summaryTableOptions}
						 >
						 </DataTable>


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
	const layoutPanelA = {bargroupgap: 0.05};
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
		 plotData={plotDataPanelA}
		 plotTitle="Distribution of the number of alleles"
		 xaxisTitle="Number of Alleles"
		 yaxisTitle="Number of Loci"
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
		 plotData={plotDataPanelB}
		 plotTitle="Distribution of the allele mode size"
		 xaxisTitle="Allele Mode Size"
		 yaxisTitle="Number of Loci"
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
	const layoutPanelC = {
		xaxis: {
			zeroline: true,
			showline: true
		},
		yaxis : {
			zeroline: true,
			showline: true
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
	// Component for Plotly Scatter with loci statistics
	const LociStatsScatter = (
		<PlotlyPlot
		 plotData={plotDataPanelC}
		 plotTitle="Locus Length Statistics"
		 xaxisTitle="Allele size (bp)"
		 yaxisTitle="Number of alleles"
		 layoutProps={layoutPanelC}
		 configOptions={configPanelC}
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
		x: Array(yDataPanelB.length).fill().map((element, index) => index),
		q1: q1PanelD,
		median: xDataMedianPanelC,
		q3: q3PanelD,
		lowerfence: xDataMinPanelC,
		upperfence: xDataMaxPanelC,
		showlegend: false,
		marker: {color: "#0570b0"}
	  }
	];
	const layoutPanelD = {
		xaxis: {
			rangeslider: {
				thickness: 0.15, 
			},
			showticklabels: false,
			range: [0, 100],
			title: {
				standoff: 0
			}
		}
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
	// Component for Plotly Boxplots with loci allele size variation
	const LociSizeBoxplots = (
		<PlotlyPlot
		 plotData={plotDataPanelD}
		 plotTitle="Locus Size Variation"
		 xaxisTitle="Loci"
		 yaxisTitle="Allele Size"
		 layoutProps={layoutPanelD}
		 configOptions={configPanelD}
		>
		</PlotlyPlot>
	);

	const TabPanelTitles = [
		"Total Alleles", 
		"Allele Mode Size", 
		"Locus Statistics", 
		"Allele Size Variation"
	];
	const TabPanelData = [
		TotalAllelesHistogram,
		AlleleModeHistogram,
		LociStatsScatter,
		LociSizeBoxplots
	];

	const locusColumnFormatting = { 
		0: {
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
		8: {
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
	const annotationsTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: "loci_annotations.tsv",
			separator: "\t",
			filterOptions: {
				useDisplayedColumnsOnly: true,
				useDisplayedRowsOnly: true
			}
		},
		filter: true,
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

	// Select conditional formatting function to apply
	// Only add href to loci names if the loci reports were generated
	let LociAnnotationsFormatting = undefined;
	if (data.lociReports === 1) {
		LociAnnotationsFormatting = {...locusColumnFormatting, ...uniprotLinkColumnFormatting};
	} else {
		LociAnnotationsFormatting = {...uniprotLinkColumnFormatting};
	}
	// DataTable component to display loci annotations
	const LociAnnotationsTable = (
		<DataTable 
		 tableData={annotationsData} 
		 tableTitle="Loci Annotations" 
		 tableOptions={annotationsTableOptions}
		 tableConditionalFormatting={LociAnnotationsFormatting}
		>
		</DataTable>
	);

	// data for CDS analysis
	const analysisData = data.analysis;
	const analysisTableOptions = {
		responsive: "simple",
		selectableRows: "none",
		selectableRowsHeader: false,
		selectableRowsHideCheckboxes: true,
		selectableRowsOnClick: false,
		print: false,
		download: true,
		downloadOptions: {
			filename: "cds_analysis.tsv",
			separator: "\t",
			filterOptions: {
				useDisplayedColumnsOnly: true,
				useDisplayedRowsOnly: true
			}
		},
		filter: true,
		search: true,
		viewColumns: true,
		pagination: true,
		rowsPerPage: 10,
		rowsPerPageOptions: [10, 20, 40, 80, 100],
		jumpToPage: true,
		draggableColumns: {
			enabled: true
		},
	};

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
		 tableTitle="CDS Analysis" 
		 tableOptions={analysisTableOptions}
		 tableConditionalFormatting={LociAnalysisFormatting}
		>
		</DataTable>
	);

	return (
		<div style={{ marginTop: "40px"}}>
			<div style={{ marginTop: "40px"}}>
				{summaryTable}
			</div>
			<div style={{ marginTop: "40px"}}>
				<TabPanelMUI
					ContentTitles={TabPanelTitles}
					ContentData={TabPanelData}
				>
				</TabPanelMUI>
			</div>
			<div style={{ marginTop: "40px"}}>
				{LociAnnotationsTable}
			</div>
			<div style={{ marginTop: "40px"}}>
				{LociAnalysisTable}
			</div>
		</div>
	  );
};

export default SchemaPage; 
