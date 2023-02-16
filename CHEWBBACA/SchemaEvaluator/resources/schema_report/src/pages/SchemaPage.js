import { useState } from 'react';
import DataTable from '../components/DataTable';
import PlotlyPlot from '../components/PlotlyPlot';
import Box from '@mui/material/Box';
import Paper from '@mui/material/Paper';
import Tabs from '@mui/material/Tabs';
import Tab from '@mui/material/Tab';


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


const SchemaPage = () => {
	const [panel, setPanel] = useState(0);

	const handleChange = (event, newValue) => {
		// where is this value coming from?
		setPanel(newValue);
	};

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
	};

	// data for Panel A
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

	// data for Panel B
	const xDataPanelB = data.mode;
	const yDataPanelB = data.loci;
	const plotDataPanelB = [
		{x: xDataPanelB,
		 y: yDataPanelB,
		 type: "histogram",
		 name: "Distribution of allele mode sizes per gene",
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

	// data for Panel C
	const xDataMinPanelC = data.min
	const xDataMaxPanelC = data.max
	const xDataMedianPanelC = data.median

	const plotDataPanelC = [
		{
		  x: xDataMinPanelC,
		  y: xDataPanelA,
		  type: "scatter",
		  name: "Min",
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
		  name: "Max",
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
			zeroline: false,
			showline: true
		},
		yaxis : {
			zeroline: false,
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
	};

	// data for CDS analysis
	const analysisData = data.analysis;
	const analysisTableOptions = {
		responsive: "simple",
		selectableRowsHeader: false,
		selectableRows: "none",
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
	};

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
								<Tab label="Locus Statistics" wrapped {...a11yProps(2)} />
								<Tab label="Allele Size Variation" wrapped {...a11yProps(3)} />
							</Tabs>
						</Box>
						<TabPanel value={panel} index={0}>
							<PlotlyPlot 
								plotData={plotDataPanelA}
								plotTitle="Number of Loci with given Number of Alleles"
								xaxisTitle="Number of Different Alleles"
								yaxisTitle="Number of Loci"
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
						<TabPanel value={panel} index={2}>
							<PlotlyPlot 
								plotData={plotDataPanelC}
								plotTitle="Locus Statistics"
								xaxisTitle="Allele size (bp)"
								yaxisTitle="Number of alleles"
								layoutProps={layoutPanelC}
								configOptions={configPanelC}
							>
							</PlotlyPlot>
						</TabPanel>
						<TabPanel value={panel} index={3}>
							<PlotlyPlot 
								plotData={plotDataPanelD}
								plotTitle="Locus Size Variation"
								xaxisTitle="Loci"
								yaxisTitle="Allele size variation"
								layoutProps={layoutPanelD}
								configOptions={configPanelD}
							>
							</PlotlyPlot>
						</TabPanel>
					</Paper>
				</Box>
			</div>
			<div style={{ marginTop: "40px"}}>
				<DataTable 
					tableData={annotationsData} 
					tableTitle="Loci Annotations" 
					tableOptions={annotationsTableOptions}
					tableConditionalFormatting={{...locusColumnFormatting, ...uniprotLinkColumnFormatting}}
				>
				</DataTable>
			</div>
			<div style={{ marginTop: "40px"}}>
				<DataTable 
					tableData={analysisData} 
					tableTitle="CDS Analysis" 
					tableOptions={analysisTableOptions}
					tableConditionalFormatting={{...locusColumnFormatting}}
				>
				</DataTable>
			</div>
		</div>
	  );
};

export default SchemaPage; 
