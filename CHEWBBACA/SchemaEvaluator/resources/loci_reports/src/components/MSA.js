import { useState, useEffect } from "react";

// Import arrays with MSA color schemes and download format options
import { msaColorSchemes, msaDownloadOptions } from '../constants';

import Box from '@mui/material/Box';
import Switch from '@mui/material/Switch';
import Select from '@mui/material/Select';
import MenuItem from '@mui/material/MenuItem';
import TextField from '@mui/material/TextField';
import InputLabel from '@mui/material/InputLabel';
import IconButton from '@mui/material/IconButton';
import FormControl from '@mui/material/FormControl';
import FileDownload from '@mui/icons-material/FileDownload';
import FormControlLabel from '@mui/material/FormControlLabel';


var msa = require("@jlab-contrib/msa");


// For more info about features https://www.npmjs.com/package/@jlab-contrib/msa
// Code examples to filter sequences https://github.com/jupyterlab-contrib/msa/blob/master/src/menu/views/FilterMenu.js
const MSA = ({ MSADAta }) => {
	const [msaView, setMSAView] = useState(undefined);

	const [msaRenderCount, setRenderCount] = useState(0);
	const handleRenderCount = (event) => {
		setRenderCount(msaRenderCount+1);
	};

	// State used to set color scheme
	const [colorScheme, setColorScheme] = useState("lesk");
	// Event handler used to change color scheme
	const handleColorChange = (event) => {
		setColorScheme(event.target.value)
	};

	// State used to trigger MSA download
	const [msaExport, setMSAExport] = useState(false);
	// Event handler used to trigger MSA download
	const handleMSAExport = (event) => {
		setMSAExport(!msaExport)
	};

	// State used to set MSA export format
	const [downloadMSAType, setDownloadMSAType] = useState("Full MSA");
	// Event handler used to change MSA export format
	const handleMSATypeExport = (event) => {
		setDownloadMSAType(event.target.value);
	};

	// State used to determine if MSA conservation values should be displayed
	const [showConservation, setShowConservation] = useState(false);
	// Event handler used to change visibility of MSA conservation values
	const handleMSAConservation = (event) => {
		setShowConservation(event.target.checked);
	};

	// State used to determine if MSA SeqLogo should be displayed
	const [showSeqLogo, setShowSeqLogo] = useState(false);
	// Event handler used to change visibility of MSA SeqLogo
	const handleSeqLogo = (event) => {
		setShowSeqLogo(event.target.checked);
	};

	// MSA component does not handle well simultaneous row and column selection
	// const [msaRow, setMSARow] = useState(undefined);
	// const handleMSARow = (value) => {
	// 	setMSARow(value);
	// };

	// State used to set first MSA column that is visible
	const [msaColumn, setMSAColumn] = useState(undefined);
	// Event handler used to change the first MSA column that is visible
	const handleMSAColumn = (value) => {
		// Need to subtract to get correct column index
		setMSAColumn(value-1)
	};

	// State used to set the motif to search for in the MSA
	const [searchMotif, setSearchMotif] = useState(undefined);
	// Event handler used to change the motif to search for in the MSA
	const handleSearchMotif = (value) => {
		setSearchMotif(value)
	};

	const colorSchemesOptions = msaColorSchemes.map((color) => {
		return <MenuItem key={color[0]} value={color[0]}>{color[1]}</MenuItem>
	})

	const msaDownloadMenuOptions = msaDownloadOptions.map((format) => {
		return <MenuItem key={format} value={format}>{format}</MenuItem>
	})

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

	useEffect(() => {
		// parsed array of the sequences
		let opts = {};
		// set your custom properties
		// @see: https://github.com/greenify/biojs-vis-msa/tree/master/src/g 
		const sequences = msa.io.fasta.parse(MSADAta);
		opts.seqs = sequences;
		opts.el = document.getElementById("MSA");
		opts.vis = {
			conserv: showConservation, 
			overviewbox: false,
			labels: true,
			labelName: true,
			labelId: false,
			leftHeader: false,
			seqlogo: showSeqLogo
		}

		const msaLength = (MSADAta.match(/>/g) || []).length;
		let MSAHeight = 450;
		if (msaLength < 30) {
			MSAHeight = msaLength*15;
		}

		opts.zoomer = {
			// need to subtract labelNameLength or the first time it renders it will
			//alignmentWidth: MSAwidth-50,
			alignmentHeight: MSAHeight,
			columnWidth: 15,
			rowHeight: 15,
			autoResize: true,
			labelIdLength: 50,
			labelNameLength: 50,
			labelPartLength: 15,
			labelCheckLength: 15,
			labelFontsize: "13px",
			labelLineHeight: "13px",
			markerFontsize: "10px",
		}
		opts.colorscheme = {
			scheme: colorScheme
		}
		opts.conf = {
			// awesome but several options are broken
			bootstrapMenu: false,
			hasRef: false,
			alphabetSize: 20
		}

		// init msa
		setMSAView(new msa.msa(opts));
	}, []);

	useEffect(() => {
		if (msaView) {
			msaView.g.colorscheme.set("scheme", colorScheme);
			msaView.g.vis.set("conserv", showConservation);
			msaView.g.vis.set("seqlogo", showSeqLogo);
		}
	}, [colorScheme, showConservation, showSeqLogo])

	if (msaView) {
		msaView.render();
		// Redefine width in first render
		if (msaRenderCount === 0) {
			// get current alignment width
			const currentWidth = msaView.g.zoomer.attributes.alignmentWidth
			// Subtract 12% of current value to avoid overflow
			msaView.g.zoomer.set("alignmentWidth", currentWidth*(1-0.12))
			handleRenderCount()
		}
	}

	useEffect(() => {
		if (msaView) {
			if (downloadMSAType === "Full MSA") {
				msa.utils.exporter.saveAsFile(msaView, "msa.fasta")
			} else if (downloadMSAType === "Row Selection") {
				msa.utils.exporter.saveSelection(msaView, "msaSelected.fasta")
			} else if (downloadMSAType === "Image (PNG)") {
				msa.utils.exporter.saveAsImg(msaView, "msa.png")
			}
		}
	}, [msaExport])

	useEffect(() => {
		if (msaView) {
			if (msaColumn !== -1) {
				msaView.g.zoomer.setLeftOffset(msaColumn);
			}
		}
		// MSA component does not handle well simultaneous row and column selection
		// 	if (msaRow) {
		// 		// get sequence indices
		// 		const records = (msaView.seqs._byId)
		// 		let seqids = []
		// 		for (const key in records) {
		// 			seqids.push(records[key].attributes.name)
		// 		}
		// 		const seqidIndex = seqids.indexOf(msaRow);
		// 		console.log(msaRow, seqidIndex)
		// 		msaView.g.zoomer.setLeftOffset(msaColumn);
		// 		if (seqidIndex) {
		// 			msaView.g.zoomer.setTopOffset(seqidIndex);
		// 			console.log("row")
		// 		} else {
		// 			console.log("No such seqid.")
		// 		}
		// 	}
		// }
	}, [msaColumn])

	useEffect(() => {
		if (msaView) {
			if (searchMotif) {
				msaView.g.user.set("searchText", searchMotif);
				msaView.g.selcol.reset();
			}
		}
	}, [searchMotif])

	return (
		<div>
			{msaMenu}
			<div id="MSA"></div>
		</div>
	)
};


export default MSA;
