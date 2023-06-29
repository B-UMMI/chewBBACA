import { useState, useEffect } from "react";

// Import array with tree format download options
import { treeDownloadOptions } from '../constants';

// Phylocanvas imports
import PhylocanvasGL, { plugins } from "@phylocanvas/phylocanvas.gl";
import { TreeTypes, Shapes } from "@phylocanvas/phylocanvas.gl";

// Package for animating vertical scrolling in the NJ Tree component
import { Element } from 'react-scroll';

import Box from '@mui/material/Box';
import Switch from '@mui/material/Switch';
import Select from '@mui/material/Select';
import Checkbox from '@mui/material/Checkbox';
import MenuItem from '@mui/material/MenuItem';
import TextField from '@mui/material/TextField';
import InputLabel from '@mui/material/InputLabel';
import IconButton from '@mui/material/IconButton';
import FormControl from '@mui/material/FormControl';
import Autocomplete from '@mui/material/Autocomplete';
import CheckBoxIcon from '@mui/icons-material/CheckBox';
import FileDownload from '@mui/icons-material/FileDownload';
import FormControlLabel from '@mui/material/FormControlLabel';
import CheckBoxOutlineBlankIcon from '@mui/icons-material/CheckBoxOutlineBlank';


const icon = <CheckBoxOutlineBlankIcon fontSize="small" />;
const checkedIcon = <CheckBoxIcon fontSize="small" />;


const PhylogeneticTree = ({ treeSource, validIDs }) => {
	const [phyloTree, setPhyloTree] = useState(undefined);

	// State used to control the type/shape of the NJ tree
	const [treeType, setTreeType] = useState(TreeTypes.Circular);
	// Event handler used to change NJ tree type
	const handleTreeTypeSelect = (event) => {
		setTreeType(event.target.value);
	};

	// State used to determine if NJ tree node labels should be displayed
	const [showNodeLabels, setShowNodeLabels] = useState(true);
	// Event handler used to change visibility of NJ tree node labels
	const handleTreeNodeLabelsCheck = (event) => {
		setShowNodeLabels(event.target.checked);
	};

	// State used to determine if NJ tree branch lenghts should be displayed
	const [showBranchLengths, setShowBranchLengths] = useState(false);
	// Event handler used to change visibility of NJ tree branch lengths
	const handleBranchLabelsCheck = (event) => {
		setShowBranchLengths(event.target.checked);
	};

	// State used to determine which NJ tree nodes should be highlighted
	const [selectedIds, setSelectedIds] = useState([]);
	// Event handler used to change NJ tree node selection
	const handleTreeNodeSelectChange = (event, newValue) => {
		setSelectedIds(newValue);
	};

	// State used to set NJ tree node size
	const [treeNodeSize, setTreeNodeSize] = useState(8);
	// Event handler used to change NJ tree node size
	const handleTreeNodeSize = (value) => {
		setTreeNodeSize(+value)
	};

	// State used to set NJ tree font size
	const [treeFontSize, setTreeFontSize] = useState(8);
	// Event handler used to change NJ tree font size
	const handleTreeFontSize = (value) => {
		setTreeFontSize(+value)
	};

	// State used to set NJ tree node shape
	const [treeNodeShape, setTreeNodeShape] = useState(Shapes.Circle);
	// Event handler used to change NJ tree node shape
	const handleTreeNodeShape = (event) => {
		setTreeNodeShape(event.target.value)
	}

	// Piece of state used to trigger tree download
	const [download, setDownload] = useState(false);
	// Event handler used to change download state and trigger download
	const handleClick = (event) => {
		setDownload(!download);
	}

	// Piece of state used to set download file format
	const [downloadType, setDownloadType] = useState("SVG");
	// Event handler used to change download file format
	const handleExport = (event) => {
		setDownloadType(event.target.value);
	};

	const treeDownloadMenuOptions = treeDownloadOptions.map((format) => {
		return <MenuItem key={format} value={format}>{format}</MenuItem>
	})

	const treeShapeOptions = Object.entries(TreeTypes).map((type) => {
		return <MenuItem key={type[0]} value={type[1]}>{type[0]}</MenuItem>
	})

	const treeNodeShapeOptions = Object.entries(Shapes).map((shape) => {
		return <MenuItem key={shape[0]} value={shape[1]}>{shape[0]}</MenuItem>
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
				limitTags={2}
				multiple
				disablePortal
				id="select-tree-id"
				options={validIDs}
				size="small"
				sx={{ height: 40, width: 300 }}
				renderOption={(props, option, { selected }) => (
					<li {...props}>
						<Checkbox
							icon={icon}
							checkedIcon={checkedIcon}
							style={{ marginRight: 8 }}
							checked={selected}
						/>
						{option}
					</li>
				)}
				renderInput={(params) => <TextField {...params} label="Select Ids" />}
				onChange={handleTreeNodeSelectChange}
				clearOnEscape
				disableCloseOnSelect
			/>
		</Box>
	);

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
					onKeyDown={(e) => {
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
					onChange={handleTreeNodeLabelsCheck}
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

	// Hook to check if tree is not undefined
	// and destroy previous tree before rerendering
	useEffect(() => {
		if (phyloTree) {
			phyloTree.destroy();
		}
	}, [treeSource, treeType, showNodeLabels, showBranchLengths, selectedIds, treeNodeSize, treeFontSize, treeNodeShape])

	useEffect(() => {
		if (phyloTree) {
			// need to set mime type based on file type
			let treeFile = undefined;
			if (downloadType === "SVG") {
				treeFile = new File([ phyloTree.exportSVG() ], "tree.svg", { type: "image/svg+xml" });
			} else if (downloadType === "Newick") {
				treeFile = new File([ phyloTree.exportNewick() ], "tree.nwk", { type: "text/plain" });
			} else if (downloadType === "JSON") {
				treeFile = new File([ phyloTree.exportJSON() ], "tree.json", { type: "text/plain" });
			}

			const anchor = document.createElement("a");

      		anchor.setAttribute("href", URL.createObjectURL(treeFile))
      		anchor.setAttribute("download", treeFile.name)
			anchor.click();
		}
	}, [download])

	useEffect(() => {
		const treeView = document.querySelector("#phyloTree");
		setPhyloTree(new PhylocanvasGL(
			treeView,
			// Default prop values at https://gitlab.com/cgps/phylocanvas/phylocanvas.gl/-/blob/main/defaults.js
			{type: treeType,
			 size: { width: 1440, height: 700 },
			 source: treeSource,
			 showLabels: true,
			 showLeafLabels: showNodeLabels,
			 interactive: true,
			 padding: 15,
			 treeToCanvasRatio: 0.25,
			 selectedIds: selectedIds,
			 showBranchLengths: showBranchLengths,
			 scaleLineAlpha: false,
			 nodeShape: treeNodeShape,
			 nodeSize: treeNodeSize,
			 fontSize: treeFontSize,
			 fontFamily: "monospace",
			 collapsedIds: [],
			 centre: [ 0.5, 0.5 ],
			// alignLabels: true
			},
			[plugins.scalebar]
		))
	}, [treeSource, treeType, showNodeLabels, showBranchLengths, selectedIds, treeNodeSize, treeFontSize, treeNodeShape])

	return (
		<div>
			{phyloTreeMenu}
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
				<div id="phyloTree" style={{ margin: "auto" }}></div>
			</Element>
		</div>
	)
};

export default PhylogeneticTree;
