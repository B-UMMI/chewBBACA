import { useState, useEffect } from "react";

// Phylocanvas imports
import PhylocanvasGL, { plugins } from "@phylocanvas/phylocanvas.gl";


const PhylogeneticTree = ({ treeSource, treeType, showLabels, showBranchLengths, selectedIds, download, downloadType, treeNodeSize, fontSize, nodeShape }) => {
	const [phyloTree, setPhyloTree] = useState(undefined);

	// Hook to check if tree is not undefined
	// and destroy previous tree before rerendering
	useEffect(() => {
		if (phyloTree) {
			phyloTree.destroy();
		}
	}, [treeSource, treeType, showLabels, showBranchLengths, selectedIds, treeNodeSize, fontSize, nodeShape])

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
			 showLeafLabels: showLabels,
			 interactive: true,
			 padding: 15,
			 treeToCanvasRatio: 0.25,
			 selectedIds: selectedIds,
			 showBranchLengths: showBranchLengths,
			 scaleLineAlpha: false,
			 nodeShape: nodeShape,
			 nodeSize: treeNodeSize,
			 fontSize: fontSize,
			 fontFamily: "monospace",
			 collapsedIds: [],
			 centre: [ 0.5, 0.5 ],
			// alignLabels: true
			},
			[plugins.scalebar]
		))
	}, [treeSource, treeType, showLabels, showBranchLengths, selectedIds, treeNodeSize, fontSize, nodeShape])

	return (
		<div id="phyloTree" style={{ margin: "auto" }}></div>
	)
};

export default PhylogeneticTree;
