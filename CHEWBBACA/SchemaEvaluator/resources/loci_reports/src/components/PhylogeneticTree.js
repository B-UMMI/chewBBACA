import { useState, useEffect } from "react";

// Phylocanvas imports
import PhylocanvasGL, { plugins } from "@phylocanvas/phylocanvas.gl";


const PhylogeneticTree = ({ treeSource, treeType, showLabels }) => {
	const [phyloTree, setPhyloTree] = useState(undefined);

	// Hook to check if tree is not undefined
	// and destroy previous tree before rerendering
	useEffect(() => {
		if (phyloTree) {
			phyloTree.destroy();
		}
	}, [treeSource, treeType, showLabels])

	useEffect(() => {
		const treeView = document.querySelector("#phyloTree");
		setPhyloTree(new PhylocanvasGL(
			treeView,
			{type: treeType,
			 size: { width: 600, height: 700 },
			 source: treeSource,
			 showLabels: showLabels,
			 showLeafLabels: true,
			 interactive: true,
			 padding: 15,
			// alignLabels: true
			},
			[plugins.scalebar]
		))
	}, [treeSource, treeType, showLabels])

	return (
		<div id="phyloTree" style={{ margin: "auto" }}></div>
	)
};

export default PhylogeneticTree;
