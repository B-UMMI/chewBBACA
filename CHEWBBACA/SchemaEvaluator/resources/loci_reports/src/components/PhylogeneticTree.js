import { useEffect, useState } from "react";

// Phylocanvas imports
import PhylocanvasGL, { TreeTypes } from "@phylocanvas/phylocanvas.gl";


const PhylogeneticTree = ({ source, treeWidth, treeHeight, showLabels, showLeafLabels, interactive }) => {
	console.log('Phylo');

	useEffect(() => {
		const tree = new PhylocanvasGL(
			document.querySelector("#demo"),
			{type: TreeTypes.Square,
			 size: { width: treeWidth, height: treeHeight },
			 source: source,
			 showLabels,
			 showLeafLabels,
			 interactive,
			 padding: 10
			}
		)
		// can I use this to destroy the tree before each rerender?
		// return tree.destroy();
	}, [])
};

export default PhylogeneticTree;
