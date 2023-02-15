import { useEffect } from "react";

// Phylocanvas imports
import PhylocanvasGL, { TreeTypes, plugins } from "@phylocanvas/phylocanvas.gl";


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
			 padding: 15,
			// alignLabels: true
			},
			[plugins.scalebar,
			 ]
		)
	}, [])
};

export default PhylogeneticTree;
