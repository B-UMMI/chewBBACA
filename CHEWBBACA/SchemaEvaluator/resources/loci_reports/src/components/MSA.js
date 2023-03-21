import { useState, useEffect } from "react";

var msa = require("@jlab-contrib/msa");


// For more info about features https://www.npmjs.com/package/@jlab-contrib/msa
// Code examples to filter sequences https://github.com/jupyterlab-contrib/msa/blob/master/src/menu/views/FilterMenu.js
const MSA = ({ MSADAta, colorScheme, msaExport, downloadMSAType, conservation, seqLogo, msaColumn, searchMotif }) => {
	const [msaView, setMSAView] = useState(undefined);

	const [msaRenderCount, setRenderCount] = useState(0);
	const handleRenderCount = (event) => {
		setRenderCount(msaRenderCount+1);
	};

	useEffect(() => {
		// parsed array of the sequences
		let opts = {};
		// set your custom properties
		// @see: https://github.com/greenify/biojs-vis-msa/tree/master/src/g 
		const sequences = msa.io.fasta.parse(MSADAta);
		opts.seqs = sequences;
		opts.el = document.getElementById("MSA");
		opts.vis = {
			conserv: conservation, 
			overviewbox: false,
			labels: true,
			labelName: true,
			labelId: false,
			leftHeader: false,
			seqlogo: seqLogo
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
			msaView.g.vis.set("conserv", conservation);
			msaView.g.vis.set("seqlogo", seqLogo);
			console.log(msaView.g.selcol.where({type: "row"}))
		}
	}, [colorScheme, conservation, seqLogo])

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
		<div id="MSA"></div>
	)
};


export default MSA;
