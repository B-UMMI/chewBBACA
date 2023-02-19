import { useEffect } from "react";

var msa = require("@jlab-contrib/msa");


const MSA = ({ MSAwidth }) => {

	useEffect(() => {
		// how to get this outside useEffect for first render?
		const msaData = window.preComputedDataInd.msa.sequences;
		// parsed array of the sequences
		var opts = {};
		// set your custom properties
		// @see: https://github.com/greenify/biojs-vis-msa/tree/master/src/g 
		//opts.seqs = msa.utils.seqgen.getDummySequences(1000,300);
		const sequences = msa.io.fasta.parse(msaData);
		opts.seqs = sequences;
		opts.el = document.getElementById("MSA");
		opts.vis = {
			conserv: false, 
			overviewbox: false,
			labels: true,
			labelName: true,
			labelId: false,
			leftHeader: false
			}

		const msaLength = (msaData.match(/>/g) || []).length;
		let MSAHeight = 450;
		if (msaLength < 30) {
			MSAHeight = msaLength*15;
		}

		opts.zoomer = {
			// need to subtract labelNameLength or the first time it renders it will
			alignmentWidth: MSAwidth-50,
			alignmentHeight: MSAHeight,
			columnWidth: 15,
			rowHeight: 15,
			autoResize: false,
			labelIdLength: 50,
			labelNameLength: 50,
			labelPartLength: 15,
			labelCheckLength: 15,
			labelFontsize: "13px",
			labelLineHeight: "13px",
			}
		opts.colorscheme = {
			scheme: "clustal"
		}
		opts.conf = {
			// awesome but several options are broken
			bootstrapMenu: false,
			hasRef: false
		}

		// init msa
		var m = new msa.msa(opts);

		// call render at the end to display the whole MSA
		m.render();

	}, [MSAwidth]);
};


export default MSA;
