import { useState, useEffect } from "react";
import MSA from '../components/MSA';


const Resized = () => {
	const [width, setWidth] = useState(0);
	//const [height, setHeight] = useState(0);

	useEffect(() => {
		const resizeObserver = new ResizeObserver((event) => {
			setWidth(event[0].contentBoxSize[0].inlineSize);
			//setHeight(event[0].contentBoxSize[0].blockSize);
		});
		resizeObserver.observe(document.getElementById("MSA"));

		// using this removes the markers scrollbar in the MSA...
		//return resizeObserver.unobserve(document.getElementById("MSA"))
	});

	console.log(width)

	return (
		<div id="MSA">
			<MSA 
				MSAwidth={width}
			>
			</MSA>
		</div>
	)
};


export default Resized;