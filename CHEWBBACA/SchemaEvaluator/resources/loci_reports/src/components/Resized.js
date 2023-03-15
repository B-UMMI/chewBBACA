import { useState, useEffect } from "react";


const Resized = ({ divID, component, colorScheme }) => {
	const [width, setWidth] = useState(0);
	//const [height, setHeight] = useState(0);

	useEffect(() => {
		const resizeObserver = new ResizeObserver((event) => {
			setWidth(event[0].contentBoxSize[0].inlineSize);
			//setHeight(event[0].contentBoxSize[0].blockSize);
		});

		resizeObserver.observe(document.getElementById(divID));

		// using this removes the markers scrollbar in the MSA...
		//return resizeObserver.unobserve(document.getElementById("MSA"))
	});

	//console.log(width, divID)

	// get constructor of the component that will be displayed and automatically resized
	const ComponentDisplay = component;

	return (
		<div id={divID}>
			<ComponentDisplay MSAwidth={width} colorScheme={colorScheme}></ComponentDisplay>
		</div>
	)
};


export default Resized;
