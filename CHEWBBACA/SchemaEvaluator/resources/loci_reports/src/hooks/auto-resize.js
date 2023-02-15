import { useEffect, useState } from "react";


const PageSize = () => {

	const [windowWidth, setWindowWidth] = useState(window.innerWidth);
	const [windowHeight, setWindowHeight] = useState(window.innerHeight);

	const setWindowDimensions = () => {
		setWindowWidth(window.innerWidth)
		setWindowHeight(window.innerHeight)
	};

	useEffect(() => {
		window.addEventListener('resize', setWindowDimensions);
		return () => {
			window.removeEventListener('resize', setWindowDimensions)
		}
	}, [])

};

export default PageSize;