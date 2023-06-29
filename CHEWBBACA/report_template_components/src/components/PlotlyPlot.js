import Plot from "react-plotly.js";


// use "...rest" notation to enable passing any number of custom event handlers
const PlotlyPlot = ({ plotData, layoutProps, configOptions, ...rest }) => {
	return (
		<Plot
            data={plotData}
            layout={layoutProps}
            config={configOptions}
            useResizeHandler={true}
            style={{ width: "100%", height: "100%" }}
            {...rest}
        />
	)
};

export default PlotlyPlot;
