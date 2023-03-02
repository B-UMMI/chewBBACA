import Plot from "react-plotly.js";


const PlotlyPlot = ({ plotData, layoutProps, configOptions }) => {

	return (
		<Plot
        data={plotData}
        layout={layoutProps}
        config={configOptions}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
      />
	)
};

export default PlotlyPlot;
