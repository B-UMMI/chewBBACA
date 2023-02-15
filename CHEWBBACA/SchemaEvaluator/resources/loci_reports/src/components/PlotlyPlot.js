import Plot from "react-plotly.js";


const PlotlyPlot = ({ plotData, plotTitle, xaxisTitle, yaxisTitle, layoutProps, configOptions }) => {

	return (
		<Plot
        data={plotData}
        layout={{
          title: {
            text: plotTitle,
          },
          xaxis: {
            title: { text: xaxisTitle },
            showgrid: true,
          },
          yaxis: {
            title: { text: yaxisTitle },
          },
          ...layoutProps
        }}
        config={configOptions}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        line={{
          width: 1,
        }}
      />
	)
};

export default PlotlyPlot;
