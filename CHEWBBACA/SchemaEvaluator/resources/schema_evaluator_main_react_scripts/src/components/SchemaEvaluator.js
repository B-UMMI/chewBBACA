import React, { Component } from "react";

import Aux from "../hoc/Aux";
import classes from "./SchemaEvaluator.css";

// Material-UI components
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
// import FormControlLabel from "@material-ui/core/FormControlLabel";
import { createMuiTheme, MuiThemeProvider } from "@material-ui/core/styles";

// Material-UI ExpansionPanel components
import Accordion from "@material-ui/core/Accordion";
import ExpandMoreIcon from "@material-ui/icons/ExpandMore";
import AccordionSummary from "@material-ui/core/AccordionSummary";
import AccordionDetails from "@material-ui/core/AccordionDetails";

// react-select import
// import Select from "react-select";

// Material-UI Datatables
import MUIDataTable from "mui-datatables";

// Plotly.js
import Plot from "react-plotly.js";


class SchemaEvaluator extends Component {
  state = {
    pre_computed_data: _preComputedData,
    locus_ind_data: _preComputedDataInd,
    pre_computed_data_boxplot: _preComputedDataBoxplot,
    cds_df_data: _cdsDf,
    cds_scatter_data: _cdsScatter,
    tabValue: 0,
    indTabValue: 0,
    scatterSelectOption: "",
    isIndOption: false,
  };

  getMuiTheme = () =>
    createMuiTheme({
      overrides: {
        MUIDataTableToolbar: {
          titleText: {
            color: "#bb7944",
          },
        },
      },
    });

  plotChangeHandler = (value) => {
    this.setState({ tabValue: value });
  };

  clickScatterPlotHandler = (event) => {
    const locus_id = event.points[0].text;
    console.log(locus_id);

    const anchor = document.createElement("a");
    anchor.href = `../html_files/${locus_id.split(".")[0]}_individual_report.html`;
    anchor.target = "_blank";
    anchor.rel = "noopener noreferrer";
    anchor.click();

  };

  clickBoxPlotHandler = (event) => {
    const locus_id = event.points[0].x;
    console.log(locus_id);

    this.setState({
      scatterSelectOption: locus_id,
      isIndOption: true,
    });
  };

  render() {
    const style = {
      buttonBar: {
        overflowX: "auto",
        display: "flex",
        justifyContent: "center",
        marginBottom: "20px",
      },
      button: {
        minWidth: "150px",
      },
      formControl: {
        margin: 0,
        display: "flex",
        wrap: "nowrap",
      },
      root: {
        width: "100%",
      },
      toolbar: {
        display: "flex",
      },
    };

    // Create panel A
    let locus_name2 = [];
    let nr_alleles = [];
    let total_al_data = [];
    for (let key in this.state.pre_computed_data.total_alleles) {
      locus_name2.push(
        this.state.pre_computed_data.total_alleles[key].locus_name
      );
      nr_alleles.push(
        this.state.pre_computed_data.total_alleles[key].nr_alleles
      );
    }

    total_al_data.push({
      x: nr_alleles,
      y: locus_name2,
      type: "histogram",
      name: "Total alleles",
    });

    let total_allele_plot = (
      <Plot
        data={total_al_data}
        layout={{
          title: {
            text: "Number of Loci with given Number of Alleles",
          },
          xaxis: {
            title: { text: "Number of Different Alleles" },
            showgrid: true,
          },
          yaxis: {
            title: { text: "Number of Loci" },
          },
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        line={{
          width: 1,
        }}
      />
    );

    // Create panel B
    let mode_data_arr = [];
    let allele_mode = [];
    let locus_name = [];

    for (let key in this.state.pre_computed_data.mode) {
      allele_mode.push(this.state.pre_computed_data.mode[key].alleles_mode);
      locus_name.push(this.state.pre_computed_data.mode[key].locus_name);
    }

    mode_data_arr.push({
      x: allele_mode,
      y: locus_name,
      type: "histogram",
      name: "Distribution of allele mode sizes per gene",
    });

    let mode_data_plot = (
      <Plot
        data={mode_data_arr}
        layout={{
          title: {
            text: "Distribution of allele mode sizes",
          },
          xaxis: {
            title: { text: "Allele Mode Size" },
            showgrid: true,
          },
          yaxis: {
            title: { text: "Number of Loci" },
          },
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        line={{
          width: 1,
        }}
      />
    );

    // Create panel C
    let locus_id = [];
    let nr_alleles_scatter = [];
    let scatter_data = [];
    let scatter_data_median = [];
    let scatter_data_min = [];
    let scatter_data_max = [];

    for (let key in this.state.pre_computed_data.scatter_data) {
      locus_id.push(this.state.pre_computed_data.scatter_data[key].locus_id);
      nr_alleles_scatter.push(
        this.state.pre_computed_data.scatter_data[key].nr_alleles
      );
      scatter_data_median.push(
        this.state.pre_computed_data.scatter_data[key].alleles_median
      );
      scatter_data_min.push(
        this.state.pre_computed_data.scatter_data[key].alleles_min
      );
      scatter_data_max.push(
        this.state.pre_computed_data.scatter_data[key].alleles_max
      );
    }

    scatter_data.push(
      {
        x: scatter_data_min,
        y: nr_alleles_scatter,
        type: "scatter",
        name: "Min",
        mode: "markers",
        marker: {
          opacity: 0.7,
          size: 4,
        },
        hovertemplate: "<b>ID</b>: %{text}",
        text: locus_id,
      },
      {
        x: scatter_data_max,
        y: nr_alleles_scatter,
        type: "scatter",
        name: "Max",
        mode: "markers",
        marker: {
          opacity: 0.7,
          size: 4,
        },
        hovertemplate: "<b>ID</b>: %{text}",
        text: locus_id,
      },
      {
        x: scatter_data_median,
        y: nr_alleles_scatter,
        type: "scatter",
        name: "Median",
        mode: "markers",
        marker: {
          opacity: 0.7,
          size: 4,
        },
        hovertemplate: "<b>ID</b>: %{text}",
        text: locus_id,
      }
    );

    let scatter_plot = (
      <Plot
        data={scatter_data}
        layout={{
          title: {
            text: "Locus Statistics",
          },
          xaxis: {
            title: { text: "Allele size in bp" },
            showgrid: true,
            zeroline: false,
          },
          yaxis: {
            title: { text: "Number of alleles" },
            zeroline: false,
          },
          hovermode: "closest",
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        line={{
          width: 1,
        }}
        onClick={(e) => this.clickScatterPlotHandler(e)}
      />
    );

    // Panel D
    let boxplotData = [];

    boxplotData.push({
      type: "box",
      name: "Locus Size Variation",
      x: this.state.pre_computed_data_boxplot.loci,
      q1: this.state.pre_computed_data_boxplot.q1,
      median: this.state.pre_computed_data_boxplot.median,
      q3: this.state.pre_computed_data_boxplot.q3,
      lowerfence: this.state.pre_computed_data_boxplot.min,
      upperfence: this.state.pre_computed_data_boxplot.max,
      showlegend: false,
    });

    let boxplot = (
      <Plot
        data={boxplotData}
        layout={{
          title: {
            text: "Locus Size Variation",
          },
          xaxis: {
            title: { text: "Loci" },
            showticklabels: false,
          },
          yaxis: {
            title: { text: "Allele size variation" },
          },
          boxgap: 0.05,
          boxgapgroup: 0.05,
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        onClick={(e) => this.clickBoxPlotHandler(e)}
      />
    );

    // CDS Analysis
    const columns = [
      {
        name: "Gene",
        label: "Locus",
        options: {
          filter: true,
          sort: true,
          display: true,
          setCellHeaderProps: (value) => {
            return {
              style: {
                fontWeight: "bold",
              },
            };
          },
        },
      },
      {
        name: "Alleles not multiple of 3",
        label: "Alleles not multiple of 3",
        options: {
          filter: true,
          sort: true,
          display: true,
          setCellHeaderProps: (value) => {
            return {
              style: {
                fontWeight: "bold",
              },
            };
          },
        },
      },
      {
        name: "Alleles w/ >1 stop codons",
        label: "Alleles w/ >1 stop codons",
        options: {
          filter: true,
          sort: true,
          display: true,
          setCellHeaderProps: (value) => {
            return {
              style: {
                fontWeight: "bold",
              },
            };
          },
        },
      },
      {
        name: "Alleles wo/ Start/Stop Codon",
        label: "Alleles wo/ Start/Stop Codon",
        options: {
          filter: true,
          sort: true,
          display: true,
          setCellHeaderProps: (value) => {
            return {
              style: {
                fontWeight: "bold",
              },
            };
          },
        },
      },
      {
        name: "CDS",
        label: "CDS",
        options: {
          filter: true,
          sort: true,
          display: true,
          setCellHeaderProps: (value) => {
            return {
              style: {
                fontWeight: "bold",
              },
            };
          },
        },
      },
    ];

    const options = {
      responsive: "vertical",
      selectableRowsHeader: false,
      selectableRows: "none",
      selectableRowsOnClick: false,
      print: false,
      download: true,
      filter: true,
      search: true,
      viewColumns: true,
      pagination: true,
      onCellClick: (cellData, cellMeta) => {
        console.log(cellData, cellMeta);

        if (cellData.includes(".fasta")) {
          this.setState({
            scatterSelectOption: cellData,
            isIndOption: true,
          });
        }
      },
    };

    const cds_data = this.state.cds_df_data;

    const cds_analysis_table = (
      <MuiThemeProvider theme={this.getMuiTheme()}>
        <MUIDataTable
          title={""}
          data={cds_data}
          columns={columns}
          options={options}
        />
      </MuiThemeProvider>
    );

    // CDS scatter
    let scatterCDS = [];

    const barWidth = Array(
      this.state.cds_scatter_data["CDS_Alleles"].length
    ).fill(0.2);

    scatterCDS.push(
      {
        type: "bar",
        name: "CDS Alleles",
        x: this.state.cds_scatter_data["CDS_Alleles"],
        orientation: "h",
        text: this.state.cds_scatter_data["genes"],
        marker_color: "#045a8d",
        width: barWidth,
      },
      {
        type: "bar",
        name: "Non multiple 3",
        x: this.state.cds_scatter_data["mult3"],
        orientation: "h",
        text: this.state.cds_scatter_data["genes"],
        marker_color: "#006d2c",
        width: barWidth,
      },
      {
        type: "bar",
        name: "> 1 stop codon",
        x: this.state.cds_scatter_data["stopC"],
        orientation: "h",
        text: this.state.cds_scatter_data["genes"],
        marker_color: "#7b3294",
        width: barWidth,
      },
      {
        type: "bar",
        name: "No Start/Stop codon",
        x: this.state.cds_scatter_data["notStart"],
        orientation: "h",
        text: this.state.cds_scatter_data["genes"],
        marker_color: "#ec7014",
        width: barWidth,
      }
    );

    let cds_scatter_plot = (
      <Plot
        data={scatterCDS}
        layout={{
          barmode: "stack",
          bargap: 0.5,
          bargroupgap: 0.0,
          title: {
            text: "Summary of problematic alleles per locus",
          },
          xaxis: {
            title: { text: "Number of occurrences" },
            gridcolor: "#eee",
          },
          yaxis: { showgrid: false, showticklabels: false },
          plot_bgcolor: "rgba(0,0,0,0)",
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
      />
    );

    return (
      <Aux>
        <div>
          <div style={{ marginTop: "40px" }}>
            <Accordion defaultExpanded>
              <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                <Typography variant="h5" className={classes.title}>
                  Schema Evaluation
                </Typography>
              </AccordionSummary>
              <AccordionDetails>
                <div
                  className={classes.mainPaper}
                  style={{ width: "100%", height: "100%" }}
                >
                  <div style={style.buttonBar}>
                    <Button
                      style={style.button}
                      className={`${
                        this.state.tabValue === 0 && classes.tabButton
                      }`}
                      onClick={() => {
                        this.plotChangeHandler(0);
                      }}
                    >
                      Allele Numbers Analysis
                    </Button>
                    <Button
                      style={style.button}
                      className={`${
                        this.state.tabValue === 1 && classes.tabButton
                      }`}
                      onClick={() => {
                        this.plotChangeHandler(1);
                      }}
                    >
                      Allele Length Analysis
                    </Button>
                    <Button
                      style={style.button}
                      className={`${
                        this.state.tabValue === 2 && classes.tabButton
                      }`}
                      onClick={() => {
                        this.plotChangeHandler(2);
                      }}
                    >
                      Locus Statistics
                    </Button>
                    <Button
                      style={style.button}
                      className={`${
                        this.state.tabValue === 3 && classes.tabButton
                      }`}
                      onClick={() => {
                        this.plotChangeHandler(3);
                      }}
                    >
                      Locus Size Variation
                    </Button>
                  </div>
                  {this.state.tabValue === 0 && total_allele_plot}
                  {this.state.tabValue === 1 && mode_data_plot}
                  {this.state.tabValue === 2 && scatter_plot}
                  {this.state.tabValue === 3 && boxplot}
                </div>
              </AccordionDetails>
            </Accordion>
          </div>
          <div id="cds-analysis" style={{ marginTop: "40px" }}>
            <Typography variant="h5" className={classes.title}>
              CDS Analysis
            </Typography>
            <div style={{ marginTop: "20px" }}>{cds_analysis_table}</div>
            <div style={{ marginTop: "20px" }}>{cds_scatter_plot}</div>
          </div>
          <div id="locus-ind" style={{ marginTop: "40px" }}>
            <div style={{ marginTop: "20px", marginBottom: "20px" }}>
              <Typography variant="h5" className={classes.title}>
                Please choose the locus to analyse in the Locus Statistics or
                Locus Size Variation plot.
              </Typography>
            </div>
          </div>
        </div>
      </Aux>
    );
  }
}

export default SchemaEvaluator;