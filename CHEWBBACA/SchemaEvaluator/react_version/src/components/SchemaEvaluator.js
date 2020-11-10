import React, { Component } from "react";
import pre_computed_data from "../data/pre_computed_data.json";
import pre_computed_data_boxplot from "../data/pre_computed_data_boxplot.json";
import cds_df_data from "../data/cds_df.json";
import cds_scatter_data from "../data/cds_scatter.json";
import locus_ind_data from "../data/pre_computed_data_ind.json";
// import alignment from "../data/AAATEST.fasta";

import Aux from "../hoc/Aux";
import classes from "./SchemaEvaluator.css";
// import classNames from "classnames";

// Material-UI components
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import {
  createMuiTheme,
  MuiThemeProvider,
  withStyles,
} from "@material-ui/core/styles";

import InputLabel from "@material-ui/core/InputLabel";
import MenuItem from "@material-ui/core/MenuItem";
import FormHelperText from "@material-ui/core/FormHelperText";
import FormControl from "@material-ui/core/FormControl";
import Select from "@material-ui/core/Select";

// Material-UI ExpansionPanel components
import Accordion from "@material-ui/core/Accordion";
import ExpandMoreIcon from "@material-ui/icons/ExpandMore";
import AccordionSummary from "@material-ui/core/AccordionSummary";
import AccordionDetails from "@material-ui/core/AccordionDetails";

// Material-UI Datatables
import MUIDataTable from "mui-datatables";

// Plotly.js
import Plot from "react-plotly.js";

// MSAViewer
// import AlignmentViewer from "react-alignment-viewer/lib/components/AlignmentViewer";
// import AlignmentChart from "react-alignment-viewer/lib/components/AlignmentChart";

// const styles = (theme) => ({
//   formControl: {
//     margin: 0,
//     display: "flex",
//     wrap: "nowrap",
//   },
// });

class SchemaEvaluator extends Component {
  state = {
    tabValue: 0,
    dropdownOption: '',
  };

  // componentDidMount() {
  //   console.log(alignment);
  // }

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

    // console.log(Object.keys(locus_ind_data));

    const hist_allele_sizes_x = locus_ind_data[locus_id]["allele_sizes"];

    const hist_allele_ids_y = locus_ind_data[locus_id]["locus_ids"];

    let ind_hist_data = [];

    ind_hist_data.push({
      x: hist_allele_sizes_x,
      y: hist_allele_ids_y,
      type: "histogram",
      name: "Locus Details",
    });

    let ind_hist = (
      <Plot
        data={ind_hist_data}
        layout={{
          xaxis: {
            title: { text: "Sequence size in bp" },
          },
          yaxis: {
            title: { text: "Number of Alleles" },
          },
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        line={{
          width: 1,
        }}
      />
    );

    ind_locus = <div>{ind_hist}</div>;
  };

  onSelectChange = (event) => {
    console.log(event.target.value);
    this.setState({ dropdownOption: event.target.value });
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
    };

    // Create panel A
    let locus_name2 = [];
    let nr_alleles = [];
    let total_al_data = [];
    for (let key in pre_computed_data.total_alleles) {
      locus_name2.push(pre_computed_data.total_alleles[key].locus_name);
      nr_alleles.push(pre_computed_data.total_alleles[key].nr_alleles);
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

    for (let key in pre_computed_data.mode) {
      allele_mode.push(pre_computed_data.mode[key].alleles_mode);
      locus_name.push(pre_computed_data.mode[key].locus_name);
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

    for (let key in pre_computed_data.scatter_data) {
      locus_id.push(pre_computed_data.scatter_data[key].locus_id);
      nr_alleles_scatter.push(pre_computed_data.scatter_data[key].nr_alleles);
      scatter_data_median.push(
        pre_computed_data.scatter_data[key].alleles_median
      );
      scatter_data_min.push(pre_computed_data.scatter_data[key].alleles_min);
      scatter_data_max.push(pre_computed_data.scatter_data[key].alleles_max);
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
      x: pre_computed_data_boxplot.loci,
      q1: pre_computed_data_boxplot.q1,
      median: pre_computed_data_boxplot.median,
      q3: pre_computed_data_boxplot.q3,
      lowerfence: pre_computed_data_boxplot.min,
      upperfence: pre_computed_data_boxplot.max,
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
    };

    const cds_data = cds_df_data;

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

    scatterCDS.push(
      {
        type: "bar",
        name: "CDS Alleles",
        x: cds_scatter_data["CDS_Alleles"],
        orientation: "h",
        text: cds_scatter_data["genes"],
        marker_color: "#045a8d",
      },
      {
        type: "bar",
        name: "Non multiple 3",
        x: cds_scatter_data["mult3"],
        orientation: "h",
        text: cds_scatter_data["genes"],
        marker_color: "#006d2c",
      },
      {
        type: "bar",
        name: "> 1 stop codon",
        x: cds_scatter_data["stopC"],
        orientation: "h",
        text: cds_scatter_data["genes"],
        marker_color: "#7b3294",
      },
      {
        type: "bar",
        name: "No Start/Stop codon",
        x: cds_scatter_data["notStart"],
        orientation: "h",
        text: cds_scatter_data["genes"],
        marker_color: "#ec7014",
      }
    );

    let cds_scatter_plot = (
      <Plot
        data={scatterCDS}
        layout={{
          barmode: "stack",
          bargap: 0.0,
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

    const dropdownOptions = Object.keys(locus_ind_data);

    let ind_locus = (
      <FormControl style={style.formControl}>
        <InputLabel id="locus-ind-select-label">
          Please select the locus to analyse.
        </InputLabel>
        <Select
          labelId="locus-ind-select-label"
          id="locus-ind-select"
          autoWidth={true}
          value={this.state.dropdownOption}
          onChange={(e) => this.onSelectChange(e)}
        >
          {dropdownOptions.map((option) => (
            <MenuItem key={option} value={option}>{option}</MenuItem>
          ))}
        </Select>
      </FormControl>
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
            {/* <div id="msa-viewer" style={{ marginTop: "20px" }}>
              <AlignmentChart
                data={alignment}
                height={4000}
                showconservation={false}
                showgap={false}
              />
            </div> */}
          </div>
          <div id="locus-ind" style={{ marginTop: "40px" }}>
            <div style={{ marginTop: "20px" }}>
              <Typography variant="h5" className={classes.title}>
                Locus Individual Analysis
              </Typography>
            </div>
            {/* <div style={{ marginTop: "20px" }}>
              <Typography variant="subtitle1">
                Please select the locus to analyse.
              </Typography>
            </div> */}
            <div style={{ marginTop: "20px" }}>{ind_locus}</div>
          </div>
        </div>
      </Aux>
    );
  }
}

export default SchemaEvaluator;
