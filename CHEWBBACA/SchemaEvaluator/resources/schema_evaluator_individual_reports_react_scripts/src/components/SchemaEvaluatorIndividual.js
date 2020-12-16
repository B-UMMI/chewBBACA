import React, { Component } from "react";

import Aux from "../hoc/Aux";
import classes from "./SchemaEvaluatorIndividual.css";

// Material-UI components
import Grid from "@material-ui/core/Grid";
import Switch from "@material-ui/core/Switch";
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import { createMuiTheme, MuiThemeProvider } from "@material-ui/core/styles";

// Material-UI ExpansionPanel components
import Accordion from "@material-ui/core/Accordion";
import ExpandMoreIcon from "@material-ui/icons/ExpandMore";
import AccordionSummary from "@material-ui/core/AccordionSummary";
import AccordionDetails from "@material-ui/core/AccordionDetails";

// react-select import
import Select from "react-select";

// Material-UI Datatables
import MUIDataTable from "mui-datatables";

// Plotly.js
import Plot from "react-plotly.js";

// MSAViewer
import MSAViewer from "react-msa-viewer";

// Phylocanvas
import { PhylogeneticTree } from "./Phylocanvas";

class SchemaEvaluator extends Component {
  state = {
    locus_ind_data: _preComputedDataInd,
    cds_df_data: _cdsDf,
    exceptions: _exceptions,
    testMSA: _msaData,
    phyloData: _phyloData,
    indTabValue: 0,
    scatterSelectOption: "",
    isIndOption: false,
    showSnack: false,
    treeType: 0,
    zoom: false,
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

  indPlotChangeHandler = (value) => {
    this.setState({ indTabValue: value });
  };

  handleSelectChange = (selectType, val) => {
    const newState = {};
    newState[selectType] = val.value;

    this.setState(newState);
  };

  render() {
    // Tree type options for Phylocanvas
    const treeTypeOption = [
      "circular",
      "rectangular",
      "diagonal",
      "hierarchical",
      "radial",
    ].map((v, i) => {
      return { label: v, value: i };
    });

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
      select: {
        width: "170px",
        marginRight: "10px",
        zIndex: 100,
      },
      zoom: {
        marginLeft: "auto",
      },
    };

    const phylocanvasComponent = (
      <div style={{ marginTop: "40px" }}>
        <Accordion defaultExpanded>
          <AccordionSummary>
            <Typography variant="h5" className={classes.title}>
              NJ tree
            </Typography>
          </AccordionSummary>
          <AccordionDetails>
            <div style={style.root}>
              <Grid container style={style.toolbar}>
                <div style={style.select}>
                  <Typography>Tree type:</Typography>
                  <Select
                    closeMenuOnSelect={false}
                    options={treeTypeOption}
                    onChange={(val) => {
                      this.handleSelectChange("treeType", val);
                    }}
                    value={treeTypeOption[this.state.treeType]}
                  />
                </div>
                <FormControlLabel
                  control={
                    <Switch
                      color={"primary"}
                      onChange={() => {
                        this.setState({ zoom: !this.state.zoom });
                      }}
                      checked={this.state.zoom}
                    />
                  }
                  label={`Zoom is ${this.state.zoom ? "enabled" : "disabled"}`}
                  style={style.zoom}
                />
              </Grid>
              <PhylogeneticTree
                zoom={this.state.zoom}
                treeType={treeTypeOption[this.state.treeType].label}
                newickString={this.state.phyloData.phylo_data}
              />
            </div>
          </AccordionDetails>
        </Accordion>
      </div>
    );

    const hist_allele_sizes_x = this.state.locus_ind_data.data.allele_sizes;

    const hist_allele_ids_y = this.state.locus_ind_data.data.locus_ids;

    let ind_hist_data = [];

    ind_hist_data.push({
      x: hist_allele_sizes_x,
      y: hist_allele_ids_y,
      type: "histogram",
      name: "Locus Details",
    });

    const ind_hist = (
      <Plot
        data={ind_hist_data}
        layout={{
          title: {
            text: this.state.locus_ind_data.locus_name,
          },
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

    // Build scatter
    const scatter_allele_ids_x = this.state.locus_ind_data.data.locus_ids.map(
      (id) => id.split("_")[id.split("_").length - 1]
    );

    const scatter_allele_sizes_y = this.state.locus_ind_data.data.allele_sizes;

    let ind_scatter_data = [];

    ind_scatter_data.push({
      x: scatter_allele_ids_x,
      y: scatter_allele_sizes_y,
      type: "scatter",
      name: "Locus Details",
      mode: "markers",
    });

    const ind_scatter = (
      <Plot
        data={ind_scatter_data}
        layout={{
          title: {
            text: this.state.locus_ind_data.locus_name,
          },
          xaxis: {
            title: { text: "Allele ID", tick0: 0, dtick: 1 },
          },
          yaxis: {
            title: { text: "Sequence size in bp", tick0: 0, dtick: 1 },
          },
          hovermode: "closest",
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "100%" }}
        line={{
          width: 1,
        }}
      />
    );

    // Build table with additional information
    const cdsTableInfo = this.state.cds_df_data;

    const indData = {
      size_range: this.state.locus_ind_data.data.size_range,
      alleles_median: this.state.locus_ind_data.data.alleles_median,
      alleles_mode: this.state.locus_ind_data.data.alleles_mode,
    };

    const indTableData = [{ ...cdsTableInfo, ...indData }];

    const indTableColumns = [
      {
        name: "Gene",
        label: "Locus",
        options: {
          filter: false,
          sort: false,
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
        name: "Number of alleles",
        label: "Number of Alleles",
        options: {
          filter: false,
          sort: false,
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
          filter: false,
          sort: false,
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
          filter: false,
          sort: false,
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
          filter: false,
          sort: false,
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
          filter: false,
          sort: false,
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
        name: "size_range",
        label: "Size Range (bp)",
        options: {
          filter: false,
          sort: false,
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
        name: "alleles_median",
        label: "Alleles Median",
        options: {
          filter: false,
          sort: false,
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
        name: "alleles_mode",
        label: "Alleles Mode",
        options: {
          filter: false,
          sort: false,
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

    const indTableOptions = {
      responsive: "vertical",
      selectableRowsHeader: false,
      selectableRows: "none",
      selectableRowsOnClick: false,
      print: false,
      download: false,
      filter: false,
      search: false,
      viewColumns: false,
      pagination: false,
    };

    const indTable = (
      <MuiThemeProvider theme={this.getMuiTheme()}>
        <MUIDataTable
          title={"Locus information"}
          data={indTableData}
          columns={indTableColumns}
          options={indTableOptions}
        />
      </MuiThemeProvider>
    );

    const excTableColumns = [
      {
        name: "allele",
        label: "Allele",
        options: {
          filter: false,
          sort: false,
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
        name: "exception",
        label: "Exceptions",
        options: {
          filter: false,
          sort: false,
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

    const excTableOptions = {
      responsive: "vertical",
      selectableRowsHeader: false,
      selectableRows: "none",
      selectableRowsOnClick: false,
      print: false,
      download: false,
      filter: false,
      search: false,
      viewColumns: false,
      pagination: false,
    };

    const excTable = (
      <MuiThemeProvider theme={this.getMuiTheme()}>
        <MUIDataTable
          title={"Exceptions"}
          data={this.state.exceptions}
          columns={excTableColumns}
          options={excTableOptions}
        />
      </MuiThemeProvider>
    );

    const msa_options = {
      ...this.state.testMSA,
      colorScheme: "clustal2",
      height: 1000,
      width: 1500,
      barMethod: "information-content",
    };

    const msa_component = (
      <div style={{ marginTop: "40px" }}>
        <Accordion defaultExpanded>
          <AccordionSummary expandIcon={<ExpandMoreIcon />}>
            <Typography variant="h5" className={classes.title}>
              Multiple Sequence Analysis
            </Typography>
          </AccordionSummary>
          <AccordionDetails>
            <div id="msa-viewer" style={{ padding: "16px" }}>
              <Typography variant="subtitle1">
                To analyse the alignment please hover over the alignment until
                the hand cursor appears and then drag it to see the remaining
                data.
              </Typography>
              <Typography variant="subtitle1">
                Note: the gray graph above the alignment represents the
                information entropy after Shannon of a column (scaled).
              </Typography>
              <br />
              <MSAViewer {...msa_options} />
            </div>
          </AccordionDetails>
        </Accordion>
      </div>
    );

    let locusIndHist = (
      <div
        style={{
          marginBottom: "2%",
          marginTop: "3%",
        }}
      >
        <div>
          <Accordion defaultExpanded>
            <AccordionDetails>
              <div
                className={classes.mainPaper}
                style={{ width: "100%", height: "100%" }}
              >
                <div style={style.buttonBar}>
                  <Button
                    style={style.button}
                    className={`${
                      this.state.indTabValue === 0 && classes.tabButton
                    }`}
                    onClick={() => {
                      this.indPlotChangeHandler(0);
                    }}
                  >
                    Histogram
                  </Button>
                  <Button
                    style={style.button}
                    className={`${
                      this.state.indTabValue === 1 && classes.tabButton
                    }`}
                    onClick={() => {
                      this.indPlotChangeHandler(1);
                    }}
                  >
                    Scatter
                  </Button>
                </div>
                {this.state.indTabValue === 0 && ind_hist}
                {this.state.indTabValue === 1 && ind_scatter}
              </div>
            </AccordionDetails>
          </Accordion>
        </div>
        <div style={{ marginTop: "20px" }}>{indTable}</div>
        <div style={{ marginTop: "20px" }}>{excTable}</div>
      </div>
    );

    return (
      <Aux>
        <div>
          <div id="locus-ind" style={{ marginTop: "40px" }}>
            <div style={{ marginTop: "20px" }}>
              <Typography variant="h5" className={classes.title}>
                Locus Individual Analysis
              </Typography>
            </div>
            {locusIndHist}
            {phylocanvasComponent}
            {msa_component}
          </div>
        </div>
      </Aux>
    );
  }
}

export default SchemaEvaluator;
