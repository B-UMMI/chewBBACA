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
import AccordionSummary from "@material-ui/core/AccordionSummary";
import AccordionDetails from "@material-ui/core/AccordionDetails";

// Material-UI Lab components
import Alert from "@material-ui/lab/Alert";

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

// SeqLogo component
import { ProteinLogo } from "logojs-react";

class SchemaEvaluator extends Component {
  state = {
    locus_ind_data: _preComputedDataInd,
    cds_df_data: _cdsDf,
    exceptions: _exceptions,
    testMSA: _msaData,
    phyloData: _phyloData,
    minLen: _minLen,
    fasta: _fasta,
    indTabValue: 0,
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

    const phylocanvasComponent =
      this.state.phyloData === "undefined" ? (
        <div style={{ marginTop: "20px", width: "100%" }}>
          <Alert variant="outlined" severity="warning">
            <Typography variant="subtitle1">
              The NJ tree and MSA were not generated because chewBBACA
              determined that this locus does not have valid alleles or has only
              one valid allele.
            </Typography>
          </Alert>
        </div>
      ) : (
        <Aux>
          <div style={{ marginTop: "40px" }}>
            <div style={{ width: "100%" }}>
              <Alert variant="outlined" severity="info">
                <Typography variant="subtitle1">
                  Note: the first number of the name of each leaf corresponds to
                  the order the sequence appears in the input sequence file. For
                  example: '1_5' means that the first sequence of the input file
                  is named '5'.
                </Typography>
              </Alert>
            </div>
            <div style={{ marginTop: "10px", width: "100%" }}>
              <Alert variant="outlined" severity="info">
                <Typography variant="subtitle1">
                  If a tree is not being displayed it means that the distance
                  between alleles is 0 and Phylocanvas won't render the tree.
                </Typography>
              </Alert>
            </div>
          </div>
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
                      label={`Zoom is ${
                        this.state.zoom ? "enabled" : "disabled"
                      }`}
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
        </Aux>
      );

    const seqLogo =
      this.state.fasta === "undefined" ? (
        <div style={{ marginTop: "20px", width: "100%" }}>
          <Alert variant="outlined" severity="warning">
            <Typography variant="subtitle1">
              No valid alignment.
            </Typography>
          </Alert>
        </div>
      ) : (
        <Aux>
          <div style={{ marginTop: "40px" }}>
            <Accordion defaultExpanded>
              <AccordionSummary>
                <Typography variant="h5" className={classes.title}>
                  Sequence Logo
                </Typography>
              </AccordionSummary>
              <AccordionDetails>
                <div id="logo" style={{ width: "100%", height: "100%" }}>
                  <ProteinLogo
                    fasta={this.state.fasta}
                    noFastaNames
                    mode={"FREQUENCY"}
                  />
                </div>
              </AccordionDetails>
            </Accordion>
          </div>
        </Aux>
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
            range: [
              Math.min(...hist_allele_sizes_x) -
                (this.state.locus_ind_data.data.alleles_mode / 100) * 20,
              Math.max(...hist_allele_sizes_x) +
                (this.state.locus_ind_data.data.alleles_mode / 100) * 20,
            ],
          },
          yaxis: {
            title: { text: "Number of Alleles" },
          },
          shapes: [
            {
              line: { color: "orange", width: 3 },
              type: "line",
              x0: this.state.minLen,
              x1: this.state.minLen,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
            {
              line: { color: "red", width: 3 },
              type: "line",
              x0: this.state.minLen - (this.state.minLen / 100) * 20,
              x1: this.state.minLen - (this.state.minLen / 100) * 20,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
            {
              fillcolor: "orange",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: this.state.minLen,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
            {
              fillcolor: "red",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: this.state.minLen - (this.state.minLen / 100) * 20,
              xref: "x",
              y0: 0,
              y1: 1,
              yref: "y domain",
            },
          ],
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
          shapes: [
            {
              line: { color: "orange", width: 3 },
              type: "line",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: this.state.minLen,
              y1: this.state.minLen,
              yref: "y",
            },
            {
              line: { color: "red", width: 3 },
              type: "line",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: this.state.minLen - (this.state.minLen / 100) * 20,
              y1: this.state.minLen - (this.state.minLen / 100) * 20,
              yref: "y",
            },
            {
              fillcolor: "orange",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: 0,
              y1: this.state.minLen,
              yref: "y",
            },
            {
              fillcolor: "red",
              line: { width: 0 },
              opacity: 0.1,
              type: "rect",
              x0: 0,
              x1: 1,
              xref: "x domain",
              y0: 0,
              y1: this.state.minLen - (this.state.minLen / 100) * 20,
              yref: "y",
            },
          ],
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

    const alleleShorterColumn = Object.keys(this.state.cds_df_data).find((k) =>
      k.includes("Alleles shorter")
    );

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
        name: alleleShorterColumn,
        label: alleleShorterColumn,
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
        name: "Missing Allele IDs",
        label: "Missing Allele IDs",
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
          customBodyRender: (value, tableMeta, updateValue) => (
            <div>{value.join(", ")}</div>
          ),
        },
      },
      {
        name: "CDS",
        label: "Valid Alleles",
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
      pagination: true,
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

    const msa_component =
      this.state.testMSA === "undefined" ? (
        <div />
      ) : (
        <div style={{ marginTop: "40px" }}>
          <Accordion defaultExpanded>
            <AccordionSummary>
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
            {seqLogo}
            {msa_component}
          </div>
        </div>
      </Aux>
    );
  }
}

export default SchemaEvaluator;
