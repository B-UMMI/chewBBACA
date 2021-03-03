import React, { Component } from "react";

import Aux from "../hoc/Aux";
import classes from "./SchemaEvaluator.css";

// Material-UI components
import Button from "@material-ui/core/Button";
import Typography from "@material-ui/core/Typography";
import { createMuiTheme, MuiThemeProvider } from "@material-ui/core/styles";

// Material-UI ExpansionPanel components
import Accordion from "@material-ui/core/Accordion";
import ExpandMoreIcon from "@material-ui/icons/ExpandMore";
import AccordionSummary from "@material-ui/core/AccordionSummary";
import AccordionDetails from "@material-ui/core/AccordionDetails";

// Material-UI Lab components
import Alert from "@material-ui/lab/Alert";

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
    total_data: _totalData,
    notConserved: _notConserved,
    notConservedMessage: _notConservedMessage,
    oneAlleleOnly: _oneAlleleOnly,
    message: _message,
    tabValue: 0,
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
    anchor.href = `../html_files/${
      locus_id.split(".")[0]
    }_individual_report.html`;
    anchor.target = "_blank";
    anchor.rel = "noopener noreferrer";
    anchor.click();
  };

  clickBoxPlotHandler = (event) => {
    const locus_id = event.points[0].x;
    console.log(locus_id);

    const anchor = document.createElement("a");
    anchor.href = `../html_files/${
      locus_id.split(".")[0]
    }_individual_report.html`;
    anchor.target = "_blank";
    anchor.rel = "noopener noreferrer";
    anchor.click();
  };

  // clickCdsTableHandler = (cellData) => {
  //   const locus_id = cellData;
  //   console.log(locus_id);

  //   const anchor = document.createElement("a");
  //   anchor.href = `../html_files/${
  //     locus_id.split(".")[0]
  //   }_individual_report.html`;
  //   anchor.target = "_blank";
  //   anchor.rel = "noopener noreferrer";
  //   anchor.click();
  // };

  totalDataColumnsHandler = (columnData, alleleShorterColumn) => {
    if ("chewBBACA_version" in columnData[0]) {
      const totalDataColumns = [
        {
          name: "chewBBACA_version",
          label: "chewBBACA Version",
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
          name: "bsr",
          label: "Blast Score Ratio",
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
          name: "total_loci",
          label: "Total Loci",
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
          name: "total_alleles",
          label: "Total Alleles",
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
          name: "total_alleles_mult3",
          label: "Total Alleles not multiple of 3",
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
          name: "total_alleles_stopC",
          label: "Total Alleles w/ >1 stop codons",
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
          name: "total_alleles_notStart",
          label: "Total Alleles wo/ Start/Stop Codon",
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
          name: "total_alleles_shorter",
          label: `Total ${alleleShorterColumn}`,
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
          name: "total_invalid_alleles",
          label: "Total Invalid Alleles",
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
      return totalDataColumns;
    } else {
      const totalDataColumns = [
        {
          name: "total_loci",
          label: "Total Loci",
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
          name: "total_alleles",
          label: "Total Alleles",
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
          name: "total_alleles_mult3",
          label: "Total Alleles not multiple of 3",
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
          name: "total_alleles_stopC",
          label: "Total Alleles w/ >1 stop codons",
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
          name: "total_alleles_notStart",
          label: "Total Alleles wo/ Start/Stop Codon",
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
          name: "total_alleles_shorter",
          label: `Total ${alleleShorterColumn}`,
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
          name: "total_invalid_alleles",
          label: "Total Invalid Alleles",
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
      return totalDataColumns;
    }
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
    let boxLoci = [];
    let boxQ1 = [];
    let boxMedian = [];
    let boxQ3 = [];
    let boxLowerfence = [];
    let boxUpperfence = [];

    for (let boxKey in this.state.pre_computed_data_boxplot["boxplot_data"]) {
      boxLoci.push(
        this.state.pre_computed_data_boxplot["boxplot_data"][boxKey].locus_name
      );
      boxQ1.push(
        this.state.pre_computed_data_boxplot["boxplot_data"][boxKey].q1
      );
      boxMedian.push(
        this.state.pre_computed_data_boxplot["boxplot_data"][boxKey].median
      );
      boxQ3.push(
        this.state.pre_computed_data_boxplot["boxplot_data"][boxKey].q3
      );
      boxLowerfence.push(
        this.state.pre_computed_data_boxplot["boxplot_data"][boxKey].min
      );
      boxUpperfence.push(
        this.state.pre_computed_data_boxplot["boxplot_data"][boxKey].max
      );
    }

    let boxplotData = [];

    boxplotData.push({
      type: "box",
      name: "Locus Size Variation",
      x: boxLoci,
      q1: boxQ1,
      median: boxMedian,
      q3: boxQ3,
      lowerfence: boxLowerfence,
      upperfence: boxUpperfence,
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
    const alleleShorterColumn = Object.keys(
      this.state.cds_df_data[0]
    ).find((k) => k.includes("Alleles shorter"));

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
          customBodyRender: (value, tableMeta, updateValue) => (
            <a
              href={`../html_files/${
                value.split(".")[0]
              }_individual_report.html`}
              target={"_blank"}
              rel={"noopener noreferrer"}
            >
              {value}
            </a>
          ),
        },
      },
      {
        name: "name",
        label: "Uniprot Annotation",
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
        name: "url",
        label: "Uniprot URL",
        options: {
          filter: true,
          sort: true,
          display: false,
          setCellHeaderProps: (value) => {
            return {
              style: {
                fontWeight: "bold",
              },
            };
          },
          customBodyRender: (value, tableMeta, updateValue) => (
            <a href={value} target={"_blank"} rel={"noopener noreferrer"}>
              {value}
            </a>
          ),
        },
      },
      {
        name: "genome",
        label: "Origin Genome",
        options: {
          filter: true,
          sort: true,
          display: false,
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
        name: "contig",
        label: "Origin Genome Contig",
        options: {
          filter: true,
          sort: true,
          display: false,
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
        name: "start",
        label: "Original Genome Start Position",
        options: {
          filter: true,
          sort: true,
          display: false,
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
        name: "stop",
        label: "Original Genome Stop Position",
        options: {
          filter: true,
          sort: true,
          display: false,
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
        name: "coding_strand",
        label: "Coding Strand",
        options: {
          filter: true,
          sort: true,
          display: false,
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
        name: "Total Invalid Alleles",
        label: "Total Invalid Alleles",
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
      // onCellClick: (cellData, cellMeta) => {
      //   if (cellData.includes(".fasta")) {
      //     this.clickCdsTableHandler(cellData);
      //   }
      // },
    };

    const cds_data = this.state.cds_df_data;

    const cds_analysis_table = (
      <MuiThemeProvider theme={this.getMuiTheme()}>
        <MUIDataTable
          title={"Loci Analysis"}
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
      },
      {
        type: "bar",
        name: alleleShorterColumn,
        x: this.state.cds_scatter_data["shorter"],
        orientation: "h",
        text: this.state.cds_scatter_data["genes"],
        marker_color: "#6a3d9a",
        width: barWidth,
      }
    );

    let cds_scatter_plot = (
      <Plot
        data={scatterCDS}
        layout={{
          autosize: true,
          barmode: "stack",
          bargap: 0.5,
          bargroupgap: 0.0,
          title: {
            text: "Summary of problematic alleles per locus",
          },
          xaxis: {
            title: { text: "Number of occurrences" },
            gridcolor: "#eee",
            type: "log",
          },
          yaxis: { showgrid: false, showticklabels: false },
          plot_bgcolor: "rgba(0,0,0,0)",
          hovermode: "closest",
        }}
        useResizeHandler={true}
        style={{ width: "100%", height: "700px" }}
      />
    );

    // Schema Summary Statistics table
    const totalDataColumns = this.totalDataColumnsHandler(
      this.state.total_data,
      alleleShorterColumn
    );

    const totalDataOptions = {
      responsive: "vertical",
      selectableRowsHeader: false,
      selectableRows: "none",
      selectableRowsOnClick: false,
      print: false,
      download: false,
      filter: false,
      search: false,
      viewColumns: true,
      pagination: false,
    };

    const total_data = this.state.total_data;

    const total_data_table = (
      <MuiThemeProvider theme={this.getMuiTheme()}>
        <MUIDataTable
          title={"Schema Summary Statistics"}
          data={total_data}
          columns={totalDataColumns}
          options={totalDataOptions}
        />
      </MuiThemeProvider>
    );

    let notConservedList = <div />;

    if (this.state.notConserved === "undefined") {
      notConservedList = (
        <Alert variant="outlined" severity="info">
          No loci with high length variability detected.
        </Alert>
      );
    } else {
      const notConservedColumns = [
        {
          name: "gene",
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
            customBodyRender: (value, tableMeta, updateValue) => (
              <a
                href={`../html_files/${
                  value.split(".")[0]
                }_individual_report.html`}
                target={"_blank"}
                rel={"noopener noreferrer"}
              >
                {value}
              </a>
            ),
          },
        },
      ];

      const notConservedOptions = {
        responsive: "vertical",
        selectableRowsHeader: false,
        selectableRows: "none",
        selectableRowsOnClick: false,
        print: false,
        download: false,
        filter: false,
        search: false,
        viewColumns: true,
        pagination: true,
        // onCellClick: (cellData, cellMeta) => {
        //   if (cellData.includes(".fasta")) {
        //     this.clickCdsTableHandler(cellData);
        //   }
        // },
      };

      notConservedList = (
        <Aux>
          <div style={{ marginTop: "40px", width: "100%" }}>
            <Alert variant="outlined" severity="info">
              {this.state.notConservedMessage}
            </Alert>
          </div>
          <div style={{ marginTop: "20px" }}>
            <MuiThemeProvider theme={this.getMuiTheme()}>
              <MUIDataTable
                title={"Loci with high variability"}
                data={this.state.notConserved}
                columns={notConservedColumns}
                options={notConservedOptions}
              />
            </MuiThemeProvider>
          </div>
        </Aux>
      );
    }

    let oneAlleleOnlyList = <div />;

    if (this.state.oneAlleleOnly === "undefined") {
      oneAlleleOnlyList = (
        <Alert variant="outlined" severity="info">
          No loci with only 1 allele were detected.
        </Alert>
      );
    } else {
      const oneAlleleOnlyColumns = [
        {
          name: "gene",
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
            customBodyRender: (value, tableMeta, updateValue) => (
              <a
                href={`../html_files/${
                  value.split(".")[0]
                }_individual_report.html`}
                target={"_blank"}
                rel={"noopener noreferrer"}
              >
                {value}
              </a>
            ),
          },
        },
      ];

      const oneAlleleOnlyOptions = {
        responsive: "vertical",
        selectableRowsHeader: false,
        selectableRows: "none",
        selectableRowsOnClick: false,
        print: false,
        download: false,
        filter: false,
        search: false,
        viewColumns: true,
        pagination: true,
        // onCellClick: (cellData, cellMeta) => {
        //   if (cellData.includes(".fasta")) {
        //     this.clickCdsTableHandler(cellData);
        //   }
        // },
      };

      oneAlleleOnlyList = (
        <Aux>
          <div>
            <MuiThemeProvider theme={this.getMuiTheme()}>
              <MUIDataTable
                title={"Loci with only 1 allele"}
                data={this.state.oneAlleleOnly}
                columns={oneAlleleOnlyColumns}
                options={oneAlleleOnlyOptions}
              />
            </MuiThemeProvider>
          </div>
        </Aux>
      );
    }

    let minLenTableColumns = [
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
          customBodyRender: (value, tableMeta, updateValue) => (
            <a
              href={`../html_files/${
                value.split(".")[0]
              }_individual_report.html`}
              target={"_blank"}
              rel={"noopener noreferrer"}
            >
              {value}
            </a>
          ),
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
    ];

    let minLenTableOptions = {
      responsive: "vertical",
      selectableRowsHeader: false,
      selectableRows: "none",
      selectableRowsOnClick: false,
      print: false,
      download: false,
      filter: false,
      search: false,
      viewColumns: true,
      pagination: true,
      // onCellClick: (cellData, cellMeta) => {
      //   if (cellData.includes(".fasta")) {
      //     this.clickCdsTableHandler(cellData);
      //   }
      // },
    };

    let minLenTable = (
      <Aux>
        <div>
          <MuiThemeProvider theme={this.getMuiTheme()}>
            <MUIDataTable
              title={`Loci with ${alleleShorterColumn}`}
              data={cds_data}
              columns={minLenTableColumns}
              options={minLenTableOptions}
            />
          </MuiThemeProvider>
        </div>
      </Aux>
    );

    return (
      <Aux>
        <div>
          <div style={{ marginTop: "40px", width: "100%" }}>
            <Alert variant="outlined" severity="info">
              {this.state.message}
            </Alert>
          </div>
          <div style={{ marginTop: "40px" }}>{total_data_table}</div>
          <div style={{ marginTop: "40px" }}>{notConservedList}</div>
          <div style={{ marginTop: "40px" }}>{oneAlleleOnlyList}</div>
          <div style={{ marginTop: "40px" }}>{minLenTable}</div>
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
                  <div style={{ marginBottom: "10px" }}>
                    <Typography variant="subtitle1">
                      In the Locus Statistics and Locus Size Variation graphs,
                      each point (locus) is clickable and will open a page with
                      more details about it.
                    </Typography>
                  </div>
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
              Loci Analysis
            </Typography>
            <div style={{ marginTop: "20px" }}>{cds_analysis_table}</div>
            <div style={{ marginTop: "20px" }}>
              <Aux>
                <div style={{ marginTop: "20px" }}>
                  <Alert variant="outlined" severity="info">
                    <Typography variant="subtitle1">
                      The x-axis of the following plot is in logarithmic scale
                      for a more compact visualization.
                    </Typography>
                  </Alert>
                </div>
                {cds_scatter_plot}
              </Aux>
            </div>
          </div>
        </div>
      </Aux>
    );
  }
}

export default SchemaEvaluator;
