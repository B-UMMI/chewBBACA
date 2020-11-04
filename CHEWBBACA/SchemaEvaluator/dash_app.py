# import six.moves.urllib.request as urlreq
# from six import PY3

import os
import json

import plotly.graph_objects as go

import dash
import dash_table
import dash_bio as dashbio
import dash_core_components as dcc
import dash_html_components as html
# import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate

import SchemaEvaluator


# Functions to read pre-computed-data
def get_data(pre_computed_data_path):
    """Get pre-computed data for panels A-C."""
    json_file = os.path.join(pre_computed_data_path, "pre_computed_data.json")
    with open(json_file) as data_file:
        return json.load(data_file)


def get_boxplot_data(pre_computed_data_path):
    """Get pre-computed data for panel D."""
    boxplot_json_file = os.path.join(
        pre_computed_data_path, "pre_computed_data_boxplot.json")
    with open(boxplot_json_file) as data_file:
        return json.load(data_file)


def get_cds_df(schema_dir):
    """Creates the dataframe for panel E."""
    data, hist_data = SchemaEvaluator.create_cds_df(schema_dir)
    return data, hist_data


def get_locus_ind_dict(pre_computed_data_path):
    """Get pre-computed data for the Locus Invidual Analysis section."""
    locus_ind_file = os.path.join(
        pre_computed_data_path, "pre_computed_data_ind.json")
    with open(locus_ind_file) as data_file:
        return json.load(data_file)


def get_locus_ind_df(locus_ind_dict, option):
    """Creates a dataframe for the Locus Invidual Analysis section."""
    df = SchemaEvaluator.create_locus_ind_df(locus_ind_dict, option)
    return df


# Functions to create Plotly charts and Dash components
def build_locus_options(pre_computed_data_path):
    """Creates the options for the dropdown menu of the Locus Invidual Analysis section."""

    data_ind = get_locus_ind_dict(pre_computed_data_path)

    # Options for selectbox
    locus_options = list(data_ind.keys())

    options = []
    for k in locus_options:
        options.append({"label": k, "value": k})

    return options, data_ind


def build_locus_histogram(data_ind, locus_id):
    """Creates a plotly Histogram for the Locus Invidual Analysis section."""

    allele_sizes_x = data_ind[locus_id]["allele_sizes"]

    allele_ids_y = data_ind[locus_id]["locus_ids"]

    dict_of_fig = dict(
        {
            "data": [
                {
                    "type": "histogram",
                    "x": allele_sizes_x,
                    "y": allele_ids_y,
                    "name": "Locus Details",
                }
            ],
            "layout": {
                "xaxis": {
                    "title": {"text": "Sequence size in bp"},
                },
                "yaxis": {
                    "title": {"text": "Number of Alleles"},
                },
            },
        }
    )

    fig = go.Figure(dict_of_fig)

    return fig


def build_locus_scatter(data_ind, locus_id):
    """Creates a plotly Scatterplot for the Locus Invidual Analysis section."""

    allele_ids_x = [ids.split("_")[-1]
                    for ids in data_ind[locus_id]["locus_ids"]]

    allele_sizes_y = data_ind[locus_id]["allele_sizes"]

    dict_of_fig = dict(
        {
            "data": [
                {
                    "type": "scatter",
                    "x": allele_ids_x,
                    "y": allele_sizes_y,
                    "name": "Locus Details",
                    "mode": "markers",
                }
            ],
            "layout": {
                "xaxis": {
                    "title": {"text": "Allele ID"},
                },
                "yaxis": {
                    "title": {"text": "Sequence size in bp"},
                },
                "hovermode": "closest",
            },
        }
    )

    fig = go.Figure(dict_of_fig)

    return fig


def build_allele_length_histogram(pre_computed_data_path):
    """Builds the Allele Length Analysis plotly histogram Figure."""

    data_mode = get_data(pre_computed_data_path)

    data_alleles_mode = data_mode["mode"]

    allele_mode_x = [mo["alleles_mode"] for mo in data_alleles_mode]

    ln_y = [ln["locus_name"] for ln in data_alleles_mode]

    dict_of_fig = dict(
        {
            "data": [
                {
                    "type": "histogram",
                    "x": allele_mode_x,
                    "y": ln_y,
                }
            ],
            "layout": {
                "title": {"text": "Distribution of allele mode sizes"},
                "xaxis": {
                    "title": {"text": "Allele Mode Size"},
                    "showgrid": True,
                },
                "yaxis": {
                    "title": {"text": "Number of Loci"},
                },
            },
        }
    )

    fig = go.Figure(dict_of_fig)

    return fig


def build_allele_number_histogram(pre_computed_data_path):
    """Builds the Allele Numbers Analysis plotly histogram Figure."""

    data_total = get_data(pre_computed_data_path)

    data_total_all = data_total["total_alleles"]

    locus_name_y = [locus_name["locus_name"] for locus_name in data_total_all]

    nr_alleles_x = [nr["nr_alleles"] for nr in data_total_all]

    dict_of_fig2 = dict(
        {
            "data": [
                {
                    "type": "histogram",
                    "x": nr_alleles_x,
                    "y": locus_name_y,
                }
            ],
            "layout": {
                "title": {"text": "Number of Loci with given Number of Alleles"},
                "xaxis": {
                    "title": {"text": "Number of Different Alleles"},
                    "showgrid": True,
                },
                "yaxis": {
                    "title": {"text": "Number of Loci"},
                },
            },
        }
    )

    fig2 = go.Figure(dict_of_fig2)

    return fig2


def build_locus_stats_scatter(pre_computed_data_path):
    """Builds the Locus Statistics Plotly scatterplot Figure."""

    data_scatter = get_data(pre_computed_data_path)

    data_total_scatter = data_scatter["scatter_data"]

    locus_id = [lid["locus_id"] for lid in data_total_scatter]

    nr_alleles_scatter = [nrs["nr_alleles"] for nrs in data_total_scatter]

    scatter_data_median = [
        s_median["alleles_median"] for s_median in data_total_scatter
    ]

    scatter_data_min = [s_min["alleles_min"] for s_min in data_total_scatter]

    scatter_data_max = [s_max["alleles_max"] for s_max in data_total_scatter]

    dict_of_fig3 = dict(
        {
            "data": [
                {
                    "type": "scatter",
                    "x": scatter_data_min,
                    "y": nr_alleles_scatter,
                    "name": "Min",
                    "mode": "markers",
                    "marker": {
                        "opacity": 0.7,
                        "size": 4,
                    },
                    "hovertemplate": "<b>ID</b>: %{text}",
                    "text": locus_id,
                },
                {
                    "type": "scatter",
                    "x": scatter_data_max,
                    "y": nr_alleles_scatter,
                    "name": "Max",
                    "mode": "markers",
                    "marker": {
                        "opacity": 0.7,
                        "size": 4,
                    },
                    "hovertemplate": "<b>ID</b>: %{text}",
                    "text": locus_id,
                },
                {
                    "type": "scatter",
                    "x": scatter_data_median,
                    "y": nr_alleles_scatter,
                    "name": "Median",
                    "mode": "markers",
                    "marker": {
                        "opacity": 0.7,
                        "size": 4,
                    },
                    "hovertemplate": "<b>ID</b>: %{text}",
                    "text": locus_id,
                },
            ],
            "layout": {
                "title": {"text": "Locus Statistics"},
                "xaxis": {
                    "title": {"text": "Allele size in bp"},
                    "showgrid": True,
                    "zeroline": False,
                },
                "yaxis": {"title": {"text": "Number of alleles"}, "zeroline": False},
                "hovermode": "closest",
            },
        }
    )

    fig3 = go.Figure(dict_of_fig3)

    return fig3


def build_boxplot(pre_computed_data_path):
    """Builds the loci size distribution plotly Boxplot figure."""

    data_boxplot = get_boxplot_data(pre_computed_data_path)

    loci_x = data_boxplot["loci"]

    q1 = data_boxplot["q1"]
    median = data_boxplot["median"]
    q3 = data_boxplot["q3"]

    lowerfence = data_boxplot["min"]
    upperfence = data_boxplot["max"]

    dict_to_fig_box = dict(
        {
            "data": [
                {
                    "type": "box",
                    "x": loci_x,
                    "q1": q1,
                    "median": median,
                    "q3": q3,
                    "lowerfence": lowerfence,
                    "upperfence": upperfence,
                    "showlegend": False,
                }
            ],
            "layout": {
                "title": {"text": "Loci Size Variation"},
                "xaxis": {
                    "title": {"text": "Loci"},
                    "showticklabels": False,
                },
                "yaxis": {
                    "title": {"text": "Allele size variation"},
                },
                "boxgap": 0.05,
            },
        }
    )

    fig_box = go.Figure(dict_to_fig_box)

    return fig_box


def build_cds_hist(cds_dict):
    """Builds the CDS Analysis scatterplot plotly Figure."""

    mult3 = cds_dict["mult3"]
    stopC = cds_dict["stopC"]
    notStart = cds_dict["notStart"]
    cds = cds_dict["CDS_Alleles"]
    names = cds_dict["genes"]

    dict_of_fig = dict(
        {
            "data": [
                {
                    "type": "bar",
                    "x": cds,
                    "name": "CDS Alleles",
                    "orientation": "h",
                    "text": names,
                    "marker_color": "#045a8d",
                },
                {
                    "type": "bar",
                    "x": mult3,
                    "name": "Non multiple 3",
                    "orientation": "h",
                    "text": names,
                    "marker_color": "#006d2c",
                },
                {
                    "type": "bar",
                    "x": stopC,
                    "name": "> 1 stop codon",
                    "orientation": "h",
                    "text": names,
                    "marker_color": "#7b3294",
                },
                {
                    "type": "bar",
                    "x": notStart,
                    "name": "No Start/Stop codon",
                    "orientation": "h",
                    "text": names,
                    "marker_color": "#ec7014",
                },
            ],
            "layout": {
                "barmode": "stack",
                "bargap": 0.0,
                "bargroupgap": 0.0,
                "title": {"text": "Summary of problematic alleles per locus"},
                "xaxis": {
                    "title": {"text": "Number of occurrences"},
                    "gridcolor": "#eee",
                },
                "yaxis": {"showgrid": False, "showticklabels": False},
                # "paper_bgcolor": 'rgba(0,0,0,0)',
                "plot_bgcolor": "rgba(0,0,0,0)",
            },
        }
    )

    fig = go.Figure(dict_of_fig)

    return fig


# Main function where all components are assembled in a layout
def main(pre_computed_data_path, schema_dir, light_mode):

    # Instatiate the app, loading the external stylesheet
    external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]

    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    # Initial Explanation text
    markdown_initial_text = dcc.Markdown(
        """
        # Schema Evaluator

        Provides summary charts that allows users to explore:
        - The diversity (number of alleles) at each locus (**Panel A**);
        - The variation of allele mode sizes per locus (**Panel B**);
        - Summary statistics (minimum allele size in blue, minimum allele size in orange and median allele size in green) for each locus (**Panel C**);
        - The loci size distribution (**Panel D**);
        - The presence of alleles that are not CDSs (when evaluating schemas called by other algorithms) (**Panel E**).

        After the panels, a **dropdown menu** allows users to select an **individual locus** to analysed.
        By selecting a locus the following will be displayed:
        - 2 charts (histogram and scatterplot) containing an **analysis of the allele sizes**;
        - a table with **summary statistics** of the alleles;
        - a **multiple sequence alignment** of the alleles.
        """
    )

    # Build the plotly Figure object for panel A - Allele Numbers Analysis.
    allele_number_fig = build_allele_number_histogram(pre_computed_data_path)

    dash_allele_number_fig = dcc.Graph(
        id="allele-number", figure=allele_number_fig)

    # Build the plotly Figure object for panel B - Allele Length Analysis.
    allele_length_fig = build_allele_length_histogram(pre_computed_data_path)

    dash_allele_length_fig = dcc.Graph(
        id="allele-length", figure=allele_length_fig)

    # Build the plotly Figure object for panel C - Locus Statistics.
    locus_statictics_fig = build_locus_stats_scatter(pre_computed_data_path)

    dash_locus_statictics_fig = dcc.Graph(
        id="locus-statistics", figure=locus_statictics_fig
    )

    # Build the plotly Figure object for panel D - Loci Size Distribution.
    boxplot = build_boxplot(pre_computed_data_path)

    dash_boxplot = dcc.Graph(id="boxplot", figure=boxplot)
    # Build the plotly Figure object for panel E - CDS Analysis.

    # Get the CDS Analysis dataframe.
    cds_analysis_df, histogram_cds_data = get_cds_df(schema_dir)

    # Create Dash's DataTable component.
    table_df = dash_table.DataTable(
        id="table-cds",
        columns=[{"name": i, "id": i} for i in cds_analysis_df.columns],
        data=cds_analysis_df.to_dict("records"),
        page_size=20,
        style_table={"height": "300px", "overflowY": "auto"},
    )

    # Build the plotly Figure object for the CDS Analysis' histogram.
    cds_bar_fig = build_cds_hist(histogram_cds_data)

    dash_cds_bar_fig = dcc.Graph(id="cds-bar-fig", figure=cds_bar_fig)

    # Locus Individual analysis
    options, data_ind = build_locus_options(pre_computed_data_path)

    # Create Dash's Dropdown component.
    dropdown = dcc.Dropdown(
        id="locus-ind-dropdown",
        options=options,
        placeholder="Please select a locus to analyse.",
    )

    # MSA
    # with open(
    #     "/home/pcerqueira/github_repos/SCHEMA_EVALUATOR_REFACTOR/GCF-000007265-protein1_aligned_prot.fasta",
    #     "r",
    # ) as m:
    #     msa_file = m.read()

    # msa_viewer = html.Div(
    #     [
    #         html.Div(
    #             dashbio.AlignmentChart(
    #                 id="my-alignment-viewer",
    #                 data=msa_file,
    #                 extension="fasta",
    #                 height=2000,
    #                 groupbars=False,
    #                 overview="heatmap",
    #             )
    #         ),
    #         html.Div(id="alignment-viewer-output"),
    #     ]
    # )

    # App layout
    app.layout = html.Div(
        id="main-div",
        className="app-body",
        children=[
            html.Div(
                children=[markdown_initial_text],
                style=dict(display="flex", justifyContent="center"),
            ),
            html.Br(),
            html.Div(
                [
                    html.Div(
                        [html.H3("A - Allele Numbers Analysis"),
                         dash_allele_number_fig],
                        className="six columns",
                    ),
                    html.Div(
                        [html.H3("B - Allele Length Analysis"),
                         dash_allele_length_fig],
                        className="six columns",
                    ),
                ],
                className="row",
            ),
            html.Div(
                [
                    html.Div(
                        [html.H3("C - Locus Statistics"),
                         dash_locus_statictics_fig],
                        className="six columns",
                    ),
                    html.Div(
                        [html.H3("D - Locus Size Variation"), dash_boxplot],
                        className="six columns",
                    ),
                ],
                className="row",
            ),
            html.Br(),
            html.H3("E - CDS Analysis"),
            html.Div(table_df),
            html.Br(),
            html.Div(dash_cds_bar_fig),
            html.Br(),
            html.H3("Locus Individual Analysis"),
            dropdown,
            html.Div(
                [
                    html.Div([dcc.Graph(id="locus-hist")],
                             className="six columns"),
                    html.Div([dcc.Graph(id="locus-scatter")],
                             className="six columns"),
                ],
                className="row",
            ),
            html.Div(
                [
                    dash_table.DataTable(
                        id="table-locus-ind",
                        page_size=20,
                        style_table={"height": "300px", "overflowY": "auto"},
                    )
                ]
            ),
            # html.Br(),
            html.H3("Multiple Sequence alignment of an Individual Locus"),
            # msa_viewer,
        ]
    )

    # App Callback functions.

    # Dropdown callback
    @app.callback(
        [
            dash.dependencies.Output("locus-hist", "figure"),
            dash.dependencies.Output("locus-scatter", "figure"),
            dash.dependencies.Output("table-locus-ind", "columns"),
            dash.dependencies.Output("table-locus-ind", "data"),
        ],
        [dash.dependencies.Input("locus-ind-dropdown", "value")],
    )
    def update_output(value):
        """Update dropdown output."""
        if value is None:
            raise PreventUpdate

        locus_ind_hist = build_locus_histogram(data_ind, value)

        locus_ind_scatter = build_locus_scatter(data_ind, value)

        locus_ind_df = get_locus_ind_df(data_ind, value)

        locus_ind_table_df_columns = [
            {"name": i, "id": i} for i in locus_ind_df.columns]

        locus_ind_table_df_data = locus_ind_df.to_dict("records")

        return (
            locus_ind_hist,
            locus_ind_scatter,
            locus_ind_table_df_columns,
            locus_ind_table_df_data,
        )

    # MSA Callback

    @app.callback(
        dash.dependencies.Output("alignment-viewer-output", "children"),
        [dash.dependencies.Input("my-alignment-viewer", "eventDatum")],
    )
    def update_output_msa(value):
        if value is None:
            return "No data."
        return str(value)

    app.run_server(debug=True)


# if __name__ == "__main__":
#     main()
#     # app.run_server(debug=True)
