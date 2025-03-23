import dash
from dash import Dash, html, dcc, callback, Output, Input, dash_table
import plotly.express as px
from dash.exceptions import PreventUpdate
import pandas as pd
import montu
import numpy as np
import dash_bootstrap_components as dbc
from utils.theme import egyptian_palette

# Preliminary data
columns = [{'name': col, 'id': col} for col in ['MN', 'HD', 'Name', 'Bayer',
                                                'RAJ2000', 'DecJ2000', 'Constellation',
                                                'Vmag', 'Distance']]

# Estilos para la tabla
style_cell = {
    'backgroundColor': 'white',
    'color': egyptian_palette["text"],
    'fontSize': 14,
    'font-family': 'Arial, sans-serif',
    'textAlign': 'center',
    'padding': '10px',
    'border': '1px solid ' + egyptian_palette["header"]
}

style_header = {
    'backgroundColor': egyptian_palette["header"],
    'color': 'white',
    'fontWeight': 'bold',
    'border': '1px solid ' + egyptian_palette["header"]
}

style_data_conditional = [
    {
        'if': {'row_index': 'odd'},
        'backgroundColor': '#faf0d7',
    }
]

style_table = {
    'maxWidth': '100%',
    'overflowX': 'auto',
    'margin': 'auto',
    'border': '2px solid ' + egyptian_palette["header"],
    'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)',
    'borderRadius': '8px',
    'marginTop': '20px',
    'marginBottom': '20px'
}

# Layout
dash.register_page(__name__)
layout = html.Div([
    dbc.Card([
        dbc.CardHeader(
            html.H3(
                "Naked Eye Stars",
                className="text-center",
                style={"color": egyptian_palette["text"]},
            )
        ),
        dbc.CardBody([
            dcc.Input(id='dummy-input', type='hidden', value=0),
            dcc.Loading(
                id="table-loading",
                children=[
                    dash_table.DataTable(
                        id='stars-table',
                        data=[],
                        columns=columns,
                        page_size=10,
                        style_cell=style_cell,
                        style_header=style_header,
                        style_data_conditional=style_data_conditional,
                        style_table=style_table
                    )
                ],
                type="default",
                fullscreen=False,
            ),
        ])
    ], style={"backgroundColor": egyptian_palette["background"]}),
], className="p-4")

# Callback
class stars(object):
    stars_visible = None

@callback(
    Output('stars-table', 'data'),
    Input('dummy-input', 'value'),
)
def update_stars(dummy):
    global stars_visible
    stars_visible = montu.Stars(filename=montu.Util._data_path('montu_stellar_catalogue_v37_visible.csv'))
    stars_visible.data['RAJ2000'] = stars_visible.data.apply(lambda row: montu.D2H(row['RAJ2000']), axis=1)
    stars_visible.data['DecJ2000'] = stars_visible.data.apply(lambda row: montu.D2H(row['DecJ2000']), axis=1)
    return stars_visible.data.to_dict('records')
