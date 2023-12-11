import dash
from dash import Dash, html, dcc, callback, Output, Input, dash_table
import plotly.express as px
from dash.exceptions import PreventUpdate
import pandas as pd
import montu
import numpy as np

# Preliminary data
columns = [{'name': col, 'id': col} for col in ['MN', 'HD', 'Name', 'Bayer',
                                                'RAJ2000', 'DecJ2000', 'Constellation',
                                                'Vmag', 'Distance']]

# Estilos para la tabla
style_cell = {
    'backgroundColor': 'white',
    'color': 'black',
    'fontSize': 14,
    'font-family': 'Arial, sans-serif',
    'textAlign': 'center',
    'padding': '10px',
    'border': '1px solid lightgrey'
}

style_header = {
    'backgroundColor': 'lightgrey',
    'fontWeight': 'bold',
    'border': '1px solid black'
}

style_data_conditional = [
    {
        'if': {'row_index': 'odd'},
        'backgroundColor': 'rgb(220, 220, 220)',
    }
]

style_table = {
    'maxWidth': '100%',
    'overflowX': 'auto',
    'margin': 'auto',
    'border': '1px solid black',
    'boxShadow': '4px 4px 4px lightgrey'
}

# Layout
dash.register_page(__name__)  # Uncomment in production
layout = html.Div([
    html.H3(children='Naked eye stars ordered by brightness', 
            style={'textAlign': 'center'}),
    dcc.Input(id='dummy-input', type='hidden', value=0),
    dcc.Loading(
        id="table-loading",
        children=[dash_table.DataTable(
            id='stars-table', 
            data=[], 
            columns=columns, 
            page_size=10, 
            style_cell=style_cell,
            style_header=style_header,
            style_data_conditional=style_data_conditional,
            style_table=style_table
        )],
        type="default", 
        fullscreen=False,
    ),
], style={'backgroundColor': '#f5e2a1'})

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
