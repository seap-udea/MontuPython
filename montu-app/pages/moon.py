import dash
import montu
from dash import Dash, html, dcc, callback, Output, Input, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import numpy as np
from utils.theme import egyptian_palette

################################################################
# Preliminary data
################################################################
module_quickstart_doc = """
This module allows you to calculate the date of the lunar quarters (lunar phases)
starting at an arbitrary date (input `Initial date`). If you leave this filed blank the app
will assume that the **initial date is the present date and time**. You may select which is the first 
quarter to be detected since the initial date (dropdown menu `Starting at quarter`). If
you select '--' the app will select the first quarter just after the date.
"""

module_field_doc = """
- **Quarter**: Type of quarter.
- **Datetime**: Date and time of quarter in UTC.
- **Delta_t ($\Delta t$)**: Time since last quarter in days ($\Delta t$).
- **Delta_0 ($\Delta t_0$)**: Time since initial date in days.
- **Caniucular**: Date in *caniucular* calendar, ie. civil egyptian calendar. In this calendar the datum, namely *horus* year 0 or `hrw 0-I-Akhet-1`, corresponds to bce 2782-07-20 00:00:00.
"""

################################################################
# Layout
################################################################
dash.register_page(__name__)

## Style for the table
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
        'backgroundColor': '#faf0d7',  # Un tono más claro del background
    }
]

style_table = {
    'maxWidth': '100%',
    'overflowX': 'auto',
    'margin': 'auto',
    'border': '1px solid ' + egyptian_palette["header"],
    'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)',
    'borderRadius': '8px',
    'marginTop': '20px',
    'marginBottom': '20px'
}

layout = html.Div([
    dbc.Card([
        dbc.CardHeader(
            html.H3(
                "Lunar Phases",
                className="text-center",
                style={"color": egyptian_palette["text"]},
            )
        ),
        dbc.CardBody([
            dcc.Markdown(module_quickstart_doc),
            
            dbc.Row([
                dbc.Col([
                    dbc.Label("Initial date (format [-]CCYY-MM-DD):", className="me-2"),
                    dbc.Input(
                        id='since',
                        placeholder='Right now',
                        value='-236-08-23',
                        type='text',
                        className="mb-3",
                        style={'border': f'1px solid {egyptian_palette["accent"]}'}
                    ),
                    
                    dbc.RadioItems(
                        id="calendar",
                        options=[
                            dict(label='Proleptic Gregorian', value='proleptic'),
                            dict(label='Mixed Calendar', value='mixed'),
                        ],
                        value='mixed',
                        inline=True,
                        className="mb-3"
                    ),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Starting at quarter:", className="me-2"),
                    dbc.Select(
                        id='starting_at',
                        options=[
                            {'label': q, 'value': q} for q in ['new', 'first', 'full', 'last', '--']
                        ],
                        value='--',
                        style={'width': '200px', 'border': f'1px solid {egyptian_palette["accent"]}'}
                    ),
                ], width=12, className="mb-3"),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Number of synodic months:", className="me-2"),
                    dbc.Input(
                        id='nummonths',
                        value='12',
                        type='number',
                        style={'width': '100px', 'border': f'1px solid {egyptian_palette["accent"]}'}
                    ),
                ], width=12, className="mb-3"),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Options:", className="me-2"),
                    dbc.Checklist(
                        id='options',
                        options=[
                            {'label': 'Show caniucular', 'value': 'show_caniucular'},
                        ],
                        value=['show_caniucular'],
                        inline=True
                    ),
                ], width=12, className="mb-3"),
            ]),

            dbc.Row([
                dbc.Col([
                    dcc.Loading(
                        id="table-loading",
                        children=[
                            dash_table.DataTable(
                                id='table',
                                data=[],
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
                ], width=12),
            ]),

            html.Hr(),

            dbc.Row([
                dbc.Col([
                    html.H3("Field explanation"),
                    dcc.Markdown(module_field_doc, mathjax=True),
                ], width=12),
            ]),
        ])
    ], style={"backgroundColor": egyptian_palette["background"]}),
], className="p-4")

################################################################
# Callbacks
################################################################
@callback(
    Output('table', 'data'),
    Input('since', 'value'),
    Input('calendar', 'value'),
    Input('starting_at', 'value'),
    Input('nummonths', 'value'),
    Input('options', 'value'),
)
def update_table(since, calendar, starting_at, nummonths, options):
    # Set time
    if since == '':
        mtime = montu.Time()
    else:
        mtime = montu.Time(since, calendar=calendar)
    
    # Set input options
    starting_at = None if starting_at == '--' else starting_at
    output = 'datepro' if calendar=='proleptic' else 'datemix'
    
    # Get lunar phases
    quarter_dates = montu.Moon.next_moon_quarters(
        since=mtime,
        starting_at=starting_at,
        numquarters=4*int(nummonths),
        output=output,
        format='columns'
    )

    # Show caniucular
    if 'show_caniucular' in options:
        for date in quarter_dates:
            mtime = montu.Time(date['Datetime'], calendar=calendar)
            date['Caniucular'] = mtime.readable.datecan

    return quarter_dates