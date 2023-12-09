import dash
import montu
from dash import Dash, html, dcc, callback, Output, Input,dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import numpy as np

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
- **delta-t**: Time since last quarter in days.
- **delta-0**: Time since initial date in days.
- **Caniucular**: Date in *caniucular* calendar, ie. civil egyptian calendar. In this calendar the datum, namely *horus* year 0 or `hrw 0-I-Akhet-1`, corresponds to bce 2782-07-20 00:00:00.
"""

################################################################
# Layout
################################################################
dash.register_page(__name__) # Uncomment in production
layout = html.Div([
    html.H3(children=f'Lunar phases', style={'textAlign':'center'}),
    html.Div([
        dcc.Markdown(module_quickstart_doc),
    ], style={'padding':'1%'}),
    html.Div([
        "Initial date (format [-]CCYY-MM-DD):",
        dcc.Input(id='since', placeholder='Right now',
                  value='-236-08-23', type='text', style={
                      'margin-left':'1%', 'border-radius': '15px', 'border': '2px solid gold'
                  }),
        dcc.RadioItems(id="calendar",
                       options=[
                           dict(label='Proleptic Gregorian', value='proleptic'),
                           dict(label='Mixed Calendar', value='mixed'),
                       ],
                       value='mixed',
                       inline=True,
                       style={'margin-left':'5%'}),
    ], style={'padding':'1%'}),
    html.Div([
        "Starting at quarter:",
        dcc.Dropdown(['new', 'first', 'full', 'last', '--'], value='--', id='starting_at', style={
            'width': '10em', 'border-radius': '15px', 'border': '2px solid gold'
        }),
    ], style={'padding':'1%'}),
    html.Div([
        "Number of synodic months:",
        dcc.Input(id='nummonths', value='12', type='number', style={
            'margin-left':'1%', 'width':'3em', 'border-radius': '15px', 'border': '2px solid gold'
        }),
    ], style={'padding':'1%'}),
    html.Div([
        "Options:",
        dcc.Checklist(
            id='options', 
            options=[
                {'label': 'Show caniucular', 'value': 'show_caniucular'},
            ],
            value=['show_caniucular'],
            inline=True
        ),
    ], style={'padding':'1%'}),
    dcc.Loading(
        id="table-loading",
        children=[
            dash_table.DataTable(
                id='table', data=[], page_size=10,
                style_cell={
                    'padding-right': '30px',
                    'padding-left': '30px',
                    'text-align': 'center',
                    'marginLeft': 'auto',
                    'marginRight': 'auto'
                }
            )],
        type="default", fullscreen=False,
    ),
    html.Hr(),
    html.Div([
        html.H3("Field explanation"),
        dcc.Markdown(module_field_doc),
    ], style={'padding':'1%'}),
], style={'backgroundColor': '#f5e2a1'})


################################################################
# Routines
################################################################

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
def update_table(since,calendar,starting_at,nummonths,options):

    # Set time
    if since == '':
        mtime = montu.Time()
    else:
        mtime = montu.Time(since,calendar=calendar)
    
    # Set input options
    starting_at = None if starting_at == '--' else starting_at
    output = 'datepro' if calendar=='proleptic' else 'datemix'
    
    # Get lunar phases
    quarter_dates = montu.Moon.next_moon_quarters(since=mtime,
                                                  starting_at=starting_at,
                                                  numquarters=4*int(nummonths),
                                                  output=output,
                                                  format='columns')

    # Show caniucular
    if 'show_caniucular' in options:
        for date in quarter_dates:
            mtime = montu.Time(date['Datetime'],calendar=calendar)
            date['Caniucular'] = mtime.readable.datecan

    # Create table
    return quarter_dates