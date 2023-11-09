import dash
import montu
from dash import Dash, html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import numpy as np

################################################################
# Preliminary data
################################################################
select_input = html.Div([
        dbc.RadioItems(
            id="input-type",
            className="btn-group",
            inputClassName="btn-check",
            labelClassName="btn btn-outline-primary",
            labelCheckedClassName="active",
            options=[
                {"label": "Calendar", "value": "calendar"},
                {"label": "Julian day", "value": "julian"},
            ],
            value="calendar",
        ),
        html.Br(),
        html.Div(id="output")
    ],className="radio-group",
    )

input_date = html.Div(
    [
        html.Div([
            html.H5('Date',{'display':'inline-block', 'border': '1px solid black'})
        ],style={'display':'inline-block', 'margin-right':20, 'margin-left':20, 'margin-top':10}),
        
        dcc.Input(id="year", type="number", value = -44, placeholder="Year", 
                  style={'marginRight':'10px', 'width':'10%', 'marginTop':'10px'}),
        dcc.Input(id="month", type="number", min=1, max=12, value = 1, placeholder="Month", 
                  style={'marginRight':'10px', 'width':'10%','marginTop':'10px'}),
        dcc.Input(id="day", type="number", placeholder="Day", value = 1, min=1, max=31, step=1, 
                  style={'marginTop':'10px'}),
        
        html.Div([
        html.H5('Time',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
        dcc.Input(id="hour", type="number", value = 0, min=0, max=24, placeholder="Hour", style={'width':'7%'}),
        
        html.Div([
        html.H6(':',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':5, 'margin-left':5}),
        dcc.Input(id="min", type="number", value = 0, min=0, max=59, debounce=True, placeholder="Min", style={'width':'7%'}),
        
        html.Div([
        html.H6(':',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':5, 'margin-left':5}),
        dcc.Input(id="sec", type="number", value = 0, min=0, max=59, placeholder="Sec", step=1, style={'width':'7%'}),
        dcc.Input(id="julian-input", type="hidden",value=0),
    ]
)

julian_date = html.Div([
    html.Div([
    html.H5('Julian Date',{'display':'inline-block', 'border': '1px solid black'})], 
    style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
    dcc.Input(id="julian-input", type="number", value=2451545.0, placeholder="Julian Day", 
              style={'width':'10%', 'marginTop':'10px'}),
])

################################################################
# Layout
################################################################
dash.register_page(__name__) # Uncomment in production
layout = html.Div([
    html.H3(children=f'Date converter', style={'textAlign':'center'}),
    html.Center(html.Div([
        dbc.RadioItems(
            id="input-type",
            className="btn-group",
            inputClassName="btn-check",
            labelClassName="btn btn-outline-primary",
            labelCheckedClassName="active",
            options=[
                {"label": "Input calendar", "value": "calendar"},
                #{"label": "Input julian day", "value": "julian"},
            ],
            value="calendar",
        ),
        html.Div(id="output")
    ],className="radio-group",
    )),
    html.Center(html.Div([
        dbc.Label("Calendar: "),
        dbc.RadioItems(
            options=[
                {"label": "Gregorian proleptic", "value": "proleptic"},
                {"label": "Mixed", "value": "mixed"}
            ],
            value="proleptic",
            id="calendar-radio-input",
            inline=True
        ),
        html.Br()
    ])), 
    html.Hr(),
    dcc.Markdown(id="date-output",style={'white-space': 'pre-line', 'font-size': '1.5em'}),
    html.Hr(),
    dcc.Input(id="year", type="hidden",value=2000),
    dcc.Input(id="month", type="hidden",value=2000),
    dcc.Input(id="day", type="hidden",value=2000),
    dcc.Input(id="hour", type="hidden",value=2000),
    dcc.Input(id="min", type="hidden",value=2000),
    dcc.Input(id="sec", type="hidden",value=2000),
    ])

################################################################
# Routines
################################################################

################################################################
# Callbacks
################################################################
@callback(
    Output("date-output","children"),
    Input("year", "value"),
    Input("month", "value"),
    Input("day", "value"),
    Input("hour", "value"),
    Input("min", "value"),
    Input("sec", "value"),
    Input("calendar-radio-input", "value"),
)
def convert_date_calendar(year, month, day, hour, min, sec, calendar):
    date = f'{int(year or 1)}-{int(month or 1):02d}-{int(day or 1):02d} {int(hour or 0):02d}:{int(min or 0):02d}:{int(sec or 0):02d}'  
    mtime = montu.Time(date, format='iso',scale='utc',calendar=calendar)    
    return get_date_output(mtime)

def get_date_output(mtime):
    # Extract variables
    spice = mtime.readable.datespice # Date in gregorian proleptic
    proleptic = mtime.readable.datepro # Date in gregorian proleptic
    mixed = mtime.readable.datemix # Date in mixed
    jd_utc = mtime.jed # Date in Julian Day (utc scale)
    jd = mtime.tt # Date in ephemerides time (tt scale)
    es = mtime.et # Date in ephemerides time (utc scale)
    delta_t = mtime.deltat

    # Express result
    result = f'''
- **Date string**:
    - **Gregorian proleptic** (human readable) : {spice} 
    - **Gregorian Proleptic** (astronomical convention): {proleptic} 
    - **Mixed** (gregorian or Julian) : {mixed} 
- **Julian day (UTC)** : {jd_utc} days
- **Julian day (TT)**: {jd} days
- **Ephemerides seconds (TT)**: {es} seconds
- **Delta-t**: {delta_t} seconds
'''
    return result

#///////////////////////////////
# Change input
#///////////////////////////////
@callback(Output("output", "children"), Input("input-type", "value"))
def display_value(value):
    return input_date 