import dash
import montu
from dash import Dash, html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
import numpy as np

# Paleta de colores egipcios
egyptian_palette = {
    'background': '#f5e2a1',  # Fondo amarillo claro
    'text': '#000000',  # Texto en negro
    'header': '#cda434',  # Encabezado en dorado egipcio
    'accent': '#d97824',  # Acento en naranja oscuro
}

dash.register_page(__name__)

# Estilos Egipcios
egyptian_style = {
    'backgroundColor': egyptian_palette['background'],  # Fondo amarillo claro
    'padding': '10px',
    'borderRadius': '5px',
    'margin': '10px 0',  # Añadido para espaciar entre las secciones
}

input_style = {
    'width': '10%',
    'marginRight': '20px',  # Aumentado para más espacio
    'marginTop': '10px',
    'marginBottom': '10px',  # Añadido para espaciar entre los inputs
    'border': '1px solid #cda434',  # Borde dorado
    'padding': '5px',
    'borderRadius': '5px'
}

label_style = {
    'marginRight': '20px',
    'marginLeft': '20px',
    'marginTop': '10px'
}

# Preliminary data
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
    ], className="radio-group",
)

input_date = html.Div(
    [
        html.Div([html.H5('Date', style=label_style)], style={'display': 'inline-block'}),
        
        dcc.Input(id="year", type="number", value=-1321, placeholder="Year", style=input_style),
        dcc.Input(id="month", type="number", min=1, max=12, value=7, placeholder="Month", style=input_style),
        dcc.Input(id="day", type="number", placeholder="Day", value=20, min=1, max=31, step=1, style=input_style),
        
        html.Div([html.H5('Time', style=label_style)], style={'display': 'inline-block'}),
        dcc.Input(id="hour", type="number", value=0, min=0, max=24, placeholder="Hour", style=input_style),
        dcc.Input(id="min", type="number", value=0, min=0, max=59, debounce=True, placeholder="Min", style=input_style),
        dcc.Input(id="sec", type="number", value=0, min=0, max=59, placeholder="Sec", step=1, style=input_style),
    ], style={'display': 'flex', 'justifyContent': 'center', 'flexWrap': 'wrap'}
)

# Layout
layout = html.Div(style={'background-color': egyptian_palette['background']}, children=[
    html.H3(children='Date converter', style={'textAlign': 'center', 'color': egyptian_palette['text']}),
    html.Center(select_input),
    html.Center(html.Div([
        dbc.Label("Calendar: ", style=label_style),
        dbc.RadioItems(
            options=[
                {"label": "Gregorian proleptic", "value": "proleptic"},
                {"label": "Mixed", "value": "mixed"}
            ],
            value="mixed",
            id="calendar-radio-input",
            inline=True,
            style=egyptian_style
        ),
    ], style={'textAlign': 'center'})),  # Estilo modificado para centralizar
    html.Hr(),
    dcc.Markdown(id="date-output", style={'white-space': 'pre-line', 'font-size': '1.5em', 'color': egyptian_palette['text']}),
    html.Hr(),
    dcc.Input(id="year", type="hidden", value=2000),
    dcc.Input(id="month", type="hidden", value=1),
    dcc.Input(id="day", type="hidden", value=1),
    dcc.Input(id="hour", type="hidden", value=0),
    dcc.Input(id="min", type="hidden", value=0),
    dcc.Input(id="sec", type="hidden", value=0),
])


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
    caniucular = mtime.readable.datecan # Date in caniucular (civil egyptian)
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
    - **Caniucular** (civil egyptian) : {caniucular} 
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
