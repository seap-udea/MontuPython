import dash
from dash import Dash, html, dcc, callback, Output, Input
from dash.exceptions import PreventUpdate
import plotly.express as px
import pandas as pd
import montu
import numpy as np
import dash_bootstrap_components as dbc
from utils.theme import egyptian_palette

################################################################
# Preliminary data
################################################################
# Planets
planets = [montu.Planet(value) for value in montu.PLANETARY_NAMES.values() if value not in ['SUN','MOON','EARTH']]
planet_names = [planet.name for planet in planets]

# Condition columns
properties = ['RAJ2000', 'DecJ2000', 'RAEpoch', 'DecEpoch',
       'RAGeo', 'DecGeo', 'el', 'az', 'ha', 'Vmag', 'rise_time', 'rise_az',
       'set_time', 'set_az', 'transit_time', 'transit_el', 'elongation',
       'earth_distance', 'sun_distance', 'is_circumpolar', 'is_neverup',
       'angsize', 'phase', 'hlat', 'hlon', 'hlong', 'datestr']

################################################################
# Layout
################################################################
dash.register_page(__name__)

layout = html.Div([
    dbc.Card([
        dbc.CardHeader(
            html.H3(
                "Planetary Ephemerides",
                className="text-center",
                style={"color": egyptian_palette["text"]},
            )
        ),
        dbc.CardBody([
            dbc.Row([
                dbc.Col([
                    dbc.Label("Initial date (format [-]CCYY-MM-DD):", className="me-2"),
                    dbc.Input(
                        id='initial-date',
                        value='-1500-01-01',
                        type='text',
                        className="mb-3",
                        style={'border': f'1px solid {egyptian_palette["accent"]}'}
                    ),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Time span (in years):", className="me-2"),
                    dbc.Input(
                        id='time-span',
                        value='10',
                        type='text',
                        className="mb-3",
                        style={'border': f'1px solid {egyptian_palette["accent"]}'}
                    ),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Number of points:", className="me-2"),
                    dbc.Input(
                        id='number-points',
                        value='120',
                        type='text',
                        className="mb-3",
                        style={'border': f'1px solid {egyptian_palette["accent"]}'}
                    ),
                ], width=12),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Planet:", className="me-2"),
                    dcc.Dropdown(
                        id='dropdown-planet',
                        options=[{'label': name, 'value': name} for name in planet_names],
                        value='Mercury',
                        multi=True,
                        style={
                            'width': '100%',
                        }
                    ),
                ], width=12, className="mb-3"),
            ]),

            dbc.Row([
                dbc.Col([
                    dbc.Label("Property:", className="me-2"),
                    dbc.Select(
                        id='dropdown-property',
                        options=[{'label': prop, 'value': prop} for prop in properties],
                        value='DecEpoch',
                        style={
                            'width': '100%',
                            'border': f'1px solid {egyptian_palette["accent"]}'
                        }
                    ),
                ], width=12, className="mb-3"),
            ]),

            dbc.Row([
                dbc.Col([
                    dcc.Loading(
                        id="graph-loading",
                        children=[dcc.Graph(id='graph-content')],
                        type="default",
                        fullscreen=False,
                    ),
                ], width=12),
            ]),
        ])
    ], style={"backgroundColor": egyptian_palette["background"]}),
], className="p-4")

################################################################
# Compute ephemerides
################################################################
def compute_ephemerides(initial='-2500-01-01',timespan=1,numpoints=52):
    global planetary_ephemerides
    # Time and place
    mtime = montu.Time(initial)
    Tebas = montu.Observer(lon=33, lat=24)

    # Range times
    mts = []
    dates = []
    for dt in np.linspace(0,float(timespan or 1)*montu.YEAR,int(numpoints or 10)):
        mt = (mtime+dt).get_readable()
        mts += [mt]
        dates += [f'{mt.readable.year}-{mt.readable.month}-{mt.readable.day}']

    # Compute planetary positions
    planetary_ephemerides = pd.DataFrame()
    for planet in planets:
        planet.reset_store()
        for mt in mts:
            planet.conditions_in_sky(at=mt,observer=Tebas,store=True)
        planet.tabulate_ephemerides()
        planet.ephemerides['datestr'] = dates
        planetary_ephemerides = pd.concat([planetary_ephemerides,planet.ephemerides],ignore_index=True)
    return planetary_ephemerides

################################################################
# Change plot
################################################################
@callback(
    Output('graph-content', 'figure'),
    Input('dropdown-planet', 'value'),
    Input('dropdown-property', 'value'),
    Input('initial-date', 'value'),
    Input('time-span', 'value'),
    Input('number-points', 'value'),
)
def update_graph(planets=['Mercury'],property='DecEpoch',
                 initial='-2500-01-01',timespan=1,numpoints=100):

    planetary_ephemerides = compute_ephemerides(initial,timespan,numpoints)

    if isinstance(planets,str):
        planets=[planets]

    if planetary_ephemerides is not None:
        mtime_initial = montu.Time(initial)
        mtime_final = (mtime_initial + float(timespan)*montu.YEAR).get_readable()
        mask = planetary_ephemerides.Name.isin(planets)
        fig = px.line(planetary_ephemerides[mask], x='datestr', y=property, color='Name')
        fig.update_layout(
            title=dict(text=f"Property '{property}' of {', '.join(planets)} starting at {mtime_initial.strftime('%Y-%m-%d')} until {mtime_final.strftime('%Y-%m-%d')}"),
            xaxis_title="Date [Month & Year]",
            xaxis=dict(
                rangeslider=dict(
                    visible=True
                ),
            )
        )
    else:
        raise PreventUpdate
    return fig