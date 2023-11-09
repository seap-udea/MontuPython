import dash
from dash import Dash, html, dcc, callback, Output, Input, dash_table
import plotly.express as px
import pandas as pd
import montu
import numpy as np

################################################################
# Preliminary data
################################################################
stars_visible = montu.Stars(filename=montu.Util._data_path('montu_stellar_catalogue_v37_visible.csv'))
stars_visible.data['RAJ2000'] = stars_visible.data.apply(lambda row:montu.D2H(row['RAJ2000']),axis=1)
stars_visible.data['DecJ2000'] = stars_visible.data.apply(lambda row:montu.D2H(row['DecJ2000']),axis=1)
columns=[{'name': col, 'id': col} for col in ['MN','HD','Name','Bayer',
                                              'RAJ2000','DecJ2000','Constellation',
                                              'Vmag','Distance'
                                              ]]

################################################################
# Layout
################################################################
dash.register_page(__name__) # Uncomment in production
layout = html.Div([
    html.H3(children=f'{len(stars_visible.data)} naked eye stars ordered by brightness', 
            style={'textAlign':'center'}),
    dash_table.DataTable(data=stars_visible.data.to_dict('records'),columns=columns,page_size=10),
    dcc.Graph(figure=px.histogram(stars_visible.data, x='Constellation', y='Distance', histfunc='count'))
    ])