import dash
from dash import Dash, html, dcc, callback, Output, Input, dash_table
import plotly.express as px
from dash.exceptions import PreventUpdate
import pandas as pd
import montu
import numpy as np

################################################################
# Preliminary data
################################################################
columns=[{'name': col, 'id': col} for col in ['MN','HD','Name','Bayer',
                                              'RAJ2000','DecJ2000','Constellation',
                                              'Vmag','Distance'
                                              ]]

################################################################
# Layout
################################################################
#dash.register_page(__name__) # Uncomment in production
layout = html.Div([
    html.H3(children=f'Naked eye stars ordered by brightness', 
            style={'textAlign':'center'}),
    dcc.Input(id='dummy-input',type='hidden',value=0),
    dcc.Loading(
        id="table-loading",
        children=[dash_table.DataTable(id='stars-table',data=[],columns=columns,page_size=10)],
        type="default",fullscreen=False,
    ),
])

################################################################
# Callback
################################################################
class stars(object):
    stars_visible = None

@callback(
    Output('stars-table', 'data'),
    Input('dummy-input', 'value'),
)
def update_stars(dummy):
    global stars_visible
    stars_visible = montu.Stars(filename=montu.Util._data_path('montu_stellar_catalogue_v37_visible.csv'))
    stars_visible.data['RAJ2000'] = stars_visible.data.apply(lambda row:montu.D2H(row['RAJ2000']),axis=1)
    stars_visible.data['DecJ2000'] = stars_visible.data.apply(lambda row:montu.D2H(row['DecJ2000']),axis=1)
    return stars_visible.data.to_dict('records')

################################################################
# Independent app code
################################################################
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = layout
if __name__ == '__main__':
    app.run(debug=True,port=8002)
