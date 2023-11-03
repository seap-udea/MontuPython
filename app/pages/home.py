import dash 
from dash import dcc, html

dash.register_page(__name__, path='/')

layout = html.Div([
    html.H5('Atronomical Ephemerides for the ancient world', className='h5'), 
    html.Img(src='assets/Montu.png', className='center', width='200px')
])