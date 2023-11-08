import dash 
from dash import dcc, html

dash.register_page(__name__, path='/')

layout = html.Div(children=[
        dcc.Markdown("""
`MontuPython` (transileterated mnṯw ꜥꜣpp(y)) is a Python package intended to compute astronomical ephemerides in the ancient world, 
thousands of years before present. It was initially designed to compute ephemerides for the ancient Egypt, 
but it can also be used to study astronomical phenomena in other sites of interest for cultural astronomy (archeoastronomy).
                     
For details see the `GitHub` [repository of the package](https://github.com/seap-udea/MontuPython).

*This package has been designed and written mostly by Jorge I. Zuluaga with the historical advise of egyptologist 
Francisco "Tito" Vivas, and contributions by Juanita Agudelo.*
"""),
        html.Img(src='assets/Montu.png',className='center'),
])
