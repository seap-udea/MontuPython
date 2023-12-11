import dash 
from dash import dcc, html

# Paleta de colores egipcios
egyptian_palette = {
    'background': '#f5e2a1',  # Fondo amarillo claro
    'text': '#000000',  # Texto en negro
    'header': '#cda434',  # Encabezado en dorado egipcio
    'accent': '#d97824',  # Acento en naranja oscuro
}

dash.register_page(__name__, path='/')

layout = html.Div(style={'background-color': egyptian_palette['background']}, children=[
    html.Img(src='assets/Montu.png', className='center', style={'width': '300px'}),
    html.Header(style={'background-color': egyptian_palette['header'], 'padding': '20px'}, children=[
        html.H1("MontuPython - Ancient Astronomy", style={'color': egyptian_palette['text']}),
    ]),
    dcc.Markdown("""

    Welcome to `MontuPython` (transliterated mnṯw ꜥꜣpp(y)), your gateway to ancient astronomy.

    *MontuPython* is a Python package designed to compute astronomical ephemerides thousands of years before the present day. It was originally developed to calculate ephemerides for ancient Egypt, but it can also be applied to the study of astronomical phenomena in other regions of cultural interest, such as archaeoastronomy sites.

    Explore the ancient skies and unravel the mysteries of the past with MontuPython!

    For more details, visit the [GitHub repository](https://github.com/seap-udea/MontuPython).

    *This package has been primarily developed by Jorge I. Zuluaga, with historical guidance from Egyptologist Francisco "Tito" Vivas and contributions by Juanita Agudelo.*

    """, style={'color': egyptian_palette['text'], 'padding': '20px'}),
])
