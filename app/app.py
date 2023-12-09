import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, page_registry, page_container
from dash.dependencies import Input, Output
from datetime import date

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], use_pages=True)
app._favicon = ("montu.ico")

# Estilo para la barra lateral
sidebar_style = {
    'backgroundColor': '#f4f4f2', # Color de fondo claro
    'padding': '10px'
}

# Estilo para los títulos y subtítulos
header_style = {
    'fontSize': 50,  # Mantenemos el tamaño original
    'textAlign': 'center',
    'backgroundColor': '#cda434',  # Dorado egipcio
    'color': 'black',  # Texto en negro para contraste
    'font-weight': 'bold'
}

subheader_style = {
    'fontSize': 32,  # Mantenemos el tamaño original
    'textAlign': 'center',
    'backgroundColor': '#cda434',  # Dorado egipcio
    'color': 'black'
}

# Creación de la barra lateral
sidebar = html.Div(
    [
        dbc.Nav(
            [dbc.NavLink(page['name'], href=page['path']) for page in page_registry.values()],
            vertical=True,
            pills=True,
            style=sidebar_style,
        )
    ]
)

# Layout principal de la aplicación
app.layout = dbc.Container([
    dbc.Row([
        dbc.Row(html.Div("")),
        dbc.Row(html.Div("Montu App", style=header_style)),
        dbc.Row(html.Div("Astronomical ephemerides for the ancient world", style=subheader_style))
    ]),
    dbc.Row(
        [
            dbc.Col(sidebar, xs=4, sm=4, md=2, lg=2, xl=2, xxl=2),
            dbc.Col(dash.page_container, xs=8, sm=8, md=10, lg=10, xl=10, xxl=10)
        ]
    )
], fluid=True)

server = app.server

if __name__ == '__main__':
    app.run_server(debug=False, port=8080)
