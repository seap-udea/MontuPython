
import dash
import dash_bootstrap_components as dbc
from dash import html, dcc, callback, page_registry, page_container
from dash.dependencies import Input, Output
from datetime import date

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    use_pages=True,
    #requests_pathname_prefix='/var/www/html/dash/MontuPython/app',
    requests_pathname_prefix='/montu-app/',
    routes_pathname_prefix='/',
)
app._favicon = ("montu.ico")

# Estilo egipcio para la barra lateral
sidebar_style = {
    'backgroundColor': '#cda434',  # Dorado egipcio
    'padding': '20px 10px',  # Aumento del padding para más espacio
    'font-family': 'Times New Roman',  # Fuente clásica
    'color': 'black',
    'border-right': '2px solid black'  # Borde para separar de la página principal
}

# Estilo para enlaces en la barra lateral
link_style = {
    'color': 'black',
    'fontWeight': 'bold',
    'padding': '8px',  # Espacio alrededor de los enlaces
    'margin': '5px 0',  # Margen entre enlaces
    'backgroundColor': 'lightgoldenrodyellow',  # Color de fondo para cada enlace
    'border-radius': '5px',  # Bordes redondeados
    'border': '1px solid black',  # Borde sutil
    'text-align': 'center',  # Centrar texto
    'text-decoration': 'none',  # Sin subrayado
    'display': 'block',  # Hacer que los enlaces llenen el espacio
}

# Estilo para los títulos y subtítulos
header_style = {
    'fontSize': 50,
    'textAlign': 'center',
    'backgroundColor': '#cda434',
    'color': 'black',
    'font-weight': 'bold'
}

subheader_style = {
    'fontSize': 32,
    'textAlign': 'center',
    'backgroundColor': '#cda434',
    'color': 'black'
}

# Creación de la barra lateral con estilo egipcio
sidebar = html.Div(
    [
        dbc.Nav(
            [dbc.NavLink(page['name'], href=page['path'], style=link_style) for page in page_registry.values()],
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
    app.run_server(debug=False, host='0.0.0.0', port=8060)

