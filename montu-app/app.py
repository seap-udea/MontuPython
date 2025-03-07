import dash
import dash_bootstrap_components as dbc
from dash import (
    html,
    dcc,
    Input,
    Output,
    State,
)
from datetime import date
from utils.theme import egyptian_palette


app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.BOOTSTRAP],
    use_pages=True,
    suppress_callback_exceptions=True,
)
app._favicon = "montu.ico"

menu_item_style = {
    "backgroundColor": "lightyellow",
    "margin": "5px 10px",
    "borderRadius": "5px",
    "color": "black",
    "fontWeight": "bold",
    "textAlign": "center",
    "padding": "10px",
}

sidebar = html.Div(
    [
        html.Div(
            [
                html.Img(
                    src="assets/Montu.png",
                    style={"width": "80px", "margin": "0 auto", "display": "block"},
                ),
                html.H4(
                    "Montu App",
                    className="text-center mt-2",
                    style={"color": "black"},
                ),
            ],
            style={"padding": "20px 0", "background-color": egyptian_palette["header"]},
        ),
        dbc.Nav(
            [
                dbc.NavLink(
                    [html.Div("Home")], href="/", active="exact", style=menu_item_style
                ),
                dbc.NavLink(
                    [html.Div("Calendar")],
                    href="/calendar",
                    active="exact",
                    style=menu_item_style,
                ),
                dbc.NavLink(
                    [html.Div("Caniucular")],
                    href="/caniucular",
                    active="exact",
                    style=menu_item_style,
                ),
                dbc.NavLink(
                    [html.Div("Moon")],
                    href="/moon",
                    active="exact",
                    style=menu_item_style,
                ),
                dbc.NavLink(
                    [html.Div("Planets")],
                    href="/planets",
                    active="exact",
                    style=menu_item_style,
                ),
                dbc.NavLink(
                    [html.Div("Stars")],
                    href="/stars",
                    active="exact",
                    style=menu_item_style,
                ),
            ],
            vertical=True,
            pills=True,
            style={"padding-bottom": "20px"},
        ),
    ],
    id="sidebar",
    style={
        "position": "fixed",
        "top": 0,
        "left": 0,
        "bottom": 0,
        "width": "250px",
        "background-color": egyptian_palette["header"],
        "transition": "all 0.3s",
        "zIndex": 1050,
        "boxShadow": "2px 0 5px rgba(0,0,0,0.2)",
        "transform": "translateX(-250px)",  # Inicialmente oculto
    },
)

navbar = html.Nav(
    style={
        "background-color": egyptian_palette["header"],
        "padding": "10px 20px",
        "display": "flex",
        "justify-content": "center",
        "align-items": "center",
        "position": "fixed",
        "top": 0,
        "left": 0,
        "right": 0,
        "zIndex": 1040,
        "boxShadow": "0 2px 5px rgba(0,0,0,0.1)",
    },
    children=[
        html.Div(
            style={
                "display": "flex",
                "align-items": "center",
                "justify-content": "center",
                "position": "relative",
                "width": "100vw",
            },
            children=[
                html.Button(
                    [
                        html.Div(
                            className="menu-icon",
                            id="menu-icon-bars",
                            children=[
                                html.Span(className="bar1"),
                                html.Span(className="bar2"),
                                html.Span(className="bar3"),
                            ],
                        )
                    ],
                    id="sidebar-toggle",
                    style={
                        "position": "absolute",
                        "zIndex": 1060,
                        "backgroundColor": egyptian_palette["secondary"],
                        "borderRadius": "8px",
                        "width": "45px",
                        "height": "45px",
                        "border": "none",
                        "cursor": "pointer",
                        "display": "flex",
                        "alignItems": "center",
                        "justifyContent": "center",
                        "boxShadow": "0 2px 5px rgba(0,0,0,0.2)",
                        "transition": "all 0.3s ease",
                    },
                ),
                html.Div(
                    style={
                        "display": "flex",
                        "alignItems": "center",
                        "justifyContent": "center",
                        "flexDirection": "column",
                    },
                    children=[
                        html.H3(
                            "Montu App",
                            style={
                                "color": egyptian_palette["text"],
                                "margin": "0",
                            },
                        ),
                        html.P("Astronomical ephemerides for the ancient world"),
                    ],
                ),
            ],
        ),
    ],
)

app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            /* Estilos para el icono de menú hamburguesa */
            .menu-icon {
                width: 25px;
                height: 20px;
                position: relative;
            }
            
            .menu-icon .bar1, .menu-icon .bar2, .menu-icon .bar3 {
                width: 100%;
                height: 3px;
                background-color: white;
                position: absolute;
                left: 0;
                border-radius: 3px;
                transition: all 0.3s ease;
            }
            
            .menu-icon .bar1 {
                top: 0;
            }
            
            .menu-icon .bar2 {
                top: 8px;
            }
            
            .menu-icon .bar3 {
                top: 16px;
            }
            
            /* Animación cuando el menú está abierto */
            .menu-icon.open .bar1 {
                transform: rotate(-45deg) translate(-5px, 6px);
            }
            
            .menu-icon.open .bar2 {
                opacity: 0;
            }
            
            .menu-icon.open .bar3 {
                transform: rotate(45deg) translate(-5px, -6px);
            }
            
            /* Iconos para los estados abierto y cerrado */
            .menu-icon:not(.open) .bar1, 
            .menu-icon:not(.open) .bar2, 
            .menu-icon:not(.open) .bar3 {
                /* Icono de menú hamburguesa (tres barras) cuando está cerrado */
                display: block;
            }
            
            .menu-icon.open .bar1, 
            .menu-icon.open .bar2, 
            .menu-icon.open .bar3 {
                /* Icono X cuando está abierto */
                display: block;
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""

app.layout = html.Div(
    [
        sidebar,
        navbar,
        dbc.Container(
            [
                html.Div(style={"height": "70px"}),
                dbc.Row(
                    html.Div(
                        style={
                            "background-color": egyptian_palette["background"],
                            "min-height": "100vh",
                            "font-family": "Roboto, sans-serif",
                            "margin-top": "20px",
                        },
                        children=[dash.page_container],
                    ),
                    id="page-content",
                ),
                dbc.Row(
                    html.Footer(
                        style={
                            "background-color": egyptian_palette["secondary"],
                            "color": "#ffffff",
                            "text-align": "center",
                            "padding": "20px",
                        },
                        children=[
                            html.P(
                                "© 2023 MontuPython - Universidad de Antioquia",
                                style={"font-weight": "bold", "margin-bottom": "10px"},
                            ),
                            html.P(
                                [
                                    "License: ",
                                    html.A(
                                        "MIT",
                                        href="https://github.com/seap-udea/MontuPython/blob/main/LICENSE",
                                        target="_blank",
                                        style={
                                            "color": egyptian_palette["highlight"],
                                            "text-decoration": "underline",
                                            "font-weight": "bold",
                                        },
                                    ),
                                ]
                            ),
                        ],
                    ),
                ),
            ],
            className="p-0",
            fluid=True,
            id="main-content",
            style={"transition": "all 0.3s"},
        ),
        dcc.Store(id="sidebar-state", data={"visible": False}),
    ]
)


@app.callback(
    [
        Output("sidebar", "style"),
        Output("main-content", "style"),
        Output("sidebar-state", "data"),
        Output("menu-icon-bars", "className"),
        Output("sidebar-toggle", "style"),
    ],
    [Input("sidebar-toggle", "n_clicks")],
    [State("sidebar-state", "data")],
)
def toggle_sidebar(n_clicks, sidebar_state):
    if n_clicks is None:
        visible = False
    else:
        visible = not sidebar_state.get("visible", False)

    sidebar_style = {
        "position": "fixed",
        "top": 0,
        "left": 0,
        "bottom": 0,
        "width": "250px",
        "background-color": egyptian_palette["header"],
        "transition": "all 0.3s",
        "zIndex": 1050,
        "boxShadow": "2px 0 5px rgba(0,0,0,0.2)",
        "transform": "translateX(0px)" if visible else "translateX(-250px)",
    }

    main_content_style = {
        "transition": "all 0.3s",
        "marginLeft": "250px" if visible else "0px",
    }

    # Clase para la animación del icono de menú
    menu_icon_class = "menu-icon open" if visible else "menu-icon"

    # Estilo del botón (se mueve cuando el sidebar está abierto)
    button_style = {
        "position": "absolute",
        "left": (
            "250px" if visible else "0"
        ),  # Se mueve a la derecha cuando el menú está abierto
        "zIndex": 1060,
        "backgroundColor": egyptian_palette["secondary"],
        "borderRadius": "8px",
        "width": "45px",
        "height": "45px",
        "border": "none",
        "cursor": "pointer",
        "display": "flex",
        "alignItems": "center",
        "justifyContent": "center",
        "boxShadow": "0 2px 5px rgba(0,0,0,0.2)",
        "transition": "all 0.3s ease",
    }

    return (
        sidebar_style,
        main_content_style,
        {"visible": visible},
        menu_icon_class,
        button_style,
    )


server = app.server
if __name__ == "__main__":
    app.run_server(debug=True, host="0.0.0.0", port=8060)
