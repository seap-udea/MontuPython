import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from utils.theme import egyptian_palette

external_stylesheets = [dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME]

dash.register_page(__name__, path="/")

layout = (
    html.Div(
        style={"max-width": "1200px", "margin": "0 auto"},
        children=[
            html.Div(
                style={
                    "display": "flex",
                    "flex-direction": "column",
                    "align-items": "center",
                    "text-align": "center",
                    "margin-bottom": "40px",
                    "padding": "40px 20px",
                    "background": f"linear-gradient(rgba(205, 164, 52, 0.2), rgba(205, 164, 52, 0.1))",
                    "border-radius": "10px",
                },
                children=[
                    html.Img(
                        src="assets/Montu.png",
                        style={"width": "200px", "margin-bottom": "20px"},
                    ),
                    html.H1(
                        "MontuPython - Ancient Astronomy",
                        style={
                            "color": egyptian_palette["text"],
                            "margin-bottom": "20px",
                            "font-size": "2.5rem",
                        },
                    ),
                    html.H3(
                        "Travel through time and explore the ancient skies",
                        style={
                            "color": egyptian_palette["secondary"],
                            "margin-bottom": "30px",
                            "font-weight": "normal",
                        },
                    ),
                    dbc.Button(
                        "Start exploration",
                        style={
                            "background-color": egyptian_palette["accent"],
                            "border": "none",
                            "padding": "10px 20px",
                        },
                        href="/moon",
                    ),
                ],
            ),
            html.Div(
                style={
                    "display": "grid",
                    "grid-template-columns": "repeat(auto-fit, minmax(300px, 1fr))",
                    "gap": "30px",
                    "margin-bottom": "40px",
                },
                children=[
                    dbc.Card(
                        style={"background-color": "white", "border": "none"},
                        children=[
                            dbc.CardBody(
                                [
                                    html.H4(
                                        "Ancient Astronomy",
                                        style={
                                            "color": egyptian_palette["accent"],
                                            "margin-bottom": "15px",
                                        },
                                    ),
                                    html.P(
                                        "Calculate precise astronomical positions for dates going back thousands of years, ideal for archaeological and cultural studies.",
                                        style={"color": egyptian_palette["text"]},
                                    ),
                                ]
                            )
                        ],
                    ),
                    dbc.Card(
                        style={"background-color": "white", "border": "none"},
                        children=[
                            dbc.CardBody(
                                [
                                    html.H4(
                                        "Egyptian Focus",
                                        style={
                                            "color": egyptian_palette["accent"],
                                            "margin-bottom": "15px",
                                        },
                                    ),
                                    html.P(
                                        "Originally developed for the study of astronomy in ancient Egypt, with specific features for this civilization.",
                                        style={"color": egyptian_palette["text"]},
                                    ),
                                ]
                            )
                        ],
                    ),
                    dbc.Card(
                        style={"background-color": "white", "border": "none"},
                        children=[
                            dbc.CardBody(
                                [
                                    html.H4(
                                        "Universal Application",
                                        style={
                                            "color": egyptian_palette["accent"],
                                            "margin-bottom": "15px",
                                        },
                                    ),
                                    html.P(
                                        "Applicable to any archaeoastronomical site or cultural study that requires precise historical astronomical calculations.",
                                        style={"color": egyptian_palette["text"]},
                                    ),
                                ]
                            )
                        ],
                    ),
                ],
            ),
            html.Div(
                style={
                    "background-color": "white",
                    "padding": "30px",
                    "border-radius": "10px",
                    "margin-bottom": "40px",
                },
                children=[
                    dcc.Markdown(
                        """
            ## About MontuPython

            **MontuPython** (transliterated mnṯw ꜥꜣpp(y)) is your gateway to ancient astronomy.

            *MontuPython* is a Python package designed to compute astronomical ephemerides thousands of years before the present day. It was originally developed to calculate ephemerides for ancient Egypt, but it can also be applied to the study of astronomical phenomena in other regions of cultural interest, such as archaeoastronomy sites.

            Explore the ancient skies and unravel the mysteries of the past with MontuPython!
            """,
                        style={"color": egyptian_palette["text"]},
                    ),
                    html.Div(
                        style={
                            "display": "flex",
                            "justify-content": "center",
                            "margin-top": "20px",
                        },
                        children=[
                            dbc.Button(
                                [
                                    "View on GitHub",
                                    html.I(className="bi bi-github"),
                                ],
                                href="https://github.com/seap-udea/MontuPython",
                                target="_blank",
                                style={
                                    "background-color": egyptian_palette["secondary"],
                                    "border": "none",
                                    "margin-right": "10px",
                                    "display": "flex",
                                    "gap": "8px",
                                    "align-items": "center",
                                },
                            ),
                        ],
                    ),
                ],
            ),
            html.Div(
                style={
                    "background-color": "rgba(205, 164, 52, 0.1)",
                    "padding": "20px",
                    "border-radius": "10px",
                    "text-align": "center",
                },
                children=[
                    html.H4(
                        "Developers",
                        style={
                            "color": egyptian_palette["secondary"],
                            "margin-bottom": "15px",
                        },
                    ),
                    html.P(
                        [
                            "This package has been primarily developed by ",
                            html.Strong("Jorge I. Zuluaga"),
                            ", with historical guidance from Egyptologist ",
                            html.Strong('Francisco "Tito" Vivas'),
                            " and contributions by ",
                            html.Strong("Juanita Agudelo"),
                            ".",
                        ],
                        style={"color": egyptian_palette["text"]},
                    ),
                ],
            ),
        ],
    ),
)
