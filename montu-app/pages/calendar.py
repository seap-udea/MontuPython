import dash
import montu
from dash import html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
from utils.theme import egyptian_palette

dash.register_page(__name__)

# Agregar estilos consistentes al inicio del archivo
style_input = {
    'border': f'1px solid {egyptian_palette["accent"]}',
    'width': '100%',
    'borderRadius': '4px'
}

button_style = {
    'backgroundColor': egyptian_palette["accent"],
    'border': 'none',
    'color': 'white',
    'borderRadius': '4px',
    'padding': '8px 16px'
}

table_style = {
    'borderColor': egyptian_palette["accent"],
    'backgroundColor': 'white',
}

layout = html.Div(
    children=[
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H3(
                        "Date Converter",
                        className="text-center",
                        style={"color": egyptian_palette["text"]},
                    )
                ),
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.ButtonGroup(
                                            [
                                                dbc.Button(
                                                    "Calendar",
                                                    id="calendar-btn",
                                                    color="primary",
                                                    active=True,
                                                    className="me-1",
                                                    n_clicks=0,
                                                    style=button_style
                                                ),
                                                dbc.Button(
                                                    "Julian day",
                                                    id="julian-btn",
                                                    color="primary",
                                                    outline=True,
                                                    className="ms-1",
                                                    n_clicks=0,
                                                    style=button_style
                                                ),
                                            ],
                                            className="d-flex justify-content-center mb-4",
                                        ),
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Calendar Type:", className="fw-bold me-2"),
                                                        dbc.RadioItems(
                                                            options=[
                                                                {"label": "Gregorian proleptic", "value": "proleptic"},
                                                                {"label": "Mixed", "value": "mixed"},
                                                            ],
                                                            value="mixed",
                                                            id="calendar-radio-input",
                                                            inline=True,
                                                        ),
                                                    ],
                                                    className="d-flex align-items-center justify-content-center mb-4",
                                                ),
                                            ]
                                        ),
                                        html.Div(id="output"),
                                    ],
                                    width=12,
                                ),
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Card(
                                            [
                                                dbc.CardHeader(
                                                    html.H5("Date Information", className="text-center")
                                                ),
                                                dbc.CardBody(
                                                    html.Div(id="date-output-container")
                                                ),
                                            ],
                                        )
                                    ],
                                    width=12,
                                ),
                            ],
                            className="mt-4",
                        ),
                        dcc.Store(id="input-type-store", data={"value": "calendar"}),
                        dcc.Input(id="year", type="hidden", value=2000),
                        dcc.Input(id="month", type="hidden", value=1),
                        dcc.Input(id="day", type="hidden", value=1),
                        dcc.Input(id="hour", type="hidden", value=0),
                        dcc.Input(id="min", type="hidden", value=0),
                        dcc.Input(id="sec", type="hidden", value=0),
                    ]
                ),
            ],
            style={"backgroundColor": egyptian_palette["background"]},
        )
    ],
    className="p-4",
)

# Actualizar el componente input_date
input_date = html.Div(
    [
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.H5("Date", className="mb-3 text-center fw-bold"),
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Year", html_for="year-input", className="fw-bold"),
                                                        dbc.Input(
                                                            id="year",
                                                            type="number",
                                                            value=-1321,
                                                            placeholder="Year",
                                                            className="mb-2",
                                                            style=style_input
                                                        ),
                                                    ],
                                                    width=4,
                                                ),
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Month", html_for="month-input", className="fw-bold"),
                                                        dbc.Input(
                                                            id="month",
                                                            type="number",
                                                            min=1,
                                                            max=12,
                                                            value=7,
                                                            placeholder="Month",
                                                            className="mb-2",
                                                            style=style_input
                                                        ),
                                                    ],
                                                    width=4,
                                                ),
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Day", html_for="day-input", className="fw-bold"),
                                                        dbc.Input(
                                                            id="day",
                                                            type="number",
                                                            placeholder="Day",
                                                            value=20,
                                                            min=1,
                                                            max=31,
                                                            step=1,
                                                            className="mb-2",
                                                            style=style_input
                                                        ),
                                                    ],
                                                    width=4,
                                                ),
                                            ],
                                        ),
                                    ],
                                    width=12,
                                    lg=6,
                                    className="mb-3 pe-lg-4 border-end-lg",
                                ),
                                dbc.Col(
                                    [
                                        html.H5("Time", className="mb-3 text-center fw-bold"),
                                        dbc.Row(
                                            [
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Hour", html_for="hour-input", className="fw-bold"),
                                                        dbc.Input(
                                                            id="hour",
                                                            type="number",
                                                            value=0,
                                                            min=0,
                                                            max=24,
                                                            placeholder="Hour",
                                                            className="mb-2",
                                                            style=style_input
                                                        ),
                                                    ],
                                                    width=4,
                                                ),
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Minute", html_for="min-input", className="fw-bold"),
                                                        dbc.Input(
                                                            id="min",
                                                            type="number",
                                                            value=0,
                                                            min=0,
                                                            max=59,
                                                            debounce=True,
                                                            placeholder="Min",
                                                            className="mb-2",
                                                            style=style_input
                                                        ),
                                                    ],
                                                    width=4,
                                                ),
                                                dbc.Col(
                                                    [
                                                        dbc.Label("Second", html_for="sec-input", className="fw-bold"),
                                                        dbc.Input(
                                                            id="sec",
                                                            type="number",
                                                            value=0,
                                                            min=0,
                                                            max=59,
                                                            placeholder="Sec",
                                                            step=1,
                                                            className="mb-2",
                                                            style=style_input
                                                        ),
                                                    ],
                                                    width=4,
                                                ),
                                            ],
                                        ),
                                    ],
                                    width=12,
                                    lg=6,
                                    className="ps-lg-4",
                                ),
                            ],
                            className="g-3",
                        ),
                    ]
                ),
            ],
            style={"backgroundColor": egyptian_palette["background"]},
        ),
    ],
)

@callback(
    Output("date-output-container", "children"),
    Input("year", "value"),
    Input("month", "value"),
    Input("day", "value"),
    Input("hour", "value"),
    Input("min", "value"),
    Input("sec", "value"),
    Input("calendar-radio-input", "value"),
)
def convert_date_calendar(year, month, day, hour, min, sec, calendar):
    date = f"{int(year or 1)}-{int(month or 1):02d}-{int(day or 1):02d} {int(hour or 0):02d}:{int(min or 0):02d}:{int(sec or 0):02d}"
    mtime = montu.Time(date, format="iso", scale="utc", calendar=calendar)
    return get_date_output_table(mtime)

def get_date_output_table(mtime):
    spice = mtime.readable.datespice
    proleptic = mtime.readable.datepro
    mixed = mtime.readable.datemix
    caniucular = mtime.readable.datecan
    jd_utc = mtime.jed
    jd = mtime.tt
    es = mtime.et
    delta_t = mtime.deltat

    table_props = {
        'bordered': True,
        'hover': True,
        'responsive': True,
        'striped': True,
        'className': "mb-4",
        'style': table_style
    }

    date_formats_table = dbc.Table(
        [
            html.Thead(
                html.Tr([
                    html.Th("Format", style={"width": "40%"}),
                    html.Th("Value"),
                ])
            ),
            html.Tbody([
                html.Tr([
                    html.Td(html.Strong("Gregorian proleptic (human readable)")),
                    html.Td(spice),
                ]),
                html.Tr([
                    html.Td(html.Strong("Gregorian Proleptic (astronomical)")),
                    html.Td(proleptic),
                ]),
                html.Tr([
                    html.Td(html.Strong("Mixed (gregorian or Julian)")),
                    html.Td(mixed),
                ]),
                html.Tr([
                    html.Td(html.Strong("Caniucular (civil egyptian)")),
                    html.Td(caniucular),
                ]),
            ]),
        ],
        **table_props
    )

    julian_info_table = dbc.Table(
        [
            html.Thead(
                html.Tr([
                    html.Th("Measurement", style={"width": "40%"}),
                    html.Th("Value"),
                ])
            ),
            html.Tbody([
                html.Tr([
                    html.Td(html.Strong("Julian day (UTC)")),
                    html.Td(f"{jd_utc} days"),
                ]),
                html.Tr([
                    html.Td(html.Strong("Julian day (TT)")),
                    html.Td(f"{jd} days"),
                ]),
                html.Tr([
                    html.Td(html.Strong("Ephemerides seconds (TT)")),
                    html.Td(f"{es} seconds"),
                ]),
                html.Tr([
                    html.Td(html.Strong("Delta-t")),
                    html.Td(f"{delta_t} seconds"),
                ]),
            ]),
        ],
        **table_props
    )

    return html.Div([
        html.H6("Date Formats", className="mt-2 mb-3"),
        date_formats_table,
        html.H6("Julian Day Information", className="mt-4 mb-3"),
        julian_info_table,
    ])

@callback(
    [Output("output", "children"), 
     Output("calendar-btn", "active"),
     Output("julian-btn", "active"),
     Output("calendar-btn", "outline"),
     Output("julian-btn", "outline"),
     Output("input-type-store", "data")],
    [Input("calendar-btn", "n_clicks"),
     Input("julian-btn", "n_clicks")],
    [Input("input-type-store", "data")]
)
def update_input_type(calendar_clicks, julian_clicks, stored_data):
    ctx = dash.callback_context
    
    if not ctx.triggered:
        input_type = "calendar"
    else:
        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == "calendar-btn":
            input_type = "calendar"
        elif button_id == "julian-btn":
            input_type = "julian"
        else:
            input_type = stored_data["value"]
    
    calendar_active = input_type == "calendar"
    julian_active = input_type == "julian"
    
    return (
        input_date,
        calendar_active,
        julian_active,
        not calendar_active,
        not julian_active,
        {"value": input_type}
    )
