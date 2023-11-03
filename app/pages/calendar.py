import dash
from dash import html, dcc, callback, Input, Output
import dash_bootstrap_components as dbc
from datetime import date
import montu as mn

dash.register_page(__name__)

def convert(n):
    if n < 10:
        return '0'+str(n)
    else: 
        return  str(n)

button_group = html.Div(
    [
        dbc.RadioItems(
            id="radios",
            className="btn-group",
            inputClassName="btn-check",
            labelClassName="btn btn-outline-primary",
            labelCheckedClassName="active",
            options=[
                {"label": "Calendar", "value": "calendar"},
                {"label": "Julian day", "value": "julian day"},
                {"label": "Sothic", "value": "sothic"}
            ],
            value="calendar",
        ),
        html.Br(),
        html.Div(id="output")
    ],
    className="radio-group",
)

radioitems = html.Div(
    [
        #dbc.Label("Choose one"),
        dbc.RadioItems(
            options=[
                {"label": "Proleptic Date", "value": "proleptic"},
                {"label": "Mixed Date", "value": "mixed"}
            ],
            value="proleptic",
            id="radioitems-input",
            inline=True
        ),
        html.Br()
    ]
)

calendar_mixed = html.Div([
        dcc.DatePickerSingle(
        id='my-date-picker-single',
        initial_visible_month=date.today(),
        date=date.today()
    ),
    html.Div([
    html.H5('Hour1',{'display':'inline-block', 'border': '1px solid black'}),
    ], style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
    dcc.Input(id="hour1", type="number", min=1, max=24, placeholder="Hour", style={'marginRight':'10px', 'width':'5%'}),
    html.Div([
    html.H6(':',{'display':'inline-block', 'border': '1px solid black'}),
    ], style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
    dcc.Input(id="min1", type="number", min=1, max=60, debounce=True, placeholder="Min", style={'marginRight':'10px', 'width':'5%'}),
    html.Div([
    html.H6(':',{'display':'inline-block', 'border': '1px solid black'}),
    ], style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
    dcc.Input(id="sec1", type="number", min=1, max=60, placeholder="Sec", step=1),
    html.Div(id='output-container-date-picker-single')
])


input_date = html.Div(
    [
        html.Div([
        html.H5('Date',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':20, 'margin-left':20, 'margin-top':10}),
        
        dcc.Input(id="year", type="number", placeholder="Year", style={'marginRight':'10px', 'width':'10%', 'marginTop':'10px'}),
        dcc.Input(id="month", type="number", min=1, max=12, debounce=True, placeholder="Month", style={'marginRight':'10px', 'width':'10%','marginTop':'10px'}),
        dcc.Input(id="day", type="number", placeholder="Day", min=1, max=31, step=1, style={'marginTop':'10px'}),
        
        html.Div([
        html.H5('Time',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
        dcc.Input(id="hour", type="number", min=1, max=24, placeholder="Hour", style={'width':'7%'}),
        
        html.Div([
        html.H6(':',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':5, 'margin-left':5}),
        dcc.Input(id="min", type="number", min=0, max=59, debounce=True, placeholder="Min", style={'width':'7%'}),
        
        html.Div([
        html.H6(':',{'display':'inline-block', 'border': '1px solid black'}),
        ], style={'display':'inline-block', 'margin-right':5, 'margin-left':5}),
        dcc.Input(id="sec", type="number", min=0, max=59, placeholder="Sec", step=1, style={'width':'7%'}),
        #html.Div(id="date-out")
    ]
)

julian_date = html.Div([
    html.Div([
    html.H5('Julian Date',{'display':'inline-block', 'border': '1px solid black'})], style={'display':'inline-block', 'margin-right':20, 'margin-left':20}),
    dcc.Input(id="julian", type="number", placeholder="Julian Day", style={'width':'20%', 'marginTop':'10px'}),
])

button = dbc.Button("Submit", id="button", n_clicks=0)

result = html.Div(id="final_output", style={'white-space': 'pre-line'})

#----------------------------------------------------------------------------------------------------------------------

layout = html.Div([
    html.H5('Date Converter', className='h5'),
    html.P('Use this date converter to convert proleptic dates or mixed gregorian dates to ...', className='p'),
    html.Center(button_group), 
    html.Center(radioitems), 
    html.Center(button), 
    html.Hr(),
    result,
    html.Hr()
    ])

#----------------------------------------------------------------------------------------------------------------------
@callback(
    Output("radioitems-checklist-output", "children"),
    [Input("radioitems-input", "value")],
)

def on_form_change(radio_items_value):
    return radio_items_value

@callback(Output("output", "children"), [Input("radios", "value")])

def display_value(value):
    if value == "calendar":
        return input_date 
    if value == "julian day":
        return julian_date
    else :
        return calendar_mixed

@callback(
    Output('output-container-date-picker-single', 'children'),
    Input('my-date-picker-single', 'date'),
    Input("hour1", "value"),
    Input("min1", "value"),
    Input("sec1", "value"),
    )

def update_output(date_value, hour1, min1, sec1):
    string_prefix = 'You have selected: '
    if date_value is not None:
        date_object = date.fromisoformat(date_value)
        date_string = date_object.strftime('%B %d, %Y')
        time = str(hour1) + ':' + str(min1) + ':' + str(sec1)
        return string_prefix + date_string + time

@callback(
    Output("day", "max"),
    [Input("month", "value")]
)
def update_max_day(month):
    if month in [4,6,9,11]:
        return 30
    elif month == 2:
        return 28
    else:
        return 31
    
@callback(
    Output("date-out", "children"),
    Input("year", "value"),
    Input("month", "value"),
    Input("day", "value"),
    Input("hour", "value"),
    Input("min", "value"),
    Input("sec", "value"),
    Input("julian", "value"),
)
def number_render(year, month, day, hour, min, sec):
    return "year: {}, month: {}, day: {}, hour: {}, min: {}, sec: {}".format(year, month, day, hour, min, sec)

@callback(
    Output("final_output", "children"),
    Input("button", "n_clicks"),
    Input("year", "value"),
    Input("month", "value"),
    Input("day", "value"),
    Input("hour", "value"),
    Input("min", "value"),
    Input("sec", "value"),
    Input("radioitems-input", "value")
)

def subbit_button(button, year, month, day, hour, min, sec, radioitems):

    year = str(year)
    month = convert(month)
    day = convert(day)

    hour = convert(hour)
    min = convert(min)
    sec = convert(sec)

    date = year + '-' + month + '-' + day + ' ' + hour + ':' + min + ':' + sec  
    mtime = mn.Time(date, format='iso',scale='utc',calendar=radioitems)

    proleptic = mtime.readable.datepro, # Date in gregorian proleptic
    mixed = mtime.readable.datemix, # Date in gregorian mixed
    jd_utc = mtime.jed, # Date in Julian Day (utc scale)
    jd = mtime.tt, # Date in ephemerides time (tt scale)
    es = mtime.et, # Date in ephemerides time (utc scale)
    delta_t = mtime.deltat

    result = f'''Date Proleptic : {proleptic} 
    Date Mixed : {mixed} 
    JD (UTC) : {jd_utc} 
    JD : {jd} 
    Ephemerides seconds : {es}
    Delta-t : {delta_t}'''
    return result
