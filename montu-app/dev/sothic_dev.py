import dash
import montu
from dash import Dash, html, dcc, callback, Output, Input, State
import plotly.express as px
import pandas as pd
import numpy as np

################################################################
# Preliminary data
################################################################
module_quickstart_doc = """
This module allows you to convert between the gregorian or julian calendar and the Egyptian civil calendar, which was called by greeks the [*Sothic*](https://books.google.com.co/books?id=xKKPUpDOTKAC&pg=PA334&lpg=PA334&dq=kunikon+annus+magnus&source=bl&ots=_OUZw0Rg8m&sig=ACfU3U3M917f-R6xil1uziu4P0NMMwnsrg&hl=es-419&sa=X&ved=2ahUKEwiAsPeG18GCAxUfRTABHa3qDzQQ6AF6BAgIEAM#v=onepage&q=kunikon%20annus%20magnus&f=false) calendar since it is based on the heliacal rise of Sirus (*sothis* or *sepedet* for the Egyptian), the Dog-star.
Please provide a julian/gregorian date and the module will convert it to the corresponding sothic date. The reverse operation is also available.
"""

module_field_doc = """
- **Quarter**: Type of quarter.
- **Datetime**: Date and time of quarter in UTC.
- **delta-t**: Time since last quarter in days.
- **delta-0**: Time since initial date in days.
- **Sothic**: Date in *Sothic* calendar, ie. civil egyptian calendar. In this calendar the datum, namely *horus* year 0 or `[hrw 0] I *akhet* 1`, corresponds to bce 2782-07-20 00:00:00.
"""

font_text = '1.2em'
font_input = '1em'


def _historical_markdown(entry: dict) -> str:
    parts = []
    for key in ("description", "details"):
        text = entry.get(key, "").strip()
        if text:
            parts.append(text)
    source = entry.get("source", "").strip()
    if source:
        parts.append(f"*Source: {source}*")
    return "\n\n".join(parts)


historical_dates = montu.load_historical_dates()
historical_dates_options = [
    dict(label=item.get("label", key), value=key)
    for key, item in historical_dates.items()
]

################################################################
# Layout
################################################################
#dash.register_page(__name__) # Uncomment in production
layout = html.Div([
    html.H3(children=f'Egyptian civil calendar', style={'textAlign':'center'}),
    html.Div([
        dcc.Markdown(module_quickstart_doc),
    ],style={'padding':'1%'}),
    
    html.Div([
        html.Em("Julian/Gregorian date:",style={'margin':'0.1em'}),

        dcc.RadioItems(id="gdate-era",options=['bce','ce'],value='bce',inline=True,
                       style={'width':'3em','display':'inline-block','align-items':'center','margin':'1%'}),
        dcc.Input(id='gdate-year', value='2782', 
                  type='number',style={'margin':'1%','font-size':font_input,'width':'4em'}),"-",
        dcc.Input(id='gdate-month', value='7', 
                  type='number',style={'margin':'1%','font-size':font_input,'width':'3em'}),"-",
        dcc.Input(id='gdate-day', value='20', 
                  type='number',style={'margin':'1%','font-size':font_input,'width':'3em'}),
        
        html.Em("And add:",style={'margin':'0.1em'}),
        
        dcc.Input(id='gdate-add', value='0', 
                  type='number',
                  style={'margin':'1%','font-size':font_input,'width':'5em','display':'inline-block'}),
        dcc.Dropdown(id='gdate-add-units',options=['days','weeks','months','years'],value='days',
                     style={'width':'5em','display':'inline-block','margin':'0%'}),

        # Button 
        html.Button('Convert to sothic',id='button-to-sothic',value='Click',n_clicks=0,
                    style={'margin-left':'1%','font-size':font_text}),
        
        # Output
        html.Br(),
        html.Em("Calendar date:",style={'margin-left':'5em'}),
        dcc.Input(id="gdate-output",style={'margin':'1%','border':'0px','font-size':font_text,'width':'20em'}),
        
    ],style={'padding':'1%','font-size':font_text}),

    html.Div([
        html.Em("Sothic date:",style={'margin':'0.1em'}),
        
        dcc.Input(id='cdate-hyear',  value = 0,
                  type='number',style={'margin':'1%','font-size':font_input,'width':'4em'}),"-",
        dcc.Dropdown(id='cdate-month',options=['I','II','III','IV'], value = 'I', 
                     style={'width':'3em','display':'inline-block','margin':'1%'}),"-",
        dcc.Dropdown(id='cdate-season',options=['akhet','peret','shemu','mesut'], value = 'akhet', 
                     style={'width':'6em','display':'inline-block','margin':'1%'}),"-",
        dcc.Input(id='cdate-day', value='20', 
                  type='number',style={'margin':'1%','font-size':font_input,'width':'3em'}),
        
        html.Em("And add:",style={'margin':'0.1em'}),
        
        dcc.Input(id='cdate-add', value='0', 
                  type='number',
                  style={'margin':'1%','font-size':font_input,'width':'5em','display':'inline-block'}),
        dcc.Dropdown(id='cdate-add-units',options=['days','weeks','months','years'],value='days',
                     style={'width':'5em','display':'inline-block','margin':'0%'}),
        
        html.Button('Convert to julian',id='button-to-julian',value='Click',n_clicks=0,
                    style={'margin-left':'1%','font-size':font_text}),

        # Output
        html.Br(),
        html.Em("Sothic date:",style={'margin-left':'5em'}),
        dcc.Input(id="cdate-output",style={'margin':'1%','border':'0px','font-size':font_text,'width':'20em'}),

    ],style={'padding':'1%','alignment':'center','font-size':font_text,'vertical-align': 'middle'}),

    html.Div([
        html.Em("Historical dates:",style={'margin':'0.1em'}),
        
        dcc.Dropdown(
            id='hdate',
            options = historical_dates_options,
            value = 'bce 2782-07-20', 
            style={'width':'10em','display':'inline-block','margin':'1%'}),

        html.Button('Convert',id='button-to-historical',value='Click',n_clicks=0,
                    style={'margin-left':'1%','font-size':font_text}),

        html.Div(id='hdate-explanation'),

    ],style={'padding':'1%','alignment':'center','font-size':font_text,'vertical-align': 'middle'}),
])

################################################################
# Routines
################################################################

################################################################
# Callbacks
################################################################
@callback(
    # Outputs
    Output('gdate-output', 'value', allow_duplicate=True),
    Output('cdate-output', 'value', allow_duplicate=True),
    Output('cdate-hyear', 'value', allow_duplicate=True),
    Output('cdate-month', 'value', allow_duplicate=True),
    Output('cdate-season', 'value', allow_duplicate=True),
    Output('cdate-day', 'value', allow_duplicate=True),
    Output('gdate-era', 'value', allow_duplicate=True),
    Output('gdate-year', 'value', allow_duplicate=True),
    Output('gdate-month', 'value', allow_duplicate=True),
    Output('gdate-day', 'value', allow_duplicate=True),
    Output('hdate-explanation', 'children'),

    # Inputs
    Input('button-to-historical', 'n_clicks'),
    State('hdate', 'value'),

    prevent_initial_call=True
)
def convert_historical(button,hdate):

    # Convert date
    mtime = montu.Time(hdate,calendar='mixed')

    # Convert datemix
    comps = mtime.readable.datemix.split(' ')
    comps = comps[0].strip('-').split('-')
    if mtime.bce:
        bce = 'B.C.E.'
        era = 'bce'
        comps[0] = int(comps[0]) + 1
    else:
        bce = 'C.E.'
        era = 'ce'
    datemix = f'{bce} {comps[0]} - {comps[1]} - {comps[2]}'

    # Get components of sothic
    hyear, month, season, day = montu.Time.parse_datesot(mtime.readable.datesot)
    
    return_tuple = (
        datemix, mtime.readable.datesot,
        hyear, month, season, day,
        era, comps[0], comps[1], comps[2], 
        dcc.Markdown(_historical_markdown(historical_dates[hdate]))
    )
    return return_tuple

@callback(
    # Outputs
    Output('gdate-output', 'value', allow_duplicate=True),
    Output('cdate-output', 'value', allow_duplicate=True),
    Output('cdate-hyear', 'value'),
    Output('cdate-month', 'value'),
    Output('cdate-season', 'value'),
    Output('cdate-day', 'value'),
    

    # Inputs
    Input('button-to-sothic', 'n_clicks'),
    State('gdate-era', 'value'),
    State('gdate-year', 'value'),
    State('gdate-month', 'value'),
    State('gdate-day', 'value'),
    State('gdate-add', 'value'),
    State('gdate-add-units', 'value'),

    prevent_initial_call=True
)
def julian_to_sothic(button,era,year,month,day,add=0,add_units='days'):

    # Is bce?
    bce='bce' if era == 'bce' else ''

    # Prepare date
    datemix = f"{bce} {year}-{int(month):02d}-{int(day):02d}"

    # Convert date
    mtime = montu.Time(datemix,calendar='mixed')

    if int(add)>0:
        kwargs = {'days': float(add)}
        if add_units == 'weeks':
            kwargs = {'days': float(add) * 7}
        elif add_units == 'months':
            kwargs = {'months': int(add)}
        elif add_units == 'years':
            kwargs = {'years': int(add)}
        mtime = mtime.add(**kwargs).get_readable()

    # Convert datemix
    comps = mtime.readable.datemix.split(' ')
    comps = comps[0].strip('-').split('-')
    if bce == 'bce':
        bce = 'B.C.E.'
        comps[0] = int(comps[0]) + 1
    else:
        bce = 'C.E.'
    datemix = f'{bce} {comps[0]} - {comps[1]} - {comps[2]}'

    # Get components of sothic
    hyear, month, season, day = montu.Time.parse_datesot(mtime.readable.datesot)
    
    return datemix, mtime.readable.datesot, hyear, month, season, day, 

@callback(
    # Outputs
    Output('gdate-output', 'value', allow_duplicate=True),
    Output('cdate-output', 'value', allow_duplicate=True),
    Output('gdate-era', 'value'),
    Output('gdate-year', 'value'),
    Output('gdate-month', 'value'),
    Output('gdate-day', 'value'),

    # Inputs
    Input('button-to-julian', 'n_clicks'),
    State('cdate-hyear', 'value'),
    State('cdate-month', 'value'),
    State('cdate-season', 'value'),
    State('cdate-day', 'value'),
    State('cdate-add', 'value'),
    State('cdate-add-units', 'value'),
    prevent_initial_call=True
)
def sothic_to_julian(button,hyear,month,season,day,add=0,add_units='days'):

    cdate = f'[hrw {hyear}] {month} {season.lower()} {int(day)}'
    mtime = montu.Time(cdate,calendar='sothic')

    if int(add)>0:
        kwargs = {'days': float(add)}
        if add_units == 'weeks':
            kwargs = {'days': float(add) * 7}
        elif add_units == 'months':
            kwargs = {'months': int(add)}
        elif add_units == 'years':
            kwargs = {'years': int(add)}
        mtime = mtime.add(**kwargs).get_readable()

    # Convert datemix
    comps = mtime.readable.datemix.split(' ')
    comps = comps[0].strip('-').split('-')
    if mtime.bce:
        bce = 'B.C.E.'
        era = 'bce'
        comps[0] = int(comps[0]) + 1
    else:
        bce = 'C.E.'
        era = 'ce'
        
    datemix = f'{bce} {comps[0]} - {comps[1]} - {comps[2]}'

    return datemix, mtime.readable.datesot, era, comps[0], comps[1], comps[2] 

################################################################
# Independent app code
################################################################
app = Dash(__name__)
app.layout = layout
if __name__ == '__main__':
    app.run(debug=True,port=8005)
