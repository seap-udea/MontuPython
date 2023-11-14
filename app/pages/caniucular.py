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
This module allows you to convert between the gregorian or julian calendar and the Egyptian civil calendar, which was called by greeks the [*caniucular*](https://books.google.com.co/books?id=xKKPUpDOTKAC&pg=PA334&lpg=PA334&dq=kunikon+annus+magnus&source=bl&ots=_OUZw0Rg8m&sig=ACfU3U3M917f-R6xil1uziu4P0NMMwnsrg&hl=es-419&sa=X&ved=2ahUKEwiAsPeG18GCAxUfRTABHa3qDzQQ6AF6BAgIEAM#v=onepage&q=kunikon%20annus%20magnus&f=false) calendar since it is based on the heliacal rise of Sirus (*sothis* or *sepedet* for the Egyptian), the Dog-star.
Please provide a julian/gregorian date and the module will convert it to the corresponding caniucular date. The reverse operation is also available.
"""

module_field_doc = """
- **Quarter**: Type of quarter.
- **Datetime**: Date and time of quarter in UTC.
- **delta-t**: Time since last quarter in days.
- **delta-0**: Time since initial date in days.
- **Caniucular**: Date in *caniucular* calendar, ie. civil egyptian calendar. In this calendar the datum, namely *horus* year 0 or `hrw 0-I-Akhet-1`, corresponds to bce 2782-07-20 00:00:00.
"""

font_text = '1.2em'
font_input = '1em'

historical_dates = {
'bce 2782-07-20':'''
This is the first *apokatastais*, ie. the date when I-Akhet-1 coincides with the heliacal rise of sopedet (Sirius).
''',
'bce 1322-07-20':'''
The second *apokatastais*, ie. the date when I-Akhet-1 coincides with the heliacal rise of sopedet (Sirius), happens 
in the middle of the new Reign.
''',
'139-07-20':'''
This is the third *apokatastais*, ie. the date when I-Akhet-1 coincides with the heliacal rise of sopedet (Sirius), and it 
was identified in the writings by Censorino (Lull, p.95)
''',
'bce 559-10-19':'''
This is a date appearing in the Papyrus Louvre 7848 where we have both, the civil date and the lunar date
(see Lull, p. 94, date there III-Shemu-13)
''',
'bce 237-08-23':'''
This date comes from an inscription in the Edfu ptolemaic temple where it is evident that egyptians know the 25 years lunar cycle
(see Lull, p.92, date there III-Shemu-7).
''',
'bce 212-08-17':'''
This date comes from an inscription in the Edfu ptolemaic temple where it is evident that egyptians know the 25 years lunar cycle
(see Lull, p.92, date there III-Shemu-7).
''',
'bce 238-03-07':'''
This is the date of the Canopus decree (see Lull, p.76, date there I-Peret-17).
''',
'384-07-23':'''
Heliacal rise of sopedet (Sirius) according to computation of Theon (see Lull, p.98, date there I-Akhet-1).
''',
'bce 688-06-11':'''
This is the date appearing in the first document with a perfect chronology, namely the hieratical papyrus Louvre E3228d.
The date falls in the third year of king Taharqa (see Lull, p.105, date there I-Peret-10).
''',
}
historical_dates_options = []
for key,item in historical_dates.items():
    historical_dates_options += [dict(label=key,value=key)]

################################################################
# Layout
################################################################
dash.register_page(__name__) # Uncomment in production
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
        html.Button('Convert to caniucular',id='button-to-caniucular',value='Click',n_clicks=0,
                    style={'margin-left':'1%','font-size':font_text}),
        
        # Output
        html.Br(),
        html.Em("Calendar date:",style={'margin-left':'5em'}),
        dcc.Input(id="gdate-output",style={'margin':'1%','border':'0px','font-size':font_text,'width':'20em'}),
        
    ],style={'padding':'1%','font-size':font_text}),

    html.Div([
        html.Em("Caniucular date:",style={'margin':'0.1em'}),
        
        dcc.Input(id='cdate-hyear',  value = 0,
                  type='number',style={'margin':'1%','font-size':font_input,'width':'4em'}),"-",
        dcc.Dropdown(id='cdate-month',options=['I','II','III','IV'], value = 'I', 
                     style={'width':'5em','display':'inline-block','margin':'1%'}),"-",
        dcc.Dropdown(id='cdate-season',options=['Akhet','Peret','Shemu','Mesut'], value = 'Akhet', 
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
        html.Em("Caniucular date:",style={'margin-left':'5em'}),
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

    # Get components of caniucular
    comps_can = mtime.readable.datecan.split(' ')
    comps_can = comps_can[1].split('-')
    
    return_tuple = (
        datemix, mtime.readable.datecan,
        comps_can[0], comps_can[1], comps_can[2], comps_can[3], 
        era, comps[0], comps[1], comps[2], 
        dcc.Markdown(historical_dates[hdate])
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
    Input('button-to-caniucular', 'n_clicks'),
    State('gdate-era', 'value'),
    State('gdate-year', 'value'),
    State('gdate-month', 'value'),
    State('gdate-day', 'value'),
    State('gdate-add', 'value'),
    State('gdate-add-units', 'value'),

    prevent_initial_call=True
)
def julian_to_caniucular(button,era,year,month,day,add=0,add_units='days'):

    # Is bce?
    bce='bce' if era == 'bce' else ''

    # Prepare date
    datemix = f"{bce} {year}-{int(month):02d}-{int(day):02d}"

    # Convert date
    mtime = montu.Time(datemix,calendar='mixed')

    if int(add)>0:
        # Add 
        if add_units == 'weeks':
            factor = 7*montu.DAY
        elif add_units == 'months':
            factor = 29.5*montu.DAY
        elif add_units == 'years':
            factor = montu.JULYEAR
        else:
            factor = montu.DAY
        mtime = mtime.add(add*factor).get_readable()

    # Convert datemix
    comps = mtime.readable.datemix.split(' ')
    comps = comps[0].strip('-').split('-')
    if bce == 'bce':
        bce = 'B.C.E.'
        comps[0] = int(comps[0]) + 1
    else:
        bce = 'C.E.'
    datemix = f'{bce} {comps[0]} - {comps[1]} - {comps[2]}'

    # Get components of caniucular
    comps = mtime.readable.datecan.split(' ')
    comps = comps[1].split('-')
    
    return datemix, mtime.readable.datecan, comps[0], comps[1], comps[2], comps[3], 

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
def caniucular_to_julian(button,hyear,month,season,day,add=0,add_units='days'):

    cdate = f'hrw {hyear}-{month}-{season}-{int(day)}'
    mtime = montu.Time(cdate,calendar='caniucular')

    if int(add)>0:
        # Add 
        if add_units == 'weeks':
            factor = 10*montu.DAY
        elif add_units == 'months':
            factor = 30*montu.DAY
        elif add_units == 'years':
            factor = montu.CALYEAR
        else:
            factor = montu.DAY
        mtime = mtime.add(add*factor).get_readable()

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

    return datemix, mtime.readable.datecan, era, comps[0], comps[1], comps[2] 