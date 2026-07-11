# %pip install montu plotly pandas numpy

import montu
import numpy as np
import pandas as pd
import plotly.express as px

INITIAL_DATE = "-1500-01-01"
TIME_SPAN_YEARS = 10
NUM_POINTS = 120
PLANETS = ["Mercury", "Venus"]
PROPERTY = "DecEpoch"

## Observer at Thebes
observer = montu.Observer(lon=33, lat=24)

# Build planet list (all except Sun, Moon, Earth)
all_planets = []
for value in montu.PLANETARY_NAMES.values():
    if value not in ("SUN", "MOON", "EARTH"):
        all_planets.append(montu.Planet(value))

# Sample dates across the time span
## Starting epoch for the ephemeris loop
mtime = montu.Time(INITIAL_DATE)
mts = []
dates = []
for dt in np.linspace(0, TIME_SPAN_YEARS * montu.YEAR, NUM_POINTS):
    ## Advance time by dt seconds and fill readable date strings (year, month, day, …)
    mt = (mtime + dt).get_readable()
    mts.append(mt)
    dates.append(f"{mt.readable.year}-{mt.readable.month}-{mt.readable.day}")

# Ephemeris for each planet
ephemerides = pd.DataFrame()
for planet in all_planets:
    planet.reset_store()
    for mt in mts:
        ## Store sky position at each sampled date
        planet.conditions_in_sky(at=mt, observer=observer, store=True)
    ## Build a table of stored ephemeris rows
    planet.tabulate_ephemerides()
    planet.ephemerides["datestr"] = dates
    ephemerides = pd.concat([ephemerides, planet.ephemerides], ignore_index=True)

mask = ephemerides.Name.isin(PLANETS)
fig = px.line(ephemerides[mask], x="datestr", y=PROPERTY, color="Name")
## End date of the span, with readable calendar strings filled in
mtime_final = (mtime + TIME_SPAN_YEARS * montu.YEAR).get_readable()
fig.update_layout(
    title=f"{PROPERTY}  ·  {', '.join(PLANETS)}  ·  {INITIAL_DATE} to {mtime_final.strftime('%Y-%m-%d')}",
    xaxis_title="Date",
    xaxis=dict(rangeslider=dict(visible=True)),
)
fig.show()
