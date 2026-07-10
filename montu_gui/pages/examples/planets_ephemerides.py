"""
MontuPython planetary ephemerides example.

Reproduces the chart logic of the Planetary Ephemerides module:
sample one ephemeris property for one or more planets over a time span.

Run:
    python montu_planets_ephemerides.py
"""

import montu
import numpy as np
import pandas as pd
import plotly.express as px

# ─── Parameters (edit as needed) ─────────────────────────────────────────────

INITIAL_DATE = "-1500-01-01"
TIME_SPAN_YEARS = 10
NUM_POINTS = 120
PLANETS = ["Mercury", "Venus"]
PROPERTY = "DecEpoch"

OBSERVER = montu.Observer(lon=33, lat=24)  # Thebes

# ─── Compute ephemerides ─────────────────────────────────────────────────────

all_planets = [
    montu.Planet(value)
    for value in montu.PLANETARY_NAMES.values()
    if value not in ("SUN", "MOON", "EARTH")
]

mtime = montu.Time(INITIAL_DATE)
mts = []
dates = []
for dt in np.linspace(0, TIME_SPAN_YEARS * montu.YEAR, NUM_POINTS):
    mt = (mtime + dt).get_readable()
    mts.append(mt)
    dates.append(f"{mt.readable.year}-{mt.readable.month}-{mt.readable.day}")

ephemerides = pd.DataFrame()
for planet in all_planets:
    planet.reset_store()
    for mt in mts:
        planet.conditions_in_sky(at=mt, observer=OBSERVER, store=True)
    planet.tabulate_ephemerides()
    planet.ephemerides["datestr"] = dates
    ephemerides = pd.concat([ephemerides, planet.ephemerides], ignore_index=True)

# ─── Plot ────────────────────────────────────────────────────────────────────

mask = ephemerides.Name.isin(PLANETS)
fig = px.line(ephemerides[mask], x="datestr", y=PROPERTY, color="Name")
mtime_final = (mtime + TIME_SPAN_YEARS * montu.YEAR).get_readable()
fig.update_layout(
    title=(
        f"Property '{PROPERTY}' of {', '.join(PLANETS)} "
        f"from {mtime.strftime('%Y-%m-%d')} to {mtime_final.strftime('%Y-%m-%d')}"
    ),
    xaxis_title="Date [Month & Year]",
    xaxis=dict(rangeslider=dict(visible=True)),
)
fig.show()
