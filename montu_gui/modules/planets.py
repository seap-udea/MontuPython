"""
Planetary ephemerides logic for MontuPython Desktop.

Computes sky conditions for the classical planets from an observer at Thebes
(lon 33°, lat 24°) and builds Plotly line charts — no Qt dependency.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import plotly.express as px

from montu_gui.utils.debug import timed_block
from montu_gui.utils.plotly_html import figure_to_html

EPHEMERIS_PROPERTIES = [
    "RAJ2000", "DecJ2000", "RAEpoch", "DecEpoch",
    "RAGeo", "DecGeo", "el", "az", "ha", "Vmag", "rise_time", "rise_az",
    "set_time", "set_az", "transit_time", "transit_el", "elongation",
    "earth_distance", "sun_distance", "is_circumpolar", "is_neverup",
    "angsize", "phase", "hlat", "hlon", "hlong", "datestr",
]

DEFAULT_INITIAL_DATE = "-1500-01-01"
DEFAULT_TIME_SPAN = 10.0
DEFAULT_NUM_POINTS = 120
DEFAULT_PLANETS = ["Mercury"]
DEFAULT_PROPERTY = "DecEpoch"

OBSERVER_LON = 33.0
OBSERVER_LAT = 24.0


def parse_montu_date(text: str) -> tuple[str, int, int, int]:
    """Parse ``-YYYY-MM-DD`` (BCE) or ``YYYY-MM-DD`` (CE) into era + components."""
    raw = (text or DEFAULT_INITIAL_DATE).strip()
    if raw.startswith("-"):
        body = raw[1:]
        era = "bce"
    else:
        body = raw
        era = "ce"
    parts = body.split("-")
    if len(parts) < 3:
        raise ValueError(f"Invalid date: {text!r}")
    year, month, day = int(parts[0]), int(parts[1]), int(parts[2])
    return era, year, month, day


def format_montu_date(era: str, year: int, month: int, day: int) -> str:
    """Build the mixed-date string expected by ``montu.Time``."""
    if era == "bce":
        return f"-{year:04d}-{month:02d}-{day:02d}"
    return f"{year:04d}-{month:02d}-{day:02d}"


@dataclass
class PlanetsPlotResult:
    ok: bool
    html: str = ""
    title: str = ""
    n_rows: int = 0
    error: str = ""


def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def get_planet_names() -> list[str]:
    """Return display names for the classical planets (no Sun, Moon, Earth)."""
    montu = _import_montu()
    planets = [
        montu.Planet(value)
        for value in montu.PLANETARY_NAMES.values()
        if value not in ("SUN", "MOON", "EARTH")
    ]
    return [planet.name for planet in planets]


def compute_ephemerides(
    initial: str = DEFAULT_INITIAL_DATE,
    timespan: float = DEFAULT_TIME_SPAN,
    numpoints: int = DEFAULT_NUM_POINTS,
) -> pd.DataFrame:
    """Sample ephemerides for all classical planets over a time span."""
    montu = _import_montu()
    planets = [
        montu.Planet(value)
        for value in montu.PLANETARY_NAMES.values()
        if value not in ("SUN", "MOON", "EARTH")
    ]

    mtime = montu.Time(initial)
    observer = montu.Observer(lon=OBSERVER_LON, lat=OBSERVER_LAT)

    mts = []
    dates = []
    span = float(timespan or 1) * montu.YEAR
    count = max(2, int(numpoints or 10))
    for dt in np.linspace(0, span, count):
        mt = (mtime + dt).get_readable()
        mts.append(mt)
        dates.append(f"{mt.readable.year}-{mt.readable.month}-{mt.readable.day}")

    planetary_ephemerides = pd.DataFrame()
    for planet in planets:
        planet.reset_store()
        for mt in mts:
            planet.conditions_in_sky(at=mt, observer=observer, store=True)
        planet.tabulate_ephemerides()
        planet.ephemerides["datestr"] = dates
        planetary_ephemerides = pd.concat(
            [planetary_ephemerides, planet.ephemerides],
            ignore_index=True,
        )
    return planetary_ephemerides


def build_planets_plot(
    initial: str = DEFAULT_INITIAL_DATE,
    timespan: float = DEFAULT_TIME_SPAN,
    numpoints: int = DEFAULT_NUM_POINTS,
    planets: list[str] | None = None,
    property: str = DEFAULT_PROPERTY,
) -> PlanetsPlotResult:
    """Compute ephemerides and return a Plotly figure as embeddable HTML."""
    selected = list(planets or DEFAULT_PLANETS)
    if not selected:
        return PlanetsPlotResult(ok=False, error="Select at least one planet.")

    prop = property or DEFAULT_PROPERTY
    if prop not in EPHEMERIS_PROPERTIES:
        return PlanetsPlotResult(ok=False, error=f"Unknown property: {prop}")

    try:
        with timed_block("planets ephemerides"):
            df = compute_ephemerides(initial, timespan, numpoints)
        montu = _import_montu()
        mtime_initial = montu.Time(initial)
        mtime_final = (mtime_initial + float(timespan) * montu.YEAR).get_readable()

        mask = df.Name.isin(selected)
        subset = df.loc[mask]
        if subset.empty:
            return PlanetsPlotResult(
                ok=False,
                error=f"No data for planet(s): {', '.join(selected)}",
            )

        fig = px.line(subset, x="datestr", y=prop, color="Name")
        title = (
            f"Property '{prop}' of {', '.join(selected)} "
            f"starting at {mtime_initial.strftime('%Y-%m-%d')} "
            f"until {mtime_final.strftime('%Y-%m-%d')}"
        )
        fig.update_layout(
            title=dict(text=title, x=0.5, xanchor="center"),
            xaxis_title="Date [Month & Year]",
            xaxis=dict(rangeslider=dict(visible=True)),
            margin=dict(l=48, r=24, t=72, b=48),
        )

        html = figure_to_html(fig)
        return PlanetsPlotResult(
            ok=True,
            html=html,
            title=title,
            n_rows=len(subset),
        )
    except Exception as exc:
        return PlanetsPlotResult(ok=False, error=str(exc))
