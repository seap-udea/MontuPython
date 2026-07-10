"""
Planetary ephemerides logic for MontuPython Desktop.

Computes sky conditions for the classical planets from a configurable observer
site and builds Plotly line charts — no Qt dependency.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px

from montu_gui.utils.debug import timed_block
from montu_gui.utils.plotly_html import figure_to_html

_PROPERTIES_FILE = Path(__file__).parent.parent / "assets" / "planet_properties.json"

_FALLBACK_PROPERTIES = [
    {"varname": v, "quantname": v, "explanation": "", "unit": "N/A"}
    for v in [
        "RAJ2000", "DecJ2000", "RAEpoch", "DecEpoch",
        "RAGeo", "DecGeo", "el", "az", "ha", "Vmag", "rise_time", "rise_az",
        "set_time", "set_az", "transit_time", "transit_el", "elongation",
        "earth_distance", "sun_distance", "is_circumpolar", "is_neverup",
        "angsize", "phase", "hlat", "hlon", "hlong", "datestr",
    ]
]


def load_planet_properties() -> list[dict]:
    """Load property metadata from ``assets/planet_properties.json``."""
    try:
        with open(_PROPERTIES_FILE, encoding="utf-8") as fh:
            data = json.load(fh)
        props = data.get("properties", [])
        return props if props else _FALLBACK_PROPERTIES
    except (FileNotFoundError, json.JSONDecodeError):
        return _FALLBACK_PROPERTIES


def get_property_catalog() -> dict[str, dict]:
    """Map ``varname`` → property record."""
    return {p["varname"]: p for p in load_planet_properties()}


EPHEMERIS_PROPERTIES = [p["varname"] for p in load_planet_properties()]

JD_TIME_PROPERTIES = frozenset({"rise_time", "set_time", "transit_time"})

DEFAULT_INITIAL_DATE = "-1500-01-01"
DEFAULT_TIME_SPAN = 10.0
DEFAULT_NUM_POINTS = 120
DEFAULT_PLANETS = ["Mercury"]
DEFAULT_PROPERTY = "DecEpoch"

# Legacy defaults (Thebes / Luxor) — callers should pass observer coordinates.
OBSERVER_LON = 32.6422
OBSERVER_LAT = 25.6967
OBSERVER_HEIGHT_KM = 0.076

VISIBLE_PLANET_NAMES = (
    "Sun",
    "Moon",
    "Mercury",
    "Venus",
    "Mars",
    "Jupiter",
    "Saturn",
)


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


def format_proleptic_date(era: str, year: int, month: int, day: int) -> str:
    """Build a proleptic Gregorian input string (historical BCE/CE years)."""
    if era == "bce":
        return f"bce {year:04d}-{month:02d}-{day:02d} 00:00:00"
    return f"{year:04d}-{month:02d}-{day:02d} 00:00:00"


def format_montu_date(era: str, year: int, month: int, day: int) -> str:
    """Alias kept for the date picker — proleptic Gregorian historical years."""
    return format_proleptic_date(era, year, month, day)


def format_datepro_label(mtime) -> str:
    """Astronomical proleptic Gregorian date for chart axes (``datepro``)."""
    mtime.get_readable()
    return mtime.readable.datepro.split()[0]


def format_human_proleptic(mtime) -> str:
    """Human proleptic Gregorian label for titles (SPICE-style ``datespice``)."""
    mtime.get_readable()
    parts = mtime.readable.datespice.split()
    if len(parts) >= 3:
        return f"{parts[0]} {parts[1]} {parts[2]}"
    return mtime.readable.datespice


def get_property_meta(varname: str) -> dict:
    """Return metadata for a property ``varname``."""
    catalog = get_property_catalog()
    return catalog.get(varname, {
        "varname": varname,
        "quantname": varname,
        "explanation": "",
        "unit": "N/A",
    })


def is_jd_time_property(varname: str) -> bool:
    """True when ``varname`` stores an instant as Julian day (UTC)."""
    if varname in JD_TIME_PROPERTIES:
        return True
    return get_property_meta(varname).get("unit", "").upper() == "JD"


def jd_time_to_parts(jd: float) -> tuple[float, str, str]:
    """Convert a Julian-day instant to decimal hours, ``hh:mm:ss``, and event date."""
    if pd.isna(jd):
        return float("nan"), "—", "—"
    montu = _import_montu()
    readable = (
        montu.Time(float(jd), format="jd", calendar="proleptic")
        .get_readable()
        .readable
    )
    decimal_hours = (
        readable.hour
        + readable.minute / 60.0
        + readable.second / 3600.0
        + readable.usecond / 3.6e9
    )
    hms = (
        f"{readable.hour:02d}:{readable.minute:02d}:{readable.second:02d}"
    )
    event_date = readable.datepro.split()[0]
    return decimal_hours, hms, event_date


def property_axis_label(varname: str) -> str:
    """Y-axis label: ``<quantname> [unit]``."""
    meta = get_property_meta(varname)
    unit = meta.get("unit", "N/A")
    quant = meta.get("quantname", varname)
    return f"{quant} [{unit}]"


def property_help_entry() -> dict:
    """Build dynamic help text listing all properties from JSON."""
    items = []
    for prop in load_planet_properties():
        quant = prop.get("quantname", prop["varname"])
        expl = prop.get("explanation", "")
        unit = prop.get("unit", "N/A")
        varname = prop["varname"]
        items.append(
            f"<li><b>{quant}</b> (<code>{varname}</code>) — {expl} "
            f"<i>[{unit}]</i></li>"
        )
    body = (
        "<p>Select the ephemeris quantity to plot on the Y-axis. "
        "Names below match the dropdown labels.</p><ul>"
        + "".join(items)
        + "</ul>"
    )
    return {"title": "Property", "body": body}


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


def _make_body(name: str):
    """Return a ``Sun``, ``Moon``, or ``Planet`` instance by display name."""
    montu = _import_montu()
    if name == "Sun":
        return montu.Sun()
    if name == "Moon":
        return montu.Moon()
    return montu.Planet(name)


def _iter_planets():
    """Yield Sun, Moon, and classical planets offered in the desktop module."""
    for name in VISIBLE_PLANET_NAMES:
        yield _make_body(name)


def get_planet_names() -> list[str]:
    """Return display names for planets shown in the desktop module."""
    return list(VISIBLE_PLANET_NAMES)


def compute_ephemerides(
    initial: str = DEFAULT_INITIAL_DATE,
    timespan: float = DEFAULT_TIME_SPAN,
    numpoints: int = DEFAULT_NUM_POINTS,
    lon: float = OBSERVER_LON,
    lat: float = OBSERVER_LAT,
    height_km: float = OBSERVER_HEIGHT_KM,
) -> pd.DataFrame:
    """Sample ephemerides for module planets over a time span."""
    planets = list(_iter_planets())

    montu = _import_montu()

    mtime = montu.Time(initial, calendar="proleptic")
    observer = montu.Observer(lon=lon, lat=lat, height=height_km)

    mts = []
    dates = []
    span = float(timespan or 1) * montu.YEAR
    count = max(2, int(numpoints or 10))
    for dt in np.linspace(0, span, count):
        mt = (mtime + dt).get_readable()
        mts.append(mt)
        dates.append(format_datepro_label(mt))

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
    lon: float = OBSERVER_LON,
    lat: float = OBSERVER_LAT,
    height_km: float = OBSERVER_HEIGHT_KM,
    observer_name: str = "",
) -> PlanetsPlotResult:
    """Compute ephemerides and return a Plotly figure as embeddable HTML."""
    selected = list(planets or DEFAULT_PLANETS)
    if not selected:
        return PlanetsPlotResult(ok=False, error="Select at least one planet.")

    prop = property or DEFAULT_PROPERTY
    if prop not in EPHEMERIS_PROPERTIES:
        return PlanetsPlotResult(ok=False, error=f"Unknown property: {prop}")

    meta = get_property_meta(prop)
    quantname = meta.get("quantname", prop)
    y_label = property_axis_label(prop)

    try:
        with timed_block("planets ephemerides"):
            df = compute_ephemerides(
                initial, timespan, numpoints, lon=lon, lat=lat, height_km=height_km,
            )
        montu = _import_montu()
        mtime_initial = montu.Time(initial, calendar="proleptic").get_readable()
        mtime_final = (
            montu.Time(initial, calendar="proleptic")
            + float(timespan) * montu.YEAR
        ).get_readable()

        mask = df.Name.isin(selected)
        subset = df.loc[mask]
        if subset.empty:
            return PlanetsPlotResult(
                ok=False,
                error=f"No data for planet(s): {', '.join(selected)}",
            )

        plot_df = subset.copy()
        if is_jd_time_property(prop):
            hours, hms_labels, event_dates = zip(
                *(jd_time_to_parts(jd) for jd in plot_df[prop]),
                strict=True,
            )
            plot_df["_plot_y"] = hours
            plot_df["_event_hms"] = hms_labels
            plot_df["_event_date"] = event_dates
            y_col = "_plot_y"
            custom_cols = ["jed", "_event_date", "_event_hms"]
            value_hover = (
                f"{quantname}=%{{customdata[2]}} "
                f"(%{{customdata[1]}}, proleptic Gregorian)<extra></extra>"
            )
        else:
            y_col = prop
            custom_cols = ["jed"]
            value_hover = f"{quantname}=%{{y}}<extra></extra>"

        fig = px.line(
            plot_df, x="datestr", y=y_col, color="Name", custom_data=custom_cols,
        )
        fig.update_traces(
            hovertemplate=(
                "Name=%{fullData.name}<br>"
                "datestr=%{x} (proleptic Gregorian, astronomical)<br>"
                "JD=%{customdata[0]:.6f}<br>"
                + value_hover
            ),
        )
        site = observer_name or f"lat {lat:.2f}°, lon {lon:.2f}°"
        title = (
            f"{quantname} of {', '.join(selected)} "
            f"from {site} "
            f"starting at {format_human_proleptic(mtime_initial)} "
            f"until {format_human_proleptic(mtime_final)}"
        )
        fig.update_layout(
            title=dict(text=title, x=0.5, xanchor="center"),
            xaxis_title="Date (proleptic Gregorian, astronomical)",
            yaxis_title=y_label,
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
