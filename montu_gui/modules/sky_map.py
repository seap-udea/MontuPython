"""
Sky map module for MontuPython Desktop (🌌).

Thin wrapper around :func:`montu.maps.polar_sky_map` with catalogue caching and
Plotly HTML embedding for the Qt UI.
"""

from __future__ import annotations

from dataclasses import dataclass

import montu
import pandas as pd

from montu.maps import (
    DEFAULT_CONSTELLATION_SET,
    DEFAULT_MAG_LIMIT,
    LINE_ECLIPTIC,
    LINE_HORIZON,
    local_solar_to_utc_datepro,
    polar_sky_map,
)
from montu_gui.utils.debug import dbg, timed_block
from montu_gui.utils.plotly_html import figure_to_html

DEFAULT_DATE = "bce 2500-06-01 00:00:00"
DEFAULT_LOCAL_HOUR = 18
DEFAULT_LOCAL_MINUTE = 0
DEFAULT_LOCAL_SECOND = 0
DEFAULT_BODIES = ["Sun", "Moon", "Mercury", "Venus"]
DEFAULT_LINES = ["Ecliptic", "Horizon"]
CONSTELLATION_SETS: tuple[tuple[str, str], ...] = (
    ("iau", "IAU Constellations"),
    ("egyptian_ancient", "Ancient Egyptian"),
    ("egyptian_dendera", "Dendera Egyptian"),
)
PRECESS_THRESHOLD_YEARS = 100.0

_CATALOG: pd.DataFrame | None = None
_PRECESSED_DATE: str | None = None
_PRECESSED_JED: float | None = None
_PRECESSED_DATA: pd.DataFrame | None = None


@dataclass
class SkyMapPlotResult:
    """HTML output ready to embed in PlotlyView (one page per hemisphere)."""

    ok: bool
    html_north: str = ""
    html_south: str = ""
    n_north: int = 0
    n_south: int = 0
    error: str = ""


def _calendar_date_only(date_str: str) -> str:
    tokens = date_str.strip().split()
    if tokens and tokens[0].lower() == "bce":
        return f"bce {tokens[1][:10]}"
    return tokens[0][:10]


def _precession_epoch_string(calendar_date: str) -> str:
    return f"{_calendar_date_only(calendar_date)} 12:00:00"


def _date_jed(date_str: str) -> float:
    return float(montu.Time(date_str, calendar="proleptic").jed)


def _years_apart(jed_a: float, jed_b: float) -> float:
    return abs(jed_a - jed_b) / 365.25


def _load_catalog() -> pd.DataFrame:
    global _CATALOG
    if _CATALOG is not None:
        dbg("sky_map: catalogue already in memory (skip reload)")
        return _CATALOG

    with timed_block("sky_map: load stellar catalogue"):
        _CATALOG = montu.Stars().data.copy()
    return _CATALOG


def _get_precessed(calendar_date: str) -> pd.DataFrame:
    """Return precessed catalogue rows, reusing cache within threshold."""
    global _PRECESSED_DATE, _PRECESSED_JED, _PRECESSED_DATA

    epoch_str = _precession_epoch_string(calendar_date)
    jed = _date_jed(epoch_str)

    if _PRECESSED_DATA is not None and _PRECESSED_JED is not None:
        delta_years = _years_apart(jed, _PRECESSED_JED)
        if delta_years <= PRECESS_THRESHOLD_YEARS:
            dbg(
                f"sky_map: reuse precession from {_PRECESSED_DATE!r} "
                f"(Δt={delta_years:.1f} yr ≤ {PRECESS_THRESHOLD_YEARS:.0f} yr, "
                f"requested {epoch_str!r})"
            )
            return _PRECESSED_DATA

    with timed_block(f"sky_map: precess stars → {epoch_str}"):
        catalog = _load_catalog()
        mtime = montu.Time(epoch_str, calendar="proleptic")
        precessed = montu.Stars(data=catalog).where_in_space(at=mtime).data.copy()
        precessed["ra_deg"] = precessed["RAEpoch"] * 15.0
        precessed["display_name"] = precessed.apply(_star_name, axis=1)

    _PRECESSED_DATE = epoch_str
    _PRECESSED_JED = jed
    _PRECESSED_DATA = precessed
    return _PRECESSED_DATA


def _star_name(row: pd.Series) -> str:
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


def clear_sky_map_cache() -> None:
    """Drop cached catalogue and precessed tables (tests / manual reset)."""
    global _CATALOG, _PRECESSED_DATE, _PRECESSED_JED, _PRECESSED_DATA
    _CATALOG = None
    _PRECESSED_DATE = None
    _PRECESSED_JED = None
    _PRECESSED_DATA = None


def _wrap_sky_map_html(html: str) -> str:
    """Inject styles for compact title and emoji body legends on dark UI."""
    extras = """
  <style>
    .main-svg .gtitle, .main-svg .gtitle text {
      font-size: 11px !important;
    }
    .main-svg .legend .legendpoints path,
    .main-svg .legend .legendpoints circle {
      display: none !important;
    }
    .main-svg .legend .legendtext {
      font-size: 12px !important;
    }
  </style>
"""
    return html.replace("</head>", extras + "</head>", 1)


def build_sky_map_plot(
    *,
    date_str: str = DEFAULT_DATE,
    mag_limit: float = DEFAULT_MAG_LIMIT,
    local_hour: int = DEFAULT_LOCAL_HOUR,
    local_minute: int = DEFAULT_LOCAL_MINUTE,
    local_second: int = DEFAULT_LOCAL_SECOND,
    bodies: list[str] | None = None,
    lines: list[str] | None = None,
    meridian_view: bool = False,
    constellation_set: str = DEFAULT_CONSTELLATION_SET,
    observer_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    height_km: float = 0.0,
    min_height: int = 420,
) -> SkyMapPlotResult:
    """Build sky map HTML with precessed stars and asterism lines."""
    try:
        local_time = (
            f"{local_hour:02d}:{local_minute:02d}:{local_second:02d} local"
        )
        obs_utc = local_solar_to_utc_datepro(
            date_str, local_hour, local_minute, local_second, lon,
        )
        dbg(f"sky_map: observation {local_time} → UTC {obs_utc}")
        dbg(f"sky_map: bodies on map: {bodies if bodies is not None else DEFAULT_BODIES}")
        dbg(f"sky_map: lines on map: {lines if lines is not None else DEFAULT_LINES}")
        dbg(f"sky_map: meridian view: {meridian_view}")
        dbg(f"sky_map: constellation set: {constellation_set}")

        observer = montu.Observer(lon=lon, lat=lat, height=height_km)
        precessed = _get_precessed(date_str)

        with timed_block("sky_map: build polar figures"):
            fig_north, fig_south = polar_sky_map(
                date_str,
                local_hour=local_hour,
                local_minute=local_minute,
                local_second=local_second,
                observer=observer,
                mag_limit=float(mag_limit),
                bodies=bodies,
                lines=lines,
                meridian_view=meridian_view,
                constellation_set=constellation_set,
                observer_name=observer_name,
                precessed_star_data=precessed,
            )

        stars = precessed[
            (precessed["Vmag"] <= float(mag_limit))
            & precessed["RAEpoch"].notna()
            & precessed["DecEpoch"].notna()
        ]
        n_north = int((stars["DecEpoch"] >= 0.0).sum())
        n_south = int((stars["DecEpoch"] <= 0.0).sum())

        html_north = _wrap_sky_map_html(
            figure_to_html(
                fig_north,
                min_height=min_height,
                page_background="#0d1117",
            ),
        )
        html_south = _wrap_sky_map_html(
            figure_to_html(
                fig_south,
                min_height=min_height,
                page_background="#0d1117",
            ),
        )

        return SkyMapPlotResult(
            ok=True,
            html_north=html_north,
            html_south=html_south,
            n_north=n_north,
            n_south=n_south,
        )
    except Exception as exc:
        return SkyMapPlotResult(ok=False, error=str(exc))
