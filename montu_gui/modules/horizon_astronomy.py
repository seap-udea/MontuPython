"""
Horizon Astronomy module for MontuPython Desktop (🌄).

Logic for rendering the horizon profile and celestial bodies near the horizon.
"""

from __future__ import annotations

from dataclasses import dataclass

import montu
import pandas as pd

from montu.maps import local_solar_to_utc_datepro
from montu_gui.utils.debug import dbg, timed_block
from montu_gui.utils.plotly_html import figure_to_html

DEFAULT_DATE = "bce 2500-06-01 00:00:00"
DEFAULT_LOCAL_HOUR = 6
DEFAULT_LOCAL_MINUTE = 0
DEFAULT_LOCAL_SECOND = 0


@dataclass
class HorizonPlotResult:
    """HTML output ready to embed in PlotlyView."""
    ok: bool
    html: str = ""
    error: str = ""
    warning: str = ""


def _wrap_horizon_html(html: str) -> str:
    """Inject styles for Plotly HTML."""
    extras = """
  <style>
    .main-svg .gtitle, .main-svg .gtitle text {
      font-size: 12px !important;
    }
  </style>
"""
    return html.replace("</head>", extras + "</head>", 1)


def build_horizon_plot(
    *,
    date_str: str = DEFAULT_DATE,
    local_hour: int = DEFAULT_LOCAL_HOUR,
    local_minute: int = DEFAULT_LOCAL_MINUTE,
    local_second: int = DEFAULT_LOCAL_SECOND,
    az_center: float = 180.0,
    az_delta: float = 180.0,
    elev_view: float = 45.0,
    show_constellations: bool = True,
    show_starnames: bool = True,
    show_asterisms: bool = False,
    show_galaxy: bool = False,
    constellation_set: str = "iau",
    bodies: list[str] | None = None,
    lat: float = 0.0,
    lon: float = 0.0,
    height_km: float = 0.0,
    location_id: str = "",
    min_height: int = 400,
    cached_horizon=None,
    calendar: str = "mixed",
) -> HorizonPlotResult:
    """Build horizon plot HTML."""
    try:
        if location_id:
            observer = montu.Observer(site=location_id)
        else:
            observer = montu.Observer(lon=lon, lat=lat, height=height_km)
            
        warning_msg = ""
        # Check if the observer has a horizon profile
        if cached_horizon is not None:
            observer.horizon = cached_horizon
        elif not hasattr(observer, "horizon") or observer.horizon is None or observer.horizon.profile is None:
            # Create a flat horizon
            if not hasattr(observer, "horizon") or observer.horizon is None:
                observer.horizon = montu.Horizon(
                    lat=float(observer.lat), 
                    lon=float(observer.lon), 
                    alt_m=float(getattr(observer, 'elevation', 0.0)),
                    observer=observer
                )
            import numpy as np
            observer.horizon.azimuths = np.linspace(0, 360, 361)
            observer.horizon.elevations = np.zeros(361)
            observer.horizon.distances = np.zeros(361)
            observer.horizon.lathorizon = np.zeros(361)
            observer.horizon.longhorizon = np.zeros(361)
            observer.horizon.profile = pd.DataFrame({
                "azimuth": observer.horizon.azimuths,
                "elevation": observer.horizon.elevations,
                "distance": observer.horizon.distances
            })
            observer.horizon._interp = lambda az: np.zeros_like(az, dtype=float) if isinstance(az, (list, np.ndarray)) else 0.0
            observer.horizon.is_computed = True
            
        # Time conversion
        obs_utc = local_solar_to_utc_datepro(
            date_str, local_hour, local_minute, local_second, lon, calendar=calendar
        )
        at = montu.Time(obs_utc, calendar="proleptic")

        with timed_block("horizon_astronomy: build horizon plot"):
            fig = observer.horizon.plot_horizon(
                at=at,
                az_center=az_center,
                az_delta=az_delta,
                elev_view=elev_view,
                show=False,
                show_planets=bodies if bodies is not None else [],
                show_title=False,
                show_stars=True,
                show_star_names=show_starnames,
                mag_limit=6.5,
                constellation_set=constellation_set,
                show_constellation_lines=show_asterisms,
                show_constellation_labels=show_constellations,
                show_constellation_boundaries=False,
                show_galaxy_equator=show_galaxy,
                show_galaxy_contours=show_galaxy,
            )

        html = _wrap_horizon_html(
            figure_to_html(
                fig,
                min_height=min_height,
                page_background="#0d1117",
            ),
        )

        return HorizonPlotResult(ok=True, html=html, warning=warning_msg)
    except Exception as exc:
        return HorizonPlotResult(ok=False, error=str(exc))
