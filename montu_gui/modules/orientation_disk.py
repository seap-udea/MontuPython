"""
Orientation Disk — extreme rise/set azimuths for celestial bodies.

For each selected body (Sun, Moon, planets, stars) the module computes,
over a 3-year span starting from a reference historical year, the
northernmost and southernmost azimuths at which the body rises (East
hemisphere) and sets (West hemisphere).

Results are displayed as an azimuth disk: a polar chart with N at the
top and colored arrows for each body's extreme orientations.

No Qt dependency — pure computation + Plotly output.
"""

from __future__ import annotations

import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dataclasses import dataclass, field
from typing import Literal

from montu_gui.utils.debug import timed_block

# ── defaults ─────────────────────────────────────────────────────────────────
DEFAULT_YEAR       = 2560
DEFAULT_ERA        = "bce"
DEFAULT_SPAN_YEARS = 3          # how many years to sweep for extremes
DEFAULT_STEP_DAYS  = 5          # sampling cadence in days
DEFAULT_MAG_LIMIT  = 1.0        # star catalogue filter for the UI dropdown

# Bodies available as planets in montu / pyephem
SOLAR_SYSTEM_BODIES = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn"]

BODY_EMOJIS: dict[str, str] = {
    "Sun": "☀",
    "Moon": "🌙",
    "Mercury": "☿",
    "Venus": "♀",
    "Mars": "♂",
    "Jupiter": "♃",
    "Saturn": "♄",
}

# Dark saturated palette — blues, greens, reds
BODY_DEFAULT_COLORS: dict[str, str] = {
    "Sun":     "#B71C1C",   # dark red
    "Moon":    "#1565C0",   # dark blue
    "Mercury": "#2E7D32",   # dark green
    "Venus":   "#C62828",   # crimson
    "Mars":    "#8E0000",   # deep red
    "Jupiter": "#1B5E20",   # forest green
    "Saturn":  "#0D47A1",   # navy blue
}

STAR_COLOR_CYCLE = [
    "#1565C0", "#1B5E20", "#B71C1C",   # blue, green, red
    "#0D47A1", "#2E7D32", "#C62828",
    "#1976D2", "#388E3C", "#D32F2F",
    "#283593",
]

DEFAULT_BODY_COLOR = "#1565C0"


# ── dataclasses ───────────────────────────────────────────────────────────────

@dataclass
class BodyConfig:
    """Configuration for a single body on the orientation disk."""
    name: str                             # display name ("Sun", "Sirius", …)
    body_type: Literal["planet", "star"]  # "planet" covers Sun/Moon/planets
    horizon_el: float = 0.0              # effective horizon altitude [degrees]
    color: str = "#FFD700"
    hip: int | None = None               # HIP catalogue number (stars only)


@dataclass
class BodyExtreme:
    """Extreme rise/set azimuths found for one body."""
    name: str
    color: str
    rise_north: float | None = None    # northernmost rise az [0-180°]
    rise_south: float | None = None    # southernmost  rise az [0-180°]
    set_north:  float | None = None    # northernmost set  az [180-360°]
    set_south:  float | None = None    # southernmost  set  az [180-360°]
    is_circumpolar: bool = False
    is_neverup:     bool = False
    error: str = ""


@dataclass
class DiskResult:
    """Raw output from :func:`compute_disk`."""
    ok: bool
    bodies: list[BodyExtreme] = field(default_factory=list)
    year: int = DEFAULT_YEAR
    era: str  = DEFAULT_ERA
    observer_name: str = ""
    lat: float = 0.0
    error: str = ""


@dataclass
class DiskPlotResult:
    """Plotly HTML ready to embed in :class:`PlotlyView`."""
    ok: bool
    html: str = ""
    error: str = ""


# ── helpers ───────────────────────────────────────────────────────────────────

def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def _year_to_jed(year: int, era: str) -> float:
    montu = _import_montu()
    y = max(1, int(year))
    if era.lower() == "bce":
        date_str = f"bce {y:04d}-01-01 00:00:00"
    else:
        date_str = f"{y:04d}-01-01 00:00:00"
    return montu.Time(date_str, calendar="proleptic").jed


def _copy_pyephem_site(site, date=None, horizon_el: float = 0.0):
    """Clone a pyephem Observer, optionally setting date and horizon."""
    import ephem as pyephem
    s = pyephem.Observer()
    s.lon       = site.lon
    s.lat       = site.lat
    s.elevation = site.elevation
    s.pressure  = site.pressure
    s.temp      = site.temp
    s.horizon   = str(horizon_el)
    s.date      = date if date is not None else site.date
    return s


def _rise_set_from_pyephem(seba, site, horizon_el: float = 0.0):
    """Return (rise_az°, set_az°, circumpolar, neverup) for one epoch.

    *seba* is a pyephem body; *site* is a pyephem Observer already set
    to the desired date.
    """
    import ephem as pyephem

    obs = _copy_pyephem_site(site, horizon_el=horizon_el)

    rise_az   = None
    set_az    = None
    circumpolar = False
    neverup     = False

    try:
        rd = obs.next_rising(seba)
        rs = _copy_pyephem_site(site, date=rd, horizon_el=horizon_el)
        seba.compute(rs)
        rise_az = math.degrees(float(seba.az))
    except pyephem.AlwaysUpError:
        circumpolar = True
    except pyephem.NeverUpError:
        neverup = True
    except Exception:
        pass

    try:
        obs2 = _copy_pyephem_site(site, horizon_el=horizon_el)
        sd = obs2.next_setting(seba)
        ss = _copy_pyephem_site(site, date=sd, horizon_el=horizon_el)
        seba.compute(ss)
        set_az = math.degrees(float(seba.az))
    except pyephem.AlwaysUpError:
        circumpolar = True
    except pyephem.NeverUpError:
        neverup = True
    except Exception:
        pass

    # Restore seba to original epoch
    seba.compute(site)
    return rise_az, set_az, circumpolar, neverup


def _star_rise_set_az(dec_deg: float, lat_deg: float, horizon_el: float = 0.0):
    """Analytical rise/set azimuths for a fixed-position celestial object.

    Returns ``(rise_az, set_az, circumpolar, neverup)`` in degrees.
    Rise az is in [0, 180]; set az is in [180, 360].
    """
    dec = math.radians(dec_deg)
    lat = math.radians(lat_deg)
    h   = math.radians(horizon_el)

    denom = math.cos(lat) * math.cos(h)
    if abs(denom) < 1e-12:
        return None, None, True, False

    cos_az = (math.sin(dec) - math.sin(lat) * math.sin(h)) / denom

    if cos_az >= 1.0:
        return None, None, True, False   # circumpolar
    if cos_az <= -1.0:
        return None, None, False, True   # never rises above h

    rise_az = math.degrees(math.acos(cos_az))   # 0-180 (east sector)
    set_az  = 360.0 - rise_az                    # 180-360 (west sector)
    return rise_az, set_az, False, False


# ── main computation ──────────────────────────────────────────────────────────

def compute_disk(
    year: int,
    era: str,
    lat: float,
    lon: float,
    height: float,
    bodies: list[BodyConfig],
    span_years: float = DEFAULT_SPAN_YEARS,
    step_days:  float = DEFAULT_STEP_DAYS,
    observer_name: str = "",
) -> DiskResult:
    """Compute extreme rise/set azimuths for each body.

    Sweeps ``span_years`` years (default 3) from Jan 1 of the reference
    year in increments of ``step_days`` days.  Planets and the Sun/Moon
    are sampled via PyEphem; stars are computed analytically from their
    precessed declination at the start epoch.

    Parameters
    ----------
    year, era : int, str
        Reference year (historical count, positive integer) and era
        (``"bce"`` or ``"ce"``).
    lat, lon, height : float
        Observer coordinates [degrees, degrees, km].
    bodies : list of BodyConfig
        Bodies to include.  Must not be empty.
    span_years : float
        Length of the search window [years].
    step_days : float
        Sampling cadence [days].
    observer_name : str
        Label used in plot titles.

    Returns
    -------
    DiskResult
    """
    if not bodies:
        return DiskResult(ok=False, error="No bodies selected.")

    try:
        montu = _import_montu()

        jed_start = _year_to_jed(year, era)
        jed_end   = jed_start + span_years * 365.25
        jed_steps = np.arange(jed_start, jed_end, float(step_days))

        # Build montu observer for PyEphem pass-through
        obs = montu.Observer(lon=lon, lat=lat, height=height / 1000.0)  # montu uses km

        result_bodies: list[BodyExtreme] = []

        for cfg in bodies:
            if cfg.body_type == "star":
                extreme = _process_star(cfg, jed_start, lat)
            else:
                extreme = _process_planet(cfg, jed_steps, obs, lat)
            result_bodies.append(extreme)

        return DiskResult(
            ok=True,
            bodies=result_bodies,
            year=year,
            era=era,
            observer_name=observer_name,
            lat=lat,
        )

    except Exception as exc:
        return DiskResult(ok=False, error=str(exc))


def _process_planet(
    cfg: BodyConfig,
    jed_steps: np.ndarray,
    obs,
    lat: float,
) -> BodyExtreme:
    """Sample rise/set azimuths for a solar-system body across all JED steps."""
    try:
        montu = _import_montu()

        # Instantiate the montu body
        name_lower = cfg.name.lower()
        if name_lower == "sun":
            body = montu.Sun()
        elif name_lower == "moon":
            body = montu.Moon()
        else:
            body = montu.Planet(cfg.name)

        rise_azs: list[float] = []
        set_azs:  list[float] = []
        circumpolar = False
        neverup     = False

        for jed in jed_steps:
            try:
                # Update pyephem date
                obs.site.date = float(jed) - montu.PYEPHEM_JD_REF
                body.seba.compute(obs.site)

                r_az, s_az, cp, nu = _rise_set_from_pyephem(
                    body.seba, obs.site, horizon_el=cfg.horizon_el
                )
                if cp:
                    circumpolar = True
                if nu:
                    neverup = True
                if r_az is not None:
                    rise_azs.append(r_az)
                if s_az is not None:
                    set_azs.append(s_az)
            except Exception:
                continue

        if not rise_azs and not set_azs:
            return BodyExtreme(
                name=cfg.name, color=cfg.color,
                is_circumpolar=circumpolar,
                is_neverup=neverup,
            )

        # Northernmost rise = smallest az (closest to 0°)
        # Southernmost rise = largest az (closest to 180°)
        # Northernmost set = largest az (closest to 360°)
        # Southernmost set = smallest az (closest to 180°)
        return BodyExtreme(
            name=cfg.name,
            color=cfg.color,
            rise_north=min(rise_azs) if rise_azs else None,
            rise_south=max(rise_azs) if rise_azs else None,
            set_north =max(set_azs)  if set_azs  else None,
            set_south =min(set_azs)  if set_azs  else None,
            is_circumpolar=circumpolar,
            is_neverup=neverup,
        )

    except Exception as exc:
        return BodyExtreme(name=cfg.name, color=cfg.color, error=str(exc))


def _process_star(cfg: BodyConfig, jed_start: float, lat: float) -> BodyExtreme:
    """Compute rise/set azimuths for a star analytically via precessed Dec."""
    try:
        montu = _import_montu()

        t_start = montu.Time(jed_start, format="jd", calendar="proleptic")

        if cfg.hip is not None:
            stars = montu.Stars()
            subset = stars.get_stars(HIP=cfg.hip)
        else:
            stars = montu.Stars()
            subset = stars.get_stars(ProperName=cfg.name)
            if subset.data.empty:
                subset = stars.get_stars(Name=cfg.name)

        if subset.data.empty:
            return BodyExtreme(
                name=cfg.name, color=cfg.color,
                error=f"Star '{cfg.name}' not found in catalogue.",
            )

        precessed = subset.where_in_space(at=t_start)
        dec_epoch = float(precessed.data["DecEpoch"].iloc[0])

        r_az, s_az, cp, nu = _star_rise_set_az(dec_epoch, lat, cfg.horizon_el)

        return BodyExtreme(
            name=cfg.name,
            color=cfg.color,
            rise_north=r_az,   # stars have a single fixed azimuth (N==S extreme)
            rise_south=r_az,
            set_north =s_az,
            set_south =s_az,
            is_circumpolar=cp,
            is_neverup=nu,
        )

    except Exception as exc:
        return BodyExtreme(name=cfg.name, color=cfg.color, error=str(exc))


# ── star catalogue helper ──────────────────────────────────────────────────────

def get_available_stars(mag_limit: float = DEFAULT_MAG_LIMIT) -> pd.DataFrame:
    """Return catalogue stars brighter than *mag_limit* sorted by Vmag.

    Columns: ``ProperName``, ``Name``, ``Vmag``, ``HIP``.
    Rows without a ProperName are excluded — only named stars are offered
    in the UI dropdown.
    """
    try:
        montu = _import_montu()
        stars = montu.Stars()
        bright = stars.get_stars(Vmag=[-2.0, float(mag_limit)])
        df = bright.data[["ProperName", "Name", "Vmag", "HIP"]].copy()
        df = df[
            df["ProperName"].notna()
            & (df["ProperName"].astype(str).str.strip() != "")
            & (~df["ProperName"].astype(str).isin(["nan", "None"]))
        ]
        df = df.sort_values("Vmag").reset_index(drop=True)
        return df
    except Exception:
        return pd.DataFrame(columns=["ProperName", "Name", "Vmag", "HIP"])


# ── Plotly disk ───────────────────────────────────────────────────────────────

def build_disk_plot(result: DiskResult, *, plot_height: int = 640) -> DiskPlotResult:
    """Build the orientation disk as embeddable Plotly HTML.

    Parameters
    ----------
    result : DiskResult
        Output from :func:`compute_disk`.
    plot_height : int
        Pixel height for the embedded Plotly figure.

    Returns
    -------
    DiskPlotResult
    """
    if not result.ok:
        return DiskPlotResult(ok=False, error=result.error)
    try:
        fig = _build_disk_figure(result, plot_height=plot_height)
        return DiskPlotResult(ok=True, html=_figure_to_disk_html(fig, plot_height))
    except Exception as exc:
        return DiskPlotResult(ok=False, error=str(exc))


def _figure_to_disk_html(fig: go.Figure, height_px: int) -> str:
    """Embed Plotly HTML sized to fill the available panel height."""
    from montu_gui.utils.plotly_html import plotly_js_path, _PLOTLY_RESIZE_SCRIPT

    div = fig.to_html(include_plotlyjs=False, full_html=False, config={
        "responsive": True,
        "displayModeBar": True,
    })
    js_url = plotly_js_path().as_uri()
    h = max(320, int(height_px))
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <script src="{js_url}"></script>
  <style>
    html, body {{
      margin: 0;
      padding: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      background: #ffffff;
    }}
    .plotly-graph-div {{
      width: 100% !important;
      min-height: {h}px !important;
      height: {h}px !important;
    }}
    .js-plotly-plot .plotly .modebar {{
      top: auto !important;
      left: auto !important;
      right: 12px !important;
      bottom: 12px !important;
    }}
    .js-plotly-plot .plotly .modebar-group {{
      background: rgba(255, 255, 255, 0.85) !important;
      border-radius: 4px;
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.12);
    }}
  </style>
</head>
<body>
{div}
{_PLOTLY_RESIZE_SCRIPT}
</body>
</html>"""


def _az_label(name: str, suffix: str) -> str:
    """Short label for a body (up to 3 chars) + suffix."""
    abbrev = name[:3] if len(name) > 3 else name
    return f"{abbrev}{suffix}"


def _build_disk_figure(result: DiskResult, *, plot_height: int = 640) -> go.Figure:
    """Internal: assemble the Plotly polar figure."""

    fig = go.Figure()

    # ── disk background ───────────────────────────────────────────────────
    # Faint ring
    theta_ring = list(range(361))
    fig.add_trace(go.Scatterpolar(
        r=[1.0] * 361,
        theta=theta_ring,
        mode="lines",
        line=dict(color="rgba(80,80,80,0.4)", width=1.5),
        hoverinfo="skip",
        showlegend=False,
    ))

    # Inner reference rings at 0.5 (subtle)
    fig.add_trace(go.Scatterpolar(
        r=[0.5] * 361,
        theta=theta_ring,
        mode="lines",
        line=dict(color="rgba(180,180,180,0.3)", width=0.8, dash="dot"),
        hoverinfo="skip",
        showlegend=False,
    ))

    # Cardinal cross lines (N-S and E-W)
    for az_pair in [(0, 180), (90, 270)]:
        fig.add_trace(go.Scatterpolar(
            r=[0, 1.1, None, 0, 1.1],
            theta=[az_pair[0], az_pair[0], None, az_pair[1], az_pair[1]],
            mode="lines",
            line=dict(color="rgba(100,100,100,0.35)", width=1),
            hoverinfo="skip",
            showlegend=False,
        ))

    # ── body arrows ───────────────────────────────────────────────────────
    has_legend = False
    for body in result.bodies:
        if body.error:
            continue
        if body.is_circumpolar or body.is_neverup:
            continue

        color = body.color
        name  = body.name

        # Collect the four extremes to draw
        arrows: list[tuple[float, bool, str]] = []  # (az, is_rise, label)
        if body.rise_north is not None:
            arrows.append((body.rise_north, True,  "△N"))
        if body.rise_south is not None and body.rise_south != body.rise_north:
            arrows.append((body.rise_south, True,  "△S"))
        elif body.rise_south is not None and body.rise_south == body.rise_north:
            # Single rise point (star): only one marker but label as "△"
            pass  # already added above as △N but rename
        if body.set_north is not None:
            arrows.append((body.set_north,  False, "▽N"))
        if body.set_south is not None and body.set_south != body.set_north:
            arrows.append((body.set_south,  False, "▽S"))

        # For stars (rise_north == rise_south), relabel
        if body.rise_north == body.rise_south and body.rise_north is not None:
            arrows = [
                (body.rise_north, True,  "△"),
                (body.set_north,  False, "▽"),
            ]

        first = True
        for az, is_rise, suffix in arrows:
            symbol = "triangle-up" if is_rise else "triangle-down"
            label_text = f"<b>{_az_label(name, suffix)}</b>"

            # Line from center to 0.90
            fig.add_trace(go.Scatterpolar(
                r=[0.0, 0.90],
                theta=[az, az],
                mode="lines",
                line=dict(color=color, width=2.2),
                showlegend=first and not has_legend,
                name=name,
                hoverinfo="skip",
            ))

            # Arrowhead marker at 0.95
            fig.add_trace(go.Scatterpolar(
                r=[0.95],
                theta=[az],
                mode="markers",
                marker=dict(
                    symbol=symbol,
                    size=12,
                    color=color,
                    line=dict(width=1.2, color="white"),
                ),
                showlegend=first,
                name=name if first else "",
                legendgroup=name,
                hovertemplate=(
                    f"<b>{name}</b><br>"
                    f"{'Rise' if is_rise else 'Set'} az: {az:.1f}°"
                    "<extra></extra>"
                ),
            ))

            # Label at 1.12
            fig.add_trace(go.Scatterpolar(
                r=[1.12],
                theta=[az],
                mode="text",
                text=[label_text],
                textfont=dict(size=9, color=color),
                showlegend=False,
                hoverinfo="skip",
            ))

            if first:
                has_legend = True
                first = False

    # ── compass labels ────────────────────────────────────────────────────
    compass = [
        (0,   "<b>N</b>"), (45,  "NE"), (90,  "<b>E</b>"),
        (135, "SE"),       (180, "<b>S</b>"), (225, "SW"),
        (270, "<b>W</b>"), (315, "NW"),
    ]
    for az, lbl in compass:
        fig.add_trace(go.Scatterpolar(
            r=[1.25],
            theta=[az],
            mode="text",
            text=[lbl],
            textfont=dict(
                size=15 if "<b>" in lbl else 11,
                color="rgba(60,60,60,0.85)",
            ),
            showlegend=False,
            hoverinfo="skip",
        ))

    # ── title & layout ────────────────────────────────────────────────────
    era_label = f"{result.year} BCE" if result.era.lower() == "bce" else f"{result.year} CE"
    title_text = (
        f"⭕  Orientation Disk  ·  {era_label}  ·  "
        f"{result.observer_name} (lat {result.lat:.2f}°)"
        f"<br><sup>Arrows show extreme rise (△) and set (▽) azimuths  ·  "
        f"N = northernmost  ·  S = southernmost</sup>"
    )

    fig.update_layout(
        polar=dict(
            angularaxis=dict(
                direction="clockwise",
                rotation=90,            # N at top
                tickmode="array",
                tickvals=list(range(0, 360, 10)),
                ticktext=[
                    f"{d}°" if d % 30 == 0 else ""
                    for d in range(0, 360, 10)
                ],
                tickfont=dict(size=9, color="rgba(80,80,80,0.8)"),
                linecolor="rgba(80,80,80,0.4)",
                gridcolor="rgba(180,180,180,0.3)",
                linewidth=1.5,
                gridwidth=0.8,
                ticks="outside",
                ticklen=6,
                showline=True,
            ),
            radialaxis=dict(
                visible=False,
                range=[0, 1.4],
            ),
            bgcolor="rgba(240, 245, 255, 0.4)",
            domain=dict(x=[0.05, 0.95], y=[0.0, 0.95]),
        ),
        title=dict(
            text=title_text,
            x=0.5,
            xanchor="center",
            font=dict(size=13),
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.08,
            xanchor="center",
            x=0.5,
            font=dict(size=12),
            bgcolor="rgba(255,255,255,0.7)",
            bordercolor="rgba(180,180,180,0.5)",
            borderwidth=1,
        ),
        margin=dict(l=60, r=60, t=110, b=80),
        height=max(320, int(plot_height)),
        autosize=True,
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig
