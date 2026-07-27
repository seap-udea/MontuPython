"""
Sky map plotting for MontuPython.

Mercator (equatorial) and polar (azimuthal) Plotly sky maps with constellation
asterisms, optional horizon and ecliptic overlays, and solar-system bodies.
"""

from __future__ import annotations

from functools import lru_cache

import montu
import numpy as np
import pandas as pd
import plotly.graph_objects as go

from montu.stars import (
    CONSTELLATION_SET_IDS,
    STELLAR_CATALOGUE,
    parse_constellation_boundaries,
    parse_constellation_lines,
    parse_constellation_names,
)


# ── Mercator (equatorial) sky map ─────────────────────────────────────────

def _import_plotly():
    try:
        import plotly.graph_objects as go
        return go
    except ImportError as exc:
        raise ImportError(
            "Plotly is required for sky maps. "
            "Install with: pip install plotly"
        ) from exc


def _mag_to_marker_size(vmag: float) -> float:
    return float(np.clip(13.0 - 2.0 * vmag, 3.0, 22.0))


def ra_deg_about_center(ra_deg, center_ra_deg):
    """Express right ascension on a branch continuous about *center_ra_deg*."""
    delta = (float(ra_deg) - float(center_ra_deg) + 180.0) % 360.0 - 180.0
    return float(center_ra_deg) + delta


def _angular_separation_ra_dec_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg):
    """Great-circle separation [deg] between two equatorial points."""
    return float(np.rad2deg(
        montu.Util.haversine_distance(
            np.deg2rad(dec1_deg), np.deg2rad(ra1_deg),
            np.deg2rad(dec2_deg), np.deg2rad(ra2_deg),
        )
    ))


def unwrap_figure_ra_deg(fig, center_ra_deg):
    """Re-express every trace ``x`` coordinate about *center_ra_deg* in place."""
    for trace in fig.data:
        if trace.x is None:
            continue
        trace.x = [
            ra_deg_about_center(x, center_ra_deg)
            if x is not None and not (isinstance(x, float) and np.isnan(x))
            else x
            for x in trace.x
        ]


def _star_display_name(row) -> str:
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


def _build_hip_lookup(
    star_data: pd.DataFrame,
    ra_col: str = "RAEpoch",
    dec_col: str = "DecEpoch",
) -> dict:
    lookup = {}
    for _, row in star_data.iterrows():
        hip = row.get("HIP", np.nan)
        if pd.isna(hip):
            continue
        ra = row.get(ra_col, np.nan)
        dec = row.get(dec_col, np.nan)
        if pd.isna(ra) or pd.isna(dec):
            continue
        ra_deg = float(ra) * 15.0
        lookup[int(hip)] = (ra_deg, float(dec))
    return lookup


def _fab_hip_ids(path=None, *, set_id: str = 'iau'):
    """HIP catalogue numbers referenced in a stick-figure file."""
    hips = set()
    for entry in parse_constellation_lines(path=path, set_id=set_id):
        for hip_a, hip_b in entry["segments"]:
            hips.add(hip_a)
            hips.add(hip_b)
    return hips


def _complete_hip_lookup(
    star_data: pd.DataFrame,
    ra_col: str = "RAEpoch",
    dec_col: str = "DecEpoch",
    at=None,
) -> dict:
    """Build HIP→(RA°, Dec°) lookup, supplementing asterism stars from the catalogue."""
    lookup = _build_hip_lookup(star_data, ra_col=ra_col, dec_col=dec_col)
    missing = _fab_hip_ids() - set(lookup.keys())
    if not missing:
        return lookup
    cat = pd.read_csv(
        montu.Util._data_path(STELLAR_CATALOGUE, check=True),
        low_memory=False,
    )
    subset = cat[cat["HIP"].isin(list(missing))].copy()
    if subset.empty:
        return lookup
    use_epoch = ra_col in ("RAEpoch",) and at is not None
    if use_epoch:
        from montu.stars import Stars
        extra = Stars(data=subset)
        extra = extra.where_in_space(at=at)
        lookup.update(_build_hip_lookup(extra.data, "RAEpoch", "DecEpoch"))
    else:
        lookup.update(_build_hip_lookup(subset, "RAJ2000", "DecJ2000"))
    return lookup


def _polyline_ra_dec(points, *, split_wrap=True):
    """Expand polygon vertices into x/y lists with ``None`` breaks at RA wraps."""
    if not points:
        return [], []
    xs, ys = [], []
    for k, (ra, dec) in enumerate(points):
        if k > 0 and split_wrap:
            ra_prev = points[k - 1][0]
            if abs(ra - ra_prev) > 180.0:
                xs.append(None)
                ys.append(None)
        xs.append(ra)
        ys.append(dec)
    return xs, ys


def mercator_sky_axes():
    """Default Plotly axis dicts for an equatorial Mercator sky map."""
    xaxis = dict(
        title="Right Ascension [h]",
        autorange="reversed",
        range=[360, 0],
        gridcolor="#1a2740",
        tickvals=list(range(0, 361, 30)),
        ticktext=[f"{v // 15}h" for v in range(0, 361, 30)],
        color="#8899aa",
        showgrid=True,
        zeroline=False,
    )
    yaxis = dict(
        title="Declination [°]",
        range=[-90, 90],
        gridcolor="#1a2740",
        tickvals=list(range(-90, 91, 30)),
        color="#8899aa",
        showgrid=True,
        zeroline=False,
    )
    return xaxis, yaxis


def mercator_sky_map(
    star_data: pd.DataFrame,
    *,
    ra_col: str = "RAEpoch",
    dec_col: str = "DecEpoch",
    mag_col: str = "Vmag",
    mag_limit: float = 6.5,
    label_bright_mag: float = 2.5,
    show_stars: bool = True,
    show_constellation_lines: bool = True,
    show_constellation_boundaries: bool = True,
    show_constellation_labels: bool = True,
    constellation_full_names: bool = False,
    label_bounds: tuple | None = None,
    label_center: tuple | None = None,
    label_radius_deg: float | None = None,
    constellation_label_font: dict | None = None,
    at=None,
    layout=None,
):
    """Build a base equatorial Mercator sky map (Plotly).

    Draws IAU constellation boundaries, asterism lines, constellation
    names (abbreviations or full names), and background stars.  Alignment
    overlays (target declination, circumpolar limit, highlighted stars,
    title, etc.) should be added by the caller on the returned figure.

    Parameters
    ----------
    star_data : pandas.DataFrame
        Stellar catalogue rows at the map epoch (must include ``HIP`` for
        constellation geometry).
    ra_col, dec_col : str
        Right ascension [hours] and declination [degrees] column names.
    mag_limit : float
        Faint limit for background stars.
    label_bright_mag : float
        Annotate stars brighter than this V magnitude.
    constellation_full_names : bool, optional
        When ``True``, label constellations with their full names (e.g.
        ``Taurus``) instead of three-letter abbreviations.
    label_bounds : tuple of float, optional
        ``(ra_min_deg, ra_max_deg, dec_min_deg, dec_max_deg)`` window used
        to keep only constellation labels whose centroid falls inside the
        field of view. Prefer :paramref:`label_center` and
        :paramref:`label_radius_deg` near RA = 0 h / 24 h.
    label_center : tuple of float, optional
        ``(ra_deg, dec_deg)`` used with :paramref:`label_radius_deg` to
        filter labels by angular distance (handles RA wrap).
    label_radius_deg : float, optional
        Angular radius [deg] for :paramref:`label_center` label filtering.
    constellation_label_font : dict, optional
        Plotly ``textfont`` dict for constellation labels.
    at : montu.Time, optional
        Epoch for precessing asterism stars missing from *star_data*.
    layout : dict, optional
        Extra keys merged into ``fig.update_layout`` (no title by default).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    go = _import_plotly()
    fig = go.Figure()
    hip_lookup = _complete_hip_lookup(
        star_data, ra_col=ra_col, dec_col=dec_col, at=at,
    )

    if show_constellation_boundaries:
        bx, by = [], []
        for poly in parse_constellation_boundaries(at=at):
            px, py = _polyline_ra_dec(poly["points"])
            if px:
                bx.extend(px + [None])
                by.extend(py + [None])
        if bx:
            fig.add_trace(go.Scatter(
                x=bx, y=by, mode="lines",
                line=dict(color="rgba(230, 120, 170, 0.65)", width=0.8),
                hoverinfo="skip", showlegend=False, name="boundaries",
            ))

    label_positions: dict[str, list[tuple[float, float]]] = {}

    if show_constellation_lines:
        lx, ly = [], []
        for entry in parse_constellation_lines():
            abbrev = entry["abbrev"]
            for hip_a, hip_b in entry["segments"]:
                pa = hip_lookup.get(hip_a)
                pb = hip_lookup.get(hip_b)
                if pa is None or pb is None:
                    continue
                ra1, dec1 = pa
                ra2, dec2 = pb
                if abs(ra2 - ra1) > 180.0:
                    lx.extend([ra1, None])
                    ly.extend([dec1, None])
                lx.extend([ra1, ra2, None])
                ly.extend([dec1, dec2, None])
                label_positions.setdefault(abbrev, []).append(pa)
                label_positions.setdefault(abbrev, []).append(pb)
        if lx:
            fig.add_trace(go.Scatter(
                x=lx, y=ly, mode="lines",
                line=dict(color="rgba(110, 125, 145, 0.55)", width=1.0),
                hoverinfo="skip", showlegend=False, name="asterisms",
            ))

    if show_constellation_labels and label_positions:
        name_labels = (
            parse_constellation_names() if constellation_full_names else {}
        )
        default_label_font = dict(
            size=9, color='rgba(130, 140, 155, 0.42)',
        )
        if constellation_label_font:
            default_label_font.update(constellation_label_font)
        label_x, label_y, label_text = [], [], []
        for abbrev, coords in label_positions.items():
            ra_mean = float(np.mean([c[0] for c in coords]))
            dec_mean = float(np.mean([c[1] for c in coords]))
            if label_center is not None and label_radius_deg is not None:
                center_ra, center_dec = label_center
                if _angular_separation_ra_dec_deg(
                    center_ra, center_dec, ra_mean, dec_mean,
                ) > float(label_radius_deg):
                    continue
            elif label_bounds is not None:
                ra_min, ra_max, dec_min, dec_max = label_bounds
                if not (
                    ra_min <= ra_mean <= ra_max
                    and dec_min <= dec_mean <= dec_max
                ):
                    continue
            label_x.append(ra_mean)
            label_y.append(dec_mean)
            label_text.append(name_labels.get(abbrev, abbrev))
        if label_x:
            fig.add_trace(go.Scatter(
                x=label_x, y=label_y, mode='text', text=label_text,
                textfont=default_label_font,
                hoverinfo='skip', showlegend=False, name='constellation labels',
            ))

    if show_stars and not star_data.empty:
        data = star_data[star_data[mag_col] <= float(mag_limit)].copy()
        if not data.empty:
            if ra_col.endswith("J2000") or ra_col == "RAJ2000":
                data["ra_deg"] = data[ra_col] * 15.0
            elif ra_col.startswith("RA"):
                data["ra_deg"] = data[ra_col] * 15.0
            else:
                data["ra_deg"] = data[ra_col]
            data["msize"] = data[mag_col].apply(_mag_to_marker_size)
            data["display_name"] = data.apply(_star_display_name, axis=1)
            label_col = data.apply(
                lambda r: r["display_name"] if r[mag_col] <= label_bright_mag else "",
                axis=1,
            )
            fig.add_trace(go.Scatter(
                x=data["ra_deg"],
                y=data[dec_col],
                mode="markers+text",
                marker=dict(
                    size=data["msize"],
                    color="white",
                    opacity=0.65,
                    symbol="circle",
                    line=dict(width=0),
                ),
                text=label_col,
                textposition="top center",
                textfont=dict(size=9, color="#8899aa"),
                name="Stars",
                customdata=np.stack([data[mag_col], data["display_name"]], axis=1),
                hovertemplate=(
                    "<b>%{customdata[1]}</b><br>"
                    "RA: %{x:.2f}°<br>"
                    "Dec: %{y:.2f}°<br>"
                    "V mag: %{customdata[0]:.2f}"
                    "<extra></extra>"
                ),
                showlegend=True,
            ))

    xaxis, yaxis = mercator_sky_axes()
    base_layout = dict(
        paper_bgcolor="#0d1117",
        plot_bgcolor="#0d1117",
        font=dict(color="white"),
        xaxis=xaxis,
        yaxis=yaxis,
        legend=dict(
            bgcolor="rgba(10,16,26,0.7)",
            bordercolor="#2c4060",
            borderwidth=1,
            font=dict(size=11),
            x=0.01, y=0.99,
            xanchor="left", yanchor="top",
        ),
        margin=dict(l=60, r=40, t=60, b=60),
        height=520,
        autosize=True,
    )
    if layout:
        for key, val in layout.items():
            if key in ("xaxis", "yaxis") and isinstance(val, dict):
                base_layout[key] = {**base_layout.get(key, {}), **val}
            else:
                base_layout[key] = val
    fig.update_layout(**base_layout)
    return fig


# ── Polar (azimuthal) sky map ─────────────────────────────────────────────

DEFAULT_MAG_LIMIT = 3.5
DEFAULT_LOCAL_HOUR = 18
DEFAULT_LOCAL_MINUTE = 0
DEFAULT_LOCAL_SECOND = 0
DEFAULT_BODIES = ["Sun"]
DEFAULT_CONSTELLATION_SET = "iau"

LINE_ECLIPTIC = "Ecliptic"
LINE_HORIZON = "Horizon"

BODY_EMOJIS: dict[str, str] = {
    "Sun": "☀",
    "Moon": "🌙",
    "Mercury": "☿",
    "Venus": "♀",
    "Mars": "♂",
    "Jupiter": "♃",
    "Saturn": "♄",
}

BODY_MAP_COLORS: dict[str, str] = {
    "Sun": "#FFD700",
    "Moon": "#1565C0",
    "Mercury": "#2E7D32",
    "Venus": "#C62828",
    "Mars": "#8E0000",
    "Jupiter": "#1B5E20",
    "Saturn": "#0D47A1",
}

STAR_COLOR_ABOVE_HORIZON = "rgba(255, 255, 255, 0.72)"
STAR_COLOR_BELOW_HORIZON = "rgba(144, 238, 180, 0.88)"
STAR_LABEL_COLOR_ABOVE = "#95aac5"
STAR_LABEL_COLOR_BELOW = "rgba(144, 238, 180, 0.75)"
HORIZON_AZ_MARK_STEP = 30
HORIZON_AZ_LABEL_COLOR = "rgba(144, 238, 180, 0.95)"

_CATALOG_SOLAR_SYSTEM_NAMES = frozenset({"Sol"})


def _make_sky_body(name: str):
    """Return a ``Sun``, ``Moon``, or ``Planet`` instance by display name."""
    if name == "Sun":
        return montu.Sun()
    if name == "Moon":
        return montu.Moon()
    return montu.Planet(name)


@lru_cache(maxsize=len(CONSTELLATION_SET_IDS))
def _constellation_entries(set_id: str = DEFAULT_CONSTELLATION_SET) -> list[dict]:
    """Return parsed asterism entries for one constellation set."""
    return parse_constellation_lines(set_id=set_id)


@lru_cache(maxsize=len(CONSTELLATION_SET_IDS))
def _constellation_name_labels(set_id: str = DEFAULT_CONSTELLATION_SET) -> dict[str, str]:
    """Return display labels for constellation abbreviations."""
    return parse_constellation_names(set_id=set_id)


def _star_name(row: pd.Series) -> str:
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


def _star_size(vmag: float) -> float:
    return float(np.clip(12.5 - 1.8 * float(vmag), 2.5, 18.0))


def _north_radius(dec_deg: float) -> float:
    """North-polar azimuthal radius: NCP at 0, equator at 90."""
    return 90.0 - float(dec_deg)


def _south_radius(dec_deg: float) -> float:
    """South-polar azimuthal radius: SCP at 0, equator at 90."""
    return 90.0 + float(dec_deg)


def _observer_sidereal_time_hours(
    obs_utc,
    *,
    lat: float,
    lon: float,
    height_km: float,
) -> float:
    """Local apparent sidereal time [hours] at the observation instant."""
    observer = montu.Observer(lon=lon, lat=lat, height=height_km)
    mtime = obs_utc if isinstance(obs_utc, montu.Time) else montu.Time(obs_utc, calendar="proleptic")
    observer.site.date = mtime.jed - montu.PYEPHEM_JD_REF
    return float(observer.site.sidereal_time() * montu.RAD / 15.0)


def _equatorial_to_horizontal(
    ra_h: np.ndarray,
    dec_deg: np.ndarray,
    *,
    lat: float,
    lst_hours: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised equatorial → horizontal (azimuth°, elevation°)."""
    ha = np.mod(lst_hours - ra_h, 24.0)
    el_rad = np.arcsin(
        np.sin(dec_deg * montu.DEG) * np.sin(lat * montu.DEG)
        + np.cos(dec_deg * montu.DEG) * np.cos(lat * montu.DEG) * np.cos(ha * 15.0 * montu.DEG)
    )
    cos_el = np.cos(el_rad)
    cos_el_safe = np.where(np.abs(cos_el) < 1e-12, 1e-12, cos_el)
    az_rad = np.arctan2(
        -np.sin(ha * 15.0 * montu.DEG) * np.cos(dec_deg * montu.DEG) / cos_el_safe,
        (
            np.sin(dec_deg * montu.DEG) - np.sin(lat * montu.DEG) * np.sin(el_rad)
        ) / (np.cos(lat * montu.DEG) * cos_el_safe),
    )
    az = np.mod(az_rad * montu.RAD, 360.0)
    el = el_rad * montu.RAD
    return az, el


def _map_theta(
    ra_deg: float,
    *,
    lst_deg: float,
    meridian_view: bool,
) -> float:
    """Map equatorial RA [deg] to polar theta, optionally meridian-aligned."""
    if meridian_view:
        return (float(ra_deg) - lst_deg) % 360.0
    return float(ra_deg)


def _map_theta_array(
    ra_deg: np.ndarray,
    *,
    lst_deg: float,
    meridian_view: bool,
) -> np.ndarray:
    """Vectorised :func:`_map_theta`."""
    values = ra_deg.astype(float)
    if meridian_view:
        return np.mod(values - lst_deg, 360.0)
    return values


def _angular_axis_ticks(
    *,
    lst_hours: float,
    meridian_view: bool,
) -> tuple[list[float], list[str]]:
    """Polar angular grid: absolute RA meridians, rotated when meridian-aligned."""
    ra_hours = list(range(0, 24, 2))
    ticktext = [f"{h}h" for h in ra_hours]
    if meridian_view:
        lst_deg = lst_hours * 15.0
        pairs = sorted(
            (
                ((h * 15.0) - lst_deg) % 360.0,
                label,
            )
            for h, label in zip(ra_hours, ticktext)
        )
        tickvals = [theta for theta, _ in pairs]
        ticktext = [label for _, label in pairs]
    else:
        tickvals = [float(h * 15) for h in ra_hours]
    return tickvals, ticktext


def _without_catalog_solar_bodies(df: pd.DataFrame) -> pd.DataFrame:
    """Drop catalogue placeholders for the Sun (and similar) from the star layer."""
    names = df["Name"].astype(str)
    return df[~names.isin(_CATALOG_SOLAR_SYSTEM_NAMES)].copy()


def _polar_hip_lookup(precessed: pd.DataFrame) -> dict[int, tuple[float, float]]:
    hip_lookup: dict[int, tuple[float, float]] = {}
    for _, row in precessed.iterrows():
        hip = row.get("HIP", np.nan)
        if pd.isna(hip):
            continue
        hip_lookup[int(hip)] = (float(row["ra_deg"]), float(row["DecEpoch"]))
    return hip_lookup


def _precession_epoch_string(calendar_date: str) -> str:
    """Date-only epoch for star precession and ecliptic (noon, proleptic)."""
    return f"{_calendar_date_only(calendar_date)} 12:00:00"


def _resolve_line_flags(
    lines: list[str] | None,
    show_horizon: bool | None,
    show_ecliptic: bool | None,
) -> tuple[bool, bool]:
    """Resolve horizon/ecliptic flags from *lines* or explicit booleans."""
    if lines is not None:
        enabled = set(lines)
        return LINE_HORIZON in enabled, LINE_ECLIPTIC in enabled
    if show_horizon is not None or show_ecliptic is not None:
        horizon = bool(show_horizon) if show_horizon is not None else False
        ecliptic = bool(show_ecliptic) if show_ecliptic is not None else True
        return horizon, ecliptic
    return False, True


def _resolve_bodies(bodies: list[str] | None) -> list[str]:
    """Return the body list to plot; empty list means plot none."""
    if bodies is None:
        return list(DEFAULT_BODIES)
    return list(bodies)


def _calendar_date_only(date_str: str) -> str:
    """Return a proleptic date string without time.

    BCE input is converted from historical year numbering (``bce 2026``)
    to astronomical year numbering (``-2025``), which is what
    ``montu.Time(..., calendar='proleptic')`` expects for stable JDs.
    """
    tokens = date_str.strip().split()
    if tokens and tokens[0].lower() == "bce":
        year_s, month_s, day_s = tokens[1][:10].split("-")
        astro_year = 1 - int(year_s)
        return f"{astro_year:04d}-{month_s}-{day_s}"
    return tokens[0][:10]


def local_solar_to_utc_time(
    calendar_date: str,
    hour: int,
    minute: int,
    second: int,
    lon: float,
    calendar: str = "mixed",
) -> montu.Time:
    """Convert local solar time at *lon* to a UTC ``montu.Time`` object."""
    date_only = _calendar_date_only(calendar_date)
    midnight = f"{date_only} 00:00:00"
    local_decimal = hour + minute / 60.0 + second / 3600.0
    utc_decimal = local_decimal - lon / 15.0
    t0 = montu.Time(midnight, calendar=calendar)
    return t0 + utc_decimal * montu.HOUR


def local_solar_to_utc_datepro(
    calendar_date: str,
    hour: int,
    minute: int,
    second: int,
    lon: float,
    calendar: str = "mixed",
) -> str:
    """Convert local solar time at *lon* to a proleptic UTC ``datepro`` string."""
    t_obs = local_solar_to_utc_time(calendar_date, hour, minute, second, lon, calendar=calendar).get_readable()
    return str(t_obs.readable.datepro)


def _local_time_label(hour: int, minute: int, second: int) -> str:
    return f"{hour:02d}:{minute:02d}:{second:02d} local"


def _sidereal_time_label(lst_hours: float) -> str:
    """Format local apparent sidereal time as ``LST HH:MM:SS``."""
    total_seconds = int(round((lst_hours % 24.0) * 3600.0))
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60
    return f"LST {hours:02d}:{minutes:02d}:{seconds:02d}"


def _compute_body_positions(
    obs_utc,
    *,
    bodies: list[str],
    lat: float,
    lon: float,
    height_km: float,
) -> dict[str, tuple[float, float, float, float, float, float]]:
    """Return (RA°, Dec°, RA h, Dec°, az°, el°) for each selected body."""
    observer = montu.Observer(lon=lon, lat=lat, height=height_km)
    mtime = obs_utc if isinstance(obs_utc, montu.Time) else montu.Time(obs_utc, calendar="proleptic")
    positions: dict[str, tuple[float, float, float, float, float, float]] = {}
    for name in bodies:
        body = _make_sky_body(name)
        body.where_in_sky(at=mtime, observer=observer)
        ra_h = float(body.position.RAEpoch)
        ra_deg = ra_h * 15.0
        dec_deg = float(body.position.DecEpoch)
        az = float(body.position.az)
        el = float(body.position.el)
        positions[name] = (ra_deg, dec_deg, ra_h, dec_deg, az, el)
    return positions


def _body_trace(
    name: str,
    ra_deg: float,
    dec_deg: float,
    hemisphere: str,
    *,
    ra_h: float = 0.0,
    az_deg: float = 0.0,
    el_deg: float = 0.0,
    show_horizontal: bool = False,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Plot one solar-system body with emoji marker on one hemisphere."""
    if hemisphere == "north":
        if dec_deg < 0.0:
            return go.Scatterpolar(
                theta=[], r=[], mode="markers", showlegend=False, name=name,
            )
        radius = _north_radius(dec_deg)
    else:
        if dec_deg > 0.0:
            return go.Scatterpolar(
                theta=[], r=[], mode="markers", showlegend=False, name=name,
            )
        radius = _south_radius(dec_deg)

    emoji = BODY_EMOJIS.get(name, "★")
    color = BODY_MAP_COLORS.get(name, "#FFD700")
    if show_horizontal:
        hovertemplate = (
            f"<b>{name}</b> {emoji}<br>"
            f"RA(epoch): {ra_h:.4f} h<br>"
            f"Dec(epoch): {dec_deg:.4f}°<br>"
            f"Azimuth: {az_deg:.4f}°<br>"
            f"Elevation: {el_deg:.4f}°"
            "<extra></extra>"
        )
    else:
        hovertemplate = (
            f"<b>{name}</b> {emoji}<br>"
            f"RA(epoch): {ra_h:.4f} h<br>"
            f"Dec(epoch): {dec_deg:.4f}°"
            "<extra></extra>"
        )
    return go.Scatterpolar(
        theta=[_map_theta(ra_deg, lst_deg=lst_deg, meridian_view=meridian_view)],
        r=[radius],
        mode="text",
        text=[emoji],
        textposition="middle center",
        textfont=dict(size=22, color=color),
        hovertemplate=hovertemplate,
        name=name,
        legendgroup=name,
        showlegend=False,
    )


def _compute_ecliptic_curve(precession_date: str) -> tuple[np.ndarray, np.ndarray]:
    """Equatorial (RA°, Dec°) samples along the ecliptic at *precession_date*."""
    from pyplanets.core.coordinates import true_obliquity

    mtime = montu.Time(precession_date, calendar="mixed")
    eps = np.radians(float(true_obliquity(mtime.obj_pyplanet)))
    lam = np.linspace(0.0, 2.0 * np.pi, 361)
    dec_rad = np.arcsin(np.sin(eps) * np.sin(lam))
    ra_rad = np.arctan2(np.sin(lam) * np.cos(eps), np.cos(lam))
    ra_deg = np.degrees(ra_rad) % 360.0
    dec_deg = np.degrees(dec_rad)
    return ra_deg, dec_deg


def _sky_curve_trace(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    hemisphere: str,
    *,
    name: str,
    color: str,
    width: float = 1.6,
    dash: str | None = None,
    showlegend: bool = True,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Plot an equatorial great-circle segment on one polar hemisphere."""
    theta: list[float | None] = []
    radius: list[float | None] = []
    segment_theta: list[float] = []
    segment_r: list[float] = []

    def _flush() -> None:
        nonlocal segment_theta, segment_r
        if len(segment_theta) >= 2:
            theta.extend(segment_theta)
            radius.extend(segment_r)
            theta.append(None)
            radius.append(None)
        segment_theta = []
        segment_r = []

    for ra, dec in zip(ra_deg, dec_deg):
        if hemisphere == "north":
            if dec < 0.0:
                _flush()
                continue
            r = _north_radius(dec)
        else:
            if dec > 0.0:
                _flush()
                continue
            r = _south_radius(dec)

        if segment_theta and abs(
            _map_theta(float(ra), lst_deg=lst_deg, meridian_view=meridian_view)
            - segment_theta[-1]
        ) > 180.0:
            _flush()

        segment_theta.append(
            _map_theta(float(ra), lst_deg=lst_deg, meridian_view=meridian_view),
        )
        segment_r.append(r)

    _flush()
    if theta and theta[-1] is None:
        theta.pop()
        radius.pop()

    line = dict(color=color, width=width)
    if dash:
        line["dash"] = dash

    return go.Scatterpolar(
        theta=theta,
        r=radius,
        mode="lines",
        line=line,
        hoverinfo="skip",
        name=name,
        showlegend=showlegend,
    )


def _compute_horizon_curve(
    obs_utc,
    *,
    lat: float,
    lon: float,
    height_km: float,
    n_samples: int = 361,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Equatorial samples and azimuth [°] along the el=0° horizon."""
    observer = montu.Observer(lon=lon, lat=lat, height=height_km)
    mtime = obs_utc if isinstance(obs_utc, montu.Time) else montu.Time(obs_utc, calendar="proleptic")
    lst = _observer_sidereal_time_hours(
        obs_utc, lat=lat, lon=lon, height_km=height_km,
    )
    lat_rad = lat * montu.DEG

    if abs(lat) >= 89.5:
        ra_deg = np.linspace(0.0, 360.0, n_samples)
        dec_deg = np.zeros(n_samples)
        az_deg = np.linspace(0.0, 360.0, n_samples, endpoint=False)
        return ra_deg, dec_deg, az_deg

    records: list[tuple[float, float, float]] = []
    for ha_h in np.linspace(-12.0, 12.0, n_samples, endpoint=False):
        ha_rad = ha_h * 15.0 * montu.DEG
        dec_rad = np.arctan2(
            -np.cos(lat_rad) * np.cos(ha_rad),
            np.sin(lat_rad),
        )
        dec = dec_rad * montu.RAD
        ra_h = np.mod(lst - ha_h, 24.0)
        az, el = montu.Astro.where_in_sky(
            RA=ra_h, Dec=dec, at=mtime, observer=observer,
        )
        if abs(el) > 0.05:
            continue
        records.append((az, ra_h * 15.0, dec))

    if not records:
        return np.array([]), np.array([]), np.array([])

    records.sort(key=lambda row: row[0])
    ra_deg = np.array([row[1] for row in records])
    dec_deg = np.array([row[2] for row in records])
    az_deg = np.array([row[0] for row in records])
    return ra_deg, dec_deg, az_deg


def _horizon_azimuth_marks_trace(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    az_deg: np.ndarray,
    hemisphere: str,
    *,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Green azimuth labels on the horizon curve."""
    if len(ra_deg) == 0:
        return go.Scatterpolar(theta=[], r=[], mode="text", showlegend=False)

    theta: list[float] = []
    radius: list[float] = []
    text: list[str] = []
    for mark_az in range(0, 360, HORIZON_AZ_MARK_STEP):
        delta = np.abs((az_deg - mark_az + 180.0) % 360.0 - 180.0)
        idx = int(np.argmin(delta))
        ra = float(ra_deg[idx])
        dec = float(dec_deg[idx])
        if hemisphere == "north":
            if dec < 0.0:
                continue
            r = _north_radius(dec)
        else:
            if dec > 0.0:
                continue
            r = _south_radius(dec)
        theta.append(_map_theta(ra, lst_deg=lst_deg, meridian_view=meridian_view))
        radius.append(r)
        text.append(f"{mark_az}°")

    return go.Scatterpolar(
        theta=theta,
        r=radius,
        mode="text",
        text=text,
        textfont=dict(size=9, color=HORIZON_AZ_LABEL_COLOR),
        hoverinfo="skip",
        showlegend=False,
        name="Azimuth marks",
    )


def _annotate_star_elevations(
    stars: pd.DataFrame,
    obs_utc,
    *,
    lat: float,
    lon: float,
    height_km: float,
) -> pd.DataFrame:
    """Add observer azimuth and elevation [deg] for each precessed star."""
    lst_hours = _observer_sidereal_time_hours(
        obs_utc, lat=lat, lon=lon, height_km=height_km,
    )
    ra_h = stars["RAEpoch"].astype(float).to_numpy()
    dec_deg = stars["DecEpoch"].astype(float).to_numpy()
    az, el = _equatorial_to_horizontal(ra_h, dec_deg, lat=lat, lst_hours=lst_hours)

    annotated = stars.copy()
    annotated["az"] = az
    annotated["el"] = el
    return annotated


def _horizon_trace(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    hemisphere: str,
    *,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Plot the observer horizon (elevation = 0°) on one polar hemisphere."""
    return _sky_curve_trace(
        ra_deg,
        dec_deg,
        hemisphere,
        name="Horizon",
        color="rgba(46, 204, 113, 0.92)",
        width=2.0,
        showlegend=True,
        lst_deg=lst_deg,
        meridian_view=meridian_view,
    )


def _ecliptic_trace(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    hemisphere: str,
    *,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Plot the ecliptic great circle on one polar hemisphere."""
    return _sky_curve_trace(
        ra_deg,
        dec_deg,
        hemisphere,
        name="Ecliptic",
        color="rgba(255, 193, 80, 0.78)",
        width=1.6,
        dash="dash",
        showlegend=True,
        lst_deg=lst_deg,
        meridian_view=meridian_view,
    )


def _prepare_precessed_data(
    calendar_date: str,
    precessed_star_data: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[int, tuple[float, float]], tuple[np.ndarray, np.ndarray]]:
    """Return precessed catalogue, HIP lookup, and ecliptic samples."""
    epoch_str = _precession_epoch_string(calendar_date)
    ecliptic_ra, ecliptic_dec = _compute_ecliptic_curve(epoch_str)

    if precessed_star_data is not None:
        precessed = precessed_star_data.copy()
    else:
        mtime = montu.Time(epoch_str, calendar="proleptic")
        stars = montu.Stars()
        precessed = stars.where_in_space(at=mtime).data.copy()

    if "ra_deg" not in precessed.columns:
        precessed["ra_deg"] = precessed["RAEpoch"] * 15.0
    if "display_name" not in precessed.columns:
        precessed["display_name"] = precessed.apply(_star_name, axis=1)

    hip_lookup = _polar_hip_lookup(precessed)
    return precessed, hip_lookup, (ecliptic_ra, ecliptic_dec)


def _constellation_boundaries_trace(
    *,
    hemisphere: str,
    at=None,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Build constellation boundaries for one hemisphere."""
    from montu.stars import parse_constellation_boundaries
    theta: list[float | None] = []
    radius: list[float | None] = []

    for poly in parse_constellation_boundaries(at=at):
        poly_theta = []
        poly_r = []
        for ra_deg, dec_deg in poly["points"]:
            if hemisphere == "north":
                if dec_deg < -5.0:
                    continue
                r_a = _north_radius(dec_deg)
            else:
                if dec_deg > 5.0:
                    continue
                r_a = _south_radius(dec_deg)
            poly_theta.append(_map_theta(ra_deg, lst_deg=lst_deg, meridian_view=meridian_view))
            poly_r.append(r_a)
            
        if poly_theta:
            theta.extend(poly_theta + [None])
            radius.extend(poly_r + [None])

    return go.Scatterpolar(
        theta=theta,
        r=radius,
        mode="lines",
        line=dict(color="rgba(230, 120, 170, 0.45)", width=0.8),
        hoverinfo="skip",
        name="Boundaries",
        showlegend=False,
    )


def _asterism_trace(
    *,
    hip_lookup: dict[int, tuple[float, float]],
    hemisphere: str,
    label_positions: dict[str, list[tuple[float, float]]],
    constellation_set: str = DEFAULT_CONSTELLATION_SET,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Build asterism lines for one hemisphere; fill *label_positions*."""
    theta: list[float | None] = []
    radius: list[float | None] = []

    for entry in _constellation_entries(constellation_set):
        abbrev = entry["abbrev"]
        for hip_a, hip_b in entry.get("segments", []):
            pa = hip_lookup.get(hip_a)
            pb = hip_lookup.get(hip_b)
            if pa is None or pb is None:
                continue
            ra_a, dec_a = pa
            ra_b, dec_b = pb

            if hemisphere == "north":
                if dec_a < 0.0 or dec_b < 0.0:
                    continue
                r_a = _north_radius(dec_a)
                r_b = _north_radius(dec_b)
            else:
                if dec_a > 0.0 or dec_b > 0.0:
                    continue
                r_a = _south_radius(dec_a)
                r_b = _south_radius(dec_b)

            theta.extend([
                _map_theta(ra_a, lst_deg=lst_deg, meridian_view=meridian_view),
                _map_theta(ra_b, lst_deg=lst_deg, meridian_view=meridian_view),
                None,
            ])
            radius.extend([r_a, r_b, None])
            label_positions.setdefault(abbrev, []).append(pa)
            label_positions.setdefault(abbrev, []).append(pb)

    return go.Scatterpolar(
        theta=theta,
        r=radius,
        mode="lines",
        line=dict(color="rgba(125,145,170,0.52)", width=1.0),
        hoverinfo="skip",
        name="Asterisms",
        showlegend=False,
    )


def _constellation_labels_trace(
    *,
    label_positions: dict[str, list[tuple[float, float]]],
    hemisphere: str,
    constellation_set: str = DEFAULT_CONSTELLATION_SET,
    constellation_full_names: bool = False,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Soft grey constellation labels (same style as mercator_sky_map)."""
    theta: list[float] = []
    radius: list[float] = []
    text: list[str] = []
    name_labels = _constellation_name_labels(constellation_set) if constellation_full_names else {}

    for abbrev, coords in label_positions.items():
        ra_mean = float(np.mean([c[0] for c in coords]))
        dec_mean = float(np.mean([c[1] for c in coords]))
        if hemisphere == "north":
            if dec_mean < 0.0:
                continue
            r = _north_radius(dec_mean)
        else:
            if dec_mean > 0.0:
                continue
            r = _south_radius(dec_mean)
        theta.append(_map_theta(ra_mean, lst_deg=lst_deg, meridian_view=meridian_view))
        radius.append(r)
        text.append(name_labels.get(abbrev, abbrev))

    return go.Scatterpolar(
        theta=theta,
        r=radius,
        mode="text",
        text=text,
        textfont=dict(size=9, color="rgba(130, 140, 155, 0.42)"),
        hoverinfo="skip",
        name="Constellation labels",
        showlegend=False,
    )


def _stars_trace(
    data: pd.DataFrame,
    hemisphere: str,
    *,
    shade_below_horizon: bool = False,
    lst_deg: float = 0.0,
    meridian_view: bool = False,
) -> go.Scatterpolar:
    """Build stars trace for one hemisphere."""
    if hemisphere == "north":
        subset = data[data["DecEpoch"] >= 0.0].copy()
        subset["r"] = subset["DecEpoch"].apply(_north_radius)
    else:
        subset = data[data["DecEpoch"] <= 0.0].copy()
        subset["r"] = subset["DecEpoch"].apply(_south_radius)

    subset["marker_size"] = subset["Vmag"].apply(_star_size)
    subset["label"] = subset.apply(
        lambda row: row["display_name"] if float(row["Vmag"]) <= 2.0 else "",
        axis=1,
    )
    if shade_below_horizon and "el" in subset.columns:
        below = subset["el"].astype(float) < 0.0
        subset["marker_color"] = np.where(
            below, STAR_COLOR_BELOW_HORIZON, STAR_COLOR_ABOVE_HORIZON,
        )
        subset["label_color"] = np.where(
            below, STAR_LABEL_COLOR_BELOW, STAR_LABEL_COLOR_ABOVE,
        )
        customdata = np.stack(
            [
                subset["display_name"],
                subset["Vmag"],
                subset["RAEpoch"],
                subset["DecEpoch"],
                subset["el"],
                subset["az"],
            ],
            axis=1,
        )
        hovertemplate = (
            "<b>%{customdata[0]}</b><br>"
            "V mag: %{customdata[1]:.2f}<br>"
            "RA(epoch): %{customdata[2]:.4f} h<br>"
            "Dec(epoch): %{customdata[3]:.4f}°<br>"
            "Azimuth: %{customdata[5]:.4f}°<br>"
            "Elevation: %{customdata[4]:.4f}°"
            "<extra></extra>"
        )
    else:
        subset["marker_color"] = STAR_COLOR_ABOVE_HORIZON
        subset["label_color"] = STAR_LABEL_COLOR_ABOVE
        customdata = np.stack(
            [
                subset["display_name"],
                subset["Vmag"],
                subset["RAEpoch"],
                subset["DecEpoch"],
            ],
            axis=1,
        )
        hovertemplate = (
            "<b>%{customdata[0]}</b><br>"
            "V mag: %{customdata[1]:.2f}<br>"
            "RA(epoch): %{customdata[2]:.4f} h<br>"
            "Dec(epoch): %{customdata[3]:.4f}°"
            "<extra></extra>"
        )

    return go.Scatterpolar(
        theta=_map_theta_array(
            subset["ra_deg"].to_numpy(),
            lst_deg=lst_deg,
            meridian_view=meridian_view,
        ),
        r=subset["r"],
        mode="markers+text",
        text=subset["label"],
        textposition="top center",
        textfont=dict(size=9, color=subset["label_color"]),
        marker=dict(
            size=subset["marker_size"],
            color=subset["marker_color"],
            line=dict(width=0),
        ),
        customdata=customdata,
        hovertemplate=hovertemplate,
        name="Stars",
        showlegend=False,
    )


def _map_title(
    *,
    observer_name: str,
    lat: float,
    lon: float,
    calendar_date: str,
    local_time: str,
    lst_hours: float,
    mag_limit: float,
) -> str:
    obs = (
        f"{observer_name} (lat {lat:.2f}°, lon {lon:.2f}°)"
        if observer_name
        else f"lat {lat:.2f}°, lon {lon:.2f}°"
    )
    date_only = _calendar_date_only(calendar_date)
    lst_label = _sidereal_time_label(lst_hours)
    return (
        f"🌌 Sky Map · {obs} · "
        f"{date_only} {local_time} · {lst_label} · V ≤ {mag_limit:.1f}"
    )


def _title_layout(text: str) -> dict:
    """Compact left-aligned title for the sky-map header row."""
    return dict(
        text=text,
        font=dict(size=11, color="#b8c5d6"),
        x=0.0,
        xanchor="left",
        y=1.0,
        yanchor="top",
        pad=dict(t=6, b=2, l=10, r=0),
    )


def polar_sky_map_figure(
    *,
    hemisphere: str,
    calendar_date: str,
    local_time: str,
    mag_limit: float,
    stars: pd.DataFrame,
    hip_lookup: dict[int, tuple[float, float]],
    ecliptic_ra: np.ndarray,
    ecliptic_dec: np.ndarray,
    horizon_ra: np.ndarray,
    horizon_dec: np.ndarray,
    horizon_az: np.ndarray,
    body_positions: dict[str, tuple[float, float, float, float, float, float]],
    selected_bodies: list[str],
    show_ecliptic: bool = True,
    show_horizon: bool = False,
    shade_below_horizon: bool = False,
    meridian_view: bool = False,
    constellation_set: str = DEFAULT_CONSTELLATION_SET,
    lst_deg: float = 0.0,
    lst_hours: float = 0.0,
    observer_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    show_constellation_lines: bool = True,
    show_constellation_boundaries: bool = False,
    show_constellation_labels: bool = True,
    constellation_full_names: bool = False,
    at=None,
) -> go.Figure:
    """Build one polar sky map for *hemisphere* (``north`` or ``south``)."""
    fig = go.Figure()
    label_positions: dict[str, list[tuple[float, float]]] = {}

    map_traces: list[go.Scatterpolar] = []
    if show_horizon:
        map_traces.append(
            _horizon_trace(
                horizon_ra, horizon_dec, hemisphere,
                lst_deg=lst_deg, meridian_view=meridian_view,
            ),
        )
        map_traces.append(
            _horizon_azimuth_marks_trace(
                horizon_ra, horizon_dec, horizon_az, hemisphere,
                lst_deg=lst_deg, meridian_view=meridian_view,
            ),
        )
    if show_ecliptic:
        map_traces.append(
            _ecliptic_trace(
                ecliptic_ra, ecliptic_dec, hemisphere,
                lst_deg=lst_deg, meridian_view=meridian_view,
            ),
        )
    asterism_trace = _asterism_trace(
        hip_lookup=hip_lookup,
        hemisphere=hemisphere,
        label_positions=label_positions,
        constellation_set=constellation_set,
        lst_deg=lst_deg,
        meridian_view=meridian_view,
    )
    if show_constellation_lines:
        map_traces.append(asterism_trace)

    if show_constellation_labels and label_positions:
        map_traces.append(
            _constellation_labels_trace(
                label_positions=label_positions,
                hemisphere=hemisphere,
                constellation_set=constellation_set,
                constellation_full_names=constellation_full_names,
                lst_deg=lst_deg,
                meridian_view=meridian_view,
            )
        )

    if show_constellation_boundaries:
        map_traces.append(
            _constellation_boundaries_trace(
                hemisphere=hemisphere,
                at=at,
                lst_deg=lst_deg,
                meridian_view=meridian_view,
            )
        )

    map_traces.append(
        _stars_trace(
            stars,
            hemisphere=hemisphere,
            shade_below_horizon=shade_below_horizon,
            lst_deg=lst_deg,
            meridian_view=meridian_view,
        )
    )
    for name in selected_bodies:
        if name not in body_positions:
            continue
        ra_deg, dec_deg, ra_h, _, az_deg, el_deg = body_positions[name]
        map_traces.append(
            _body_trace(
                name, ra_deg, dec_deg, hemisphere,
                ra_h=ra_h,
                az_deg=az_deg,
                el_deg=el_deg,
                show_horizontal=show_horizon,
                lst_deg=lst_deg,
                meridian_view=meridian_view,
            ),
        )

    for trace in map_traces:
        fig.add_trace(trace)

    title = _map_title(
        observer_name=observer_name,
        lat=lat,
        lon=lon,
        calendar_date=calendar_date,
        local_time=local_time,
        lst_hours=lst_hours,
        mag_limit=mag_limit,
    )

    tickvals, ticktext = _angular_axis_ticks(
        lst_hours=lst_hours,
        meridian_view=meridian_view,
    )

    fig.update_layout(
        title=_title_layout(title),
        paper_bgcolor="#0d1117",
        plot_bgcolor="#0d1117",
        font=dict(color="#e6edf7"),
        margin=dict(l=40, r=12, t=40, b=40),
        legend=dict(
            x=0.01,
            y=0.99,
            xanchor="left",
            yanchor="top",
            bgcolor="rgba(10, 16, 26, 0.75)",
            bordercolor="#2c4060",
            borderwidth=1,
            font=dict(size=11, color="#e6edf7"),
        ),
        polar=dict(
            bgcolor="#0d1117",
            domain=dict(x=[0.0, 1.0], y=[0.0, 1.0]),
            radialaxis=dict(
                range=[0, 90],
                tickmode="array",
                tickvals=[0, 15, 30, 45, 60, 75, 90],
                ticktext=["90°", "75°", "60°", "45°", "30°", "15°", "0°"],
                gridcolor="#2a3648",
                linecolor="#2a3648",
                angle=90,
            ),
            angularaxis=dict(
                direction="clockwise",
                rotation=90,
                tickmode="array",
                tickvals=tickvals,
                ticktext=ticktext,
                gridcolor="#1f2a38",
                linecolor="#2a3648",
            ),
        ),
        autosize=True,
    )
    return fig


def polar_sky_map(
    calendar_date: str,
    *,
    local_hour: int = DEFAULT_LOCAL_HOUR,
    local_minute: int = DEFAULT_LOCAL_MINUTE,
    local_second: int = DEFAULT_LOCAL_SECOND,
    observer: montu.Observer,
    mag_limit: float = DEFAULT_MAG_LIMIT,
    bodies: list[str] | None = None,
    lines: list[str] | None = None,
    show_horizon: bool | None = None,
    show_ecliptic: bool | None = None,
    meridian_view: bool = False,
    constellation_set: str = DEFAULT_CONSTELLATION_SET,
    show_constellation_lines: bool = True,
    show_constellation_boundaries: bool = False,
    show_constellation_labels: bool = True,
    constellation_full_names: bool = False,
    observer_name: str = "",
    precessed_star_data: pd.DataFrame | None = None,
) -> tuple[go.Figure, go.Figure]:
    """Build north and south polar sky map figures.

    Parameters
    ----------
    calendar_date : str
        Proleptic calendar date (``bce YYYY-MM-DD`` or ``YYYY-MM-DD``).
    local_hour, local_minute, local_second : int
        Local solar time at the observer longitude.
    observer : montu.Observer
        Observing site (latitude, longitude, height).
    mag_limit : float
        Faint limit for background stars [V magnitude].
    bodies : list[str], optional
        Solar-system bodies to plot (``Sun``, ``Moon``, planet names).
        ``None`` plots the Sun; ``[]`` plots none.
    lines : list[str], optional
        Overlay lines to draw (``LINE_ECLIPTIC``, ``LINE_HORIZON``).
        ``None`` defaults to ecliptic only; ``[]`` plots none.
    show_horizon, show_ecliptic : bool, optional
        Alternative to *lines* when explicit booleans are preferred.
    meridian_view : bool
        Rotate the map so the meridian aligns with the vertical axis.
    constellation_set : str
        Constellation asterism set id (see ``CONSTELLATION_SET_IDS``).
    observer_name : str
        Optional label for the map title.
    precessed_star_data : pandas.DataFrame, optional
        Pre-precessed stellar catalogue.  When ``None``, loads
        ``montu.Stars()`` and precesses to noon on *calendar_date*.

    Returns
    -------
    tuple[plotly.graph_objects.Figure, plotly.graph_objects.Figure]
        North- and south-hemisphere polar sky maps.
    """
    selected = _resolve_bodies(bodies)
    show_horizon_flag, show_ecliptic_flag = _resolve_line_flags(
        lines, show_horizon, show_ecliptic,
    )
    local_time = _local_time_label(local_hour, local_minute, local_second)

    lat = float(observer.lat)
    lon = float(observer.lon)
    height_km = float(observer.height)

    obs_time = local_solar_to_utc_time(
        calendar_date, local_hour, local_minute, local_second, lon,
    )
    obs_utc = str(obs_time.get_readable().readable.datepro)
    lst_hours = _observer_sidereal_time_hours(
        obs_time, lat=lat, lon=lon, height_km=height_km,
    )
    lst_deg = lst_hours * 15.0

    precessed, hip_lookup, (ecliptic_ra, ecliptic_dec) = _prepare_precessed_data(
        calendar_date, precessed_star_data,
    )

    if show_horizon_flag:
        horizon_ra, horizon_dec, horizon_az = _compute_horizon_curve(
            obs_time,
            lat=lat,
            lon=lon,
            height_km=height_km,
        )
    else:
        horizon_ra = np.array([])
        horizon_dec = np.array([])
        horizon_az = np.array([])

    body_positions = (
        _compute_body_positions(
            obs_time,
            bodies=selected,
            lat=lat,
            lon=lon,
            height_km=height_km,
        )
        if selected
        else {}
    )

    stars = _without_catalog_solar_bodies(precessed)
    stars = stars[
        (stars["Vmag"] <= float(mag_limit))
        & stars["RAEpoch"].notna()
        & stars["DecEpoch"].notna()
    ].copy()

    if show_horizon_flag:
        stars = _annotate_star_elevations(
            stars,
            obs_time,
            lat=lat,
            lon=lon,
            height_km=height_km,
        )

    fig_kwargs = dict(
        calendar_date=calendar_date,
        local_time=local_time,
        mag_limit=float(mag_limit),
        stars=stars,
        hip_lookup=hip_lookup,
        ecliptic_ra=ecliptic_ra,
        ecliptic_dec=ecliptic_dec,
        horizon_ra=horizon_ra,
        horizon_dec=horizon_dec,
        horizon_az=horizon_az,
        body_positions=body_positions,
        selected_bodies=selected,
        show_ecliptic=show_ecliptic_flag,
        show_horizon=show_horizon_flag,
        shade_below_horizon=show_horizon_flag,
        meridian_view=meridian_view,
        constellation_set=constellation_set,
        show_constellation_lines=show_constellation_lines,
        show_constellation_boundaries=show_constellation_boundaries,
        show_constellation_labels=show_constellation_labels,
        constellation_full_names=constellation_full_names,
        at=obs_time,
        lst_deg=lst_deg,
        lst_hours=lst_hours,
        observer_name=observer_name,
        lat=lat,
        lon=lon,
    )

    fig_north = polar_sky_map_figure(hemisphere="north", **fig_kwargs)
    fig_south = polar_sky_map_figure(hemisphere="south", **fig_kwargs)
    return fig_north, fig_south
