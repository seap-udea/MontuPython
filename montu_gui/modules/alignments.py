"""
Star alignment analysis for MontuPython Desktop.

Computes which bright stars pass through a given azimuth/elevation direction
at a configurable observer location and date range, using stellar precession.
No Qt dependency — pure computation and Plotly output.

The core insight: for a fixed direction on the horizon (azimuth + elevation),
the declination that passes through that direction is:

    dec = arcsin( sin(lat)·sin(el) + cos(lat)·cos(el)·cos(az) )

This is constant for a fixed structure at a fixed latitude. As the Earth's
axis precesses (~26 000-year cycle), different stars sweep past that
declination, creating alignment windows that last centuries.

Famous example: the northern shaft of the King's Chamber of the Great Pyramid
of Khufu (Giza) points at az≈0°, el≈31.7°, placing the target declination
at ~88–89°. Around 2 560 BCE the star Thuban (α Draconis) — then the North
Pole Star — passed through that region of the sky.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from dataclasses import dataclass, field
from pathlib import Path

from montu_gui.utils.debug import timed_block
from montu_gui.utils.plotly_html import figure_to_html
from montu.stars import mercator_sky_map

# ── Great Pyramid of Khufu (Cheops) defaults ─────────────────────────────────
# King's Chamber northern shaft — "imperishable stars" / Thuban alignment.
# Elevation ~31.7° (Gantenbrink 1992), azimuth ~0° (simplified from ~2° W of N).
DEFAULT_AZ            = 0.0      # degrees from North, clockwise  (N=0°, E=90°)
DEFAULT_EL            = 31.7     # degrees above the horizon
DEFAULT_YEAR_START    = 2620     # BCE
DEFAULT_YEAR_END      = 2420     # BCE
DEFAULT_ERA_START     = "bce"
DEFAULT_ERA_END       = "bce"
DEFAULT_MAG_LIMIT     = 4.0      # visual (Johnson V) magnitude
DEFAULT_DEC_TOL       = 1.0      # degrees — tolerance around target declination
DEFAULT_N_EPOCHS      = 5        # sample epochs spread across the year range

# Giza plateau — Great Pyramid site coordinates
DEFAULT_LAT           = 29.9792  # °N
DEFAULT_LON           = 31.1342  # °E
DEFAULT_ALT_M         = 75.0     # metres above sea level
DEFAULT_OBSERVER_NAME = "Giza"


# ── helpers ───────────────────────────────────────────────────────────────────

def _import_montu():
    try:
        import montu
        return montu
    except Exception as exc:
        raise ImportError(f"Cannot import montu: {exc}") from exc


def compute_target_declination(az: float, el: float, lat: float) -> float:
    """Declination of the sky direction (az, el) at observer latitude lat.

    For a fixed architectural direction (a shaft, a doorway, an alignment
    stone), all stars that pass through it share this epoch-of-date declination.
    The formula is exact for a spherical Earth model.

    Parameters
    ----------
    az : float
        Azimuth [degrees], measured from **North** clockwise (N=0°, E=90°,
        S=180°, W=270°).
    el : float
        Elevation (altitude) above the horizon [degrees].
    lat : float
        Geodetic latitude of the observer [degrees].

    Returns
    -------
    float
        Declination [degrees].
    """
    az_r  = np.radians(az)
    el_r  = np.radians(el)
    lat_r = np.radians(lat)
    sin_dec = (
        np.sin(lat_r) * np.sin(el_r)
        + np.cos(lat_r) * np.cos(el_r) * np.cos(az_r)
    )
    return float(np.degrees(np.arcsin(np.clip(sin_dec, -1.0, 1.0))))



def circumpolar_declination(lat: float) -> tuple[float, str]:
    """Return (declination limit, pole) for circumpolar stars.

    Northern sites: stars with dec ≥ limit never set (circumpolar around NCP).
    Southern sites: stars with dec ≤ limit never set (circumpolar around SCP).
    """
    if lat >= 0:
        return 90.0 - lat, "north"
    return -90.0 - lat, "south"


def _circumpolar_label(lat: float, limit: float, pole: str) -> str:
    if pole == "north":
        return f"circumpolar  dec ≥ {limit:+.1f}°"
    return f"circumpolar  dec ≤ {limit:+.1f}°"


def _mag_to_marker_size(vmag: float) -> float:
    return float(np.clip(13.0 - 2.0 * vmag, 3.0, 22.0))


def _star_display_name_from_row(row) -> str:
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


def _prepare_sky_frames(
    sky_df: pd.DataFrame,
    stars_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Add display columns; split background vs aligned stars."""
    bg = sky_df.copy()
    bg["ra_deg"] = bg["RAEpoch"] * 15.0
    bg["msize"] = bg["Vmag"].apply(_mag_to_marker_size)
    bg["display_name"] = bg.apply(_star_display_name_from_row, axis=1)

    aligned_names: set[str] = set()
    if not stars_df.empty:
        al = stars_df.copy()
        al["ra_deg"] = al["RAEpoch"] * 15.0
        al["msize"] = al["Vmag"].apply(_mag_to_marker_size) * 1.9
        al["display_name"] = al.apply(_star_display_name_from_row, axis=1)
        aligned_names = set(al["display_name"])
        for _, row in stars_df.iterrows():
            pn = str(row.get("ProperName", ""))
            n = str(row.get("Name", ""))
            aligned_names.add(pn if pn not in ("", "nan", "None") else n)
    else:
        al = stars_df.copy()

    bg_bg = bg[~bg["display_name"].isin(aligned_names)]
    return bg_bg, al, aligned_names


def _year_to_montu_str(year: int, era: str) -> str:
    """Return a montu proleptic date string for mid-year at a given historical year."""
    y = max(1, int(year))
    if era.lower() == "bce":
        return f"bce {y:04d}-07-01 00:00:00"
    return f"{y:04d}-07-01 00:00:00"


def _era_label(era: str, year: int) -> str:
    return f"{year} BCE" if era.lower() == "bce" else f"{year} CE"


def _format_ra_hms(ra_hours: float) -> str:
    """Decimal hours → ``HH h MM m SS.s s`` string."""
    if pd.isna(ra_hours):
        return "—"
    h = int(ra_hours)
    rem = (ra_hours - h) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{h:02d}h {m:02d}m {s:04.1f}s"


def _format_dec_dms(dec_deg: float) -> str:
    """Decimal degrees → ``±DD° MM′ SS.s″`` string."""
    if pd.isna(dec_deg):
        return "—"
    sign = "+" if dec_deg >= 0 else "−"
    abs_dec = abs(dec_deg)
    d = int(abs_dec)
    rem = (abs_dec - d) * 60
    m = int(rem)
    s = (rem - m) * 60
    return f"{sign}{d:02d}° {m:02d}′ {s:04.1f}″"


def star_display_name(row) -> str:
    """Return proper name when available, else catalogue name."""
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


RESULT_TABLE_COLUMNS = (
    "Star",
    "V mag",
    "RA at best epoch",
    "Dec at best epoch",
    "Δ Dec",
    "Best year",
)


def alignment_table_row(row) -> list[str]:
    """One results-table row for the desktop GUI."""
    return [
        star_display_name(row),
        f"{row['Vmag']:.2f}",
        _format_ra_hms(row["RAEpoch"]),
        _format_dec_dms(row["DecEpoch"]),
        f"{row['delta_dec']:.2f}°",
        str(row["epoch_label"]),
    ]


def _spice_year(readable) -> str:
    """Extract the year+era component from a datespice string."""
    parts = str(readable.datespice).split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return str(readable.datespice)


# ── result dataclasses ────────────────────────────────────────────────────────

@dataclass
class AlignmentResult:
    """Raw output from :func:`find_alignment_stars`."""
    ok: bool
    dec_target: float = 0.0
    stars_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    sky_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    midpoint_label: str = ""
    midpoint_jed: float = 0.0
    error: str = ""


@dataclass
class AlignmentPlotResult:
    """Plotly sky-map HTML ready to embed in PlotlyView."""
    ok: bool
    map_html: str = ""
    title: str = ""
    n_stars: int = 0
    dec_target: float = 0.0
    error: str = ""


# ── main computation ──────────────────────────────────────────────────────────

def find_alignment_stars(
    az: float          = DEFAULT_AZ,
    el: float          = DEFAULT_EL,
    lat: float         = DEFAULT_LAT,
    year_start: int    = DEFAULT_YEAR_START,
    era_start: str     = DEFAULT_ERA_START,
    year_end: int      = DEFAULT_YEAR_END,
    era_end: str       = DEFAULT_ERA_END,
    mag_limit: float   = DEFAULT_MAG_LIMIT,
    dec_tolerance: float = DEFAULT_DEC_TOL,
    n_epochs: int      = DEFAULT_N_EPOCHS,
) -> AlignmentResult:
    """Find stars that pass through direction (az, el) over a given year range.

    For each sampled epoch the stellar catalogue is precessed; stars whose
    epoch declination falls within ``dec_tolerance`` degrees of the target
    are collected. The table reports the epoch at which each star's declination
    was closest to the target.

    Parameters
    ----------
    az, el : float
        Azimuth and elevation of the alignment direction [degrees].
    lat : float
        Observer geographic latitude [degrees N].
    year_start, year_end : int
        Year range (historical BCE/CE counts, positive integers).
    era_start, era_end : str
        ``"bce"`` or ``"ce"``.
    mag_limit : float
        Faint magnitude limit (only stars with Vmag ≤ mag_limit are searched).
    dec_tolerance : float
        Half-width of the declination band [degrees].
    n_epochs : int
        Number of sample epochs uniformly distributed across the range.

    Returns
    -------
    AlignmentResult
    """
    try:
        montu = _import_montu()

        dec_target = compute_target_declination(az, el, lat)

        # Filter catalogue at J2000 for speed
        with timed_block("alignments: load catalogue"):
            allstars = montu.Stars()
            bright = allstars.get_stars(Vmag=[-3.0, float(mag_limit)])

        # Convert year bounds to JED
        t_start = montu.Time(
            _year_to_montu_str(year_start, era_start), calendar="proleptic"
        )
        t_end = montu.Time(
            _year_to_montu_str(year_end, era_end), calendar="proleptic"
        )
        jed_a = min(t_start.jed, t_end.jed)
        jed_b = max(t_start.jed, t_end.jed)

        n = max(3, int(n_epochs))
        jed_samples = np.linspace(jed_a, jed_b, n)
        jed_mid = 0.5 * (jed_a + jed_b)

        # Midpoint epoch object (for sky map background)
        t_mid = montu.Time(jed_mid, format="jd", calendar="proleptic")
        t_mid.get_readable()
        midpoint_label = _spice_year(t_mid.readable)

        # Precess at midpoint for the sky map background stars
        with timed_block("alignments: precess sky map"):
            sky_stars = bright.where_in_space(at=t_mid)
        sky_df = sky_stars.data.copy()

        # Search across all sample epochs
        # For each star keep the epoch where Δdec is smallest
        star_best: dict = {}   # original index → best record

        with timed_block("alignments: epoch search"):
            for jed in jed_samples:
                t = montu.Time(jed, format="jd", calendar="proleptic")
                t.get_readable()
                epoch_label = _spice_year(t.readable)

                precessed = bright.where_in_space(at=t)
                data = precessed.data

                for idx, row in data.iterrows():
                    dec_epoch = row.get("DecEpoch", np.nan)
                    if pd.isna(dec_epoch):
                        continue
                    delta = abs(float(dec_epoch) - dec_target)
                    prev = star_best.get(idx)
                    if prev is None or delta < prev["delta_dec"]:
                        star_best[idx] = {
                            "Name": str(row.get("Name", "")),
                            "ProperName": str(row.get("ProperName", "")),
                            "Vmag": float(row.get("Vmag", np.nan)),
                            "RAEpoch": float(row.get("RAEpoch", np.nan)),
                            "DecEpoch": float(dec_epoch),
                            "delta_dec": float(delta),
                            "epoch_label": epoch_label,
                            "jed": float(jed),
                        }

        aligned = sorted(
            [v for v in star_best.values() if v["delta_dec"] <= float(dec_tolerance)],
            key=lambda r: (r["Vmag"], r["delta_dec"]),
        )

        stars_df = (
            pd.DataFrame(aligned)
            if aligned
            else pd.DataFrame(
                columns=[
                    "Name", "ProperName", "Vmag",
                    "RAEpoch", "DecEpoch", "delta_dec",
                    "epoch_label", "jed",
                ]
            )
        )

        return AlignmentResult(
            ok=True,
            dec_target=dec_target,
            stars_df=stars_df,
            sky_df=sky_df,
            midpoint_label=midpoint_label,
            midpoint_jed=float(jed_mid),
        )

    except Exception as exc:
        return AlignmentResult(ok=False, error=str(exc))


# ── Plotly output ─────────────────────────────────────────────────────────────

def build_alignment_plots(
    result: AlignmentResult,
    az: float          = DEFAULT_AZ,
    el: float          = DEFAULT_EL,
    lat: float         = DEFAULT_LAT,
    observer_name: str = DEFAULT_OBSERVER_NAME,
    dec_tolerance: float = DEFAULT_DEC_TOL,
    mag_limit: float   = DEFAULT_MAG_LIMIT,
    year_start: int    = DEFAULT_YEAR_START,
    era_start: str     = DEFAULT_ERA_START,
    year_end: int      = DEFAULT_YEAR_END,
    era_end: str       = DEFAULT_ERA_END,
) -> AlignmentPlotResult:
    """Build the sky-map Plotly figure as embeddable HTML.

    Parameters
    ----------
    result : AlignmentResult
        Output from :func:`find_alignment_stars`.
    (rest)
        Original search parameters — used for chart titles and labels.

    Returns
    -------
    AlignmentPlotResult
    """
    if not result.ok:
        return AlignmentPlotResult(ok=False, error=result.error)

    try:
        dec_target  = result.dec_target
        stars_df    = result.stars_df
        sky_df      = result.sky_df

        era_range = (
            f"{_era_label(era_start, year_start)} – {_era_label(era_end, year_end)}"
        )
        title = (
            f"Star alignments  ·  az = {az:.1f}°  el = {el:.1f}°  "
            f"from {observer_name} (lat {lat:.2f}°)  ·  {era_range}"
        )

        map_fig = _build_alignment_sky_map(
            sky_df=sky_df,
            stars_df=stars_df,
            dec_target=dec_target,
            dec_tolerance=dec_tolerance,
            mag_limit=mag_limit,
            era_range=era_range,
            midpoint_label=result.midpoint_label,
            midpoint_jed=result.midpoint_jed,
            observer_name=observer_name,
            az=az, el=el, lat=lat,
            title=title,
        )
        map_html = figure_to_html(map_fig)

        return AlignmentPlotResult(
            ok=True,
            map_html=map_html,
            title=title,
            n_stars=len(stars_df),
            dec_target=dec_target,
        )

    except Exception as exc:
        return AlignmentPlotResult(ok=False, error=str(exc))


def _build_alignment_sky_map(
    *,
    sky_df: pd.DataFrame,
    stars_df: pd.DataFrame,
    dec_target: float,
    dec_tolerance: float,
    mag_limit: float,
    era_range: str,
    midpoint_label: str,
    midpoint_jed: float,
    observer_name: str,
    az: float,
    el: float,
    lat: float,
    title: str,
) -> go.Figure:
    """Base IAU sky map from :mod:`montu.stars` plus alignment overlays."""
    circumpolar_dec, pole = circumpolar_declination(lat)
    bg_bg, al, _ = _prepare_sky_frames(sky_df, stars_df)
    montu = _import_montu()
    mtime = montu.Time(midpoint_jed, format="jd", calendar="proleptic")

    fig = mercator_sky_map(
        bg_bg,
        mag_limit=mag_limit,
        label_bright_mag=min(2.5, mag_limit),
        at=mtime,
        show_constellation_boundaries=False,
    )

    fig.add_hrect(
        y0=dec_target - dec_tolerance,
        y1=dec_target + dec_tolerance,
        fillcolor="rgba(212, 175, 55, 0.13)",
        line_width=0,
    )
    fig.add_hline(
        y=dec_target,
        line_dash="dash",
        line_color="#d4af37",
        line_width=1.5,
        annotation_text=f"  target  {dec_target:+.2f}°",
        annotation_position="top left",
        annotation_font_color="#d4af37",
        annotation_font_size=11,
    )
    fig.add_hline(
        y=circumpolar_dec,
        line_dash="dot",
        line_color="#5eb3ff",
        line_width=1.5,
        annotation_text=f"  {_circumpolar_label(lat, circumpolar_dec, pole)}",
        annotation_position="bottom left",
        annotation_font_color="#5eb3ff",
        annotation_font_size=11,
    )

    if not al.empty:
        fig.add_trace(go.Scatter(
            x=al["ra_deg"],
            y=al["DecEpoch"],
            mode="markers+text",
            marker=dict(
                size=al["msize"],
                color="#ffd700",
                symbol="star",
                line=dict(width=1.2, color="#cc8800"),
            ),
            text=al["display_name"],
            textposition="top center",
            textfont=dict(size=12, color="#ffd700"),
            name="Aligned ★",
            customdata=np.stack(
                [al["Vmag"], al["delta_dec"], al["epoch_label"]], axis=1,
            ),
            hovertemplate=(
                "<b>%{text}</b><br>"
                "RA: %{x:.2f}°<br>"
                "Dec: %{y:.2f}°<br>"
                "V mag: %{customdata[0]:.2f}<br>"
                "Δ Dec: %{customdata[1]:.2f}°<br>"
                "Best year: %{customdata[2]}"
                "<extra></extra>"
            ),
            showlegend=True,
        ))

    subtitle = (
        f"az = {az:.1f}°  ·  el = {el:.1f}°  ·  "
        f"{observer_name} (lat {lat:.2f}°)  ·  {era_range}"
    )
    fig.update_layout(
        title=dict(
            text=(
                f"{title}<br>"
                f"<sup>Sky map — Mercator — epoch midpoint: {midpoint_label}"
                f"  ·  {subtitle}</sup>"
            ),
            x=0.5,
            xanchor="center",
            font=dict(size=13),
        ),
        margin=dict(l=60, r=40, t=110, b=60),
        autosize=True,
    )
    return fig
