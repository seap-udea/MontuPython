"""
MontuPython orientation disk example.

Reproduces the logic of the Orientation Disk module using only the
``montu`` package: extreme rise and set azimuths for celestial bodies
over a multi-year window.

Run:
    python montu_orientation_disk.py
"""

import montu
import math
import numpy as np

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

# ─── Parameters (edit as needed) ─────────────────────────────────────────────

YEAR = 2560
ERA = "bce"          # "bce" or "ce"
SPAN_YEARS = 3
STEP_DAYS = 5
HORIZON_EL = 0.0     # effective horizon altitude [°]

OBSERVER = montu.Observer(lon=31.1342, lat=29.9792, height=0.075)  # Giza
BODIES = ["Sun", "Moon", "Venus"]


# ─── Helpers ─────────────────────────────────────────────────────────────────

def year_to_jed(year: int, era: str) -> float:
    y = max(1, int(year))
    if era.lower() == "bce":
        date_str = f"bce {y:04d}-01-01 00:00:00"
    else:
        date_str = f"{y:04d}-01-01 00:00:00"
    return montu.Time(date_str, calendar="proleptic").jed


def make_body(name: str):
    if name.lower() == "sun":
        return montu.Sun()
    if name.lower() == "moon":
        return montu.Moon()
    return montu.Planet(name)


def rise_set_az(body, observer, jed: float, horizon_el: float = 0.0):
    """Return (rise_az°, set_az°) at one epoch, or (None, None) on failure."""
    import ephem as pyephem

    site = pyephem.Observer()
    site.lon = observer.site.lon
    site.lat = observer.site.lat
    site.elevation = observer.site.elevation
    site.pressure = observer.site.pressure
    site.temp = observer.site.temp
    site.horizon = str(horizon_el)
    site.date = float(jed) - montu.PYEPHEM_JD_REF

    rise_az = set_az = None
    try:
        rd = site.next_rising(body.seba)
        rs = pyephem.Observer()
        rs.lon, rs.lat = site.lon, site.lat
        rs.elevation, rs.pressure, rs.temp = site.elevation, site.pressure, site.temp
        rs.horizon = str(horizon_el)
        rs.date = rd
        body.seba.compute(rs)
        rise_az = math.degrees(float(body.seba.az))
    except (pyephem.AlwaysUpError, pyephem.NeverUpError):
        pass

    try:
        site2 = pyephem.Observer()
        site2.lon, site2.lat = site.lon, site.lat
        site2.elevation, site2.pressure, site2.temp = site.elevation, site.pressure, site.temp
        site2.horizon = str(horizon_el)
        site2.date = site.date
        sd = site2.next_setting(body.seba)
        ss = pyephem.Observer()
        ss.lon, ss.lat = site.lon, site.lat
        ss.elevation, ss.pressure, ss.temp = site.elevation, site.pressure, site.temp
        ss.horizon = str(horizon_el)
        ss.date = sd
        body.seba.compute(ss)
        set_az = math.degrees(float(body.seba.az))
    except (pyephem.AlwaysUpError, pyephem.NeverUpError):
        pass

    return rise_az, set_az


# ─── 1. Sweep azimuths ───────────────────────────────────────────────────────

jed_start = year_to_jed(YEAR, ERA)
jed_end = jed_start + SPAN_YEARS * 365.25
jed_steps = np.arange(jed_start, jed_end, STEP_DAYS)

print(f"\nOrientation disk  ·  {YEAR} {ERA.upper()}  ·  {OBSERVER.site.lat}°N")
print(f"Window: {SPAN_YEARS} years  ·  step {STEP_DAYS} days\n")

for name in BODIES:
    body = make_body(name)
    rise_azs, set_azs = [], []

    for jed in jed_steps:
        body.seba.compute(OBSERVER.site)
        r_az, s_az = rise_set_az(body, OBSERVER, jed, HORIZON_EL)
        if r_az is not None:
            rise_azs.append(r_az)
        if s_az is not None:
            set_azs.append(s_az)

    if not rise_azs and not set_azs:
        print(f"{name:8s}  — no rise/set events in window")
        continue

    rise_n = min(rise_azs) if rise_azs else float("nan")
    rise_s = max(rise_azs) if rise_azs else float("nan")
    set_s  = min(set_azs)  if set_azs  else float("nan")
    set_n  = max(set_azs)  if set_azs  else float("nan")
    print(
        f"{name:8s}  rise {rise_n:6.1f}°–{rise_s:6.1f}°   "
        f"set {set_s:6.1f}°–{set_n:6.1f}°"
    )

# ─── 2. Optional polar plot ──────────────────────────────────────────────────

if go is not None:
    fig = go.Figure()
    colors = {"Sun": "#B71C1C", "Moon": "#1565C0", "Venus": "#C62828"}

    for name in BODIES:
        body = make_body(name)
        rise_azs, set_azs = [], []
        for jed in jed_steps:
            r_az, s_az = rise_set_az(body, OBSERVER, jed, HORIZON_EL)
            if r_az is not None:
                rise_azs.append(r_az)
            if s_az is not None:
                set_azs.append(s_az)
        if not rise_azs:
            continue
        color = colors.get(name, "#5eb3ff")
        for az in sorted(set([min(rise_azs), max(rise_azs)])):
            fig.add_trace(go.Scatterpolar(
                r=[0, 0.9], theta=[az, az], mode="lines",
                line=dict(color=color, width=2), name=f"{name} rise",
            ))
            fig.add_trace(go.Scatterpolar(
                r=[0.95], theta=[az], mode="markers",
                marker=dict(symbol="triangle-up", size=10, color=color),
                showlegend=False,
            ))

    fig.update_layout(
        polar=dict(
            angularaxis=dict(direction="clockwise", rotation=90),
            radialaxis=dict(visible=False, range=[0, 1.2]),
        ),
        title=f"Orientation disk — {YEAR} {ERA.upper()} — Giza",
        height=640,
    )
    fig.show()
