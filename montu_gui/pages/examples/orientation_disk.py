# %pip install montu plotly

import math

import montu
import plotly.graph_objects as go

YEAR, ERA = 2560, "bce"
SPAN_YEARS, STEP_DAYS, HORIZON_EL = 3, 5, 0.0
## Observer at Giza (lon/lat in degrees, height in km)
observer = montu.Observer(lon=31.1342, lat=29.9792, height=0.075)
BODIES = ["Sun", "Moon", "Venus"]


def make_body(name):
    if name.lower() == "sun":
        ## Solar disk for rise/set calculations
        return montu.Sun()
    if name.lower() == "moon":
        ## Lunar disk for rise/set calculations
        return montu.Moon()
    ## Planet instance by name
    return montu.Planet(name)


def rise_set_az(body, jed):
    """Rise and set azimuths at one Julian day (uses PyEphem via montu)."""
    import ephem as pyephem

    site = pyephem.Observer()
    site.lon = observer.site.lon
    site.lat = observer.site.lat
    site.elevation = observer.site.elevation
    site.horizon = str(HORIZON_EL)
    site.date = float(jed) - montu.PYEPHEM_JD_REF

    rise_az = set_az = None
    try:
        rd = site.next_rising(body.seba)
        rs = pyephem.Observer()
        rs.lon, rs.lat = site.lon, site.lat
        rs.elevation = site.elevation
        rs.horizon = str(HORIZON_EL)
        rs.date = rd
        body.seba.compute(rs)
        rise_az = math.degrees(float(body.seba.az))
    except (pyephem.AlwaysUpError, pyephem.NeverUpError):
        pass

    try:
        site.date = float(jed) - montu.PYEPHEM_JD_REF
        sd = site.next_setting(body.seba)
        ss = pyephem.Observer()
        ss.lon, ss.lat = site.lon, site.lat
        ss.elevation = site.elevation
        ss.horizon = str(HORIZON_EL)
        ss.date = sd
        body.seba.compute(ss)
        set_az = math.degrees(float(body.seba.az))
    except (pyephem.AlwaysUpError, pyephem.NeverUpError):
        pass

    return rise_az, set_az


date_str = f"bce {YEAR:04d}-01-01 00:00:00" if ERA == "bce" else f"{YEAR:04d}-01-01 00:00:00"
## Julian day at the start of the reference year
jed_start = montu.Time(date_str, calendar="proleptic").jed
jed_end = jed_start + SPAN_YEARS * 365.25

print(f"Orientation disk  ·  {YEAR} {ERA.upper()}  ·  Giza")

for name in BODIES:
    body = make_body(name)
    rise_azs, set_azs = [], []
    jed = jed_start
    while jed < jed_end:
        ## Update the body's ephemeris for the observer site
        body.seba.compute(observer.site)
        r_az, s_az = rise_set_az(body, jed)
        if r_az is not None:
            rise_azs.append(r_az)
        if s_az is not None:
            set_azs.append(s_az)
        jed += STEP_DAYS

    if rise_azs or set_azs:
        print(
            f"{name:8s}  rise {min(rise_azs, default=0):6.1f}°–{max(rise_azs, default=0):6.1f}°   "
            f"set {min(set_azs, default=0):6.1f}°–{max(set_azs, default=0):6.1f}°"
        )

# Polar plot of extreme azimuths
fig = go.Figure()
colors = {"Sun": "#B71C1C", "Moon": "#1565C0", "Venus": "#C62828"}
for name in BODIES:
    body = make_body(name)
    rise_azs = []
    jed = jed_start
    while jed < jed_end:
        r_az, _ = rise_set_az(body, jed)
        if r_az is not None:
            rise_azs.append(r_az)
        jed += STEP_DAYS
    if not rise_azs:
        continue
    color = colors.get(name, "#5eb3ff")
    for az in (min(rise_azs), max(rise_azs)):
        fig.add_trace(go.Scatterpolar(
            r=[0, 0.9], theta=[az, az], mode="lines",
            line=dict(color=color, width=2), name=f"{name}",
        ))

fig.update_layout(
    polar=dict(angularaxis=dict(direction="clockwise", rotation=90), radialaxis=dict(visible=False)),
    title=f"Orientation disk — {YEAR} {ERA.upper()}",
    height=640,
)
fig.show()
