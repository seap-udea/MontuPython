# %pip install montu plotly

import math

import montu
import plotly.graph_objects as go

YEAR, ERA = 2560, "bce"
SPAN_YEARS, STEP_DAYS, HORIZON_EL = 3, 5, 0.0
## Observador en Giza (lon/lat en grados, altura en km)
observer = montu.Observer(lon=31.1342, lat=29.9792, height=0.075)
BODIES = ["Sun", "Moon", "Venus"]

COLORS = {
    "Sun": "#B71C1C",
    "Moon": "#1565C0",
    "Venus": "#C62828",
}


def make_body(name):
    if name.lower() == "sun":
        return montu.Sun()
    if name.lower() == "moon":
        return montu.Moon()
    return montu.Planet(name)


def _copy_site(site, *, date=None, horizon_el=0.0):
    import ephem as pyephem

    s = pyephem.Observer()
    s.lon = site.lon
    s.lat = site.lat
    s.elevation = site.elevation
    s.pressure = site.pressure
    s.temp = site.temp
    s.horizon = str(horizon_el)
    s.date = date if date is not None else site.date
    return s


def rise_set_az(seba, site, horizon_el=0.0):
    """Azimuts de salida y puesta en la epoca del observador ya definida en *site*."""
    import ephem as pyephem

    rise_az = set_az = None
    circumpolar = neverup = False

    try:
        rd = _copy_site(site, horizon_el=horizon_el).next_rising(seba)
        rs = _copy_site(site, date=rd, horizon_el=horizon_el)
        seba.compute(rs)
        rise_az = math.degrees(float(seba.az))
    except pyephem.AlwaysUpError:
        circumpolar = True
    except pyephem.NeverUpError:
        neverup = True

    try:
        sd = _copy_site(site, horizon_el=horizon_el).next_setting(seba)
        ss = _copy_site(site, date=sd, horizon_el=horizon_el)
        seba.compute(ss)
        set_az = math.degrees(float(seba.az))
    except pyephem.AlwaysUpError:
        circumpolar = True
    except pyephem.NeverUpError:
        neverup = True

    seba.compute(site)
    return rise_az, set_az, circumpolar, neverup


def compute_extremes(name):
    """Azimuts extremos norte/sur de salida y puesta en la ventana de busqueda."""
    body = make_body(name)
    rise_azs, set_azs = [], []
    circumpolar = neverup = False

    jed = jed_start
    while jed < jed_end:
        observer.site.date = float(jed) - montu.PYEPHEM_JD_REF
        body.seba.compute(observer.site)
        r_az, s_az, cp, nu = rise_set_az(body.seba, observer.site, HORIZON_EL)
        if cp:
            circumpolar = True
        if nu:
            neverup = True
        if r_az is not None:
            rise_azs.append(r_az)
        if s_az is not None:
            set_azs.append(s_az)
        jed += STEP_DAYS

    return {
        "name": name,
        "color": COLORS.get(name, "#1565C0"),
        "rise_north": min(rise_azs) if rise_azs else None,
        "rise_south": max(rise_azs) if rise_azs else None,
        "set_north": max(set_azs) if set_azs else None,
        "set_south": min(set_azs) if set_azs else None,
        "circumpolar": circumpolar,
        "neverup": neverup,
    }


def _az_label(name, suffix):
    abbrev = name[:3] if len(name) > 3 else name
    return f"{abbrev}{suffix}"


def build_disk_figure(extremes, *, year, era, observer_name="Giza", lat=29.9792):
    """Disco polar de orientacion con la misma disposicion de MontuPython Desktop."""
    fig = go.Figure()
    theta_ring = list(range(361))

    fig.add_trace(go.Scatterpolar(
        r=[1.0] * 361, theta=theta_ring, mode="lines",
        line=dict(color="rgba(80,80,80,0.4)", width=1.5),
        hoverinfo="skip", showlegend=False,
    ))
    fig.add_trace(go.Scatterpolar(
        r=[0.5] * 361, theta=theta_ring, mode="lines",
        line=dict(color="rgba(180,180,180,0.3)", width=0.8, dash="dot"),
        hoverinfo="skip", showlegend=False,
    ))
    for az_pair in ((0, 180), (90, 270)):
        fig.add_trace(go.Scatterpolar(
            r=[0, 1.1, None, 0, 1.1],
            theta=[az_pair[0], az_pair[0], None, az_pair[1], az_pair[1]],
            mode="lines",
            line=dict(color="rgba(100,100,100,0.35)", width=1),
            hoverinfo="skip", showlegend=False,
        ))

    has_legend = False
    for body in extremes:
        if body["circumpolar"] or body["neverup"]:
            continue

        color = body["color"]
        name = body["name"]
        arrows = []
        if body["rise_north"] is not None:
            arrows.append((body["rise_north"], True, "△N"))
        if body["rise_south"] is not None and body["rise_south"] != body["rise_north"]:
            arrows.append((body["rise_south"], True, "△S"))
        if body["set_north"] is not None:
            arrows.append((body["set_north"], False, "▽N"))
        if body["set_south"] is not None and body["set_south"] != body["set_north"]:
            arrows.append((body["set_south"], False, "▽S"))

        first = True
        for az, is_rise, suffix in arrows:
            symbol = "triangle-up" if is_rise else "triangle-down"
            fig.add_trace(go.Scatterpolar(
                r=[0.0, 0.90], theta=[az, az], mode="lines",
                line=dict(color=color, width=2.2),
                showlegend=first and not has_legend,
                name=name, hoverinfo="skip",
            ))
            fig.add_trace(go.Scatterpolar(
                r=[0.95], theta=[az], mode="markers",
                marker=dict(symbol=symbol, size=12, color=color,
                            line=dict(width=1.2, color="white")),
                showlegend=first, name=name if first else "",
                legendgroup=name,
                hovertemplate=(
                    f"<b>{name}</b><br>"
                    f"{'Rise' if is_rise else 'Set'} az: {az:.1f}°<extra></extra>"
                ),
            ))
            fig.add_trace(go.Scatterpolar(
                r=[1.12], theta=[az], mode="text",
                text=[f"<b>{_az_label(name, suffix)}</b>"],
                textfont=dict(size=9, color=color),
                showlegend=False, hoverinfo="skip",
            ))
            if first:
                has_legend = True
                first = False

    for az, lbl in (
        (0, "<b>N</b>"), (45, "NE"), (90, "<b>E</b>"), (135, "SE"),
        (180, "<b>S</b>"), (225, "SW"), (270, "<b>W</b>"), (315, "NW"),
    ):
        fig.add_trace(go.Scatterpolar(
            r=[1.25], theta=[az], mode="text", text=[lbl],
            textfont=dict(size=15 if "<b>" in lbl else 11, color="rgba(60,60,60,0.85)"),
            showlegend=False, hoverinfo="skip",
        ))

    era_label = f"{year} BCE" if era.lower() == "bce" else f"{year} CE"
    fig.update_layout(
        polar=dict(
            angularaxis=dict(
                direction="clockwise",
                rotation=90,
                tickmode="array",
                tickvals=list(range(0, 360, 10)),
                ticktext=[f"{d}°" if d % 30 == 0 else "" for d in range(0, 360, 10)],
                tickfont=dict(size=9, color="rgba(80,80,80,0.8)"),
                linecolor="rgba(80,80,80,0.4)",
                gridcolor="rgba(180,180,180,0.3)",
            ),
            radialaxis=dict(visible=False, range=[0, 1.4]),
            bgcolor="rgba(240, 245, 255, 0.4)",
        ),
        title=dict(
            text=(
                f"⭕  Orientation Disk  ·  {era_label}  ·  "
                f"{observer_name} (lat {lat:.2f}°)<br><sup>"
                "Arrows show extreme rise (△) and set (▽) azimuths  ·  "
                "N = northernmost  ·  S = southernmost</sup>"
            ),
            x=0.5, xanchor="center", font=dict(size=13),
        ),
        legend=dict(orientation="h", yanchor="bottom", y=-0.08,
                    xanchor="center", x=0.5),
        margin=dict(l=60, r=60, t=110, b=80),
        height=640,
    )
    return fig


date_str = (
    f"bce {YEAR:04d}-01-01 00:00:00" if ERA == "bce"
    else f"{YEAR:04d}-01-01 00:00:00"
)
jed_start = montu.Time(date_str, calendar="proleptic").jed
jed_end = jed_start + SPAN_YEARS * 365.25

print(f"Orientation disk  ·  {YEAR} {ERA.upper()}  ·  Giza")

extremes = []
for name in BODIES:
    body = compute_extremes(name)
    extremes.append(body)
    if body["rise_north"] is not None or body["set_north"] is not None:
        print(
            f"{name:8s}  rise "
            f"{body['rise_north'] or 0:6.1f}°–{body['rise_south'] or 0:6.1f}°   "
            f"set {body['set_south'] or 0:6.1f}°–{body['set_north'] or 0:6.1f}°"
        )

fig = build_disk_figure(
    extremes,
    year=YEAR,
    era=ERA,
    observer_name="Giza",
    lat=observer.lat,
)
fig.show()
