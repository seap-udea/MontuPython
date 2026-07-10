"""
MontuPython star-alignment example.

Reproduces the search logic of the Star Alignments module using only the
``montu`` package.  The final section shows how to build a base IAU sky map
with :func:`montu.mercator_sky_map` and complement it with alignment overlays.

The King's Chamber northern shaft of the Great Pyramid of Khufu (Giza) points
almost due North at an elevation of ~31.7°. For an observer at Giza this
direction corresponds to a high northern declination. Around 2560 BCE the star
Thuban (α Draconis) passed through that declination window as precession
swept it toward the pole.

Run:
    python montu_star_alignments.py
"""

import montu
import numpy as np

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

# ─── Parameters (edit as needed) ─────────────────────────────────────────────

AZ = 0.0          # azimuth [°] from North, clockwise
EL = 31.7         # elevation [°] — King's Chamber north shaft
LAT = 29.9792     # observer latitude [°N] — Giza

YEAR_START = 2620
YEAR_END = 2420   # historical BCE years
MAG_LIMIT = 4.0   # Johnson V limiting magnitude
DEC_TOL = 1.0     # half-width of declination band [°]
N_EPOCHS = 7      # sample epochs across the year range


# ─── Target declination from azimuth, elevation, and latitude ────────────────

def compute_target_declination(az: float, el: float, lat: float) -> float:
    """Declination of the sky direction (az, el) at geographic latitude lat."""
    az_r = np.radians(az)
    el_r = np.radians(el)
    lat_r = np.radians(lat)
    sin_dec = (
        np.sin(lat_r) * np.sin(el_r)
        + np.cos(lat_r) * np.cos(el_r) * np.cos(az_r)
    )
    return float(np.degrees(np.arcsin(np.clip(sin_dec, -1.0, 1.0))))


def year_to_proleptic(year: int, era: str = "bce") -> str:
    """Build a proleptic Gregorian date string for mid-year."""
    y = max(1, int(year))
    if era.lower() == "bce":
        return f"bce {y:04d}-07-01 00:00:00"
    return f"{y:04d}-07-01 00:00:00"


def spice_year(readable) -> str:
    """Year label from ``Time.readable.datespice`` (e.g. ``2487 B.C.``)."""
    parts = str(readable.datespice).split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return str(readable.datespice)


def star_name(row) -> str:
    pn = str(row.get("ProperName", ""))
    if pn not in ("", "nan", "None"):
        return pn
    return str(row.get("Name", ""))


# ─── 1. Target declination ───────────────────────────────────────────────────

dec_target = compute_target_declination(AZ, EL, LAT)
print(f"\nAlignment direction:  az = {AZ:.1f}°  el = {EL:.1f}°")
print(f"Observer latitude:    {LAT:.4f}°N  (Giza)")
print(f"Target declination:   {dec_target:+.3f}°\n")

# ─── 2. Load bright stars and sample epochs ──────────────────────────────────

print(
    f"Searching {YEAR_START}–{YEAR_END} BCE for stars "
    f"brighter than V = {MAG_LIMIT} within ±{DEC_TOL}° of target dec …\n"
)

allstars = montu.Stars()
bright = allstars.get_stars(Vmag=[-3.0, MAG_LIMIT])

t_start = montu.Time(year_to_proleptic(YEAR_START, "bce"), calendar="proleptic")
t_end = montu.Time(year_to_proleptic(YEAR_END, "bce"), calendar="proleptic")
jed_a = min(t_start.jed, t_end.jed)
jed_b = max(t_start.jed, t_end.jed)
jed_samples = np.linspace(jed_a, jed_b, max(3, N_EPOCHS))

# For each star keep the epoch where |Dec − target| is smallest
star_best: dict = {}

for jed in jed_samples:
    mtime = montu.Time(jed, format="jd", calendar="proleptic")
    mtime.get_readable()
    epoch_label = spice_year(mtime.readable)

    precessed = bright.where_in_space(at=mtime)
    for idx, row in precessed.data.iterrows():
        dec_epoch = row.get("DecEpoch", np.nan)
        if np.isnan(dec_epoch):
            continue
        delta = abs(float(dec_epoch) - dec_target)
        prev = star_best.get(idx)
        if prev is None or delta < prev["delta_dec"]:
            star_best[idx] = {
                "name": star_name(row),
                "Vmag": float(row["Vmag"]),
                "RAEpoch": float(row["RAEpoch"]),
                "DecEpoch": float(dec_epoch),
                "delta_dec": float(delta),
                "epoch_label": epoch_label,
            }

aligned = sorted(
    [v for v in star_best.values() if v["delta_dec"] <= DEC_TOL],
    key=lambda r: (r["Vmag"], r["delta_dec"]),
)

# ─── 3. Print results (brightest first) ─────────────────────────────────────

if not aligned:
    print("No stars found within the tolerance band.")
else:
    print(f"Found {len(aligned)} aligned star(s):\n")
    for row in aligned:
        print(
            f"  {row['name']:<20}  V = {row['Vmag']:+.2f}  "
            f"RA = {row['RAEpoch']:.4f} h  "
            f"Dec = {row['DecEpoch']:+.3f}°  "
            f"Δ Dec = {row['delta_dec']:.3f}°  "
            f"({row['epoch_label']})"
        )
    print()

# ─── 4. Complement a base sky map with alignment overlays ────────────────────

if go is None:
    print("Install plotly to generate the sky map (pip install plotly).")
else:
    jed_mid = 0.5 * (jed_a + jed_b)
    t_mid = montu.Time(jed_mid, format="jd", calendar="proleptic")
    t_mid.get_readable()
    mid_label = spice_year(t_mid.readable)
    mid_sky = bright.where_in_space(at=t_mid).data

    # Base map from the package: IAU boundaries, asterisms, labels, stars
    fig = montu.mercator_sky_map(mid_sky, mag_limit=MAG_LIMIT, at=t_mid)

    # Caller adds alignment-specific layers
    fig.add_hrect(
        y0=dec_target - DEC_TOL,
        y1=dec_target + DEC_TOL,
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
    )
    circumpolar_dec = 90.0 - LAT
    fig.add_hline(
        y=circumpolar_dec,
        line_dash="dot",
        line_color="#5eb3ff",
        line_width=1.5,
        annotation_text=f"  circumpolar  dec ≥ {circumpolar_dec:+.1f}°",
        annotation_position="bottom left",
        annotation_font_color="#5eb3ff",
    )

    if aligned:
        ra_deg = [r["RAEpoch"] * 15.0 for r in aligned]
        decs = [r["DecEpoch"] for r in aligned]
        names = [r["name"] for r in aligned]
        vmags = [r["Vmag"] for r in aligned]
        deltas = [r["delta_dec"] for r in aligned]
        epochs = [r["epoch_label"] for r in aligned]
        sizes = [float(np.clip(13.0 - 2.0 * v, 3.0, 22.0)) * 1.9 for v in vmags]
        fig.add_trace(go.Scatter(
            x=ra_deg,
            y=decs,
            mode="markers+text",
            marker=dict(
                size=sizes,
                color="#ffd700",
                symbol="star",
                line=dict(width=1.2, color="#cc8800"),
            ),
            text=names,
            textposition="top center",
            textfont=dict(size=12, color="#ffd700"),
            name="Aligned ★",
            customdata=np.stack([vmags, deltas, epochs], axis=1),
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

    fig.update_layout(
        title=dict(
            text=(
                f"Khufu north shaft alignments  ·  epoch {mid_label}<br>"
                f"<sup>az = {AZ:.1f}°  ·  el = {EL:.1f}°  ·  "
                f"Giza (lat {LAT:.2f}°)  ·  "
                f"{YEAR_START}–{YEAR_END} BCE  ·  "
                f"found {len(aligned)} star(s)</sup>"
            ),
            x=0.5,
            xanchor="center",
            font=dict(size=13),
        ),
        margin=dict(l=60, r=40, t=110, b=60),
    )

    out_html = "star_alignment_map.html"
    fig.write_html(out_html)
    print(f"Sky map written to {out_html}\n")
