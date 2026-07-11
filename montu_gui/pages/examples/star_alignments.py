# %pip install montu plotly

import math

import montu
import plotly.graph_objects as go

# Khufu north shaft: az 0°, el 31.7°, observer at Giza
AZ, EL, LAT = 0.0, 31.7, 29.9792
YEAR_START, YEAR_END = 2620, 2420
MAG_LIMIT, DEC_TOL, N_EPOCHS = 4.0, 1.0, 7

# Target declination of the sight line
az_r, el_r, lat_r = map(math.radians, (AZ, EL, LAT))
sin_dec = math.sin(lat_r) * math.sin(el_r) + math.cos(lat_r) * math.cos(el_r) * math.cos(az_r)
dec_target = math.degrees(math.asin(max(-1.0, min(1.0, sin_dec))))
print(f"Target declination: {dec_target:+.3f}°")

# Bright stars and sample epochs
## Get all stars
allstars = montu.Stars()
## From all stars filter those with magnitudes below MAG_LIMIT
bright = allstars.get_stars(Vmag=[-3.0, MAG_LIMIT])

## Start and end of the search window
t_start = montu.Time(f"bce {YEAR_START:04d}-07-01 00:00:00", calendar="proleptic")
t_end = montu.Time(f"bce {YEAR_END:04d}-07-01 00:00:00", calendar="proleptic")
jed_a, jed_b = min(t_start.jed, t_end.jed), max(t_start.jed, t_end.jed)
jed_step = (jed_b - jed_a) / max(1, N_EPOCHS - 1)

star_best = {}
for i in range(N_EPOCHS):
    ## Build a time object for each sampled Julian day
    mtime = montu.Time(jed_a + i * jed_step, format="jd", calendar="proleptic")
    ## Fill mtime.readable with dates in every calendar format (for labels)
    mtime.get_readable()
    epoch_label = mtime.readable.datespice.split()[:2]
    epoch_label = f"{epoch_label[0]} {epoch_label[1]}" if len(epoch_label) >= 2 else str(mtime.readable.datespice)

    ## Precess the star catalogue to this epoch
    precessed = bright.where_in_space(at=mtime)
    for idx, row in precessed.data.iterrows():
        dec_epoch = row.get("DecEpoch")
        if dec_epoch is None or (isinstance(dec_epoch, float) and math.isnan(dec_epoch)):
            continue
        delta = abs(float(dec_epoch) - dec_target)
        prev = star_best.get(idx)
        if prev is None or delta < prev["delta_dec"]:
            name = row.get("ProperName") or row.get("Name")
            star_best[idx] = {
                "name": str(name),
                "Vmag": float(row["Vmag"]),
                "RAEpoch": float(row["RAEpoch"]),
                "DecEpoch": float(dec_epoch),
                "delta_dec": delta,
                "epoch_label": epoch_label,
            }

aligned = [v for v in star_best.values() if v["delta_dec"] <= DEC_TOL]
aligned.sort(key=lambda r: (r["Vmag"], r["delta_dec"]))

print(f"\nAligned stars ({YEAR_START}–{YEAR_END} BCE):")
for row in aligned:
    print(
        f"  {row['name']:<16}  V={row['Vmag']:+.1f}  "
        f"Dec={row['DecEpoch']:+.2f}°  Δ={row['delta_dec']:.2f}°  ({row['epoch_label']})"
    )

# Sky map at mid-epoch
jed_mid = 0.5 * (jed_a + jed_b)
t_mid = montu.Time(jed_mid, format="jd", calendar="proleptic")
## Fill t_mid.readable with dates in every calendar format (for the map title)
t_mid.get_readable()
## Star positions at the mid-point of the search window
mid_sky = bright.where_in_space(at=t_mid).data

## Base Mercator sky map with IAU boundaries and stars
fig = montu.mercator_sky_map(mid_sky, mag_limit=MAG_LIMIT, at=t_mid)
fig.add_hrect(
    y0=dec_target - DEC_TOL, y1=dec_target + DEC_TOL,
    fillcolor="rgba(212, 175, 55, 0.13)", line_width=0,
)
fig.add_hline(y=dec_target, line_dash="dash", line_color="#d4af37")

if aligned:
    fig.add_trace(go.Scatter(
        x=[r["RAEpoch"] * 15.0 for r in aligned],
        y=[r["DecEpoch"] for r in aligned],
        mode="markers+text",
        text=[r["name"] for r in aligned],
        textposition="top center",
        marker=dict(size=12, color="#ffd700", symbol="star"),
        name="Aligned",
    ))

fig.show()
