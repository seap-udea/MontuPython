# %pip install montu
"""Solar eclipse catalogue search and local circumstances with MontuPython."""

import montu

# ---------------------------------------------------------------------------
# 1. Find eclipses in the NASA Five Millennium catalogue with get_eclipses()
# ---------------------------------------------------------------------------
# Historical years use astronomical numbering in the catalogue:
# 585 BCE = year -584; 600 BCE = year -599.

eclipses = montu.SolarEclipses()

# Thales' eclipse: 28 May 585 BCE (NASA catalogue year -584-05-28).
thales_day = eclipses.get_eclipses(year=-584, month=5, day=28)
print(f"Eclipses on 28 May 585 BCE: {thales_day.number}")

eclipse = montu.SolarEclipse(thales_day.data.iloc[0])
print(eclipse)

# Filter a date range and eclipse types (tuple = OR of catalogue codes).
window = eclipses.get_eclipses(
    year=[-599, -499],
    eclipse_type=("T", "Tm", "Ts", "Tn"),
)
print(f"Total eclipses in catalogue, 600–500 BCE: {window.number}")

# ---------------------------------------------------------------------------
# 2. Local circumstances at an observing site with conditions_eclipse()
# ---------------------------------------------------------------------------
troy = montu.Observer(site="troy")
cond = eclipse.conditions_eclipse(troy)

print(f"\nThales' eclipse at Troy ({troy.lat:.4f}°N, {troy.lon:.4f}°E)")
print(f"  Local kind      : {cond.kind}")
print(f"  Visible         : {cond.visible}")
print(f"  Magnitude       : {cond.magnitude:.3f}")
print(f"  Sun altitude    : {cond.sun_altitude_deg:.1f}°")
print(f"  Maximum (UT)    : {montu.Time(cond.jed_max, format='jd', scale='utc').readable.datepro}")
print(f"  Maximum (local) : {troy.get_local_time(cond.time_max)}")

# Scan the window and keep eclipses visible at Troy.
visible_at_troy = []
for _, row in window.data.iterrows():
    candidate = montu.SolarEclipse(row)
    local = candidate.conditions_eclipse(troy)
    if local.visible:
        y = int(row.year)
        visible_at_troy.append(
            (
                1 - y,  # historical BCE year
                int(row.month),
                int(row.day),
                local.kind,
                round(local.magnitude, 3),
            )
        )

print(f"\nTotal eclipses visible at Troy, 600–500 BCE: {len(visible_at_troy)}")
for entry in visible_at_troy:
    print(f"  {entry[0]} BCE {entry[1]:02d}-{entry[2]:02d}  {entry[3]:8s}  mag={entry[4]:.3f}")
