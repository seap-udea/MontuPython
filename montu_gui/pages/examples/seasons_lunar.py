"""
MontuPython astronomical seasons & lunar phases example.

Reproduces the calculations shown in the Seasons & Lunar Phases module:
  - Exact dates of the four astronomical seasons for this year
  - All lunar quarter phases (new moon, first quarter, full moon, last quarter)

Dates are shown in:
  - Gregorian proleptic (SPICE style)
  - Mixed Julian/Gregorian
  - Caniucular (Egyptian civil)

Run:
    python montu_seasons_lunar.py
"""

from datetime import date

import montu

YEAR = date.today().year   # change to any year, e.g. YEAR = -1000

# ─── Astronomical Seasons ────────────────────────────────────────────────────

# Start from Jan 1 of the requested year; next_seasons() returns the
# Julian Days of the four upcoming astronomical events.

t_start = montu.Time(f"{YEAR}-01-01 12:00:00", calendar="mixed")
vernal, summer, autumnal, _ = montu.Sun.next_seasons(at=t_start)
_, _, _, winter = montu.Sun.next_seasons(
    at=montu.Time(autumnal, format="jd", calendar="mixed")
)

SEASON_NAMES = [
    "Northward equinox (Spring)",
    "Northern solstice (Summer)",
    "Southward equinox (Autumn)",
    "Southern solstice (Winter)",
]

print(f"─── Astronomical Seasons for year {YEAR} ───")
print(f"{'Season':<32} {'Gregorian proleptic':<32} {'Mixed Julian/Greg.':<26} Caniucular")
print("─" * 120)
for name, jed in zip(SEASON_NAMES, [vernal, summer, autumnal, winter]):
    mt = montu.Time(jed, format="jd", calendar="mixed", full=True)
    print(
        f"{name:<32} "
        f"{mt.readable.datespice:<32} "
        f"{mt.readable.datemix:<26} "
        f"{mt.readable.datecan}"
    )

# ─── Lunar Phases ────────────────────────────────────────────────────────────

# Fetch ~56 quarters (≈14 months) starting at the first new moon of the year,
# then keep only those whose date falls within the requested year.

PHASE_ICONS = {"new": "🌑", "first": "🌓", "full": "🌕", "last": "🌗"}
PHASE_NAMES = {"new": "New Moon", "first": "First Quarter",
               "full": "Full Moon", "last": "Last Quarter"}

t_end = montu.Time(f"{YEAR + 1}-01-01 12:00:00", calendar="mixed")

quarters = montu.Moon.next_moon_quarters(
    since=t_start,
    starting_at="new",
    numquarters=56,
    output="mtime",
    format="columns",
)

print(f"\n─── Lunar Phases for year {YEAR} ───")
print(f"{'':3} {'Phase':<16} {'Mixed Julian/Greg.':<26} Caniucular")
print("─" * 75)
for item in quarters:
    mt = item["Datetime"]
    if mt.jed >= t_end.jed:
        break
    q = item["Quarter"]
    print(
        f"{PHASE_ICONS.get(q, '')} "
        f"{PHASE_NAMES.get(q, q):<16} "
        f"{mt.readable.datemix:<26} "
        f"{mt.readable.datecan}"
    )
