# %pip install montu

import montu

YEAR = -1500

# Seasons of the year
## Reference instant at the start of the year
t_start = montu.Time(f"{YEAR}-01-01 12:00:00", calendar="mixed")
## Julian days of the next equinoxes and solstices after t_start
vernal, summer, autumnal, _ = montu.Sun.next_seasons(at=t_start)
_, _, _, winter = montu.Sun.next_seasons(
    at=montu.Time(autumnal, format="jd", calendar="mixed")
)

print(f"Seasons for {YEAR}")
for label, jed in (
    ("Spring equinox", vernal),
    ("Summer solstice", summer),
    ("Autumn equinox", autumnal),
    ("Winter solstice", winter),
):
    ## Rebuild a full Time object from each Julian day
    mt = montu.Time(jed, format="jd", calendar="mixed", full=True)
    print(f"\n{label}")
    print(mt)

# Lunar quarters in the same year
t_end = montu.Time(f"{YEAR + 1}-01-01 12:00:00", calendar="mixed")
## List of new/first/full/last moon times starting at t_start
quarters = montu.Moon.next_moon_quarters(
    since=t_start,
    starting_at="new",
    numquarters=56,
    output="mtime",
    format="columns",
)

print(f"\nLunar phases for {YEAR}")
for item in quarters:
    mt = item["Datetime"]
    if mt.jed >= t_end.jed:
        break
    print(f"\n{item['Quarter']}")
    print(mt)
