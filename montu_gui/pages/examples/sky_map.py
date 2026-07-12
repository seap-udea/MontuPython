# %pip install montu plotly

import montu
from montu.maps import LINE_ECLIPTIC, LINE_HORIZON, polar_sky_map

# Observer at Giza; local evening on 1 January 2500 BCE
DATE = "bce 2500-01-01"
LOCAL_HOUR, LOCAL_MINUTE, LOCAL_SECOND = 18, 0, 0
MAG_LIMIT = 3.5
CONSTELLATION_SET = "egyptian_ancient"  # iau | egyptian_ancient | egyptian_dendera
BODIES = ["Sun", "Moon"]
LINES = [LINE_ECLIPTIC, LINE_HORIZON]

observer = montu.Observer(lon=31.1342, lat=29.9792, height=0.075)

## Precess the stellar catalogue to the map epoch (noon on DATE)
mtime = montu.Time(f"{DATE} 12:00:00", calendar="proleptic")
stars = montu.Stars().where_in_space(at=mtime)

fig_north, fig_south = polar_sky_map(
    DATE,
    local_hour=LOCAL_HOUR,
    local_minute=LOCAL_MINUTE,
    local_second=LOCAL_SECOND,
    observer=observer,
    mag_limit=MAG_LIMIT,
    bodies=BODIES,
    lines=LINES,
    meridian_view=False,
    constellation_set=CONSTELLATION_SET,
    observer_name="Giza",
    precessed_star_data=stars.data,
)

print(
    f"Sky map · Giza · {DATE} "
    f"{LOCAL_HOUR:02d}:{LOCAL_MINUTE:02d}:{LOCAL_SECOND:02d} local · "
    f"V ≤ {MAG_LIMIT} · {CONSTELLATION_SET}"
)

fig_north.show()
fig_south.show()
