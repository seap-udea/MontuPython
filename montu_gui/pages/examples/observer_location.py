"""
MontuPython observer location example.

Shows how to define an observing site with ``montu.Observer`` and use it
in sky calculations (local time, planet altitude/azimuth, etc.).

Run:
    python montu_observer_location.py
"""

import json
from pathlib import Path

import montu

# ─── Method 1: explicit coordinates (decimal degrees + height in km) ────────

SITE_NAME = "Thebes (Luxor)"
LAT = 25.6967          # degrees north (positive)
LON = 32.6422          # degrees east (positive)
ALT_M = 76.0           # metres above sea level

observer = montu.Observer(
    lon=LON,
    lat=LAT,
    height=ALT_M / 1000.0,   # MontuPython expects kilometres
)

print(f"Observer: {SITE_NAME}")
print(f"  Latitude:  {LAT}°")
print(f"  Longitude: {LON}°")
print(f"  Height:    {ALT_M} m  ({ALT_M / 1000:.3f} km)")

# ─── Method 2: load a predefined ancient site from locations.json ───────────

# Path when this file lives in montu_gui/pages/examples/ inside the repo:
LOCATIONS_FILE = (
    Path(__file__).resolve().parent.parent.parent / "assets" / "locations.json"
)

if LOCATIONS_FILE.exists():
    with open(LOCATIONS_FILE, encoding="utf-8") as fh:
        catalog = json.load(fh)
    site_id = catalog.get("default", "thebes")
    entry = next(
        loc for loc in catalog["locations"] if loc["id"] == site_id
    )
    observer = montu.Observer(
        lon=entry["lon"],
        lat=entry["lat"],
        height=entry["alt_m"] / 1000.0,
    )
    SITE_NAME = entry["name"]
    LAT, LON, ALT_M = entry["lat"], entry["lon"], entry["alt_m"]
    print(f"\nLoaded from JSON ({LOCATIONS_FILE.name}): {SITE_NAME}")

# ─── Sexagesimal (DMS) display ───────────────────────────────────────────────

print("\nCoordinates in sexagesimal notation:")
print(f"  Lat {montu.Util.dec2hex(LAT)}")
print(f"  Lon {montu.Util.dec2hex(LON)}")

# ─── Use the observer in sky calculations ───────────────────────────────────

DATE = "-1500-06-21 12:00:00"   # near summer solstice, 1500 BCE
mtime = montu.Time(DATE)

print(f"\n─── Sky at {SITE_NAME} on {DATE} ───")
print(f"Local solar time: {observer.get_local_time(mtime)}")

jupiter = montu.Planet("JUPITER")
jupiter.conditions_in_sky(at=mtime, observer=observer)
p = jupiter.position
c = jupiter.condition
print(
    f"Jupiter: altitude {montu.Util.dec2hex(p.el)}°, "
    f"azimuth {montu.Util.dec2hex(p.az)}°, "
    f"magnitude {c.Vmag:.1f}"
)

venus = montu.Planet("VENUS")
venus.conditions_in_sky(at=mtime, observer=observer)
p = venus.position
c = venus.condition
print(
    f"Venus:   altitude {montu.Util.dec2hex(p.el)}°, "
    f"azimuth {montu.Util.dec2hex(p.az)}°, "
    f"magnitude {c.Vmag:.1f}"
)

print(
    "\nPass the same ``observer`` to ``planet.conditions_in_sky``, "
    "``Sun.next_seasons``, ephemeris loops, etc."
)
