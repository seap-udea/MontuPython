# %pip install montu

import montu

# Observer at Thebes (Luxor)
## Define an observing site (lon/lat in degrees, height in km)
observer = montu.Observer(lon=32.6422, lat=25.6967, height=0.076)

print(f"Site: Thebes  ·  lat {observer.site.lat}°  ·  lon {observer.site.lon}°")
## Convert decimal degrees to sexagesimal notation
print(f"Sexagesimal: {montu.Util.dec2sex(observer.site.lat)}, {montu.Util.dec2sex(observer.site.lon)}")

# Sky at summer solstice, 1500 BCE
## Parse a date string into a Montu time object
mtime = montu.Time("-1500-06-21 12:00:00")
## Local solar time at the observer for that instant
print(f"\nLocal time: {observer.get_local_time(mtime)}")

for name in ("JUPITER", "VENUS"):
    ## Create a planet instance
    planet = montu.Planet(name)
    ## Compute altitude, azimuth and magnitude as seen from the observer
    planet.conditions_in_sky(at=mtime, observer=observer)
    print(planet)
