# %pip install montu

import montu

# Observador en Tebas (Luxor)
## Define un sitio de observacion (lon/lat en grados, altura en km)
observer = montu.Observer(lon=32.6422, lat=25.6967, height=0.076)

print(f"Sitio: Tebas  ·  lat {observer.site.lat}°  ·  lon {observer.site.lon}°")
## Convertir grados decimales a notacion sexagesimal
print(f"Sexagesimal: {montu.Util.dec2sex(observer.site.lat)}, {montu.Util.dec2sex(observer.site.lon)}")

# Cielo en el solsticio de verano, 1500 AEC
## Analizar una fecha a un objeto de tiempo de Montu
mtime = montu.Time("-1500-06-21 12:00:00")
## Hora solar local del observador en ese instante
print(f"\nHora local: {observer.get_local_time(mtime)}")

for name in ("JUPITER", "VENUS"):
    ## Crear una instancia de planeta
    planet = montu.Planet(name)
    ## Calcular altitud, azimut y magnitud vistos por el observador
    planet.conditions_in_sky(at=mtime, observer=observer)
    print(planet)
