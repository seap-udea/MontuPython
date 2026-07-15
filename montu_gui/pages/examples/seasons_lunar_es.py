# %pip install montu

import montu

YEAR = -1500

# Estaciones del año
## Instante de referencia al inicio del año
t_start = montu.Time(f"{YEAR}-01-01 12:00:00", calendar="mixed")
## Dias julianos de los proximos equinoccios y solsticios despues de t_start
vernal, summer, autumnal, _ = montu.Sun.next_seasons(at=t_start)
_, _, _, winter = montu.Sun.next_seasons(
    at=montu.Time(autumnal, format="jd", calendar="mixed")
)

print(f"Estaciones para {YEAR}")
for label, jed in (
    ("Equinoccio de primavera", vernal),
    ("Solsticio de verano", summer),
    ("Equinoccio de otono", autumnal),
    ("Solsticio de invierno", winter),
):
    ## Reconstruir un objeto Time completo a partir de cada dia juliano
    mt = montu.Time(jed, format="jd", calendar="mixed", full=True)
    print(f"\n{label}")
    print(mt)

# Cuartos lunares en el mismo año
t_end = montu.Time(f"{YEAR + 1}-01-01 12:00:00", calendar="mixed")
## Lista de instantes de luna nueva/creciente/llena/menguante desde t_start
quarters = montu.Moon.next_moon_quarters(
    since=t_start,
    starting_at="new",
    numquarters=56,
    output="mtime",
    format="columns",
)

print(f"\nFases lunares para {YEAR}")
for item in quarters:
    mt = item["Datetime"]
    if mt.jed >= t_end.jed:
        break
    print(f"\n{item['Quarter']}")
    print(mt)
