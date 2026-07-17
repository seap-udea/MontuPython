# %pip install montu

import montu

# Fecha gregoriana -> todos los formatos de calendario
## Analizar una fecha juliana/gregoriana mixta
t = montu.Time("2026-07-10 00:00:00", calendar="mixed")
print("2026-07-10")
print(t)

# Otra fecha historica
## Las fechas AEC usan la convencion "bce YYYY-MM-DD"
t = montu.Time("bce 1500-06-21 12:00:00", calendar="mixed")
print(f"\nbce 1500-06-21")
print(t)

# Epoca sothic (primera apocatastasis)
## Fecha del calendario civil egipcio (año de Horus, mes, dia)
t = montu.Time("[hrw 0] I akhet 1", calendar="sothic")
print(f"\n[hrw 0] I akhet 1")
print(t)

# Fecha historica (Decreto de Canopo)
## MontuPython incluye fechas predefinidas (mismo catalogo de Calendario -> Fechas historicas)
historical_dates = montu.load_historical_dates()
date_key = "bce 238-03-07"
entry = historical_dates[date_key]
t = montu.Time(date_key, calendar="mixed")
print(f"\n{entry['label']}")
print(f"  Juliano/Gregoriano (mixto): {date_key}")
print(f"  Fecha civil conocida: {entry['egyptian_date']}")
print(f"  Fecha civil calculada: {t.readable.datesot}")
print(t)
