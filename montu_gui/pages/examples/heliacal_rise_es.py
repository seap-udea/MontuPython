# %pip install montu
"""Calculo de levantamientos heliacos con valores por defecto de MontuPython Desktop."""

import montu

# Ubicación por defecto del observador en MontuPython Desktop: Tebas (Luxor).
observer = montu.Observer(lon=32.6422, lat=25.6967, height=0.076)

# Cuerpo e intervalo por defecto: Sirio, ventana de primera apocatastasis, calendario mixto.
sirius = montu.Stars(subset="bright", ProperName="Sirius")
start = montu.Time("133-06-01", calendar="mixed")
end = start + 10 * 365 * montu.DAY

# Parametros por defecto del modelo de visibilidad Schaefer (1987).
calculator = montu.HeliacalRise(
    model="schaefer1987",
    k=0.25,
    limiting_mag_zenith=6.0,
    sun_depression=-10.0,
)

events = calculator.compute(sirius, observer, start, end)
calculator.print_rises(
    events,
    title="Levantamientos heliacos de Sirio",
    body_label="Sirius",
)

# Cada fila reporta el instante de deteccion, hora local, altitud/azimut
# del cuerpo y del Sol, y las cantidades de visibilidad del modelo.
print(events)
