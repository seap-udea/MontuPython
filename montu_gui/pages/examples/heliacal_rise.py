# %pip install montu
"""Heliacal-rise calculation with MontuPython Desktop default values."""

import montu

# Default observer location in MontuPython Desktop: Thebes (Luxor).
observer = montu.Observer(lon=32.6422, lat=25.6967, height=0.076)

# Default body and interval: Sirius, third-apokatastasis window, mixed calendar.
sirius = montu.Stars(subset="bright", ProperName="Sirius")
start = montu.Time("133-06-01", calendar="mixed")
end = start + 10 * 365 * montu.DAY

# Default parameters for the Schaefer (1987) visibility model.
calculator = montu.HeliacalRise(
    model="schaefer1987",
    k=0.25,
    limiting_mag_zenith=6.0,
    sun_depression=-10.0,
)

events = calculator.compute(sirius, observer, start, end)
calculator.print_rises(
    events,
    title="Heliacal rises of Sirius",
    body_label="Sirius",
)

# Each row reports the detection instant, local time, body/Sun
# altitude and azimuth, and the model visibility quantities.
print(events)
