"""
MontuPython calendar conversion example.

Reproduces the conversions shown in the Calendar Calculator module:
  - Gregorian / Julian ↔ Caniucular (Egyptian civil)
  - Julian Day, ephemeris time, Delta-T

Run:
    python montu_calendar_conversion.py
"""

from datetime import date

import montu

# ─── Example 1: Any date → Caniucular ──────────────────────────────────────
# Convert today's date to all calendar representations,
# including the Egyptian civil (caniucular) calendar.

t = montu.Time("2026-07-10 00:00:00", calendar="mixed")

print(f"─── Example 1: {date} → Caniucular ───")
print("Gregorian proleptic (SPICE): ", t.readable.datespice)
print("Gregorian proleptic (astron):", t.readable.datepro)
print("Mixed Julian/Gregorian:      ", t.readable.datemix)
print("Caniucular (Egyptian civil): ", t.readable.datecan)
print("Julian Day (UTC):            ", t.jed)
print("TT ephemeris seconds J2000:  ", t.tt)
print("UTC ephemeris seconds J2000: ", t.et)
print("Delta-T (seconds):           ", t.deltat)

# ─── Example 2: today → Caniucular ─────────────────────────────────────────
# Convert today's date to all calendar representations,
# including the Egyptian civil (caniucular) calendar.

today = date.today()
today_montu = f"{today.year}-{today.month:02d}-{today.day:02d} 00:00:00"
t = montu.Time(today_montu, calendar="mixed")

print(f"─── Example 1: {today.strftime('%d %B %Y')} (today) → Caniucular ───")
print("Gregorian proleptic (SPICE): ", t.readable.datespice)
print("Gregorian proleptic (astron):", t.readable.datepro)
print("Mixed Julian/Gregorian:      ", t.readable.datemix)
print("Caniucular (Egyptian civil): ", t.readable.datecan)
print("Julian Day (UTC):            ", t.jed)
print("TT ephemeris seconds J2000:  ", t.tt)
print("UTC ephemeris seconds J2000: ", t.et)
print("Delta-T (seconds):           ", t.deltat)


# ─── Example 3: First Apokatastasis → Gregorian ──────────────────────────────
# The First Apokatastasis is the epoch of the caniucular (Egyptian civil)
# calendar. It marks Horus year 0, I-Akhet-1 — the moment when the heliacal
# rising of Sirius (Sopedet) coincided with New Year's Day of the civil year.
# Historical equivalent: BCE 2782-07-20 (Mixed Julian/Gregorian).

t_epoch = montu.Time("hrw 0-I-Akhet-1", calendar="caniucular")

print("\n─── Example 2: First Apokatastasis (epoch of the caniucular calendar) ───")
print("Caniucular:                  ", t_epoch.readable.datecan)
print("Mixed Julian/Gregorian:      ", t_epoch.readable.datemix)
print("Gregorian proleptic (SPICE): ", t_epoch.readable.datespice)
print("Julian Day (UTC):            ", t_epoch.jed)
print("Delta-T (seconds):           ", t_epoch.deltat)
