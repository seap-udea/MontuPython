# %pip install montu

import montu

# Gregorian date → all calendar formats
## Parse a mixed Julian/Gregorian date
t = montu.Time("2026-07-10 00:00:00", calendar="mixed")
print("2026-07-10")
print(t)

# Another historical date
## BCE dates use the "bce YYYY-MM-DD" convention
t = montu.Time("bce 1500-06-21 12:00:00", calendar="mixed")
print(f"\nbce 1500-06-21")
print(t)

# Caniucular epoch (First Apokatastasis)
## Egyptian civil calendar date (Horus year, month, day)
t = montu.Time("hrw 0-I-Akhet-1", calendar="caniucular")
print(f"\nhrw 0-I-Akhet-1")
print(t)
