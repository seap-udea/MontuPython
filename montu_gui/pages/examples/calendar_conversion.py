# %pip install montu

import montu

# Gregorian date -> all calendar formats
## Analyze a mixed Julian/Gregorian date
t = montu.Time("2026-07-10 00:00:00", calendar="mixed")
print("2026-07-10")
print(t)

# Another historical date
## BCE dates use the convention "bce YYYY-MM-DD"
t = montu.Time("bce 1500-06-21 12:00:00", calendar="mixed")
print(f"\nbce 1500-06-21")
print(t)

# Sothic epoch (first apokatastasis)
## Egyptian civil-calendar date (Horus year, month, day)
t = montu.Time("[hrw 0] I akhet 1", calendar="sothic")
print(f"\n[hrw 0] I akhet 1")
print(t)

# Historical date (Canopus Decree)
## MontuPython includes predefined dates (same catalogue as Calendar -> Historical dates)
historical_dates = montu.load_historical_dates()
date_key = "bce 238-03-07"
entry = historical_dates[date_key]
t = montu.Time(date_key, calendar="mixed")
print(f"\n{entry['label']}")
print(f"  Julian/Gregorian (mixed): {date_key}")
print(f"  Known civil date: {entry['egyptian_date']}")
print(f"  Computed civil date: {t.readable.datesot}")
print(t)
