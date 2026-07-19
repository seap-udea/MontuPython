"""Shared helpers for MontuPython organic regression snapshots.

Organic CSV files freeze MontuPython outputs at a given package version so
future releases can be compared against the same reference epoch grid.
Regenerate with ``bin/generate_organic_planetary_ephemeris.py`` and
``bin/generate_organic_stellar_positions.py`` when catalogue coordinates or
ephemeris algorithms change intentionally.
"""

from __future__ import annotations

import csv
from datetime import date
from io import StringIO
from pathlib import Path
from typing import Iterable

import montu
from montu.stars import Star

TESTS_DIR = Path(__file__).resolve().parent
PLANETARY_ORGANIC_CSV = TESTS_DIR / "test-planetary-ephemeris-organic.csv"
STELLAR_ORGANIC_CSV = TESTS_DIR / "test-stellar-positions-organic.csv"

OBSERVER_LON = 32.642
OBSERVER_LAT = 25.696699
OBSERVER_HEIGHT_KM = 0.076

PLANETARY_BODIES = ("Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn")
STELLAR_NAMES = ("Sirius", "Canopus")
# Quarterly samples within each Saros-selected year (every 3 months).
PLANETARY_SAMPLES_PER_YEAR = (1, 4, 7, 10)
SAROS_DAYS = 6585.321
STELLAR_YEAR_STEP = 100
YEAR_START_BCE = 3000
PLANETARY_YEAR_END_CE = 1599
YEAR_END_CE = 1600

METADATA_KEYS = ("generated", "montu_version", "description")


def package_version() -> str:
    return montu.version


def snapshot_metadata(description: str, *, row_count: int | None = None) -> dict[str, str]:
    import ephem

    metadata = {
        "generated": date.today().isoformat(),
        "montu_version": package_version(),
        "pyephem_version": ephem.__version__,
        "description": description,
    }
    if row_count is not None:
        metadata["row_count"] = str(row_count)
    return metadata


def historical_era_year(at: montu.Time) -> tuple[str, int]:
    """Return historical BCE/CE era and 1-based year from a proleptic Time."""
    at.get_readable()
    astro_year = int(at.readable.year)
    if astro_year <= 0:
        return "bce", 1 - astro_year
    return "ce", astro_year


def planetary_saros_sample_years() -> list[tuple[str, int]]:
    """Historical years sampled one per Saros cycle from 3000 BCE to 1599 CE."""
    start = montu.Time("bce 3000-01-01 12:00:00", calendar="proleptic")
    end_jed = montu.Time(
        f"{PLANETARY_YEAR_END_CE:04d}-12-31 12:00:00",
        calendar="proleptic",
    ).jed

    years: list[tuple[str, int]] = []
    seen: set[tuple[str, int]] = set()
    jed = start.jed
    while jed <= end_jed:
        era, year = historical_era_year(montu.Time(jed, format="jd", calendar="proleptic"))
        key = (era, year)
        if key not in seen:
            seen.add(key)
            years.append(key)
        jed += SAROS_DAYS

    final = ("ce", PLANETARY_YEAR_END_CE)
    if final not in seen:
        years.append(final)
    return years


def stellar_sample_years() -> list[tuple[str, int]]:
    years = [( "bce", year) for year in range(YEAR_START_BCE, 0, -STELLAR_YEAR_STEP)]
    years += [("ce", year) for year in range(STELLAR_YEAR_STEP, YEAR_END_CE + 1, STELLAR_YEAR_STEP)]
    return years


def format_proleptic_datetime(era: str, year: int, month: int, day: int, hour: int = 12) -> str:
    if era == "bce":
        astro_year = 1 - year
    else:
        astro_year = year
    return f"{astro_year:04d}-{month:02d}-{day:02d} {hour:02d}:00:00"


def calendar_date_only(era: str, year: int, month: int, day: int) -> str:
    if era == "bce":
        astro_year = 1 - year
    else:
        astro_year = year
    return f"{astro_year:04d}-{month:02d}-{day:02d}"


def jed_to_proleptic_datetime(jed: float) -> str:
    text = montu.Time(jed, format="jd", calendar="proleptic").readable.datepro
    return text.split(".")[0]


def geocentric_observer() -> montu.Observer:
    return montu.Observer(lon=0, lat=0, height=0)


def thebes_like_observer() -> montu.Observer:
    return montu.Observer(
        lon=OBSERVER_LON,
        lat=OBSERVER_LAT,
        height=OBSERVER_HEIGHT_KM,
        pressure=0,
        temperature=0,
    )


def make_planetary_body(name: str):
    if name == "Sun":
        return montu.Sun()
    if name == "Moon":
        return montu.Moon()
    return montu.Planet(name)


def load_star(name: str) -> tuple[Star, montu.Stars]:
    catalogue = montu.Stars(subset="bright", ProperName=name)
    if catalogue.number != 1:
        raise RuntimeError(f"Expected one catalogue row for {name}, got {catalogue.number}")
    row = catalogue.data.iloc[0]
    return Star(row), catalogue


def proper_motion_j2000(star_row, at: montu.Time) -> tuple[float, float]:
    dt = (at.jed - montu.JED_2000) * montu.MARCSEC / 365.25
    raj2000t = float(star_row["RAJ2000"]) + float(star_row["pmRA"]) * dt / 15.0
    decj2000t = float(star_row["DecJ2000"]) + float(star_row["pmDec"]) * dt
    return raj2000t, decj2000t


def compute_planetary_row(name: str, date_time: str, observer: montu.Observer) -> dict[str, str]:
    mtime = montu.Time(date_time, calendar="proleptic")
    body = make_planetary_body(name)
    body.conditions_in_sky(at=mtime, observer=observer)
    return {
        "Name": name,
        "Date and Time": date_time,
        "RAJ2000": f"{body.position.RAJ2000:.8f}",
        "DecJ2000": f"{body.position.DecJ2000:.8f}",
        "Mag": f"{body.condition.Vmag:.5f}",
        "Phase": f"{body.condition.phase:.5f}",
        "Dist_AU": f"{body.condition.earth_distance:.8f}",
    }


def compute_stellar_row(
    name: str,
    era: str,
    year: int,
    star: Star,
    star_row,
    observer: montu.Observer,
) -> dict[str, str]:
    sample_time = format_proleptic_datetime(era, year, 1, 1, hour=12)
    mtime = montu.Time(sample_time, calendar="proleptic")
    star.conditions_in_sky(at=mtime, observer=observer)
    raj2000t, decj2000t = proper_motion_j2000(star_row, mtime)
    return {
        "Name": name,
        "Historical_Year": str(year),
        "Era": era,
        "Calendar_Date": calendar_date_only(era, year, 1, 1),
        "Sample_Time": sample_time,
        "RAJ2000": f"{star.position.RAJ2000:.8f}",
        "DecJ2000": f"{star.position.DecJ2000:.8f}",
        "RAJ2000t": f"{raj2000t:.8f}",
        "DecJ2000t": f"{decj2000t:.8f}",
        "RAEpoch": f"{star.position.RAEpoch:.8f}",
        "DecEpoch": f"{star.position.DecEpoch:.8f}",
        "az": f"{star.position.az:.8f}",
        "el": f"{star.position.el:.8f}",
        "Rise": jed_to_proleptic_datetime(star.condition.rise_time),
        "Transit": jed_to_proleptic_datetime(star.condition.transit_time),
        "Set": jed_to_proleptic_datetime(star.condition.set_time),
        "Transit_el_deg": f"{star.condition.transit_el:.8f}",
        "Lon_deg": f"{OBSERVER_LON:.6f}",
        "Lat_deg": f"{OBSERVER_LAT:.6f}",
        "Height_m": f"{OBSERVER_HEIGHT_KM * 1000:.1f}",
    }


def write_metadata_lines(metadata: dict[str, str]) -> list[str]:
    return [f"# {key}={value}" for key, value in metadata.items()]


def write_organic_csv(
    path: Path,
    fieldnames: list[str],
    rows: Iterable[dict[str, str]],
    metadata: dict[str, str],
) -> None:
    buffer = StringIO()
    for line in write_metadata_lines(metadata):
        buffer.write(line + "\n")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path.write_text(buffer.getvalue(), encoding="utf-8")


def read_organic_csv(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    metadata: dict[str, str] = {}
    data_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, _, value = line[1:].partition("=")
            metadata[key.strip()] = value.strip()
            continue
        data_lines.append(line)
    reader = csv.DictReader(data_lines)
    rows = list(reader)
    return metadata, rows


def generate_planetary_rows() -> list[dict[str, str]]:
    observer = geocentric_observer()
    rows: list[dict[str, str]] = []
    for era, year in planetary_saros_sample_years():
        for month in PLANETARY_SAMPLES_PER_YEAR:
            date_time = format_proleptic_datetime(era, year, month, 1, hour=12)
            for name in PLANETARY_BODIES:
                rows.append(compute_planetary_row(name, date_time, observer))
    return rows


def generate_stellar_rows() -> list[dict[str, str]]:
    observer = thebes_like_observer()
    stars = {name: load_star(name) for name in STELLAR_NAMES}
    rows: list[dict[str, str]] = []
    for era, year in stellar_sample_years():
        for name in STELLAR_NAMES:
            star, catalogue = stars[name]
            rows.append(
                compute_stellar_row(
                    name,
                    era,
                    year,
                    star,
                    catalogue.data.iloc[0],
                    observer,
                )
            )
    return rows
