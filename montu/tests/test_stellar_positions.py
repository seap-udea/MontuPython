"""Compare Spica rise/transit/set times against Stellarium reference data.

Reference rows live in ``test-stellar-positions.csv``: the first year of
each century from 1600 to 2000, sampled on the 1st day of every third
month (January, April, July, October).  Site: 25.696699°N, 32.642°E,
76 m (Thebes-like).  Star: Spica (Stellarium export).

MontuPython precesses the catalogue star with PyEphem and searches for the
nearest rise, transit, and set around each calendar date.

Tolerances
----------
- Rise, transit, and set times: <= 1 minute vs Stellarium
- Transit altitude: <= 1 arcminute

Stellarium rise/set in the reference export use the geometric horizon (no
atmospheric refraction).  The observer therefore uses ``pressure=0`` and
``temperature=0`` so PyEphem matches that convention.
"""

from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import ephem
import pytest

import montu
from montu.stars import Star

STELLAR_CSV = Path(__file__).with_name("test-stellar-positions.csv")
TIME_TOL_MIN = 1.0
ALT_TOL_ARCMIN = 1.0
ARCSEC_PER_ARCMIN = 60.0

OBSERVER_LON = 32.642
OBSERVER_LAT = 25.696699
OBSERVER_HEIGHT_KM = 0.076
STAR_NAME = "Spica"


def _parse_time(value: str) -> datetime:
    return datetime.strptime(value.strip(), "%Y-%m-%d %H:%M:%S")


def _time_diff_minutes(computed: datetime, reference: datetime) -> float:
    return abs((computed - reference).total_seconds()) / 60.0


def _to_datetime(jed: float) -> datetime:
    text = montu.Time(jed, format="jd", calendar="proleptic").readable.datepro
    return _parse_time(text.split(".")[0])


def _load_rows() -> list[dict[str, str]]:
    with STELLAR_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _row_id(row: dict[str, str]) -> str:
    return f"{row['Name'].strip()}-{row['Calendar_Date'].strip()}"


def _collect_horizon_events(star: Star, observer: montu.Observer, calendar_date: str) -> dict[str, list[float]]:
    year, month, day = map(int, calendar_date.split("-"))
    site = ephem.Observer()
    site.lon = observer.site.lon
    site.lat = observer.site.lat
    site.elevation = observer.site.elevation
    site.pressure = observer.site.pressure
    site.temp = observer.site.temp

    body = star.seba
    window_start = ephem.Date(f"{year}/{month}/{day}") - 1
    window_end = ephem.Date(f"{year}/{month}/{day}") + 2

    events: dict[str, list[float]] = {"rise": [], "transit": [], "set": []}
    mapping = {
        "rise": site.next_rising,
        "transit": site.next_transit,
        "set": site.next_setting,
    }
    for kind, finder in mapping.items():
        site.date = window_start
        try:
            instant = finder(body)
            while float(instant) <= float(window_end):
                events[kind].append(float(instant) + montu.PYEPHEM_JD_REF)
                site.date = instant + 1e-8
                instant = finder(body)
        except (ephem.AlwaysUpError, ephem.NeverUpError):
            continue
    return events


def _nearest_event(reference: datetime, event_jeds: list[float]) -> float:
    ref_jed = montu.Time(
        reference.strftime("%Y-%m-%d %H:%M:%S"),
        calendar="proleptic",
    ).jed
    return min(event_jeds, key=lambda jed: abs(jed - ref_jed))


def _transit_altitude_deg(star: Star, observer: montu.Observer, transit_jed: float) -> float:
    site = ephem.Observer()
    site.lon = observer.site.lon
    site.lat = observer.site.lat
    site.elevation = observer.site.elevation
    site.pressure = observer.site.pressure
    site.temp = observer.site.temp
    site.date = transit_jed - montu.PYEPHEM_JD_REF
    star.seba.compute(site)
    return float(star.seba.alt) * montu.RAD


@pytest.fixture(scope="module")
def thebes_like_observer():
    return montu.Observer(
        lon=OBSERVER_LON,
        lat=OBSERVER_LAT,
        height=OBSERVER_HEIGHT_KM,
        pressure=0,
        temperature=0,
    )


@pytest.fixture(scope="module")
def spica_star():
    catalogue = montu.Stars(subset="bright", ProperName=STAR_NAME)
    assert catalogue.number == 1
    return Star(catalogue.data.iloc[0])


@pytest.mark.parametrize("row", _load_rows(), ids=_row_id)
def test_spica_horizon_events_match_stellarium(row, thebes_like_observer, spica_star):
    calendar_date = row["Calendar_Date"].strip()
    ref_rise = _parse_time(row["Rise"])
    ref_transit = _parse_time(row["Transit"])
    ref_set = _parse_time(row["Set"])
    ref_alt = float(row["Altitude_deg"])

    events = _collect_horizon_events(spica_star, thebes_like_observer, calendar_date)
    for kind, reference in (("rise", ref_rise), ("transit", ref_transit), ("set", ref_set)):
        assert events[kind], f"No MontuPython {kind} events near {calendar_date}"
        nearest_jed = _nearest_event(reference, events[kind])
        computed = _to_datetime(nearest_jed)
        diff = _time_diff_minutes(computed, reference)
        assert diff <= TIME_TOL_MIN, (
            f"{STAR_NAME} {kind} on {calendar_date}: "
            f"Stellarium={reference}, Montu={computed}, diff={diff:.3f} min "
            f"(tolerance {TIME_TOL_MIN} min)"
        )

    transit_jed = _nearest_event(ref_transit, events["transit"])
    altitude = _transit_altitude_deg(spica_star, thebes_like_observer, transit_jed)
    alt_diff_arcmin = abs(altitude - ref_alt) * ARCSEC_PER_ARCMIN
    assert alt_diff_arcmin <= ALT_TOL_ARCMIN, (
        f"{STAR_NAME} transit altitude on {calendar_date}: "
        f"Stellarium={ref_alt:.5f}°, Montu={altitude:.5f}°, "
        f"diff={alt_diff_arcmin:.3f} arcmin"
    )
