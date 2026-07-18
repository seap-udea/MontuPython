"""Compare MontuPython planet ephemerides against Stellarium reference values.

Reference data live in ``test-planetary-ephemeris.csv`` (geocentric J2000 equatorial
coordinates, illuminated phase, geocentric distance, heliocentric elongation,
and lunar elongation exported from Stellarium).

Solar elongation checks (planet vs Sun) validate the planet ephemeris together
with the Sun.  Lunar elongation checks (Moon vs Sun) validate the Moon at the
same epoch.

Tolerances
----------
- Equatorial position: angular separation <= 1 arcmin
- Solar and lunar elongation: <= 1 arcmin (Stellarium reports positive degrees;
  PyEphem may use the opposite sign, so comparisons use ``abs``)
- Phase: relative difference <= 10 %
- Geocentric distance: relative difference <= 1 %
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

import montu

EPHEMERIS_CSV = Path(__file__).with_name("test-planetary-ephemeris.csv")
ARCSEC_PER_ARCMIN = 60.0
POSITION_TOL_ARCMIN = 1.0
ELONGATION_TOL_ARCMIN = 1.0
REL_TOL_PHASE = 0.10
REL_TOL_DISTANCE = 0.01


def _parse_deg(value: str) -> float:
    return float(str(value).replace("°", "").strip())


def _parse_percent(value: str) -> float:
    return float(str(value).replace("%", "").strip())


def _relative_difference(computed: float, reference: float) -> float:
    if reference == 0:
        return abs(computed - reference)
    return abs(computed - reference) / abs(reference)


def _angular_separation_deg(
    ra_deg_a: float,
    dec_deg_a: float,
    ra_deg_b: float,
    dec_deg_b: float,
) -> float:
    ra_a, dec_a, ra_b, dec_b = map(
        math.radians,
        (ra_deg_a, dec_deg_a, ra_deg_b, dec_deg_b),
    )
    cos_delta = (
        math.sin(dec_a) * math.sin(dec_b)
        + math.cos(dec_a) * math.cos(dec_b) * math.cos(ra_b - ra_a)
    )
    cos_delta = min(1.0, max(-1.0, cos_delta))
    return math.degrees(math.acos(cos_delta))


def _angular_separation_arcmin(
    ra_deg_a: float,
    dec_deg_a: float,
    ra_deg_b: float,
    dec_deg_b: float,
) -> float:
    return _angular_separation_deg(ra_deg_a, dec_deg_a, ra_deg_b, dec_deg_b) * ARCSEC_PER_ARCMIN


def _elongation_difference_arcmin(computed_deg: float, reference_deg: float) -> float:
    """Compare elongations regardless of PyEphem/Stellarium sign convention."""
    return abs(abs(computed_deg) - abs(reference_deg)) * ARCSEC_PER_ARCMIN


def _load_ephemeris_rows() -> list[dict[str, str]]:
    with EPHEMERIS_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _row_id(row: dict[str, str]) -> str:
    return f"{row['Name'].strip()}-{row['Date and Time'].strip()}"


@pytest.fixture(scope="module")
def geocentric_observer():
    """Geocentric reference frame (Stellarium exports geocentric ephemerides)."""
    return montu.Observer(lon=0, lat=0, height=0)


@pytest.mark.parametrize("row", _load_ephemeris_rows(), ids=_row_id)
def test_planet_ephemeris_matches_stellarium(row, geocentric_observer):
    name = row["Name"].strip()
    date_time = row["Date and Time"].strip()

    mtime = montu.Time(date_time, calendar="proleptic")
    planet = montu.Planet(name)
    planet.conditions_in_sky(at=mtime, observer=geocentric_observer)

    sun = montu.Sun()
    sun.where_in_sky(at=mtime, observer=geocentric_observer)

    moon = montu.Moon()
    moon.conditions_in_sky(at=mtime, observer=geocentric_observer)

    ra_ref = _parse_deg(row["RA (J2000)"])
    dec_ref = _parse_deg(row["Dec (J2000)"])
    phase_ref = _parse_percent(row["Phase"])
    dist_ref = float(row["Dist_AU"].strip())
    solar_elong_ref = _parse_deg(row["Elong."])
    lunar_elong_ref = _parse_deg(row["Lunar_Elong."])

    ra_deg = planet.position.RAJ2000 * 15.0
    dec_deg = planet.position.DecJ2000
    phase = planet.condition.phase
    dist_au = planet.condition.earth_distance

    sun_ra_deg = sun.position.RAJ2000 * 15.0
    sun_dec_deg = sun.position.DecJ2000
    moon_ra_deg = moon.position.RAJ2000 * 15.0
    moon_dec_deg = moon.position.DecJ2000

    separation_arcmin = _angular_separation_arcmin(ra_ref, dec_ref, ra_deg, dec_deg)
    assert separation_arcmin <= POSITION_TOL_ARCMIN, (
        f"{name} @ {date_time}: separation {separation_arcmin:.3f} arcmin "
        f"(RA/Dec ref={ra_ref:.5f}°/{dec_ref:.5f}°, "
        f"montu={ra_deg:.5f}°/{dec_deg:.5f}°)"
    )

    assert _relative_difference(phase, phase_ref) <= REL_TOL_PHASE, (
        f"{name} @ {date_time}: phase ref={phase_ref:.2f}%, montu={phase:.2f}%, "
        f"relative diff={_relative_difference(phase, phase_ref):.3%}"
    )
    assert _relative_difference(dist_au, dist_ref) <= REL_TOL_DISTANCE, (
        f"{name} @ {date_time}: distance ref={dist_ref:.6f} AU, montu={dist_au:.6f} AU, "
        f"relative diff={_relative_difference(dist_au, dist_ref):.3%}"
    )

    solar_elong_planet = _elongation_difference_arcmin(
        planet.condition.elongation, solar_elong_ref,
    )
    assert solar_elong_planet <= ELONGATION_TOL_ARCMIN, (
        f"{name} @ {date_time}: solar elongation ref={solar_elong_ref:.5f}°, "
        f"montu={planet.condition.elongation:.5f}°, diff={solar_elong_planet:.3f} arcmin"
    )

    solar_elong_from_coords = _elongation_difference_arcmin(
        _angular_separation_deg(sun_ra_deg, sun_dec_deg, ra_deg, dec_deg),
        solar_elong_ref,
    )
    assert solar_elong_from_coords <= ELONGATION_TOL_ARCMIN, (
        f"{name} @ {date_time}: Sun-planet separation ref={solar_elong_ref:.5f}°, "
        f"from Montu coords={_angular_separation_deg(sun_ra_deg, sun_dec_deg, ra_deg, dec_deg):.5f}°, "
        f"diff={solar_elong_from_coords:.3f} arcmin"
    )

    solar_elong_sun_vs_ref_planet = _elongation_difference_arcmin(
        _angular_separation_deg(sun_ra_deg, sun_dec_deg, ra_ref, dec_ref),
        solar_elong_ref,
    )
    assert solar_elong_sun_vs_ref_planet <= ELONGATION_TOL_ARCMIN, (
        f"{name} @ {date_time}: Sun vs Stellarium planet ref elongation "
        f"ref={solar_elong_ref:.5f}°, "
        f"sep(Sun, ref planet)="
        f"{_angular_separation_deg(sun_ra_deg, sun_dec_deg, ra_ref, dec_ref):.5f}°, "
        f"diff={solar_elong_sun_vs_ref_planet:.3f} arcmin"
    )

    lunar_elong_moon = _elongation_difference_arcmin(
        moon.condition.elongation, lunar_elong_ref,
    )
    assert lunar_elong_moon <= ELONGATION_TOL_ARCMIN, (
        f"{name} @ {date_time}: lunar elongation ref={lunar_elong_ref:.5f}°, "
        f"montu={moon.condition.elongation:.5f}°, diff={lunar_elong_moon:.3f} arcmin"
    )

    lunar_elong_from_coords = _elongation_difference_arcmin(
        _angular_separation_deg(moon_ra_deg, moon_dec_deg, sun_ra_deg, sun_dec_deg),
        lunar_elong_ref,
    )
    assert lunar_elong_from_coords <= ELONGATION_TOL_ARCMIN, (
        f"{name} @ {date_time}: Moon-Sun separation ref={lunar_elong_ref:.5f}°, "
        f"from Montu coords="
        f"{_angular_separation_deg(moon_ra_deg, moon_dec_deg, sun_ra_deg, sun_dec_deg):.5f}°, "
        f"diff={lunar_elong_from_coords:.3f} arcmin"
    )
