import math

import pytest

import montu


def test_catalogue_site_pressure_temperature_override():
    thebes_geo = montu.Observer(site="thebes", pressure=0, temperature=0)
    assert thebes_geo.pressure == 0
    assert thebes_geo.temperature == 0
    assert thebes_geo.site.pressure == 0
    assert thebes_geo.site.temp == 0
    assert thebes_geo.site_id == "thebes"
    entry = next(
        item for item in montu.Observer.list(details=True) if item["id"] == "thebes"
    )
    assert thebes_geo.lat == pytest.approx(entry["lat"])
    assert thebes_geo.lon == pytest.approx(entry["lon"])


def test_catalogue_site_coordinate_override():
    entry = next(
        item for item in montu.Observer.list(details=True) if item["id"] == "thebes"
    )
    thebes = montu.Observer(site="thebes", lon=30.0, height=0.1)
    assert thebes.lon == pytest.approx(30.0)
    assert thebes.lat == pytest.approx(entry["lat"])
    assert thebes.height == pytest.approx(0.1)


def test_catalogue_site_exposes_all_location_metadata():
    thebes = montu.Observer(site="thebes")
    entry = next(
        item for item in montu.Observer.list(details=True) if item["id"] == "thebes"
    )
    for key, value in entry.items():
        assert getattr(thebes, key) == value
    assert thebes.site_id == "thebes"
    assert thebes.site_name == entry["name"]
    assert thebes.height == pytest.approx(entry["alt_m"] / 1000.0)
    assert thebes.pressure == pytest.approx(entry["pressure_mbar"])
    assert thebes.temperature == pytest.approx(entry["temperature_c"])


def test_manual_observer_has_no_catalogue_metadata():
    obs = montu.Observer(lon=0, lat=0)
    assert obs.site_id is None
    assert obs.site_name is None
    assert not hasattr(obs, "region")


def test_distance_to_between_catalogue_sites():
    alexandria = montu.Observer(site="alexandria")
    aswan = montu.Observer(site="aswan")
    distance_km = aswan.distance_to(alexandria)
    assert math.isclose(distance_km, 843.37, rel_tol=0, abs_tol=0.5)


def test_distance_to_matches_haversine():
    a = montu.Observer(lon=31.2001, lat=31.2001, height=7 / 1000)
    b = montu.Observer(site="aswan")
    arc = montu.Util.haversine_distance(
        a.lat * montu.DEG, a.lon * montu.DEG,
        b.lat * montu.DEG, b.lon * montu.DEG,
    )
    radius = float(montu.load_planets().loc["Earth", "Rmean"])
    assert math.isclose(a.distance_to(b), arc * radius, rel_tol=0, abs_tol=1e-6)
    assert math.isclose(a.distance_to(b, units="deg"), arc * montu.RAD, rel_tol=0, abs_tol=1e-12)


def test_sidereal_time_at_thebes():
    thebes = montu.Observer(site="thebes")
    mtime = montu.Time("bce 1500-06-21 12:00:00", calendar="mixed")
    assert thebes.sidereal_time(mtime) == "07:14:54.925"
    assert thebes.sidereal_time(mtime, hms=False) == pytest.approx(
        7.248590315909614, abs=1e-9
    )


def test_observer_repr_for_catalogue_and_manual_sites():
    entry = next(
        item for item in montu.Observer.list(details=True) if item["id"] == "thebes"
    )
    thebes = montu.Observer(site="thebes")
    assert repr(thebes) == (
        f"Observer('thebes'/'{entry['name']}'/"
        f"{entry['lat']:.6f}°, {entry['lon']:.6f}°, {entry['alt_m']:.0f} m/"
        f"P={entry['pressure_mbar']} mbar, T={entry['temperature_c']} °C)"
    )

    manual = montu.Observer(lon=31.134, lat=29.979, height=0.075)
    assert repr(manual) == (
        "Observer('29.979000°, 31.134000°, 75 m/P=1013.25 mbar, T=15.0 °C)"
    )


def test_observer_str_includes_catalogue_metadata():
    thebes = montu.Observer(site="thebes")
    text = str(thebes)
    assert text.startswith("Observer\n")
    assert "Site: Thebes (Luxor) [thebes]" in text
    assert "Region: Egypt · Ancient Egypt" in text
    assert "lat 25.696700°, lon 32.642200°, elevation 76 m" in text
    assert "Ancient Waset" in text

    manual = montu.Observer(lon=31.134, lat=29.979, height=0.075)
    manual_text = str(manual)
    assert "Site:" not in manual_text
    assert "Description:" not in manual_text
    assert "lat 29.979000°, lon 31.134000°, elevation 75 m" in manual_text
