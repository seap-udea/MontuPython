import math

import montu


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
