import montu


def test_load_planets_exposes_synodic_orbit():
    planets = montu.load_planets()
    assert "SynodicOrbit" in planets.columns
    assert planets.loc["Mars", "SiderealOrbit"] > 0


def test_astro_where_in_sky(andes_observer):
    mtime = montu.Time("2024-05-01 19:00:00")
    az, el = montu.Astro.where_in_sky(
        RA=6.770358,
        Dec=-16.751203,
        at=mtime,
        observer=andes_observer,
    )
    assert 0 <= az <= 360
    assert -90 <= el <= 90
