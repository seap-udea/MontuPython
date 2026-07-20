import pytest

import montu


def test_show_position_without_where_in_sky(capsys):
    spica = montu.Star(
        montu.Stars(subset="bright", ProperName="Spica").data.iloc[0]
    )
    assert spica.show_position() is None
    out = capsys.readouterr().out
    assert "no sky position stored" in out
    assert "where_in_sky" in out


def test_show_conditions_without_conditions_in_sky(capsys):
    mars = montu.Planet("Mars")
    assert mars.show_conditions() is None
    out = capsys.readouterr().out
    assert "no sky conditions stored" in out
    assert "conditions_in_sky" in out


def test_show_position_for_star(thebes_observer, capsys):
    spica = montu.Star(
        montu.Stars(subset="bright", ProperName="Spica").data.iloc[0]
    )
    mtime = montu.Time("bce 1500-06-21 12:00:00", calendar="mixed")
    spica.where_in_sky(at=mtime, observer=thebes_observer)
    assert hasattr(spica.position, "RAJ2000t")
    assert hasattr(spica.position, "DecJ2000t")
    expected_ra, expected_dec = montu.Star._propagated_j2000_coordinates(
        spica._star_data, mtime
    )
    assert spica.position.RAJ2000t == pytest.approx(expected_ra, abs=1e-9)
    assert spica.position.DecJ2000t == pytest.approx(expected_dec, abs=1e-9)
    assert spica.show_position() is None
    out = capsys.readouterr().out
    assert "Spica — sky position" in out
    assert "Epoch:" in out
    assert "Site: Thebes (Luxor) [thebes]" in out
    assert "RA (J2000, proper motion):" in out
    assert "Dec (J2000, proper motion):" in out
    assert "RA (J2000):" in out
    assert "Azimuth:" in out
    assert "Elevation:" in out


def test_show_conditions_star_omits_planet_fields(thebes_observer, capsys):
    spica = montu.Star(
        montu.Stars(subset="bright", ProperName="Spica").data.iloc[0]
    )
    mtime = montu.Time("bce 1500-06-21 12:00:00", calendar="mixed")
    spica.conditions_in_sky(at=mtime, observer=thebes_observer)
    assert spica.show_conditions() is None
    out = capsys.readouterr().out
    assert "Rise time (UTC):" in out
    assert "Distance from Earth" not in out
    assert "Illuminated fraction" not in out


def test_show_conditions_planet_includes_phase(giza_observer, capsys):
    mars = montu.Planet("Mars")
    mtime = montu.Time("-1000-03-21 20:00:00")
    mars.conditions_in_sky(at=mtime, observer=giza_observer)
    assert mars.show_conditions() is None
    out = capsys.readouterr().out
    assert "Distance from Earth:" in out
    assert "Illuminated fraction:" in out
    assert "Heliocentric latitude:" in out
