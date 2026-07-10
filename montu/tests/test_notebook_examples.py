import numpy as np
import pandas as pd
import pytest

import montu


pytestmark = pytest.mark.notebooks


def test_basic_functions_notebook_smoke(andes_observer):
    visible = montu.Stars().get_stars(Vmag=[-2, 4])
    mtime = montu.Time("2024-05-01 19:00:00")

    visible.where_in_sky(at=mtime, observer=andes_observer, inplace=True)
    assert {"az", "el", "HA"}.issubset(visible.data.columns)

    az, el = montu.Util.where_in_sky(RA=6.770358, Dec=-16.751203, at=mtime, observer=andes_observer)
    assert 0 <= az <= 360
    assert -90 <= el <= 90

    allstars = montu.Stars()
    bright = allstars.get_stars(Vmag=[-2, 6.5])
    assert len(bright.data) >= len(visible.data)


def test_montime_notebook_smoke():
    proleptic = montu.Time("-2500-01-01 12:00:00.00", format="iso", scale="utc", calendar="proleptic")
    mixed = montu.Time("-2500-01-22 12:00:00.00", format="iso", scale="utc", calendar="mixed")
    from_jd = montu.Time(807954, format="jd", scale="utc").get_readable()
    from_tt = montu.Time(0, format="tt", scale="utc").get_readable()

    assert proleptic.readable.datepro.startswith("-2500-01-01")
    assert mixed.readable.datemix.startswith("-2500-01-22")
    assert from_jd.jed == 807954
    assert isinstance(from_tt.readable.datepro, str)
    assert "bce" not in from_tt.readable.datepro.lower()


def test_planesticios_notebook_smoke(egypt_observer):
    mars = montu.Planet("Mars")
    now = montu.Time()

    mars.where_in_sky(at=now, observer=egypt_observer, store=True)
    mars.conditions_in_sky(at=now, observer=egypt_observer, store=True)
    mars.tabulate_ephemerides()

    assert isinstance(mars.ephemerides, pd.DataFrame)
    assert not mars.ephemerides.empty
    assert {"az", "el", "Vmag"}.issubset(mars.ephemerides.columns)


def test_titovivas_notebook_smoke():
    allstars = montu.Stars()
    aldebaran = allstars.get_stars(ProperName="Aldebaran")
    mtime = montu.Time("BCE1472-01-01 00:00:00")
    senenmut = montu.Observer(lon=32, lat=26, height=0.089)

    aldebaran.where_in_sky(at=mtime, observer=senenmut, inplace=True)
    az, el = montu.Util.where_in_sky(
        aldebaran.data.RAEpoch.iloc[0],
        aldebaran.data.DecEpoch.iloc[0],
        at=mtime,
        observer=senenmut,
    )

    assert np.isfinite(az)
    assert np.isfinite(el)
    assert -90 <= el <= 90


def test_venus_azimuths_notebook_smoke(egypt_observer):
    mtime = montu.Time("-300-12-04 12:00:00")
    venus = montu.Planet("Venus")
    sun = montu.Sun()

    venus.conditions_in_sky(at=mtime, observer=egypt_observer)
    sun.conditions_in_sky(at=mtime, observer=egypt_observer)

    assert np.isfinite(venus.condition.elongation)
    assert np.isfinite(venus.position.az)
    assert np.isfinite(sun.position.az)


def test_montunctions_notebook_smoke(egypt_observer):
    mars = montu.Planet("Mars")
    mtime = montu.Time("-700-01-01 00:00:00.00")
    allstars = montu.Stars()
    aldebaran = allstars.get_stars(ProperName="Aldebaran")

    mars.conditions_in_sky(at=mtime, observer=egypt_observer)
    aldebaran.where_in_sky(at=mtime, observer=egypt_observer, inplace=True)

    assert np.isfinite(mars.condition.Vmag)
    assert np.isfinite(aldebaran.data.az.iloc[0])


def test_get_stars_around_accepts_series_center():
    allstars = montu.Stars()
    aldebaran = allstars.get_stars(ProperName="Aldebaran")
    hyades = allstars.get_stars_around(
        center=[aldebaran.data.RAJ2000, aldebaran.data.DecJ2000],
        radius=5.5,
        Vmag=[-1, 5],
    )

    assert len(hyades.data) > 0
    assert (hyades.data.Vmag <= 5).all()


def test_where_in_space_value_for_polar_stars():
    star_names = ("Polaris", "Vega", "Thuban", "Deneb", "Alderamin", "Kochab")
    stars = montu.Stars().get_stars(ProperName=star_names)
    past = montu.Time() + (-5000 * montu.YEAR)
    pstars = stars.where_in_space(at=past)

    assert isinstance(pstars, montu.Stars)
    for star in star_names:
        dec = pstars.value_for(star, "DecEpoch")
        assert np.isfinite(dec)
        assert -90 <= dec <= 90


def test_stars_scalar_after_where_in_space():
    aldebaran = montu.Stars().get_stars(ProperName="Aldebaran")
    mtime = montu.Time("-700-01-01 00:00:00")
    aldebaran.where_in_space(at=mtime, inplace=True)
    dec = aldebaran.scalar("DecJ2000t")
    ra = aldebaran.scalar("RAJ2000t")
    assert np.isfinite(dec)
    assert np.isfinite(ra)
