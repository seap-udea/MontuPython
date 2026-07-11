import math

import numpy as np
import pandas as pd
import pytest

import montu


pytestmark = pytest.mark.docstrings


def test_util_examples_match_docstrings():
    assert montu.Util.dec2sex(15.5) == "15:30:00.000"
    assert montu.Util.dec2sex(15.5, string=False) == (15.0, 30, 0.0)
    assert montu.Util.sex2dec("15:30:00.000") == 15.5
    assert montu.Util.sex2dec((-7, 15, 0.0)) == -7.25
    assert np.allclose(montu.Util.arange(0, 1, 0.25), [0.0, 0.25, 0.5, 0.75, 1.0])


def test_time_examples_cover_supported_input_formats():
    proleptic = montu.Time("-2500-01-01 12:00:00.00", format="iso", scale="utc", calendar="proleptic")
    mixed = montu.Time("-2500-01-22 12:00:00.00", format="iso", scale="utc", calendar="mixed")
    julian_day = montu.Time(807954, format="jd", scale="utc").get_readable()
    terrestrial = montu.Time(0, format="tt", scale="utc").get_readable()

    assert proleptic.readable.datepro.startswith("-2500-01-01")
    assert mixed.readable.datemix.startswith("-2500-01-22")
    assert julian_day.jed == 807954
    assert isinstance(terrestrial.tt, float)
    assert terrestrial.readable.datepro

    shifted = terrestrial.add(365 * montu.DAY)
    assert math.isclose(shifted.diff(terrestrial), 365.0, rel_tol=0, abs_tol=1e-6)


def test_stars_examples_filter_and_project_to_sky(andes_observer):
    allstars = montu.Stars()
    bright = allstars.get_stars(Vmag=[-2, 4])
    assert bright.number > 0

    mtime = montu.Time("2024-05-01 19:00:00")
    sky = bright.where_in_sky(at=mtime, observer=andes_observer)

    assert isinstance(sky, pd.DataFrame)
    assert {"az", "el", "HA"}.issubset(sky.columns)
    assert np.isfinite(sky["az"]).all()
    assert np.isfinite(sky["el"]).all()


def test_sun_and_planet_examples_compute_conditions(giza_observer):
    mtime = montu.Time("-1000-03-21 18:00:00")

    sun = montu.Sun()
    sun.where_in_sky(at=mtime, observer=giza_observer)
    assert 0 <= sun.position.az <= 360
    assert -90 <= sun.position.el <= 90

    mars = montu.Planet("Mars")
    mars.conditions_in_sky(at=mtime, observer=giza_observer)
    assert mars.condition.rise_time > 0
    assert np.isfinite(mars.condition.elongation)

    seasons = montu.Sun.next_seasons(at=montu.Time("-1000-01-01"))
    assert len(seasons) == 4
    assert list(seasons) == sorted(seasons)
