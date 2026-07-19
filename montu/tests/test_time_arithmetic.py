import math

import pytest

import montu


def test_tt_add_and_subtract_with_constants():
    t = montu.Time("2000-01-01 00:00:00")
    shifted = t + 100 * montu.YEAR + 120 * montu.DAY
    back = shifted - 100 * montu.YEAR - 120 * montu.DAY
    assert math.isclose(back.tt, t.tt, rel_tol=0, abs_tol=1e-6)


def test_calendar_add_and_subs():
    t = montu.Time("2000-01-01 00:00:00")
    forward = t.add(years=1, days=10)
    back = forward.subs(years=1, days=10)
    assert forward.readable.datepro.startswith("2001-01-11")
    assert back.readable.datepro.startswith("2000-01-01")


def test_time_minus_time_returns_tt_seconds():
    t1 = montu.Time("2000-01-01 00:00:00")
    t2 = montu.Time("2001-01-01 00:00:00")
    diff_tt = t2 - t1
    assert isinstance(diff_tt, float)
    assert math.isclose(diff_tt, t2.tt - t1.tt, rel_tol=0, abs_tol=1e-6)


def test_diff_returns_calendar_delta():
    t1 = montu.Time("2000-01-01 00:00:00")
    t2 = montu.Time("2001-01-01 00:00:00")
    diff_cal = t2.diff(t1)
    assert isinstance(diff_cal, montu.CalendarDelta)
    years, days, hours, minutes, seconds = diff_cal
    assert years == 1
    assert days == hours == minutes == seconds == 0
    assert math.isclose(diff_cal.to_days(), 366.0, rel_tol=0, abs_tol=1e-6)


def test_sothic_apokatastasis_calendar_diff():
    t_apo2 = montu.Time("[bce 1322] I akhet 1", calendar="sothic")
    t_apo3 = montu.Time("[139] I akhet 1", calendar="sothic")
    diff_cal = t_apo3.diff(t_apo2)
    assert diff_cal.years == 1461
    assert math.isclose(diff_cal.to_years(), 1461.0, rel_tol=0, abs_tol=1e-6)


def test_mixed_gregorian_reform_diff():
    t1 = montu.Time("1582-10-04", calendar="mixed")
    t2 = montu.Time("1582-10-15", calendar="mixed")
    assert math.isclose(t2.diff(t1).to_days(), 1.0, rel_tol=0, abs_tol=1e-6)


def test_mixed_gregorian_reform_add():
    t1 = montu.Time("1582-10-04", calendar="mixed")
    t2 = t1.add(days=1)
    assert t2.readable.datemix.startswith("1582-10-15")
    t3 = t2.subs(days=1)
    assert t3.readable.datemix.startswith("1582-10-04")


def test_mixed_and_proleptic_agree_after_reform():
    t_pro = montu.Time("2000-01-01 00:00:00", calendar="proleptic")
    t_mix = montu.Time("2000-01-01 00:00:00", calendar="mixed")
    pro = t_pro.add(years=100, days=120)
    mix = t_mix.add(years=100, days=120)
    assert pro.readable.datemix == mix.readable.datemix


def test_mixed_calendar_add_and_subs():
    t = montu.Time("2000-01-01 00:00:00", calendar="mixed")
    forward = t.add(years=100, days=120)
    assert forward.readable.datemix.startswith("2100-")
    back = forward.subs(years=100, days=120)
    assert abs(back.jtd - t.jtd) <= 1.0


def test_sub_deprecated_seconds_still_works():
    t = montu.Time("2000-01-01 00:00:00")
    back = t.sub(montu.DAY)
    assert math.isclose(back.tt, t.tt - montu.DAY, rel_tol=0, abs_tol=1e-6)
