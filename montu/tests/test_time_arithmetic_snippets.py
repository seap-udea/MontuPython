"""Time arithmetic examples from ``MontuPython-CodeSnippets.ipynb``.

These cases mirror the notebook cells on ``diff``, ``add``, ``subs``, and TT
``+`` / ``-``. Expected values were cross-checked with
`timeanddate.com <https://www.timeanddate.com/date/duration.html>`_ date
difference and date calculators (Gregorian reform, Columbus quincentenary,
century shifts, etc.).
"""

import math

import pytest

import montu


# --- diff(): reform, civil span, Sothic cycle --------------------------------


def test_gregorian_reform_mixed_calendar_diff():
    """1582-10-04 → 1582-10-15 is one civil day in the mixed calendar."""
    t1 = montu.Time("1582-10-04", calendar="mixed")
    t2 = montu.Time("1582-10-15", calendar="mixed")
    assert math.isclose(t2.diff(t1).to_days(), 1.0, rel_tol=0, abs_tol=1e-6)


def test_gregorian_reform_proleptic_calendar_diff():
    """Same civil labels span eleven days in the proleptic Gregorian calendar."""
    t1 = montu.Time("1582-10-04", calendar="proleptic")
    t2 = montu.Time("1582-10-15", calendar="proleptic")
    assert math.isclose(t2.diff(t1).to_days(), 11.0, rel_tol=0, abs_tol=1e-6)


def test_proleptic_vs_mixed_same_civil_date():
    """1582-10-04 proleptic is ten Julian days after the mixed-calendar label."""
    t_mixed = montu.Time("1582-10-04", calendar="mixed")
    t_proleptic = montu.Time("1582-10-04", calendar="proleptic")
    assert t_proleptic.diff(t_mixed).to_days() == pytest.approx(-10.0, abs=1e-5)


def test_columbus_quincentenary_calendar_diff():
    """1492-10-12 to 1992-10-12 (mixed): timeanddate.com component breakdown."""
    t1 = montu.Time("1492-10-12", calendar="mixed")
    t2 = montu.Time("1992-10-12", calendar="mixed")
    delta = t2.diff(t1)

    assert delta.years == 499
    assert delta.days == 29
    assert delta.hours == 23
    assert delta.minutes == 57
    assert delta.seconds == 32
    assert math.isclose(delta.to_days(), 182612.0, rel_tol=0, abs_tol=1e-6)
    assert math.isclose(delta.to_years(), 500.306849, rel_tol=0, abs_tol=1e-6)


def test_sothic_apokatastasis_calendar_diff_and_subs():
    """Second to third Sothic apokatastasis spans 1461 calendar years."""
    t_apo3 = montu.Time("[139] I akhet 1 00:00:00", calendar="sothic")
    t_apo2 = montu.Time("[bce 1322] I akhet 1", calendar="sothic")

    delta_apo = t_apo3.diff(t_apo2)
    assert delta_apo.years == 1461
    assert math.isclose(delta_apo.to_years(), 1461.0, rel_tol=0, abs_tol=1e-6)

    t_apo1 = t_apo2.subs(years=delta_apo.to_years())
    assert t_apo1.readable.datemix.startswith("-2781-07-20")
    assert t_apo1.readable.datesot == "[hrw 0] I akhet 1"


# --- add(): mixed reform step and calendar sum --------------------------------


def test_mixed_reform_add_one_day():
    """1582-10-04 + 1 civil day skips the historical gap to 1582-10-15."""
    t1 = montu.Time("1582-10-04", calendar="mixed")
    tcal = t1.add(days=1)
    assert tcal.readable.datemix.startswith("1582-10-15 00:00:00")
    assert "1582-10-15" in repr(tcal)


def test_calendar_add_century_and_days():
    """2000-01-01 + 100 years + 120 days → 2100-05-01 (calendar units)."""
    t = montu.Time("2000-01-01 00:00:00")
    tcal = t.add(years=100, days=120)
    assert tcal.readable.datepro.startswith("2100-05-01 00:00:00")


# --- TT ``+`` / ``-``: ephemeris seconds --------------------------------------


def test_tt_add_century_and_days():
    """2000-01-01 + 100×YEAR + 120×DAY (TT seconds, not civil months)."""
    t = montu.Time("2000-01-01 00:00:00")
    tcal = t + 100 * montu.YEAR + 120 * montu.DAY
    assert tcal.readable.datepro.startswith("2100-04-30 23:57:40")


def test_tt_subtract_century_and_days():
    """2000-01-01 − 100×YEAR − 120×DAY (TT seconds)."""
    t = montu.Time("2000-01-01 00:00:00")
    tcal = t - 100 * montu.YEAR - 120 * montu.DAY
    assert tcal.readable.datepro.startswith("1899-09-02 00:01:07")


def test_tt_add_subtract_round_trip():
    t = montu.Time("2000-01-01 00:00:00")
    back = (t + 100 * montu.YEAR + 120 * montu.DAY) - 100 * montu.YEAR - 120 * montu.DAY
    assert back.readable.datepro.startswith("2000-01-01 00:00:00")


# --- subs(): calendar units ---------------------------------------------------


def test_calendar_subs_century_and_days():
    """2000-01-01 − 100 years − 120 days (calendar units; differs from TT ``-``)."""
    t = montu.Time("2000-01-01 00:00:00")
    tcal = t.subs(years=100, days=120)
    assert tcal.readable.datepro.startswith("1899-09-03 00:00:00")
