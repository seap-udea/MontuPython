import montu
from montu_gui.modules.date_converter import (
    SOTHIC_YEAR_BCE,
    SOTHIC_YEAR_CE,
    SOTHIC_YEAR_HORUS,
    SOTHIC_YEAR_MIXED,
    _format_sothic,
    sothic_display_year,
    sothic_era_year_mode,
    sothic_horus_from_year,
    sothic_to_julian_gregorian,
)


def test_sothic_format_mixed_bce():
    assert _format_sothic(2026, "I", "akhet", 2, year_mode=SOTHIC_YEAR_BCE) == (
        "[bce 2026] I akhet 2"
    )
    assert _format_sothic(-2025, "I", "akhet", 2, year_mode=SOTHIC_YEAR_MIXED) == (
        "[-2025] I akhet 2"
    )


def test_sothic_horus_from_mixed_year_modes():
    assert sothic_horus_from_year(756, SOTHIC_YEAR_HORUS) == 756
    assert sothic_horus_from_year(2026, SOTHIC_YEAR_BCE) == 756
    assert sothic_horus_from_year(-2025, SOTHIC_YEAR_MIXED) == 756


def test_sothic_display_year_modes():
    assert sothic_display_year(756, SOTHIC_YEAR_HORUS) == 756
    assert sothic_display_year(756, SOTHIC_YEAR_BCE) == 2026
    assert sothic_display_year(756, SOTHIC_YEAR_MIXED) == -2025


def test_sothic_era_year_mode():
    assert sothic_era_year_mode(0) == SOTHIC_YEAR_BCE
    assert sothic_era_year_mode(756) == SOTHIC_YEAR_BCE
    assert sothic_era_year_mode(4810) == SOTHIC_YEAR_CE


def test_sothic_to_julian_gregorian_mixed_year_modes():
    ref = sothic_to_julian_gregorian(
        756, "I", "akhet", 2, year_mode=SOTHIC_YEAR_HORUS
    )
    for mode, year in (
        (SOTHIC_YEAR_BCE, 2026),
        (SOTHIC_YEAR_MIXED, -2025),
    ):
        result = sothic_to_julian_gregorian(
            year, "I", "akhet", 2, year_mode=mode
        )
        assert result.ok
        assert result.jd_utc == ref.jd_utc
        assert result.can_hyear == 756

    modern = sothic_to_julian_gregorian(
        2026, "I", "akhet", 2, year_mode=SOTHIC_YEAR_CE
    )
    assert modern.ok
    assert modern.can_hyear == 4810
