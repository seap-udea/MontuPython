import pytest

from montu_gui.modules.date_converter import parse_sothic_launch_arg


def test_parse_full_sothic_date_highlights_day():
    req = parse_sothic_launch_arg("[hrw 0] I akhet 1")
    assert req.horus_year == 0
    assert req.month == "I"
    assert req.season == "akhet"
    assert req.day == 1
    assert req.highlight_day is True


def test_parse_full_sothic_date_without_brackets():
    req = parse_sothic_launch_arg("hrw 0 I akhet 1")
    assert req.horus_year == 0
    assert req.highlight_day is True


def test_parse_horus_year_only():
    req = parse_sothic_launch_arg("hrw 0")
    assert req.horus_year == 0
    assert req.highlight_day is False


def test_parse_horus_year_only_accepts_hwr_typo():
    req = parse_sothic_launch_arg("[hwr 0]")
    assert req.horus_year == 0
    assert req.highlight_day is False


def test_parse_bce_year_only():
    req = parse_sothic_launch_arg("bce 1341")
    assert req.horus_year > 0
    assert req.highlight_day is False


def test_parse_astronomical_year_only():
    req = parse_sothic_launch_arg("-1341")
    assert req.horus_year > 0
    assert req.highlight_day is False


def test_parse_empty_raises():
    with pytest.raises(ValueError, match="empty"):
        parse_sothic_launch_arg("   ")
