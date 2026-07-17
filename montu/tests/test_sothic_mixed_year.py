import pytest
import montu


def test_sothic_mixed_bce_year():
    t = montu.Time('[bce 2026] I akhet 2', calendar='sothic')
    ref = montu.Time('[hrw 756] I akhet 2', calendar='sothic')
    assert t.jed == ref.jed
    assert t.hed == ref.hed


def test_sothic_mixed_astronomical_year():
    t = montu.Time('[-2025] I akhet 2', calendar='sothic')
    ref = montu.Time('[hrw 756] I akhet 2', calendar='sothic')
    assert t.jed == ref.jed


def test_sothic_mixed_bce_and_astronomical_equivalent():
    t_bce = montu.Time('[bce 2026] I akhet 2', calendar='sothic')
    t_astro = montu.Time('[-2025] I akhet 2', calendar='sothic')
    assert t_bce.jed == t_astro.jed


def test_sothic_mixed_apokatastasis():
    t = montu.Time('[bce 2782] I akhet 1', calendar='sothic')
    ref = montu.Time('[hrw 0] I akhet 1', calendar='sothic')
    assert t.jed == ref.jed


def test_sothic_mixed_ce_year():
    t = montu.Time('[2026] I akhet 1', calendar='sothic')
    ref = montu.Time('[hrw 4810] I akhet 1', calendar='sothic')
    assert t.jed == ref.jed


def test_sothic_mixed_ce_tag():
    t = montu.Time('[ce 2026] I akhet 1', calendar='sothic')
    ref = montu.Time('[2026] I akhet 1', calendar='sothic')
    assert t.jed == ref.jed
