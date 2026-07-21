import montu
from montu.stars import (
    clear_stellar_catalogue_cache,
    load_stellar_catalogue,
    stellar_catalogue_cache_status,
)


def setup_function():
    clear_stellar_catalogue_cache()


def test_cache_reuses_same_tier_without_second_disk_read(capsys):
    montu.Stars(subset='bright')
    out_first = capsys.readouterr().out
    assert out_first.count("Loading stellar catalogue") == 1

    montu.Stars(subset='bright')
    out_second = capsys.readouterr().out
    assert "Loading stellar catalogue" not in out_second
    assert stellar_catalogue_cache_status() == {'bright': 291}


def test_visible_after_bright_reads_visible_file(capsys):
    montu.Stars(subset='bright')
    capsys.readouterr()

    visible = montu.Stars(subset='visible')
    out = capsys.readouterr().out
    assert "montu_stellar_catalogue_v38_visible.csv" in out
    assert visible.number > 5000
    assert set(stellar_catalogue_cache_status()) == {'bright', 'visible'}


def test_bright_after_visible_is_derived_not_reread(capsys):
    montu.Stars(subset='visible')
    capsys.readouterr()

    bright = montu.Stars(subset='bright')
    out = capsys.readouterr().out
    assert "Loading stellar catalogue" not in out
    assert bright.number < 500
    assert set(stellar_catalogue_cache_status()) == {'visible', 'bright'}


def test_visible_after_full_is_derived(capsys):
    montu.Stars()
    capsys.readouterr()

    visible = montu.Stars(subset='visible')
    out = capsys.readouterr().out
    assert "Loading stellar catalogue" not in out
    assert visible.number > 5000
    assert 'full' in stellar_catalogue_cache_status()
    assert 'visible' in stellar_catalogue_cache_status()


def test_derived_bright_matches_packaged_subset():
    direct = load_stellar_catalogue('bright')
    clear_stellar_catalogue_cache()
    load_stellar_catalogue('visible')
    derived = load_stellar_catalogue('bright', verbose=False)
    assert len(derived) == len(direct)
    assert set(direct['HIP']) == set(derived['HIP'])
