import montu
import pytest


@pytest.fixture(scope='module')
def mars_aldebaran():
    mars = montu.Planet('Mars')
    aldebaran = montu.Stars(
        subset='bright', ProperName='Aldebaran', return_as='Star',
    )
    return mars, aldebaran


def test_print_dict_scalar_mapping(capsys):
    montu.Util.print_dict({'visible': True, 'separation': 4.275, 'note': None})
    out = capsys.readouterr().out
    assert 'visible' in out
    assert 'yes' in out
    assert '4.28' in out
    assert '—' in out


def test_print_dict_dictobj(capsys):
    data = montu.Dictobj(dict={'in_conjunction': False, 'maxseparation': 5})
    montu.Util.print_dict(data, title='Visibility')
    out = capsys.readouterr().out
    assert out.startswith('Visibility')
    assert 'in_conjunction' in out
    assert 'no' in out


def test_print_dict_expands_list_of_dicts(capsys):
    montu.Util.print_dict({
        'visible': True,
        'body_conditions': [
            {'name': 'Mars', 'el': 72.26, 'above_horizon': True},
            {'name': 'Aldebaran', 'el': 76.05, 'above_horizon': True},
        ],
    })
    out = capsys.readouterr().out
    assert 'body_conditions' in out
    assert 'see below' in out
    assert 'Mars' in out
    assert 'Aldebaran' in out


def test_print_dict_conjunction_is_visible(mars_aldebaran, capsys):
    mars, aldebaran = mars_aldebaran
    site = montu.Observer(lat=6, lon=-75)
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07 09:00:00'),
        observer=site,
    )
    conditions = conj.is_visible(from_site=site, verbose=False)
    montu.Util.print_dict(conditions, title='Conjunction visibility')
    out = capsys.readouterr().out
    assert 'Conjunction visibility' in out
    assert 'visible_from_site' in out
    assert 'body_conditions' in out
    assert 'Mars' in out


def test_print_dict_nested_dict(capsys):
    montu.Util.print_dict({'meta': {'a': 1, 'b': 2}}, expand_tables=False)
    out = capsys.readouterr().out
    assert 'meta' in out
    assert 'a' in out
