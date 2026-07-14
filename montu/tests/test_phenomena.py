import pytest

import montu


@pytest.fixture(scope="module")
def thebes():
    return montu.Observer(lon=33, lat=24, height=0.075)


@pytest.fixture(scope="module")
def sirius():
    return montu.Stars(subset='bright', ProperName='Sirius')


def test_heliacal_rise_class_models(thebes, sirius):
    start = montu.Time('139-07-01 00:00:00', calendar='mixed')
    end = montu.Time(start.jed + 44, format='jd', scale='utc')

    for model in montu.HELIACAL_RISE_MODELS:
        calculator = montu.HeliacalRise(model=model)
        result = calculator.compute(sirius, thebes, start, end)
        assert calculator.source in montu.HELIACAL_RISE_SOURCES[model]
        assert 'jed' in result.columns
        assert result.model.iloc[0] == model


def test_schaefer1987_sirius_third_apokatastasis(thebes, sirius):
    start = montu.Time('139-07-01 00:00:00', calendar='mixed')
    end = montu.Time(start.jed + 44, format='jd', scale='utc')

    result = montu.HeliacalRise(
        model='schaefer1987', k=0.25, sun_depression=-11.0,
    ).compute(sirius, thebes, start, end)
    assert len(result) == 1
    assert int(result.day_jed.iloc[0] - start.jed) == 19


def test_heliacal_rise_function_wrapper(thebes, sirius, capsys):
    start = montu.Time('139-07-01 00:00:00', calendar='mixed')
    end = montu.Time(start.jed + 44, format='jd', scale='utc')
    result = montu.heliacal_rise(sirius, thebes, start, end, model='schaefer1987')
    assert not result.empty
    captured = capsys.readouterr().out
    assert 'schaefer1987' in captured
    assert 'source:' in captured


def test_print_rises_returns_dataframe(thebes, sirius):
    start = montu.Time('139-07-01 00:00:00', calendar='mixed')
    end = montu.Time(start.jed + 44, format='jd', scale='utc')
    calculator = montu.HeliacalRise(model='schaefer1987')
    result = calculator.compute(sirius, thebes, start, end)
    printed = calculator.print_rises(result, title='Test')
    assert printed is result
    assert not printed.empty


def test_sebau_body_uses_ephemeris_magnitude(thebes):
    venus = montu.Planet('Venus')
    start = montu.Time('2020-01-01 00:00:00')
    end = montu.Time(start.jed + 30, format='jd', scale='utc')

    result = montu.HeliacalRise(model='schaefer1987', sun_depression=-6.0).compute(
        venus, thebes, start, end,
    )
    if not result.empty:
        assert result.vmag.notna().all()
