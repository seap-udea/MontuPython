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


@pytest.mark.parametrize(
    (
        'model', 'kwargs', 'expected_parameters',
        'expected_criterion', 'expected_quantities',
    ),
    [
        (
            'ptolemy',
            {'arcus_visionis_crit': 15.0, 'ptolemy_refraction_deg': 34.0 / 60.0},
            ['arcus_visionis_crit=15.0', f'ptolemy_refraction_deg={34.0 / 60.0}'],
            'criterion: AV_calc >= AV_crit and h_sun < 0°',
            'AV: solar depression at object rise',
        ),
        (
            'schaefer1987',
            {'k': 0.25, 'limiting_mag_zenith': 6.0, 'sun_depression': -10.0},
            ['k=0.25', 'limiting_mag_zenith=6.0', 'sun_depression=-10.0'],
            'criterion: h_star > 0° and V_observed <= V_limit(local)',
            'V_observed: extinguished object magnitude',
        ),
        (
            'schaefer1985',
            {'k': 0.25, 'limiting_mag_zenith': 6.0, 'step_minutes': 2.0},
            ['k=0.25', 'limiting_mag_zenith=6.0', 'step_minutes=2.0'],
            'criterion: h_star > 0° and V <= V_limit',
            'V: object magnitude; V_limit: limiting magnitude at object altitude',
        ),
        (
            'belokrylov2011',
            {'k': 0.25, 'reference_extinction': 0.25, 'step_minutes': 2.0},
            ['k=0.25', 'reference_extinction=0.25', 'step_minutes=2.0'],
            'criterion: h_star > 0° and h_sun <= h_theor',
            "m': extinction-corrected magnitude",
        ),
    ],
)
def test_compute_verbose_reports_model_parameters(
    thebes, sirius, capsys, model, kwargs, expected_parameters,
    expected_criterion, expected_quantities,
):
    day = montu.Time('139-07-20 00:00:00', calendar='mixed')
    calculator = montu.HeliacalRise(model=model, **kwargs)

    calculator.compute(sirius, thebes, day, day, verbose=True)

    captured = capsys.readouterr().out
    assert f'HeliacalRise verbose — model={model}' in captured
    assert 'day 001' in captured
    assert 'visible=' in captured
    for parameter in expected_parameters:
        assert parameter in captured
    assert expected_criterion in captured
    assert expected_quantities in captured
    assert 'criterion satisfied:' in captured
    assert 'body_azimuth_deg' not in captured
    assert 'body_ra_hours' not in captured
    assert 'body_dec_deg' not in captured
    if model == 'ptolemy':
        assert 't_rise: local object rise time' in captured
        assert 't_rise=' in captured
    if model in ('schaefer1985', 'belokrylov2011'):
        assert 'morning scan: h_sun=-18.0° to sunrise, every 2.0 minutes' in captured
    else:
        assert 'morning scan:' not in captured


@pytest.mark.parametrize(
    ('refraction_deg', 'expected_horizon_deg'),
    [(0.0, 0.0), (0.5, -0.5)],
)
def test_ptolemy_refraction_sets_horizon_directly(
    thebes, sirius, refraction_deg, expected_horizon_deg,
):
    start = montu.Time('139-07-01 00:00:00', calendar='mixed')
    end = montu.Time(start.jed + 44, format='jd', scale='utc')
    result = montu.HeliacalRise(
        model='ptolemy',
        ptolemy_refraction_deg=refraction_deg,
    ).compute(sirius, thebes, start, end)

    assert result.target_horizon_deg.iloc[0] == pytest.approx(expected_horizon_deg)


def test_heliacal_rise_wrapper_accepts_verbose(thebes, sirius, capsys):
    day = montu.Time('139-07-20 00:00:00', calendar='mixed')

    montu.heliacal_rise(
        sirius, thebes, day, day, model='schaefer1987', verbose=True,
    )

    captured = capsys.readouterr().out
    assert 'HeliacalRise verbose — model=schaefer1987' in captured
    assert 'day 001' in captured


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
