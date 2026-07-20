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


def test_print_rises_is_display_only(thebes, sirius, capsys):
    start = montu.Time('139-07-01 00:00:00', calendar='mixed')
    end = montu.Time(start.jed + 44, format='jd', scale='utc')
    calculator = montu.HeliacalRise(model='schaefer1987')
    result = calculator.compute(sirius, thebes, start, end)
    assert calculator.print_rises(result, title='Test') is None
    assert not result.empty
    captured = capsys.readouterr().out
    assert 'Test —' in captured


def test_sebau_body_uses_ephemeris_magnitude(thebes):
    venus = montu.Planet('Venus')
    start = montu.Time('2020-01-01 00:00:00')
    end = montu.Time(start.jed + 30, format='jd', scale='utc')

    result = montu.HeliacalRise(model='schaefer1987', sun_depression=-6.0).compute(
        venus, thebes, start, end,
    )
    if not result.empty:
        assert result.vmag.notna().all()


@pytest.fixture(scope="module")
def solar_eclipses():
    return montu.SolarEclipses()


def test_solar_eclipses_catalogue_loads(solar_eclipses):
    assert solar_eclipses.number == 11898
    subset = solar_eclipses.get_eclipses(year=[-1400, -1200])
    assert 400 < subset.number < 600
    assert subset.get_eclipse(eclipse_type='T').number > 0


def test_solar_eclipse_dallas_2024_totality(solar_eclipses):
    """Dallas was inside the 2024-04-08 path of totality (~3m51s)."""
    eclipse = solar_eclipses.get_eclipses(year=2024, month=4, day=8).eclipse(0)
    dallas = montu.Observer(lon=-96.7970, lat=32.7767, height=0.14)
    cond = eclipse.conditions_eclipse(dallas)

    assert cond.kind == 'total'
    assert cond.visible is True
    assert cond.magnitude == pytest.approx(1.015, abs=0.01)
    assert cond.duration_umbra_seconds == pytest.approx(231, abs=5)
    assert cond.jed_c1 is not None and cond.jed_c4 is not None
    assert cond.jed_c2 < cond.jed_max < cond.jed_c3


def test_solar_eclipse_thebes_2024_below_horizon(solar_eclipses, thebes):
    """Same eclipse is daytime in Texas but night in Egypt."""
    eclipse = solar_eclipses.get_eclipses(year=2024, month=4, day=8).eclipse(0)
    cond = eclipse.conditions_eclipse(thebes)
    assert cond.kind == 'partial'
    assert cond.visible is False
    assert cond.sun_altitude_deg < 0


def test_solar_eclipse_from_series_and_indexing(solar_eclipses):
    subset = solar_eclipses.get_eclipses(year=2024, month=4, day=8)
    by_index = subset[0]
    by_series = montu.SolarEclipse(subset.data.iloc[0])
    assert repr(by_index).startswith('<SolarEclipse')
    assert by_index.data.year == by_series.data.year == 2024


def test_solar_eclipse_str_is_compact(solar_eclipses):
    eclipse = solar_eclipses.get_eclipses(year=-584, month=5, day=28).eclipse(0)
    text = str(eclipse)
    assert 'Catalogue fields' not in text
    assert 'Calendar year of the eclipse' not in text
    assert 'Eclipse type         : T (total)' in text
    assert 'Greatest eclipse' in text
    assert 'Central path' in text
    assert 'lat_ge, lng_ge       : 38.2N, 45.0W' in text
    assert 'sun_alt, sun_azm     : 71.1°, 158.5°' in text
    assert 'ΔT assumed           : 18383.9 s' in text
    assert '271.5 km' in text
    assert 'Besselian elements' not in text
    assert 'path_map             :' in text
    assert 'Ecl=-05840528' in text
    assert eclipse.path_map.endswith('Mag=0')


def test_solar_eclipse_path_map_and_cond_map(solar_eclipses):
    eclipse = solar_eclipses.get_eclipses(year=-584, month=5, day=28).eclipse(0)
    assert 'Ecl=-05840528' in eclipse.path_map
    troy = montu.Observer(site='troy')
    cond = eclipse.conditions_eclipse(troy)
    assert isinstance(cond, montu.EclipseConditions)
    assert 'cond_map' in cond.__dict__
    assert 'Lat=' in cond.cond_map
    assert 'Lng=' in cond.cond_map
    assert 'LC=1' in cond.cond_map
    assert f'Lat={troy.lat}' in cond.cond_map
    assert f'Lng={troy.lon}' in cond.cond_map


def test_eclipse_conditions_show_details(solar_eclipses, capsys):
    eclipse = solar_eclipses.get_eclipses(year=-584, month=5, day=28).eclipse(0)
    troy = montu.Observer(site='troy')
    cond = eclipse.conditions_eclipse(troy)
    assert cond.show_details() is None
    out = capsys.readouterr().out
    assert 'Eclipse local circumstances' in out
    assert 'Kind                 : total' in out
    assert 'C1 (first contact)' in out
    assert 'cond_map             :' in out
    assert 'LC=1' in out
