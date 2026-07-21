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
    assert '(alt ' in out and '°, az ' in out
    assert 'cond_map             :' in out
    assert 'LC=1' in out


def test_solar_eclipses_list_heclipses(solar_eclipses):
    rows = solar_eclipses.list_heclipses()
    assert len(rows) == 29
    ids = {row['heclipseid'] for row in rows}
    assert 'amarna-1338bce' in ids
    assert 'thales-585bce' in ids
    assert 'legendary_chinese-2137bce' in ids
    amarna = next(r for r in rows if r['heclipseid'] == 'amarna-1338bce')
    assert amarna['date'] == 'bce 1338-05-14'
    assert 'Akhenaten' in amarna['description']


def test_solar_eclipse_from_heclipseid_catalogue():
    amarna = montu.SolarEclipse('amarna-1338bce')
    assert amarna.heclipseid == 'amarna-1338bce'
    assert amarna.location_id == 'amarna'
    assert amarna.date_key == 'bce 1338-05-14'
    assert amarna.data is not None
    assert int(amarna.data.year) == -1337
    assert 'heclipseid' in str(amarna)


def test_solar_eclipse_from_heclipseid_keyword():
    thales = montu.SolarEclipse(heclipseid='thales-585bce')
    assert thales.location_id == 'miletus'
    cond = thales.conditions_eclipse(montu.Observer(site='miletus'))
    assert cond.kind == 'partial'
    assert cond.magnitude == pytest.approx(0.971, abs=0.01)


def test_solar_eclipse_historical_only_no_conditions():
    legendary = montu.SolarEclipse('legendary_chinese-2137bce')
    assert legendary.data is None
    assert legendary.in_catalogue is False
    assert 'legendary annular eclipse' in legendary.description.lower()
    with pytest.raises(ValueError, match='no NASA catalogue row'):
        legendary.conditions_eclipse(montu.Observer(lat=35, lon=105))


def test_solar_eclipse_unknown_heclipseid():
    with pytest.raises(ValueError, match='unknown historical eclipse id'):
        montu.SolarEclipse('not-a-real-eclipse')


@pytest.fixture(scope='module')
def mars_aldebaran():
    mars = montu.Planet('Mars')
    aldebaran = montu.Stars(
        subset='bright', ProperName='Aldebaran', return_as='Star',
    )
    return mars, aldebaran


def test_conjunction_mars_aldebaran_september_2022(mars_aldebaran):
    mars, aldebaran = mars_aldebaran
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07 12:00:00'),
        observer='geocentric',
    )
    assert conj.is_geocentric
    assert conj.in_conjunction
    assert conj.visible_from_site is None
    assert conj.above_horizon is False  # below geocentric horizon at noon
    assert conj.separation == pytest.approx(4.275, abs=0.05)
    assert conj.position_angle == pytest.approx(169.5, abs=2.0)
    assert len(conj.pairs) == 1
    assert conj.body_conditions[0]['name'] == 'Mars'
    assert conj.body_conditions[0]['phase'] is not None
    assert conj.body_conditions[1]['name'] == 'Aldebaran'


def test_conjunction_is_visible_modes(mars_aldebaran, capsys):
    mars, aldebaran = mars_aldebaran
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07 12:00:00'),
        observer='geocentric',
    )
    site = montu.Observer(lat=10, lon=-75)

    only_time = conj.is_visible(at=montu.Time('2022-09-07 12:00:00'))
    assert only_time.in_conjunction is True
    assert only_time.visible is None
    assert only_time.visible_from_site is None

    only_site = conj.is_visible(from_site=site)
    assert only_site.in_conjunction is True
    # ~07:00 local: bodies above horizon, but Sun not deep enough (< -5°)
    assert only_site.above_horizon is True
    assert only_site.sun_altitude > montu.CONJUNCTION_SUN_MAX_ALTITUDE_DEG
    assert only_site.visible_from_site is False
    assert only_site.visible is False

    both = conj.is_visible(
        from_site=site, at=montu.Time('2022-09-07 12:00:00'),
    )
    assert both.in_conjunction is True
    assert both.visible_from_site is False
    assert both.visible is False

    # Predawn: bodies above horizon and Sun below -5°
    night = conj.is_visible(
        from_site=site, at=montu.Time('2022-09-07 09:00:00'),
    )
    assert night.above_horizon is True
    assert night.sun_altitude < montu.CONJUNCTION_SUN_MAX_ALTITUDE_DEG
    assert night.visible_from_site is True
    assert night.visible is True

    out = capsys.readouterr().out
    assert 'in conjunction' in out
    assert 'Is visible from site=' in out


def test_conjunction_show_details(mars_aldebaran, capsys):
    mars, aldebaran = mars_aldebaran
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07 12:00:00'),
        observer='geocentric',
    )
    conj.show_details()
    out = capsys.readouterr().out
    assert 'Conjunction: Mars–Aldebaran' in out
    assert 'Observer             : geocentric' in out
    assert 'Is visible from site : n/a (geocentric)' in out
    assert 'Angular separation' in out
    assert 'Position angle' in out
    assert 'Phase' in out
    assert 'Angular size' in out
    assert 'Earth distance' not in out
    assert 'Rise (UTC)' not in out
    assert 'Set (UTC)' not in out


def test_conjunction_explorer_finds_september_2022(mars_aldebaran):
    mars, aldebaran = mars_aldebaran
    explorer = montu.ConjunctionExplorer(
        bodies=[mars, aldebaran], maxseparation=5,
    )
    conjs = explorer.search(
        start=montu.Time('2022-09-01'),
        end=montu.Time('2022-10-01'),
        observer='geocentric',
        verbose=False,
    )
    assert len(conjs) == 1
    conj = conjs[0]
    assert isinstance(conj, montu.Conjunction)
    assert conj.is_geocentric
    assert conj.in_conjunction
    assert conj.visible_from_site is None
    assert abs(conj.mtime.jed - montu.Time('2022-09-07').jed) < 1.5
    assert conj.separation == pytest.approx(4.275, abs=0.05)


def test_conjunction_explorer_topocentric(mars_aldebaran, capsys):
    mars, aldebaran = mars_aldebaran
    site = montu.Observer(lat=6, lon=-75)
    explorer = montu.ConjunctionExplorer(
        bodies=[mars, aldebaran], maxseparation=5,
    )
    conjs = explorer.search(
        start=montu.Time('2022-09-01'),
        end=montu.Time('2022-10-01'),
        observer=site,
        verbose=False,
    )
    assert len(conjs) == 1
    conj = conjs[0]
    assert conj.is_geocentric is False
    assert conj.body_conditions[0]['el'] is not None
    assert conj.above_horizon is True
    assert all(bc['above_horizon'] for bc in conj.body_conditions)
    # Closest approach is near local sunrise → Sun above -5°
    assert conj.sun_altitude is not None
    assert conj.sun_altitude > montu.CONJUNCTION_SUN_MAX_ALTITUDE_DEG
    assert conj.visible_from_site is False
    conj.show_details()
    out = capsys.readouterr().out
    assert 'Is visible from site : no' in out
    assert 'Sun < -5°' in out
    assert 'Sun altitude' in out
    assert 'above horizon: yes' in out
    assert 'Rise (UTC)' in out
    assert 'Set (UTC)' in out


def test_conjunction_alias_typo():
    assert montu.Conjuntion is montu.Conjunction


def test_conjunction_explore_lapse_mars_aldebaran(mars_aldebaran, capsys):
    mars, aldebaran = mars_aldebaran
    site = montu.Observer(lat=6, lon=-75)
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07'),
        observer=site,
    )
    lapse = conj.explore_lapse(verbose=False)
    assert lapse is not None
    start, end = lapse
    assert isinstance(start, montu.Time)
    assert end.jed > start.jed
    assert 8.0 <= (end.jed - start.jed) <= 11.0
    assert conj.lapse.duration_days == pytest.approx(end.jed - start.jed)
    assert conj.lapse.start_separation == pytest.approx(5.0, abs=0.05)
    assert conj.lapse.end_separation == pytest.approx(5.0, abs=0.05)
    out = capsys.readouterr().out
    assert out == ''


def test_conjunction_explore_lapse_no_conjunction(mars_aldebaran, capsys):
    mars, aldebaran = mars_aldebaran
    site = montu.Observer(lat=6, lon=-75)
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=1,
        mtime=montu.Time('2022-09-07'),
        observer=site,
    )
    assert conj.explore_lapse(verbose=True) is None
    out = capsys.readouterr().out
    assert 'No hay conjunción en esas condiciones.' in out


def test_conjunction_plot_lapse_builds_figure(mars_aldebaran):
    pytest.importorskip('plotly')
    mars, aldebaran = mars_aldebaran
    site = montu.Observer(lat=6, lon=-75)
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07'),
        observer=site,
    )
    lapse = conj.explore_lapse(verbose=False)
    fig = conj.plot_lapse(
        lapse[0], lapse[1], step_hours=6, show=False, return_fig=True,
    )
    assert fig is not None
    assert len(fig.data) >= 3


def test_conjunction_plot_map_builds_figure(mars_aldebaran):
    pytest.importorskip('plotly')
    mars, aldebaran = mars_aldebaran
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=5,
        mtime=montu.Time('2022-09-07'),
        observer='geocentric',
    )
    fig = conj.plot_map(show=False, return_fig=True)
    assert fig is not None
    assert len(fig.data) >= 2
    assert 'Conjunction map' in fig.layout.title.text


def test_conjunction_plot_map_skips_when_not_in_conjunction(mars_aldebaran):
    pytest.importorskip('plotly')
    mars, aldebaran = mars_aldebaran
    conj = montu.Conjunction(
        bodies=[mars, aldebaran],
        maxseparation=1,
        mtime=montu.Time('2022-09-07'),
        observer='geocentric',
    )
    assert conj.plot_map(show=False, return_fig=True) is None


def test_conjunction_plot_map_ra_wrap_near_zero():
    pytest.importorskip('plotly')
    bodies = [
        montu.Planet('Mars'),
        montu.Planet('Jupiter'),
        montu.Planet('Saturn'),
    ]
    conj = montu.Conjunction(
        bodies=bodies,
        maxseparation=10,
        mtime=montu.Time('-0005-02-20 07:00:00', calendar='mixed'),
        observer='geocentric',
    )
    fig = conj.plot_map(show=False, return_fig=True)
    assert fig is not None
    stars = next(t for t in fig.data if t.name == 'Stars')
    assert len(stars.x) > 10
    ra_lo, ra_hi = sorted(fig.layout.xaxis.range)
    assert all(ra_lo <= float(x) <= ra_hi for x in stars.x)
    assert max(stars.x) - min(stars.x) < 90.0
