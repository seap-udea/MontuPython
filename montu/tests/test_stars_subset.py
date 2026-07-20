import pytest
import montu

def test_stars_all():
    allstars = montu.Stars()
    # It should have a large number of stars (around 100,000+)
    assert allstars.number > 100000

def test_stars_visible():
    visible = montu.Stars(subset='visible')
    # It should have fewer stars than the full catalogue, but still a reasonable number (around 9,000+)
    assert visible.number < 100000
    assert visible.number > 5000
    # Let's verify maximum magnitude is indeed around 6.5
    assert visible.data.Vmag.max() <= 7.0

def test_stars_bright():
    bright = montu.Stars(subset='bright')
    # It should have even fewer stars
    assert bright.number < 5000
    assert bright.number > 100
    # Let's verify maximum magnitude is around 3.5
    assert bright.data.Vmag.max() <= 4.0

def test_stars_constructor_filtering():
    # Test instantiating single star by proper name
    sirius = montu.Stars(ProperName='Sirius')
    assert sirius.number == 1
    assert sirius.data.iloc[0].ProperName == 'Sirius'

    # Test multiple criteria (e.g. Orion stars with Vmag in [-3, 1])
    orion = montu.Stars(Vmag=[-3, 1], Constellation='Ori')
    assert orion.number > 0
    assert (orion.data.Constellation == 'Ori').all()
    assert orion.data.Vmag.max() <= 1.0

    # Test combining subset with filtering options
    orion_visible = montu.Stars(subset='visible', Vmag=[-3, 1], Constellation='Ori')
    assert orion_visible.number == orion.number


def test_stars_return_as_star_returns_single_object():
    spica = montu.Stars(subset="bright", ProperName="Spica", return_as="Star")
    assert isinstance(spica, montu.Star)
    assert spica.name == "Spica"


def test_get_stars_return_as_star_returns_single_object():
    spica = montu.Stars(subset="bright").get_stars(
        ProperName="Spica", return_as="Star"
    )
    assert isinstance(spica, montu.Star)
    assert spica.name == "Spica"


def test_get_stars_return_as_star_returns_name_dict_for_multiple():
    stars = montu.Stars(subset="bright").get_stars(
        Constellation="Ori", Vmag=[-3, 1], return_as="Star"
    )
    assert isinstance(stars, dict)
    assert stars
    assert all(isinstance(star, montu.Star) for star in stars.values())
    assert "Betelgeuse" in stars or "Rigel" in stars


def test_stars_return_as_star_empty_filter_raises():
    with pytest.raises(ValueError, match="No stars matched"):
        montu.Stars(subset="bright").get_stars(
            ProperName="NotARealStar", return_as="Star"
        )


def test_stars_conditions_in_sky_returns_dataframe(thebes_observer):
    epoch = montu.Time("bce 1500-01-01 12:00:00", calendar="proleptic")
    orion = montu.Stars(subset="visible", Vmag=[-2, 4], Constellation="Ori")
    conditions = orion.conditions_in_sky(at=epoch, observer=thebes_observer)
    assert conditions is not orion.data
    assert len(conditions) == orion.number
    for col in (
        "ha",
        "rise_time",
        "set_time",
        "transit_time",
        "transit_el",
        "elongation",
        "is_circumpolar",
        "is_neverup",
    ):
        assert col in conditions.columns


def test_stars_conditions_in_sky_matches_single_star(thebes_observer):
    epoch = montu.Time("bce 1500-01-01 12:00:00", calendar="proleptic")
    row = montu.Stars(subset="visible", ProperName="Rigel").data.iloc[0]
    star = montu.Star(row)
    star.conditions_in_sky(at=epoch, observer=thebes_observer)

    batch = montu.Stars(subset="visible", Vmag=[-2, 4], Constellation="Ori")
    conditions = batch.conditions_in_sky(at=epoch, observer=thebes_observer)
    rigel = conditions[conditions.ProperName == "Rigel"].iloc[0]

    assert rigel.rise_time == pytest.approx(star.condition.rise_time)
    assert rigel.set_time == pytest.approx(star.condition.set_time)
    assert rigel.transit_time == pytest.approx(star.condition.transit_time)
    assert rigel.transit_el == pytest.approx(star.condition.transit_el)


def test_stars_conditions_in_sky_inplace(thebes_observer):
    epoch = montu.Time("bce 1500-01-01 12:00:00", calendar="proleptic")
    orion = montu.Stars(subset="visible", Vmag=[-2, 4], Constellation="Ori")
    result = orion.conditions_in_sky(at=epoch, observer=thebes_observer, inplace=True)
    assert result is None
    assert "transit_el" in orion.data.columns


def test_stars_conditions_in_sky_requires_observer():
    epoch = montu.Time("bce 1500-01-01 12:00:00", calendar="proleptic")
    orion = montu.Stars(subset="visible", Vmag=[-2, 4], Constellation="Ori")
    with pytest.raises(ValueError, match="valid montu.Observer"):
        orion.conditions_in_sky(at=epoch, observer=None)
