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
