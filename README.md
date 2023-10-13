# Montu Python /mnṯw ꜥꜣpp(y)/
## Astronomical ephemerides for the ancient world

<!-- This are visual tags that you may add to your package at the beginning with useful information on your package --> 
[![version](https://img.shields.io/pypi/v/montu?color=blue)](https://pypi.org/project/montu/)
[![downloads](https://img.shields.io/pypi/dw/montu)](https://pypi.org/project/montu/)

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/montu/data/montu-python-logo-complete.png?raw=true" alt="Logo""/></p>
<!-- Fuente: https://symbolikon.com/downloads/montu-egyptian-god/-->

`MontuPython` (transileterated is a Python package intended to compute astronomical ephemerides in the ancient world, thousands of years before present.
It was initially designed to compute ephemerides for the ancient Egypt, but it can also be used to study astronomical 
phenomena in other sites of interest for cultural astronomy (archeoastronomy).

## Download and install

Describe here how the package can be downloaded and install it in
different arquitectures.

If you are using `PyPI` installation it's as simple as:

```
pip install montu
```

You can also test the unstable version of the package with:

```
pip install -i https://test.pypi.org/simple/ montu
```

## Quick start

In this section you should provide the most simple instructions to use
your package.

For instance:

```
import montu
print(montu.version())
```

## Code examples

For a fully-fledged working example see `examples/montunctions.ipynb`.

In order to properly use MontuPython you need to import the package and load the required
data:

```python
import montu
from montu import *
Montu.load_kernels()
allstars=Stars()
```

It may take a long time to download the planetary kernels, so be patient! 

### Position of mars at an historical date

First we need to prepare some basic objects, the Earth (where the observer is), the location on Earth, the time of 
observation and the object to be observed:

```python
earth = PlanetaryBody('Earth')
tebas = ObservingSite(planet=earth,lon=33,lat=24,height=0)
mtime = MonTime('bce2501-01-01 12:00:00')
mars = PlanetaryBody('Mars')
```

Now we can ask for the position of Mars:

```python
mars.calculate_sky_position(mtime,tebas)
```

The result will be:
```
Computing position of body 'mars' at epoch: jtd = 807954.6909685184 
Updating orientation of site (old time 2000-01-01 11:58:56.126200, new time 2501 B.C. 01-01 12:00:00.000000)
Method 'SPICE':
	Coordinates @ J2000: 
		Equatorial: (12.0, 31, 48.75374471536361) (1.0, 37, 12.183893707210274)
		Ecliptic: (186.0, 39, 46.94915035185886) (4.0, 38, 36.30772065411392)
	Coordinates @ Epoch : 
		Equatorial: (8.0, 32, 9.7958860775843) (24.0, 6, 28.554848777837947)
		Ecliptic: (124.0, 21, 21.541927979009188) (4.0, 39, 5.338625036125535)
	Local true sidereal time:  (20.0, 52, 25.32259203620299)
	Hour angle @ Epoch:  (12.0, 20, 15.526705958618692)
	Local coordinates @ Epoch:  (6.0, 11, 24.27528943044642) (-41.0, 38, 31.093764289017827)
```

There are different methods to calculate the position of a planets with MontuPython. You can try all of them:

```python
mars.calculate_sky_position(mtime,tebas,method='all')
```

> **NOTE**: Some methods depend on the availability of a network connection.

The result of the previous command will be:
```
Computing position of body 'mars' at epoch: jtd = 807954.6909685184 
Method 'Horizons':
	Coordinates @ J2000: 
		Equatorial: (12.0, 31, 49.14719999999818) (1.0, 37, 6.7080000000001405)
		Ecliptic: (186.0, 39, 54.55755751664583) (4.0, 38, 33.60945464297373)
	Coordinates @ Epoch : 
		Equatorial: (8.0, 32, 9.283199999998999) (24.0, 6, 29.555999999998903)
		Ecliptic: (124.0, 21, 14.472000000019989) (4.0, 39, 4.619159999998601)
	Local true sidereal time:  (20.0, 52, 25.33680635060051)
	Hour angle @ Epoch:  (12.0, 20, 16.053606350601513)
	Local coordinates @ Epoch:  (6.0, 11, 33.727199999998945) (-41.0, 38, 29.317200000006665)
Method 'VSOP87':
	Coordinates @ J2000: 
		Equatorial: (12.0, 31, 48.35959491216798) (1.0, 37, 33.10917010023772)
		Ecliptic: (186.0, 39, 33.204611828884936) (4.0, 38, 53.19324054798258)
	Coordinates @ Epoch : 
		Equatorial: (8.0, 32, 9.71805413547699) (24.0, 6, 40.9201083707282)
		Ecliptic: (124.0, 21, 17.39024402412724) (4.0, 39, 17.039572599530572)
	Local true sidereal time:  (20.0, 52, 25.322505532880086)
	Hour angle @ Epoch:  (12.0, 20, 15.604451397396701)
	Local coordinates @ Epoch:  (6.0, 11, 23.905564366115826) (-41.0, 38, 18.685147982311605)
Method 'SPICE':
	Coordinates @ J2000: 
		Equatorial: (12.0, 31, 48.75374471536361) (1.0, 37, 12.183893707210274)
		Ecliptic: (186.0, 39, 46.94915035185886) (4.0, 38, 36.30772065411392)
	Coordinates @ Epoch : 
		Equatorial: (8.0, 32, 9.7958860775843) (24.0, 6, 28.554848777837947)
		Ecliptic: (124.0, 21, 21.541927979009188) (4.0, 39, 5.338625036125535)
	Local true sidereal time:  (20.0, 52, 25.32259203620299)
	Hour angle @ Epoch:  (12.0, 20, 15.526705958618692)
	Local coordinates @ Epoch:  (6.0, 11, 24.27528943044642) (-41.0, 38, 31.093764289017827)
```

## What's new

Versions 0.5.*:
- First fully-fledge working script, yay!
- MonTime class is now able to produce dates in all relevant calendars and formats.
- Verified ephemerides calculations for planets.
- Major improvements of all functionalities.

Version 0.1.*:

- First classes created and tested with study case.
- A proper identifying image of the project has been found.
- The project is started!

------------

This package has been designed and written by Jorge I. Zuluaga with the historical advise of egyptologist Francisco "Tito" Vivas (C) 2023