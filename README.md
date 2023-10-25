# Montu Python /mnṯw ꜥꜣpp(y)/
## Astronomical ephemerides for the ancient world

<!-- This are visual tags that you may add to your package at the beginning with useful information on your package --> 
[![version](https://img.shields.io/pypi/v/montu?color=blue)](https://pypi.org/project/montu/)
[![downloads](https://img.shields.io/pypi/dw/montu)](https://pypi.org/project/montu/)

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/montu/data/montu-python-logo-complete.png?raw=true" alt="Logo""/></p>
<!-- Fuente: https://symbolikon.com/downloads/montu-egyptian-god/-->

`MontuPython` (transileterated mnṯw ꜥꜣpp(y)) is a Python package intended to compute astronomical ephemerides in the ancient world, thousands of years before present. It was initially designed to compute ephemerides for the ancient Egypt, but it can also be used to study astronomical phenomena in other sites of interest for cultural astronomy (archeoastronomy).

## Download and install

Describe here how the package can be downloaded and install it in
different arquitectures.

If you are using `PyPI` installation it's as simple as:

```
pip install montu
```

## Quick start

In this section you should provide the most simple instructions to use
your package.

For instance:

```
from montu import *
```

or for a safe import:

```
import montu as mn
```

## Basic code examples

A `Google Colab` notebook with some basic commands are available [here](https://colab.research.google.com/drive/1AbvT5u3yf40UPV9QldbsroeG2teUo52p?usp=sharing).

In order to properly use MontuPython you need to import the package:

```python
from montu import *
```

Importing the package will also load some basic data required for calculations, including a complete database with the brightest stars.

### Get an historical date in all calendars and scales

One of the most interesting and basic functionalities of MontuPython is to convert date among 
different type of calendars and astronomical scales.  You may taste these functionalities using:

```python
mtime = MonTime('2501 b.c.e. 01-01 12:00:00')
```

other alternative formats for the same date are:

```python
mtime = MonTime('bce2501-01-01 12:00:00')
mtime = MonTime('-2500-01-01 12:00:00')
print(mtime)
```

whose output will be:

```
Montu Time Object:
--------------------------
Date in proleptic UTC: -2500-01-01 12:00:00.0000
Date in mixed UTC: -2500-01-22 12:00:00:
Date in SPICE format: 2501 B.C. 01-01 12:00:00.00
General:
    Components: [-1, 2500, 1, 1, 12, 0, 0, 0]
    Is bce: True
    Is Julian: True
Uniform scales:
    Terrestrial time:
        tt: -142006202700.32004
        jtd: 807954.6909685181
    UTC time:
        et: -142006262400.00003
        jed: 807953.9999999997
    Delta-t = TT - UTC = 59699.68000000001
Objects:
    Date in datetime64 format: -2500-01-01T12:00:00.000000
    Date in PyPlanet Epoch: 807953.9999999997
    Date in PyEphem Epoch: -2501/1/22 12:00:00
    Date in AstroPy Time: 807954.6909685181
Astronomical properties at Epoch:
    True obliquity of ecliptic: 23:58:33.587
    True nutation longitude: 00:00:10.214
    Greenwhich Meridian Sidereal Time: 18:40:25.323
```

Notice that the date in Gregorian proleptic will be bce2501-01-01 but in mixed calendar will be bce 2501-01-22.

### Position of mars at an historical date and place

For this we need to load precision information about the planets using:

```python
Montu.load_kernels(PRECISION_KERNELS)
```

> **NOTE**: *Precision kernels* are binary files used by the NASA NAIF SPICE information system which is used by `MontuPython` to calculate with the highest precision the position of the planets in the present and in the past. Some of them are large files (in the order of Giga bytes) so you need enough space to download them.

First we need to prepare some basic objects, the Earth (where the observer is), the location on Earth, the time of observation and the object to be observed:

```python
earth = PlanetaryBody('Earth')
tebas = ObservingSite(planet=earth,lon=33,lat=24,height=0)
mtime = MonTime('-2500-01-01 12:00:00')
mars = PlanetaryBody('Mars')
```

Now we can ask for the position of Mars:

```python
mars.calculate_sky_position(mtime,tebas)
```

The result will be:
```
Computing position of body 'mars' at epoch: jtd = 807954.6909685181 
Method 'SPICE':
	Position Epoch: prolectic gregorian -2500-01-01 12:00:00.0000, JED = 807953.9999999997
	Coordinates @ J2000: 
		Equatorial: 12:31:48.754 01:37:12.184
		Ecliptic: 186:39:46.949 04:38:36.308
	Coordinates @ Epoch : 
		Equatorial: 08:32:9.796 24:06:28.555
		Ecliptic: 124:21:21.542 04:39:5.339
	Observing conditions: 
		Distance to site [au]:  0.6604503488431318
		Distance to sun [au]:  1.6261149729988635
		Solar elongation [deg]:  157:49:18.876
		Phase angle [deg]:  13:21:31.981
		Magnitude:  -1.1
	Other properties: 
		Local true sidereal time:  20:52:25.323
		Hour angle @ Epoch:  12:20:15.527
		Local coordinates @ Epoch:  06:11:24.275 -41:38:31.094
```

There are different methods to calculate the position of a planets with MontuPython. You can try all of them:

```python
mars.calculate_sky_position(mtime,tebas,method='all')
```

> **NOTE**: Some methods depend on the availability of a network connection.

The result of the previous command will be:
```
Computing position of body 'mars' at epoch: jtd = 807954.6909685181 
Method 'Horizons':
	Position Epoch: prolectic gregorian -2500-01-01 12:00:00.0000, JED = 807953.9999999997
	Coordinates @ J2000: 
		Equatorial: 12:31:49.147 01:37:6.708
		Ecliptic: 186:39:54.558 04:38:33.609
	Coordinates @ Epoch : 
		Equatorial: 08:32:9.283 24:06:29.556
		Ecliptic: 124:21:14.472 04:39:4.619
	Observing conditions: 
		Distance to site [au]:  0.66052182424896
		Distance to sun [au]:  1.626124866723
		Solar elongation [deg]:  157:47:51.000
		Phase angle [deg]:  13:22:14.880
		Magnitude:  -1.1
	Other properties: 
		Local true sidereal time:  20:52:25.323
		Hour angle @ Epoch:  12:20:16.054
		Local coordinates @ Epoch:  06:11:33.727 -41:38:29.317
Method 'VSOP87':
	Position Epoch: prolectic gregorian -2500-01-01 12:00:00.0000, JED = 807953.9999999997
	Coordinates @ J2000: 
		Equatorial: 12:31:48.360 01:37:33.109
		Ecliptic: 186:39:33.205 04:38:53.193
	Coordinates @ Epoch : 
		Equatorial: 08:32:9.718 24:06:40.920
		Ecliptic: 124:21:17.390 04:39:17.040
	Observing conditions: 
		Distance to site [au]:  0.6604883074760437
		Distance to sun [au]:  1.626107096672058
		Solar elongation [deg]:  157:49:6.094
		Phase angle [deg]:  13:21:48.632
		Magnitude:  -1.13
	Other properties: 
		Local true sidereal time:  20:52:25.323
		Hour angle @ Epoch:  12:20:15.604
		Local coordinates @ Epoch:  06:11:23.906 -41:38:18.685
Method 'SPICE':
	Position Epoch: prolectic gregorian -2500-01-01 12:00:00.0000, JED = 807953.9999999997
	Coordinates @ J2000: 
		Equatorial: 12:31:48.754 01:37:12.184
		Ecliptic: 186:39:46.949 04:38:36.308
	Coordinates @ Epoch : 
		Equatorial: 08:32:9.796 24:06:28.555
		Ecliptic: 124:21:21.542 04:39:5.339
	Observing conditions: 
		Distance to site [au]:  0.6604503488431318
		Distance to sun [au]:  1.6261149729988635
		Solar elongation [deg]:  157:49:18.876
		Phase angle [deg]:  13:21:31.981
		Magnitude:  -1.1
	Other properties: 
		Local true sidereal time:  20:52:25.323
		Hour angle @ Epoch:  12:20:15.527
		Local coordinates @ Epoch:  06:11:24.275 -41:38:31.094
```

### Create a map of the sky

Now we want to create a map of the sky surrounding Aldebaran. For that purpose we get first 
the information about Aldebaran:

```python
aldebaran = ALL_STARS.get_stars(ProperName='Aldebaran')
```

Now we select all stars between magnitude -1 and 5 around Aldebaran in a radius of 10 degrees:

```python
hyades = ALL_STARS.get_stars_area(RA=aldebaran.data.RAJ2000,Dec=aldebaran.data.DecJ2000,
								  radius=10,Mag=[-1,5])
```

Plot the selected stars:

```python
fig,ax = hyades.plot_stars(pad=0.0,labels=False,figargs=dict(figsize=(8,8)))
```

The resulting figure will be:

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/dev/gallery/hyades.png?raw=true" alt="Logo""/></p>

### More complicated example: evolution of pole stars

Choose from database all bright stars that according to [wikipedia](https://en.wikipedia.org/wiki/Pole_star#Precession_of_the_equinoxes) were or will be close to the celestial North pole:

```python
star_names = ('Polaris','Vega','Thuban','Deneb','Alderamin','Kochab')
stars = ALL_STARS.get_stars(ProperName=star_names)
```

Now precess the position of all stars from -20 000 to 20 000 years from 2000:

```python
now = MonTime()
df = pd.DataFrame()
for dt in tqdm.tqdm(np.linspace(-20000*YEAR,20000*YEAR,1000)):
    past = now + dt
    pstars = SkyCoordinates.precess_coordinates(stars.data,past)

    row = dict(tt = past.tt)
    for star in star_names:
        row.update({star:float(pstars[pstars.ProperName == star].DecEpoch)})
    df = pd.concat([df,pd.DataFrame([row])])
```

Now plot declinations as a function of time:

```python
fig,ax = plt.subplots(figsize=(10,6))
for star in star_names:
    ax.plot(df['tt'],df[star],label=star)

ax.legend(loc='lower center',ncol=len(star_names))
ax.set_xlabel("Time [year]")
ax.set_ylabel("Declination [deg]")
ax.axvline(MonTime().tt,color='k',lw=3)
ax.text(0.5,1.01,'Now',ha='center',transform=ax.transAxes)
ax.margins(0)
ax.set_xticks(np.linspace(df['tt'].min(),df['tt'].max(),14))
ax.grid()
MonTime.set_time_ticks(ax)
Montu.montu_mark(ax)
```

The resulting figure will be:

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/dev/gallery/pole-stars.png?raw=true" alt="Logo""/></p>

Check the date when star is close to the pole:

```python
for star in star_names:
    imax = df[star].argmax()
    mtime = MonTime(df.iloc[imax].tt)
    print(f"Star {star} will be the closest to the pole at {mtime.datespice} (declination {D2H(df.iloc[imax][star])})")
```

The output will be:
```
Star Polaris will be the closest to the pole at 2083-11-07 23:25:59.5813 (declination 89:31:51.354)
Star Vega will be the closest to the pole at 11692 B.C. 10-12 08:20:37.619000 (declination 87:24:10.508)
Star Thuban will be the closest to the pole at 2803 B.C. 11-12 23:35:40.709700 (declination 89:54:19.213)
Star Deneb will be the closest to the pole at 14735 B.C. 08-30 11:40:55.71900 (declination 86:57:22.016)
Star Alderamin will be the closest to the pole at 7569-06-13 07:52:18.7355 (declination 87:55:56.040)
Star Kochab will be the closest to the pole at 1041 B.C. 08-30 22:55:00.192000 (declination 83:30:4.833)
```

### Other example scripts

You will find a complete set of runable example notebooks –`Colab` notebooks– in the [Google Drive public repository of `MontuPython`](https://drive.google.com/drive/folders/11L59yZ3A1g1ZT7v_dLDPwLMnRMR-tFgE?usp=sharing).

## Advanced examples

For a fully-fledged working example see `examples/montunctions.ipynb`.

## What's new

Versions 0.7.*:
- A Google Drive repository with examples was created.
- More interesting examples added.
- A new class SkyCoordinates allows to precess coordinates of any set of objects.
- Star catalogue coordinates are now given in J2000.
- Some improvements in the `plot_stars`` method.
- A new `sid_to_tai` for `MonTime` class for calculating the ratio of sidereal to TAI time.

Versions 0.6.*:
- Refactor of MonTime class.
- Corrected a problem with tt during the years 300-1582.
- General cleaning of the package.
- Kernels are separated into basic kernels (which are automatically loaded).
- Julian day dates are rounded-up to 7 figures to avoid representation artifacts.
- Package now load ALL_STARS from import.

Versions 0.5.*:
- New properties (distances, phase angle, magnituded) added to PlanetaryBody class method.
- Solved the DEBUG problem.
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