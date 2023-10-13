
# ToDo

Science:

- Calculate Montunctions and draw trajectories of the planet for Tito's PhD thesis.

Urgent: 

- All classes include `verbose` option.

- Correct for time-of-flight and aberration the coordinates calculated with SPICE.

- Create __str__ method for:    
  * MonTime ((SOLVED))
  * Site
  * Planet
  * SkyObject

- Create a comparison matrix between different coordinates calculation methods.

- Study differences between methods along time.

- Precess using FK5 routine of astropy and compare.

- Write routines to transform to local equatorial to local altazimutal and viceversa.

- When computing coordinates of a planet create a DataFrame storing coordinates at different dates.

Others:

- Study in depth manuals for:
  * PyEphem
  * PyPlanets/PyMeeus
  * AstroPy, especially coordinates and time.

- Produce plots with `plotly`.

- Create special Class to compute interesting events for archeoastronomy:
  * Sunrise, Sunset
  * Lunar phases
  * Equinoxes and solstices.
  
- Include a modulus to calculate lunar and solar eclipses. See what is the best package in the market.

- See if we can get information on the topography of a given site and convert it on information about the horizon.

- Create a SkyObject Class to deal with coordinates, aberration correction, etc.
  * Store coordinates as pyplanets_Angle
  * Precess coordinates
  * Correct for aberration

- Test: fecha = -54-06-03, salida de Alkaid por azimuth 18 grados, 09:16:00, Dendera

- Test calculations with historical solar eclipse.

- The package uses four different tools:
  * AstroPy
  * PyEphem
  * PyPlanets
  * SPICE
  * Sky
  * PyMeeus

- Create maps of the sky in different projections
  * Especially maps of the circumpolar region.

- Create maps of the sky inlcuding the Horizon.

- Get a more complete database with stars, including all proper names and Bayer and Flamsteed names with proper characters.

## Useful documentation:

- Reference frames: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/Tutorials/pdf/individual_docs/23_lunar-earth_pck-fk.pdf
- astroquery Horizons documentation: https://astroquery.readthedocs.io/en/latest/api/astroquery.jplhorizons.HorizonsClass.html
- PyEphem quick reference: https://rhodesmill.org/pyephem/quick.html
- PyEphem detailed tutorial: https://rhodesmill.org/pyephem/tutorial
- PyPlanets information in PyPI: https://pypi.org/project/pyplanets/
- PyPlanets repository in GitHub: https://github.com/martin5f/pyplanets
- PyMeeus documentation: https://pymeeus.readthedocs.io/en/latest/
- API reference guide of SPICE: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/index.html
- SOFA, Standards of Fundamental Astronomy: http://www.iausofa.org/current.html
- Useful routines of PyMeeus: https://github.com/architest/pymeeus/blob/master/pymeeus/Coordinates.py
- Examples of PyPlanets: https://github.com/martin5f/pyplanets/blob/master/examples/pymeeus/usage_coordinates.py
- Coordinates transformation: http://star-www.st-and.ac.uk/~fv/webnotes/chapter7.htm
- Specification of times in Horizons: https://ssd.jpl.nasa.gov/horizons/manual.html#time

## Test case

```python
# Planet
mars = PlanetaryBody('Mars')
earth = PlanetaryBody('Earth')

# Site
senenmut = ObservingSite(planet=earth,lon=33,lat=24,height=0)

# Initial epoch
mtime_initial = MonTime('-2500-01-01 12:00:00.00',scale='utc',calendar='proleptic')
print(mtime_initial)

# Initialize quantities for ephemerides calculation
mars.calculate_sky_position(mtime_initial,senenmut,method='SPICE',verbose=1)
```

It should produce:

```
Montu Time Object:
--------------------------
General:
    Calendar: proleptic
    Is bce: True
    Components UTC: [-1, 2500, 1, 1, 12, 0, 0, 0]
Uniform scales:
    Delta-t = TT - UTC = 59699.68000000001
    Terrestrial time:
        tt: -142006202700.32
        jtd: 807954.6909685184
    UTC time:
        et: -142006262400.0
        jed: 807953.9999999999
Strings:
    Date in native format: -2500-01-01 12:00:00.0
    Date in SPICE format: 2501 B.C. 01-01 12:00:00.000000
    Date in mixed calendar: -2501-1-22 12:00:00
Objects:
    Date in datetime64 format: -2500-01-01T12:00:00.000
    Date in datetime format: 2500-01-01 12:00:00
    Date in PyPlanet Epoch: 807953.9999999999
    Date in PyEphem Epoch: -2501/1/22 12:00:00
    Date in AstroPy Time: 807954.6909685184
Astronomical properties at Epoch:
    True obliquity of ecliptic: 23:58:33.587
    True nutation longitude: 00:00:10.214
    Greenwhich Meridian Sidereal Time: 18:40:25.323
Hash: 8781012572531953949

Computing position of body 'mars' at epoch: jtd = 807954.6909685184 
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

## Solved

- Check the problems in MonTime:
  * Parsing error for years between 0 and 100 ((SOLVED))
  * Dates between 500 and 1600 ((SOLVED))

- Update date strings when changing the tt or jd variables ((SOLVED))
