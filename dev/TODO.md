
# ToDo

Science:

- Calculate Montunctions and draw trajectories of the planet for Tito's PhD thesis.

Urgent: 

- All classes include `verbose` option.

- Update date strings when changing the tt or jd variables.

- Correct for time-of-flight and aberration the coordinates calculated with SPICE.

- Create __str__ method for:    
  * MonTime
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

## Solved

- Check the problems in MonTime:
  * Parsing error for years between 0 and 100 ((SOLVED))
  * Dates between 500 and 1600 ((SOLVED))

