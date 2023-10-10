
# ToDo

- Create a SkyObject Class to deal with coordinates, aberration correction, etc.
  * Store coordinates as pyplanets_Angle
  * Precess coordinates
  * Correct for aberration

- Correct for time-of-flight and aberration the coordinates calculated with SPICE.

- Check the problems in MonTime:
  * Parsing error for years between 0 and 100
  * Dates between 500 and 1600

- Create __str__ method for:    
  * MonTime
  * Site
  * Planet
  * SkyObject

- Test: fecha = -54-06-03, salida de Alkaid por azimuth 18 grados, 09:16:00, Dendera

- Create a comparison matrix between different coordinates calculation methods.

- Study differences between methods along time.

- Test calculations with historical solar eclipse.

- Precess using FK5 routine of astropy.

- Write routines to transform to local equatorial to local altazimutal and viceversa.

- When computing coordinates of a planet create a DataFrame storing coordinates at different dates.

- The package uses four different tools:
  * AstroPy
  * PyEphem
  * PyPlanets
  * SPICE
  * Sky
  * PyMeeus

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