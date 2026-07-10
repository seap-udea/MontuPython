# What's New

This file collects the release notes and the main changes in MontuPython.

## Version 0.10.*

- Expanded and cleaned the API documentation for the core modules with improved NumPy-style docstrings.
- Stabilized the Read the Docs and Sphinx build, including autodoc imports and rendered examples.
- Standardized the notebooks in examples with a shared header, Colab badge, logo, setup cells, footer, and normalized first-level titles.
- Improved how example notebooks are exposed in the documentation and removed duplicated-looking entries in the examples navigation.
- Added a pytest-based test suite derived from docstring examples and workflows from the example notebooks.
- Added dedicated make targets for the test suite and separated smoke, docstring, notebook, and structure checks.
- Moved the tests into montu/tests and configured them to be distributed with the package.
- Fixed issues uncovered by the new tests in util, stars, sebau, and time-related workflows.
- Removed warnings during test execution, including catalogue-loading dtype warnings and deprecated PyEphem event access.
- Added the imontu package script to verify an installation and run the packaged tests.

## Versions 0.9.*

- Major refactoring of inner machinery.
- After many tests we decided to work with the algorithms in PyEphem, with some support from PyMeeus and PyPlanets, and minor usage of spiceypy.
- Improvements on class Time to make it faster.
- Generic class Sebau for celestial objects.
- New classes Sun and Moon.
- Complete example of Montunctions adapted to the new package architecture.
- The general design of the package is maintained.
- New functions for calculating twilight time.
- Juanita Agudelo created a basic Dash app for interacting with the package; it is now included in the repository, although not yet in the public package.
- MontuApp has been developed in depth.
- Lunar phases added.
- Conversion gregorian/julian <-> civil egyptian implemented.

## Versions 0.8.*

- Major refactoring of classes.
- New class Heka, intended to perform calculations.
- New plot of positions in a stereographic projection with respect to the horizon.
- Horizontal coordinates of planets and stars can now be computed.
- New stellar catalogue incorporated into the package.
- Classes were divided into modules such as stars, planets, and observer.
- README was converted into a notebook.

## Versions 0.7.*

- A Google Drive repository with examples was created.
- More interesting examples were added.
- A new class SkyCoordinates allows precession of coordinates for arbitrary sets of objects.
- Star catalogue coordinates are now given in J2000.
- Improvements in the plot_stars method.
- A new sid_to_tai method was added to MonTime for calculating the ratio of sidereal to TAI time.

## Versions 0.6.*

- Refactor of MonTime class.
- Corrected a problem with tt during the years 300-1582.
- General cleaning of the package.
- Kernels were separated into basic kernels, which are automatically loaded.
- Julian day dates are rounded to 7 figures to avoid representation artifacts.
- The package now loads ALL_STARS on import.

## Versions 0.5.*

- New properties such as distances, phase angle, and magnitude were added to the PlanetaryBody class method.
- Solved the DEBUG problem.
- First fully fledged working script.
- MonTime can now produce dates in all relevant calendars and formats.
- Verified ephemerides calculations for planets.
- Major improvements of all functionalities.

## Version 0.1.*

- First classes created and tested with the study case.
- A proper identifying image of the project was found.
- The project started.
