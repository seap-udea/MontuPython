# What's New

This file collects the release notes and the main changes in MontuPython.

## Version 0.30.0 (breakthrough)

- **Sothic (Egyptian civil) calendar** — terminology and API updated from *caniucular* to **sothic** across the library, examples, and documentation (`calendar='sothic'`, `Time.readable.datesot`, `Time.parse_datesot()`).
- **Mixed-year sothic input** — civil dates accept Horus year (`[hrw N]`), mixed historical BCE/CE (`[bce YYYY]`), astronomical year (`[-YYYY]`), and CE tags (`[ce YYYY]` / `[YYYY]`); see `montu/tests/test_sothic_mixed_year.py`.
- **`Time.parse_datesot()`** — parse `[hrw YEAR] MONTH SEASON DAY` strings into Horus year, month, season, and day components.
- **Historical dates catalogue** — updates in `montu/data/historical_dates.json` for sothic-era entries and civil-calendar cross-references.
- **Example notebooks** — Egyptian calendar, MonTime, Heliacal Rises, and Basic Functions refreshed for sothic naming and mixed-year conventions.
- **Tests** — `montu/tests/test_sothic_mixed_year.py` and docstring-example coverage extended for the new civil-date formats.

## Version 0.22.1 (major)

- **Solar eclipses** — new classes in `montu/phenomena.py` for local circumstances from polynomial Besselian elements:
  - **`SolarEclipses`** — catalogue loader and filter (same conventions as `Stars.get_stars` / `get_eclipses`); alias `get_eclipse`.
  - **`SolarEclipse`** — single-eclipse object with `conditions_eclipse(observer)` (kind, visibility, magnitude, obscuration, contacts C1–C4, umbral duration, solar altitude) for any `Observer` / city.
  - Bundled NASA Five Millennium Catalog of Solar Eclipses (Espenak & Meeus): **11 898** eclipses from **−1999 to +3000**, file `montu/data/nasa_5mcse_besselian.csv`.
  - **`SolarEclipse.__str__`** — human-readable catalogue summary with field explanations.
- **Example notebooks**
  - **`examples/MontuPython-HeliacalRises.ipynb`** — full step-by-step tutorial of heliacal-rise calculation (Ptolemy, Schaefer 1985/1987, Belokrylov 2011) from Thebes / Sirius.
  - **`examples/MontuPython-SolarEclipses.ipynb`** — catalogue search, modern validation (Dallas 2024), and historical survey of eclipses visible from Thebes (case study −1257 July 27).
- **Tests** — `montu/tests/test_phenomena.py` extended for catalogue loading and Dallas 2024 local totality.

## Version 0.21.6 (minor)

- Fixed time calculation errors in the azimuthal projection map algorithm.
- Introduced Ptolemy's arcus visionis algorithm for heliacal rising calculation.

## Version 0.21.5 (minor)

- Localization and consistency updates to support the Desktop ES/EN bilingual experience.
- Terminology alignment for seasons, heliacal-rise wording, and related Spanish text conventions.
- Supporting notes/content refresh to keep library and Desktop documentation aligned.

## Version 0.21.4 (major)

- **`montu/observer.py`** — predefined ancient-world observing sites bundled in `montu/data/locations.json` (same catalogue as MontuPython Desktop).
- **`Observer(site=…)`** — construct an observer from a catalogue id (e.g. `montu.Observer(site='memphis')`); coordinates and altitude are taken from the JSON entry.
- **`Observer.list()`** — return the list of available site ids (`'thebes'`, `'memphis'`, …).
- **`Observer.list(details=True)`** — return full metadata dicts (`id`, `name`, `lat`, `lon`, `alt_m`, `region`, `era`, `description`).
- New attributes **`site_id`** and **`site_name`** when the observer is created via `site=…`.

## Version 0.21.3 (major)

- **`montu/phenomena.py`** — new module for celestial phenomena predicted from observational visibility models.
- **`HeliacalRise`** — class to search for heliacal-rise mornings over a date interval at an observer site; supports **Schaefer 1985**, **Schaefer 1987**, and **Belokrylov et al. 2011** models with configurable extinction, limiting magnitude, solar depression, and twilight-scan parameters.
- **`heliacal_rise()`** — convenience wrapper that runs a search and prints a formatted table (model source, local time, body/Sun altitude and azimuth, and model visibility quantities).
- **Stars and planets** — accepts a precessed `Stars` row or a `Sebau` planet; catalogue V magnitude for stars and ephemeris magnitude for planets.
- **`HELIACAL_RISE_MODELS`** / **`HELIACAL_RISE_SOURCES`** — model identifiers and bibliographic source strings exported at package level.
- **Tests** — `montu/tests/test_phenomena.py` covers all three models, the function wrapper, and the third *apokatastasis* Sirius case (Thebes, Schaefer 1987).

## Version 0.21.2 (minor)

- **Example notebooks** — `%matplotlib inline` moved to the `import montu` cell; `%load_ext autoreload` / `%autoreload 2` commented with a note that they fail on Google Colab (Python 3.12+, removed `imp` module). Applied in `examples/` and `docs/examples/`.
- **`README.ipynb`** — pole-star precession loop prints a short start/finish message instead of a tqdm progress bar (cleaner output in Jupyter and Colab).
- **`polar_sky_map`** — figure layout uses `autosize` rather than a fixed pixel height so azimuthal maps resize better when embedded in MontuPython Desktop.

## Version 0.21.1 (major)

- **`montu/maps.py`** — new module that unifies Plotly sky-map plotting: equatorial Mercator and azimuthal (polar) projections live together instead of being split across `stars.py` and a separate helper module.
- **`polar_sky_map`** — builds north- and south-hemisphere azimuthal maps from a precessed catalogue: limiting magnitude, optional **horizon** (elevation 0° with azimuth marks) and **ecliptic**, solar-system bodies, **meridian view**, and **LST** in the title. High-level API: `polar_sky_map(...)`; per-hemisphere: `polar_sky_map_figure(...)`.
- **Constellation asterism sets** — choose among `iau`, `egyptian_ancient`, and `egyptian_dendera` via `constellation_set` (Stellarium sky-culture stick figures and name files shipped in `montu/data/`). New helpers: `CONSTELLATION_SET_IDS`, `constellation_data_files()`, `parse_constellation_names()`.
- **`mercator_sky_map`** — moved into `montu/maps.py`; still exported from `montu`, `montu.maps`, and `montu.stars` so existing notebooks and `from montu.stars import mercator_sky_map` keep working.
- **`Stars.polar_sky_map()`** — convenience wrapper on a precessed `Stars` catalogue (mirrors `Stars.mercator_sky_map()`).
- **`local_solar_to_utc_datepro()`** — convert local solar time at an observer longitude to a proleptic UTC `datepro` string (used by polar sky maps and MontuPython Desktop).

## Version 0.21.0 (breakthrough)

- **SPICE removed from the active package** — kernels are no longer loaded on import; time conversions use NumPy and PyMeeus instead of `spiceypy`. `spiceypy` was dropped from install dependencies.
- **Archived SPICE code** — earlier cycle modules (`__cycle_*.py`) and NAIF kernels moved to `contrib/montu-deprecated/` (outside the installable package).
- **`load_historical_dates()`** — historical event catalogue shipped in `montu/data/historical_dates.json`; shared by the library, Desktop, and MontuApp.
- **New example** — `examples/MontuPython-EgyptianCalendar.ipynb` compares known Egyptian civil dates from historical events against MontuPython.
- **`Time.__str__`** — clearer readable output with attribute paths (e.g. `.readable.datepro`) and weekday merged into the components line.
- **CLI module renamed** — `montu/imontu_cli.py` → `montu/cli.py` (entry points `imontu` and `montu-gui` unchanged).
- **Release pipeline** — `release-pipeline.sh` at the repository root automates version bump, notebook execution, docs, PyPI upload, and optional Desktop tagging.
- **Developer environment** — standard virtualenv directory is now `.venv` (was `.montuenv`).
- **Documentation** — corrected Read the Docs URL (`montupython.readthedocs.io`); README examples split into tutorial and advanced sections; logo images use WebP.
- **Metadata** — CITATION.cff and `.zenodo.json` authorship updates.

## Version 0.20.x (breakthrough)

- Renamed **`Util.dec2hex`** to **`Util.dec2sex`** (decimal → sexagesimal DMS/HMS); added inverse **`Util.sex2dec`**.
- Removed the old `dec2hex` / `hex2dec` aliases; updated tests, notebooks, and examples to the new names.
- New top-level alias **`S2D`** for sexagesimal-to-decimal conversion; **`D2H`** now points to `dec2sex`.
- **`mercator_sky_map`** exported at package level (`montu.mercator_sky_map` and `from montu.stars import mercator_sky_map`) for notebooks and Colab workflows.
- Version bump to **0.20.0** with updated README, badges, AI assistance disclosure, and Zenodo/Citation metadata.
- **MontuPython Desktop** (`montu_gui/`) — PySide6 GUI with Calendar Calculator and Seasons & Lunar Phases modules.
- Refactored and standardized example notebooks; cleaned outdated development artefacts.
- Expanded pytest coverage for notebook workflows and docstring examples.
- Fixed `Stars.get_stars_around` to accept one-element pandas Series as centre coordinates.
- Fixed `Util._linear_map` (missing slope) affecting `Stars.plot_stars`.
- Fixed `Observer.get_local_time` and Julian-day inputs from `Sun.when_is_twilight`.
- Added `Stars.value_for` and `Stars.scalar` for robust scalar access from star catalogues.
- `where_in_space` now returns a `Stars` object for API consistency.
- `Stars.mercator_sky_map` layout now sets `autosize=True` so sky maps resize correctly when embedded in the Desktop GUI.

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
