# What's New

This file collects the release notes and the main changes in MontuPython.

## Version 0.42.0 (major)

- **Cached stellar catalogue** — when loading the stellar catalogue in a session, the package check if the catalogue was already loaded. If so, it will return the cached catalogue instead of loading it again.
- **Senenmut's "Montunction"** — added the BCE 1466-12-06 Mars–Aldebaran historical conjunction to `montu/data/historical-conjunctions.json`, featuring its occurrence near a marsticio and winter solstice according to Francisco Vivas' research.
- **`montu.Horizon` & `Observer.horizon_profile()`** — computes the real topographic horizon profile (elevation vs. azimuth) for any observer by automatically downloading, mosaicking, and caching Copernicus GLO-30 DEM tiles.
- **Smooth horizon profiles** — implemented bilinear interpolation (`scipy.ndimage.map_coordinates`) when sampling elevations from the DEM, eliminating stair-step artifacts at high resolutions.
- **`Horizon.plot_horizon()`** — interactive Plotly chart of the local topography. Constrain the view using `az_center`, `az_delta`, and `elev_view`.
- **`Horizon.plot_horizon(at=...)`** — new functionality that superimposes the computed terrain profile with the visible sky (stars scaled by magnitude and colored by B-V index, asterism lines, and names) at a specific `montu.Time` when `at` is provided. Correctly occludes stars behind local mountains.
- **Horizon plot UI** — dynamic chart titles showing the observer's local time and the UTC date in SPICE format, with automatically adjusted top margins to prevent layout overlaps.
- **Horizon class representation** — clean `__repr__` and detailed `__str__` showing calculation status, elevation range, and the parameters used during the grid scan (`max_dist`, `az_step`, `coarse_step`).
- **`Horizon.plot_map()`** — interactive 2D map method that plots the topographical contour of the observer's horizon on top of geographic map tiles, prominently featuring the observer's location marker.
- **Detailed celestial overlay in `plot_horizon`** — extended the sky plot to dynamically include the Sun (represented with ☀️), IAU constellation boundaries, asterism connections, and star/constellation name labels. Added dedicated toggles (`show_boundaries`, `show_asterism`, `show_starnames`, `show_constname`) to selectively filter elements.
- **Enhanced horizon visual realism** — updated the horizon profile in `plot_horizon` to use an earth ochre color and a highly opaque fill (`alpha=0.85`), effectively and realistically occluding celestial bodies that are positioned behind the mountains.
- **Observer Location Bugfixes** — resolved an issue where manually calculated horizons in Montu Desktop omitted the `site_name` and produced generic titles. Additionally, updated the observer horizon example notebook with map plotting demonstrations and precise sunrise calculation routines.
- **MontuPython Desktop UI Improvements** — adjusted the compact sidebar layout and scaled the logo to prevent clipping on small screens. The sidebar now correctly accommodates longer translated strings (like Spanish) without overflowing, and the configuration buttons' labels were shortened. The home page modules list was updated to a space-efficient 2-column layout, and the module order was reorganized.
- **DEM Cache Management** — centralized the DEM cache folder to always reside next to the installation directory (`dem_cache_dir()`). Added a dedicated "Clear cache" button to the home page's Configuration section to allow users to easily delete downloaded topographic tiles and free up disk space.

## Version 0.41.0 (major)
- **Pinned astronomy stack** — `ephem==4.2.1`, `pymeeus==0.5.12`, `pyplanets==0.4.2` (reproducible ephemerides; stored in `requirements-astronomy.txt`).
- **Organic regression snapshots** — `montu/tests/test-planetary-ephemeris-organic.csv` and `montu/tests/test-stellar-positions-organic.csv` regenerated with this stack (`make organic-snapshots`).
- **`Conjunction`** — new class in `montu/phenomena.py` to evaluate angular conjunctions of two or more sky bodies (planets, stars, Moon, Sun) at one epoch; reports maximum pairwise separation, per-pair geometry, rise/set, phase, and angular size (planets).
- **`ConjunctionExplorer`** — scan a date interval for local separation minima below a threshold; returns fully computed `Conjunction` objects sorted by epoch.
- **`Conjuntion`** — backward-compatible alias for the common misspelling of `Conjunction`.
- **Site visibility** — `visible_from_site`, `is_visible(from_site=…, at=…)`, and `CONJUNCTION_SUN_MAX_ALTITUDE_DEG` (−5°): conjunction is *visible* when all bodies are above the horizon and the Sun is below −5°; geocentric runs mark visibility as `n/a`.
- **`explore_lapse()`** — find the UTC interval during which the group stays within `maxseparation` around the reference day; reports start/end separations and local times.
- **`plot_lapse()`** — Plotly three-panel chart (separation, body elevations, Sun altitude) with green bands for site-visible intervals.
- **`plot_map()`** — Plotly equatorial sky map centred on the geometric mean of the body directions, with stars from the visible catalogue (same spirit as `Stars.plot_stars` / `mercator_sky_map`); options `mag_plotlimit` (default 5.0) and `mag_namelimit` (default 3.5); only plots when `in_conjunction` is true.
- **`Conjunction.show_details()`** — formatted report (epoch, observer, separation, visibility, per-body sky conditions, pairwise table).
- **`Util.print_dict()`** — print a mapping as a Key | Value table with nested sub-tables for lists of dicts (used in conjunction visibility examples); formats `Time`, `Observer`, booleans, and Julian dates.
- **`ConjunctionExplorer.search(..., verbose=False)`** — progress bar off by default (quiet output for notebooks and Read the Docs).
- **Tests** — `montu/tests/test_phenomena.py` extended with Mars–Aldebarán (September 2022) reference cases, explorer search, lapse edges, `plot_lapse`, and `plot_map`; `montu/tests/test_util_print_dict.py` for `Util.print_dict`.
- **`examples/MontuPython-Conjunctions.ipynb`** — tutorial for `Conjunction`, `ConjunctionExplorer`, lapse, visibility, sky maps, and ground-truth validation cases (planet–star pairs, planetary trios, mixed groupings, retrograde triple crossings).

## Version 0.40.0 (breaktrhough)

- **Pinned astronomy stack** — `ephem==4.2.1`, `pymeeus==0.5.12`, `pyplanets==0.4.2` (reproducible ephemerides; stored in `requirements-astronomy.txt`).
- **Organic regression snapshots** — `montu/tests/test-planetary-ephemeris-organic.csv` and `montu/tests/test-stellar-positions-organic.csv` regenerated with this stack (`make organic-snapshots`).
- **`montu[desktop]` optional extra** — `pip install montu[desktop]` installs PySide6, Pygments, and plotly for MontuPython Desktop and `imontu --gui` / `imontu --sothic` (see `requirements-desktop.txt`).
- **Single `locations.json`** — removed the duplicate under `montu_gui/assets/`; MontuPython Desktop and the library both read `montu/data/locations.json`.
- **Catalogue metadata on `Observer(site=…)`** — every field from the JSON entry (`region`, `era`, `description`, `pressure_mbar`, `temperature_c`, …) is exposed as an instance attribute.
- **Climatology in `locations.json`** — mean surface pressure (from `alt_m`) and mean annual 2 m air temperature (NASA POWER 1991–2020) for all 72 bundled sites; used as default `pressure` and `temperature` unless overridden.
- **Explicit overrides with `site=…`** — constructor arguments (`lon`, `lat`, `height`, `pressure`, `temperature`) take precedence over catalogue values when passed explicitly (e.g. `Observer(site='thebes', pressure=0, temperature=0)` for a geometric horizon).
- **`Observer.sidereal_time(mtime)`** — local apparent sidereal time at the site (PyEphem convention, consistent with `Stars.where_in_sky`).
- **`Observer.__repr__` / `Observer.__str__`** — clearer site label, coordinates (6 decimals), elevation, and atmospheric summary.
- **`Stars.return_as='Star'`** — filter the catalogue and get a single `Star` or a name-keyed `dict` of `Star` objects instead of a `Stars` subset (`Stars(subset='bright', ProperName='Spica', return_as='Star')`).
- **`Stars.conditions_in_sky(at, observer)`** — batch rise/set/transit and related columns for every row in a subset (PyEphem per star; same engine as `Star.conditions_in_sky`); returns a DataFrame or updates in place with `inplace=True`.
- **`Stars.__new__` fix** — constructing `Stars(subset=…, **filters)` no longer loads the stellar CSV twice.
- **`Star.where_in_sky`** — adds `RAJ2000t` and `DecJ2000t` (J2000 coordinates with proper motion at epoch), included in `show_position()`.
- **`show_position()` / `show_conditions()`** — on `Sebau` and subclasses (`Star`, `Sun`, `Moon`, `Planet`); print epoch, site, and labelled quantities with units; return `None` (stdout only). Call `where_in_sky` or `conditions_in_sky` first.
- **`Planet('Sun')` / `Planet('Moon')`** — homogeneous factory call; returns a `Sun` or `Moon` instance (same construction as `Sun()` / `Moon()`). Numeric ids `10` and `301` also dispatch. Other names still return `Planet`.
- **`montu/tests/test_time_arithmetic_snippets.py`** — arithmetic on `Time` (`diff`, `add`, `subs`, `+`/`-`) from the Code Snippets notebook, including Gregorian reform and Sothic cycle cases.
- **`montu/tests/test_sky_show.py`**, **`montu/tests/test_stars_subset.py`**, **`montu/tests/test_observer.py`** — coverage for sky reporting, `return_as`, batch `conditions_in_sky`, site overrides, sidereal time, **`Star.show_properties()`**, and compact **`Star.__repr__`**.
- **`Star.show_properties()`** — print catalogue identifiers, coordinates, photometry, spectral type, and related fields for a single `Star`; return `None` (stdout only).
- **`Star.__repr__`** — compact one-line summary (name, constellation, J2000 coordinates, V magnitude, distance, spectral type); replaces the verbose inherited `Sebau` representation.
- **`HeliacalRise.print_rises()`** — display-only; returns `None` (no DataFrame chaining).
- **`SolarEclipse.path_map`** — Xavier Jubier Google-Map URL for the central-eclipse path, set when the catalogue row is loaded.
- **`EclipseConditions`** and **`show_details()`** — typed return from `conditions_eclipse()`; includes **`cond_map`** (Jubier URL for local circumstances at the observer site) and a formatted local report (kind, magnitude, contacts C1–C4, umbra duration).
- **`SolarEclipse.__str__`** — compact catalogue summary reorganized into Catalogue, Greatest eclipse, and Central path blocks (field names and units only; `path_map` at the end).
- **`montu/tests/test_phenomena.py`** — coverage for compact `SolarEclipse.__str__`, `path_map`, `cond_map`, `EclipseConditions.show_details()`, and display-only `print_rises`.
- **`examples/MontuPython-CodeSnippets.ipynb`** — promoted from `examples/inprogress/`; titled runnable snippets with explanations across utilities, time, observers, stars, planets, sky maps, heliacal rises, and solar eclipses (Hyades cluster plot, package constants, calendar/JD conversions, batch precession, `return_as='Star'`, geometric vs refracting observer).
- **`README.ipynb`** — Quickstart in-depth examples: Sirius heliacal rise at the 139 CE apokatastasis (Schaefer 1987, Thebes) and Thales' eclipse (585 BCE) with Xavier Jubier cross-checks; Peret Sopedet illustration under What's new.
- **`gallery/peret-sopedet-illustration.webp`** — heliacal rise of Sirius (Peret Sopedet) artwork for README and docs.
- **MontuPython Desktop tests** — location module expectations aligned with current Thebes climatology in `locations.json` (catalogue-driven assertions instead of stale hard-coded values).

## Version 0.31.0 (improve)

- **Stellarium validation — planetary ephemerides** — new parametrized test suite `montu/tests/test_planetary_ephemeris.py` against reference rows in `montu/tests/test-planetary-ephemeris.csv` (geocentric J2000 RA/Dec, phase, distance, solar elongation, and lunar elongation exported from Stellarium for Mercury through Saturn).
- **Stellarium validation — stellar precession and horizon events** — new test suite `montu/tests/test_stellar_positions.py` against `montu/tests/test-stellar-positions.csv` (Spica at Thebes-like coordinates: rise, transit, set, and transit altitude from 1600 to 2000; sampled on the 1st day of January, April, July, and October in the first year of each century).
- **Geometric vs refracting horizon** — stellar reference rise/set times use Stellarium’s geometric horizon (`pressure=0`, `temperature=0`); default `Observer` settings still apply standard atmospheric refraction. The feature tour notebook documents both conventions.
- **Example notebook** — `examples/MontuPython-FeatureTour.ipynb` (end-to-end package tour); `examples/docignore` lets `make docs-prepare` skip selected notebooks.
- **Tests and packaging** — pytest markers registered in `montu/tests/conftest.py` for packaged runs (`imontu --tests`); fixes for NumPy `timedelta64` deprecation in `Util.dt2cal()` and a spurious `RuntimeWarning` in heliacal-rise limiting-magnitude calculations.
- **`D2H` → `D2S`** — renamed the top-level decimal-to-sexagesimal alias (`D2S = Util.dec2sex`) across the package, documentation, notebooks, MontuPython Desktop, and MontuApp.
- **`imontu --sothic`** — open the interactive Sothic year calendar from the command line (same widget as MontuPython Desktop). A full civil date highlights that day (e.g. `"[hrw 0] I akhet 1"`); year-only forms open the Horus year without a highlighted cell (`"hrw 0"`, `"bce 1341"`, `"-1341"`). Requires PySide6 (installed via `make install-dev` or MontuPython Desktop dependencies).

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
- New top-level alias **`S2D`** for sexagesimal-to-decimal conversion; **`D2H`** (later renamed **`D2S`** in 0.30.1) now points to `dec2sex`.
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
