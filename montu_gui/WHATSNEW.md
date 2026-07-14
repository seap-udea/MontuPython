# What's New — MontuPython Desktop

Release notes for the graphical front-end to MontuPython.

## Version 0.2.5 (Desktop)

- **Bilingual localization (ES/EN)** — expanded and normalized translations across labels, forms, results text, and module UI.
- **Spanish terminology fixes** — updated wording for astronomical seasons and heliacal rises; standardized visible text from "ano" to "año".
- **Spanish help content restored** — rebuilt `help_es.json` with a module/key structure aligned to `help.json` for consistent fallback behavior.

## Version 0.2.4

- **Observer sites from the library** — predefined locations now load from `montu/data/locations.json`; the Desktop uses the same catalogue as `montu.Observer.list()` / `montu.Observer(site=…)`.
- **Module order** — sidebar and Home list reordered: Calendar, Observer Location, Sky Map, Seasons & Lunar Phases, Planetary Ephemerides, Orientation Disk, Star Alignments, Heliacal Rises (Home remains first in the sidebar).
- **Orientation Disk (⭕)** — search window per body is **max(2 years, orbital period)** (Sun fixed at 2 years; Jupiter and Saturn use their sidereal periods); help text and reference-year label updated.
- **Sky Map (🌌) defaults** — date **June** 2500 BCE; **Horizon** and **Meridian view** on by default; bodies **Sun, Moon, Mercury, and Venus** selected (`default.json` updated).
- **Seasons & Lunar Phases (🎑)** — below **Caniucular** on each season card: Sun rise/set **azimuth** and **local time** for the observing site on that season’s calendar day; below the year selector: **Observing site** line from the global location; module listens to **LocationState** and recalculates on site change.
- **Seasons fix** — solstice sunrise azimuths now use that day’s sunrise/sunset (local noon anchor) and display degrees correctly (no erroneous `RAD` conversion).
- **Heliacal Rises (🌅)** — module description extended with a **this link** to open a non-modal **Historical Sirius heliacal rises** window (`montu/data/historical-heliacal-rises.json`); application stays usable while the list is open.

## Version 0.2.3

- **Heliacal Rises (🌅)** — new module to find first visible morning appearances of named bright stars and classical planets using the library’s `HeliacalRise` models (Schaefer 1985/1987, Belokrylov 2011).
- **Inputs** — shared observer location; star/planet selector; initial date with BCE/CE and year/month/day spinboxes; mixed or proleptic calendar; year range (`StepSpinBox`); model-specific parameters (extinction, limiting magnitude, solar depression, scan step).
- **Results** — summary line (interval, model source, elapsed time) and table with mixed, proleptic, and caniucular dates, local time, body/Sun altitude and azimuth, and days since the prior event.
- **Peret Sopedet illustration** — full-width image below the table with English title *Illustration of Egyptian Peret Sopedet on the Giza plateau*; contextual help on the heliacal rise of Sirius (Sopedet) and image credits (Nano Banana / Preview, based on Bob Moler’s Stellarium–GIMP simulation).
- **Defaults** — factory observer preset **Thebes (Luxor)**; Sirius, BCE 2782-06-01, 10-year range, Schaefer 1987 (`k=0.25`, `limiting_mag_zenith=6.5`, `sun_depression=−7.8`).
- **Let's Python!** — example script `montu_gui/pages/examples/heliacal_rise.py` reproduces the Desktop defaults.
- **Observer location fix** — Star Alignments no longer overwrites the global observer to Giza on startup; saved observer config is reapplied after all pages are built, and Observer Location refreshes from shared state when opened.
- **Help** — `help.json` entries for all Heliacal Rises inputs, results, and Peret Sopedet; module listed on Home.

## Version 0.2.2

- **Frozen-bundle assets** — new `montu_gui/utils/bundle_paths.py` resolves `home.json`, `help.json`, logos, and other files from the PyInstaller `_MEIPASS` tree; all `montu_gui/assets/` files are copied explicitly in `montu-desktop.spec`. Fixes empty Home page and “No help text found” dialogs on Windows installs.
- **Windows icon** — `montu-python-logo-complete.ico` used for the `.exe` and taskbar; sidebar logo path fixed in frozen builds (`main.py` → `gui_asset()`).
- **User config** — `default.json` bundled and writable `config.json` path under `%AppData%` on Windows (same pattern as macOS Application Support).
- **Module branding** — compact **Module brand** strip (sidebar emoji, title, one-line description from `home.json`) in each module’s input panel; top title bars removed to recover vertical space.
- **Plotly embed** — maps and charts fill the available panel height in tall windows and gain scroll bars when the view is shorter than the minimum map height (`plotly_html.py`; Sky Map minimum 560 px).
- **Orientation Disk example** — `montu_gui/pages/examples/orientation_disk.py` redrawn to match the Desktop azimuth disk (four rise/set extremes per body, triangle markers, compass labels).

## Version 0.2.0

- **Sky Map (🌌)** — new module with separate azimuthal maps for the northern and southern celestial hemispheres.
- Precessed stars to the selected date; limiting magnitude filter; local solar time and observer site from **Observer Location**.
- **Bodies on map** — Sun, Moon, and classical planets as emoji symbols; **Lines on map** — ecliptic (precessed) and horizon (elevation 0°).
- **Horizon** mode — azimuth marks every 30° on the horizon circle; stars below the horizon tinted green; tooltips show RA, Dec, azimuth, and elevation (four decimal places).
- **Meridian view** — rotate the map so the local meridian (LST) is at the top; RA grid labels follow absolute meridians.
- Map title includes **LST**; legend lists overlay lines only (bodies remain on the map without legend entries).
- **Constellation set** — choose **IAU Constellations**, **Ancient Egyptian**, or **Dendera Egyptian** asterisms (Stellarium sky-culture data in `montu/data`; help cites sources and licences).
- **Persistent configuration** — module parameters and observer location saved to `montu_gui/user/config.json`; factory defaults in `montu_gui/user/default.json`.
- **Home** — redesigned welcome page with module list (icon + description) and **Save configuration** / **Reset configuration** buttons.
- **Contextual help** expanded across modules; shared field text in `help.json` under `_common` (resolved via `$ref`); Sky Map and Orientation Disk fields documented; **Map** help explains window sizing (no scroll bars, no zoom in this projection).
- All main modules export and restore their parameters through `export_config` / `apply_config` (Calendar, Location, Planets, Seasons, Star Alignments, Orientation Disk, Sky Map).
- **PlotlyView** — each embedded chart uses its own temp HTML file so multiple plots on one page do not overwrite each other.

## Version 0.1.5

- **Historical dates from the library** — Calendar Calculator and examples now use `montu.load_historical_dates()` and the catalogue bundled in `montu/data/historical_dates.json` (removed duplicate data from the Desktop PyInstaller spec).
- **Calendar example** — added the Canopus Decree (`bce 238-03-07`) to `montu_gui/pages/examples/calendar_conversion.py`.
- **Aligned with MontuPython 0.20.5** — Desktop works with the SPICE-free library (no kernel loading on import).
- **Release tooling** — `release-pipeline.sh` supports Desktop releases via `--tag` (`desktop-release`, git tag, and `desktop-ci`).

## Version 0.1.2

- **Desktop distribution** — PyInstaller packaging (`montu_gui/montu-desktop.spec`), local build scripts (`bin/build-desktop.sh`, `bin/build-desktop.ps1`), and GitHub Actions CI for macOS and Windows releases (tag `desktop-v*`).
- **README** — download section for Desktop builds from [GitHub Releases](https://github.com/seap-udea/MontuPython/releases).
- **DEVELOPER.md** — `make` targets for `desktop-build`, `desktop-package`, `desktop-release`, and `desktop-ci`.
- **Let's Python!** — new **Copy and Test in Colab** button copies the example script and opens `examples/MontuPython-TestCode.ipynb`.
- Let's Python! button moved to the left sidebar column (compact layout) on Calendar, Location, Planets, Seasons, Star Alignments, and Orientation Disk pages.
- Colab-friendly example scripts for all six modules in `montu_gui/pages/examples/`; bundled in the PyInstaller build.
- Examples updated for MontuPython **0.20.1**: `dec2sex` / `sex2dec` naming and explicit `mercator_sky_map` import in Star Alignments.

## Version 0.1.1

- **Orientation Disk (⭕)** — new module for the northernmost and southernmost rise and set azimuths of celestial bodies on a polar disk.
- Reference year with BCE / CE selector; computation sweeps a **3-year window** from 1 January (long enough to capture Venus’s synodic cycle).
- Celestial bodies: Sun (added by default), Moon, Mercury, Venus, Mars, Jupiter, Saturn, and named bright stars from the Montu stellar catalogue.
- Star magnitude filter in the add-body dropdown (default limiting magnitude **V = 1.0**).
- Per-body horizon altitude (`h`, default 0°) for rise/set visibility (e.g. after twilight).
- Add and remove bodies from the disk; duplicate bodies are prevented; at least one body must remain.
- Polar azimuth disk (Plotly): N at top, graduated azimuth scale, coloured arrows per body — △ rise (east), ▽ set (west), N/S northernmost / southernmost extremes.
- Dark colour palette (blues, greens, reds) for body arrows and list colour swatches.
- Celestial bodies panel expands vertically with the sidebar; the body list grows as entries are added and scrolls when needed.
- Disk plot fills the available window and resizes with the main window.
- **Let's Python!** example script: `montu_gui/pages/examples/orientation_disk.py`.
- Pure-logic backend: `montu_gui/modules/orientation_disk.py` (no Qt dependency).
- **Star Alignments** — sky map expands horizontally when the window is maximised or the splitter is moved; removed the fixed-width scroll-area wrapper.
- **Plotly embed** (`plotly_html.py`, `plotly_view.py`) — responsive width, `autosize`, and automatic `Plotly.Plots.resize()` on load and widget resize; benefits Star Alignments, Planetary Ephemerides, and Orientation Disk.

## Version 0.1.0

- First packaged release of **MontuPython Desktop**, the PySide6 GUI for MontuPython.
- **Home** page with project overview, credits, and quick links to the repository and documentation.
- **Observer Location** — set observer latitude, longitude, and altitude; map view and named presets.
- **Calendar Calculator** — convert among Julian/Gregorian, proleptic Gregorian, caniucular (Egyptian civil), and Julian Day formats; historical date presets; runnable Python example export.
- **Seasons & Lunar Phases** — solstices, equinoxes, and lunar phase tables for a chosen year and site.
- **Planetary Ephemerides** — time-series charts of ephemeris properties for selected planets.
- **Star Alignments** — find stars passing through a fixed azimuth/elevation direction over a date range, with Mercator sky map.
- Contextual **help** tooltips on form fields (editable via `montu_gui/assets/help.json`).
- Unified visual **theme** and collapsible sidebar navigation.
- **Let's Python!** dialogs to view, copy, and download example scripts for each module.
- Optional `--debug` mode logging conversions and UI events to the terminal.
