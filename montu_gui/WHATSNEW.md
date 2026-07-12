# What's New — MontuPython Desktop

Release notes for the graphical front-end to MontuPython.

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
