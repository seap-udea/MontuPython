# Montu Python
## Astronomical ephemerides for the ancient world

<p align="center"><img src="https://github.com/seap-udea/MontuPython/raw/main/montu/data/montu-python-logo-complete.webp" width="300" alt="MontuPython logo"/></p>

[![PyPI](https://img.shields.io/pypi/v/montu?color=blue)](https://pypi.org/project/montu/) [![Downloads](https://img.shields.io/pypi/dw/montu)](https://pypi.org/project/montu/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Python](https://img.shields.io/badge/python-3.8%20|%203.9%20|%203.10%20|%203.11-blue)](https://www.python.org/) [![GitHub](https://img.shields.io/badge/GitHub-seap--udea%2FMontuPython-blue?logo=github)](https://github.com/seap-udea/MontuPython) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21298063.svg)](https://doi.org/10.5281/zenodo.21298063) [![arXiv](https://img.shields.io/badge/arXiv-preprint%20forthcoming-b31b1b)](#) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/README.ipynb)
[![Powered by PyEphem](https://img.shields.io/badge/powered%20by-PyEphem-blue)](https://pypi.org/project/ephem/) [![Powered by PyMeeus](https://img.shields.io/badge/powered%20by-PyMeeus-blue)](https://pypi.org/project/pymeeus/) [![Powered by spiceypy](https://img.shields.io/badge/powered%20by-spiceypy-blue)](https://pypi.org/project/spiceypy/) [![Powered by pyplanets](https://img.shields.io/badge/powered%20by-pyplanets-blue)](https://pypi.org/project/pyplanets/)

`MontuPython` is a Python package intended to compute astronomical ephemerides in the ancient world, thousands of years before present. It was initially designed to compute ephemerides for the ancient Egypt, but it can also be used to study astronomical phenomena in other sites of interest for cultural astronomy (archeoastronomy).

The package was originally developed to streamline calculations on archeoastronomy that were previously performed manually using astronomical software such as Stellarium. While there are already outstanding libraries for positional astronomy—developed by a vibrant community of professionals and enthusiasts—such as [PyEphem](https://pypi.org/project/ephem/), [pyplanets](https://pypi.org/project/pyplanets/), [PyMeeus](https://pypi.org/project/pymeeus/), [spiceypy](https://pypi.org/project/spiceypy/), and others (many of which make MontuPython possible)—our goal here is different. We aim to provide a tool that not only simplifies interaction with these libraries, but, more importantly, is specialized for the study of astronomical phenomena that occurred thousands of years ago.

While MontuPython places a special emphasis on ancient Egypt—understandably, since Egyptian religion, culture, and statehood were profoundly intertwined with astronomy—the tool is broadly applicable to many other contexts. In fact, MontuPython can be used to study archaeoastronomical phenomena in any culture around the world. Its features are designed to support research across a diverse range of civilizations, making it a versatile resource for all archaeologists, historians, and cultural astronomers interested in humanity's relationship with the sky.

MontuPython has benefited from the contributions and advice of several collaborators:

- Francisco "Tito" Vivas provided an actual archaeoastronomical fase of study who originally motivated the development of this software.
- Prof. José Lull of the Universitá Autónoma de Barcelona (UAB), provided scientific, archaeoastronomical, and Egyptological advice.
- Juanita Agudelo developed the first versions of the MontuPython web app.
- Luis Arroyo developed the version of the web app which precedes the desktop app.

## Useful Resources 

For more information before you begin, please refer to:

- **MontuPython documentation**: The full documentation of the package is available at [https://montupython.readthedocs.io](https://montupython.readthedocs.io/en/latest/).

- **MontuPython Desktop**: This is the graphical front-end to the library, built with PySide6. It offers pre-built tools for users who prefer not to program at the Python level, or who want quick access to core functionality — including a **Calendar Calculator** and **Seasons & Lunar Phases** module.

  Pre-built installers for **macOS** and **Windows** are published on **[GitHub Releases](https://github.com/seap-udea/MontuPython/releases)**. Open that page and look for a release whose tag starts with `desktop-v` (for example, `desktop-v0.1.1`). Download the file that matches your operating system:

    | Platform | File to download | How to run |
    |----------|------------------|------------|
    | **[Descarga versión de macOS](https://jorgezuluaga.github.io/sh/montupython-desktop-macos)** | `*-macos.dmg` or `*-macos.zip` | Open the `.dmg` and drag **MontuPython Desktop** to Applications, or unzip the `.zip` and open `MontuPython Desktop.app`. If macOS shows a security warning the first time, right-click the app and choose **Open**. |
    | **[Descarga versión de Windows](https://jorgezuluaga.github.io/sh/montupython-desktop-windows)** | `*-windows.zip` | Unzip the archive, open the `MontuPython-Desktop` folder, and double-click `MontuPython-Desktop.exe`. If SmartScreen warns you, choose **More info** → **Run anyway**. |

    No Python installation is required. Most features work offline; the observer map (OpenStreetMap) needs an internet connection.
    
    If you already use the `montu` library, install the optional **Desktop** extra and launch the GUI from a terminal:

    ```bash
    pip install montu[desktop]
    montu-gui
    # or
    imontu --gui
    ```

    `montu[desktop]` adds **PySide6**, **Pygments**, and **plotly**. On first launch, `montu-gui` / `imontu --gui` also fetches the `montu_gui` sources if they are not already present in your checkout. The same extra is required for `imontu --sothic`.

## Download and install

### Library only

Install the core package from PyPI:

```bash
pip install montu
```

This pulls in the pinned astronomy stack (`ephem`, `pymeeus`, `pyplanets`), Dash (MontuApp), and the usual scientific Python dependencies. It does **not** install PySide6 or other MontuPython Desktop GUI packages.

### Library + Desktop GUI

Add the optional **Desktop** extra when you want the PySide6 GUI or `imontu --sothic`:

```bash
pip install montu[desktop]
```

From a development checkout:

```bash
pip install -e ".[desktop]"
```

## Usage in Google Colab

Uncomment and run in the next cell (replace the branch name if needed), then **restart the runtime** before importing `montu`:


```python
# Uncomment if you are in Google Colab
# !pip install -q "git+https://github.com/seap-udea/MontuPython@maat"
# For Desktop GUI in a local Python environment (not typical on Colab):
# !pip install -q "montu[desktop] @ git+https://github.com/seap-udea/MontuPython@maat"
```

## Tutorial by example notebooks

The following notebooks introduce MontuPython step by step. Open any of them in Google Colab and start working immediately — they are ideal for training new users or for testing code generated by MontuPython Desktop.


- **Code snippets** Runnable examples — one MontuPython feature per snippet, each with a short explanation. Covers utilities, time and calendars, observers, stars, Sun/Moon/planets, seasons, sky maps, heliacal rises, solar eclipses, conjunctions, and planetary stations. A good first overview before the topic-specific notebooks below.

   **File**: `MontuPython-CodeSnippets.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-CodeSnippets.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-CodeSnippets.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-CodeSnippets.ipynb)

- **Time: Basic time functionalities** Illustrates the `montu.Time` class across its full range of calendars and time scales. Parse dates in ISO, SPICE, Julian day, and terrestrial time; convert between proleptic Gregorian, mixed Julian/Gregorian, ancient Egyptian civil, and sothic (Horus year) calendars; and add or subtract calendar periods.

   **File**: `MontuPython-MonTime.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-MonTime.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-MonTime.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-MonTime.ipynb)

- **Egyptian civil calendar: historical cross-checks** Compares MontuPython's computed Egyptian civil dates against the civil dates recorded in Egyptology sources for key historical events (apokatastasis anchors, Ptolemaic decrees, dated papyri). Uses the same event catalogue as MontuPython Desktop's Calendar Calculator.

   **File**: `MontuPython-EgyptianCalendar.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-EgyptianCalendar.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-EgyptianCalendar.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-EgyptianCalendar.ipynb)

## Advanced examples


The following notebooks present archaeoastronomical computations of broader historical and interpretive interest.


- **Heliacal rises: Sirius step by step** Reconstructs the first morning visibility of Sirius from Thebes around 2782 BCE — the era of the first *apokatastasis* of the Sothic cycle. Works through all four heliacal-rise models in MontuPython (Ptolemy, Schaefer 1985/1987, Belokrylov 2011) and explains each visibility criterion rather than treating `HeliacalRise` as a black box.

   **File**: `MontuPython-HeliacalRises.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-HeliacalRises.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-HeliacalRises.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-HeliacalRises.ipynb)

- **Solar eclipses** Predicts **local circumstances** of solar eclipses from the bundled NASA Five Millennium catalogue (Espenak & Meeus): search the catalogue, validate the 2024 Dallas totality, and survey historical eclipses visible from Thebes.

   **File**: `MontuPython-SolarEclipses.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-SolarEclipses.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-SolarEclipses.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-SolarEclipses.ipynb)

- **Conjunctions** Finds angular conjunctions of planets and stars with `Conjunction` and `ConjunctionExplorer`: separation and visibility at one epoch, lapse intervals, Plotly sky maps, and reference cases from simple pairs to planetary trios and retrograde triple crossings.

   **File**: `MontuPython-Conjunctions.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-Conjunctions.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-Conjunctions.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-Conjunctions.ipynb)

- **Montunctions** Computes *marstices* — the stationary points of Mars — and their conjunctions with Aldebaran near the winter solstice throughout ancient Egyptian history. Includes periodicity analysis and publication-ready plots. The methodology and figures supported Francisco Vivas Fernández's Ph.D. thesis on the astronomical and landscape orientations of Senenmut's monuments in the Theban necropolis (2023).

   **File**: `MontuPython-Montunctions.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-Montunctions.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-Montunctions.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-Montunctions.ipynb)

- **Example: Venus Azimutal Distributions** Computes the distribution of Venus rising and setting azimuths at dawn and dusk over millennia, as observed from Thebes. Distinguishes morning-star and evening-star apparitions and summarizes their preferred horizon directions — a classic archaeoastronomical exercise for Venus cult and calendar studies.

   **File**: `MontuPython-VenusAzimuths.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-VenusAzimuths.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-VenusAzimuths.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-VenusAzimuths.ipynb)

- **Observer Horizons** Computes and plots the real topographical horizon profile for an observer by downloading and processing Copernicus GLO-30 DEM tiles. Demonstrates interactive Plotly and Leaflet maps, celestial object occultation, and sunrise azimuth calculations.

   **File**: `MontuPython-ObserverHorizon.ipynb`. **Links**: [GitHub](https://github.com/seap-udea/MontuPython/blob/main/examples/MontuPython-ObserverHorizon.ipynb) | [ReadTheDocs](https://montupython.readthedocs.io/en/latest/examples/MontuPython-ObserverHorizon.html) | [Colab](https://colab.research.google.com/github/seap-udea/MontuPython/blob/main/examples/MontuPython-ObserverHorizon.ipynb)


## What's new


Release notes are now maintained in [WHATSNEW.md](WHATSNEW.md).

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/gallery/peret-sopedet-illustration.webp?raw=true" alt="Illustration of the heliacal rise of Sirius (Peret Sopedet) on the Giza plateau" width="800"/></p>

<p align="center"><em>Illustration of Egyptian Peret Sopedet on the Giza plateau</em> — the heliacal rise of Sirius (<i>Sopedet</i>) that marks the civil New Year (<code>[hrw 0] I <i>akhet</i> 1</code>). Image generated with Nano Banana, edited in Preview, based on Bob Moler’s Stellarium–GIMP simulation (<a href="https://bobmoler.blog/2017/02/07/02072017-ephemeris-sirius-an-important-star-in-history/">Bob Moler's Ephemeris Blog</a>, 7 February 2017).</p>


## AI assistance disclosure

Starting with version *0.20.0*, portions of the code, inline documentation, code review, and debugging in this repository were assisted by AI language models, including **Cursor Composer 2.5**, **Anthropic Claude Sonnet 4.6**.

The human authors maintain that all scientific ideas, the overall project conception, the package and notebook architecture, the design of the numerical and scientific workflows, and their interpretation remain original contributions of the human authors. AI tools were used exclusively as coding and writing assistants — analogous to a spell-checker or a compiler — and bear no authorship over the scientific content of this project. AI models also assisted with translating text from Spanish (the native language of the human authors) into English, and with English spelling and grammar review.

## Quickstart

In the rest of this document we show you basic applications of the package. Since how it is imported until how it can be used for realistic calculations.

This README is at the same time an executable notebook. A runable version of the README can be find in the [Google Drive Repo of the package](https://drive.google.com/drive/folders/11L59yZ3A1g1ZT7v_dLDPwLMnRMR-tFgE?usp=sharing) as the file [`README.ipynb`]() 

### Preparation

> **NOTE:** This notebook is autogenerated with `Jupyter`. Some of the commands you will find here are only intended to make the README executable. In no sense are mandatory when using `MontuPython`.

You may import the package using:


```python
%matplotlib inline
import matplotlib.pyplot as plt
!mkdir -p gallery

import montu
```

    MontuPython version 0.43.3. 𓇍𓇋𓇋𓏏𓅓𓊵 𓎛𓎡𓄿𓀭𓎛𓈖𓂝𓎡 (ii-ti m Htp, HkAx Hn'-k)


It is important that before using the most interesting commands of the package, load relevant data:


```python
# Useful aliases
from montu import D2S, PRINTDF, TABLEDF
# Load stars
allstars = montu.Stars()
```

    Loading stellar catalogue montu_stellar_catalogue_v38.csv


Two very basic codes that do something non-trivial albeit simple calculations in `MontuPython are:`

1. Compute the position of mars at a given time and while observing from a given site on Earth:


```python
mtime = montu.Time('-2500-01-01 12:00:00.00')
tebas = montu.Observer(lon=33,lat=24)
mars = montu.Planet('Mars')
mars.conditions_in_sky(at=mtime,observer=tebas)
mars.show_conditions()
```

    Mars — sky conditions
      Epoch: -2500-01-22 12:00:00 / -2500-01-01 12:00:00.000000  (JED 807954.000000)
      Site: Custom site — lat 24.000000°, lon 33.000000°, 0 m  (P=1013.25 mbar, T=15.0 °C)
      Name: Mars
      Hour angle: 12:20:15.604 h
      Visual magnitude: -1.13 mag
      Rise time (UTC): -2500-01-22 16:49:49
      Rise azimuth: 63.133665°
      Set time (UTC): -2500-01-23 06:24:24
      Set azimuth: 296.924743°
      Transit time (UTC): -2500-01-22 23:37:37
      Transit elevation: 89.841412°
      Elongation from Sun: -157.818359°
      Distance from Earth: 0.660488 AU
      Distance from Sun: 1.626107 AU
      Angular diameter: 14.171"
      Illuminated fraction: 98.65 %
      Heliocentric latitude: 1.888422°
      Heliocentric longitude: 111.258645°
      Heliocentric longitude (alt.): 111.258645°
      Circumpolar: no
      Never rises: no


2. Obtain the information about a star from the stellar catalogue and, as in the case of the planet, obtain the position of the star in the sky.


```python
mtime = montu.Time('-2500-01-01 12:00:00.00')
tebas = montu.Observer(lon=33,lat=24,height=0)
aldebaran = allstars.get_stars(ProperName='Aldebaran', return_as='Star')
aldebaran.where_in_sky(at=mtime,observer=tebas)
aldebaran.show_position()
```

    Aldebaran — sky position
      Epoch: -2500-01-22 12:00:00 / -2500-01-01 12:00:00.000000  (JED 807954.000000)
      Site: Custom site — lat 24.000000°, lon 33.000000°, 0 m  (P=1013.25 mbar, T=15.0 °C)
      Name: Aldebaran
      RA (J2000): 04:35:35.599 h
      Dec (J2000): 16.745934°
      RA (J2000, proper motion): 04:35:36.404 h
      Dec (J2000, proper motion): 16.745996°
      RA (epoch): 00:36:32.065 h
      Dec (epoch): -2.261722°
      RA (geocentric): 00:36:32.065 h
      Dec (geocentric): -2.261722°
      Azimuth: 107.617963°
      Elevation: 29.633599°


### Working with time

One of the most interesting and basic functionalities of MontuPython is to convert date among 
different type of calendars and astronomical scales.  You may taste these functionalities using:


```python
mtime = montu.Time('bce2501-01-01 12:00:00')
mtime
```




    Time('-2500-01-01 12:00:00.000000'/'-2500-01-22 12:00:00'/'[hrw 280] I shemu 17'/JED 807954.0/JTD 807954.6909688)



other alternative formats for the same date are:


```python
mtime = montu.Time('2501 b.c.e. 01-01 12:00:00')
mtime = montu.Time('-2500-01-01 12:00:00')
```

If you print this time object you will get:


```python
print(mtime)
```

    2501 B.C. 01-01 12:00:00.00 / [hrw 280] I shemu 17:
        Date in ISO format: 2501 B.C. 01-01 12:00:00.00
        Date in proleptic UTC: -2500-01-01 12:00:00.000000
        Date in mixed UTC: -2500-01-22 12:00:00
        Weekday: 2 (monday)
        Date in sothic format: [hrw 280] I shemu 17
        Terrestrial time: tt [seconds]: -142006202700.3
        UTC time: jed [days]: 807954.0
        Delta-t = TT - UTC [seconds]: 59699.7
    


For more details:


```python
mtime.details()
```

    Montu Time Object:
    -------------------------- 
    Readable:
        Date in proleptic UTC (.readable.datepro): -2500-01-01 12:00:00.000000
        Date in mixed UTC (.readable.datemix): -2500-01-22 12:00:00
        Date in SPICE format (.readable.datespice): 2501 B.C. 01-01 12:00:00.00
        Date in sothic format (.readable.datesot): [hrw 280] I shemu 17
        Weekday (.readable.weekday): 2 (monday)
        Components (.readable.comps): [-1, 2500, 1, 1, 12, 0, 0, 0]
    Objects:
        Date in datetime64 format (.readable.obj_datetime64): -2500-01-01T12:00:00.000000
        Date in PyPlanet Epoch (.obj_pyplanet): 807954.0
        Date in PyEphem Epoch (.obj_pyephem): -2501/1/22 12:00:00
    General:
        Is bce (.bce): True
        Is Julian (.isjulian): True
    Uniform scales:
        Terrestrial time:
            tt (.tt): -142006202700.3
            jtd (.jtd): 807954.6909688
            htd (.htd): 102457.19096879999
        UTC time:
            et (.et): -142006262400.0
            jed (.jed): 807954.0
            hed (.hed): 102456.5
        Delta-t = TT - UTC (.deltat): 59699.7
    


Notice that the date in Gregorian proleptic will be **bce 2501-01-01** but in the mixed calendar that uses Julian calendar before its adoption at 1582-10-04, will be **bce 2501-01-22**.

#### Egyptian civil calendar

In the previous output you may also notice that the class `montu.Time` automatically convert the gregorian/julian date into the *civil egyptian calendar* or the *Sothic* calendar. In the previous example, the date bce 2501-01-22 (julian) correspond to the civil egyptian date of **I *shemu* 17** (MontuPython: `[hrw 280] I shemu 17`). Since the years in the *Sothic calendar* were regularly referred to specific king's reigns, we have introduced in `MontuPython` the so-called *Horus years* (abreviated *hrw*). The zero *Horus year* is bce 2782 which corresponds to the first *apokatastasis*, namely the year when the day **I *akhet* 1** coincide with the heliacal rise of Sirius. 

You may also provide a date in the *Sothic calendar* and obtain the corresponding julian date:


```python
mtime = montu.Time('[hrw 1461] I akhet 1',calendar='sothic')
print(mtime)
```

    1322 B.C. 07-20 08:45:50.800000 / [hrw 1461] I akhet 1:
        Date in ISO format: 1322 B.C. 07-20 08:45:50.800000
        Date in proleptic UTC: -1321-07-20 08:45:50.8
        Date in mixed UTC: -1321-07-20 00:00:00
        Weekday: 3 (tuesday)
        Date in sothic format: [hrw 1461] I akhet 1
        Terrestrial time: tt [seconds]: -104784376449.2
        UTC time: jed [days]: 1238762.5
        Delta-t = TT - UTC [seconds]: 31550.8
    


As you may see, this is the second *apokatastasis*, ie. in bce 1322-07-20, the heliacal rise of Sirius happened again at I *akhet* 1.

#### Operations with dates

You may add or substract time to a given date. This is done by adding or substracting seconds to the reference time:


```python
mtime = montu.Time('2001-01-01 12:00:00',format='iso')
(mtime,
 mtime - 12*montu.HOUR, 
 mtime + 1*montu.DAY, 
 mtime - 3*montu.CALYEAR, 
 mtime + 20*montu.JULYEAR)
```




    (Time('2001-01-01 12:00:00.000000'/'2001-01-01 12:00:00'/'[hrw 4784] I shemu 14'/JED 2451911.0/JTD 2451911.0007419),
     Time('2001-01-01 00:00:00.000000'/'2001-01-01 00:00:00'/'[hrw 4784] I shemu 14'/JED 2451910.5/JTD 2451910.5007419),
     Time('2001-01-02 12:00:00.000000'/'2001-01-02 12:00:00'/'[hrw 4784] I shemu 15'/JED 2451912.0/JTD 2451912.0007419),
     Time('1998-01-02 12:00:01.097278'/'1998-01-02 12:00:00'/'[hrw 4781] I shemu 14'/JED 2450816.0000127/JTD 2450816.0007419),
     Time('2021-01-01 11:59:51.904328'/'2021-01-01 11:59:59'/'[hrw 4804] I shemu 19'/JED 2459215.9999063/JTD 2459216.0007419))



As you may notice, adding or substracting and integer number of seconds not necesarily correspond to adding or sutracting days or years to the calendar. This is because of the difference in UT and TT: UT is always behind TT by a certain amount of seconds.  Normally leapseconds are included every once in a while. However to calculate ephemerides in the ancient world, `MontuPython` uses a continuous model of deltat that small discrepancies from year to year.

To add calendar periods of time you may use the `add` method of `Montu.Time`:


```python
mtime = montu.Time('2001-01-01 12:00:00')
mtime2 = mtime.add(years=1)
print(mtime2)
```

    2002-01-01 12:00:00.000000 / [hrw 4785] I shemu 14:
        Date in ISO format: 2002-01-01 12:00:00.000000
        Date in proleptic UTC: 2002-01-01 12:00:00.000000
        Date in mixed UTC: 2002-01-01 12:00:00
        Weekday: 3 (tuesday)
        Date in sothic format: [hrw 4785] I shemu 14
        Terrestrial time: tt [seconds]: 63158464.3
        UTC time: jed [days]: 2452276.0
        Delta-t = TT - UTC [seconds]: 64.3
    


### Working with stars

The stellar catalogue included with `MontuPython` contains more than 119 000 stars, including almost 9 000 visible to the naked eye.  Information about stars is stored in a `pandas Data Frame` whose columns are:


```python
allstars.data.columns
```




    Index(['MN', 'HD', 'HR', 'HIP', 'Gl', 'Name', 'OtherDesignations',
           'ProperName', 'Bayer', 'Flamsteed', 'Constellation', 'RAJ2000',
           'DecJ2000', 'GalLonJ2000', 'GalLatJ2000', 'pmRA', 'pmDec', 'RadVel',
           'Distance', 'Vmag', 'Vmag_min', 'Vmag_max', 'B-V', 'SpType',
           'Luminosity', 'XJ2000', 'YJ2000', 'ZJ2000', 'VXJ2000', 'VYJ2000',
           'VZJ2000', 'Primary', 'MultipleID', 'IsMultiple', 'IsVariable'],
          dtype='str')



Altough you may manipulate this DataFrame using the conventional commands in pandas, we have designed several useful methods to obtain subsets of the catalogue. For instance if you want to extract the stars visible to naked eye, the command would be:


```python
stars = allstars.get_stars(Vmag=[-2,6.5])
print(f"There is {stars.number} visible to the naked eye in the catalogue")
```

    There is 8920 visible to the naked eye in the catalogue


You can use any of the properties of the stars to filter them. A common filter is to look for single stars:


```python
aldebaran = stars.get_stars(ProperName='Aldebaran', return_as='Star')
aldebaran
```




    Star('Aldebaran'/'Tau'/04:35:55.237/16:30:33.484/V=0.87 mag/d=20.43 pc/Sp=K5III)



To get the whole properties of a Star:


```python
aldebaran.show_properties()
```

    Aldebaran — catalogue properties
      Proper name: Aldebaran
      Name: Aldebaran
      Bayer designation: α Tau
      Flamsteed number: 87 Tau
      Constellation: Tau
      HIP: 21421
      HD: 29139
      HR: 1457
      Gliese: Gl 171.1A
      Other designations: 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau
      RA (J2000): 04:35:55.237 h
      Dec (J2000): 16:30:33.484°
      Galactic longitude: 180:58:18.868°
      Galactic latitude: -20:14:53.892°
      Proper motion RA: 62.78 mas/yr
      Proper motion Dec: -189.36 mas/yr
      Radial velocity: 54.5 km/s
      Distance: 20.4332 pc
      Visual magnitude: 0.87 mag
      B−V colour index: 1.538
      Spectral type: K5III
      Luminosity: 163.23 L☉
      Multiple system: yes
      Variable star: yes


All information about a star is stored in the `data` attribute:


```python
aldebaran = stars.get_stars(ProperName='Aldebaran')
TABLEDF(aldebaran.data)
```

    |    |   MN |    HD |   HR |   HIP | Gl        | Name      | OtherDesignations                                                                     | ProperName   | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   GalLonJ2000 |   GalLatJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |   VXJ2000 |   VYJ2000 |   VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |
    |----|------|-------|------|-------|-----------|-----------|---------------------------------------------------------------------------------------|--------------|---------|-------------|-----------------|-----------|------------|---------------|---------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|-----------|-----------|-----------|-----------|--------------|--------------|--------------|
    | 14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran    | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |       180.972 |      -20.2483 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |       163.23 |  7.02722 |  18.2876 |  5.80666 | 1.528e-05 | 5.709e-05 | -2.14e-06 |     21368 | Gl 171.1     |            1 |            1 |


> **NOTE:** The `TABLEDF` is a special function of `MontuPython` that allows you to produce a readable table out of a pandas DataFrame, a list of dictionaries and other iterable python objects. See `tabulate` package for details.

Another useful method included with the class `Stars` is that of filtering the getting the stars close to a given point in the sky. For illustrare, below is the command to obtain all stars in the sky with magnitudes less than 5 and that are at 5.5 degrees or less than Aldebaran:


```python
hyades = stars.get_stars_around(center=[aldebaran.data.RAJ2000,aldebaran.data.DecJ2000],radius=5.5,Vmag=[-1,5])
hyades
```




    18 star(s):
    |      |   MN |    HD |   HR |   HIP | Gl        | Name           | OtherDesignations                                                                     | ProperName     | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   GalLonJ2000 |   GalLatJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |    VXJ2000 |    VYJ2000 |    VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |
    |------|------|-------|------|-------|-----------|----------------|---------------------------------------------------------------------------------------|----------------|---------|-------------|-----------------|-----------|------------|---------------|---------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|------------|------------|------------|-----------|--------------|--------------|--------------|
    |   14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran      | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran      | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |       180.972 |      -20.2483 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |    163.23    |  7.02722 |  18.2876 |  5.80666 |  1.528e-05 |  5.709e-05 | -2.14e-06  |     21368 | Gl 171.1     |            1 |            1 |
    |  259 |  260 | 28319 | 1412 | 20894 | nan       | Chamukuy       | 78 Tau/78The2Tau/Chamukuy/HD 28319/HIP 20894/HR 1412/HYG 20842/MN 260/θ2 Tau          | Chamukuy       | θ2 Tau  | 78 Tau      | Tau             |   4.47771 |    15.8709 |       180.348 |      -22.0104 | 108.66 |  -26.39 |     40   |    46.1042 |   3.4  |      3.411 |      3.391 | 0.179 | A7III    |     80.8351  | 17.2097  |  40.8716 | 12.6082  | -6.49e-06  |  4.718e-05 |  5.51e-06  |     20842 | nan          |            0 |            1 |
    |  301 |  302 | 28305 | 1409 | 20889 | nan       | Ain            | 74 Tau/74Eps Tau/Ain/HD 28305/HIP 20889/HR 1409/HYG 20837/MN 302/ε Tau                | Ain            | ε Tau   | 74 Tau      | Tau             |   4.47694 |    19.1804 |       177.596 |      -19.9237 | 107.23 |  -36.77 |     39   |    44.964  |   3.53 |    nan     |    nan     | 1.014 | K0III    |     68.1711  | 16.4884  |  39.1368 | 14.7728  | -5.89e-06  |  4.622e-05 |  5.53e-06  |     20837 | nan          |            0 |            0 |
    |  344 |  345 | 27371 | 1346 | 20205 | nan       | Prima Hyadum   | 54 Tau/54Gam Tau/HD 27371/HIP 20205/HR 1346/HYG 20155/MN 345/Prima Hyadum/γ Tau       | Prima Hyadum   | γ Tau   | 54 Tau      | Tau             |   4.32989 |    15.6276 |       179.084 |      -23.8155 | 115.29 |  -23.86 |     39   |    49.5295 |   3.65 |    nan     |    nan     | 0.981 | G8III    |     74.0628  | 20.1974  |  43.2117 | 13.3426  | -8.16e-06  |  4.792e-05 |  5.23e-06  |     20155 | nan          |            0 |            0 |
    |  394 |  395 | 27697 | 1373 | 20455 | nan       | Secunda Hyadum | 61 Tau/61Del1Tau/HD 27697/HIP 20455/HR 1373/HYG 20405/MN 395/Secunda Hyadum/δ1 Tau    | Secunda Hyadum | δ1 Tau  | 61 Tau      | Tau             |   4.38225 |    17.5425 |       178.015 |      -22.0086 | 107.75 |  -28.84 |     39   |    47.7099 |   3.77 |    nan     |    nan     | 0.983 | G8III    |     61.546   | 18.696   |  41.472  | 14.3805  | -6.26e-06  |  4.675e-05 |  5.66e-06  |     20405 | nan          |            0 |            0 |
    |  429 |  430 | 28307 | 1411 | 20885 | nan       | θ1 Tau         | 77 Tau/77The1Tau/HD 28307/HIP 20885/HR 1411/HYG 20833/MN 430/θ1 Tau                   | nan            | θ1 Tau  | 77 Tau      | Tau             |   4.47625 |    15.9622 |       180.257 |      -21.9696 | 104.76 |  -15.01 |     40   |    47.3261 |   3.84 |    nan     |    nan     | 0.952 | G7III    |     56.8068  | 17.6738  |  41.929  | 13.0149  | -6.5e-06   |  4.645e-05 |  7.94e-06  |     20833 | nan          |            0 |            0 |
    |  556 |  557 | 31421 | 1580 | 22957 | nan       | ο2 Ori         | 9 Ori/9Omi2Ori/HD 31421/HIP 22957/HR 1580/HYG 22904/MN 557/ο2 Ori                     | nan            | ο2 Ori  | 9 Ori       | Ori             |   4.93952 |    13.5145 |       186.626 |      -18.0522 | -77.77 |  -45.97 |      1   |    57.0125 |   4.06 |    nan     |    nan     | 1.158 | K2III    |     67.2977  | 15.1933  |  53.3112 | 13.3233  |  2.176e-05 | -2.08e-06  | -1.211e-05 |     22904 | nan          |            0 |            0 |
    |  693 |  694 | 29388 | 1473 | 21589 | nan       | 90 Tau         | 90 Tau/HD 29388/HIP 21589/HR 1473/HYG 21536/MN 694                                    | nan            | nan     | 90 Tau      | Tau             |   4.63596 |    12.5108 |       184.737 |      -22.2407 | 101.73 |  -14.9  |     45   |    47.081  |   4.27 |    nan     |    nan     | 0.122 | A6V      |     37.8094  | 16.0671  |  43.0638 | 10.199   | -5.79e-06  |  5.09e-05  |  6.65e-06  |     21536 | nan          |            0 |            0 |
    |  721 |  722 | 27962 | 1389 | 20648 | nan       | δ3 Tau         | 68 Tau/68Del3Tau/HD 27962/HIP 20648/HR 1389/HYG 20597/MN 722/δ3 Tau                   | nan            | δ3 Tau  | 68 Tau      | Tau             |   4.42483 |    17.9279 |       178.118 |      -21.2948 | 108.26 |  -32.47 |     35   |    45.5373 |   4.3  |      4.309 |      4.299 | 0.049 | A2IV     |     34.4191  | 17.3648  |  39.6945 | 14.0174  | -7.36e-06  |  4.28e-05  |  4.2e-06   |     20597 | nan          |            0 |            1 |
    |  894 |  895 | 28052 | 1394 | 20713 | nan       | 71 Tau         | 71 Tau/HD 28052/HIP 20713/HR 1394/HYG 20661/MN 895                                    | nan            | nan     | 71 Tau      | Tau             |   4.43909 |    15.6183 |       180.182 |      -22.6029 | 114.66 |  -33.3  |     38   |    49.0918 |   4.48 |      4.488 |      4.468 | 0.262 | F0V...   |     33.8844  | 18.7872  |  43.3865 | 13.2169  | -9.32e-06  |  4.715e-05 |  2.83e-06  |     20661 | nan          |            0 |            1 |
    | 1072 | 1073 | 28910 | 1444 | 21273 | nan       | ρ Tau          | 86 Tau/86Rho Tau/HD 28910/HIP 21273/HR 1444/HYG 21220/MN 1073/ρ Tau                   | nan            | ρ Tau   | 86 Tau      | Tau             |   4.56414 |    14.8444 |       182.05  |      -21.665  | 103.69 |  -25.94 |     40   |    48.5201 |   4.65 |      4.662 |      4.642 | 0.255 | A8V      |     28.3139  | 17.2181  |  43.6262 | 12.4307  | -7.59e-06  |  4.719e-05 |  4.58e-06  |     21220 | nan          |            0 |            1 |
    | 1100 | 1101 | 29488 | 1479 | 21683 | nan       | σ2 Tau         | 92 Tau/92Sig2Tau/HD 29488/HIP 21683/HR 1479/HYG 21630/MN 1101/σ2 Tau                  | nan            | σ2 Tau  | 92 Tau      | Tau             |   4.65458 |    15.918  |       181.995 |      -19.9736 |  82.4  |  -19.53 |     36   |    47.6872 |   4.67 |    nan     |    nan     | 0.147 | A5Vn     |     26.8411  | 15.8209  |  43.0435 | 13.0788  | -5.24e-06  |  4.097e-05 |  5.76e-06  |     21630 | nan          |            0 |            0 |
    | 1123 | 1124 | 28100 | 1396 | 20732 | nan       | π Tau          | 73 Tau/73Pi Tau/HD 28100/HIP 20732/HR 1396/HYG 20680/MN 1124/π Tau                    | nan            | π Tau   | 73 Tau      | Tau             |   4.44344 |    14.7138 |       180.989 |      -23.1211 |  -7.57 |  -31.15 |     32   |   127.714  |   4.69 |    nan     |    nan     | 0.979 | G8III    |    188.973   | 48.9561  | 113.411  | 32.4382  |  1.879e-05 |  3.17e-05  | -1.034e-05 |     20680 | nan          |            0 |            0 |
    | 1144 | 1145 | 30959 | 1556 | 22667 | nan       | ο1 Ori         | 4 Ori/4Omi1Ori/HD 30959/HIP 22667/HR 1556/HYG 22614/MN 1145/ο1 Ori                    | nan            | ο1 Ori  | 4 Ori       | Ori             |   4.87554 |    14.2506 |       185.428 |      -18.3898 |  -2.62 |  -56.13 |     -8   |   199.601  |   4.71 |      4.758 |      4.668 | 1.773 | M3Sv     |    453.315   | 56.1317  | 185.136  | 49.1345  |  4e-06     |  4.47e-06  | -5.465e-05 |     22614 | nan          |            0 |            1 |
    | 1242 | 1243 | 28527 | 1427 | 21029 | Gl 170.1  | HD 28527       | Gl 170.1/HD 28527/HIP 21029/HR 1427/HYG 20977/MN 1243                                 | nan            | nan     | nan         | Tau             |   4.50934 |    16.194  |       180.383 |      -21.4524 | 104.98 |  -25.14 |     41   |    43.1965 |   4.78 |    nan     |    nan     | 0.17  | A6IV     |     19.8976  | 15.781   |  38.364  | 12.0472  | -4.45e-06  |  4.696e-05 |  6.64e-06  |     20977 | nan          |            0 |            0 |
    | 1275 | 1276 | 27819 | 1380 | 20542 | nan       | δ2 Tau         | 64 Tau/64Del2Tau/HD 27819/HIP 20542/HR 1380/HYG 20490/MN 1276/δ2 Tau                  | nan            | δ2 Tau  | 64 Tau      | Tau             |   4.4016  |    17.4441 |       178.288 |      -21.8596 | 109.99 |  -33.47 |     39   |    49.4805 |   4.8  |    nan     |    nan     | 0.154 | A7V      |     25.633   | 19.182   |  43.1321 | 14.8332  | -7.67e-06  |  4.769e-05 |  4.3e-06   |     20490 | nan          |            0 |            0 |
    | 1485 | 1486 | 27045 | 1329 | 19990 | nan       | ω2 Tau         | 50 Tau/50Ome2Tau/HD 27045/HIP 19990/HR 1329/HYG 19940/MN 1486/ω2 Tau                  | nan            | ω2 Tau  | 50 Tau      | Tau             |   4.28768 |    20.5786 |       174.605 |      -21.0331 | -40.88 |  -61.45 |     16   |    28.9436 |   4.93 |    nan     |    nan     | 0.259 | A3m      |      7.78395 | 11.7443  |  24.4195 | 10.1735  |  1.312e-05 |  1.405e-05 | -2.32e-06  |     19940 | nan          |            0 |            0 |
    | 1544 | 1545 | 28292 | 1407 | 20877 | nan       | 75 Tau         | 75 Tau/HD 28292/HIP 20877/HR 1407/HYG 20825/MN 1545                                   | nan            | nan     | 75 Tau      | Tau             |   4.47399 |    16.3597 |       179.901 |      -21.7454 |   8.12 |   17.79 |     18   |    57.241  |   4.96 |    nan     |    nan     | 1.137 | K2IIIvar |     29.621   | 21.3634  |  50.5985 | 16.1229  |  4.25e-06  |  1.587e-05 |  9.92e-06  |     20825 | nan          |            0 |            0 |



We can map the stars:


```python
fig,axs = hyades.plot_stars()
fig.savefig('gallery/hyades.png')
plt.close(fig) # Used only for README generation

```

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/gallery/hyades.png?raw=true" alt="Logo""/></p>

Now you can precess the stars:


```python
mtime = montu.Time('-2500-01-01 12:00:00.00')
hyades.where_in_space(at=mtime,inplace=True)
hyades
```




    18 star(s):
    |      |   MN |    HD |   HR |   HIP | Gl        | Name           | OtherDesignations                                                                     | ProperName     | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   GalLonJ2000 |   GalLatJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |    VXJ2000 |    VYJ2000 |    VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |           tt |    jed |   RAJ2000t |   DecJ2000t |   RAEpoch |   DecEpoch |
    |------|------|-------|------|-------|-----------|----------------|---------------------------------------------------------------------------------------|----------------|---------|-------------|-----------------|-----------|------------|---------------|---------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|------------|------------|------------|-----------|--------------|--------------|--------------|--------------|--------|------------|-------------|-----------|------------|
    |   14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran      | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran      | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |       180.972 |      -20.2483 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |    163.23    |  7.02722 |  18.2876 |  5.80666 |  1.528e-05 |  5.709e-05 | -2.14e-06  |     21368 | Gl 171.1     |            1 |            1 | -1.42006e+11 | 807954 |    4.59345 |     16.746  |  0.610393 |  -2.25214  |
    |  259 |  260 | 28319 | 1412 | 20894 | nan       | Chamukuy       | 78 Tau/78The2Tau/Chamukuy/HD 28319/HIP 20894/HR 1412/HYG 20842/MN 260/θ2 Tau          | Chamukuy       | θ2 Tau  | 78 Tau      | Tau             |   4.47771 |    15.8709 |       180.348 |      -22.0104 | 108.66 |  -26.39 |     40   |    46.1042 |   3.4  |      3.411 |      3.391 | 0.179 | A7III    |     80.8351  | 17.2097  |  40.8716 | 12.6082  | -6.49e-06  |  4.718e-05 |  5.51e-06  |     20842 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.46865 |     15.9039 |  0.509023 |  -3.52894  |
    |  301 |  302 | 28305 | 1409 | 20889 | nan       | Ain            | 74 Tau/74Eps Tau/Ain/HD 28305/HIP 20889/HR 1409/HYG 20837/MN 302/ε Tau                | Ain            | ε Tau   | 74 Tau      | Tau             |   4.47694 |    19.1804 |       177.596 |      -19.9237 | 107.23 |  -36.77 |     39   |    44.964  |   3.53 |    nan     |    nan     | 1.014 | K0III    |     68.1711  | 16.4884  |  39.1368 | 14.7728  | -5.89e-06  |  4.622e-05 |  5.53e-06  |     20837 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.46801 |     19.2264 |  0.452486 |  -0.316303 |
    |  344 |  345 | 27371 | 1346 | 20205 | nan       | Prima Hyadum   | 54 Tau/54Gam Tau/HD 27371/HIP 20205/HR 1346/HYG 20155/MN 345/Prima Hyadum/γ Tau       | Prima Hyadum   | γ Tau   | 54 Tau      | Tau             |   4.32989 |    15.6276 |       179.084 |      -23.8155 | 115.29 |  -23.86 |     39   |    49.5295 |   3.65 |    nan     |    nan     | 0.981 | G8III    |     74.0628  | 20.1974  |  43.2117 | 13.3426  | -8.16e-06  |  4.792e-05 |  5.23e-06  |     20155 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.32028 |     15.6575 |  0.374386 |  -4.29558  |
    |  394 |  395 | 27697 | 1373 | 20455 | nan       | Secunda Hyadum | 61 Tau/61Del1Tau/HD 27697/HIP 20455/HR 1373/HYG 20405/MN 395/Secunda Hyadum/δ1 Tau    | Secunda Hyadum | δ1 Tau  | 61 Tau      | Tau             |   4.38225 |    17.5425 |       178.015 |      -22.0086 | 107.75 |  -28.84 |     39   |    47.7099 |   3.77 |    nan     |    nan     | 0.983 | G8III    |     61.546   | 18.696   |  41.472  | 14.3805  | -6.26e-06  |  4.675e-05 |  5.66e-06  |     20405 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.37327 |     17.5786 |  0.392658 |  -2.24724  |
    |  429 |  430 | 28307 | 1411 | 20885 | nan       | θ1 Tau         | 77 Tau/77The1Tau/HD 28307/HIP 20885/HR 1411/HYG 20833/MN 430/θ1 Tau                   | nan            | θ1 Tau  | 77 Tau      | Tau             |   4.47625 |    15.9622 |       180.257 |      -21.9696 | 104.76 |  -15.01 |     40   |    47.3261 |   3.84 |    nan     |    nan     | 0.952 | G7III    |     56.8068  | 17.6738  |  41.929  | 13.0149  | -6.5e-06   |  4.645e-05 |  7.94e-06  |     20833 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.46752 |     15.9809 |  0.506667 |  -3.45849  |
    |  556 |  557 | 31421 | 1580 | 22957 | nan       | ο2 Ori         | 9 Ori/9Omi2Ori/HD 31421/HIP 22957/HR 1580/HYG 22904/MN 557/ο2 Ori                     | nan            | ο2 Ori  | 9 Ori       | Ori             |   4.93952 |    13.5145 |       186.626 |      -18.0522 | -77.77 |  -45.97 |      1   |    57.0125 |   4.06 |    nan     |    nan     | 1.158 | K2III    |     67.2977  | 15.1933  |  53.3112 | 13.3233  |  2.176e-05 | -2.08e-06  | -1.211e-05 |     22904 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.946   |     13.5719 |  0.996528 |  -3.88634  |
    |  693 |  694 | 29388 | 1473 | 21589 | nan       | 90 Tau         | 90 Tau/HD 29388/HIP 21589/HR 1473/HYG 21536/MN 694                                    | nan            | nan     | 90 Tau      | Tau             |   4.63596 |    12.5108 |       184.737 |      -22.2407 | 101.73 |  -14.9  |     45   |    47.081  |   4.27 |    nan     |    nan     | 0.122 | A6V      |     37.8094  | 16.0671  |  43.0638 | 10.199   | -5.79e-06  |  5.09e-05  |  6.65e-06  |     21536 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.62748 |     12.5295 |  0.716907 |  -6.18727  |
    |  721 |  722 | 27962 | 1389 | 20648 | nan       | δ3 Tau         | 68 Tau/68Del3Tau/HD 27962/HIP 20648/HR 1389/HYG 20597/MN 722/δ3 Tau                   | nan            | δ3 Tau  | 68 Tau      | Tau             |   4.42483 |    17.9279 |       178.118 |      -21.2948 | 108.26 |  -32.47 |     35   |    45.5373 |   4.3  |      4.309 |      4.299 | 0.049 | A2IV     |     34.4191  | 17.3648  |  39.6945 | 14.0174  | -7.36e-06  |  4.28e-05  |  4.2e-06   |     20597 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.41581 |     17.9685 |  0.425544 |  -1.71979  |
    |  894 |  895 | 28052 | 1394 | 20713 | nan       | 71 Tau         | 71 Tau/HD 28052/HIP 20713/HR 1394/HYG 20661/MN 895                                    | nan            | nan     | 71 Tau      | Tau             |   4.43909 |    15.6183 |       180.182 |      -22.6029 | 114.66 |  -33.3  |     38   |    49.0918 |   4.48 |      4.488 |      4.468 | 0.262 | F0V...   |     33.8844  | 18.7872  |  43.3865 | 13.2169  | -9.32e-06  |  4.715e-05 |  2.83e-06  |     20661 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.42954 |     15.6599 |  0.476608 |  -3.90687  |
    | 1072 | 1073 | 28910 | 1444 | 21273 | nan       | ρ Tau          | 86 Tau/86Rho Tau/HD 28910/HIP 21273/HR 1444/HYG 21220/MN 1073/ρ Tau                   | nan            | ρ Tau   | 86 Tau      | Tau             |   4.56414 |    14.8444 |       182.05  |      -21.665  | 103.69 |  -25.94 |     40   |    48.5201 |   4.65 |      4.662 |      4.642 | 0.255 | A8V      |     28.3139  | 17.2181  |  43.6262 | 12.4307  | -7.59e-06  |  4.719e-05 |  4.58e-06  |     21220 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.5555  |     14.8768 |  0.607744 |  -4.19944  |
    | 1100 | 1101 | 29488 | 1479 | 21683 | nan       | σ2 Tau         | 92 Tau/92Sig2Tau/HD 29488/HIP 21683/HR 1479/HYG 21630/MN 1101/σ2 Tau                  | nan            | σ2 Tau  | 92 Tau      | Tau             |   4.65458 |    15.918  |       181.995 |      -19.9736 |  82.4  |  -19.53 |     36   |    47.6872 |   4.67 |    nan     |    nan     | 0.147 | A5Vn     |     26.8411  | 15.8209  |  43.0435 | 13.0788  | -5.24e-06  |  4.097e-05 |  5.76e-06  |     21630 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.64772 |     15.9424 |  0.674876 |  -2.81982  |
    | 1123 | 1124 | 28100 | 1396 | 20732 | nan       | π Tau          | 73 Tau/73Pi Tau/HD 28100/HIP 20732/HR 1396/HYG 20680/MN 1124/π Tau                    | nan            | π Tau   | 73 Tau      | Tau             |   4.44344 |    14.7138 |       180.989 |      -23.1211 |  -7.57 |  -31.15 |     32   |   127.714  |   4.69 |    nan     |    nan     | 0.979 | G8III    |    188.973   | 48.9561  | 113.411  | 32.4382  |  1.879e-05 |  3.17e-05  | -1.034e-05 |     20680 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.44407 |     14.7527 |  0.505399 |  -4.73255  |
    | 1144 | 1145 | 30959 | 1556 | 22667 | nan       | ο1 Ori         | 4 Ori/4Omi1Ori/HD 30959/HIP 22667/HR 1556/HYG 22614/MN 1145/ο1 Ori                    | nan            | ο1 Ori  | 4 Ori       | Ori             |   4.87554 |    14.2506 |       185.428 |      -18.3898 |  -2.62 |  -56.13 |     -8   |   199.601  |   4.71 |      4.758 |      4.668 | 1.773 | M3Sv     |    453.315   | 56.1317  | 185.136  | 49.1345  |  4e-06     |  4.47e-06  | -5.465e-05 |     22614 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.87576 |     14.3208 |  0.91669  |  -3.46436  |
    | 1242 | 1243 | 28527 | 1427 | 21029 | Gl 170.1  | HD 28527       | Gl 170.1/HD 28527/HIP 21029/HR 1427/HYG 20977/MN 1243                                 | nan            | nan     | nan         | Tau             |   4.50934 |    16.194  |       180.383 |      -21.4524 | 104.98 |  -25.14 |     41   |    43.1965 |   4.78 |    nan     |    nan     | 0.17  | A6IV     |     19.8976  | 15.781   |  38.364  | 12.0472  | -4.45e-06  |  4.696e-05 |  6.64e-06  |     20977 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.50059 |     16.2254 |  0.533307 |  -3.10096  |
    | 1275 | 1276 | 27819 | 1380 | 20542 | nan       | δ2 Tau         | 64 Tau/64Del2Tau/HD 27819/HIP 20542/HR 1380/HYG 20490/MN 1276/δ2 Tau                  | nan            | δ2 Tau  | 64 Tau      | Tau             |   4.4016  |    17.4441 |       178.288 |      -21.8596 | 109.99 |  -33.47 |     39   |    49.4805 |   4.8  |    nan     |    nan     | 0.154 | A7V      |     25.633   | 19.182   |  43.1321 | 14.8332  | -7.67e-06  |  4.769e-05 |  4.3e-06   |     20490 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.39243 |     17.486  |  0.411904 |  -2.26983  |
    | 1485 | 1486 | 27045 | 1329 | 19990 | nan       | ω2 Tau         | 50 Tau/50Ome2Tau/HD 27045/HIP 19990/HR 1329/HYG 19940/MN 1486/ω2 Tau                  | nan            | ω2 Tau  | 50 Tau      | Tau             |   4.28768 |    20.5786 |       174.605 |      -21.0331 | -40.88 |  -61.45 |     16   |    28.9436 |   4.93 |    nan     |    nan     | 0.259 | A3m      |      7.78395 | 11.7443  |  24.4195 | 10.1735  |  1.312e-05 |  1.405e-05 | -2.32e-06  |     19940 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.29109 |     20.6554 |  0.26793  |   0.459077 |
    | 1544 | 1545 | 28292 | 1407 | 20877 | nan       | 75 Tau         | 75 Tau/HD 28292/HIP 20877/HR 1407/HYG 20825/MN 1545                                   | nan            | nan     | 75 Tau      | Tau             |   4.47399 |    16.3597 |       179.901 |      -21.7454 |   8.12 |   17.79 |     18   |    57.241  |   4.96 |    nan     |    nan     | 1.137 | K2IIIvar |     29.621   | 21.3634  |  50.5985 | 16.1229  |  4.25e-06  |  1.587e-05 |  9.92e-06  |     20825 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.47332 |     16.3374 |  0.506045 |  -3.09247  |



Or compute their horizontal positions:


```python
hyades.where_in_sky(at=mtime,observer=tebas,inplace=True)
hyades
```




    18 star(s):
    |      |   MN |    HD |   HR |   HIP | Gl        | Name           | OtherDesignations                                                                     | ProperName     | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   GalLonJ2000 |   GalLatJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |   B-V | SpType   |   Luminosity |   XJ2000 |   YJ2000 |   ZJ2000 |    VXJ2000 |    VYJ2000 |    VZJ2000 |   Primary | MultipleID   |   IsMultiple |   IsVariable |           tt |    jed |   RAJ2000t |   DecJ2000t |   RAEpoch |   DecEpoch |      HA |      az |      el |     zen |
    |------|------|-------|------|-------|-----------|----------------|---------------------------------------------------------------------------------------|----------------|---------|-------------|-----------------|-----------|------------|---------------|---------------|--------|---------|----------|------------|--------|------------|------------|-------|----------|--------------|----------|----------|----------|------------|------------|------------|-----------|--------------|--------------|--------------|--------------|--------|------------|-------------|-----------|------------|---------|---------|---------|---------|
    |   14 |   15 | 29139 | 1457 | 21421 | Gl 171.1A | Aldebaran      | 87 Tau/87Alp Tau/Aldebaran/Gl 171.1A/HD 29139/HIP 21421/HR 1457/HYG 21368/MN 15/α Tau | Aldebaran      | α Tau   | 87 Tau      | Tau             |   4.59868 |    16.5093 |       180.972 |      -20.2483 |  62.78 | -189.36 |     54.5 |    20.4332 |   0.87 |      0.888 |      0.858 | 1.538 | K5III    |    163.23    |  7.02722 |  18.2876 |  5.80666 |  1.528e-05 |  5.709e-05 | -2.14e-06  |     21368 | Gl 171.1     |            1 |            1 | -1.42006e+11 | 807954 |    4.59345 |     16.746  |  0.610393 |  -2.25214  | 20.2633 | 107.596 | 29.5909 | 60.4091 |
    |  259 |  260 | 28319 | 1412 | 20894 | nan       | Chamukuy       | 78 Tau/78The2Tau/Chamukuy/HD 28319/HIP 20894/HR 1412/HYG 20842/MN 260/θ2 Tau          | Chamukuy       | θ2 Tau  | 78 Tau      | Tau             |   4.47771 |    15.8709 |       180.348 |      -22.0104 | 108.66 |  -26.39 |     40   |    46.1042 |   3.4  |      3.411 |      3.391 | 0.179 | A7III    |     80.8351  | 17.2097  |  40.8716 | 12.6082  | -6.49e-06  |  4.718e-05 |  5.51e-06  |     20842 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.46865 |     15.9039 |  0.509023 |  -3.52894  | 20.3647 | 109.747 | 30.2697 | 59.7303 |
    |  301 |  302 | 28305 | 1409 | 20889 | nan       | Ain            | 74 Tau/74Eps Tau/Ain/HD 28305/HIP 20889/HR 1409/HYG 20837/MN 302/ε Tau                | Ain            | ε Tau   | 74 Tau      | Tau             |   4.47694 |    19.1804 |       177.596 |      -19.9237 | 107.23 |  -36.77 |     39   |    44.964  |   3.53 |    nan     |    nan     | 1.014 | K0III    |     68.1711  | 16.4884  |  39.1368 | 14.7728  | -5.89e-06  |  4.622e-05 |  5.53e-06  |     20837 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.46801 |     19.2264 |  0.452486 |  -0.316303 | 20.4212 | 106.974 | 32.6028 | 57.3972 |
    |  344 |  345 | 27371 | 1346 | 20205 | nan       | Prima Hyadum   | 54 Tau/54Gam Tau/HD 27371/HIP 20205/HR 1346/HYG 20155/MN 345/Prima Hyadum/γ Tau       | Prima Hyadum   | γ Tau   | 54 Tau      | Tau             |   4.32989 |    15.6276 |       179.084 |      -23.8155 | 115.29 |  -23.86 |     39   |    49.5295 |   3.65 |    nan     |    nan     | 0.981 | G8III    |     74.0628  | 20.1974  |  43.2117 | 13.3426  | -8.16e-06  |  4.792e-05 |  5.23e-06  |     20155 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.32028 |     15.6575 |  0.374386 |  -4.29558  | 20.4993 | 111.726 | 31.5992 | 58.4008 |
    |  394 |  395 | 27697 | 1373 | 20455 | nan       | Secunda Hyadum | 61 Tau/61Del1Tau/HD 27697/HIP 20455/HR 1373/HYG 20405/MN 395/Secunda Hyadum/δ1 Tau    | Secunda Hyadum | δ1 Tau  | 61 Tau      | Tau             |   4.38225 |    17.5425 |       178.015 |      -22.0086 | 107.75 |  -28.84 |     39   |    47.7099 |   3.77 |    nan     |    nan     | 0.983 | G8III    |     61.546   | 18.696   |  41.472  | 14.3805  | -6.26e-06  |  4.675e-05 |  5.66e-06  |     20405 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.37327 |     17.5786 |  0.392658 |  -2.24724  | 20.481  | 109.49  | 32.4222 | 57.5778 |
    |  429 |  430 | 28307 | 1411 | 20885 | nan       | θ1 Tau         | 77 Tau/77The1Tau/HD 28307/HIP 20885/HR 1411/HYG 20833/MN 430/θ1 Tau                   | nan            | θ1 Tau  | 77 Tau      | Tau             |   4.47625 |    15.9622 |       180.257 |      -21.9696 | 104.76 |  -15.01 |     40   |    47.3261 |   3.84 |    nan     |    nan     | 0.952 | G7III    |     56.8068  | 17.6738  |  41.929  | 13.0149  | -6.5e-06   |  4.645e-05 |  7.94e-06  |     20833 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.46752 |     15.9809 |  0.506667 |  -3.45849  | 20.367  | 109.698 | 30.3358 | 59.6642 |
    |  556 |  557 | 31421 | 1580 | 22957 | nan       | ο2 Ori         | 9 Ori/9Omi2Ori/HD 31421/HIP 22957/HR 1580/HYG 22904/MN 557/ο2 Ori                     | nan            | ο2 Ori  | 9 Ori       | Ori             |   4.93952 |    13.5145 |       186.626 |      -18.0522 | -77.77 |  -45.97 |      1   |    57.0125 |   4.06 |    nan     |    nan     | 1.158 | K2III    |     67.2977  | 15.1933  |  53.3112 | 13.3233  |  2.176e-05 | -2.08e-06  | -1.211e-05 |     22904 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.946   |     13.5719 |  0.996528 |  -3.88634  | 19.8772 | 106.071 | 23.7372 | 66.2628 |
    |  693 |  694 | 29388 | 1473 | 21589 | nan       | 90 Tau         | 90 Tau/HD 29388/HIP 21589/HR 1473/HYG 21536/MN 694                                    | nan            | nan     | 90 Tau      | Tau             |   4.63596 |    12.5108 |       184.737 |      -22.2407 | 101.73 |  -14.9  |     45   |    47.081  |   4.27 |    nan     |    nan     | 0.122 | A6V      |     37.8094  | 16.0671  |  43.0638 | 10.199   | -5.79e-06  |  5.09e-05  |  6.65e-06  |     21536 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.62748 |     12.5295 |  0.716907 |  -6.18727  | 20.1568 | 110.55  | 26.2424 | 63.7576 |
    |  721 |  722 | 27962 | 1389 | 20648 | nan       | δ3 Tau         | 68 Tau/68Del3Tau/HD 27962/HIP 20648/HR 1389/HYG 20597/MN 722/δ3 Tau                   | nan            | δ3 Tau  | 68 Tau      | Tau             |   4.42483 |    17.9279 |       178.118 |      -21.2948 | 108.26 |  -32.47 |     35   |    45.5373 |   4.3  |      4.309 |      4.299 | 0.049 | A2IV     |     34.4191  | 17.3648  |  39.6945 | 14.0174  | -7.36e-06  |  4.28e-05  |  4.2e-06   |     20597 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.41581 |     17.9685 |  0.425544 |  -1.71979  | 20.4482 | 108.657 | 32.2619 | 57.7381 |
    |  894 |  895 | 28052 | 1394 | 20713 | nan       | 71 Tau         | 71 Tau/HD 28052/HIP 20713/HR 1394/HYG 20661/MN 895                                    | nan            | nan     | 71 Tau      | Tau             |   4.43909 |    15.6183 |       180.182 |      -22.6029 | 114.66 |  -33.3  |     38   |    49.0918 |   4.48 |      4.488 |      4.468 | 0.262 | F0V...   |     33.8844  | 18.7872  |  43.3865 | 13.2169  | -9.32e-06  |  4.715e-05 |  2.83e-06  |     20661 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.42954 |     15.6599 |  0.476608 |  -3.90687  | 20.3971 | 110.411 | 30.4939 | 59.5061 |
    | 1072 | 1073 | 28910 | 1444 | 21273 | nan       | ρ Tau          | 86 Tau/86Rho Tau/HD 28910/HIP 21273/HR 1444/HYG 21220/MN 1073/ρ Tau                   | nan            | ρ Tau   | 86 Tau      | Tau             |   4.56414 |    14.8444 |       182.05  |      -21.665  | 103.69 |  -25.94 |     40   |    48.5201 |   4.65 |      4.662 |      4.642 | 0.255 | A8V      |     28.3139  | 17.2181  |  43.6262 | 12.4307  | -7.59e-06  |  4.719e-05 |  4.58e-06  |     21220 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.5555  |     14.8768 |  0.607744 |  -4.19944  | 20.266  | 109.552 | 28.6561 | 61.3439 |
    | 1100 | 1101 | 29488 | 1479 | 21683 | nan       | σ2 Tau         | 92 Tau/92Sig2Tau/HD 29488/HIP 21683/HR 1479/HYG 21630/MN 1101/σ2 Tau                  | nan            | σ2 Tau  | 92 Tau      | Tau             |   4.65458 |    15.918  |       181.995 |      -19.9736 |  82.4  |  -19.53 |     36   |    47.6872 |   4.67 |    nan     |    nan     | 0.147 | A5Vn     |     26.8411  | 15.8209  |  43.0435 | 13.0788  | -5.24e-06  |  4.097e-05 |  5.76e-06  |     21630 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.64772 |     15.9424 |  0.674876 |  -2.81982  | 20.1988 | 107.62  | 28.4704 | 61.5296 |
    | 1123 | 1124 | 28100 | 1396 | 20732 | nan       | π Tau          | 73 Tau/73Pi Tau/HD 28100/HIP 20732/HR 1396/HYG 20680/MN 1124/π Tau                    | nan            | π Tau   | 73 Tau      | Tau             |   4.44344 |    14.7138 |       180.989 |      -23.1211 |  -7.57 |  -31.15 |     32   |   127.714  |   4.69 |    nan     |    nan     | 0.979 | G8III    |    188.973   | 48.9561  | 113.411  | 32.4382  |  1.879e-05 |  3.17e-05  | -1.034e-05 |     20680 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.44407 |     14.7527 |  0.505399 |  -4.73255  | 20.3683 | 110.972 | 29.6994 | 60.3006 |
    | 1144 | 1145 | 30959 | 1556 | 22667 | nan       | ο1 Ori         | 4 Ori/4Omi1Ori/HD 30959/HIP 22667/HR 1556/HYG 22614/MN 1145/ο1 Ori                    | nan            | ο1 Ori  | 4 Ori       | Ori             |   4.87554 |    14.2506 |       185.428 |      -18.3898 |  -2.62 |  -56.13 |     -8   |   199.601  |   4.71 |      4.758 |      4.668 | 1.773 | M3Sv     |    453.315   | 56.1317  | 185.136  | 49.1345  |  4e-06     |  4.47e-06  | -5.465e-05 |     22614 | nan          |            0 |            1 | -1.42006e+11 | 807954 |    4.87576 |     14.3208 |  0.91669  |  -3.46436  | 19.957  | 106.289 | 24.9889 | 65.0111 |
    | 1242 | 1243 | 28527 | 1427 | 21029 | Gl 170.1  | HD 28527       | Gl 170.1/HD 28527/HIP 21029/HR 1427/HYG 20977/MN 1243                                 | nan            | nan     | nan         | Tau             |   4.50934 |    16.194  |       180.383 |      -21.4524 | 104.98 |  -25.14 |     41   |    43.1965 |   4.78 |    nan     |    nan     | 0.17  | A6IV     |     19.8976  | 15.781   |  38.364  | 12.0472  | -4.45e-06  |  4.696e-05 |  6.64e-06  |     20977 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.50059 |     16.2254 |  0.533307 |  -3.10096  | 20.3404 | 109.107 | 30.1721 | 59.8279 |
    | 1275 | 1276 | 27819 | 1380 | 20542 | nan       | δ2 Tau         | 64 Tau/64Del2Tau/HD 27819/HIP 20542/HR 1380/HYG 20490/MN 1276/δ2 Tau                  | nan            | δ2 Tau  | 64 Tau      | Tau             |   4.4016  |    17.4441 |       178.288 |      -21.8596 | 109.99 |  -33.47 |     39   |    49.4805 |   4.8  |    nan     |    nan     | 0.154 | A7V      |     25.633   | 19.182   |  43.1321 | 14.8332  | -7.67e-06  |  4.769e-05 |  4.3e-06   |     20490 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.39243 |     17.486  |  0.411904 |  -2.26983  | 20.4618 | 109.34  | 32.162  | 57.838  |
    | 1485 | 1486 | 27045 | 1329 | 19990 | nan       | ω2 Tau         | 50 Tau/50Ome2Tau/HD 27045/HIP 19990/HR 1329/HYG 19940/MN 1486/ω2 Tau                  | nan            | ω2 Tau  | 50 Tau      | Tau             |   4.28768 |    20.5786 |       174.605 |      -21.0331 | -40.88 |  -61.45 |     16   |    28.9436 |   4.93 |    nan     |    nan     | 0.259 | A3m      |      7.78395 | 11.7443  |  24.4195 | 10.1735  |  1.312e-05 |  1.405e-05 | -2.32e-06  |     19940 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.29109 |     20.6554 |  0.26793  |   0.459077 | 20.6058 | 107.794 | 35.3961 | 54.6039 |
    | 1544 | 1545 | 28292 | 1407 | 20877 | nan       | 75 Tau         | 75 Tau/HD 28292/HIP 20877/HR 1407/HYG 20825/MN 1545                                   | nan            | nan     | 75 Tau      | Tau             |   4.47399 |    16.3597 |       179.901 |      -21.7454 |   8.12 |   17.79 |     18   |    57.241  |   4.96 |    nan     |    nan     | 1.137 | K2IIIvar |     29.621   | 21.3634  |  50.5985 | 16.1229  |  4.25e-06  |  1.587e-05 |  9.92e-06  |     20825 | nan          |            0 |            0 | -1.42006e+11 | 807954 |    4.47332 |     16.3374 |  0.506045 |  -3.09247  | 20.3677 | 109.337 | 30.5291 | 59.4709 |



And plot them again:


```python
fig,axs = hyades.plot_stars(coords=['RAEpoch','DecEpoch'])
fig.savefig('gallery/hyades-precessed.png')
plt.close(fig) # Used only for README generation

```

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/gallery/hyades-precessed.png?raw=true" alt="Logo""/></p>

### Working with planets and observing sites

`MontuPython` allows calculate the position of all planets in the solar system, including the moon:


```python
mars = montu.Planet('Mars')
neptune = montu.Planet('neptune')
jupiter = montu.Planet('JUPITER')
moon = montu.Planet('Moon')
sun = montu.Planet('Sun')
```

You may create an observing site:


```python
tebas = montu.Observer(lon=33,lat=24,height=0)
```

The straight routine method to calculate position of the planet in the sky at any date is:


```python
mtime = montu.Time('-2500-01-01 12:00:00')
mars.where_in_sky(mtime,tebas)
mars.show_position()
```

    Mars — sky position
      Epoch: -2500-01-22 12:00:00 / -2500-01-01 12:00:00.000000  (JED 807954.000000)
      Site: Custom site — lat 24.000000°, lon 33.000000°, 0 m  (P=1013.25 mbar, T=15.0 °C)
      Name: Mars
      RA (J2000): 12:31:51.677 h
      Dec (J2000): 1.620449°
      RA (epoch): 08:32:09.718 h
      Dec (epoch): 24.111367°
      RA (geocentric): 08:32:09.640 h
      Dec (geocentric): 24.114107°
      Azimuth: 6.189974°
      Elevation: -41.638524°


You may also compute other atrometric conditions (rise time, set time, elongation, etc.):


```python
mars.conditions_in_sky(mtime,tebas)
mars.show_conditions()
```

    Mars — sky conditions
      Epoch: -2500-01-22 12:00:00 / -2500-01-01 12:00:00.000000  (JED 807954.000000)
      Site: Custom site — lat 24.000000°, lon 33.000000°, 0 m  (P=1013.25 mbar, T=15.0 °C)
      Name: Mars
      Hour angle: 12:20:15.604 h
      Visual magnitude: -1.13 mag
      Rise time (UTC): -2500-01-22 16:49:49
      Rise azimuth: 63.133665°
      Set time (UTC): -2500-01-23 06:24:24
      Set azimuth: 296.924743°
      Transit time (UTC): -2500-01-22 23:37:37
      Transit elevation: 89.841412°
      Elongation from Sun: -157.818359°
      Distance from Earth: 0.660488 AU
      Distance from Sun: 1.626107 AU
      Angular diameter: 14.171"
      Illuminated fraction: 98.65 %
      Heliocentric latitude: 1.888422°
      Heliocentric longitude: 111.258645°
      Heliocentric longitude (alt.): 111.258645°
      Circumpolar: no
      Never rises: no


All times are given in Julian Days. If you want the exact date of set use:


```python
montu.Time.get_date(mars.condition.rise_time)
```




    (-2501, 1, 22, 16, 49, 57.292158)



You may get also the position at different times and store it into a data frame:


```python
import numpy as np

# Reset the data stored in the planet
mars.reset_store()

# Loop on different times
for deltat in np.arange(0,1*montu.DAY,1*montu.HOUR):
    mars.where_in_sky(mtime + deltat,tebas,store=True)

# Show results
mars
```




    Object Mars positions:
    |            tt |    jed | Name   |   RAJ2000 |   DecJ2000 |   RAEpoch |   DecEpoch |   RAGeo |   DecGeo |        el |        az |
    |---------------|--------|--------|-----------|------------|-----------|------------|---------|----------|-----------|-----------|
    | -142006202700 | 807954 | Mars   |   12.531  |    1.62045 |   8.53603 |    24.1114 | 8.53601 |  24.1141 | -41.6385  |   6.18997 |
    | -142006199100 | 807954 | Mars   |   12.5304 |    1.62593 |   8.53532 |    24.1153 | 8.53524 |  24.1179 | -38.0804  |  23.5034  |
    | -142006195500 | 807954 | Mars   |   12.5297 |    1.63142 |   8.53461 |    24.1193 | 8.53447 |  24.1217 | -31.0231  |  37.8391  |
    | -142006191900 | 807954 | Mars   |   12.5291 |    1.63692 |   8.53388 |    24.1233 | 8.53369 |  24.1256 | -21.5453  |  48.9456  |
    | -142006188300 | 807954 | Mars   |   12.5284 |    1.64242 |   8.53314 |    24.1274 | 8.53292 |  24.1294 | -10.5149  |  57.4706  |
    | -142006184700 | 807954 | Mars   |   12.5277 |    1.64793 |   8.53238 |    24.1316 | 8.53214 |  24.1332 |   1.8128  |  64.1498  |
    | -142006181100 | 807954 | Mars   |   12.5271 |    1.65344 |   8.53161 |    24.1358 | 8.53136 |  24.137  |  14.211   |  69.5518  |
    | -142006177500 | 807954 | Mars   |   12.5264 |    1.65896 |   8.53081 |    24.14   | 8.53058 |  24.1409 |  27.2472  |  74.0721  |
    | -142006173900 | 807954 | Mars   |   12.5258 |    1.66448 |   8.53    |    24.1441 | 8.5298  |  24.1447 |  40.5806  |  77.979   |
    | -142006170300 | 807954 | Mars   |   12.5251 |    1.67001 |   8.52917 |    24.1482 | 8.52902 |  24.1485 |  54.1045  |  81.4499  |
    | -142006166700 | 807954 | Mars   |   12.5244 |    1.67555 |   8.52833 |    24.1523 | 8.52823 |  24.1524 |  67.7491  |  84.572   |
    | -142006163100 | 807954 | Mars   |   12.5238 |    1.68109 |   8.52748 |    24.1562 | 8.52744 |  24.1562 |  81.4615  |  87.0426  |
    | -142006159500 | 807954 | Mars   |   12.5231 |    1.68664 |   8.52663 |    24.1601 | 8.52666 |  24.1601 |  84.7964  | 272.924   |
    | -142006155900 | 807955 | Mars   |   12.5224 |    1.69219 |   8.52578 |    24.1638 | 8.52587 |  24.1639 |  71.0758  | 274.763   |
    | -142006152300 | 807955 | Mars   |   12.5218 |    1.69775 |   8.52493 |    24.1675 | 8.52508 |  24.1677 |  57.4117  | 277.794   |
    | -142006148700 | 807955 | Mars   |   12.5211 |    1.70331 |   8.52409 |    24.1711 | 8.52429 |  24.1716 |  43.8553  | 281.178   |
    | -142006145100 | 807955 | Mars   |   12.5204 |    1.70888 |   8.52327 |    24.1746 | 8.52349 |  24.1754 |  30.4713  | 284.969   |
    | -142006141500 | 807955 | Mars   |   12.5197 |    1.71446 |   8.52245 |    24.1781 | 8.5227  |  24.1793 |  17.3547  | 289.322   |
    | -142006137900 | 807955 | Mars   |   12.5191 |    1.72004 |   8.52166 |    24.1816 | 8.5219  |  24.1831 |   4.71203 | 294.483   |
    | -142006134300 | 807955 | Mars   |   12.5184 |    1.72562 |   8.52087 |    24.1851 | 8.52111 |  24.187  |  -7.62412 | 300.812   |
    | -142006130700 | 807955 | Mars   |   12.5177 |    1.73122 |   8.52011 |    24.1886 | 8.52031 |  24.1908 | -18.9437  | 308.83    |
    | -142006127100 | 807955 | Mars   |   12.517  |    1.73681 |   8.51935 |    24.1922 | 8.51951 |  24.1947 | -28.852   | 319.24    |
    | -142006123500 | 807955 | Mars   |   12.5163 |    1.74242 |   8.51861 |    24.1959 | 8.51871 |  24.1985 | -36.5767  | 332.753   |
    | -142006119900 | 807955 | Mars   |   12.5157 |    1.74803 |   8.51787 |    24.1997 | 8.5179  |  24.2024 | -41.0646  | 349.434   |'
    Object Mars conditions:
    '



### Working with the sun and the moon

There are two special planetary objects, the Sun and the Moon. Although you can define them as if they where planets:


```python
sun = montu.Planet('Sun')
moon = montu.Planet('Moon')
```

and as we did with planets, you can calculate positions and conditions:


```python
mtime = montu.Time('-2500-01-01 12:00:00')
Tebas = montu.Observer(lon=33,lat=24)
moon.conditions_in_sky(at=mtime,observer=Tebas)
moon.show_conditions()
```

    Moon — sky conditions
      Epoch: -2500-01-22 12:00:00 / -2500-01-01 12:00:00.000000  (JED 807954.000000)
      Site: Custom site — lat 24.000000°, lon 33.000000°, 0 m  (P=1013.25 mbar, T=15.0 °C)
      Name: Moon
      Hour angle: 11:40:41.857 h
      Visual magnitude: -11.92 mag
      Rise time (UTC): -2500-01-22 17:53:53
      Rise azimuth: 67.113965°
      Set time (UTC): -2500-01-23 07:34:34
      Set azimuth: 289.880687°
      Transit time (UTC): -2500-01-23 00:47:47
      Transit elevation: 85.450065°
      Elongation from Sun: -148.449280°
      Distance from Earth: 0.002553 AU
      Distance from Sun: 0.997457 AU
      Angular diameter: 1879.512"
      Illuminated fraction: 92.65 %
      Heliocentric latitude: 5.076801°
      Heliocentric longitude: 133.820563°
      Heliocentric longitude (alt.): 133.820563°
      Circumpolar: no
      Never rises: no


There are additional methods, both proper and static, that can be used with these objects. For instance, you can:

- calculate the time of twilight:


```python
dusk_time, dawn_time = montu.Sun.when_is_twilight(day=mtime,observer=Tebas,sunbelow=-18)
Tebas.get_local_time(dusk_time), Tebas.get_local_time(dawn_time) 
```




    ('05:34:16.262', '18:55:22.138')



- Get the date of solstices and equinoxes:


```python
mtime = montu.Time('2023-01-01 12:00:00')
vernal, summer, autumn, winter = montu.Sun.next_seasons(at=mtime)

(
    montu.Time.get_date(vernal,format='mtime'),
    montu.Time.get_date(summer,format='mtime'),
    montu.Time.get_date(autumn,format='mtime'),
    montu.Time.get_date(winter,format='mtime'),
)
```




    (Time('2023-03-20 21:25:28.5'/'2023-03-20 21:24:24'/'[hrw 4806] IV shemu 7'/JED 2460024.3918415/JTD 2460024.392691),
     Time('2023-06-21 14:59:09.1'/'2023-06-21 14:57:57'/'[hrw 4807] III akhet 5'/JED 2460117.123559/JTD 2460117.1244109),
     Time('2023-09-23 06:51:14.0'/'2023-09-23 06:50:50'/'[hrw 4807] II peret 9'/JED 2460210.7847257/JTD 2460210.7855787),
     Time('2023-12-22 03:28:23.1'/'2023-12-22 03:27:27'/'[hrw 4807] I shemu 9'/JED 2460300.6438576/JTD 2460300.6447118))



- Get the lunar phases:


```python
montu.Moon.next_moon_quarters(since=mtime,output='datepro')
```




    {'full': ['2023-01-06 23:09:10.797107', 7.907907465007156, 5.46470835711807],
     'last': ['2023-01-15 02:11:35.301130', 8.126672363374382, 13.591380720492452],
     'new': ['2023-01-21 20:54:23.100479', 6.7797199985943735, 20.371100719086826],
     'first': ['2023-01-28 15:20:01.694387',
      6.767808315809816,
      27.138909034896642]}



### An indepth examples

Short case studies: **Sirius heliacal rise** at the 139 CE Sothic apokatastasis, **Mars–Aldebarán** (September 2022), **Thales' eclipse** with Xavier Jubier cross-checks, and **polar stars** over 40 kyr.

### Heliacal rise of Sirius at the 139 CE apokatastasis


The Egyptian civil year begins on **I akhet 1**, when Sirius reappears in the dawn sky. The **139 CE apokatastasis** is one Sothic New Year in that cycle (`[hrw 2922] I akhet 1`). The `HeliacalRise` class (Schaefer 1987 model) predicts first-morning visibility of Sirius during the year anchored on that date:



```python
tebas = montu.Observer(site="thebes")
sirius = montu.Stars(subset="bright", ProperName="Sirius")

calculator = montu.HeliacalRise(
    model="ptolemy",
    arcus_visionis_crit=14.0,  # assumed Arcus Visionis
)

t_start = montu.Time("139-06-01")
events = calculator.compute(
    sirius,
    tebas,
    t_start,
    t_start.add(years=1),
)

calculator.print_rises(
    events,
    title="Heliacal rises of Sirius (139 CE apokatastasis)",
    body_label="Sirius",
)

```

    Heliacal rises of Sirius (139 CE apokatastasis) — 1 date(s)
      [1] 139-07-20 00:00:00  04:07:44.893  0139-07-19 00:00:00.000000  [hrw 2922] I akhet 1  Sirius -0.57°  Sun -14.18°
      source: Toomer, G. J. (1998). Ptolemy's Almagest. Princeton University Press. Book XIII, Chapter 7: "On the heliacal risings and settings of the planets".


The first heliacal rising falls on the same civil date as **I akhet 1** — the canonical anchor of the Egyptian year. See [`MontuPython-HeliacalRises.ipynb`](examples/MontuPython-HeliacalRises.ipynb) for a step-by-step treatment and alternative visibility models.


<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/montu_gui/assets/illustrations/solar-eclipses-illustration.webp?raw=true" alt="Illustration of a hypothetical total solar eclipse observed at Stonehenge" width="800"/></p>

<p align="center"><em>Illustration of a hypothetical total solar eclipse observed at Stonehenge</em>. Image generated with Nano Banana, edited in Preview.</p>

### Mars and Aldebarán Conjunction (September 2022)


The `Conjunction` and `ConjunctionExplorer` classes evaluate angular groupings of planets and stars. Here we recover the Mars–Aldebarán conjunction of 7 September 2022:



```python
mars = montu.Planet("Mars")
aldebaran = montu.Stars(subset="bright", ProperName="Aldebaran", return_as="Star")
explorer = montu.ConjunctionExplorer(bodies=[mars, aldebaran], maxseparation=5)
conjs = explorer.search(
    start=montu.Time("2022-09-01"),
    end=montu.Time("2022-10-01"),
    observer="geocentric",
)
conjs[0].show_details()

```

    Conjunction: Mars–Aldebaran
      Epoch (UTC)          : 2022-09-07 14:28:28
      Julian Day (UTC)     : 2459830.103152
      Observer             : geocentric
      Angular separation   : 4.2746° (max allowed 5.0°)
      In conjunction       : yes
      Is visible from site : n/a (geocentric)
      Pair Mars–Aldebaran
        Separation         : 4.2746°
        Position angle     : 170.19° (N→E)
      Mars
        Phase              : 85.42%
        Angular size       : 0.169 arcmin
        V magnitude        : -0.22
      Aldebaran
        V magnitude        : 0.87


See [`MontuPython-Conjunctions.ipynb`](examples/MontuPython-Conjunctions.ipynb) for lapse intervals, visibility checks, sky maps, and mixed planet–star groupings.


### Thales' eclipse (585 BCE)


The [NASA Five Millennium Canon of Solar Eclipses](https://eclipse.gsfc.nasa.gov/SEpubs/5MCSE.html) records the **total eclipse of 28 May 585 BCE** (catalogue year −584), associated with Thales of Miletus. MontuPython loads the catalogue row and computes local circumstances for any `Observer`.

We start by looking all eclipses in a given year interval:


```python
catalogue = montu.SolarEclipses()
eclipses = catalogue.get_eclipses(year=[-590, -580], eclipse_type='T')
TABLEDF(eclipses.data)

```

    |      |   year |   month |   day | td_ge    |      dt |   luna_num |   saros | eclipse_type   |    gamma |   magnitude | lat_ge   | lng_ge   |   lat_dd_ge |   lng_dd_ge |   sun_alt |   sun_azm |   path_width | central_duration   |   duration_secs |   cat_no |   canon_plate |   julian_date |   t0 |        x0 |       x1 |        x2 |        x3 |        y0 |         y1 |         y2 |        y3 |        d0 |        d1 |     d2 |      mu0 |     mu1 |   mu2 |      l10 |       l11 |       l12 |       l20 |       l21 |       l22 |    tan_f1 |    tan_f2 |   tmin |   tmax |   etype |   PNS |   UNS |   NCN |   nSer |   nSeq |   nJLE |
    |------|--------|---------|-------|----------|---------|------------|---------|----------------|----------|-------------|----------|----------|-------------|-------------|-----------|-----------|--------------|--------------------|-----------------|----------|---------------|---------------|------|-----------|----------|-----------|-----------|-----------|------------|------------|-----------|-----------|-----------|--------|----------|---------|-------|----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|--------|--------|---------|-------|-------|-------|--------|--------|--------|
    | 3367 |   -588 |       2 |    14 | 10:48:38 | 18449.6 |     -32008 |      50 | T              | -0.02436 |     1.0251  | 16.8S    | 99.7E    |   -16.8416  |     99.6838 |      88.5 |     348.1 |         85.3 | 02m29s             |           149.1 |     3368 |           169 |   1.50633e+06 |   11 |  0.109287 | 0.548731 | -4.2e-05  | -8.38e-06 | -0.000409 |  0.123347  |  9.52e-05  | -1.81e-06 | -15.3705  |  0.012714 |  3e-06 | 340.561  | 15.0017 |     0 | 0.544127 |  9.84e-05 | -1.21e-05 | -0.002001 |  9.79e-05 | -1.2e-05  | 0.0046912 | 0.0046678 |     -3 |      3 |       1 |     0 |     0 |     0 |     73 |     35 |     44 |
    | 3369 |   -587 |       2 |     3 | 01:50:15 | 18434.7 |     -31996 |      60 | T              | -0.71738 |     1.05077 | 63.3S    | 110.7W   |   -63.3363  |   -110.722  |      43.9 |     339.4 |        242.4 | 03m21s             |           200.6 |     3370 |           169 |   1.50669e+06 |    2 |  0.213551 | 0.57455  | -3.68e-05 | -9.76e-06 | -0.691366 |  0.0976911 |  0.0001687 | -1.57e-06 | -18.4352  |  0.010715 |  4e-06 | 205.673  | 15.0002 |     0 | 0.536144 |  1.66e-05 | -1.3e-05  | -0.009944 |  1.65e-05 | -1.3e-05  | 0.0047057 | 0.0046823 |     -3 |      3 |       1 |    -1 |     0 |     0 |     72 |     25 |     20 |
    | 3376 |   -585 |       6 |     9 | 02:19:58 | 18398.7 |     -31967 |      47 | T              | -0.42678 |     1.0652  | 2.4S     | 134.2W   |    -2.36187 |   -134.223  |      64.7 |     346.8 |        235   | 06m10s             |           369.6 |     3377 |           169 |   1.50755e+06 |    2 | -0.080431 | 0.557901 |  3.09e-05 | -9.07e-06 | -0.460893 |  0.142135  | -0.000127  | -2.47e-06 |  22.2152  |  0.005645 | -5e-06 | 212.15   | 15.0001 |     0 | 0.533471 |  8.27e-05 | -1.24e-05 | -0.012604 |  8.23e-05 | -1.24e-05 | 0.0045952 | 0.0045723 |     -3 |      3 |       1 |     0 |     0 |     0 |     72 |     41 |      2 |
    | 3378 |   -584 |       5 |    28 | 19:28:50 | 18383.9 |     -31955 |      57 | T              |  0.32013 |     1.07977 | 38.2N    | 45.0W    |    38.1559  |    -45.0206 |      71.1 |     158.5 |        271.5 | 06m04s             |           364.2 |     3379 |           169 |   1.5079e+06  |   19 | -0.366269 | 0.55457  |  6.67e-05 | -9.37e-06 |  0.216789 |  0.18194   | -0.0001464 | -3.25e-06 |  20.3647  |  0.008166 | -4e-06 | 107.344  | 15.0012 |     0 | 0.53022  | -1.2e-05  | -1.28e-05 | -0.015838 | -1.19e-05 | -1.27e-05 | 0.0045938 | 0.0045709 |     -3 |      3 |       1 |     0 |     0 |     0 |     73 |     33 |     40 |
    | 3385 |   -582 |      10 |     1 | 20:59:24 | 18348   |     -31926 |      44 | T              |  0.75199 |     1.0455  | 40.6N    | 31.1W    |    40.5615  |    -31.1115 |      41   |     219.9 |        225.3 | 03m16s             |           195.6 |     3386 |           170 |   1.50876e+06 |   21 |  0.373619 | 0.507081 | -2.4e-05  | -8.46e-06 |  0.652639 | -0.285167  | -3.11e-05  |  4.98e-06 |  -0.74782 | -0.016179 |  0     | 136.9    | 15.004  |     0 | 0.537299 |  3.53e-05 | -1.3e-05  | -0.008795 |  3.51e-05 | -1.3e-05  | 0.0047218 | 0.0046983 |     -3 |      3 |       1 |     1 |     0 |     0 |     72 |     49 |      2 |
    | 3387 |   -581 |       9 |    21 | 13:01:27 | 18333.1 |     -31914 |      54 | T              |  0.06896 |     1.04501 | 7.0N     | 62.0E    |     7.00881 |     61.9918 |      86   |     209   |        150.4 | 03m48s             |           227.6 |     3388 |           170 |   1.50911e+06 |   13 |  0.021509 | 0.503022 | -3.4e-06  | -8.15e-06 |  0.066971 | -0.281069  | -2.75e-05  |  4.77e-06 |   3.51354 | -0.015938 | -1e-06 |  16.1856 | 15.004  |     0 | 0.539017 | -6.54e-05 | -1.27e-05 | -0.007086 | -6.51e-05 | -1.27e-05 | 0.0047085 | 0.004685  |     -3 |      3 |       1 |     0 |     0 |     0 |     74 |     40 |      4 |


Let's load the particular eclipse we want to study:


```python
eclipse = montu.SolarEclipse(eclipses.data.iloc[3])
print(eclipse)
```

    SolarEclipse
    Date (catalogue): -0584-05-28
    Catalogue
      Eclipse type         : T (total)
      γ                    : 0.32013 R⊕
      magnitude            : 1.07977
      julian_date          : 1507900.31200 (JD TT)
      ΔT assumed           : 18383.9 s
      saros                : 57
      luna_num             : -31955
      cat_no               : 3379
    Greatest eclipse
      td_ge (TT)           : 19:28:50
      lat_ge, lng_ge       : 38.2N, 45.0W
      lat_dd_ge            : 38.15594°
      lng_dd_ge            : -45.02063°
      sun_alt, sun_azm     : 71.1°, 158.5°
    Central path
      path_width           : 271.5 km
      central_duration     : 06m04s
      duration_secs        : 364.2 s
      path_map             : http://xjubier.free.fr/en/site_pages/solar_eclipses/xSE_GoogleMap3.php?Ecl=-05840528&Acc=2&Umb=1&Lmt=1&Mag=0



```python
troy = montu.Observer(site="troy")
conditions = eclipse.conditions_eclipse(troy)
conditions.show_details()
```

    Eclipse local circumstances
      Catalogue date       : -0584-05-28 (T, total)
      Observer             : lat 39.957500°, lon 26.238900°, 30 m
      Kind                 : total
      Visible              : yes
      Magnitude            : 1.031
      Obscuration          : 1.000
      Moon/Sun radius ratio: 1.0667
      Sun altitude at max  : 14.77°
      Maximum (UTC)        : -584-05-28 15:56:56
      Maximum (JD UT)      : 1507900.164301
      Maximum (JD TT)      : 1507900.377078
      t_max                : 2.042425 h = 122.545470 min (from catalogue t0)
    Contacts (UTC)
      C1 (first contact)   : -584-05-28 14:57:57 (alt 26.0°, az 275.5°)
      C2 (second contact)  : -584-05-28 15:54:54 (alt 15.1°, az 284.1°)
      C3 (third contact)   : -584-05-28 15:58:58 (alt 14.4°, az 284.7°)
      C4 (fourth contact)  : -584-05-28 16:51:51 (alt 4.8°, az 292.7°)
      Umbra duration       : 00:03:40
      cond_map             : http://xjubier.free.fr/en/site_pages/solar_eclipses/xSE_GoogleMap3.php?Ecl=-05840528&Acc=2&Umb=1&Lmt=1&Mag=0&Lat=39.9575&Lng=26.2389&Elv=30.0&Zoom=9&LC=1


**Xavier Jubier cross-checks** — the `print` output above lists `path_map` and `cond_map` as plain text; use these links to open the interactive maps in a browser:

- [`path_map`](http://xjubier.free.fr/en/site_pages/solar_eclipses/xSE_GoogleMap3.php?Ecl=-05840528&Acc=2&Umb=1&Lmt=1&Mag=0) — central eclipse path for 28 May 585 BCE (catalogue year −584)
- [`cond_map`](http://xjubier.free.fr/en/site_pages/solar_eclipses/xSE_GoogleMap3.php?Ecl=-05840528&Acc=2&Umb=1&Lmt=1&Mag=0&Lat=39.9575&Lng=26.2389&Elv=30.0&Zoom=9&LC=1) — local circumstances at Troy (Ilion), the observer site in this example

Further case studies: [`MontuPython-SolarEclipses.ipynb`](examples/MontuPython-SolarEclipses.ipynb).


### Evolution of polar stars

Choose from database all bright stars that according to [wikipedia](https://en.wikipedia.org/wiki/Pole_star#Precession_of_the_equinoxes) were or will be close to the celestial North pole:


```python
star_names = ('Polaris','Vega','Thuban','Deneb','Alderamin','Kochab')
stars = stars.get_stars(ProperName=star_names)
stars
```




    6 star(s):
    |     |   MN |     HD |   HR |    HIP | Gl     | Name      | OtherDesignations                                                                   | ProperName   | Bayer   | Flamsteed   | Constellation   |   RAJ2000 |   DecJ2000 |   GalLonJ2000 |   GalLatJ2000 |   pmRA |   pmDec |   RadVel |   Distance |   Vmag |   Vmag_min |   Vmag_max |    B-V | SpType       |   Luminosity |     XJ2000 |     YJ2000 |    ZJ2000 |    VXJ2000 |   VYJ2000 |    VZJ2000 |   Primary |   MultipleID |   IsMultiple |   IsVariable |
    |-----|------|--------|------|--------|--------|-----------|-------------------------------------------------------------------------------------|--------------|---------|-------------|-----------------|-----------|------------|---------------|---------------|--------|---------|----------|------------|--------|------------|------------|--------|--------------|--------------|------------|------------|-----------|------------|-----------|------------|-----------|--------------|--------------|--------------|
    |   5 |    6 | 172167 | 7001 |  91262 | Gl 721 | Vega      | 3 Lyr/3Alp Lyr/Gl 721/HD 172167/HIP 91262/HR 7001/HYG 90978/MN 6/Vega/α Lyr         | Vega         | α Lyr   | 3 Lyr       | Lyr             |  18.6156  |    38.7837 |       67.4482 |      19.2374  | 201.02 |  287.46 |    -12.1 |     7.6787 |   0.03 |    nan     |    nan     | -0.001 | A0Vvar       |      49.9344 |   0.960565 |   -5.90801 |   4.80973 |  4.76e-06  | 1.734e-05 |  5.9e-07   |     90979 |          nan |            0 |            0 |
    |  21 |   22 | 197345 | 7924 | 102098 | nan    | Deneb     | 50 Cyg/50Alp Cyg/Deneb/HD 197345/HIP 102098/HR 7924/HYG 101766/MN 22/α Cyg          | Deneb        | α Cyg   | 50 Cyg      | Cyg             |  20.6905  |    45.2803 |       84.2847 |       1.99755 |   1.56 |    1.55 |     -5   |   432.9    |   1.25 |      1.294 |      1.224 |  0.092 | A2Ia         |   51617.9    | 197.251    | -232.113   | 307.601   | -1.33e-06  | 6.62e-06  | -1.34e-06  |    101767 |          nan |            0 |            1 |
    |  48 |   49 |   8890 |  424 |  11767 | nan    | Polaris   | 1 UMi/1Alp UMi/HD 8890/HIP 11767/HR 424/HYG 11734/MN 49/Polaris/α UMi               | Polaris      | α UMi   | 1 UMi       | UMi             |   2.52975 |    89.2641 |      123.28   |      26.4614  |  44.22 |  -11.74 |    -17   |   132.626  |   1.97 |      1.993 |      1.953 |  0.636 | F7:Ib-IIv SB |    2495.74   |   1.3431   |    1.04763 | 132.615   | -1.171e-05 | 2.692e-05 | -1.748e-05 |     11734 |          nan |            0 |            1 |
    |  59 |   60 | 131873 | 5563 |  72607 | nan    | Kochab    | 7 UMi/7Bet UMi/HD 131873/HIP 72607/HR 5563/HYG 72379/Kochab/MN 60/β UMi             | Kochab       | β UMi   | 7 UMi       | UMi             |  14.8451  |    74.1555 |      112.647  |      40.5026  | -32.29 |   11.91 |     17   |    40.1445 |   2.07 |    nan     |    nan     |  1.465 | K4IIIvar     |     208.545  |  -8.05816  |   -7.42971 |  38.6194  | -6.11e-06  | 2.91e-06  |  1.736e-05 |     72380 |          nan |            0 |            0 |
    |  91 |   92 | 203280 | 8162 | 105199 | Gl 826 | Alderamin | 5 Cep/5Alp Cep/Alderamin/Gl 826/HD 203280/HIP 105199/HR 8162/HYG 104860/MN 92/α Cep | Alderamin    | α Cep   | 5 Cep       | Cep             |  21.3096  |    62.5856 |      100.999  |       9.17207 | 149.91 |   48.27 |    -11.5 |    15.0376 |   2.45 |    nan     |    nan     |  0.257 | A7IV-V       |      20.6253 |   5.27611  |   -4.4832  |  13.3488  |  5.7e-07   | 1.386e-05 | -8.82e-06  |    104861 |          nan |            0 |            0 |
    | 351 |  352 | 123299 | 5291 |  68756 | nan    | Thuban    | 11 Dra/11Alp Dra/HD 123299/HIP 68756/HR 5291/HYG 68536/MN 352/Thuban/α Dra          | Thuban       | α Dra   | 11 Dra      | Dra             |  14.0732  |    64.3758 |      110.524  |      50.9587  | -56.52 |   17.19 |    -13   |    92.9368 |   3.67 |    nan     |    nan     | -0.049 | A0III SB     |     256.094  | -34.416    |  -20.7588  |  83.7964  | -2.25e-06  | 2.838e-05 | -8.64e-06  |     68537 |          nan |            0 |            0 |



Now precess the position of all stars from -20 000 to 20 000 years from 2000:


```python
import pandas as pd
import time

now = montu.Time()
epochs = np.linspace(-20000 * montu.YEAR, 20000 * montu.YEAR, 1000)
n = len(epochs)
print(f"Precessing {len(star_names)} stars over {n:,} epochs …")

t0 = time.perf_counter()
rows = []
for dt in epochs:
    past = now + dt
    pstars = stars.where_in_space(at=past)
    row = {"tt": past.tt}
    for star in star_names:
        row[star] = pstars.value_for(star, "DecEpoch")
    rows.append(row)

df = pd.DataFrame(rows)
print(f"Done in {time.perf_counter() - t0:.1f} s.")

```

    Precessing 6 stars over 1,000 epochs …


    Done in 1.5 s.


Now plot declinations as a function of time:


```python
fig,ax = plt.subplots(figsize=(10,6))
for star in star_names:
    ax.plot(df['tt'],df[star],label=star)

ax.legend(loc='lower center',ncol=len(star_names))
ax.set_xlabel("Time [year]")
ax.set_ylabel("Declination [deg]")
ax.axvline(montu.Time().tt,lw=3)
ax.text(0.5,1.01,'Now',ha='center',transform=ax.transAxes)
ax.margins(0)
ax.set_xticks(np.linspace(df['tt'].min(),df['tt'].max(),14))
ax.grid()
montu.Time.set_time_ticks(ax)
montu.Util.montu_mark(ax)
fig.tight_layout()
fig.savefig('gallery/pole-stars.png')
plt.close(fig) # Used only for README generation

```

<p align="center"><img src="https://github.com/seap-udea/MontuPython/blob/main/gallery/pole-stars.png?raw=true" alt="Logo""/></p>

Check date when star is close to the pole:


```python
for star in star_names:
    imax = df[star].argmax()
    mtime = montu.Time(df.iloc[imax].tt).get_readable()
    print(f"Star {star} will be the closest to the pole at {mtime.readable.datespice} (declination {montu.D2S(df.iloc[imax][star])})")
```

    Star Polaris will be the closest to the pole at 2086-08-14 03:12:22.305610 (declination 89:31:55.825)
    Star Vega will be the closest to the pole at 11609 B.C. 08-16 20:03:50.106176 (declination 86:22:03.656)
    Star Thuban will be the closest to the pole at 2800 B.C. 08-19 03:23:38.299200 (declination 89:56:05.009)
    Star Deneb will be the closest to the pole at 14732 B.C. 06-05 15:32:41.798400 (declination 86:57:15.601)
    Star Alderamin will be the closest to the pole at 7532-03-04 21:01:57.002880 (declination 87:58:43.013)
    Star Kochab will be the closest to the pole at 1078 B.C. 05-23 11:30:39.902400 (declination 83:29:32.450)


### Working with Observer Horizons

MontuPython can download real terrain data (Copernicus GLO-30 DEM, 30 m resolution) and compute the true visible horizon profile for any site on Earth. Here we show how the Sun rose through the **Royal Wadi** of ancient Akhetaten (modern Amarna) — the sacred notch that Akhenaten chose as the *akhet* ('horizon') symbol of his new city.


```python
site = montu.Observer(site='amarna')
site.horizon_profile(
    max_dist=30,    # search radius in km
    az_step=0.5,    # azimuth resolution in degrees
    coarse_step=0.1 # radial scan step in km
)
print(site.horizon)
```

    Obtaining horizon profile...


    Horizon for 'Amarna (Akhetaten)'
      Coordinates: lat=27.6444, lon=30.9014, alt=90 m
      Status: computed (720 pts)
      Elevation range: [-0.22°, 1.22°]
      Parameters: max_dist=30 km, az_step=0.5°, coarse_step=0.1 km


The `horizon_profile()` call downloads and caches the DEM tiles automatically. Once computed, `site.horizon.get_elevation(az)` returns the terrain elevation angle at any azimuth.

We can now refine the flat-horizon sunrise time to account for the real terrain:


```python
sun = montu.Sun()
mtime = montu.Time('bce1341-10-21')
sun.conditions_in_sky(at=mtime, observer=site, horizon=True)

mtime_rise_hor = montu.Time(sun.condition.rise_time_hor, format='jd')
print(f"Flat-horizon rise: {site.get_local_time(sun.condition.rise_time)}")
print(f"Terrain-corrected: {site.get_local_time(sun.condition.rise_time_hor)}")
print(f"Rise azimuth (terrain): {sun.condition.rise_az_hor:.2f}°")
```

    Flat-horizon rise: 06:10:12.448
    Terrain-corrected: 06:15:42.721
    Rise azimuth (terrain): 102.94°


The terrain correction shifts the sunrise by a few minutes because the Sun must climb above the wadi walls before it becomes visible.

Finally, we can plot the sky at the horizon-corrected sunrise moment, showing the Sun emerging from the Royal Wadi:


```python
fig = site.horizon.plot_horizon(
    at=mtime_rise_hor,
    az_center=90,
    az_delta=40,
    elev_view=4,
    mag_limit=-1,
    show_asterism=False,
    show_starnames=False,
    show_planets=['Sun'],
)
```



See [`MontuPython-ObserverHorizon.ipynb`](examples/MontuPython-ObserverHorizon.ipynb) for a full walkthrough including multi-site comparisons, Valley of the Kings, the Senenmut temple, and the north shaft of the Khufu pyramid.




```python
montu.Horizon.clean_cache()
```

------------

This package has been designed and written by [Jorge I. Zuluaga](https://jorgezuluaga.github.io) (C) 2023-present
