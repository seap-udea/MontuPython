MontuPython API
===============

This is the complete API reference for MontuPython.  Classes are listed
first; standalone utility functions follow at the end of each module
section.

.. contents:: Module Index
   :local:
   :depth: 1

----

montu.time
----------

.. autoclass:: montu.Time
   :members:
   :member-order: bysource

montu.observer
--------------

.. autoclass:: montu.Observer
   :members:
   :member-order: bysource

Celestial Bodies (montu.sebau)
------------------------------

.. autoclass:: montu.Sun
   :members:
   :member-order: bysource

.. autoclass:: montu.Moon
   :members:
   :member-order: bysource

.. autoclass:: montu.Planet
   :members:
   :member-order: bysource

montu.stars
-----------

.. autoclass:: montu.Star
   :members:
   :member-order: bysource

.. autoclass:: montu.Stars
   :members:
   :member-order: bysource

**Standalone functions**

.. autofunction:: montu.stars.load_stellar_catalogue

.. autofunction:: montu.stars.clear_stellar_catalogue_cache

.. autofunction:: montu.stars.stellar_catalogue_cache_status

.. autofunction:: montu.stars.parse_constellation_boundaries

.. autofunction:: montu.stars.parse_constellation_lines

.. autofunction:: montu.stars.parse_constellation_names

.. autofunction:: montu.stars.constellation_data_files

montu.heka
----------

.. autoclass:: montu.Astro
   :members:
   :member-order: bysource

montu.horizon
-------------

.. autoclass:: montu.horizon.Horizon
   :members:
   :member-order: bysource

montu.phenomena
---------------

.. autoclass:: montu.phenomena.HeliacalRise
   :members:
   :member-order: bysource

.. autoclass:: montu.phenomena.SolarEclipses
   :members:
   :member-order: bysource

.. autoclass:: montu.phenomena.SolarEclipse
   :members:
   :member-order: bysource

.. autoclass:: montu.phenomena.EclipseConditions
   :members:
   :member-order: bysource

.. autoclass:: montu.phenomena.Conjunction
   :members:
   :member-order: bysource

.. autoclass:: montu.phenomena.ConjunctionExplorer
   :members:
   :member-order: bysource

**Standalone functions**

.. autofunction:: montu.phenomena.heliacal_rise

montu.maps
----------

.. autofunction:: montu.maps.mercator_sky_map

.. autofunction:: montu.maps.mercator_sky_axes

.. autofunction:: montu.maps.polar_sky_map

.. autofunction:: montu.maps.polar_sky_map_figure

**Utility functions**

.. autofunction:: montu.maps.ra_deg_about_center

.. autofunction:: montu.maps.unwrap_figure_ra_deg

.. autofunction:: montu.maps.local_solar_to_utc_time

.. autofunction:: montu.maps.local_solar_to_utc_datepro

montu.physics
-------------

.. autofunction:: montu.physics.load_planets

montu.util
----------

.. autoclass:: montu.Util
   :members:
   :member-order: bysource

.. autoclass:: montu.Dictobj
   :members:
   :member-order: bysource

**Standalone functions**

.. autofunction:: montu.util.load_historical_dates

.. autofunction:: montu.util.load_historical_conjunctions

.. autofunction:: montu.util.load_historical_solar_eclipses

.. autofunction:: montu.util.historical_conjunctions_by_id

.. autofunction:: montu.util.get_historical_conjunction

.. autofunction:: montu.util.list_historical_conjunctions

.. autofunction:: montu.util.historical_solar_eclipses_by_id

.. autofunction:: montu.util.get_historical_solar_eclipse

.. autofunction:: montu.util.list_historical_solar_eclipses
