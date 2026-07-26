###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import ephem as pyephem
import json
import math

###############################################################
# Module constants
###############################################################
_LOCATIONS_DATA = None
_EARTH_RMEAN_KM = None


def _earth_radius_km() -> float:
    """Mean Earth radius [km] from the bundled JPL planetary table."""
    global _EARTH_RMEAN_KM
    if _EARTH_RMEAN_KM is None:
        _EARTH_RMEAN_KM = float(montu.load_planets().loc["Earth", "Rmean"])
    return _EARTH_RMEAN_KM


def _load_locations_data() -> dict:
    """Load the bundled ancient-world sites catalogue (cached)."""
    global _LOCATIONS_DATA
    if _LOCATIONS_DATA is None:
        path = montu.Util._data_path("locations.json", check=True)
        with open(path, encoding="utf-8") as fh:
            _LOCATIONS_DATA = json.load(fh)
    return _LOCATIONS_DATA


def _find_location_entry(site_id: str) -> dict:
    """Return a location dict by ``id`` (case-insensitive)."""
    key = (site_id or "").strip().lower()
    if not key:
        raise ValueError("site id must be a non-empty string")
    for entry in _load_locations_data().get("locations", []):
        if str(entry.get("id", "")).lower() == key:
            return entry
    raise ValueError(f"Unknown observing site: {site_id!r}")


def _apply_location_entry(observer: "Observer", entry: dict) -> None:
    """Copy every catalogue field from *entry* onto *observer*."""
    for key, value in entry.items():
        setattr(observer, key, value)


def _pressure_from_altitude(alt_m: float) -> float:
    """Estimate sea-level-standard pressure [mbar] from elevation [m]."""
    return 1013.25 * math.exp(-alt_m / 8434.5)


###############################################################
# Stars Class
###############################################################
class Observer(object):
    """Observing site with geodetic coordinates and atmospheric conditions.

    Parameters
    ----------
    lon : float, optional
        Geodetic longitude of the observing site [degrees]. East positive.
        Default is 0.
    lat : float, optional
        Geodetic latitude of the observing site [degrees]. North positive.
        Default is 0.
    height : float, optional
        Elevation of the observing site above sea level [km].
        Default is 0.
    site : str, optional
        Predefined site id from :meth:`list` (e.g. ``'memphis'``). When
        given, ``lon``, ``lat``, and ``height`` are taken from the catalogue
        unless explicitly overridden.
    pressure : float, optional
        Atmospheric pressure at the observing site [mbar]. Default is 1013.25
        for a manual site, or the catalogue value when ``site=…`` is used.
    temperature : float, optional
        Air temperature at the observing site [°C]. Default is 15 for a manual
        site, or the catalogue value when ``site=…`` is used.
    relative_humidity : float, optional
        Relative humidity of the air (0–1). Default is 0.
    obswl : float, optional
        Observing wavelength [microns]. Default is 0.6.

    Attributes
    ----------
    lon : float
        Geodetic longitude [degrees].
    lat : float
        Geodetic latitude [degrees].
    height : float
        Elevation [km].
    site_id : str or None
        Catalogue id when created via ``site=…``, else ``None``.
    site_name : str or None
        Display name from the catalogue when ``site=…`` was used.
    id, name, alt_m, region, era, description, name_es, …
        When ``site=…`` is used, every field from the bundled
        ``locations.json`` entry is also exposed as an instance attribute
        (same names as in the catalogue).
    pressure_mbar, temperature_c
        Catalogue entries include mean surface pressure (from ``alt_m``) and
        mean annual air temperature (NASA POWER 1991–2020 climatology, 2 m).
        These become ``pressure`` and ``temperature`` on the observer unless
        overridden explicitly in the constructor.
    pressure : float
        Atmospheric surface pressure [mbar].
    temperature : float
        Surface temperature [°C].
    relative_humidity : float
        Relative humidity (0–1).
    obswl : float
        Observing wavelength [microns].
    site : pyephem.Observer
        Underlying PyEphem observer object used for internal computations.

    Examples
    --------
    Create an observer at Rionegro (Colombia):

    >>> import montu
    >>> rionegro = montu.Observer(lon=-75, lat=6, height=2.5)

    Create an observer at the Great Pyramid of Giza:

    >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)

    Pick a predefined ancient-world site:

    >>> memphis = montu.Observer(site='memphis')
    >>> memphis.region
    'Egypt'
    >>> repr(memphis)
    "Observer('memphis'/'Memphis'/29.845800°, 31.250800°, 25 m/P=1010.25 mbar, T=21.6 °C)"
    >>> montu.Observer.list()[:3]
    ['thebes', 'memphis', 'giza']
    """
    def __init__(self,
                 lon=None, lat=None, height=None,
                 site=None,
                 pressure=None, temperature=None,
                 relative_humidity=0, obswl=0.6,
                 horizon_profile=None):
        """Initialise the Observer; see class docstring for parameter details."""

        self.site_id = None
        self.site_name = None

        if site is not None:
            entry = _find_location_entry(site)
            _apply_location_entry(self, entry)
            self.site_id = entry.get("id", site)
            self.site_name = entry.get("name", site)
            alt_m = float(entry.get("alt_m", 0.0 if height is None else height * 1000.0))
            if lon is None:
                lon = float(entry.get("lon", 0.0))
            if lat is None:
                lat = float(entry.get("lat", 0.0))
            if height is None:
                height = alt_m / 1000.0
            if pressure is None:
                if "pressure_mbar" in entry:
                    pressure = float(entry["pressure_mbar"])
                else:
                    pressure = _pressure_from_altitude(alt_m)
            if temperature is None:
                if "temperature_c" in entry:
                    temperature = float(entry["temperature_c"])

        if lon is None:
            lon = 0.0
        if lat is None:
            lat = 0.0
        if height is None:
            height = 0.0
        if pressure is None:
            pressure = 1013.25
        if temperature is None:
            temperature = 15.0
        self.lon = float(lon)
        self.lat = float(lat)
        self.height = float(height)
        self.pressure = float(pressure)
        self.temperature = float(temperature)
        self.relative_humidity = relative_humidity
        self.obswl = obswl

        # Set pyephem observer
        self.site = pyephem.Observer()
        self.site.lon = str(self.lon)
        self.site.lat = str(self.lat)
        self.site.pressure = self.pressure
        self.site.temp = self.temperature
        self.site.elevation = self.height

        # Horizon profile (populated by horizon_profile())
        self.horizon = None
        
        if horizon_profile is not None:
            self.horizon_profile(**horizon_profile)

    def __repr__(self):
        alt_m = self.height * 1000.0
        coords = f"{self.lat:.6f}°, {self.lon:.6f}°, {alt_m:.0f} m"
        atmosphere = f"P={self.pressure} mbar, T={self.temperature} °C"
        if self.site_id:
            name = self.site_name or self.site_id
            return (
                f"Observer('{self.site_id}'/'{name}'/"
                f"{coords}/{atmosphere})"
            )
        return f"Observer('{coords}/{atmosphere})"

    def __str__(self):
        lines = ["Observer"]
        if self.site_id:
            name = self.site_name or self.site_id
            lines.append(f"  Site: {name} [{self.site_id}]")
        region = getattr(self, "region", "")
        era = getattr(self, "era", "")
        if region or era:
            lines.append(f"  Region: {' · '.join(part for part in (region, era) if part)}")
        alt_m = self.height * 1000.0
        lines.append(
            f"  Coordinates: lat {self.lat:.6f}°, lon {self.lon:.6f}°, "
            f"elevation {alt_m:.0f} m ({self.height:.3f} km)"
        )
        lines.append(
            "  Atmosphere: "
            f"P={self.pressure} mbar, T={self.temperature} °C, "
            f"RH={self.relative_humidity}, λ={self.obswl} μm"
        )
        description = getattr(self, "description", "")
        if description:
            lines.append(f"  Description: {description}")
        return "\n".join(lines)

    @classmethod
    def list(cls, details=False):
        """Return predefined observing sites from the bundled catalogue.

        Parameters
        ----------
        details : bool, optional
            If ``False`` (default), return site ids only.
            If ``True``, return a list of dicts with full metadata
            (``id``, ``name``, ``lat``, ``lon``, ``alt_m``, …).

        Returns
        -------
        list
            Site ids or site detail dicts.

        Examples
        --------
        >>> import montu
        >>> montu.Observer.list()
        ['thebes', 'memphis', ...]
        >>> montu.Observer.list(details=True)[0]['name']
        'Thebes (Luxor)'
        """
        locations = _load_locations_data().get("locations", [])
        if details:
            return [dict(entry) for entry in locations]
        return [entry.get("id", "") for entry in locations]

    def get_local_time(self,mtime,hms=True):
        """Compute local solar time at the observing site.

        Parameters
        ----------
        mtime : montu.Time
            Epoch at which the local time is computed.
        hms : bool, optional
            If ``True`` (default) return a formatted ``HH:MM:SS.sss`` string.
            If ``False`` return the local hour as a decimal float.

        Returns
        -------
        str or float
            Local solar time as a ``HH:MM:SS.sss`` string when ``hms=True``,
            or as a decimal hour when ``hms=False``.

        Examples
        --------
        >>> import montu
        >>> giza = montu.Observer(lon=31.134, lat=29.979, height=0.075)
        >>> mtime = montu.Time('-1000-03-21 06:00:00')
        >>> giza.get_local_time(mtime)
        '08:04:32.400'
        >>> giza.get_local_time(mtime, hms=False)
        8.075666...
        """
        if not isinstance(mtime, montu.Time):
            mtime = montu.Time(mtime, format='jd', scale='utc')
        # Use the UTC Julian Day fraction directly. For ancient dates,
        # readable calendar components may come from a different calendar
        # representation and can shift the clock hour by many hours.
        utc_hour = ((mtime.jed + 0.5) % 1.0) * 24.0
        hour = (utc_hour + self.lon / 15.0) % 24.0
        from montu.util import Util
        return Util.dec2sex(hour) if hms else hour

    def sidereal_time(self, mtime, hms=True):
        """Compute local apparent sidereal time at the observing site.

        Uses the underlying PyEphem observer (same convention as
        :func:`montu.maps._observer_sidereal_time_hours` and
        :meth:`montu.Stars.where_in_sky`).

        Parameters
        ----------
        mtime : montu.Time
            Epoch at which the sidereal time is computed (UTC Julian Day).
        hms : bool, optional
            If ``True`` (default) return a formatted ``HH:MM:SS.sss`` string.
            If ``False`` return the sidereal hour as a decimal float.

        Returns
        -------
        str or float
            Local apparent sidereal time as ``HH:MM:SS.sss`` when ``hms=True``,
            or as decimal hours when ``hms=False``.

        Examples
        --------
        >>> import montu
        >>> thebes = montu.Observer(site='thebes')
        >>> mtime = montu.Time('bce 1500-06-21 12:00:00', calendar='mixed')
        >>> thebes.sidereal_time(mtime)
        '07:14:54.925'
        >>> thebes.sidereal_time(mtime, hms=False)
        7.248590...
        """
        if not isinstance(mtime, montu.Time):
            mtime = montu.Time(mtime, format='jd', scale='utc')
        self.site.date = mtime.jed - montu.PYEPHEM_JD_REF
        hour = float(self.site.sidereal_time() * montu.RAD / 15.0)
        return montu.D2S(hour) if hms else hour

    def distance_to(self, other, units="km"):
        """Great-circle distance to another observing site.

        Parameters
        ----------
        other : Observer
            Second site.
        units : {'km', 'm', 'rad', 'deg'}, optional
            ``'km'`` (default) returns surface distance in kilometres.
            ``'m'`` returns metres.  ``'rad'`` and ``'deg'`` return the
            angular separation on the celestial sphere.

        Returns
        -------
        float
            Distance between the two sites.

        Notes
        -----
        Uses a spherical Earth with mean radius from the bundled planetary
        table (same convention as :func:`montu.Util.haversine_distance`).
        Site elevations are not included.

        Examples
        --------
        >>> import montu
        >>> alexandria = montu.Observer(site='alexandria')
        >>> aswan = montu.Observer(site='aswan')
        >>> round(aswan.distance_to(alexandria))
        843
        """
        if not isinstance(other, Observer):
            raise TypeError("distance_to() expects another Observer instance.")
        arc = montu.Util.haversine_distance(
            self.lat * montu.DEG, self.lon * montu.DEG,
            other.lat * montu.DEG, other.lon * montu.DEG,
        )
        if units == "rad":
            return arc
        if units == "deg":
            return arc * montu.RAD
        radius_km = _earth_radius_km()
        if units == "km":
            return arc * radius_km
        if units == "m":
            return arc * radius_km * 1000.0
        raise ValueError(
            f"Unknown units {units!r}; use 'km', 'm', 'rad', or 'deg'."
        )

    def horizon_profile(self,
                        max_dist:    float = 30.0,
                        az_step:     float = 1.0,
                        coarse_step: float = 3.0,
                        tmpdir:      "str | None" = None,
                        site_name:   str   = None,
                        verbose:     bool  = False) -> "montu.Horizon":
        """Compute the real visible horizon profile for this observing site.

        Downloads the Copernicus GLO-30 DEM (30 m/pixel) automatically and
        runs a two-phase radial scan (coarse + fine refinement) to find the
        maximum elevation angle at each azimuth direction.

        Parameters
        ----------
        max_dist : float
            Maximum distance to search for peaks [km]. Default: 30.
        az_step : float
            Step size in azimuth for the profile [degrees]. Default: 1.
        coarse_step : float
            Radial step size for the initial coarse scan [km]. Default: 3.
        tmpdir : str or None
            Directory to store downloaded DEM tiles.
            If None, uses a 'montu_dem' folder inside the system's temporary directory.
            Default: None.
        site_name : str, optional
            Overrides the site name attached to the Horizon instance.
            If None, uses ``self.name`` if available, else ``'User Site'``.
        verbose : bool
            If True, prints detailed progress messages during the scan.

        Returns
        -------
        self.horizon : montu.Horizon
            The computed horizon object (also stored in ``self.horizon``).

        Examples
        --------
        >>> import montu
        >>> udea = montu.Observer(lat=6.266152, lon=-75.569335, height=1.468)
        >>> udea.horizon_profile()                 # download DEM & compute
        >>> udea.horizon.get_elevation(90)         # elevation looking East
        >>> udea.horizon.plot_horizon()            # interactive Plotly chart

        >>> # Higher resolution, larger area
        >>> udea.horizon_profile(max_dist=40, az_step=0.5, coarse_step=2)
        """
        alt_m = self.height * 1000.0
        
        if site_name is None:
            site_name = self.site_name or self.site_id
            if not site_name:
                site_name = f"MontuSite (lat. {self.lat}, lon. {self.lon}, alt. {alt_m})"
                
        h = montu.Horizon(lat=self.lat, lon=self.lon, alt_m=alt_m,
                          site_name=site_name, observer=self)
        h.get_profile(max_dist=max_dist, az_step=az_step,
                      coarse_step=coarse_step, tmpdir=tmpdir, verbose=verbose)
        self.horizon = h
        return h
