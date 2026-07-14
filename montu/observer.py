###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import ephem as pyephem
import json

###############################################################
# Module constants
###############################################################
_LOCATIONS_DATA = None


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
        Atmospheric pressure at the observing site [mbar]. Default is 1013.25.
    temperature : float, optional
        Air temperature at the observing site [°C]. Default is 15.
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
    >>> montu.Observer.list()[:3]
    ['thebes', 'memphis', 'giza']
    """
    def __init__(self,
                 lon=0, lat=0, height=0,
                 site=None,
                 pressure=1013.25, temperature=15,
                 relative_humidity=0, obswl=0.6):
        """Initialise the Observer; see class docstring for parameter details."""

        self.site_id = None
        self.site_name = None

        if site is not None:
            entry = _find_location_entry(site)
            self.site_id = entry.get("id", site)
            self.site_name = entry.get("name", site)
            lon = float(entry.get("lon", lon))
            lat = float(entry.get("lat", lat))
            height = float(entry.get("alt_m", height * 1000.0)) / 1000.0

        # Properties of the site
        self.lon = lon
        self.lat = lat
        self.height = height

        # Atmospheric properties
        self.pressure = pressure
        self.temperature = temperature
        self.relative_humidity = relative_humidity
        self.obswl = obswl

        # Set pyephem observer
        self.site = pyephem.Observer()
        self.site.lon = str(self.lon)
        self.site.lat = str(self.lat)
        self.site.pressure = self.pressure
        self.site.temp = self.temperature
        self.site.elevation = self.height

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
        mtime.get_readable()
        comps = mtime.readable.comps
        hour = (comps[4] + comps[5] / 60.0 + comps[6] / 3600.0 + comps[7] / (1e6*3600.0)) + self.lon / 15
        return montu.D2H(hour) if hms else hour
