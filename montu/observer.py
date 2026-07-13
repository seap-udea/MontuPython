###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import ephem as pyephem

###############################################################
# Module constants
###############################################################

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
    """
    def __init__(self,
                 lon=0,lat=0,height=0,
                 pressure=1013.25,temperature=15,
                 relative_humidity=0,obswl=0.6):
        """Initialise the Observer; see class docstring for parameter details."""
            
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
