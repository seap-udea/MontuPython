###############################################################
# Montu interdependencies
###############################################################
import montu

###############################################################
# Required packages
###############################################################
import numpy as np


class Astro(object):
    """Low-level astronomical coordinate transforms."""

    @staticmethod
    def where_in_sky(RA=0, Dec=0, at=None, observer=None):
        """Compute horizontal coordinates (azimuth, elevation) for RA/Dec.

        Parameters
        ----------
        RA : float, optional
            Right ascension [hours]. Default is 0.
        Dec : float, optional
            Declination [degrees]. Default is 0.
        at : montu.Time, optional
            Epoch of the observation. Defaults to the current time.
        observer : montu.Observer
            Observing site.

        Returns
        -------
        az : float
            Azimuth [degrees], measured from North through East.
        el : float
            Elevation above the horizon [degrees].

        Raises
        ------
        ValueError
            If *observer* is not a :class:`montu.Observer` instance.

        Examples
        --------
        >>> import montu
        >>> rionegro = montu.Observer(lon=-75, lat=6, height=2.5)
        >>> mtime = montu.Time("2024-05-01 19:00:00")
        >>> az, el = montu.Astro.where_in_sky(
        ...     RA=6.770358, Dec=-16.751203,
        ...     at=mtime, observer=rionegro)
        >>> montu.Util.dec2sex(az), montu.Util.dec2sex(el)
        ('...', '...')
        """
        if at is None:
            at = montu.Time()

        if not isinstance(observer, montu.Observer):
            raise ValueError("You must provide a valid montu.Observer")

        observer.site.date = at.jed - montu.PYEPHEM_JD_REF
        ltst = observer.site.sidereal_time() * montu.RAD / 15

        HA = np.mod(ltst - RA, 24)

        lat = observer.lat
        el = np.arcsin(
            np.sin(Dec * montu.DEG) * np.sin(lat * montu.DEG)
            + np.cos(Dec * montu.DEG) * np.cos(lat * montu.DEG) * np.cos(HA * 15 * montu.DEG)
        ) * montu.RAD
        az = np.arctan2(
            -np.sin(HA * 15 * montu.DEG) * np.cos(Dec * montu.DEG) / np.cos(el * montu.DEG),
            (np.sin(Dec * montu.DEG) - np.sin(lat * montu.DEG) * np.sin(el * montu.DEG))
            / (np.cos(lat * montu.DEG) * np.cos(el * montu.DEG)),
        ) * montu.RAD
        az = np.mod(az, 360)

        return az, el
