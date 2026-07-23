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

    @staticmethod
    def true_obliquity(at: montu.Time) -> float:
        """True obliquity of the ecliptic [degrees] including nutation."""
        from pyplanets.core.coordinates import true_obliquity
        return float(true_obliquity(at.obj_pyplanet))

    @staticmethod
    def mean_obliquity(at: montu.Time) -> float:
        """Mean obliquity of the ecliptic [degrees] (secular variation only)."""
        from pyplanets.core.coordinates import mean_obliquity
        return float(mean_obliquity(at.obj_pyplanet))

    @staticmethod
    def nutation_longitude(at: montu.Time) -> float:
        """Nutation in longitude (Delta psi) [degrees]."""
        from pyplanets.core.coordinates import nutation_in_longitude
        return float(nutation_in_longitude(at.obj_pyplanet))

    @staticmethod
    def equation_of_time(at: montu.Time) -> float:
        """Equation of Time [minutes]. 
        Apparent Solar Time - Mean Solar Time."""
        import math
        from pyplanets.core.coordinates import ecliptical2equatorial, true_obliquity, nutation_longitude
        from pyplanets.sun import Sun
        from pyplanets.core.epoch import JDE2000
        from pyplanets.core.angle import Angle
        
        epoch = at.obj_pyplanet
        t = (epoch - JDE2000) / 365250
        l0 = (280.4664567 +
              t * (360007.6982779 +
                   t * (0.03032028 +
                        t * (1.0 / 49931.0 +
                             t * (-1.0 / 15300.0 - t * 1.0 / 2000000.0)))))
        l0 = Angle(l0).to_positive()
        
        lon, lat, r = Sun.apparent_geocentric_position(epoch)
        epsilon = true_obliquity(epoch)
        alpha, dec = ecliptical2equatorial(lon, lat, epsilon)
        alpha = alpha.to_positive()
        deltapsi = nutation_longitude(epoch)
        
        # PyPlanets returns Mean - Apparent (l0 - alpha). We want Apparent - Mean.
        e_deg = alpha() - deltapsi() * math.cos(epsilon.rad()) - (l0() - 0.0057183)
        # Wrap to [-180, 180] degrees to avoid wrap-around discontinuities
        e_deg = (e_deg + 180) % 360 - 180
        
        # Convert degrees to minutes of time (1 deg = 4 min)
        return e_deg * 4.0

    @staticmethod
    def greenwich_sidereal_time(at: montu.Time, apparent: bool = True) -> float:
        """Greenwich Sidereal Time [hours].
        If apparent=True (GAST), includes nutation. If False (GMST), does not.
        """
        if apparent:
            from pyplanets.core.coordinates import true_obliquity, nutation_in_longitude
            eps = true_obliquity(at.obj_pyplanet)
            dpsi = nutation_in_longitude(at.obj_pyplanet)
            return float(at.obj_pyplanet.apparent_sidereal_time(eps, dpsi)) * 24.0
        else:
            return float(at.obj_pyplanet.mean_sidereal_time()) * 24.0
